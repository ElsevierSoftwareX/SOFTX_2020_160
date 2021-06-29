#!/usr/bin/env python

# Copyright (C) 2019 Marco Drago, Igor Yakushin, Sergey Klimenko
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import cWB_conf
import commands, sys, os, glob
import os.path
if (os.environ['SITE_CLUSTER']=="CASCINA"):
  from glue.segments import *
  from glue.segmentsUtils import *
else:
  from ligo.segments import *
  from ligo.segments.utils import *
from run_utils import *
from collections import deque
import getpass

class JOB:
    def __init__(self,start,end,basedir):
        self.start=start
        self.end=end
        self.basedir=basedir
        self.dir="%s/%s/%d-%d"%(self.basedir, str(self.start)[:5],self.start, self.end)
        a=commands.getstatusoutput("mkdir -p %s"%(self.dir))
        self.missing_data=[]
    def framelist(self):
        input_dir="%s/input"%(self.dir)
        #commands.getstatusoutput("mkdir -p %s"%(input_dir))
        os.chdir(input_dir)
        #f=open("burst.in","w")
        #print >>f,"%d %d"%(int(self.start)-cWB_conf.job_offset, int(self.end)+cWB_conf.job_offset)
        #f.close()
        #for ifo in conf.ifos:
        for ifo in range(len(cWB_conf.ifos)):
            count_qm=cWB_conf.bkg_dir[ifo].count("?")
            string_qm="?"*count_qm
            fn="%s_burst.in"%(cWB_conf.ifos[ifo])
            f=open(fn,"w")
            print >>f,"%d %d"%(int(self.start)-cWB_conf.job_offset, int(self.end)+cWB_conf.job_offset)
            f.close()
            fn="%s.frames"%(cWB_conf.ifos[ifo])
            f=open(fn,"w")
            #pattern=glob.glob("%s*.gwf"%(cWB_conf.bkg_dir[ifo].replace("?????","%s"%(str(self.start)[:5]))))
	    #print repr(pattern)
            #if (int(str(self.start)[:5])<int(str(self.end)[:5]) and cWB_conf.bkg_dir[ifo].find("?????")!=-1):
            #	pattern+=glob.glob("%s*.gwf"%(cWB_conf.bkg_dir[ifo].replace("?????","%s"%(str(self.end)[:5]))))
            if (count_qm==0):
              pattern=glob.glob("%s*.gwf"%cWB_conf.bkg_dir[ifo])
            else:
              pattern=glob.glob("%s*.gwf"%(cWB_conf.bkg_dir[ifo].replace(string_qm,"%s"%(str(self.start)[:count_qm]))))
              if (int(str(self.start)[:count_qm])<int(str(self.end)[:count_qm])):
                pattern+=glob.glob("%s*.gwf"%(cWB_conf.bkg_dir[ifo].replace(string_qm,"%s"%(str(self.end)[:count_qm]))))
	    #print len(pattern)
            print >>f,"\n".join(pattern)
            f.close()
        fs=cWB_conf.bkg_framesize
        self.fstart=(self.start - cWB_conf.job_offset)/fs*fs
        self.fend=(self.end+cWB_conf.job_offset)/fs*fs
        if(self.fend<self.end+cWB_conf.job_offset):
            self.fend+=fs
        #for i,ifo in zip(range(len(conf.ifos)),conf.ifos):
        for i,ifo in zip(range(len(cWB_conf.ifos)),cWB_conf.ifos):
            count_qm=cWB_conf.bkg_dir[i].count("?")
            string_qm="?"*count_qm
            starts=range(self.fstart,self.fend,fs)
            for j,s in zip(range(len(starts)),starts):
                #fn="%s%d-%d.gwf"%(cWB_conf.bkg_dir[i].replace("?????","%s"%(str(s)[:5])),s,fs)
                if (count_qm==0):
                  fn="%s%d-%d.gwf"%(cWB_conf.bkg_dir[i],s,fs)
                else:
                  fn="%s%d-%d.gwf"%(cWB_conf.bkg_dir[i].replace(string_qm,"%s"%(str(s)[:count_qm])),s,fs)
                print ">> %s"%fn
                if(not os.path.exists(fn)):
                    print "missing?"
                    self.missing_data.append(fn)
        if(len(self.missing_data)>0):
            f=open("%s/input/missing_data.txt"%self.dir,"w")
            print >>f,"\n".join(self.missing_data)
            f.close()

    def setup(self,superlag=""):
        try:
            condor_requirements_lines=open(cWB_conf.condor_requirements_file).readlines()
            condor_requirements="".join(condor_requirements_lines)
        except:
            condor_requirements=""
        try:
            accounting_group_user="accounting_group_user=%s"%(cWB_conf.accounting_group_user)
        except:
            accounting_group_user=""
        os.chdir(self.dir)
        #a=commands.getstatusoutput("mkdir -p input condor OUTPUT OUTPUT_CED tmp tmp_ced")
        a=commands.getstatusoutput("mkdir -p input condor output tmp log")
        #if (os.environ['SITE_CLUSTER']=="CASCINA"):
        #  commands.getstatusoutput("ln -s %s .rootlogon.C"%os.environ['CWB_ROOTLOGON_FILE'])
        if (superlag==""):
            self.framelist()
            superlag="0"
            for j in range(1,len(cWB_conf.ifos)):
                superlag+=",0"
        os.chdir(self.dir)
        f=open("condor/%s.sub"%(cWB_conf.label),"w")
        lagsize=cWB_conf.bkg_nlags
        lagoff=1;
        tmpfile=cWB_conf.bkg_par
        tmpfile=tmpfile.replace(".C","_split.C")
        if (cWB_conf.bkg_split>1):
          environment="CWB_JOBID=$(JOB);CWB_UFILE=$(CWB_UFILE);CWB_STAGE=$(CWB_STAGE);SLAG_SHIFT=$(SLAG_SHIFT)"
        else:
          environment="PID=1;CWB_JOBID=1;CWB_UFILE=%s;CWB_STAGE=CWB_STAGE_FULL;SLAG_SHIFT=%s"%(tmpfile,superlag)
        sub="""universe = vanilla
getenv = true
request_memory = 3000
executable = %s/tools/online/bin/net_bkg.csh
environment = %s
output = %s/%s_$(PID)_%s_%d_%d.out
error = %s/%s_$(PID)_%s_%d_%d.err
log = %s/%s_%s_%d_%d.log
notification = never
rank=memory
accounting_group = %s
%s
%s
queue
""" % (os.environ['HOME_WAT'], environment, cWB_conf.log_path, cWB_conf.label, cWB_conf.label, self.start, self.end, cWB_conf.log_path, cWB_conf.label, cWB_conf.label, self.start, self.end, cWB_conf.log_path, cWB_conf.label, cWB_conf.label, self.start, self.end,cWB_conf.accounting_group,accounting_group_user,condor_requirements)
        if (os.environ['SITE_CLUSTER']=="CASCINA"):
          sub="""universe = vanilla
getenv = true
request_memory = 3000
executable = %s/tools/online/bin/net_bkg.csh
environment = %s
Initialdir = %s/log
output = %s_$(PID)_%s_%d_%d.out
error = %s_$(PID)_%s_%d_%d.err
log = %s_%s_%d_%d.log
notification = never
rank=memory
accounting_group = %s
%s
%s
queue
""" % (os.environ['HOME_WAT'], environment, self.dir, cWB_conf.label, cWB_conf.label, self.start, self.end, cWB_conf.label, cWB_conf.label, self.start, self.end, cWB_conf.label, cWB_conf.label, self.start, self.end,cWB_conf.accounting_group,accounting_group_user,condor_requirements)
        print >>f,sub
        f.flush()
        f.close()
        if (cWB_conf.bkg_split>1):
          f=open("condor/%s.dag"%(cWB_conf.label),"w")
          for l in range(1,cWB_conf.bkg_split+1):
            lagsize=cWB_conf.bkg_nlags/cWB_conf.bkg_split
            lagoff=lagsize*(l-1)+1
            print >>f,"""JOB A%d %s/condor/%s.sub
VARS A%d PID="%d" JOB="%d" CWB_UFILE="%s" CWB_STAGE="CWB_STAGE_FULL" SLAG_SHIFT="%s"
RETRY A%d 3000"""%(l,self.dir,cWB_conf.label,l,l,1,tmpfile,superlag,l)
          f.flush()
          f.close()

    def launch(self):
        #print "launch"
        if(os.path.exists("%s/input/missing_data.txt"%self.dir)):
           return
        os.chdir("%s/condor"%self.dir)
        if (cWB_conf.bkg_split>1):
          com="condor_submit_dag %s.dag"%(cWB_conf.label)
        else:
          com="condor_submit %s.sub"%(cWB_conf.label)
        #print com
        a=commands.getstatusoutput(com)
        #print a
        com="touch ../start"
        a=commands.getstatusoutput(com)
class TS:
    def __init__(self):
        self.considered_segments=segmentlist([])
        self.processed_segments=segmentlist([])
        self.running_segments=segmentlist([])
        self.missing_segments=segmentlist([])
        self.job_segments=segmentlist([])
        self.run_segments=segmentlist([])
        self.job_queue=[]
    def update_considered(self):
        #print "In update_considered()"
        found=0
        try:
            f=open(considered_segments_file)
            #self.considered_segments=fromsegwizard(f)
            old=fromsegwizard(f)
            f.close()
            found+=1
        except Exception, e:
            #print "cannot read considered segments file %s\n%s"%(considered_segments_file,e)
            #self.considered_segments=segmentlist([])                        
            old=segmentlist([])                        
        try:
            f=open("%s/%s/seg_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir))
            extra=fromsegwizard(f)
            f.close()
            found+=1
        except:
            extra=segmentlist([])
        if (found>0):
          considered=old | extra; considered.coalesce()
          f=open(considered_segments_file,"w")
          tosegwizard(f,considered)
          f.close()
          self.considered_segments=considered

    def update_processed(self):
        a=glob.glob("%s/%s*/?????/*"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir))
        processed=[]
        #print "In update_processed()"
        for b in a:
            if(os.path.exists("%s/start"%b) and not os.path.exists("%s/finished"%b)):
              c=glob.glob("%s/output/*.root"%b)
              if(len(c)==cWB_conf.bkg_split):
                commands.getstatusoutput("touch %s/finished"%b)
                if (not os.path.exists("%s/log.tgz"%b)):
                  commands.getstatusoutput("tar -czf %s/log.tgz %s/log; rm -rf %s/log"%(b,b,b))
            if(os.path.exists("%s/start"%b) and os.path.exists("%s/finished"%b)):
                processed.append(b)
        try:
            self.processed_segments=segmentlist(map(lambda x: segment(map(int,x.split("/")[-1].split("-"))),processed))
            self.processed_segments.coalesce()
        except Exception,e:
            print e
            #print len(processed)
            self.processed_segments=segmentlist([])
        f=open(processed_segments_file,"w")
        tosegwizard(f,self.processed_segments)
        f.close()
    def update_missing(self):
        #print "In update_missing()"
        try:
            f=open(missing_segments_file)
            self.missing_segments=fromsegwizard(f)
            f.close()
        except Exception, e:
            print "cannot read missing segments file %s\n%s"%(missing_segments_file,e)
            self.missing_segments=segmentlist([])                                
    def update_running(self):
        #print "In update_running()"
        a=glob.glob("%s/%s*/?????/*"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir))
        running=[]
        for b in a:
            #print b
            if(os.path.exists("%s/start"%b) and not os.path.exists("%s/finished"%b)):
                running.append(b)
        #print running
        try:
            self.running_segments=segmentlist(map(lambda x: segment(map(int,x.split("/")[-1].split("-"))),running))
            self.running_segments.coalesce()
        except Exception,e:
            print e
            print len(running)
            self.running_segments=segmentlist([])
        f=open(running_segments_file,"w")
        tosegwizard(f,self.running_segments)
        f.close()
    def break_seg(self):
        #segments_dir="%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.seg_dir)
        #run_segments_file="%s/%s"%(segments_dir,cWB_conf.run_segments_file)
        f=open(run_segments_file)
        lines = f.readlines()
        f.close()
        f=open(job_segments_file,"w")
        for line in lines:
           #print line,
           if(line.find("#")==-1):
               seg=line.split("\t")
               start=int(seg[1])+cWB_conf.job_offset
               stop=int(seg[2])-cWB_conf.job_offset
               n = int(stop-start)/cWB_conf.bkg_job_duration
               #print "%d %d %d"%(start,stop,n)
               s=start
               e=stop
               if (n>1):
                   for i in range(0,n-1):
                       s = start+ i*cWB_conf.bkg_job_duration
                       e = start+(i+1)*cWB_conf.bkg_job_duration
                       print>>f, "%d %d"%(s,e)
                   s = e
                   e = stop
               n = int(e-s)/cWB_conf.bkg_job_duration
               if (n>1):
                   print>>f, "error"
               if (n==0):
                   if (e-s>cWB_conf.bkg_job_minimum):
                       print>>f, "%d %d"%(s,e)
               if (n==1):
                   if (cWB_conf.bkg_job_minimum==cWB_conf.bkg_job_duration):
                       e = s+ cWB_conf.bkg_job_duration
                       print>>f, "%d %d"%(s,e)
                   else:
                       half = (e-s)/2
                       if (half>cWB_conf.bkg_job_minimum):
                            e = s+ half
                            print>>f, "%d %d"%(s,e)
                            s = e
                            e = s+ half
                            print>>f, "%d %d"%(s,e)
                       else:
                            e = s+cWB_conf.bkg_job_duration
                            print>>f, "%d %d"%(s,e)
        f.close()
    def compute_job_segments(self):
        #print "In compute_job_segments"
        rs = self.considered_segments - self.processed_segments; rs.coalesce();
        rs = rs - self.running_segments; rs.coalesce()
        rs = rs - self.missing_segments; rs.coalesce()
        if(abs(rs)==0):
            return
        current_time=get_time()
        self.run_segments=segmentlist(rs)
        #self.run_segments=segmentlist(rs[:-1])
        #n=(abs(rs[-1])-2*cWB_conf.job_offset)/cWB_conf.bkg_job_duration
        #if(current_time - rs[-1][1] > cWB_conf.bkg_delay):
        #    self.run_segments.append(rs[-1])
        #elif(n-2>0):
        #    self.run_segments.append(segment(rs[-1][0],rs[-1][0]+(n-2)*cWB_conf.bkg_job_duration+2*cWB_conf.job_offset))
        self.run_segments.coalesce()
        self.run_segments = segmentlist(filter(lambda x: abs(x)>=cWB_conf.bkg_job_duration+2*cWB_conf.job_offset, self.run_segments))
        f=open(run_segments_file,"w")
        tosegwizard(f,self.run_segments)
        f.close()
        self.break_seg()
        #com="./break.py > %s"%(job_segments_file)
        #print com
        #a=commands.getstatusoutput(com)
        f=open(job_segments_file)
        self.job_segments=map(lambda x: map(int, x.strip().split(" ")), f.readlines())
        f.close()
        self.job_queue=map(lambda x:JOB(x[0],x[1],job_dir),self.job_segments)
        missing_dirs=glob.glob("%s/%s/*/*/input/missing_data.txt"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir))
        for m in missing_dirs:
           #print m
           mis_files=open(m).readlines()
           leng=len(mis_files)
           test=0
           for l in mis_files:
             if (os.path.exists(l.replace("\n",""))):
               #print "%s exists"%l
               test+=1
           if (len(mis_files)==test):
              self.job_queue.append(JOB(int(m.split("/")[-3].split("-")[0]),int(m.split("/")[-3].split("-")[1]),job_dir))
              #commands.getstatusoutput("rm -rf %s"%("/".join(m.split("/")[:-2])))
              commands.getstatusoutput("rm %s"%(m))
        #print repr(self.job_queue)
    def launch(self):
        #print "launch"
        missing_jobs=[]
        for job in self.job_queue:
            job.setup()
            if(len(job.missing_data)>0):
                missing_jobs.append(job)
            else:
                job.launch()
        if(len(missing_jobs)>0):
            missing_segments=segmentlist(map(lambda x: segment(x.start - cWB_conf.job_offset, x.end + cWB_conf.job_offset), missing_jobs)); self.missing_segments.coalesce()
            try:
                f=open(missing_segments_file)
                seg=fromsegwizard(f)
                f.close()
            except Exception,e:
                print e
                seg=segmentlist([])
            missing_segments = missing_segments | seg; missing_segments.coalesce()
            f=open(missing_segments_file,"w")
            tosegwizard(f,missing_segments)
            f.close()
    def superlag(self):
        #print "In superlag"
        started_dir=glob.glob("%s/%s/?????/*/start"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir))
        started_dir.sort()
        #print len(started_dir)
        for d in range(0,len(started_dir)):
          #print "="*30
          #print "%i %s"%(d,started_dir[d][:-6])
          #print "="*30
          for i in range(1,len(bkg_superlag_index)):
            #print "%i %i %i"%(bkg_superlag_index[i],bkg_superlag_list[0][i], bkg_superlag_list[1][i])
            ok=True
            for j in range(0,len(cWB_conf.ifos)):
               if (d+bkg_superlag_list[j][i]<0 or d+bkg_superlag_list[j][i]>=len(started_dir)):
                 ok=False
            if (ok==True):
               #print "-"*30
               #print "%i %i %i"%(bkg_superlag_index[i],bkg_superlag_list[0][i], bkg_superlag_list[1][i])    
               second_start_job=int("/".join(started_dir[d+bkg_superlag_list[0][i]].split("/")[-2:-1]).split("-")[0])
               second_end_job=int("/".join(started_dir[d+bkg_superlag_list[0][i]].split("/")[-2:-1]).split("-")[1])
               second_dir="%s"%("/".join(started_dir[d+bkg_superlag_list[0][i]].split("/")[-3:-1]))
               new_dir="%s/%s_%i"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir,i)
               #print new_dir
               superlag_string="0"
               slag_string="%s"%(bkg_superlag_list[0][i])
               for j in range(1,len(cWB_conf.ifos)):
                  temp_start_job=int("/".join(started_dir[d+bkg_superlag_list[j][i]].split("/")[-2:-1]).split("-")[0])
                  superlag_string+=",%i"%(second_start_job-temp_start_job)
                  slag_string+=",%i"%(bkg_superlag_list[j][i])
               slag_string+=",%i"%(bkg_superlag_index[i])
               if (not os.path.exists("%s/%s/start"%(new_dir,second_dir))):
                  #print superlag_string
                  sjob=JOB(second_start_job,second_end_job,new_dir)
                  sjob.setup("%s"%(superlag_string))
                  for j in range(0,len(cWB_conf.ifos)):
                     temp_dir="/".join(started_dir[d+bkg_superlag_list[j][i]].split("/")[:-1])
                     com="ln -s %s/input/%s.frames %s/%s/input/."%(temp_dir,cWB_conf.ifos[j],new_dir,second_dir)
                     commands.getstatusoutput(com)
                     com="ln -s %s/input/%s_burst.in %s/%s/input/."%(temp_dir,cWB_conf.ifos[j],new_dir,second_dir)
                     commands.getstatusoutput(com)
                  sjob.launch()
               if (os.path.exists("%s/%s/finished"%(new_dir,second_dir)) and not os.path.exists("%s/%s/copied"%(new_dir,second_dir))):
                  files=glob.glob("%s/%s/output/*.root"%(new_dir,second_dir))
                  if(len(files)==cWB_conf.bkg_split):
                    for file in files:
                       startfile=file.split("/")[len(file.split("/"))-1]
                       subfile=startfile.replace("slag0","slag%i"%i)
                       com=("echo %s/%s/output/%s %s/%s/output/%s %s %s | root -b -l %s/tools/online/bin/AdjustSlag.C"%(new_dir,second_dir,startfile,job_dir,second_dir,subfile,slag_string,superlag_string,os.environ['HOME_WAT']))
                       if (os.environ['SITE_CLUSTER']=="CASCINA"):
                         com=("echo %s/%s/output/%s %s/%s/output/%s %s %s | root -b -l -n %s %s/tools/online/bin/AdjustSlag.C"%(new_dir,second_dir,startfile,job_dir,second_dir,subfile,slag_string,superlag_string,os.environ['CWB_ROOTLOGON_FILE'],os.environ['HOME_WAT']))
                       print com
                       commands.getstatusoutput(com)
                       com="rm  %s/%s/output/%s"%(new_dir,second_dir,startfile)
                       commands.getstatusoutput(com)
                    commands.getstatusoutput("touch %s/%s/copied"%(new_dir,second_dir))
    def cycle(self):
        print "================================================ Cycle %d ================================================="%get_time()
        #commands.getstatusoutput("condor_release %s"%cWB_conf.user)
        self.release()
        self.update_running()
        self.update_processed()
        self.update_missing()
        self.update_considered()
        self.compute_job_segments()
        self.launch()
        if (cWB_conf.bkg_njobs>1):
           self.superlag()
    def release(self):
        user=getpass.getuser()
        com="condor_q %s -nobatch | grep 'H ' "%user
        a=commands.getstatusoutput(com)
        jobs=a[1].split("\n")
	if jobs[0]!='':
         for job in jobs:
           j_id=job.split(" ")[0]
           com="condor_q -long %s | grep -w 'HoldReason ='"%j_id
           b=commands.getstatusoutput(com)
           if (b[1].find("Job has gone over memory limit of")!=-1):
              memorylimit=b[1].split(" ")[len(b[1].split(" "))-2]
              mem=int(memorylimit)+1000
              com="condor_qedit %s RequestMemory %s"%(j_id,mem)
              commands.getstatusoutput(com)
           else:
              print "Job %s in held not for memory limit"%j_id
           com="condor_release %s"%j_id
           commands.getstatusoutput(com)
    def run(self):
        stop=False
        while(True and not stop):
            try:
                self.cycle()
            except Exception,e:
                print "Crashed cycle"
                print e
                sys.stdout.flush()
                sys.stderr.flush()
            stop=os.path.exists("%s/sHuTdOwN"%cWB_conf.bkg_run_dir)
        
if(__name__=="__main__"):
    segments_dir="%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.seg_dir)
    job_dir="%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir)
    considered_segments_file="%s/%s"%(segments_dir,cWB_conf.considered_segments_file)
    processed_segments_file="%s/%s"%(segments_dir,cWB_conf.processed_segments_file)
    running_segments_file="%s/%s"%(segments_dir,cWB_conf.running_segments_file)
    missing_segments_file="%s/%s"%(segments_dir,cWB_conf.missing_segments_file)
    run_segments_file="%s/%s"%(segments_dir,cWB_conf.run_segments_file)
    job_segments_file="%s/%s"%(segments_dir,cWB_conf.job_segments_file)
    f=open("%s/superlaglist.txt"%cWB_conf.config_dir)
    lines=map(lambda x: x.replace("\n",""),filter(lambda x: x.find("#")!=0, f.readlines()))
    f.close()
    #print repr(lines)
    bkg_superlag_index=deque([])
    bkg_superlag_list={}
    for i in range(len(cWB_conf.ifos)):
        bkg_superlag_list[i]=deque([])
    for line in lines:
      #print line
      arrays=line.split(" ")
      #print repr(arrays)
      first=False
      index=0
      for array in arrays:
        if (array!=""):
           if (first==False):
             bkg_superlag_index.append(int(array))
             first=True
           else:
             bkg_superlag_list[index].append(int(array))
             index=index+1
    print repr(bkg_superlag_index)
    for i in range(len(cWB_conf.ifos)):
      print repr(bkg_superlag_list[i])
   
    ts=TS()
    ts.run()

