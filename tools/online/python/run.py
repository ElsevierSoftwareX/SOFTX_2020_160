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

import commands, os, sys, glob
if (os.environ['SITE_CLUSTER']=="CASCINA"):
  from glue.segments import *
  from glue.segmentsUtils import *
else:
  from ligo.segments import *
  from ligo.segments.utils import *
from collections import deque
import pickle
import cWB_conf
from Trigger import *
from run_utils import *
from ligo.gracedb.rest import GraceDb, HTTPError
from lal import LIGOTimeGPS

class JOB:
    def __init__(self,start,end,launch_time=-1,status=0):
        self.dir="%s/%s/%s/%d-%d"%(cWB_conf.run_dir,cWB_conf.jobs_dir,str(start)[:6],start,end)
        self.link="%s/%s/%s/%d-%d"%(cWB_conf.web_link,cWB_conf.jobs_dir,str(start)[:6],start,end)
        print self.dir
        if (status!=2):
            self.make_dirs()
        self.start=start
        self.end=end
        self.back=start
        self.launch_time=launch_time
        self.queue_time=-1
        self.completion_time=-1
        self.triggers_dump_time=-1
        self.triggers_send_time=-1        
        self.status=status # 0 - not launched yet, 1 running, 2 finished successfully, 3 - failed, 4 - missing frames
        #self.frames={}
        self.delay_start_ready=-1
        self.delay_start_launch=-1
        self.delay_launch_completion=-1
        self.delay_start_send=-1
        self.triggers_all=[]
        if (status!=2):
           self.dump_job()
        #com="top -b -n 1"
        #a=commands.getstatusoutput(com)
        #self.top=a[1]
        self.cat_segs=[-1]*5
        for cat in range(5):
            self.cat_segs[cat]={}

    def __repr__(self):
        msg="""
        dir=%s, start=%d, end=%d, launch_time=%d, completion_time=%d, status=%s
        """ % (self.dir,self.start,self.end,self.launch_time,self.completion_time,self.interpret_status())
        return msg
    def interpret_status(self):
        if(self.status==0):
            return "Not launched yet"
        if(self.status==1):
            return "Running since %d"%(self.launch_time)
        if(self.status==2):
            return "Completed successfully at %d"%(self.completion_time)
        if(self.status==3):
            return "Failed at %d"%(self.completion_time)

    def make_dirs(self):
        #subdirs=["input","OUTPUT","OUTPUT_CED","tmp","tmp_ced","config","OUTPUT_PE","tmp_pe"]
        subdirs=["input","OUTPUT","tmp","config"]
        for subdir in subdirs:
           com="mkdir -p %s/%s"%(self.dir,subdir)
           a=commands.getstatusoutput(com)
           #print a
        #com="ln -s %s/plugins %s/plugins"%(os.environ['HOME_CWB'],self.dir)
        #a=commands.getstatusoutput(com)
        #print a
        
    def generate_inputs(self):
        ifo_to_index={}
        for i in range(len(cWB_conf.ifos)):
           ifo_to_index[cWB_conf.ifos[i]]=i
        for ifo in cWB_conf.ifos:
            infile="%s/input/%s.frames"%(self.dir,ifo)

            #ifiles=glob.glob("%s/%s/*"%(cWB_conf.frames_dir[cWB_conf.ifo_to_index[ifo]],ifo))
	    ifiles=glob.glob("%s/*.gwf"%(cWB_conf.frames_dir[ifo_to_index[ifo]]))

            ifiles.sort()
            f=open(infile,"w")
            print >>f,"\n".join(ifiles)
            f.flush()
            f.close()
            #temp_frames=map(lambda x: FRAME(x), ifiles)
            #self.frames[ifo]=map(lambda x: not(x.end<self.start-cWB_conf.job_offset and x.start>self.end+cWB_conf.job_offset),temp_frames)
            tmp_frames=map(lambda x: FRAME(x), ifiles)
            n=len(tmp_frames)
            if(n<=0):
                print "In job.configure(): n=%d ifo=%s"%(n,ifo)
                self.status=4
                return
            coverage=segmentlist(map(lambda x: segment(x.start,x.end),tmp_frames));coverage.coalesce()
            m=segmentlist([segment(self.start-cWB_conf.job_offset,self.end+cWB_conf.job_offset)])-coverage;m.coalesce()
            if(abs(m)>0):
                print "In job.configure(): m=%s"%(repr(m))
                self.status=4
                return            
            cat1=segmentlist([segment(int(self.back)-cWB_conf.job_offset,int(self.end)+cWB_conf.job_offset)]);cat1.coalesce()
            f=open("%s/input/%s_cat2.in"%(self.dir,ifo),"w")
            try:
               fin=open("%s/%s/%s_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,ifo,"CAT2"))
               cat2_temp_seg=fromsegwizard(fin,LIGOTimeGPS)
               fin.close()
               if (abs(cat2_temp_seg)>0):
                 cat2_temp_seg&=cat1
                 for s in cat2_temp_seg:
                   print>>f,"%s %s"%(s[0],s[1])
               else:
                 print >>f,"%s %s"%(int(self.back)-cWB_conf.job_offset,int(self.end)+cWB_conf.job_offset)
            except:
               print >>f,"%s %s"%(int(self.back)-cWB_conf.job_offset,int(self.end)+cWB_conf.job_offset)
            f.close()
        f=open("%s/input/burst.in"%(self.dir),"w")
        print >>f,"%s %s"%(int(self.back)-cWB_conf.job_offset,int(self.end)+cWB_conf.job_offset)
        f.close()

    def configure(self,allframes):
        self.generate_inputs()

    def launch(self):
        print "In launch: dir=%s"%(self.dir)
        os.chdir(self.dir)
        com="%s/tools/online/bin/net_zero.csh %s"%(os.environ['HOME_WAT'],cWB_conf.zerolag_par)
        #print "com=%s"%(com)
        #print "run_dir=%s"%(cWB_conf.run_dir)
        #print "run_dir=%s"%(cWB_conf.online)
        a=commands.getstatusoutput(com)
        if(cWB_conf.debug==1):
            print "Job launched"
            print a
    def publish(self,segs,jobs_published):
        #print "Before merge_triggers()"; sys.stdout.flush()
        #self.merge_triggers()
        print "Before dump_triggers()"; sys.stdout.flush()
        self.dump_triggers(segs)
        self.check_triggers(jobs_published)
        self.send_to_external_collaborations()
        print "before web_page()"; sys.stdout.flush()        
        self.extract_triggers()
        print "before web_page()"; sys.stdout.flush()        
        self.web_page()
        print "after web_page()"; sys.stdout.flush()        
        self.dump_job()
    def merge_triggers(self):
        com="mkdir -p %s/OUTPUT.merged"%(self.dir)
        a=commands.getstatusoutput(com)
        print a
        a=glob.glob("%s/OUTPUT/*.root"%self.dir)
        f=open("%s/OUTPUT.merged/list"%self.dir,"w")
        print >>f,"\n".join(a)
        f.close()
        os.chdir("%s/OUTPUT.merged"%(self.dir))
        #com="%s/bin/merge.csh %s true"%(cWB_conf.run_dir,cWB_conf.run_dir)
        com="echo list cWB_online.root true | root -l -b -q %s/tools/online/RUN_cWB/bin/mergePROD.C"%(os.environ['HOME_WAT'])
        if (os.environ['SITE_CLUSTER']=="CASCINA"):
          com="echo list cWB_online.root true | root -l -b -q -n %s %s/tools/online/RUN_cWB/bin/mergePROD.C"%(os.envinron['CWB_ROOTLOGON_FILE'],os.environ['HOME_WAT'])
        a=commands.getstatusoutput(com)
        #print a
        count=3
        while(a[1].find("Error")!=-1):
            a=commands.getstatusoutput(com)
            print a
            count-=1
            if(count<=0):
                print "Error: giving up on com"
                break
            os.system("sleep %d"%cWB_conf.sleep)
            
    def extract_triggers(self):
        root_files=glob.glob("%s/OUTPUT/*.root"%self.dir)
        print "%s/OUTPUT/*.root root_files: %i"%(self.dir,len(root_files))
        os.chdir("%s/OUTPUT.merged"%(self.dir))
        live_file="live.root"
        print "live"
        if (not os.path.exists(live_file)):
             #com="%s/bin/getLive.csh %s %s %s %s %s"%(cWB_conf.run_dir,cWB_conf.run_dir,root_files[0],live_file,self.start,self.end)
             com="echo %s %s %s %s | root -l -b -q %s/tools/online/bin/ExtractLive.C"%(root_files[0],live_file,self.start,self.end,os.environ['HOME_WAT'])
             if (os.environ['SITE_CLUSTER']=="CASCINA"):
               com="echo %s %s %s %s | root -l -b -q -n %s %s/tools/online/bin/ExtractLive.C"%(root_files[0],live_file,self.start,self.end,os.environ['CWB_ROOTLOGON_FILE'],os.environ['HOME_WAT'])
             commands.getstatusoutput(com)
        nop_triggers=filter(lambda x: x.topublish==False, self.triggers_all)
        print "nop: %i"%len(nop_triggers)
        for trigger in nop_triggers:
           file="trigger_%s.root"%trigger.time[0]
           if (os.path.exists(file)):
              #print "delete file %s"%file
              commands.getstatusoutput("rm %s"%file)
        p_triggers=filter(lambda x: x.topublish==True, self.triggers_all)
        print "p: %i"%len(p_triggers)
        for trigger in p_triggers:
          for root_file in root_files: 
              #com="%s/bin/select.csh %s %s %s %s %s"%(cWB_conf.run_dir,cWB_conf.run_dir,root_file,trigger.time[0],trigger.start[0],trigger.stop[0])
              com="echo %s %s %s %s | root -l -b -q %s/tools/online/bin/SelectTriggers.C"%(root_file,trigger.time[0],trigger.start[0],trigger.stop[0],os.environ['HOME_WAT'])
              if (os.environ['SITE_CLUSTER']=="CASCINA"):
                com="echo %s %s %s %s | root -l -b -q -n %s %s/tools/online/bin/SelectTriggers.C"%(root_file,trigger.time[0],trigger.start[0],trigger.stop[0],os.environ['CWB_ROOTLOGON_FILE'],os.environ['HOME_WAT'])
              #print com
              res=commands.getstatusoutput(com)

    def dump_triggers(self,segs):
        commands.getstatusoutput("mkdir -p %s/OUTPUT.merged"%(self.dir))
        os.chdir("%s/OUTPUT.merged"%(self.dir))
        commands.getstatusoutput("mkdir -p TRIGGERS")
        tfiles=glob.glob("%s/OUTPUT/*.txt"%(self.dir))
        #a=commands.getstatusoutput("cat %s > TRIGGERS/triggers.txt"%(" ".join(tfiles)))
        print "1";sys.stdout.flush()
        for tfile in tfiles:
            print "tfile=%s"%tfile
            self.triggers_all+=ascii_dump_2_triggers(tfile,segs,self)
            print "-->len(triggers_all)=%d"%len(self.triggers_all);sys.stdout.flush()
        print "2";sys.stdout.flush()            
        self.triggers_all.sort(cmp=t_cmp_time)
        print "3";sys.stdout.flush()
        self.selected_triggers=filter(lambda x: x.selected>=3, self.triggers_all)
        #if(cWB_conf.apply_veto and len(self.selected_triggers)>0):
        #    try:
        #        self.dq()
        #    except Exception,e:
        #        print "Ignoring errors in job.dq() for now"
        #        print e
        #        print "="*30
        #for t in self.triggers_all:
            #if(t.selected>=3):
            #    t.qscan()
            #t.dump_error_region("%s/OUTPUT.merged/TRIGGERS/error_region_%s"%(self.dir,t.time[0]),"%s/OUTPUT.merged/TRIGGERS/error_region_%s"%(self.link,t.time[0]))
        print "4";sys.stdout.flush()
        self.triggers_dump_time=get_time_nooff()
        self.delay_start_dump=self.triggers_dump_time-self.start #time between the start time of the processed segment and the moment triggers are dumped
        self.delay_start_launch=self.launch_time-self.start #time between the start time of the processed segment and the moment cWB job was launched: h(t) generation + h(t) transfer + h(t) discovery
        self.delay_launch_completion=self.completion_time - self.launch_time #time it took to run cWB trigger production part

    def select_triggers(self):
        pass
    def send_to_external_collaborations(self):
        print "In send_to_external_collaborations"
        self.triggers_send_time=get_time_nooff()
        #print "Before significance"
        #for trigger in filter(lambda x: x.selected>=3 and x.topublish==True, self.triggers_all):
        for trigger in self.triggers_all:
          #trigger.significance()
          #print trigger.path
          #f=open(trigger.path,"a")
          #print >>f,trigger.path
          #f.close()
          if (trigger.selected>=3 and trigger.topublish==True):
            if(trigger.send_time>0):
                print "trigger %s has been already sent"%(repr(trigger.time[0]))
                continue
            self.trigger_to_external_collaborations(trigger)
    def check_triggers(self,jobs_list):
      for trigger in self.triggers_all:
        #print trigger.path
        coinc_jobs=filter(lambda x: float(trigger.time[0])>int(x.back) and float(trigger.time[0])<int(x.end) and x!=self, jobs_list)
        #print repr(coinc_jobs)
        tseg=segment(float(trigger.start[0]),float(trigger.stop[0]))
        tband=segment(float(trigger.left[0]),float(trigger.right[0]))
        final_trig=trigger
        for jj in coinc_jobs:
           for trigger2 in jj.triggers_all:
               #if (tseg.intersects(segment(float(trigger2.start[0]),float(trigger2.stop[0]))) and trigger2.topublish==True):
               if (tseg.intersects(segment(float(trigger2.start[0]),float(trigger2.stop[0]))) and tband.intersects(segment(float(trigger2.left[0]),float(trigger2.right[0]))) and trigger2.topublish==True):
                  #if ((float(trigger2.netcc[cWB_conf.id_cc])>cWB_conf.th_cc and float(trigger2.rho[cWB_conf.id_rho])>float(trigger.rho[cWB_conf.id_rho])) or (float(trigger2.netcc[cWB_conf.id_cc])<cWB_conf.th_cc and float(trigger.netcc[cWB_conf.id_cc])<cWB_conf.th_cc and float(trigger2.rho[cWB_conf.id_rho])>float(trigger.rho[cWB_conf.id_rho]))):
                  if ((trigger2.selected>0 and float(trigger2.rho[cWB_conf.id_rho])>float(trigger.rho[cWB_conf.id_rho])) or (trigger2.selected==0 and trigger.selected==0 and float(trigger2.rho[cWB_conf.id_rho])>float(trigger.rho[cWB_conf.id_rho]))):
                      final_trig=trigger2
                      trigger.topublish=False
                  else:
                      trigger2.topublish=False
                      #trigger.graceid=trigger2.graceid
                      trigger.first_graceid="%s%s"%(trigger2.graceid,trigger2.first_graceid)
                      trigger.gracedb_link="%s"%(trigger2.gracedb_link)
                      trigger.comments=trigger2.comments
                      jj.extract_triggers()
                      jj.web_page()
                      jj.dump_job()
        self.dump_job()

    def trigger_to_external_collaborations(self,trigger):
        #print "Before record_passed"
        #trigger.record_passed()
        #print "Before record_dq"
        #trigger.record_dq()
        #trigger.path=trigger.link.replace(cWB_conf.web_link,cWB_conf.online_dir)
        #ced_link=trigger.compute_ced()
        #print >>f,ced_link
        #print >>f,ced_link+"skyprobcc.fits"
        print "Before send_time"
        trigger.send_time=self.triggers_send_time
        if (cWB_conf.sendtogracedb==True and float(trigger.sgnf[2].split(" ")[1])>0.):
         print "Before gracedb"
         try:
            client = GraceDb(cWB_conf.gracedb_client)
         except:
            client = GraceDb()
         try:
            if (trigger.first_graceid==""):
              r = client.createEvent(cWB_conf.gracedb_group, cWB_conf.gracedb_analysis, trigger.path, search=cWB_conf.gracedb_search)
              r = r.json()
              graceid = r["graceid"]
              trigger.graceid=graceid
              trigger.gracedb_link="https://gracedb.ligo.org/events/view/%s"%(trigger.graceid)
              trigger.comments="Sent to GraceDB at %s, delay=%.3f"%(repr(self.triggers_send_time),float(self.triggers_send_time) - float(trigger.time[0]))
            else:
              trigger.comments+="<br>Uploaded at %s (delay=%.3f)"%(repr(self.triggers_send_time),float(self.triggers_send_time) - float(trigger.time[0]))
         except HTTPError, e:
            print "Error on sending to gracedb"
            print e
            #f = open("tmp.html", "w")
            #f.write(str(e))
            #f.close()
        else:
         trigger.graceid="T00000"
         trigger.gracedb_link="https://gracedb.ligo.org/events/view/%s"%(trigger.graceid)
     
        rate=float(trigger.sgnf[2].split(" ")[1])
        num=int(trigger.sgnf[2].split(" ")[0])
        if(num>0):
          if (rate*3600*24*365>1):
            ifar1_f="%.1f per year"%(rate*3600*24*365)
          else:
            ifar1_f="1 per %.1f year"%(1./(rate*3600*24*365))
        else:
          if (rate*3600*24*365>1):
            ifar1_f="<%.1f per year"%(rate*3600*24*365)
          else:
            ifar1_f=">1 per %.1f year"%(1./(rate*3600*24*365))
        if (rate<0):
         ifar1_f="WARNING NEG IFAR %s"%ifar1_f
        msg=\
"""Subject:  %s %s: Event GPS %s %s

%s
network = %s, injection: %s, search: %s
rho = %s, netcc = %s, Far: %s Hz, %s
frequency band: [%.2f,%.2f] Hz, %s

CED: %s

Gracedb entry: %s
-------------------------------------------------------------------
File on %s: %s
File on Web: %s/OUTPUT.merged/TRIGGERS/trigger_%s.txt
-------------------------------------------------------------------
Main online page: %s
        """%("".join(cWB_conf.ifos), cWB_conf.gracedb_search, trigger.time[0], ifar1_f,
             cWB_conf.title,
             "".join(cWB_conf.ifos), trigger.inj_found, cWB_conf.gracedb_search,
             trigger.rho[cWB_conf.id_rho], trigger.netcc[cWB_conf.id_cc], trigger.sgnf[2].split(" ")[1], ifar1_f, float(trigger.low[0]) ,float(trigger.high[0]), trigger.cutclass,
             trigger.compute_ced(), trigger.gracedb_link, os.environ['SITE_CLUSTER'],
             trigger.path, trigger.job.link,trigger.time[0],cWB_conf.web_link)

        if (cWB_conf.sendmail==True and trigger.first_graceid=="" and (float(trigger.rho[cWB_conf.id_rho])>th_rho_mail or rate<0)):
         try:
            com="echo \"%s\" | /usr/sbin/sendmail -F %sonline-%s-%s %s"%(msg,cWB_conf.gracedb_analysis,cWB_conf.gracedb_search,"".join(cWB_conf.ifos),",".join(cWB_conf.emails))
            a=commands.getstatusoutput(com)
            print "="*30
            print msg
            print "="*30
            print a
            print "="*30
         except HTTPError, e:
            print "Error on sending mail"
            print e
         #try:
         #  rate=float(trigger.sgnf[1].split(" ")[1])
         #  if (rate<cWB_conf.phone_alert):
         #   com=""" cat %s | mailx -s "Event passed criteria for follow up"  %s """%(cWB_conf.phone_par,cWB_conf.phone_mail)
         #   print com
         #   commands.getstatusoutput(com)
         #except ValueError:
         #   print "could not convert data to an integer"
         #except:
         #   print "Error on sending phone alert"
        if (cWB_conf.sendtogracedb==True and trigger.first_graceid==""):
          if (trigger.inj_found==True):
            try:
              r = client.writeLabel(trigger.graceid,"INJ")
            except:
              print "No INJ info sent"
          #client = GraceDb()
          ced_dir="%s/OUTPUT/%s"%(trigger.job.dir,trigger.compute_ced_dir())
          fit_file="%s/skyprobcc.fits"%(ced_dir)
          #fit_file_cwb="%s/skyprobcc_cWB.fits"%(ced_dir)
          fit_file_cwb="%s/cWB.fits"%(ced_dir)
          print "fit_file = %s"%(fit_file)
          print "fit_file_cwb = %s"%(fit_file_cwb)
          #com ="cp %s %s"%(fit_file,fit_file_cwb)
          com ="cp %s %s;gzip %s"%(fit_file,fit_file_cwb,fit_file_cwb)
          fit_file_cwb="%s/cWB.fits.gz"%(ced_dir)
          #print com
          commands.getstatusoutput(com)
          try:
            r = client.writeLog(trigger.graceid, "cWB skymap fit", filename=fit_file_cwb, tag_name=["sky_loc", "lvem","public"])
            client.writeLabel(trigger.graceid,"SKYMAP_READY")
          except:
            print "no skymap sent"
          png_file_cwb="%s/cWB.png"%(ced_dir)
          com="ligo-skymap-plot %s -o %s --contour 50 90 --annotate"%(fit_file,png_file_cwb)
          print com
          commands.getstatusoutput(com)
          try:
            r = client.writeLog(trigger.graceid, "cWB skymap png", filename=png_file_cwb, tag_name="sky_loc")
          except:
            print "no ligo skymap sent"
          try:
            r = client.writeLog(trigger.graceid, "<a href=\"%s\">cWB CED</a>"%(trigger.compute_ced()), tag_name="analyst_comments")
          except:
            print "no CED sent"
          html_code="""
<table>
<tr><th colspan=2>cWB parameter estimation</th></tr>
<tr><td>Frequency [Hz]</td><td align=right>%.3f</td></tr>
<tr><td>Bandwidth [Hz]</td><td align=right>%.3f</td></tr>
<tr><td>Duration [ms]</td><td align=right>%.2f</td></tr>
<tr><td>Hrss</td><td align=right>%.3e</td></tr>
</table>
"""%(float(trigger.frequency[0]),float(trigger.bandwidth[0]),float(trigger.duration[0])*1000.,float(trigger.strain[0]))
          print html_code
          try:
            r = client.writeLog(trigger.graceid, html_code, tag_name="pe")
          except:
            print "no pe sent"
       
    def format_for_sending(self):
        pass
    def send(self):
        pass
    def dump_job(self):
        #print "dump job()->%s"%self.dir; sys.stdout.flush()
        f=open("%s/job.pickle"%self.dir,"wb")
        pickle.dump(self,f,protocol=2)
        f.close()
    def web_page(self):
        title="%s, cWB, all-sky, unmodeled, <br>job %d - %d (%s - %s),<br>generated at %s"%(cWB_conf.label, self.start, self.end, tconvert(self.start),tconvert(self.end),tconvert(get_time_nooff()))
        statistics_table=\
        """
        <table border=1>
        <tr>
        <th>Launch time</th><td>%d<td>%s</td>
        </tr>
        <tr>
        <th>Launch time - start of the interval</th><td>%d</td><td></td>
        </tr>
        <tr>
        <th>Trigger generation completion time</th><td>%d</td><td>%s</td>
        </tr>
        <tr>
        <th>Completion time - launch time</th><td>%d</td><td></td>
        </tr>
        <tr>
        <th>Triggers send time</th><td>%d</td><td>%s</td>
        </tr>
        <tr>
        <th>Total: send time - start of the interval</th><td>%d</d>
        </tr>
        <tr>
        <th>Log (zipped)</th><td><a href=\"log.zip\">log</a></td>
        </tr>
        <tr>
        <th>Job dir</th><td><a href=\".\">dir</a></td>
        </tr>
        </table>

        """ % (self.launch_time,tconvert(self.launch_time),self.delay_start_launch,self.completion_time,tconvert(self.completion_time),\
               self.delay_launch_completion, self.triggers_send_time,tconvert(self.triggers_send_time), self.delay_start_send)

        p_triggers=filter(lambda x: x.topublish==True, self.triggers_all)
        if(len(p_triggers)>0):
            trigger_table_rows="\n".join(map(lambda x: x[1].to_html_table_row(x[0]), zip(range(1,len(p_triggers)+1),p_triggers)))
            trigger_table_header=p_triggers[0].table_header(cWB_conf.ifos,True)
        else:
            trigger_table_rows=""
            trigger_table_header="<tr><th>No triggers found</th></tr>"
        trigger_table="""
        <table border=1>
        %s
        %s
        </table>
        """%(trigger_table_header, trigger_table_rows)

        nop_triggers=filter(lambda x: x.topublish==False, self.triggers_all)
        if(len(nop_triggers)>0):
            trigger_table_rows="\n".join(map(lambda x: x[1].to_html_table_row(x[0]), zip(range(1,len(nop_triggers)+1),nop_triggers)))
            trigger_table_header=nop_triggers[0].table_header(cWB_conf.ifos,True)
            trigger_table+="""
        <h3 align="center"> Triggers in this job but not reported because more significant in other job</h3>
        <table border=1>
        %s
        %s
        </table>
            """%(trigger_table_header, trigger_table_rows)

        selection_explaination=\
        """
        <table>
        <tr><th>0</th><td>no postproduction cuts;</td></tr>
        <tr><th>1</th><td>only coherent triggers selected but no cut on rho or DQ;</td></tr>
        <tr><th>2</th><td>rho cut for internal LSC follow up with qscans, ced, etc.;</td></tr>
        <tr><th>3</th><td>rho cut for external collaborations;</td><tr>
        <tr><th>4</th><td>rho cut and DQ and vetoes for external collaborations;</td></tr>
        <tr><th>5</th><td>rho cut and DQ and vetoes that would most likely be chosen for offline papaper.</td></tr>
        </table>
        """
        #htmltop=self.top.replace("\n","<br>\n")
        msg=\
        """
        <html>
        <title>%s</title>
        <body>
        <h1 align=center>%s</h1>
        %s
        <h2 align=center>Triggers</h2>
        %s
        <h3>Explaining selection values from the trigger table:</h3>
        %s
        </body>
        </html>

        """ % (title,title,statistics_table, trigger_table, selection_explaination)
        fn="%s/job.html"%(self.dir)
        #fn="%s/job.html"%(self.link.replace(cWB_conf.web_link,cWB_conf.online_dir))
        f=open(fn,"w")
        print >>f, msg
        f.close()

    def join_triggers(self,trig_file):
        print "In join_triggers"
        #print "root_file"
        #print root_file
        #print "trig_file"
        #print trig_file
        num=deque([])
        triggers2=deque([])
        for c in trig_file:
            #print c
            trigger_dir="/".join(c.split("/")[:-3])
            f=open("%s/job.pickle"%(trigger_dir),"rb")
            c_job=pickle.load(f)
            f.close()
            segs=segmentlist([segment(int(c.split("/")[8].split("-")[0]),int(c.split("/")[8].split("-")[1]))])
            lines=filter(lambda y: len(y)>0, map(lambda x: x.strip(), open(c).readlines()))
            trigger=Trigger(lines,segs)
            trigger.job=c_job
            trigger.dir="%s/OUTPUT.merged/TRIGGERS"%(c_job.dir)
            trigger.path="%s/OUTPUT.merged/TRIGGERS/trigger_%s.txt"%(c_job.dir, trigger.time[0])
            trigger.link="%s/OUTPUT.merged/TRIGGERS/trigger_%s.txt"%(c_job.link, trigger.time[0])
            triggers2.append(trigger)
        for q in range(0,len(triggers2)):
            trigger=triggers2[q]
            t_num=q
            tseg=segment(float(trigger.start[0]),float(trigger.stop[0]))
            for l in range(0,len(triggers2)):
                trigger2=triggers2[l]
                d=trig_file[l]
                segs=segmentlist([segment(int(d.split("/")[8].split("-")[0]),int(d.split("/")[8].split("-")[1]))])
                if (tseg.intersects(segment(float(trigger2.start[0]),float(trigger2.stop[0])))):
                  if ((float(trigger2.netcc[cWB_conf.id_cc])>cWB_conf.th_cc and float(trigger2.rho[cWB_conf.id_rho])>float(trigger.rho[cWB_conf.id_rho])) or (float(trigger2.netcc[cWB_conf.id_cc])<cWB_conf.th_cc and float(trigger.netcc[cWB_conf.id_cc])<cWB_conf.th_cc and float(trigger2.rho[cWB_conf.id_rho])>float(trigger.rho[cWB_conf.id_rho]))):
                    #print "%f %f %f"%(float(trigger2.netcc[cWB_conf.id_cc]),float(trigger2.rho[cWB_conf.id_rho]),float(trigger.rho[cWB_conf.id_rho]))
                    trigger=trigger2
                    t_num=l
            if (len(num)==0):
               self.triggers_all.append(trigger)
               num.append(t_num)
               #print "t_num: %i"%t_num
            else:
               #print "t_num: %i"%t_num
               copy=1
               for i in num:
                   if(t_num==i):
                    copy=0
               if(copy==1):
                  self.triggers_all.append(trigger)
                  num.append(t_num)
        #print "job.dir: %s"%self.dir
        olddir=os.getcwd()
        os.chdir(self.dir)
        #input_file="input/root_files.txt"
        #f=open(input_file,"wb") 
        for trigger in self.triggers_all:
            trigger_dir="/".join(trigger.dir.split("/")[:-2])
            com=("head -n -3 %s > %s/OUTPUT/%s"%(trigger.path,self.dir,trigger.path.split("/")[len(trigger.path.split("/"))-1]))
            #com=("cp %s %s/OUTPUT/."%(trigger.path,self.dir))
            #print com
            commands.getstatusoutput(com)
            trigger.dir=trigger.dir.replace(trigger_dir,self.dir)
            trigger.path=trigger.path.replace(trigger_dir,self.dir)
            trigger.link=trigger.link.replace("/".join(trigger_dir.split("/")[7:9]),"/".join(self.dir.split("/")[6:8]))
            trigger.job.link=trigger.job.link.replace("/".join(trigger_dir.split("/")[7:9]),"/".join(self.dir.split("/")[6:8]))
            #print trigger.dir
            #print trigger.path
            #print trigger.link
            root_files=glob.glob("%s/OUTPUT/*.root"%trigger_dir)
            #print root_files
            for root_file in root_files:
                #print >>f, "%s %s %s %s"%(root_file,trigger.time[0],trigger.start[0],trigger.stop[0]) 
                #com="%s/bin/select.csh %s %s %s %s %s"%(cWB_conf.run_dir,cWB_conf.run_dir,root_file,trigger.time[0],trigger.start[0],trigger.stop[0]) 
                com="echo %s %s %s %s | root -l -b -q %s/tools/online/bin/SelectTriggers.C"%(root_file,trigger.time[0],trigger.start[0],trigger.stop[0],os.environ['HOME_WAT'])
                if (os.environ['SITE_CLUSTER']=="CASCINA"):
                  com="echo %s %s %s %s | root -l -b -q -n %s %s/tools/online/bin/SelectTriggers.C"%(root_file,trigger.time[0],trigger.start[0],trigger.stop[0],os.environ['CWB_ROOTLOGON_FILE'],os.environ['HOME_WAT'])
                #print com
                res=commands.getstatusoutput(com)
                #print res
        os.chdir(olddir)
    def zip_files(self):
       log_file="%s/log"%self.dir
       com="zip %s.zip %s;rm %s"%(log_file,log_file,log_file) 
       commands.getstatusoutput(com)
       #input_dir="%s/input"%self.dir
       #com="tar -czf %s.tgz %s; rm -rf %s"%(input_dir,input_dir,input_dir)
       #commands.getstatusoutput(com)

class FRAME:
    def __init__(self,path):
        self.pfn = path.strip()
        self.lfn = os.path.basename(self.pfn)
        a = self.lfn.split("-")
        self.ifos = a[0]
        self.type = a[1]
        self.start = int(a[-2])
        self.duration = int(a[-1].split(".")[0])
        self.end = self.start+self.duration
    def __cmp__(self,other):
        return cmp(self.start,other.start)
    def __repr__(self):
        msg="""
        ======= frame ======
        pfn = %s
        lfn = %s
        ifos = %s
        type = %s
        start = %d
        end = %d
        duration = %d
        ====================        
        """ % (self.pfn, self.lfn, self.ifos, self.type, self.start, self.end, self.duration)
        return msg
    
class cWB:
    def __init__(self, run_end=-1):
        self.start_loop = get_time()
        self.jobs_completed=deque([])
        self.jobs_running=deque([])
        self.jobs_failed=deque([])
        self.jobs_queued=deque([])
        self.jobs_published=deque([])
        self.frames={}
        self.frames_segs={}
        self.last_read={}
        self.run_end=run_end

    def __repr__(self):
        pass
    def schedule(self):
        frame_coverage=self.frames_segs[cWB_conf.ifos[0]]
        for ifo in cWB_conf.ifos[1:]:
            frame_coverage &= self.frames_segs[ifo]
            frame_coverage.coalesce()
        jobs_completed_segs=segmentlist(map(lambda x: segment(x.start,x.end), \
                                            list(self.jobs_completed) + \
                                            list(self.jobs_published) + \
                                            list(self.jobs_failed)
                                            )
                                        )
        jobs_completed_segs.coalesce()
        jobs_running_segs=segmentlist(map(lambda x: segment(x.start,x.end), self.jobs_running)); jobs_running_segs.coalesce()
        new_segs= frame_coverage & self.segs_global; new_segs.coalesce()
        new_segs-=jobs_completed_segs; new_segs.coalesce()
        new_segs-=jobs_running_segs; new_segs.coalesce()
        current_time=get_time()
        look_back=segmentlist([segment(current_time-cWB_conf.look_back,current_time)])
        new_segs&=look_back
        new_segs.coalesce()

        if(cWB_conf.debug==1):
            print "new_segs=%s"%(repr(new_segs))
            print "abs(new_segs)=%s"%(repr(map(abs,new_segs)))
            print "frame_coverage=%s"%(repr(frame_coverage))
            print "look_back=%s"%(repr(look_back))

            #print "GLOBAL SEGS: %s"%repr(self.segs)
            print "GLOBAL GLOBAL: %s"%repr(self.segs_global)
        self.scheduling_policy_simplest(new_segs,frame_coverage)
        
    def scheduling_policy_simplest(self,segs,frame_coverage):
        if(len(segs)==0 or max(map(lambda x:abs(x),segs))<cWB_conf.moving_step+2*cWB_conf.job_offset):
            return
        if (max(map(lambda x:abs(x),self.segs_global))<cWB_conf.seg_duration+2*cWB_conf.job_offset):
            return
        for seg in segs:
            if(abs(seg)<cWB_conf.moving_step+2*cWB_conf.job_offset):
                continue
            a=segmentlist([segment(seg[0]-cWB_conf.job_offset,seg[0])])
            a = a & self.segs_global
            if(len(a)==1 and a[0][1]==seg[0]):
                jstart=a[0][0]+cWB_conf.job_offset
            else:
                jstart=seg[0]+cWB_conf.job_offset
            jend=jstart+cWB_conf.moving_step
            back=cWB_conf.seg_duration-cWB_conf.moving_step
            temp_back=jstart-back
            #print "GLOBAL: %s"%repr(self.segs_global)
            while(jend + cWB_conf.job_offset <= seg[1]):
                whiteseg=segmentlist([segment(temp_back-cWB_conf.job_offset,jend+cWB_conf.job_offset)])
                #print repr(whiteseg)
                whiteseg = whiteseg & self.segs_global
                #print repr(whiteseg)
                #print jstart
                #print jend
                #print temp_back
                if (abs(whiteseg)<cWB_conf.seg_duration+2*cWB_conf.job_offset):
                    #print "too short"
                    jend=jend+cWB_conf.moving_step
                    temp_back += cWB_conf.moving_step
                    #print jstart
                    #print jend
                    #print temp_back
                else:
                    #print "ok"
                    #print jstart
                    #print jend
                    #print temp_back
                    job=JOB(jstart,jend)
                    job.back=temp_back
                    #if (jstart-back < seg[0]):
                    #    job.back=seg[0]+cWB_conf.job_offset
                    job.queue_time=get_time_nooff()
                    self.jobs_queued.append(job)
                    jstart=jend
                    jend=jstart+cWB_conf.moving_step
                    temp_back=jstart-back

    def launch(self):
        i=len(self.jobs_running)
#        while(len(self.jobs_queued)>0 and i<=cWB_conf.max_jobs):        
        while(len(self.jobs_queued)>0):
            if(cWB_conf.debug==1):
                print "In launch: number of running jobs=%d=%d"%(i,len(self.jobs_running))
            job=self.jobs_queued.popleft()
            if(job.status==4):
                self.jobs_failed.append(job)
                continue
            job.status=1
            if(job.launch_time<0):
                job.launch_time = get_time_nooff()
                job.launch()
                self.jobs_running.append(job)
            if(cWB_conf.debug==1):
                print "job start=%d end=%d launched=%d status=\"%s\""%(job.start,job.end,job.launch_time,job.interpret_status())
            job.dump_job()
            i+=1
    def configure(self):
        for job in self.jobs_queued:
            job.configure(self.frames)
    def check(self):
        c=get_time_nooff()
        print "In check()"
        print "c=%s"%c
        print "n=%d"%(len(self.jobs_running))
        for job in self.jobs_running:
            #print job.dir
            a=glob.glob(job.dir+"/OUTPUT/finished")
            #print a
            if(len(a)>0):
                job.status=2
                job.completion_time=c
            else:
                if(c-job.launch_time>cWB_conf.job_timeout):
                    job.status=3
                    job.completion_time=c
            try:
                job.dump_job()
            except Exception,e:
                print "job.dump failed"
                print e
        aa=deque([])
        while(len(self.jobs_running)>0):
            job=self.jobs_running.pop()
            if(job.status==1):
                aa.append(job)
            elif(job.status==2):
                self.jobs_completed.append(job)
            elif(job.status==3):
                self.jobs_failed.append(job)
        self.jobs_running=aa

        if(cWB_conf.debug==1):
            print "In check:"
            print "Jobs queued: %d"%(len(self.jobs_queued))
            print "*"*10
            print ("-"*5+"\n").join(map(repr,self.jobs_queued))
            print "*"*10            
            print "Jobs running: %d"%(len(self.jobs_running))
            print "*"*10
            print ("-"*5+"\n").join(map(repr,self.jobs_running))
            print "*"*10                                    
            print "Jobs completed: %d"%(len(self.jobs_completed))
            print "*"*10
            print ("-"*5+"\n").join(map(repr,self.jobs_completed))
            print "*"*10                        
            print "Jobs failed: %d"%(len(self.jobs_failed))
            print "*"*10
            print ("-"*5+"\n").join(map(repr,self.jobs_failed))
            print "*"*10
            print "Jobs published: %d"%(len(self.jobs_published))
            print "*"*10
            print ("-"*5+"\n").join(map(repr,self.jobs_published))
            print "*"*10                                    
        self.jobs_failed=filter(lambda x: segment(x.start-cWB_conf.job_offset, x.end+cWB_conf.job_offset) in \
                                segment(self.cycle_start-cWB_conf.look_back, self.cycle_start), self.jobs_failed)
    def publish(self):
#        end=get_time()
#        start=end-3600
        try:
            job=self.jobs_completed.popleft()
        except:
            print "Why am I here: len(jobs_completed)=%d"%(len(self.jobs_completed))
            return
        while(True):
#            self.discover_segments(start,end)
            job.publish(self.segs_global,self.jobs_published)
            self.jobs_published.append(job)
            job.dump_job()
            job.zip_files()
            try:
                job=self.jobs_completed.popleft()
            except:
                print "Why am I here: len(jobs_completed)=%d"%(len(self.jobs_completed))                
                break
        self.jobs_published=filter(lambda x: segment(x.start-cWB_conf.job_offset, x.end+cWB_conf.job_offset) in \
                            segment(self.cycle_start-5*cWB_conf.look_back, self.cycle_start), self.jobs_published)
    def run(self):
        for ifo in cWB_conf.ifos:
          self.last_read[ifo]=0#get_time()-5*cWB_conf.look_back
        while(True):
            try:
                self.cycle_start=get_time()
                if(self.run_end>0 and self.cycle_start > self.run_end):
                    print "Run is finished, run_end=%s, current_time=%s"%(self.run_end, self.cycle_start)
                    return
                self.cycle()
            except Exception,e:
                print "The cycle broke due to the exception %s"%(repr(e))
                print "Restarting ..."
                sys.exit(1)
            if(cWB_conf.debug==1):
                print "Sleepting for %d seconds"%cWB_conf.sleep
            sys.stdout.flush()
            sys.stderr.flush()
            os.system("sleep %d"%cWB_conf.sleep)        
    def copy(self):
      fn_combined_global="%s/%s/seg_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir)
      try:
        f=open(fn_combined_global)
        self.segs_global=fromsegwizard(f)
        f.close()
      except Exception,e:
        self.segs_global=segmentlist([])
      if (abs(self.segs_global)>0):
        started_dirs=glob.glob("%s/%s/??????"%(cWB_conf.run_dir.replace(cWB_conf.online_dir,cWB_conf.production_dir),cWB_conf.jobs_dir))
        for d in started_dirs:
          dirs_orig=glob.glob("%s/*"%d) 
          dirs_copied=glob.glob("%s/*"%d.replace(cWB_conf.production_dir,cWB_conf.online_dir))
          print "%s %i %i"%(d,len(dirs_orig),len(dirs_copied))
          if (len(dirs_orig)>len(dirs_copied)):
            for o in dirs_orig:
             c=o.replace(cWB_conf.production_dir,cWB_conf.online_dir)
             if (os.path.exists("%s/OUTPUT/finished"%o) and not os.path.exists("%s/OUTPUT/finished"%c)):
              try:
               f=open("%s/job.pickle"%o,"rb")
               j=pickle.load(f)
               f.close()
               if (j.status==2):
                start=j.start
                stop=j.end
                job=JOB(start,stop)
                job.back=j.back
                job.launch_time=j.launch_time
                job.completion_time=j.completion_time
                job.status=j.status # 0 - not launched yet, 1 running, 2 finished successfully, 3 - failed, 4 - missing frames
        #self.frames={}
                job.delay_start_ready=j.delay_start_ready
                job.delay_start_launch=j.delay_start_launch
                job.dump_job()
                com="cp -r %s/OUTPUT %s/."%(o,c)
                print com
                commands.getstatusoutput(com)
                com="cp -r %s/input/* %s/input/."%(o,c)
                print com
                commands.getstatusoutput(com)
                self.jobs_completed.append(job)
              except Exception,e:
               print e
      #except Exception,e:
      #  print e
      #  self.segs_global=segmentlist([])
      print repr(self.segs_global)
      print len(self.jobs_completed)
    def cycle(self):
      if(cWB_conf.debug==1):
            print "="*20 + "Cycle start time=%d"%(self.cycle_start) + "="*20
      if (pp_redo==False):
        print "In cycle before discover_segments"; sys.stdout.flush()
        self.discover_segments(self.cycle_start - cWB_conf.look_back, self.cycle_start)
        print "In cycle before schedule"; sys.stdout.flush()        
        self.schedule()
        print "In cycle before configure"; sys.stdout.flush()                
        self.configure()
        print "In cycle before launch"; sys.stdout.flush()                
        self.launch()
        print "In cycle before check"; sys.stdout.flush()                
        self.check()
      else:
        print "In cycle before copy"; sys.stdout.flush()                
        self.copy()
      print "In cycle before publish"; sys.stdout.flush()                
      self.publish()
      print "In cycle at the end"; sys.stdout.flush()                
        #f=open("%s/PICKLES/run_%d.pickle"%(cWB_conf.online_dir, self.cycle_start),"wb")
        #pickle.dump(self,f,protocol=2)
        #f.close()

    def discover_frames(self,ifo):
        ifo_to_index={}
        for i in range(len(cWB_conf.ifos)):
           ifo_to_index[cWB_conf.ifos[i]]=i
        print "In discover_frames %s"%(ifo)
        #z=glob.glob("%s/%s/*.gwf"%(cWB_conf.frames_dir[cWB_conf.ifo_to_index[ifo]], ifo))
        z=glob.glob("%s/*.gwf"%(cWB_conf.frames_dir[ifo_to_index[ifo]]))
	z.sort()
        self.frames[ifo]=map(lambda x: FRAME(x), z)
        print "Number of unfiltered frames %d"%(len(self.frames[ifo]))
        #self.frames[ifo]=filter(lambda x: segment(x.start,x.end) in segment(self.cycle_start - cWB_conf.look_back, self.cycle_start),self.frames[ifo])
        self.frames_segs[ifo]=segmentlist(map(lambda x: segment(x.start,x.end),self.frames[ifo]))
        self.frames_segs[ifo].coalesce()
        self.frames[ifo]=filter(lambda x: float(x.start)>= float(self.last_read[ifo]),self.frames[ifo])
        print "Number of filtered frames %d"%(len(self.frames[ifo]))

    def discover_segments(self,start,end):
        s={}
        for ifo in cWB_conf.ifos:
            self.discover_frames(ifo)
            s[ifo]=segmentlist([])
            #print ifo
            #print len(self.frames[ifo])
            inj_sg={}
            for inj in range(0,len(cWB_conf.inj_name)):
                inj_sg[inj]=segmentlist([])
            cat2_sg={}
            try:
                cat2veto=len(cWB_conf.cat2_channel[ifo])
            except:
                cat2veto=0
            for veto in range(0,cat2veto):
                cat2_sg[veto]=segmentlist([])
            veto_sg={}
            try:
                nveto=len(cWB_conf.veto_channel[ifo])
            except:
                nveto=0
            for veto in range(0,nveto):
                veto_sg[veto]=segmentlist([])
            #print "NVETO: %i"%(nveto)
            for frame in self.frames[ifo]:
                    #if (len(cWB_conf.DQ_channel[ifo])>1):
                    try:
                      temps=segmentlist([])
                      for idq in range(len(cWB_conf.DQ_channel[ifo])):
                        dq=len(cWB_conf.DQ_channel[ifo])-idq-1
                        bitmask_list=self.getbitmask_fromfile(frame,ifo,cWB_conf.DQ_channel[ifo][dq],cWB_conf.DQ_channel_rate[ifo][dq])
                        print "In frame_analysis_segments: %s %s"%(frame.pfn,cWB_conf.DQ_channel[ifo][dq])
                        if (idq==0):
                           temps=self.frame_analysis_segments(frame,ifo,cWB_conf.bitmask[ifo][dq],cWB_conf.DQ_channel_rate[ifo][dq],bitmask_list[1])
                        else:
                           temps&=self.frame_analysis_segments(frame,ifo,cWB_conf.bitmask[ifo][dq],cWB_conf.DQ_channel_rate[ifo][dq],bitmask_list[1])
                      s[ifo]|=temps
                    except:
                      bitmask_list=self.getbitmask_fromfile(frame,ifo,cWB_conf.DQ_channel[ifo],cWB_conf.DQ_channel_rate[ifo])
                      print "In frame_analysis_segments: %s %s"%(frame.pfn,cWB_conf.DQ_channel[ifo])
                      s[ifo]|=self.frame_analysis_segments(frame,ifo,cWB_conf.bitmask[ifo],cWB_conf.DQ_channel_rate[ifo],bitmask_list[1])
                    for inj in range(0,len(cWB_conf.inj_name)):
                       try:
                         t_inj_bitmask=cWB_conf.inj_bitmask[ifo][inj]
                       except:
                         t_inj_bitmask=cWB_conf.inj_bitmask[inj]
                       try:
                         inj_seg=self.getinj_frame(frame,ifo,bitmask_list[1],t_inj_bitmask,cWB_conf.DQ_channel_rate[ifo][0],False)
                       except:
                         inj_seg=self.getinj_frame(frame,ifo,bitmask_list[1],t_inj_bitmask,cWB_conf.DQ_channel_rate[ifo],False)
                       inj_sg[inj]|=inj_seg
                    for veto in range(0,cat2veto):
                       #print cWB_conf.veto_channel[ifo][veto]
                       cat2_bitmask=self.getbitmask_fromfile(frame,ifo,cWB_conf.cat2_channel[ifo][veto],cWB_conf.cat2_rate[ifo][veto])
                       #print cat2_bitmask[1]
                       cat2_seg=self.getinj_frame(frame,ifo,cat2_bitmask[1],cWB_conf.cat2_bitmask[ifo][veto],cWB_conf.cat2_rate[ifo][veto],True)
                       #print repr(veto_seg)
                       cat2_sg[veto]|=cat2_seg
                    for veto in range(0,nveto):
                       #print cWB_conf.veto_channel[ifo][veto]
                       veto_bitmask=self.getbitmask_fromfile(frame,ifo,cWB_conf.veto_channel[ifo][veto],cWB_conf.veto_rate[ifo][veto])
                       #print veto_bitmask[1]
                       veto_seg=self.getveto_frame(frame,ifo,veto_bitmask[1],cWB_conf.veto_rate[ifo][veto])
                       #print repr(veto_seg)
                       veto_sg[veto]|=veto_seg
                    self.last_read[ifo]=frame.start
                         
            for inj in range(0,len(cWB_conf.inj_name)):
                fn_global="%s/%s/%s_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,cWB_conf.inj_name[inj],ifo)
                try:
                    f=open(fn_global)
                    tinj_seg=fromsegwizard(f,LIGOTimeGPS)
                    f.close()
                except Exception,e:
                    #print e
                    tinj_seg=segmentlist([])
                tinj_seg|=inj_sg[inj]
                tinj_seg.coalesce()
                fn_global="%s/%s/%s_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,cWB_conf.inj_name[inj],ifo)
                try:
                    f=open(fn_global,"w")
                    tosegwizard(f,tinj_seg,True,LIGOTimeGPS)            
                    f.close()
                except Exception,e:
                    print e
            for veto in range(0,cat2veto):
                fn_global="%s/%s/%s_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,ifo,"CAT2")
                #print fn_global
                try:
                    f=open(fn_global)
                    tcat2_seg=fromsegwizard(f,LIGOTimeGPS)
                    f.close()
                except Exception,e:
                    #print e
                    tcat2_seg=segmentlist([])
                tcat2_seg|=cat2_sg[veto]
                #tcat2_seg.coalesce()
                fn_global="%s/%s/%s_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,ifo,"CAT2")
                try:
                    f=open(fn_global,"w")
                    tosegwizard(f,tcat2_seg,True,LIGOTimeGPS)
                    f.close()
                except Exception,e:
                    print e
            for veto in range(0,nveto):
                fn_global="%s/%s/%s_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,cWB_conf.veto_channel[ifo][veto].replace(":","_"),ifo)
                #print fn_global
                try:
                    f=open(fn_global)
                    tveto_seg=fromsegwizard(f,LIGOTimeGPS)
                    f.close()
                except Exception,e:
                    #print e
                    tveto_seg=segmentlist([])
                tveto_seg|=veto_sg[veto]
                #tveto_seg.coalesce()
                fn_global="%s/%s/%s_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,cWB_conf.veto_channel[ifo][veto].replace(":","_"),ifo)
                try:
                    f=open(fn_global,"w")
                    tosegwizard(f,tveto_seg,True,LIGOTimeGPS)            
                    f.close()
                except Exception,e:
                    print e
            s[ifo].coalesce()
            #fn_current="%s/%s/seg_%s_current.txt"%(cWB_conf.online_dir,"SEGMENTS",ifo)
            fn_current="%s/%s/seg_%s_current.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,ifo)
            #print "fn_current=%s"%fn_current
            try:
                f=open(fn_current,"w")
                tosegwizard(f,s[ifo])            
                f.close()
            except Exception,e:
                print e
            #fn_global="%s/%s/seg_%s_global.txt"%(cWB_conf.online_dir,"SEGMENTS",ifo)
            fn_global="%s/%s/seg_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,ifo)
            #print "fn_global=%s"%fn_global
            try:
                f=open(fn_global)
                sg=fromsegwizard(f)
                f.close()
            except Exception,e:
                print e
                sg=segmentlist([])
            sg|=s[ifo]
            sg.coalesce()
            #fn_global="%s/%s/seg_%s_global.txt"%(cWB_conf.online_dir,"SEGMENTS",ifo)
            fn_global="%s/%s/seg_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,ifo)
            try:
                f=open(fn_global,"w")
                tosegwizard(f,sg)            
                f.close()
            except Exception,e:
                print e
            #setattr(self,"segs_%s"%ifo,sg)
            s[ifo]=segmentlist([])
            s[ifo].append(segment(self.cycle_start-5*cWB_conf.look_back, self.cycle_start))
            s[ifo]&=sg
        try:
            combined=reduce(lambda x,y: (x & y).coalesce(), s.values())
        except Exception,e:
            print "Is it here: e=%s"%(repr(e))
        #combined_off=segmentlist(map(lambda x: segment(x[0], x[1]), segmentlist(filter(lambda y: y[0]< y[1], combined))));
        #combined_off.coalesce() 
        #fn_combined_current="%s/%s/seg_current.txt"%(cWB_conf.online_dir,"SEGMENTS")
        fn_combined_current="%s/%s/seg_current.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir)
        #print "fn_combined_current=%s"%fn_combined_current
        try:
            f=open(fn_combined_current,"w")
            tosegwizard(f,combined)
            f.close()
        except Exception,e:
            print e
            #combined=segmentlist([])
        #fn_combined_global="%s/%s/seg_global.txt"%(cWB_conf.online_dir,"SEGMENTS")
        fn_combined_global="%s/%s/seg_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir)
        #print "fn_combined_global=%s"%fn_combined_global
        try:
            f=open(fn_combined_global)
            combined_global=fromsegwizard(f)
            f.close()
        except Exception,e:
            print e
            combined_global=segmentlist([])
        combined_global |= combined
        combined_global.coalesce()
        try:
            f=open(fn_combined_global,"w")
            tosegwizard(f,combined_global)
            f.close()
        except Exception,e:
            print e
        #self.segs=combined_off
        #self.segs=combined_global & segmentlist([segment(self.cycle_start - cWB_conf.look_back, self.cycle_start)])
        self.segs_global=combined_global

    def getinj_frame(self,frame,ifo,bitmask_list,bitmask,DQ_channel_rate_ifo,opposite):
        start=frame.start
        states=map(lambda x: int(float(x.strip())) & bitmask, bitmask_list.split("\n"))
        if (opposite==False):
          found=map(lambda x: x!=bitmask, states)
        else:
          found=map(lambda x: x==bitmask, states)
        dt=1./DQ_channel_rate_ifo
        tsegs=segmentlist([])
        for i in range(0,len(found)):
            if (found[i]==True):
                start_seg=LIGOTimeGPS(start)+dt*i
                stop_seg=LIGOTimeGPS(start_seg)+dt
                tsegs.append(segment(start_seg,stop_seg))
        tsegs.coalesce()
        return tsegs

    def getveto_frame(self,frame,ifo,veto_bitmask,veto_rate):
        start=frame.start
        states=map(lambda x: int(float(x.strip())), veto_bitmask.split("\n"))
        print "VETO: states=%s"%(repr(states))
        found=map(lambda x: x==1, states)
        dt=1./veto_rate
        #print "DT: %f"%dt
        tsegs=segmentlist([])
        for i in range(0,len(found)):
            if (found[i]==True):
                start_seg=LIGOTimeGPS(start)+dt*i
                stop_seg=LIGOTimeGPS(start_seg)+dt
                tsegs.append(segment(start_seg,stop_seg))
        tsegs.coalesce()
        return tsegs
        
    def getbitmask_fromfile(self,frame_t,ifo,DQ_channel_ifo,DQ_channel_rate_ifo):
        #tmp_bin="%s/SEGMENTS/frame_bin"%(cWB_conf.online_dir)
        tmp_bin="%s/%s/frame_bin"%(cWB_conf.run_dir,cWB_conf.seg_dir)
        #com="%s/LIGOTOOLS/ligotools/bin/FrameDataDump -I%s -C%s -O%s"%(os.environ['HOME_LIBS'], frame_t.pfn,DQ_channel_ifo,tmp_bin)
        com="/usr/bin/env FrameDataDump -I%s -C%s -O%s"%(frame_t.pfn,DQ_channel_ifo,tmp_bin)
        #print com
        a=commands.getstatusoutput(com)
        if(ifo.find("Z1")>=0):
            com="%s/tools/online/bin/extract_float %s %d"%(os.environ['HOME_WAT'],tmp_bin, frame_t.duration*DQ_channel_rate_ifo)
        else:
            com="%s/tools/online/bin/extract_int %s %d"%(os.environ['HOME_WAT'],tmp_bin, frame_t.duration*DQ_channel_rate_ifo)
        #print com
        a=commands.getstatusoutput(com)
        #print a
        return a

    def frame_analysis_segments(self,frame_t,ifo,bitmask_ifo,DQ_channel_rate_ifo,bitmask_list):
        states=map(lambda x:int(float(x.strip())) & bitmask_ifo,bitmask_list.split("\n"))
        print "states=%s"%(repr(states))
        s=segmentlist([])
        rate=DQ_channel_rate_ifo
        nsec=len(states)/rate
        for q in range(0,nsec):
           fstart=frame_t.start+q
           bad=0
           for b in range(0,rate):
             #print "%i %i %i"%(q*rate+b,int(float(bitmask_list[q*rate+b])),states[q*rate+b])
             if (states[q*rate+b]!=bitmask_ifo):
               bad=1
           if bad==0:
              s.append(segment(fstart,fstart+1))
           #print "%i %i"%(fstart,bad)
        s.coalesce()
        return s

    def divide(self,segment):
        pass

if(__name__=="__main__"):
    cWB_conf.run_offset=0
    current_time=get_time()
    print "current_time=%d"%(current_time)
    #run_start=current_time
    #if(cWB_conf.run_start>0):
    try:
        run_start=int(cWB_conf.run_start)
    except:
        run_start=current_time
    print "run_start=%d"%(run_start)
    cWB_conf.run_offset=current_time - run_start
    print "run_offset=%d"%(cWB_conf.run_offset)
    try:
        run_end=int(cWB_conf.run_end)
    except:
        run_end=-1
    print "run_end=%d"%(run_end)

    f=open("offset.txt","w")
    print >>f,cWB_conf.run_offset
    f.close()        

    cwb=cWB(run_end)
    try:
       gracedb_group=cWB_conf.gracedb_group
       cWB_conf.sendtogracedb=True
       #try:
       #  gracedb_client=cWB_conf.gracedb_client
       #except:
       #  print "standard graceDB client"
    except:
       gracedb_group=""
       cWB_conf.sendtogracedb=False
    try:
       nmails=len(cWB_conf.emails)
       cWB_conf.sendmail=True
    except:
       cWB_conf.sendmail=False
    try:
       production_dir=cWB_conf.production_dir 
       pp_redo=True
    except:
       pp_redo=False
    try:
       th_rho_mail=cWB_conf.th_rho_mail
    except:
       th_tho_mail=cWB_conf.th_rho_lum
    cwb.run()

