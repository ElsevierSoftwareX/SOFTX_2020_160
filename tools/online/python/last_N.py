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

import glob, os, sys, commands,pickle
from optparse import OptionParser

import cWB_conf
if (os.environ['SITE_CLUSTER']=="CASCINA"):
  from glue.segments import *
  from glue.segmentsUtils import *
else:
  from ligo.segments import *
  from ligo.segments.utils import *
from rho_histogram import *

def findced(odir):
    lines=open("%s/body.html"%odir).readlines()
    for l in lines:
       if (l.find("href")!=-1 and l.find("ced")!=-1):
         cedstring=l.split("\"")[1]
         launchced(cedstring,odir)

def launchced(cedstring,odir):
   try:
      print cedstring
      line=cedstring.split("_")
      (lag,slag,time)=(int(line[6].replace("lag","")),int(line[5].replace("slag","")),line[9])
      print "%i %i %s"%(lag,slag,time)
      tmpfile=cWB_conf.bkg_par
      tmpfile=tmpfile.replace(".C","_split.C")
      olddir=os.getcwd()
      if (slag>0):
         newdir="%s/%s_%i/%s"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir,slag,time[:5])
      else:
         newdir="%s/%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir,time[:5])
      os.chdir(newdir)
      bkg_dirs=glob.glob("*")
      bkg_dir=filter(lambda x: int(x.split("-")[0])<=float(time) and int(x.split("-")[1])>float(time), bkg_dirs)
      print "%s/%s"%(newdir,bkg_dir[0])
      os.chdir(bkg_dir[0])
      if (not os.path.exists("ced")):
         commands.getstatusoutput("mkdir ced")
      if (cWB_conf.bkg_split>1):
        daglines=open("condor/%s.dag"%(cWB_conf.label)).readlines()
        slag_string=filter(lambda y: y.find("SLAG_SHIFT")!=-1,daglines[1].split(" "))[0].replace("SLAG_SHIFT=","").replace("\n","")
      else:
        daglines=open("condor/%s.sub"%(cWB_conf.label)).readlines()
        slag_string=filter(lambda y: y.find("SLAG_SHIFT")!=-1,daglines[4].split(";"))[0].replace("SLAG_SHIFT=","").replace("\n","")
      if (not os.path.exists("start_ced_%s_%i"%(time,lag))):
         #com="%s/bin/ced.csh %s %s 1 %i %i"%(cWB_conf.bkg_run_dir,tmpfile,slag_string[1:-1],int(float(time)),lag)
         #com="%s/bin/ced.csh %s %s 1 %i %i"%(cWB_conf.bkg_run_dir,tmpfile,slag_string,int(float(time)),lag)
         com="%s/tools/online/bin/ced_bkg.csh %s %s 1 %i %i"%(os.environ['HOME_WAT'],tmpfile,slag_string,int(float(time)),lag)
         print com
         commands.getstatusoutput(com)
         commands.getstatusoutput("touch start_ced_%s_%i"%(time,lag))
      else:
         ced_dir="ced/ced_%s_%s_%s_slag%i_lag%i_1_job1/%s"%(bkg_dir[0].split("-")[0],cWB_conf.bkg_job_duration,bkg_dir[0],0,lag,cedstring.split("/")[len(cedstring.split("/"))-1])
         #print "%s/%s/%s"%(newdir,bkg_dir[0],ced_dir)
         link="mkdir -p %s/%s/ced/%s;ln -s %s/%s/%s %s/%s/ced/%s"%(olddir,odir,cedstring.split("/")[1],newdir,bkg_dir[0],ced_dir,olddir,odir,cedstring.split("/")[1])
         print link
         commands.getstatusoutput(link)
      os.chdir(olddir)
   except:
      print "ced not done"
      os.chdir(olddir)

def merge(path,rootfiles):
	olddir=os.getcwd()
	os.chdir(path)
        namefile="cwb_bkg"
        commands.getstatusoutput("rm -rf %s;mkdir %s"%(namefile,namefile)) 
	os.chdir("%s"%(namefile))
        commands.getstatusoutput("mkdir -p config output condor merge report/postprod")
        version="M1"
	#f=open("merge/merge_%s.%s.lst"%(namefile,version),"w")
	#print >>f,"\n".join(rootfiles)
	#f.flush()
	#f.close()
        for rootfile in rootfiles:
            commands.getstatusoutput("ln -sf %s output/"%(rootfile))
           
        commands.getstatusoutput("cp %s config/user_parameters.C"%(cWB_conf.bkg_par))
        commands.getstatusoutput("cp %s config/user_pparameters.C"%(cWB_conf.pp_par))
        #commands.getstatusoutput("ln -sf %s/template.merged/postprod.csh"%(cWB_conf.run_dir))
        #commands.getstatusoutput("ln -sf %s/template.merged/setcuts.csh"%(cWB_conf.run_dir))
	#commands.getstatusoutput("ln -s %s/template.merged/mergePROD.C"%(cWB_conf.run_dir))	
	#commands.getstatusoutput("echo merge/merge_%s.%s.lst %s.%s | root -l -b -q %s/tools/online/RUN_cWB/template.merged/mergePROD.C"%(namefile,version,namefile,version,os.environ['HOME_WAT']))
        commands.getstatusoutput("""%s/scripts/cwb_merge.csh %s"""%(os.environ['HOME_CWB'],version))
        #RHO THRESHOLD
        #commands.getstatusoutput("./postprod.csh %s"%(version))
        b_version=version
        rho_cut="rho[%i]>%f"%(cWB_conf.id_rho,cWB_conf.th_rho_lum)
        rho_version="Rho"
        #com="""./setcuts.csh %s "%s" "%s" """%(b_version,rho_cut,rho_version)
        com="""%s/scripts/cwb_setcuts.csh %s "%s" "%s" """%(os.environ['HOME_CWB'],b_version,rho_cut,rho_version)
        #print com
        commands.getstatusoutput(com)
        version="%s.C_%s"%(b_version,rho_version)
        #commands.getstatusoutput("./postprod.csh %s"%(version))
        commands.getstatusoutput("%s/scripts/cwb_report.csh %s create"%(os.environ['HOME_CWB'],version))
        pp_dir=glob.glob("report/postprod/*")
        #print pp_dir
        commands.getstatusoutput("rm -rf ../plot")
        commands.getstatusoutput("mv %s ../plot"%(pp_dir[0]))
        commands.getstatusoutput("cp ../plot/data/live.txt ../plot_live.txt")
        commands.getstatusoutput("cp ../plot/data/events_sorted.txt ../plot_events_sorted.txt")
        try:
          production_dir=cWB_conf.production_dir
        except:
          findced("../plot")
        try:
         #T_Cuts_list=cWB_conf.Cuts_list
         #T_Cuts_list.append("OR_cut")
         if (len(cWB_conf.Cuts_list)>1):
            T_Cuts_list=["OR_cut"]
         else:
            T_Cuts_list=[]
         for n in cWB_conf.Cuts_list:
            T_Cuts_list.append(n)
         print repr(T_Cuts_list)
         for n in T_Cuts_list:
           #com="""./setcuts.csh %s '--tcuts %s --label %s' """%(version,n,n)
           com="""%s/scripts/cwb_setcuts.csh %s '--tcuts %s --label %s' """%(os.environ['HOME_CWB'],version,n,n)
           print com
           commands.getstatusoutput(com)
           if (os.path.exists("merge/live_%s.%s.C_%s.root"%(namefile,version,n))):
             #com="./postprod.csh %s.C_%s"%(version,n)
             com="%s/scripts/cwb_report.csh %s.C_%s create"%(os.environ['HOME_CWB'],version,n)
             commands.getstatusoutput(com)
             pp_dir=glob.glob("report/postprod/*")
             print pp_dir
             #commands.getstatusoutput("mv %s plot_%s"%(pp_dir[0],n))
             print commands.getstatusoutput("rm -rf ../plot%s"%n)
             #print commands.getstatusoutput("mkdir -p FOMs/plot_%s"%n)
             print commands.getstatusoutput("mv %s ../plot%s"%(pp_dir[0],n))
             commands.getstatusoutput("cp ../plot%s/data/live.txt ../plot%s_live.txt"%(n,n))
             commands.getstatusoutput("cp ../plot%s/data/events_sorted.txt ../plot%s_events_sorted.txt"%(n,n))
             try:
               production_dir=cWB_conf.production_dir
             except:
               findced("../plot%s"%(n))
           else:
             print "no event in this class"
             print commands.getstatusoutput("rm -rf ../plot%s"%n)
        except:
         print "no class"
        #try:
        #   th_qveto=cWB_conf.th_qveto
        #   qveto_cut="(frequency[0]<1024"
        #   for i in range(0,len(cWB_conf.ifos)):
        #      qveto_cut="%s&&Qveto[%i]>%f&&Qveto[%i]>%f"%(qveto_cut,i,th_qveto,i+len(cWB_conf.ifos),th_qveto)
        #   qveto_cut="%s)||frequency[0]>=1024"%(qveto_cut)
        #   #print qveto_cut
        #   qveto_version="Qveto"
        #   com="""./setcuts.csh %s "%s" "%s" """%(version,qveto_cut,qveto_version)
        #   #print com
        #   commands.getstatusoutput(com)
        #   com="./postprod.csh %s.C_%s"%(version,qveto_version)
        #   #print com
        #   commands.getstatusoutput(com)
        #   pp_dir=glob.glob("report/postprod/*")
        #   #print pp_dir
        #   commands.getstatusoutput("rm -rf ../plot_qveto")
        #   commands.getstatusoutput("mv %s ../plot_qveto"%(pp_dir[0]))
        #except Exception,e:
        #   th_veto=0.
        #   #print e
        commands.getstatusoutput("rm -rf ../merge")
        commands.getstatusoutput("mv merge ../.")
	os.chdir("../")
	commands.getstatusoutput("rm -rf %s"%(namefile))
	os.chdir(olddir)

def significance(dir, rho_th, plotdir):
     try:
        f=open("%s/segments.txt"%(dir))
        #times=filter(lambda x: len(x)>0, f.readlines()[0].strip().split(" "))
        times=filter(lambda x: len(x.strip().split(" "))>0, f.readlines())
        f.close()
        #start=float(times[0])
        #stop=float(times[1])
        start=map(lambda x: float(x.split(" ")[0]),times)
        stop=map(lambda x: float(x.split(" ")[1]),times)
        f=open("%s/%s_live.txt"%(dir,plotdir))
        lines=filter(lambda x:x.find("nonzero")==0,f.readlines())
        f.close()
        nonzero_live = float(lines[0].replace("nonzero lags live=",""))
        f=open("%s/%s_events_sorted.txt"%(dir,plotdir))
        lines=f.readlines()
        f.close()
        n=0
        for l in lines:
            line=filter(lambda x: len(x)>0, l.strip().split(" "))
            (rho,cc,neted,penalty)=(float(line[2]),float(line[3]),float(line[10]),float(line[11]))
            #print "%.3f %.3f %.3f %.3f"%(rho,cc,neted,penalty)
            if (rho>rho_th):
               n=n+1
        try:
          trials=cWB_conf.Cuts_trials
        except:
          try:
            trials=len(cWB_conf.Cuts_list)
          except:
            trials=1
        if (n==0):
            m=1*trials
        else:
            m=n*trials
        #return [n,m/float(nonzero_live),start,stop,float(stop-start)] #n, n/float(self.nonzero_LT)]
        return [n,m/float(nonzero_live),start[0],stop[-1],sum(stop)-sum(start)]
     except:
        return [-1,-1,-1,-1,-1]

if(__name__=="__main__"):

        job_dir="%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir)

        parser = OptionParser()
	parser.add_option("-r","--recompute",dest="recompute",type=int,help="1 - to recompute, 0 - use the existing ones, if available")
        parser.add_option("-n","--nhours",dest="N",type="int", help="compute/use historgram for N last hours")
	parser.add_option("-f","--offset",dest="offset",type="int",help="do not use the last offset hours",default=0)
        parser.add_option("-o","--outputdir",dest="outputdir",type="string")	
	parser.add_option("-R","--rho",dest="rho",type="float")
	#parser.add_option("-m","--mult",dest="m",type="int",default=6)	
	parser.add_option("-m","--mult",dest="mdir",type="int",default=0)	
        parser.add_option("-d","--indir",dest="plotdir",type="string",default="plot")	
	
        (options, args) = parser.parse_args()
	if(not options.N):
		parser.error("N must be specified")
	if(not options.outputdir):
		parser.error("Output directory must be specified")
	if(options.recompute==None):
		parser.error("Recompute must be specified")
	if(not options.rho):
                options.rho=cWB_conf.th_rho_lum
		#parser.error("rho must be specified")
	if(options.offset!=None and options.offset>options.N):
		parser.error("offset must be smaller than N")

        options.m=options.mdir
        if (options.mdir==0):
             options.m=int(3600./cWB_conf.bkg_job_duration+0.5)

	if(options.recompute==0):
	        dir="%s/FOM_%d_%d_%d"%(options.outputdir, options.N, options.offset, options.mdir)			
#		fn="%s/histogram_%d_%d_%d.pickle"%(options.outputdir, options.N, options.offset, options.m)
		try:
#			f=open(fn,"rb")
#			h=pickle.load(f)
#			f.close()
#			z=h.significance(options.rho)
#			print """%d %g %d %d %d"""%(z[0],z[1],h.segs[0][0],h.segs[-1][1],abs(h.segs))
			z=significance(dir,options.rho,options.plotdir)
			print """%d %g %d %d %d"""%(z[0],z[1],z[2],z[3],z[4])
		except Exception,e:
			print e
			options.recompute=1
	if(options.recompute==0):
		sys.exit(0)

	a=glob.glob("%s/?????/*/finished"%job_dir)
	a.sort()

#	while(len(a)>0):
#		tfn=a[-1].replace("finished","OUTPUT.merged/triggers.txt")
#		tf=open(tfn)
#		lines=filter(lambda x: x.find("#")==-1 and len(x.strip())>0, tf.readlines())
#		tf.close()
#		if(len(lines)>0):
#			break
#		else:
#			a.pop()

	#print a
	if(options.offset==0):
		considered_jobs=a[-options.m*options.N:]
	else:
		considered_jobs=a[-options.m*options.N:-options.m*options.offset]
        if (options.N==10000):
                considered_jobs=a[:]
	#print considered_jobs
        if (len(considered_jobs)==0):
            print "no background jobs in desired interval"
            sys.exit(0)
	considered_jobs=map(lambda x: "/".join(x.split("/")[:-1]), considered_jobs)
	considered_segments=segmentlist(map(lambda x: segment(map(int,x.split("/")[-1].split("-"))), considered_jobs));#considered_segments.coalesce()
	considered_livetime=abs(considered_segments)

	total_livetime=0

        rootfiles=[]
	for b in considered_jobs:
                rfiles=glob.glob("%s/output/w*.root"%b)
                rootfiles+=rfiles

	dir="%s/FOM_%d_%d_%d"%(options.outputdir, options.N, options.offset, options.mdir)			
        #print dir
        commands.getstatusoutput("mkdir %s"%(dir))
        merge(dir,rootfiles)
	f=open("%s/segments.txt"%(dir),"w")
	#print >>f,"%s %s"%(considered_segments[0][0],considered_segments[len(considered_segments)-1][1])
        for seg in considered_segments:
          print >>f,"%s %s"%(seg[0],seg[1])
        f.close()

	z=significance(dir,options.rho,options.plotdir)
	print """%d %g %d %d %d"""%(z[0],z[1],z[2],z[3],z[4])

        #print len(considered_jobs)
        #print len(considered_segments)
        #print considered_livetime
