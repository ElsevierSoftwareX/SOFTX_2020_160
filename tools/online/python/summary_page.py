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

from run import *
from run_utils import *
from Trigger import *
import commands, os, sys, glob
if (os.environ['SITE_CLUSTER']=="CASCINA"):
  from glue.segments import *
  from glue.segmentsUtils import *
else:
  from ligo.segments import *
  from ligo.segments.utils import *
from collections import deque
import cWB_conf
import pickle
import os.path
from ligo.gracedb.rest import GraceDb

import matplotlib
matplotlib.use('Agg')
import pylab

#import injection_search

def first_n_digits(gps, n):
    return str(gps)[0:n]

def trigger_2_segment(t,segs):
    for s in segs:
        if(t.time[0]>=s[0] and t.time[0]<=s[1]):
            return s
    return None

def make_empty_page(start,end,rootdir,output_file_name,title,roothtml,current_time,ifar_flag):
        fn_tmp="%s.tmp"%(roothtml)
        ff=open("%s.tmp"%(roothtml),"w")
        print >>ff,"<html>"
        print >>ff,"<title>%s</title>"%(title)
        print >>ff,"<body>"
        print >>ff,"""<h1 align=center><font size=10 face="Courier" color="#DF013A">%s</font></h1>"""%(title)
        t=get_time_nooff()
        print >>ff,"<center><i>The page generated at %d : %s</i></center><p>"%(t,tconvert(t))
        print >>ff,"<hr size=3 noshade color=\"blue\">"
        print >>ff,"<h2 align=center>Run time statistics</h2>"
        print >>ff,"<center>No jobs were running during the specified time interval</center><p>"

        #f = open(cWB_conf.segs1)
        f = open("%s/%s/seg_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir))
        segments = fromsegwizard(f)
        f.close()

        segments &= segmentlist([segment(start,end)]); segments.coalesce()

        segments_fnr="science_segments.txt"
        segments_fn="%s/%s"%(rootdir,segments_fnr)
        fff=open(segments_fn,"w")
        tosegwizard(fff,segments)
        fff.close()
    
        single_segments={}
        segments_fnr={}
        segments_fn={}
        for ifo in cWB_conf.ifos:
            f=open("%s/%s/seg_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,ifo))
            single_segments[ifo]=fromsegwizard(f)
            f.close()
            single_segments[ifo] &=segmentlist([segment(start,end)]);
            single_segments[ifo].coalesce()
            segments_fnr[ifo]="%s_segments.txt"%ifo
            segments_fn[ifo]="%s/%s"%(rootdir,segments_fnr[ifo])
            fff=open(segments_fn[ifo],"w")
            tosegwizard(fff,single_segments[ifo])
            fff.close()

        print >>ff,"<center><table border=1>"
        print >>ff,"<tr><th align=left>current time</th><td>%d : %s</td></tr>"%(current_time,tconvert(current_time))
        print >>ff,"<tr><th align=left>start</th><td>%d : %s</td></tr>"%(start,tconvert(start))
        print >>ff,"<tr><th align=left>end</th><td>%d : %s</td></tr>"%(end,tconvert(end))
        print >>ff,"</table></center>"
        print >>ff,"<h2 align=center>Livetime</h2>"
        print >>ff,"<center><table border=1 cellpadding=5 bgcolor=\"yellow\">"
        print >>ff,"<tr><th align=left><a href=\"%s\">Science segments</a></th><td>%d</td></tr>"%("science_segments.txt",abs(segments))
        for ifo in cWB_conf.ifos:
          print >>ff,"<tr><th align=left><a href=\"%s\">%s</a></th>"%(segments_fnr[ifo],ifo)
          print >>ff,"<td>%i</td></tr>"%(abs(single_segments[ifo]))
        print >>ff,"</table></center>"

        ff.flush()
        ff.close()
        a=commands.getstatusoutput("mv %s %s"%(fn_tmp,roothtml))


def summary_page(start,end,job_rootdir,rootdir,qPrintJobs,qPrintTriggers,output_file_name,title, ifar_flag="", segment_table=False):
    commands.getstatusoutput("mkdir -p %s"%rootdir)
    roothtml="%s/%s"%(rootdir,output_file_name)
    start6=first_n_digits(start,6)
    end6=first_n_digits(end,6)
    big_dirs_N=range(int(start6),int(end6)+1)
    job_dirs=[]
    for N in big_dirs_N:
        job_dirs+=glob.glob(cWB_conf.run_dir+"/"+cWB_conf.jobs_dir+"/%s/*"%N)
    job_dirs=filter(lambda x:int(x.split("/")[-1].split("-")[0])>=start and int(x.split("/")[-1].split("-")[1])<=end,job_dirs)
    job_dirs.sort()
    print "In summary_pages: job_dirs=%i"%(len(job_dirs))

    olddir=os.getcwd()
    bkg_postprod_dir="%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.postprod_dir)
    os.chdir("%s/python"%cWB_conf.online_dir)
    if (ifar_flag!=""):
        option_n=int(ifar_flag.split("_")[1])
        bkg_command = "./last_N.py -n %i -r 1 -o %s"%(option_n,bkg_postprod_dir)
        dodir=1
    else:
        #bkg_command= "./daily_rates.py"
        bkg_command= "./query_rate.py -s %d -e %d -r 6 -d %s/FOM_daily_%d-%d"%(start,end-1,bkg_postprod_dir,start,end-1)
        ifar_flag="FOM_daily_%d-%d"%(start,end-1)
        dodir=0
    print bkg_command
    commands.getstatusoutput(bkg_command)
    os.chdir(olddir)

    current_time=get_time_nooff()

    if(len(job_dirs)==0):
        make_empty_page(start,end,rootdir,output_file_name,title,roothtml,current_time,ifar_flag)
        return
            
    current_time=get_time()

    jobs=map(lambda x: x+"/job.pickle",job_dirs)
    job_dirs=None

    delay_launch_completion=[]
    delay_launch_start=[]
    delay_start_completion=[]
    delay_start_send=[]
    job_segments=segmentlist([])
    running_jobs=[]
    failed_jobs=[]
    jobs_str=[]

    JOBS=[]
    rootfiles=[]
    wavefiles=[]
    livefiles=[]

    n=0

    for job in jobs:
      try:
        #print "job=%s"%(repr(job))
        f=open(job,"rb")
        j=pickle.load(f)
        #if (j.completion_time<0):
        #    old_pickle="%s/job.pickle"%(j.dir)
        #    f2=open(old_pickle,"rb")
        #    j2=pickle.load(f2)
        #    f2.close()
        #    com="cp %s %s"%(old_pickle,job)
        #    commands.getstatusoutput(com)
        #    f.close()
        #    f=open(job,"rb")
        #    j=pickle.load(f)
        #print "job after pickle"
        #JOBS.append(j)
        f.close()
        n+=len(j.triggers_all)
        #if (j.triggers_all>0):
        #  JOBS.append(job)
        pick_file=job
        if(j.status==2):
            jdir=cWB_conf.web_link+"/"+"/".join(j.dir.split("/")[-3:])
            p_triggers=filter(lambda x: x.topublish==True, j.triggers_all)
            if (len(p_triggers)>0):
              JOBS.append(pick_file)
            if (len(j.triggers_all)-len(p_triggers) > 0):
               s_tr="    [%s]"%(len(j.triggers_all)-len(p_triggers))
            else:
               s_tr=""
            jobs_str.append("<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d%s</td><td><a href=\"%s/job.html\">link</a></td></tr>"%\
                  #(tconvert(job.start),job.start,job.end,job.launch_time,job.completion_time,job.completion_time - job.launch_time,job.launch_time-job.start, len(p_triggers),s_tr,jdir))
                  (tconvert(j.start),j.start,j.end,j.launch_time,j.completion_time,j.completion_time - j.launch_time,j.launch_time-j.end, len(p_triggers),s_tr,jdir))
            #delay_launch_completion.append(j.delay_launch_completion)
            delay_launch_completion.append(j.completion_time - j.launch_time)
            delay_launch_start.append(j.launch_time - j.end)#start)
            delay_start_completion.append(j.completion_time - j.end)#start - cWB_conf.seg_duration/2.)
            delay_start_send.append(j.delay_start_send)
            job_segments.append(segment(j.start,j.end))
        elif(j.status==3 or j.status==4):
            failed_jobs.append(j)
        else:
            if(j.end+cWB_conf.job_timeout > current_time):
                print "j.end=%d timeout=%d sum=%d current_time=%d"%(j.end, cWB_conf.job_timeout, j.end + cWB_conf.job_timeout, current_time)
                running_jobs.append(j)
        rfiles=glob.glob("%s/OUTPUT/w*.root"%j.dir)
        rootfiles+=rfiles
        rfiles=glob.glob("%s/OUTPUT.merged/trigger_*.root"%j.dir)
        wavefiles+=rfiles
        rfiles=glob.glob("%s/OUTPUT.merged/live.root"%j.dir)
        livefiles+=rfiles
      except:
        pass
    #jobs=None
    list_job_segments_fnr="list_job_segments.txt"
    list_job_segments_fn="%s/%s"%(rootdir,list_job_segments_fnr)
    if (not os.path.exists(list_job_segments_fn)):
       dodir=1
    if(dodir==0 and os.path.exists(list_job_segments_fn)):
           linesfile=len(open(list_job_segments_fn).readlines())
           print "linesfile:%s %i, jobs:%i"%(list_job_segments_fn,linesfile,len(job_segments))
           if (len(job_segments)>=linesfile):
             dodir=1
    if (dodir==0):
      print "dir already done"
      return
    fff=open(list_job_segments_fn,"w")
    tosegwizard(fff,job_segments)
    fff.close()
    job_segments.coalesce()

    if(qPrintTriggers and n > 0):
        print "before FOMs"
        plot_FOMs(rootfiles,wavefiles,livefiles,rootdir,start,end)
        print "after FOMs"
    del rootfiles[:]
    del wavefiles[:]
    del livefiles[:]

    fn_tmp="%s.tmp"%(roothtml)
    ff=open("%s.tmp"%(roothtml),"w")
    print >>ff,"<html>"
    print >>ff,"""<head>
<link rel="stylesheet" type="text/css" href="/~vedovato/waveburst/ced-1.0/scripts/shadowbox-3.0.3/shadowbox.css">
<script type="text/javascript" src="/~vedovato/waveburst/ced-1.0/scripts/shadowbox-3.0.3/shadowbox.js"></script>
<script type="text/javascript">
Shadowbox.init();
</script>
<script type="text/javascript" src="%s/tabber.js"></script>                                                         
<link rel="stylesheet" href="%s/tabber.css" TYPE="text/css" MEDIA="screen">
</head>"""%(cWB_conf.web_link,cWB_conf.web_link)
    print >>ff,"<title>%s</title>"%(title)
    print >>ff,"<body>"
    #print >>ff,"""<h1 align=center><font size=10 face="Courier" color="#FF8000">%s</font></h1>"""%(title)
    print >>ff,"""<h1 align=center><font size=10 face="Courier" color="#DF013A">%s</font></h1>"""%(title)
    t=get_time_nooff()
    print >>ff,"<center><i>The page generated at %d : %s</i></center><p>"%(t,tconvert(t))

    print >>ff,"<hr size=3 noshade color=\"blue\">"
    #print >>ff,"<hr size=3 noshade color=\"blue\">"

    max_delay_launch_completion=-1
    min_delay_launch_completion=-1
    avg_delay_launch_completion=-1
    jobs_completed_less_ten=-1
    if(len(delay_launch_completion)>0):
        #print "delay_launch_completion=%s"%repr(delay_launch_completion)
        max_delay_launch_completion=max(delay_launch_completion)
        tmp_delay_launch_completion=filter(lambda x:x>20, delay_launch_completion)
        jobs_completed_less_ten=len(filter(lambda x:x<=600, delay_launch_completion))/float(len(delay_launch_completion))*100
        try:
            min_delay_launch_completion=min(tmp_delay_launch_completion)
            avg_delay_launch_completion=sum(tmp_delay_launch_completion)/float(len(tmp_delay_launch_completion))
        except Exception,e:
            print "Got exception %s"%(repr(e))
            print "tmp_delay_launch_completion=%s"%(repr(tmp_delay_launch_completion))
            min_delay_launch_completion=-1
            avg_delay_launch_completion=-1

    run_statistics_fn="%s/run_statistics.html"%(rootdir) # is that the right directory?
    run_statistics_fnr="run_statistics.html"

    #bis=injection_search.injection_search(cWB_conf.injections_pickle)

    try:
      nclass=len(cWB_conf.Cuts_list)
      if (nclass>1):
         dir_bkg="OR_cut"
      else:
         dir_bkg=cWB_conf.Cuts_list[0]
    except:
      dir_bkg=""
    if (os.path.exists("%s/TIME_SHIFTS/POSTPRODUCTION/%s/plot%s/index.html"%(cWB_conf.online_dir,ifar_flag,dir_bkg))):
      bkgready=True
    else:
      bkgready=False

    if (bkgready==True):
      stat_bkg_title="<th>Background</th>"
      stat_bkg_content="""
    <td>
    <a href="%s/POSTPRODUCTION/%s/segments.txt" target=_new>Segment list</a><br>
    <a href="%s/POSTPRODUCTION/%s/plot%s/data/job_status.html" target=_new>Running time</a>
    </td>"""%(cWB_conf.web_link,ifar_flag,cWB_conf.web_link,ifar_flag,dir_bkg)
    else:
      stat_bkg_title=""
      stat_bkg_content=""
    #<h2 align=center>Run statistics</h2>
    stat_table="""
    <center>
    <table border=1 cellpadding=5 bgcolor=\"yellow\">
    <tr>
    <th>Zero lag jobs</th>
    <th>Zero lag run time (s)</th>
    %s
    </tr>
    <tr align=\"center\">
    <td><table>
    <tr><td align="left"><a href=\"%s\" target=_new>Job list</a></td><td align="right">%d</td>
    <tr><td align="left"><a href=\"%s\" target=_new>Segment list</a></td><td align="right">%d s</td>
    <tr><td align="left"><a href=\"%s\" target=_new>More details</a></td></tr>
    </table></td>
    <td><table>
    <tr><td align=\"left\">Maximum</td><td align=\"right\">%d</td></tr>
    <tr><td align=\"left\">Minimum</td><td align=\"right\">%d</td></tr>
    <tr><td align=\"left\">Average</td><td align=\"right\">%.2f</td></tr>
    <tr><td align=\"left\">Jobs completed within 10 minutes</td><td align=\"right\">%.2f</td></tr>
    </table></td>
    %s
    </tr>
    </table>
    </center><p>
    """%(stat_bkg_title,"list_job_segments.txt",len(delay_launch_completion),"job_segments.txt",abs(job_segments),run_statistics_fnr,max_delay_launch_completion,min_delay_launch_completion,avg_delay_launch_completion, jobs_completed_less_ten,stat_bkg_content)
    print >>ff, "<h3 align=center>General information</h3>"
    print >>ff, stat_table

    if(segment_table):
        link="%s/INVESTIGATIONS/SEGMENTS/%s-%s"%(cWB_conf.web_link,start,end)
        print >>ff,"""<p><a href="%s" target=_new>Segments<ga><p>"""%(link)
        print >>ff,"""<p><a href="comments.html" target=_new>Comments</a><p>"""

    run_statistics(stat_table,delay_launch_completion,delay_launch_start,delay_start_completion,rootdir,start,end, run_statistics_fn, qPrintJobs, title, output_file_name,job_segments, running_jobs, qPrintTriggers, jobs_str)
    del delay_launch_completion[:]
    del delay_launch_start[:]
    del delay_start_completion[:]
    del delay_start_send[:]
    del job_segments[:]
    del running_jobs[:]
    del failed_jobs[:]
    del jobs_str

    seg_table="""
    <table border=1>
    <tr><th>Network</th><th>Segment type</th><th>Livetime, hours</th><th>Duty cycle, \%</th></tr>
    <tr><td>L1H1</td><td>Processed by cWB</td><td><a href=\"\">N</a></td><td>N</td></tr>
    <tr><td>L1H1</td><td>Up</td><td><a href=\"\">N</a></td><td>N</td></tr>
    <tr><td>L1H1</td><td>cWB - UP</td><td><a href=\"\">N hours</a></td><td>N</td></tr>
    </table>
    """

    print "%i triggers"%(n)
    triggers=[]
    #ttlines=[]
    for j in JOBS:
        f=open(j,"rb")
        job=pickle.load(f)
        for t in job.triggers_all:
         if (t.topublish==True):
            t.job_start=job.start
            t.job_dir=job.dir
            #t.error_region=[]
            t.lines=[]
            triggers.append(t)
        #triggers+=job.triggers_all
        f.close()
        #try:
        #    ttlines+=open("%s/OUTPUT.merged/TRIGGERS/triggers.txt"%(job.dir)).readlines()
        #except:
        #    print "No %s/OUTPUT.merged/TRIGGERS/triggers.txt yet"%(job.dir)
    #triggers=filter(lambda x: x.topublish==True, triggers)
    print "qPrintTriggers=%s ntriggers=%d"%(qPrintTriggers, len(triggers))
    #tfo=open("%s/triggers.txt"%(rootdir),"w")
    #print >>tfo,"".join(ttlines)
    #tfo.flush()
    #tfo.close()

    #del ttlines
    del JOBS[:]

    injection_triggers=[]
    #for name in cWB_conf.inj_name:
         #injection_triggers+=filter(lambda x: x.cat_passed_html.find("INJECTION_%s"%name)==0,triggers)
    injection_triggers+=filter(lambda x: x.inj_found==True,triggers)
    print len(injection_triggers)
    for name in cWB_conf.inj_name:
        #triggers=filter(lambda x: x.cat_passed_html.find("INJECTION_%s"%name)==-1,triggers)
        triggers=filter(lambda x: x.inj_found==False,triggers)
    print len(triggers)

    string_miss=""
    miss_file="%s/%s/missing.txt"%(bkg_postprod_dir,ifar_flag)
    if(os.path.exists(miss_file)):
        miss=open(miss_file)
        missing_segments = fromsegwizard(miss); missing_segments.coalesce()
        miss.close()
        if(len(missing_segments)>0):
            string_miss="<br><font color=\"red\"><a href=\"%s/POSTPRODUCTION/%s/missing.txt\" target=_new>Some data have still not processed</a></font>"%(cWB_conf.web_link,ifar_flag)

    if(qPrintTriggers and (len(triggers)+len(injection_triggers)) > 0):

        gw_candidates=filter(lambda x: x.selected==4, triggers)
        lumin_triggers=filter(lambda x: x.selected>=3, triggers)
        coherent_triggers=filter(lambda x: x.selected>=1, triggers)
        for ifo in cWB_conf.ifos:
            for category in range(0,4):
                try:
                    lumin_triggers=filter(lambda x: x.cat_passed[category][ifo], lumin_triggers)
                except:
                    pass
        #for name in cWB_conf.inj_name:
        #    lumin_triggers=filter(lambda x: x.cat_passed_html.find("INJECTION_%s"%name)==-1,lumin_triggers)
        #    gw_candidates=filter(lambda x: x.cat_passed_html.find("INJECTION_%s"%name)==-1,gw_candidates)
        if(len(gw_candidates)>0):
            trigger_table_gw_fn="%s/triggers_gw.html"%(rootdir)
            trigger_table_gw_fnr="triggers_gw.html"
            trigger_table(trigger_table_gw_fn, gw_candidates, "GW event candidates",injection_triggers)        
            gw_candidates.sort(cmp=t_cmp_threshold, reverse=False)
            trigger_table_gw_r_fn="%s/triggers_gw_r.html"%(rootdir)
            trigger_table_gw_r_fnr="triggers_gw_r.html"
            trigger_table(trigger_table_gw_r_fn, gw_candidates, "GW event candidates sorted by rho",injection_triggers)

        if(len(lumin_triggers)>0):
            trigger_table_lumin_fn="%s/triggers_lumin.html"%(rootdir)
            trigger_table_lumin_fnr="triggers_lumin.html"
            trigger_table(trigger_table_lumin_fn, lumin_triggers, "GraceDB triggers after DQs and vetoes",injection_triggers)
            lumin_triggers.sort(cmp=t_cmp_threshold, reverse=False)
            trigger_table_lumin_r_fn="%s/triggers_lumin_r.html"%(rootdir)
            trigger_table_lumin_r_fnr="triggers_lumin_r.html"
            trigger_table(trigger_table_lumin_r_fn, lumin_triggers, "GraceDB triggers sorted by rho",injection_triggers)

        if(len(injection_triggers)>0):
            trigger_table_injection_fn="%s/triggers_injection.html"%(rootdir)
            trigger_table_injection_fnr="triggers_injection.html"
            trigger_table(trigger_table_injection_fn, injection_triggers, "Hardware injections",injection_triggers)

        if(len(triggers)>0):
            trigger_table_t_fn="%s/triggers_t.html"%(rootdir)
            trigger_table_t_fnr="triggers_t.html"
            trigger_table(trigger_table_t_fn, triggers, "Triggers sorted by time",injection_triggers,True)
            triggers.sort(cmp=t_cmp_threshold, reverse=False)
            trigger_table_r_fn="%s/triggers_r.html"%(rootdir)
            trigger_table_r_fnr="triggers_r.html"
            trigger_table(trigger_table_r_fn, triggers, "Triggers sorted by rho",injection_triggers)

        if(len(coherent_triggers)>0):
            coherent_table_t_fn="%s/coherent_t.html"%(rootdir)
            coherent_table_t_fnr="coherent_t.html"
            trigger_table(coherent_table_t_fn, coherent_triggers, "Triggers sorted by time",injection_triggers,True)
            coherent_triggers.sort(cmp=t_cmp_threshold, reverse=False)
            coherent_table_r_fn="%s/coherent_r.html"%(rootdir)
            coherent_table_r_fnr="coherent_r.html"
            trigger_table(coherent_table_r_fn, coherent_triggers, "Triggers sorted by rho",injection_triggers)

       # print >>ff, "<h2 align=center>%d triggers</h2>"%(len(triggers))
        print >>ff, "<center>"
        print >>ff, "<h3 align=center>Triggers</h3>"
        print >>ff, "<table border=1 cellpadding=5>"
        print >>ff,"<tr align=\"center\">"
        print >>ff,"<td><b>All triggers</b><br>No selection</td>"
        try:
          nclass=len(cWB_conf.Cuts_list)
          print >>ff,"<td bgcolor=\"lightslategray\"><b>Selected triggers</b><br>%s</td>"%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)))
        except:
          print >>ff,"<td bgcolor=\"lightslategray\"><b>Selected triggers</b><br>netcc>=%.1f</td>"%(cWB_conf.th_cc)
        print >>ff,"<td bgcolor=\"red\"><b>Gracedb triggers</b><br>(Selected & rho>=%.1f)</td>"%(cWB_conf.th_rho_lum)
        try:
          print >>ff,"<td bgcolor=\"chartreuse\"><b>GW candidates</b><br>(Selected & FAR<%.2e)</td>"%(cWB_conf.th_far_off)
        except:
          print >>ff,"<td bgcolor=\"chartreuse\"><b>GW candidates</b><br>(Selected & rho>=%.1f)</td>"%(cWB_conf.th_rho_off)
        print >>ff,"<td bgcolor=\"cyan\">Injections</td>"
        print >>ff,"</tr>"
        print >>ff,"<tr align=\"center\">"
        print >>ff,"<td><b>%d triggers</b></a><br>"%(len(triggers))
        if (len(triggers)>0):
          print >>ff,"<a href=\"%s\" target=_new>sorted by time</a><br>"%(trigger_table_t_fnr)
          print >>ff,"<a href=\"%s\" target=_new>sorted by rho</a>"%(trigger_table_r_fnr)
        #print >>ff,"<br><a href=\"triggers.txt\" target=_new>in text format</a></td>"
        if (len(coherent_triggers)>0):
          print >>ff,"<td bgcolor=\"lightslategray\"><b>%d triggers</b></a><br>"%(len(coherent_triggers))
          print >>ff,"<a href=\"%s\" target=_new>sorted by time</a><br>"%(coherent_table_t_fnr)
          print >>ff,"<a href=\"%s\" target=_new>sorted by rho</a>"%(coherent_table_r_fnr)
        #print >>ff,"<br><a href=\"triggers.txt\" target=_new>in text format</a></td>"
        else:
          print >>ff,"<td bgcolor=\"lightslategray\"><b>no triggers</b><br>on the given time interval</td>"
        if(len(lumin_triggers)==0):
             print >>ff,"<td bgcolor=\"red\"><b>no triggers</b><br>on the given time interval</td>"
        else:
            print >>ff,"<td bgcolor=\"red\"><b>%d triggers</b><br><a href=\"%s\" target=_new>sorted by time</a><br>"%(len(lumin_triggers),trigger_table_lumin_fnr)
            print >>ff,"<a href=\"%s\" target=_new>sorted by rho</a></td>"%(trigger_table_lumin_r_fnr)
        if(len(gw_candidates)==0):
             print >>ff,"<td bgcolor=\"chartreuse\"><b>no candidates</b><br>on the given time interval</td>"
        else:
            print >>ff,"<td bgcolor=\"chartreuse\"><b>%d candidates</b><br><a href=\"%s\" target=_new>sorted by time</a><br>"%(len(gw_candidates),trigger_table_gw_fnr)
            print >>ff,"<a href=\"%s\" target=_new>sorted by rho</a></td>"%(trigger_table_gw_r_fnr)        
        if(len(injection_triggers)>0):
           print >>ff,"<td bgcolor=\"cyan\"><b>%d</b><br><a href=%s target=_new>Hardware injections</a></td>"%(len(injection_triggers),trigger_table_injection_fnr)
        else:
           print >>ff,"<td bgcolor=\"cyan\">No Hardware injections</td>"
        print >>ff,"</tr>"
        print >>ff,"</table>"
        print >>ff, "</center>"

        #add loudest
        if (len(triggers)>0):
          if (len(coherent_triggers)>0):
            #coherent_triggers.sort(cmp=t_cmp_threshold, reverse=False)
            t=coherent_triggers[0]
          else:
            t=triggers[0]
          print >>ff, "<h3 align=center>Loudest event</h3>"
          print >>ff, "<center>"
          print >>ff,"<table border=1>"
          header=t.table_header(cWB_conf.ifos,False)
          print >>ff, header
          jdir=cWB_conf.web_link+"/"+"/".join(t.job_dir.split("/")[-3:])
          print >>ff, t.to_html_table_row(-1)
          print >>ff, "</table>"
          print >>ff, "</center>"

        print >>ff, "<h3 align=center>Standard cWB pages</h3>"
        print >>ff, "<center>"
        print >>ff,"<table border=1 cellpadding=5 bgcolor=\"#F5A9E1\">"
        print >>ff,"<tr align=\"center\">"
        try:
            nclass=len(cWB_conf.Cuts_list)
            print >>ff,"""<th><a href="%s/Cuts.hh.html" target=_new>Class</a></th>"""%(cWB_conf.web_link)
        except:
            print >>ff,"Class"
        print >>ff,"<th bgcolor=white></th>"
        print >>ff,"<th colspan=2>Zero lag</th>"
        if (bkgready==True):
           print >>ff,"<th bgcolor=white></th>"
           print >>ff,"<th colspan=2>Background</th>"
        print >>ff,"</tr>"
        if(len(triggers) > 0):
           try:
             nclass=len(cWB_conf.Cuts_list)
             print >>ff,"<tr bgcolor=white align=\"center\">"
             print >>ff,"<td><i>Detchar</i></td>"
             print >>ff,"<td bgcolor=white></td>"
             print >>ff,"<td><a href=\"FOMs/plot\" target=_new>Link</a></td>"
             print >>ff,"<td><a href=\"FOMs/plot/data/EVENTS.txt\" target=_new>EVENTS</a></td>"
             if (bkgready==True):
                print >>ff,"<td bgcolor=white></td>"
                print >>ff,"<td><a href=\"%s/POSTPRODUCTION/%s/plot\" target=_new>Link</a></td>"%(cWB_conf.web_link,ifar_flag)
                print >>ff,"<td><a href=\"%s/POSTPRODUCTION/%s/plot/data/EVENTS.txt\" target=_new>EVENTS</a></td>"%(cWB_conf.web_link,ifar_flag)
             print >>ff,"</tr>"
             if (nclass>1):
               print >>ff,"<tr align=\"center\">"
               #print >>ff,"""<td rowspan=%i>Zero lag</td>"""%(nclass+1)
               #print >>ff,"<td><i>Selected</i></td>"
               print >>ff,"<td><i>%s</i></td>"%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)))
               print >>ff,"<td bgcolor=white></td>"
               print >>ff,"<td><a href=\"FOMs/plotOR_cut\" target=_new>Link</a></td>"
               print >>ff,"<td><a href=\"FOMs/plotOR_cut/data/EVENTS.txt\" target=_new>EVENTS</a></td>"
               if (bkgready==True):
                 print >>ff,"<td bgcolor=white></td>"
                 print >>ff,"<td><a href=\"%s/POSTPRODUCTION/%s/plotOR_cut\" target=_new>Link</a></td>"%(cWB_conf.web_link,ifar_flag)
                 print >>ff,"<td><a href=\"%s/POSTPRODUCTION/%s/plotOR_cut/data/EVENTS.txt\" target=_new>EVENTS</a></td>"%(cWB_conf.web_link,ifar_flag)
             #print >>ff,"</tr>"
             for i in range(len(cWB_conf.Cuts_list)):
              n=cWB_conf.Cuts_list[i]
              try:
                name=" (%s)"%(cWB_conf.Cuts_name[i])
              except:
                name=""
              print >>ff,"<tr align=\"center\">"
              print >>ff,"<td><i>%s</i>%s</td>"%(n.replace("_cut",""),name)
              print >>ff,"<td bgcolor=white></td>"
              print >>ff,"<td><a href=\"FOMs/plot%s\" target=_new>Link</a></td>"%(n)
              print >>ff,"<td><a href=\"FOMs/plot%s/data/EVENTS.txt\" target=_new>EVENTS</a></td>"%(n)
              if (bkgready==True):
                 print >>ff,"<td bgcolor=white></td>"
                 print >>ff,"<td><a href=\"%s/POSTPRODUCTION/%s/plot%s\" target=_new>Link</a></td>"%(cWB_conf.web_link,ifar_flag,n)
                 print >>ff,"<td><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/EVENTS.txt\" target=_new>EVENTS</a></td>"%(cWB_conf.web_link,ifar_flag,n)
             print >>ff,"</tr>"
           except Exception,e:
             print >>ff,"<tr align=\"center\">"
             print >>ff,"<td><i>Selected</i></td>"
             print >>ff,"<td bgcolor=white></td>"
             print >>ff,"<td><a href=\"FOMs/plot\" target=_new>Link</a></td>"
             print >>ff,"<td><a href=\"FOMs/plot/data/EVENTS.txt\" target=_new>EVENTS</a></td>"
             if (bkgready==True):
                print >>ff,"<td bgcolor=white></td>"
                print >>ff,"<td><a href=\"%s/POSTPRODUCTION/%s/plot\" target=_new>Link</a></td>"%(cWB_conf.web_link,ifar_flag)
                print >>ff,"<td><a href=\"%s/POSTPRODUCTION/%s/plot/data/EVENTS.txt\" target=_new>EVENTS</a></td>"%(cWB_conf.web_link,ifar_flag)
             print >>ff,"</tr>"
        print >>ff, "</table>"
        print >>ff, "</center>"

        print >>ff,"<hr size=5 noshade color=\"white\">"
        print >>ff,"<hr width=90% size=3 noshade color=\"blue\">"
        print >>ff,"<hr width=90% size=3 noshade color=\"blue\">"

        plot_zero=True
        tabber=False
        try:
          nclass=len(cWB_conf.Cuts_list)
          if (nclass>1):
            dir_bkg="OR_cut"
            tabber=True
          else:
            dir_bkg=cWB_conf.Cuts_list[0]
          if (len(coherent_triggers)==0):
           plot_zero=False
        except:
          dir_bkg="" 
 
        print >>ff, "<h2 align=center><font size=8 face=\"Courier\"  color=\"#DF013A\">Figures of merit</font></h2>"
        print >>ff, "<center>"
        if (plot_zero):
          print >>ff, "<h3 align=center><font size=6 face=\"Courier\"  color=\"#01DF01\">Zero lag</font></h3>"
          print >>ff,"<table border=0 cellspacing = 20 width=60% bgcolor=\"#F2F2F2\">"
          print >>ff,"<tr><th colspan=\"2\">Ranking statistic vs time and frequency</th></tr>"
          print >>ff,"<tr>"
          if (tabber==False):
             print >>ff,"<td align=\"left\"><a href=\"FOMs/plot%s/data/rho_time.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/rho_time.gif\" width=400></a></td><td align=\"right\"><a href=\"FOMs/plot%s/data/rho_frequency.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/rho_frequency.gif\" width=400></a></td>"%(dir_bkg,dir_bkg,dir_bkg,dir_bkg)
          else:
             print >>ff,"<td align=\"center\">"
             print >>ff,"<div class=\"tabber\">"
             print >>ff,"<div class=\"tabbertab\"><h2>%s</h2><a href=\"FOMs/plot%s/data/rho_time.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/rho_time.gif\" width=400></a></div>"%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)),dir_bkg,dir_bkg)
             for qq in range(len(cWB_conf.Cuts_list)):
               print >>ff,"<div class=\"tabbertab\"><h2>%s</h2><a href=\"FOMs/plot%s/data/rho_time.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/rho_time.gif\" width=400></a></div>"%(cWB_conf.Cuts_list[qq].replace("_cut",""),cWB_conf.Cuts_list[qq],cWB_conf.Cuts_list[qq])
             print >>ff,"</div>"
             print >>ff,"</td>"
             print >>ff,"<td align=\"center\">"
             print >>ff,"<div class=\"tabber\">"
             print >>ff,"<div class=\"tabbertab\"><h2>%s</h2><a href=\"FOMs/plot%s/data/rho_frequency.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/rho_frequency.gif\" width=400></a></div>"%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)),dir_bkg,dir_bkg)
             for qq in range(len(cWB_conf.Cuts_list)):
               print >>ff,"<div class=\"tabbertab\"><h2>%s</h2><a href=\"FOMs/plot%s/data/rho_frequency.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/rho_frequency.gif\" width=400></a></div>"%(cWB_conf.Cuts_list[qq].replace("_cut",""),cWB_conf.Cuts_list[qq],cWB_conf.Cuts_list[qq])
             print >>ff,"</div>"
             print >>ff,"</td>"
          print >>ff,"</tr>"
          print >>ff,"<tr align=\"justify\"><td colspan=\"2\"><i>Ranking statistic as a function of gps time (left) and estimated frequency (right). The rho parameter is a approximately equal to median SNR per detector. In green GW candidates triggers, in red triggers send to gracedb, in black all the others.</i></td></tr>"
          print >>ff,"</table><p>"
          print >>ff,"<table border=0 cellspacing = 20 width=60% bgcolor=\"#F2F2F2\">"
          print >>ff,"<table border=0 cellspacing = 20 width=60% bgcolor=\"#F2F2F2\">"
          print >>ff,"<tr>"
          if (tabber==False):
              print >>ff,"<td align=\"left\"><a href=\"FOMs/plot%s/data/ra_vs_dec_online.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/ra_vs_dec_online.gif\" width=400></a></td><td align=\"right\"><a href=\"FOMs/plot%s/data/phi_vs_theta_online.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/phi_vs_theta_online.gif\" width=400></a></td>"%(dir_bkg,dir_bkg,dir_bkg,dir_bkg)
          else:
             print >>ff,"<td align=\"center\">"
             print >>ff,"<div class=\"tabber\">"
             print >>ff,"<div class=\"tabbertab\"><h2>%s</h2><a href=\"FOMs/plot%s/data/ra_vs_dec_online.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/ra_vs_dec_online.gif\" width=400></a></div>"%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)),dir_bkg,dir_bkg)
             for qq in range(len(cWB_conf.Cuts_list)):
               print >>ff,"<div class=\"tabbertab\"><h2>%s</h2><a href=\"FOMs/plot%s/data/ra_vs_dec_online.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/ra_vs_dec_online.gif\" width=400></a></div>"%(cWB_conf.Cuts_list[qq].replace("_cut",""),cWB_conf.Cuts_list[qq],cWB_conf.Cuts_list[qq])
             print >>ff,"</div>"
             print >>ff,"</td>"
             print >>ff,"<td align=\"center\">"
             print >>ff,"<div class=\"tabber\">"
             print >>ff,"<div class=\"tabbertab\"><h2>%s</h2><a href=\"FOMs/plot%s/data/phi_vs_theta_online.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/phi_vs_theta_online.gif\" width=400></a></div>"%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)),dir_bkg,dir_bkg)
             for qq in range(len(cWB_conf.Cuts_list)):
               print >>ff,"<div class=\"tabbertab\"><h2>%s</h2><a href=\"FOMs/plot%s/data/phi_vs_theta_online.gif\" rel=\"shadowbox[gallery]\"><img src=\"FOMs/plot%s/data/phi_vs_theta_online.gif\" width=400></a></div>"%(cWB_conf.Cuts_list[qq].replace("_cut",""),cWB_conf.Cuts_list[qq],cWB_conf.Cuts_list[qq])
             print >>ff,"</div>"
             print >>ff,"</td>"
          print >>ff,"</tr>"
          print >>ff,"<tr align=\"justify\"><td colspan=\"2\"><i>Estimated sky coordinates for each trigger. Left: right ascension and declination, right: latitude and longitude. Colors refer to the number of events inside each bin.</i></td></tr>"
          print >>ff,"</table>"
        if (bkgready==True):
           write_bkg(ff,ifar_flag,string_miss,False)
    else:
        print >>ff,"<h2 align=center>No zero lag triggers</h2>"
        print >>ff, "<center>"
        if (bkgready==True):
           write_bkg(ff,ifar_flag,string_miss)
    print >>ff,"</body>"
    print >>ff,"</html>"

    ff.flush()
    ff.close()

    a=commands.getstatusoutput("mv %s %s"%(fn_tmp,roothtml))

def write_bkg(ff,ifar_flag,string_miss="",putlink=True):
        tabber=False
        try:
          nclass=len(cWB_conf.Cuts_list)
          if (nclass>1):
            dir_bkg="OR_cut"
            tabber=True
          else:
            dir_bkg=cWB_conf.Cuts_list[0]
        except:
          dir_bkg=""
        print >>ff,"<center>"
        print >>ff, "<h3 align=center><font size=6 face=\"Courier\"  color=\"#01DF01\">Background</font>%s</h3>"%(string_miss)
        print >>ff,"<table border=0 cellspacing = 20 width=40% bgcolor=\"#F2F2F2\">"
        if (putlink==True):
             print >>ff,"<tr><th><a href=\"%s/POSTPRODUCTION/%s/plot%s\" target=_new>Standard web page</a></th></tr>"%(cWB_conf.web_link,ifar_flag,dir_bkg)
        print >>ff,"<tr><th align=\"center\">Cumulative rate vs ranking statistic</th></tr>"
        print >>ff,"<tr><th align=\"center\">"
        if (tabber==False):
           print >>ff,"<a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rate_threshold.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rate_threshold.gif\" width=400></a>"%(cWB_conf.web_link,ifar_flag,dir_bkg,cWB_conf.web_link,ifar_flag,dir_bkg)
        else:
           print >>ff,"""<div class="tabber">"""
           print >>ff,"""<div class="tabbertab"><h2>%s</h2><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rate_threshold.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rate_threshold.gif\" width=400></a></div>"""%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)),cWB_conf.web_link,ifar_flag,dir_bkg,cWB_conf.web_link,ifar_flag,dir_bkg)
           for qq in range(len(cWB_conf.Cuts_list)):
              print >>ff,"""<div class="tabbertab"><h2>%s</h2><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rate_threshold.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rate_threshold.gif\" width=400></a></div>"""%(cWB_conf.Cuts_list[qq].replace("_cut",""),cWB_conf.web_link,ifar_flag,cWB_conf.Cuts_list[qq],cWB_conf.web_link,ifar_flag,cWB_conf.Cuts_list[qq])
           print >>ff,"</div>"
        print >>ff,"</td></tr>"
        print >>ff,"<tr align=\"justify\"><td><i>Cumulative rate of triggers according to their amplitude distribution.</i></td></tr>"
        print >>ff,"</table><p>"
        print >>ff,"<table border=0 cellspacing = 20 width=60% bgcolor=\"#F2F2F2\">"
        print >>ff,"<tr><th colspan=\"2\">Ranking statistics vs time and frequency</th></tr>"
        print >>ff,"<tr>"
        if (tabber==False):
           print >>ff,"<td align=\"left\"><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_time.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_time.gif\" width=400></a></td><td align=\"right\"><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_frequency.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_frequency.gif\" width=400></a></td>"%(cWB_conf.web_link,ifar_flag,dir_bkg,cWB_conf.web_link,ifar_flag,dir_bkg,cWB_conf.web_link,ifar_flag,dir_bkg,cWB_conf.web_link,ifar_flag,dir_bkg)
        else:
           print >>ff,"<td align=\"center\">"
           print >>ff,"<div class=\"tabber\">"
           print >>ff,"""<div class="tabbertab"><h2>%s</h2><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_time.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_time.gif\" width=400></a></div>"""%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)),cWB_conf.web_link,ifar_flag,dir_bkg,cWB_conf.web_link,ifar_flag,dir_bkg)
           for qq in range(len(cWB_conf.Cuts_list)):
              print >>ff,"""<div class="tabbertab"><h2>%s</h2><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_time.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_time.gif\" width=400></a></div>"""%(cWB_conf.Cuts_list[qq].replace("_cut",""),cWB_conf.web_link,ifar_flag,cWB_conf.Cuts_list[qq],cWB_conf.web_link,ifar_flag,cWB_conf.Cuts_list[qq])
           print >>ff,"</div>"
           print >>ff,"</td>"
           print >>ff,"<td align=\"center\">"
           print >>ff,"<div class=\"tabber\">"
           print >>ff,"""<div class="tabbertab"><h2>%s</h2><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_frequency.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_frequency.gif\" width=400></a></div>"""%("+".join(map(lambda x: x.replace("_cut",""),cWB_conf.Cuts_list)),cWB_conf.web_link,ifar_flag,dir_bkg,cWB_conf.web_link,ifar_flag,dir_bkg)
           for qq in range(len(cWB_conf.Cuts_list)):
              print >>ff,"""<div class="tabbertab"><h2>%s</h2><a href=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_frequency.gif\" rel=\"shadowbox[gallery]\"><img src=\"%s/POSTPRODUCTION/%s/plot%s/data/rho_frequency.gif\" width=400></a></div>"""%(cWB_conf.Cuts_list[qq].replace("_cut",""),cWB_conf.web_link,ifar_flag,cWB_conf.Cuts_list[qq],cWB_conf.web_link,ifar_flag,cWB_conf.Cuts_list[qq])
           print >>ff,"</div>"
           print >>ff,"</td>"
        print >>ff,"</tr>"
        print >>ff,"<tr align=\"justify\"><td colspan=\"2\"><i>Ranking statistic as a function of gps time (left) and estimated frequency (right). The rho parameter is a approximately equal to median SNR per detector. In green GW candidates triggers, in red triggers which have same significance as ones send to gracedb, in black all the others.</i></td></tr>"
        print >>ff,"</table><p>"
        print >>ff,"</center>"

def findced(cedstring):
    #print cedstring.split("_")
    trigger_time=float(cedstring.split("_")[len(cedstring.split("_"))-1].replace("\/",""))
    #print trigger_time
    files=glob.glob("%s/%s/??????/*/job.pickle"%(cWB_conf.run_dir,cWB_conf.jobs_dir))
    files=filter(lambda x: int(x.split("/")[-2].split("-")[0])-cWB_conf.seg_duration<=trigger_time and int(x.split("/")[-2].split("-")[1])>trigger_time, files)
    for file in files:
       f=open(file,"rb")
       j=pickle.load(f)
       f.close()
       if (j.back<trigger_time and j.end>trigger_time):
          if (len(j.triggers_all)>0):
             for t in j.triggers_all:
               if (t.topublish==True and abs(float(t.start[0])-trigger_time)<.1):
                 try:
                   ced=t.ced_link.split("\"")[1]
                 except:
                   ced=t.compute_ced().replace("OUTPUT","OUTPUT_CED")
    return ced

def substitute_ced(ifile):
    tfile=ifile.replace(".html","_temp.html")
    commands.getstatusoutput("mv %s %s"%(ifile,tfile))
    lines=open(tfile).readlines()
    ff=open(ifile,"w")
    for l in lines:
       if (l.find("href")!=-1 and l.find("ced")!=-1):
          cedstring=l.split("\"")[1]
          #print cedstring
          newced=findced(cedstring)
          #print "OLD: %s NEW: %s"%(cedstring,newced)
          print >>ff,l.replace(cedstring,newced)
       else:
          print >>ff,l
    ff.close()

def plot_FOMs(rootfiles,wavefiles,livefiles,rootdir,start,end):
    olddir=os.getcwd()
    namefile="cWB_online"
    version="M1"
    commands.getstatusoutput("rm -rf %s/%s;mkdir -p %s/%s"%(rootdir,namefile,rootdir,namefile))
    os.chdir("%s/%s"%(rootdir,namefile))
    commands.getstatusoutput("mkdir -p output config condor merge report/postprod")
    f=open("merge/merge_%s.%s.lst"%(namefile,version),"w")
    print >>f, "\n".join(rootfiles)
    f.close()
    f=open("merge/wave_%s.%s.lst"%(namefile,version),"w")
    print >>f, "\n".join(wavefiles)
    f.close()
    f=open("merge/live_%s.%s.lst"%(namefile,version),"w")
    print >>f, "\n".join(livefiles)
    f.close()
    for rootfile in rootfiles:
        commands.getstatusoutput("ln -sf %s output/"%(rootfile))
    commands.getstatusoutput("cp %s config/user_parameters.C"%(cWB_conf.zerolag_par))
    commands.getstatusoutput("cp %s config/user_pparameters.C"%(cWB_conf.pp_par))
    #commands.getstatusoutput("ln -sf %s/template.merged/postprod.csh"%(cWB_conf.run_dir))
    #commands.getstatusoutput("ln -sf %s/template.merged/setcuts.csh"%(cWB_conf.run_dir))
    #commands.getstatusoutput("ln -sf %s/template.merged/mergeTREE.C"%(cWB_conf.run_dir))    
    commands.getstatusoutput("echo merge/wave_%s.%s.lst wave_%s.%s waveburst | root -l -b -q %s/tools/online/bin/mergeTREE.C"%(namefile,version,namefile,version,os.environ['HOME_WAT']))
    commands.getstatusoutput("echo merge/wave_%s.%s.root %s | root -l -b -q %s/tools/online/bin/addHistory.C"%(namefile,version,rootfiles[0],os.environ['HOME_WAT']))
    commands.getstatusoutput("echo merge/live_%s.%s.lst live_%s.%s liveTime | root -l -b -q %s/tools/online/bin/mergeTREE.C"%(namefile,version,namefile,version,os.environ['HOME_WAT']))
    print commands.getstatusoutput("mkdir -p ../FOMs")
    clock_utc=tconvert(end)
    #com="./postprod.csh %s %i %i"%(version,0,0)
    com="%s/scripts/cwb_report.csh %s create %i %i"%(os.environ['HOME_CWB'],version,0,0)
    print com
    commands.getstatusoutput(com)
    pp_dir=glob.glob("report/postprod/*")
    print pp_dir
    print commands.getstatusoutput("rm -rf ../FOMs/plot")
    print commands.getstatusoutput("mv %s ../FOMs/plot"%(pp_dir[0]))
    substitute_ced("../FOMs/plot/body.html")
    try:
       if (len(cWB_conf.Cuts_list)>1):
          T_Cuts_list=["OR_cut"]
       else:
          T_Cuts_list=[]
       for n in cWB_conf.Cuts_list:
            T_Cuts_list.append(n)
       for n in T_Cuts_list:
         #com="""./setcuts.csh %s '--tcuts %s --label %s' """%(version,n,n)
         com="""%s/scripts/cwb_setcuts.csh %s '--tcuts %s --label %s' """%(os.environ['HOME_CWB'],version,n,n)
         #print com
         commands.getstatusoutput(com)
         if (os.path.exists("merge/live_%s.%s.C_%s.root"%(namefile,version,n))):
           #com="./postprod.csh %s.C_%s %i %i"%(version,n,0,0)
           com="%s/scripts/cwb_report.csh %s.C_%s create %i %i"%(os.environ['HOME_CWB'],version,n,0,0)
           #print com
           commands.getstatusoutput(com)
           pp_dir=glob.glob("report/postprod/*")
           print pp_dir
           print commands.getstatusoutput("rm -rf ../FOMs/plot%s"%n)
           print commands.getstatusoutput("mv %s ../FOMs/plot%s"%(pp_dir[0],n))
           substitute_ced("../FOMs/plot%s/body.html"%n)
         else:
           print "no event in this class"
    except:
       print "no class"
    print commands.getstatusoutput("rm -rf ../FOMs/merge")
    print commands.getstatusoutput("mv merge ../FOMs/.")
    os.chdir("../")
    print commands.getstatusoutput("rm -rf %s"%(namefile))
    os.chdir(olddir)

def run_statistics(stat_table,delay_launch_completion,delay_launch_start,delay_start_completion,rootdir,start,end,run_statistics_fn, qPrintJobs, title, output_file_name, job_segments, running_jobs, qPrintTriggers,jobs_str):
    fn_tmp="%s.tmp"%(run_statistics_fn)
    ff=open(fn_tmp,"w")
    print >>ff,"<html>"
    print >>ff,"<title>%s</title>"%(title)
    print >>ff,"<body>"

    print >>ff, stat_table
    print >>ff,"<p>"
    print >>ff,"<center>"
    print >>ff,"<table>"
    print >>ff,"<tr align=\"center\">"
    print >>ff,"<td><img src=\"running_time.png\"></td>"
    print >>ff,"<td><img src=\"launch_time.png\"></td>"
    print >>ff,"<td><img src=\"completion_time.png\"></td>"
    print >>ff,"</tr>"
    print >>ff,"<tr>"
    print >>ff,"<td><i>Time analysis for each segment.</i></td>"
    print >>ff,"<td><i>Delay time between the starting of the analysis<br>and the end of the segment</i></td>"
    print >>ff,"<td><i>Difference between the end of the analysis job<br>and the end segment time</i></td>"
    print >>ff,"</tr>"
    print >>ff,"</table>"
    print >>ff,"</center>"
    
    pylab.figure(figsize=(4,4))
    try:
        b=pylab.hist(delay_launch_completion)
    except:
        pass
    pylab.grid(True)
    pylab.xlabel("job running time, seconds")
    pylab.ylabel("events")
    pylab.title("Job running time")
    pylab.savefig("%s/running_time.png"%(rootdir))

    pylab.figure(figsize=(4,4))
    try:
        b=pylab.hist(delay_launch_start)
    except:
        pass
    pylab.grid(True)
    pylab.xlabel("delay launch time, seconds")
    pylab.ylabel("events")
    pylab.title("Delay launch time")
    pylab.savefig("%s/launch_time.png"%(rootdir))

    pylab.figure(figsize=(4,4))
    try:
        b=pylab.hist(delay_start_completion)
    except:
        pass
    pylab.grid(True)
    pylab.xlabel("job completion time, seconds")
    pylab.ylabel("events")
    pylab.title("Job completion time")
    pylab.savefig("%s/completion_time.png"%(rootdir))

    current_time=get_time()

    #f = open(cWB_conf.segs1)
    f = open("%s/%s/seg_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir))
    segments = fromsegwizard(f)
    f.close()

    segments &= segmentlist([segment(start,end)]); segments.coalesce()

    segments_fnr="science_segments.txt"
    segments_fn="%s/%s"%(rootdir,segments_fnr)
    fff=open(segments_fn,"w")
    tosegwizard(fff,segments)
    fff.close()

    job_segments_fnr="job_segments.txt"
    job_segments_fn="%s/%s"%(rootdir,job_segments_fnr)
    fff=open(job_segments_fn,"w")
    tosegwizard(fff,job_segments)
    fff.close()

    diff_segments_job_segments = segments - job_segments; diff_segments_job_segments.coalesce()

    diff_segments_job_segments_fnr="diff_segments_job_segments.txt"
    diff_segments_job_segments_fn="%s/%s"%(rootdir,diff_segments_job_segments_fnr)
    fff=open(diff_segments_job_segments_fn,"w")
    tosegwizard(fff,diff_segments_job_segments)
    fff.close()
    
    print >>ff,"<h2 align=center>Job coverage of segments</h2>"
    print >>ff,"<table border=1>"
    print >>ff,"<tr><th align=left>Current time</th><td>%d :  %s</td></tr>"%(current_time,tconvert(current_time))
    print >>ff,"<tr><th align=left>Start of the considered interval</th><td>%d :  %s</td></tr>"%(start,tconvert(start))
    print >>ff,"<tr><th align=left>End of the considered interval</th><td>%d :  %s</td></tr>"%(end,tconvert(end))    
    print >>ff,"<tr><th align=left>The end time of the last known Science segment</th><td>%d :  %s</td></tr>"%(segments[-1][-1],tconvert(segments[-1][-1]))
    print >>ff,"<tr><th align=left>Livetime in the <a href=\"%s\">segments</a> up to the current time</th><td>%d</td></tr>"%(segments_fnr,abs(segments))
    print >>ff,"<tr><th align=left>Livetime in <a href=\"%s\">segments</a> of completed or running jobs</th><td>%d</td></tr>"%(job_segments_fnr,abs(job_segments))
    print >>ff,"<tr><th align=left>Duration of <a href=\"%s\">difference segments</a> between the above two lists</th><td>%s</td>"%(diff_segments_job_segments_fnr, map(lambda x: abs(x), diff_segments_job_segments))
    print >>ff,"</table>"

    note="""
<p>
Each cWB job takes 180 seconds of input data plus 8 second of offset at the beginning and end.
Each job overlaps the previous one of 120 seconds (excluding the first one), so the jobs time shift is 60 seconds.
This means that at least 60 seconds of unprocessed data are required to launch a new job.
Any job which length is less than 180 seconds is not considered at the moment.
<p>

"""
#In the future we are going to use a more complicated scheduling algorithm: job segments will be allowed to intersect and therefore less data should be wasted and there will be smaller chance that a signal is lost on the boundaries between adjacent job segments.
#If jobs are correctly covering Science segments on the available frames there should not be any segments with the duration greater than 76 seconds for the second from the bottom row in the table.

    print >>ff,note

    if(qPrintJobs):
        
        print >>ff,"<h2 align=center>Completed jobs</h2>"
        print >>ff,"<table border=1>"
        print >>ff,"<tr><th>start date/time</th><th>start</th><th>end</th><th>launch</th><th>completion</th><th>completion-launch</th><th>launch-start</th><th>number of triggers</th><th>link</th></tr>"
        for job_str in jobs_str:
           print  >>ff,job_str
        #for job in JOBS:
        #    jdir=cWB_conf.web_link+"/"+"/".join(job.dir.split("/")[-3:])
        #    p_triggers=filter(lambda x: x.topublish==True, job.triggers_all)
        #    if (len(job.triggers_all)-len(p_triggers) > 0):
        #       s_tr="    [%s]"%(len(job.triggers_all)-len(p_triggers))
        #    else:
        #       s_tr=""
        #    if(job.status==2):
        #        print >>ff,"<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d%s</td><td><a href=\"%s/job.html\">link</a></td></tr>"%\
        #              (tconvert(job.start),job.start,job.end,job.launch_time,job.completion_time,job.completion_time - job.launch_time,job.launch_time-job.start, len(p_triggers),s_tr,jdir)
#                      (tconvert(job.start),job.start,job.end,job.launch_time,job.completion_time,job.completion_time - job.launch_time,job.launch_time-job.start, len(job.triggers_all),jdir)
        print >>ff,"</table>"

        print >>ff,"<h2 align=center>Running jobs</h2>"
        print >>ff,"<table border=1>"
        print >>ff,"<tr><th>start date/time</th><th>start</th><th>end</th><th>launch</th><th>completion</th><th>completion-launch</th><th>launch-start</th><th>number of triggers</th><th>link</th></tr>"
        for job in running_jobs:
            jdir=cWB_conf.web_link+"/"+"/".join(job.dir.split("/")[-3:])            
            print >>ff,"<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>...</td><td>...</td><td>%d</td><td>...</td><td><a href=\"%s/log\">log</a></td></tr>"%\
                  (tconvert(job.start),job.start,job.end,job.launch_time,job.launch_time-job.start,\
                   jdir)
        print >>ff,"</table>"


    ff.flush()
    ff.close()
    commands.getstatusoutput("mv %s %s"%(fn_tmp,run_statistics_fn))

def llcache(trigger):
    start,end=map(int,os.getcwd().split("/")[-1].split("-"))
    print "%s %s"%(start,end)
    start=trigger.job.back-cWB_conf.job_offset

    ms=start/4*4-8
    me=end/4*4
    if(me!=end):
        me+=4
    me+=8

    for ifo,dir,fnames in zip(cWB_conf.ifos,cWB_conf.bkg_dir,cWB_conf.bkg_fnames):
        f=open("input/%s_scratch.frames"%ifo,"w")
        f1=open("input/%s_scratch.missed"%ifo,"w")
        for t in range(ms,me+4,4):
            five=str(t)[:5]
            p="%s/%s/%s%s/%s%s-4.gwf"%(dir,ifo,fnames,five,fnames,t)
            #print p
            if(os.path.exists(p)):
                print >>f,p
            else:
                print >>f1,p
                #sys.exit(1)
        f.close()
        f1.close()


def followup(trigger):
    print "In followup for trigger\n%s"%(repr(trigger))
    olddir=os.getcwd()
    os.chdir(trigger.job.dir)
    print "dir=%s"%(os.getcwd())
    try:
       pe_par=cWB_conf.pe_par
       do_pe=True
    except:
       do_pe=False
    try:
        print "Searching for frames"
        llcache(trigger)
        launch_ced(trigger,do_pe)
    except Exception,e:
        print "Not found frames for %f"%float(trigger.time[0])
    os.chdir(olddir)

def launch_ced(trigger,do_pe):
    #print "In launch_ced for trigger\n%s"%(repr(trigger))
    ced_dir="OUTPUT_CED/%s"%(trigger.compute_ced_dir())
    tmp_ced_dir="tmp_ced/%s"%(trigger.compute_ced_dir())
    #print "ced_dir in launch_ced = %s"%(ced_dir)
    #print "tmp_dir in launch_ced = %s"%(tmp_ced_dir)
    if((not os.path.exists("%s/start_ced"%(trigger.job.dir))) and (not os.path.exists(ced_dir)) and (not os.path.exists(tmp_ced_dir))):
        com="%s/tools/online/RUN_cWB/bin/followup.csh %.3f %s %s true"%(os.environ['HOME_WAT'], float(trigger.start[0]),"ced",cWB_conf.ced_par)
        #print com
        a=commands.getstatusoutput(com)
        #print a
        if (do_pe==True):
          com="%s/tools/online/RUN_cWB/bin/followup.csh %.3f %s %s false"%(os.environ['HOME_WAT'], float(trigger.start[0]),"pe",cWB_conf.pe_par)
          #print com
          a=commands.getstatusoutput(com)
          #print a
        try:
           num_off=len(cWB_conf.science_segment_offset)
           com="ln -s %s/OUTPUT_CED/%s %s/OUTPUT_CED/."%(trigger.job.dir,ced_dir.split("/")[1],"/".join(trigger.dir.split("/")[:-2]))
           #print com
           a=commands.getstatusoutput(com)
           #print a
        except:
           num_off=0

def trigger_table(trigger_table_fn, triggers, title, injections, firstTime=False):
    ff=open(trigger_table_fn+".tmp","w")
    print >>ff, "<html>"
    print >>ff, "<title>%s</title>"%(title)
    print >>ff, "<body>"
    print >>ff, "<h1 align=center>%s</h1>"%(title)
    print >>ff,"<table border=1>"
    header=triggers[0].table_header(cWB_conf.ifos,True)
    print >>ff, header
    counter=1

    for t in triggers:
        jdir=cWB_conf.web_link+"/"+"/".join(t.job_dir.split("/")[-3:])
        t.cat_passed_txt=""
        #try:
        #    #if(t.selected>=3 or t.cat_passed_html!=""):
        #    #    followup(t)
        #except Exception,e:
        #    print e
        #for ifo in cWB_conf.ifos:
        #    t.cat_passed_html+="<tr><th>%s</th>"%(ifo)
        #    for c in range(5):
        #        try:
        #            t.cat_passed_txt+="CAT%d_%s = %d\n"%(c+1,ifo,int(t.cat_passed[c][ifo]))
        #        except:
        #            t.cat_passed_txt+="CAT%d_%s = %d\n"%(c+1,ifo,1)
        #        try:
        #            t.cat_passed_html+="<td>%d</td>"%(int(t.cat_passed[c][ifo]))
        #        except:
        #            t.cat_passed_html+="<td>%d</td>"%(1)
        #    t.cat_passed_html+="</tr>\n"
        #t.cat_passed_html+="</table>\n"
        #if(firstTime):
        #    try:
        #        indx=injections.search(float(t.time[0]),cWB_conf.injections_match_window)
        #        if(indx>=0):
        #            t.comments+="<br>Injection %s"%(injections.injections[indx].GravEn_SimID)
        #            t.comments+="<br>(t_d - t_i) = %.3f"%(float(t.time[0]) - injections.injections[indx].sort_time)
        #        else:
        #            t.comments+="<br><blink><b><font color=\"red\">Candidate trigger!</font></b></blink>"
        #    except Exception,e:
        #        print "error while assigning injection dq to a trigger"
        #        print e
        print >>ff, t.to_html_table_row(counter)                                                                                                                            
        counter+=1
    print >>ff, "</table>"
    print >>ff, "</body>"
    print >>ff, "</html>"
    ff.flush()
    ff.close()
    commands.getstatusoutput("mv %s.tmp %s"%(trigger_table_fn, trigger_table_fn))
    
if(__name__=="__main__"):
  jobs=glob.glob("%s/%s/??????/*/OUTPUT/finished"%(cWB_conf.run_dir,cWB_conf.jobs_dir))
  jobs.sort()
  print len(jobs)
  jobs=map(lambda x: x.split("/")[len(x.split("/"))-3],jobs)
  starts=map(lambda x: int(x.split("-")[0]),jobs)
  stops=map(lambda x: int(x.split("-")[1]),jobs)
  start=starts[0]
  end=stops[-1]
  page_rootdir="%s/%s/the_whole_run"%(cWB_conf.run_dir,cWB_conf.summaries_dir)
  page_rootlink="%s/%s/the_whole_run"%(cWB_conf.web_link,cWB_conf.summaries_dir)
  summary_page(start,end,cWB_conf.web_dir,page_rootdir,True,True,"the_whole_run2.html","The whole run","FOM_10000_0_0")

