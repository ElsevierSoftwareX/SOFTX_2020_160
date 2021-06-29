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
import summary_page
#import injections

def get_start_of_the_day_gps(ctime):
    a=commands.getstatusoutput("lalapps_tconvert -f \"%D %T\" "+repr(ctime))
    t=a[1].split()
    timezero="00:00:00"
    a=commands.getstatusoutput("lalapps_tconvert "+t[0]+" "+timezero)
    day=int(a[1].strip())
    return day

try:
  choice=sys.argv[1]
except Exception,e:
  choice="all"

#f = open("%s/%s/seg_global.txt"%(cWB_conf.online_dir,cWB_conf.seg_dir))
#segments = fromsegwizard(f)
#f.close()
jobs=glob.glob("%s/%s/??????/*/OUTPUT/finished"%(cWB_conf.run_dir,cWB_conf.jobs_dir))
jobs.sort()
print len(jobs)
jobs=map(lambda x: x.split("/")[len(x.split("/"))-3],jobs)
starts=map(lambda x: int(x.split("-")[0]),jobs)
stops=map(lambda x: int(x.split("-")[1]),jobs)

#print "livetime=%d"%(abs(segments))
current_time=get_time()

if (choice.find("all")!=-1 or choice.find("daily")!=-1):
 print "-"*50
 print "Each day"
 rootdir="%s/%s/daily"%(cWB_conf.run_dir,cWB_conf.summaries_dir)
 rootlink="%s/%s/daily"%(cWB_conf.web_link,cWB_conf.summaries_dir)
 cWB_conf.run_offset=int(open("offset.txt").readlines()[0])
 try:
  #start=segments[0][0]
  #end=segments[-1][1]
  start=starts[0]
  end=current_time#stops[-1]
  if (cWB_conf.run_offset>0):
    end=stops[-1]

  day=get_start_of_the_day_gps(end-1000)
  print "day=%d"%(day)

  while( segmentlist([segment(day,day+24*3600)]).intersects(segmentlist([segment(start,end)])) ):
    page_rootdir="%s/%d-%d"%(rootdir,day,day+3600*24)
    page_rootlink="%s/%d-%d"%(rootlink,day,day+3600*24)    
    print "page_rootdir=%s"%page_rootdir
    commands.getstatusoutput("mkdir -p %s"%(page_rootdir))
    a=commands.getstatusoutput("lalapps_tconvert -f \"%D\" "+repr(day))
    cdate=a[1].strip()

    summary_page.summary_page(day,day+3600*24+1,cWB_conf.web_dir,page_rootdir,True,True,"day.html","%s"%(cdate))    

    day=get_start_of_the_day_gps(day-3600*24+1000)
    print "="*80

 except:
  print "no jobs"
    
 commands.getstatusoutput("./run_calendar.py")

#----------------------------------------------------------
#if (choice.find("all")!=-1 or choice.find("hour")!=-1):
# print "-"*50
# print "Last hour"
# end=current_time
# start=end-3600
# page_rootdir="%s/%s/last_hour"%(cWB_conf.run_dir,cWB_conf.summaries_dir)
# page_rootlink="%s/%s/last_hour"%(cWB_conf.web_link,cWB_conf.summaries_dir)
# summary_page.summary_page(start,end,cWB_conf.web_dir,page_rootdir,True,True,"last_hour.html","Last hour","FOM_1_0_0")

#----------------------------------------------------------
#if (choice.find("all")!=-1 or choice.find("mid")!=-1):
# print "-"*50
# print "Last 12 hours"
# end=current_time
# start=end-3600*12
# page_rootdir="%s/%s/last_12_hours"%(cWB_conf.run_dir,cWB_conf.summaries_dir)
# page_rootlink="%s/%s/last_12_hours"%(cWB_conf.web_link,cWB_conf.summaries_dir)
# summary_page.summary_page(start,end,cWB_conf.web_dir,page_rootdir,True,True,"last_12_hours.html","Last 12 hours","FOM_12_0_0")

#----------------------------------------------------------
if (choice.find("all")!=-1 or choice.find("day")!=-1):
 print "-"*50
 print "Last day"
 end=current_time
 start=end-3600*24
 page_rootdir="%s/%s/last_day"%(cWB_conf.run_dir,cWB_conf.summaries_dir)
 page_rootlink="%s/%s/last_day"%(cWB_conf.web_link,cWB_conf.summaries_dir)
 summary_page.summary_page(start,end,cWB_conf.web_dir,page_rootdir,True,True,"last_day.html","Last day","FOM_24_0_0")

#----------------------------------------------------------
if (choice.find("all")!=-1 or choice.find("week")!=-1):
 print "-"*50
 print "Last week"
 end=current_time
 start=end-3600*24*7
 page_rootdir="%s/%s/last_week"%(cWB_conf.run_dir,cWB_conf.summaries_dir)
 page_rootlink="%s/%s/last_week"%(cWB_conf.web_link,cWB_conf.summaries_dir)
 summary_page.summary_page(start,end,cWB_conf.web_dir,page_rootdir,True,True,"last_week.html","Last week","FOM_168_0_0")

#----------------------------------------------------------
if (choice.find("all")!=-1 or choice.find("run")!=-1):
 print "-"*50
 print "All run"
 try:
  #start=segments[0][0]
  #end=segments[-1][1]
  start=starts[0]
  end=stops[-1]
 except:
  start=0
  end=current_time
 page_rootdir="%s/%s/the_whole_run"%(cWB_conf.run_dir,cWB_conf.summaries_dir)
 page_rootlink="%s/%s/the_whole_run"%(cWB_conf.web_link,cWB_conf.summaries_dir)
 summary_page.summary_page(start,end,cWB_conf.web_dir,page_rootdir,True,True,"the_whole_run.html","The whole run","FOM_10000_0_0")

#----------------------------------------------------------
if (choice.find("all")!=-1 or choice.find("check")!=-1):
  commands.getstatusoutput("./checkstatus.py")
