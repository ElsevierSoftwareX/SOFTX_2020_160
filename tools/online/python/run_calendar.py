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
import commands, os, sys, glob, time
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
import calendar

def make_month(year,month,job_dir):
    print "year=%s month=%s"%(repr(year),repr(month))
    #print "In make_month job_dir=%s"%job_dir
    #a=commands.getstatusoutput("cal %d %d"%(month,year))
    #print a
    #b=map(lambda x:x.split(),a[1].split("\n"))
    #print b
    cal = calendar.TextCalendar(calendar.SUNDAY)
    a=cal.monthdayscalendar(year, month)
    b=[]
    namemonths=['January','February','March','April','May','June','July','August','September','October','November','December']
    monthname=[namemonths[month-1], '%i'%year]
    monthweek= ['Su', 'Mo', 'Tu', 'We', 'Th', 'Fr', 'Sa']
    b.append(monthname)
    b.append(monthweek)
    for i in a:
      u=[]
      for w in i:
         if (w!=0):
            u.append('%u'%w)
      b.append(u);
    print b
    title=" ".join(b[0])
    header="<tr><th>%s</th></tr>"%("</th><th>".join(b[1]))
    if(len(b[2])<7):
        c=['']*(7-len(b[2]))
        b[2]=c+b[2]
    if(len(b[-1])<7):
       c=['']*(7-len(b[-1]))
       b[-1]=b[-1]+c
    rows=[header]
    for i in range(2,len(b)):
        row=[]
        for d in b[i]:
            if(len(d)>0):
                g=commands.getstatusoutput("lalapps_tconvert %d/%s/%d"%(month,d,year))
                gps=g[-1].strip()
                h=glob.glob("%s/%s-*"%(job_dir,gps))
                print "In make_month: h=%s"%repr(h)
                if(len(h)==1):
                    mydir="/".join(h[0].split("/")[-3:])
                    row.append("<td><a href=\"%s/day.html\" target=\"Data\">%s</a></td>"%(mydir,d))
                elif(len(h)>1):
                    print "Something is wrong"
                else:
                    row.append("<td>%s</td>"%(d))
            else:
                row.append("<td></td>")
        srow="<tr>"+"\n".join(row)+"</tr>"
        rows.append(srow)
    table="<table border=1>\n%s\n</table>"%("\n".join(rows))
    return "<h4>%s</h4>\n%s"%(title,table)

def increment_month(year,month):
    if(month<=11):
        return [year,month+1]
    else:
        return [year+1,1]


try:
    cWB_conf.run_offset=int(open("offset.txt").readlines()[0])
except:
    pass

rpath="%s/daily"%cWB_conf.summaries_dir
rootdir="%s/%s"%(cWB_conf.web_dir,rpath)

dirs=glob.glob(rootdir+"/*")
dirs.sort()

start=int(dirs[0].split("/")[-1].split("-")[0])
end=int(dirs[-1].split("/")[-1].split("-")[1])

a=commands.getstatusoutput("lalapps_tconvert -f \"%Y\" "+repr(start))
start_year=int(a[1].strip())
a=commands.getstatusoutput("lalapps_tconvert -f \"%m\" "+repr(start))
start_month=int(a[1].strip())

a=commands.getstatusoutput("lalapps_tconvert -f \"%Y\" "+repr(end))
end_year=int(a[1].strip())
a=commands.getstatusoutput("lalapps_tconvert -f \"%m\" "+repr(end))
end_month=int(a[1].strip())

current_year=start_year
current_month=start_month


t=get_time()
t1=tconvert(t)

months=[]
months.append(make_month(current_year,current_month,rootdir))
if(start_year!=end_year or start_month!=end_month):
    [current_year,current_month]=increment_month(current_year,current_month)
    while(current_year!=end_year or current_month!=end_month):
        months.append(make_month(current_year,current_month,rootdir))
        [current_year,current_month]=increment_month(current_year,current_month)
#    months.append(make_month(end_month,end_year,rootdir))
    months.append(make_month(end_year,end_month,rootdir))
#try:
#    run_pid=open(cWB_conf.run_pid).readline().strip()
#except Exception,e:
#    print e
#    run_pid=""

#if(len(run_pid)>0):
#    com="ps -fp %s"%run_pid
#    a=commands.getstatusoutput(com)
#    out=a[1].split("\n")
#    print out
#    n=len(out)-1
#    try:
#        psout=out[1]
#    except:
#        n=0
#else:
#    n=0
#if(n==1):
#    uptime="%.2f days"%((time.time()-os.stat(cWB_conf.run_pid).st_ctime)/(3600*24.))
#    filter(lambda x:len(x)>0,psout.split(" "))[6]
#    status="<li><b><font color=green>uptime=%s<br>pid=%s</font></b></li>"%(uptime,run_pid)
#else:
#    status="<li><b><font color=red>Down</font></b></li>"

#a=commands.getstatusoutput("hostname -f")
#hostname=a[1].strip()

#f_f=open("%s/title.html"%cWB_conf.web_dir,"w")
#print >>f_f, """
#<html>
#<title>%s</title>
#<body>
#<h1 align=center>%s<h1>

#<!--<h2>Status</h2>
#<ul>
#%s
#<li><b>hostname</b><br>%s</li>
#<li><b>Page generated at</b><br>%d</br>%s</li>
#</ul>
#<hr>
#-->
#"""%(cWB_conf.title,cWB_conf.title,status, hostname, t, t1)

#print >>f_f, "</body>"
#print >>f_f, "</html>"

#f_f.flush()
#f_f.close()

ff=open("%s/calendar.html"%cWB_conf.web_dir,"w")
print >>ff, """
<html>
<body>

<h3>Calendar</h3>
"""

print >>ff,"\n".join(months)
print >>ff, "</body>"
print >>ff, "</html>"

ff.flush()
ff.close()

