#!/usr/bin/env python

# Copyright (C) 2019 Marco Drago, Andrea Miani
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

import commands, os, sys

import cWB_conf
from run import *
import Trigger
import run_utils

def send_mails(line):
    msg="""Subject:  Check online status at %s %s

%s"""%(os.environ['SITE_CLUSTER'],cWB_conf.online_dir,line)
    #print msg
    com="echo \"%s\" | /usr/sbin/sendmail -F %sonline-%s %s"%(msg,cWB_conf.gracedb_analysis,cWB_conf.gracedb_search,",".join(cWB_conf.error_emails))
    commands.getstatusoutput(com)
    print "Ok - Controllo def send_mails"
def return_error(value,time):
  timelimit=2*cWB_conf.sleep
  machine=os.environ["HOSTNAME"]
  if (value>0):
    return "OK","RUNNING AT %s"%machine
  elif (time<timelimit):
    return "WARNING","NOT RUNNING AT %s, LAST UPDATE: %i s"%(machine,time)
  else:
    return "ERROR","NOT RUNNING AT %s, LAST UPDATE: %i s"%(machine,time)
  print "Ok - Controllo def return_error"  
def def_colors(w):
    if (w=="ERROR"):
      return "red"
    if (w=="WARNING"):
      return "yellow"
    if (w=="OK"):
      return "limegreen"
def findnotfin(started,pos_s,finished,pos_f):
  started.sort()
  started_time=map(lambda x: x.split("/")[pos_s],started)
  st=segmentlist([])
  for s in started_time:
    st.append(segment(s.split("-")[0],s.split("-")[1]))
  st.coalesce()
  finished.sort()
  finished_time=map(lambda x: x.split("/")[pos_f],finished)
  fin=segmentlist([])
  for f in finished_time:
    fin.append(segment(f.split("-")[0],f.split("-")[1]))
  fin.coalesce()
  st-=fin
  #print repr(st)
  logs=[]
  for s in st:
    start=float(s[0])
    end=start+cWB_conf.moving_step
    while(end<=float(s[1])):
      dirs=filter(lambda x: x.find("%i-%i"%(start,end))!=-1,started)
      if (len(dirs)>0):
        #print repr(dirs)
        logs.append(dirs[0])
        start+=cWB_conf.moving_step
        end+=cWB_conf.moving_step
      else:
        end+=cWB_conf.moving_step
  return logs
  print "Ok - Controllo del findnotfin"

#---------- CLASS: TIME ----------#

class Time:
  def __init__(self,n,t=0):
    self.value=t
    self.name=n
  def __sub__(self,otherTime):
    return Time("(%s-%s)"%(self.name,otherTime.name),self.value-otherTime.value)
  def table_line(self):
        if self.value != 0:
            return "<tr><th align=left>%s</th><td>%d : %s</td></tr>"%(self.name,self.value,tconvert(self.value))
        else:
            return "<tr><th align=left>%s</th><td>%s</td></tr>"%(self.name,"detector is not operative")
        print "OK - Controllo class Time - def table line"    
  def last_analyzed_time(self,path):
    updirs=glob.glob("%s"%path)
    if updirs:
        downdirs=glob.glob("%s/*"%updirs[len(updirs)-1])
        last_dir=downdirs[len(downdirs)-1]
        last_analyzed_time=int(last_dir.split("/")[len(last_dir.split("/"))-1].split("-")[1])
    else:
        last_analyzed_time=0
    self.value=last_analyzed_time
    print "OK - Controllo class Time - def lat"
  def last_coincident_time(self,ifile):
    f=open(ifile)
    combined_global=fromsegwizard(f)
    f.close()
    if combined_global:
        if ((int(combined_global[-1][1])-int(combined_global[-1][0]))>(cWB_conf.seg_duration + 2*cWB_conf.job_offset)): # if instead of while
            last_coincident_time=int(combined_global[-1][1])
            print "lct = %i" %last_coincident_time
            print "Controllo - LCT - 1"
        else:
            #####
            i=1
            while ((int(combined_global[-1-i][1])-int(combined_global[-1-i][0]))<(cWB_conf.seg_duration + 2*cWB_conf.job_offset)):
                i+=1
                #####
                last_coincident_time=int(combined_global[-1-i][1])
                print "Controllo - LCT - 2"
            else:
                last_coincident_time=int(combined_global[-1-i][1])
                print "Controllo - LCT - 3"
    else:
        last_coincident_time=0
        print "Controllo - LCT - 4"
    self.value=last_coincident_time
    print "OK - Controllo class Time - def lct"
  def getlastfromfile(self,idir,ifile):
    a=commands.getstatusoutput("grep  %s %s | grep gwf | tail -n 1"%(idir,ifile))
    read_file=a[1].split(" ")[2]
    if (read_file.find(idir)==-1):
     print "Not found %s"%idir
    else:
     f=FRAME(read_file)
     self.value=f.end
  def getlastfromdir(self,idir):
    files=glob.glob("%s/*.gwf"%idir)
    files.sort()
    last_file=files[len(files)-1]
    f=FRAME(last_file)
    self.value=f.end
  print "Ok - Controllo class Time"

#---------- CLASS: CONDITION ----------#

class Condition:
  def __init__(self,Time,Value,w,line):
    bad="%s is more than %s (%i>%i), %s"%(Time.name,Value.name,Time.value,Value.value,line)
    good="%s is less than %s (%i<%i)"%(Time.name,Value.name,Time.value,Value.value)
    if (Time.value > Value.value):
      self.table_line="<tr><th colspan=\"2\" bgcolor=\"%s\">%s: %s</th></tr>"%(def_colors(w),w,bad)
      if (def_colors(w).find("red")!=-1):
        self.tmail=bad
    else:
      self.table_line="<tr><th colspan=\"2\" bgcolor=\"%s\">%s: %s</th></tr>"%(def_colors("OK"),"OK",good)
    print "Ok - Controllo class Condition"

#---------- CLASS: CHECK ----------#

class check:
  def __init__(self):
    self.times=[]
    self.conditions=[]
  def print_times(self):
    table="<center><table border=1 cellpadding=5>"
    for t in self.times:
      table+=t.table_line()
    for c in self.conditions:
      table+=c.table_line
    table+="</table></center>"
    return table
  def print_cond(self):
    res=""
    for c in self.conditions:
      try:
       res="%s\n%s"%(res,c.tmail)
      except:
       pass
    return res
  print "Ok - Controllo class check"

#-------------------------------------#
#---------- STARTING SCRIPT ----------#
#-------------------------------------#

if(__name__=="__main__"):
 current_time=Time("Current time",get_time_nooff())

 last_send_file="%s/log/sentstatus"%cWB_conf.online_dir
 send_mail=True
 text_mail=""
 if (os.path.exists(last_send_file)):
  a=commands.getstatusoutput("date -r %s | xargs tconvert"%(last_send_file))
  #print a[1]
  last_send=float(a[1].split("\n")[len(a[1].split("\n"))-1])
  #print current_time.value-last_send
  if ((current_time.value-last_send)<3600):
   send_mail=False
 print "Controllo - Starting Script"

#---------- Here I check for FRAMES ----------#

# Be avare that if some frames are missing there could be problems.
 frames_C=check()
 frames_C.times.append(current_time)
 frames_time=Time("Delay Frame Time",30)
 for i in range(0,len(cWB_conf.ifos)):
   last=Time("Last Frame %s"%cWB_conf.ifos[i])
   last.getlastfromdir(cWB_conf.frames_dir[i])
   frames_C.times.append(last)
   frames_C.conditions.append(Condition(current_time-last,frames_time,"WARNING","check if frames are updated"))
 print "Controllo Frames"

#---------- Here I check for TIMES ----------#

 current_time=Time("Current time",get_time_nooff())
 read_C=check()
 read_C.times.append(current_time)
 read_time=Time("Delay Frame Time",30)
 for i in range(0,len(cWB_conf.ifos)):
   last=Time("Last Read %s"%cWB_conf.ifos[i])
   try:
      last.getlastfromfile(cWB_conf.frames_dir[i],"%s/log/run.log"%(cWB_conf.production_dir))
   except:
     last.getlastfromfile(cWB_conf.frames_dir[i],"%s/log/run.log"%(cWB_conf.online_dir))
   read_C.times.append(last)
   read_C.conditions.append(Condition(current_time-last,read_time,"ERROR","check if pipeline is reading frames"))
 print "Controllo Times"

#---------- Here I check for ZERO LAG ANALYSIS ----------# 
 
 page_checkdir="%s/%s/check.html"%(cWB_conf.run_dir,cWB_conf.summaries_dir)
 current_time=Time("Current time",get_time_nooff())
 zero_l_an=Time("Last Analyzed Time")
 zero_l_an.last_analyzed_time("%s/%s/??????"%(cWB_conf.run_dir,cWB_conf.jobs_dir))
 zero_l_coinc=Time("Last Coincident Time")
 zero_l_coinc.last_coincident_time("%s/%s/seg_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir))
 seg_dur=Time("Segment length",cWB_conf.seg_duration+2*cWB_conf.job_offset)
 mov_step=Time("Moving step",cWB_conf.moving_step+cWB_conf.job_offset)

 print "Controllo 1 ZeroLagAnalysis"

 zero_C=check()
 zero_C.times.append(current_time)
 zero_C.times.append(zero_l_an)
 zero_C.times.append(zero_l_coinc)
 zero_C.conditions.append(Condition(zero_l_coinc-zero_l_an,seg_dur,"ERROR","check if pipeline is running"))
 zero_C.conditions.append(Condition(zero_l_coinc-zero_l_an,mov_step,"WARNING","check if data are good in the last period"))

 print "Controllo 2 ZeroLagAnalysis"

#---------- Here I check for BGK ANALYSIS ----------#

 current_time=Time("Current time",get_time_nooff())
 bkg_l_an=Time("Last Analyzed Time")
 bkg_l_an.last_analyzed_time("%s/%s/?????"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir))
 bkg_l_coinc=Time("Last Coincident Time")
 bkg_l_coinc.last_coincident_time("%s/%s/considered.txt"%(cWB_conf.bkg_run_dir,cWB_conf.seg_dir))
 bkg_seg_dur=Time("Segment length",cWB_conf.bkg_job_duration+2*cWB_conf.job_offset)

 bkg_C=check()
 bkg_C.times.append(current_time)
 bkg_C.times.append(bkg_l_an)
 bkg_C.times.append(bkg_l_coinc)
 bkg_C.conditions.append(Condition(bkg_l_coinc-bkg_l_an,bkg_seg_dur,"ERROR","check if jobs are submitted"))

 zero_notfin=glob.glob("%s/%s/??????/*/log"%(cWB_conf.run_dir,cWB_conf.jobs_dir))
 zero_st=glob.glob("%s/%s/??????/*/log.zip"%(cWB_conf.run_dir,cWB_conf.jobs_dir))
 zero_st+=zero_notfin
 zero_fin=glob.glob("%s/%s/??????/*/OUTPUT/finished"%(cWB_conf.run_dir,cWB_conf.jobs_dir))
 #zero_notfin=findnotfin(zero_st,-2,zero_fin,-3)
 bkg_st=[]
 bkg_st.append(glob.glob("%s/%s/?????/*/start"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir)))
 bkg_fin=[]
 bkg_fin.append(glob.glob("%s/%s/?????/*/finished"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir)))
 for l in range(1,cWB_conf.bkg_njobs):
  bkg_st.append(glob.glob("%s/%s_%i/?????/*/start"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir,l)))
  bkg_fin.append(glob.glob("%s/%s_%i/?????/*/finished"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir,l)))

 #logs=["run","web_pages_run","web_pages_week","web_pages_day","web_pages_mid","web_pages_hour","run_ts"]
 if hasattr(cWB_conf, 'production_dir'):
   logs=["run","web_pages_week","web_pages_day","web_pages_daily"]
 else:
   logs=["run","web_pages_week","web_pages_day","web_pages_daily","run_ts"]
 #logs=["run","web_pages_all","run_ts"]
 lens=[]
 gpss=[]
 for l in logs:
  a=commands.getstatusoutput("fuser %s/log/%s.log"%(cWB_conf.online_dir,l))
  lens.append(len(a[1]))
  a=commands.getstatusoutput("date -r %s/log/%s.log | xargs tconvert"%(cWB_conf.online_dir,l))
  gpss.append(get_time_nooff()-float(a[1].split("\n")[len(a[1].split("\n"))-1]))
  #print "%s %i"%(l,len(a[1]))

#---------- CREATION OF THE WEB PAGE ----------#

 current_time=Time("Current time",get_time_nooff())
 ff=open("%s"%(page_checkdir),"w")

#---------- Header of the Page ----------#

 print >>ff,"<html>"
 print >>ff,"<title>STATUS</title>"
 print >>ff,"<body>"
 print >>ff,"""<h1 align=center><font size=10 face="Courier" color="#DF013A">STATUS</font></h1>"""
 print >>ff,"<center><i>The page was generated by %s at %d : %s</i></center><p>"%(os.environ["HOSTNAME"],current_time.value,tconvert(current_time.value))
 print >>ff,"<hr size=3 noshade color=\"blue\">"

#---------- TABLE: Where and if scripts are running ----------#

 print >>ff, "<h3 align=center><font size=6 face=\"Courier\"  color=\"#01DF01\">Running scripts</font></h3>"
 print >>ff,"<center><table border=1 cellpadding=5>"
 for l,ll,tt in zip(logs,lens,gpss):
  s,ss=return_error(ll,tt)
  print >>ff,"<tr><td>%s</td><td bgcolor=\"%s\">%s</td></tr>"%(l,def_colors(s),ss)
  if (def_colors(s).find("red")!=-1):
    text_mail="%s\n%s %s"%(text_mail,l,ss)
 print >>ff,"</table></center>"

#---------- TABLE: Are the Frames Update? ----------#
 
 print >>ff,"<center><table border=0 cellpadding=5>"
 print >>ff,"<tr><td>"

 print >>ff, "<h3 align=center><font size=6 face=\"Courier\"  color=\"#01DF01\">Frames update</font></h3>"
 print >>ff,frames_C.print_times()
 text_mail="%s%s"%(text_mail,frames_C.print_cond())

#---------- TABLE: Are the Frames Read? ----------#

 print >>ff,"</td><td>"

 print >>ff, "<h3 align=center><font size=6 face=\"Courier\"  color=\"#01DF01\">Frames reading</font></h3>"
 print >>ff,read_C.print_times()
 text_mail="%s%s"%(text_mail,read_C.print_cond())

 print >>ff,"</td></tr>"
 print >>ff,"</table>"

#---------- TABLE: Zero Lag ----------#
 
 print >>ff,"<center><table border=0 cellpadding=5>"
 print >>ff,"<tr><td>"

 print >>ff, "<h3 align=center><font size=6 face=\"Courier\"  color=\"#01DF01\">Zero lag</font></h3>"
 print >>ff,zero_C.print_times()
 text_mail="%s%s"%(text_mail,zero_C.print_cond())

#---------- TABLE: Background ----------#
 
 print >>ff,"</td><td>"

 print >>ff, "<h3 align=center><font size=6 face=\"Courier\"  color=\"#01DF01\">Background</font></h3>"
 print >>ff,bkg_C.print_times()
 text_mail="%s%s"%(text_mail,bkg_C.print_cond())

 print >>ff,"</td></tr>"
 print >>ff,"</table>"

#---------- TABLE: Jobs Started & Finished ----------#

 print >>ff, "<h3 align=center><font size=6 face=\"Courier\"  color=\"#01DF01\">Started and finished</font></h3>"
 print >>ff,"<center><table border=1 cellpadding=5>"
 print >>ff,"<tr><td>Analysis</td><td>Started</td><td>Finished</td><td>Percentage</td><td>Log of not finished</td></tr>"
 notfin="<td>"
 lines=1
 if (len(zero_st)>0):
     for f in range(len(zero_notfin)):
       notfin+="<a href=\"%s\">%i</a>,"%(zero_notfin[f].replace(cWB_conf.run_dir,cWB_conf.web_link),f+1)
       if ((f+1)%25==0):
         notfin+="</td><tr><td>"
         lines+=1
     print >>ff,"<tr><td rowspan=\"%i\">Zero lag</td><td rowspan=\"%i\">%i</td><td rowspan=\"%i\">%i</td><td rowspan=\"%i\">%.1f</td>%s</td></tr>"%(lines,lines,len(zero_st),lines,len(zero_fin),lines,100.*len(zero_fin)/len(zero_st),notfin)
 for l in range(0,cWB_conf.bkg_njobs):
  if (len(bkg_st[l])>0):
    ratio=100.*len(bkg_fin[l])/len(bkg_st[l])
  else:
    ratio=0.
  print >>ff,"<tr><td>Superlag %i</td><td>%i</td><td>%i</td><td>%.1f</td></tr>"%(l,len(bkg_st[l]),len(bkg_fin[l]),ratio)
 print >>ff,"</table></center>"

 print >>ff,"</body>"
 print >>ff,"</html>"
 print text_mail
 if (send_mail==True and len(text_mail)>0):
  send_mails(text_mail)
  commands.getstatusoutput("rm %s;touch %s"%(last_send_file,last_send_file))
