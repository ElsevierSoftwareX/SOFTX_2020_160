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

import os, sys, glob
import math
if (os.environ['SITE_CLUSTER']=="CASCINA"):
  from glue.segments import *
  from glue.segmentsUtils import *
else:
  from ligo.segments import *
  from ligo.segments.utils import *
from run_utils import *

class Trigger:
    def __init__(self,lines,segs,job=None):
        self.job=job
        #self.lines=lines
        inmap=False
        skymap=[]
        counter=0
        notmaplines=[]
        for line in lines:
            if(not inmap):
                a=filter(lambda x: len(x)>0, line.split(" "))
                if(a[0].find("map_lenght")!=-1):
                    inmap=True
                    firsttime=True
                    maplength=int(a[1])
                else:
                    setattr(self,a[0].replace(":",""),a[1:])
                    notmaplines.append(line)
            else:
                if(firsttime):
                    firsttime=False
                    continue
                else:
                    skymap.append(filter(lambda x: len(x)>0, line.split(" ")))
                    counter=counter+1
                    if (counter==maplength):
                       inmap=False
        zeromap=skymap
        #self.error_region=map
        #try:
        #    qveto=cWB_conf.qveto_plugin
        #    #self.error_region=map[:-1]
        #    zeromap=map[:-1]
        #    Qveto=map[-1:][0]
        #    setattr(self,Qveto[0].replace(":",""),Qveto[1:])
        #except Exception,e:
        #    #self.error_region=map
        #    zeromap=map
        #self.error_region=filter(lambda x: float(x[7])>0, zeromap)
        self.DQ=[".."]
        self.comments=".."
        self.id=-1
        self.ced_link=""
        #self.qscan_links=[""]*len(cWB_conf.ifos)
        #self.qscan_link=""
        self.qscript=[","]*len(cWB_conf.ifos)*2
        self.seg=None
        #self.error_region_link=".."
        for s in segs:
            if(float(self.time[0])>=s[0] and float(self.time[0])<=s[1]):
                self.seg=s
                self.left_science="%.3f"%(float(self.time[0])-s[0])
                self.right_science="%.3f"%(s[1]-float(self.time[0]))
                break
        self.energy_disbalance="%.2e"%(float(self.netED[0]))
        if(self.job != None):
            self.dir="%s/OUTPUT.merged/TRIGGERS"%(self.job.dir)
            self.path="%s/OUTPUT.merged/TRIGGERS/trigger_%s.txt"%(self.job.dir, self.time[0])
            self.link="%s/OUTPUT.merged/TRIGGERS/trigger_%s.txt"%(self.job.link, self.time[0])
            f=open(self.path,"w")
            #f.write("\n".join(self.lines)+"\n")
            #f.write("\n".join(lines)+"\n")
            f.write("\n".join(notmaplines)+"\n")
            #print >>f,"wat version = %s\nonline version = %s\nsearch=%s"%(cWB_conf.version_wat, cWB_conf.version, cWB_conf.search)
            for ifo in cWB_conf.ifos:
              try:
                channel_reads="%s,%s"%(channel_reads,cWB_conf.channelname[ifo])
              except:
                channel_reads="%s"%(cWB_conf.channelname[ifo])
              try:
                t_lldq_flags="%s"%(",".join( map(lambda x,y: "%s:%s"%(x,y),cWB_conf.DQ_channel[ifo],cWB_conf.bitmask[ifo])))
              except:
                t_lldq_flags="%s"%("%s:%s"%(cWB_conf.DQ_channel[ifo],cWB_conf.bitmask[ifo]))
              try:
                lldq_flags="%s,%s"%(lldq_flags,t_lldq_flags)
              except:
                lldq_flags="%s"%(t_lldq_flags)
            print>>f, """
code version: %s
hoft version: %s,%s
lldq_flags: %s
"""%(cWB_conf.code_version,cWB_conf.hoft_version,channel_reads,lldq_flags)
            f.close()
        self.cutclass=""
        self.getclass()
        self.sgnf=[-1]*5
        self.significance()
        self.selected=self.select()
        self.send_time=-1
        self.cat_passed=[-1]*5
        for cat in range(0,5):
            self.cat_passed[cat]={}
        self.cat_passed_txt=""
        self.cat_passed_html=""
        self.inj_found=False
        self.getinjselect()
        self.ifar=[-1]*5
        self.creation_time=get_time()
        self.first_graceid=""
        self.graceid=""
        self.gracedb_link=""
        self.topublish=True
    #def dump_error_region(self,fn,link):
    #    self.error_region_link_txt="<a href=\"%s.txt\">txt</a>"%(link)
    #    f=open("%s.txt"%fn,"w")
    #    print >>f, "#skyID  theta   DEC     step   phi     R.A    step  probability    cumulative"
    #    error_region=add_counter_column(self.error_region)
    #    print >>f,"\n".join(map(lambda x:" ".join(x),error_region))
    #    f.close()
    #    self.error_region_link_html="<a href=\"%s.html\">html</a>"%(link)
    #    f=open("%s.html"%fn,"w")
    #    print >>f, totable(["counter","skyID","theta","DEC","step","phi","R.A","step","probability","cumulative"],error_region)
    #    f.close()
    def getclass(self):
        try:
           nclass=len(cWB_conf.Cuts_list)
           rootfile=glob.glob("%s/OUTPUT/*.root"%self.job.dir)[0]
           Cuts_file="%s/Cuts.hh"%(cWB_conf.config_dir)#,cWB_conf.Cuts_file.split("/")[len(cWB_conf.Cuts_file.split("/"))-1])
           for n in cWB_conf.Cuts_list:
              #print n
              com="echo %s %s %s %s | root -l -b -q %s/tools/online/bin/Cuts.C"%(rootfile,Cuts_file,n,self.time[0],os.environ['HOME_WAT'])
              if (os.environ['SITE_CLUSTER']=="CASCINA"):
                com="echo %s %s %s %s | root -l -b -q -n %s %s/tools/online/bin/Cuts.C"%(rootfile,Cuts_file,n,self.time[0],os.environ['CWB_ROOTLOGON_FILE'],os.environ['HOME_WAT'])
              #print com
              a=commands.getstatusoutput(com)
              result=a[1].split("\n")[len(a[1].split("\n"))-1]
              print "%s %s"%(com,result)
              if (int(result)==1):
                self.cutclass=n
        except:
           print "Only one class"
    def getinjselect(self):
        for inj in range(0,len(cWB_conf.inj_name)):
            inj_found=False
            for ifo in cWB_conf.ifos:
                fn_global="%s/%s/%s_%s_global.txt"%(cWB_conf.run_dir,cWB_conf.seg_dir,cWB_conf.inj_name[inj],ifo)
                print "fn_global=%s"%fn_global
                try:
                     f=open(fn_global)
                     inj_seg=fromsegwizard(f,float)
                     f.close()
                except Exception,e:
                     print e
                     inj_seg=segmenlist([])
                if (len(inj_seg)>0):
                     print inj_seg
                     for seg in inj_seg:
                         if (float(seg[0])<float(self.time[0]) and float(seg[1])>float(self.time[0])):
                            inj_found=True
            if (inj_found==True):
               if (self.inj_found==False):
                  self.cat_passed_html+="INJ %s"%(cWB_conf.inj_name[inj])
               else:
                  self.cat_passed_html+="<br>INJ %s"%(cWB_conf.inj_name[inj])
               self.inj_found=True
               print "INJ_%s"%(cWB_conf.inj_name[inj])
    def select(self): #0 - no extra thresholds, 1 - coherent, 2 - online publishing, 3 - external collaboration, 4 - offline threshold
        print map(type, [self.netcc[cWB_conf.id_cc], self.rho[cWB_conf.id_rho]]); sys.stdout.flush()
        try:
          nclass=len(cWB_conf.Cuts_list)
          coherent=len(self.cutclass)>0
        except:
          coherent=float(self.netcc[cWB_conf.id_cc])>cWB_conf.th_cc
        #print "%.2e %2e %s"%(float(self.sgnf[2].split(" ")[1]),cWB_conf.th_far_off,float(self.sgnf[2].split(" ")[1])<=cWB_conf.th_far_off)
        try:
          offline=float(self.rho[cWB_conf.id_rho])>=cWB_conf.th_rho_lum and float(self.sgnf[2].split(" ")[1])<=cWB_conf.th_far_off and coherent
        except:
          offline=float(self.rho[cWB_conf.id_rho])>=cWB_conf.th_rho_off and coherent
        lumin=float(self.rho[cWB_conf.id_rho])>=cWB_conf.th_rho_lum and coherent
        if(offline):
            return 4
        elif(lumin):
            return 3
        elif(coherent):
            return 1
        else:
            return 0
    def table_header(self,ifos,id=False):
        if(id):
            idh="<th>id</th>"
        else:
            idh=""
        #header="<tr>%s<th>ced<br>qscan</th><th>%s</th>"%(idh,"<br>".join(map(lambda x:"t"+x,ifos)))
        header="<tr>"
        if (id):
          header+="%s<th>Gracedb</th><th>ced</th>"%(idh)
        header+="<th>%s<br>(GPS, s)</th>"%("<br>".join(map(lambda x:"time"+x,ifos)))
        header+="<th><a href=\"%s/header_caption.html#rho\" target=_new>rho</a></th>"%(cWB_conf.web_link)
        header+="<th><a href=\"%s/header_caption.html#netCC\" target=_new>netCC</a></th>"%(cWB_conf.web_link)
        header+="<th><a href=\"%s/header_caption.html#ifar\" target=_new>IFAR</a><br>on last day (days)<br>latency (hours)</th><th>F (Hz)</th><th>%s</th>"%(cWB_conf.web_link, "<br>".join(map(lambda x:"snr"+x, ifos)))
        header+="<th><a href=\"%s/header_caption.html#netED\" target=_new>netED</a><br><a href=\"%s/header_caption.html#penalty\" target=_new>penalty</a></th>"%(cWB_conf.web_link, cWB_conf.web_link)
        header+="<th>%s</th><th>dT (s)<br>dF (Hz)</th></th>"%("<br>".join(map(lambda x:"hrss"+x, ifos)))
        header+="<th><a href=\"%s/header_caption.html#DEC\" target=_new>DEC</a> (&#176;)<br><a href=\"%s/header_caption.html#RA\" target=_new>R.A.</a> (&#176;)<br><a href=\"%s/header_caption.html#phi\" target=_new>phi</a> (&#176;)</th>"%(cWB_conf.web_link, cWB_conf.web_link, cWB_conf.web_link)
        #if (loudest==False):
        #  header+="<th><a href=\"%s/header_caption.html#error_region\">error<br>region</a></th>"%(cWB_conf.web_link)
        try:
          #th_qveto=cWB_conf.th_qveto
          #do_qveto=True
          nclass=len(cWB_conf.Cuts_list)
          do_class=True
        except:
          #do_qveto=False
          do_class=False
        if (do_class==True):
           header+="""<th><a href="%s/Cuts.hh.html" target=_new>Class</th>"""%(cWB_conf.web_link)
        #if (do_qveto==True):
        #  header+="<th>%s<br>Strain-Whitened</th>"%("<br>".join(map(lambda x:"Qveto"+x, ifos)))
        #header+="<th>[<a href=\"%s/header_caption.html#js\">js</a>,<a href=\"%s/header_caption.html#t\">t</a>)<br>[<a href=\"%s/header_caption.html#ss\">ss</a>,<a href=\"%s/header_caption.html#t\">t</a>)</th><th>(<a href=\"%s/header_caption.html#t\">t</a>,<a href=\"%s/header_caption.html#je\">je</a>]<br>(<a href=\"%s/header_caption.html#t\">t</a>,<a href=\"%s/header_caption.html#se\">se</a>]</th>"%(cWB_conf.web_link, cWB_conf.web_link, cWB_conf.web_link, cWB_conf.web_link, cWB_conf.web_link, cWB_conf.web_link, cWB_conf.web_link, cWB_conf.web_link)
        header+="<th>DQ</th>"
        if (id):
          header+="<th><a href=\"%s/header_caption.html#S\" target=_new>S</a></th>"%(cWB_conf.web_link)
        try:
          pe_par=cWB_conf.pe_par
          do_pe=True
        except:
          do_pe=False
        if (do_pe==True):
          header+="<th>Parameter Estimation<br>%s</th>"%"<br>".join(map(lambda x: "value_"+x, ifos))
        if (id):
          header+="<th>comments</th></tr>"
        return header
    
    def compute_ced(self):
        ced_link = "%s/OUTPUT/%s"%(self.job.link,self.compute_ced_dir())
        return ced_link

    def compute_ced_dir(self):
        print "In compute_ced_link"
        [seg_start,seg_stop]=self.job.dir.split("/")[-1].split("-")
        seg_white=self.job.back
        ced_dir = "ced_%s_%s_%s-%s_slag0_lag0_1_job1/"%(seg_white,cWB_conf.seg_duration, seg_start, seg_stop)
        for i in range(len(cWB_conf.ifos)):
            ced_dir += cWB_conf.ifos[i]
        for i in range(len(cWB_conf.ifos)):
            ced_dir += "_%.3f"%(float(self.start[i]))
        return ced_dir

    def get_ced_link(self):
        print "In get_ced_link"
        t="%s/OUTPUT/ced*/*%.3f*"%(self.job.dir,float(self.start[0]))
        a=glob.glob(t)
        if(len(a)>1):
            print "Warning: there are more than one CEDs: %s"%("|".join(a))
        elif(len(a)==0):
            print "Warning: No CEDs found satisfying %s"%(t)
            return "NA"
        print "Found ced: %s"%(a[0])
        ced_link="https://ldas-jobs.ligo.caltech.edu/~waveburst/" + "/".join(a[0].split("/")[3:])
        return ced_link

    def to_html_table_row(self,id=-1):
        #print "In to_html_table_row"
        #print self.start
        self.id=id
        #print self.job
        sys.stdout.flush();sys.stderr.flush()        
        #print self.job.dir
        sys.stdout.flush();sys.stderr.flush()
        t="%s/OUTPUT/ced*/*%.3f*"%(self.job.dir,float(self.start[0]))
        #print "In to_html_table_row t=%s"%(t)
        a=glob.glob(t)
        #print "a = %s"%(repr(a))
        sys.stdout.flush();sys.stderr.flush()
        gracedbid=""
        try:
            if(len(self.graceid)>0):
                #gracedbid="<a href=\"https://gracedb.ligo.org/events/view/%s\">%s</a>"%(self.graceid,self.graceid)
                gracedbid="<a href=\"%s\">%s</a>"%(self.gracedb_link,self.graceid)
            if(len(self.first_graceid)>0):
                #gracedbid="<a href=\"https://gracedb.ligo.org/events/view/%s\">%s</a>"%(self.first_graceid,self.first_graceid)
                gracedbid="<a href=\"%s\">%s</a>"%(self.gracedb_link,self.first_graceid)
        except:
            gracedbid="N/A"

        if(len(a)!=1):
            #print "No ced"
            self.ced_link="none"
        else:
            ced_path=a[0]

            sys.stdout.flush();sys.stderr.flush()
            tokens=ced_path.split("/")
            subdir="/".join(tokens[-3:])
            #print "subdir=%s"%(subdir)
            self.ced_link="<a href=\"%s/%s/\">ced</a>"%(self.job.link.replace("job.html",""),subdir)
            #print "ced_link=%s"%(self.ced_link)
            sys.stdout.flush();sys.stderr.flush()
        if(self.id==-1):
            ID=""
        else:
            ID="<td align=right><a href=\"%s\">%d</a></td><td align=right>%s</td>"%("%s/OUTPUT.merged/TRIGGERS/trigger_%s.txt"%(self.job.link.replace("job.html",""),self.time[0]), self.id, gracedbid)
            #print ID
        self.time_f=map(lambda x:"%.3f"%(float(x)), self.time)
        self.rho_f="%.1f"%float(self.rho[cWB_conf.id_rho])
        self.frequency_f=map(lambda x:"%.1f"%float(x), self.frequency)
        self.snr_f=map(lambda x:"%.1f"%math.sqrt(float(x)), self.snr)
        self.hrss_f=map(lambda x:"%.2g"%float(x), self.hrss)
        self.duration_f=map(lambda x: "%.3f"%float(x), self.duration)
        self.bandwidth_f=map(lambda x: "%d"%int(float(x)), self.bandwidth)        
        self.netcc_f=map(lambda x: "%.3f"%float(x), self.netcc)
        self.penalty_f="%.3f"%float(self.penalty[0])
        self.theta_f=map(lambda x: "%.1f"%float(x), self.theta)
        self.phi_f=map(lambda x: "%.1f"%float(x), self.phi)
        self.left_f=map(lambda x: "%.1f"%float(x),self.left)
        self.right_f=map(lambda x: "%.1f"%float(x),self.right)
        self.netED="%.3f"%float(self.energy_disbalance)
        self.fres=map(lambda x: str(int(x)/2), self.rate)
        self.erA_f=map(lambda x: "%.1f"%float(x), self.erA)
        self.utc=tconvert(self.time[0],int)
        self.erA_table()
        #self.error_region_gif()
        try:
          nclass=len(cWB_conf.Cuts_list)
          class_table="""<td align=center>%s</td comment="Class">"""%(self.cutclass)
          #print class_table
          #qveto=cWB_conf.qveto_plugin
          #th_qveto=cWB_conf.th_qveto
          #self.qveto_f="%.2f - %.2f"%(float(self.Qveto[0]),float(self.Qveto[len(cWB_conf.ifos)])) #map(lambda x:"%.2f"%float(x), self.Qveto)
          #for i in range(1,len(cWB_conf.ifos)):
          #    self.qveto_f+="<br>%.2f - %.2f"%(float(self.Qveto[i]),float(self.Qveto[i+len(cWB_conf.ifos)]))
          #qveto_table="""<td align=right>%s</td comment="Qveto">"""%(self.qveto_f)
          #qveto_table="""<td align=right>%s</td comment="Qveto">"""%("<br>".join(self.qveto_f)+"<br>")
        except:
          #self.qveto_f=""
          #qveto_table=""
          class_table=""
        try:
          pe_par=cWB_conf.pe_par
          pe_file="%s/OUTPUT_PE/pe_%.3f.html"%(self.job.dir,float(self.time[0]))
          if (os.path.exists(pe_file)):
            lines=open(pe_file).readlines()
            pe_table="<td align=\"center\">" 
            for line in lines:
              pe_table+=line
            pe_table+="</td>" 
          else:
            pe_table="<td align=\"left\">no pe</td>"
        except:
          pe_table=""
        try:
            #print "IFAR: sgnf[2]=%s"%(self.sgnf[1])
            rate=float(self.sgnf[2].split(" ")[1])
            num=int(self.sgnf[2].split(" ")[0])
            if(num>0):
                self.ifar1_f="%.1f"%(1./rate/(3600*24))
            else:
                self.ifar1_f=">%.1f"%(1./rate/(3600*24))
            self.ifar_latency_f="%.1f"%((self.creation_time-float(self.sgnf[2].split(" ")[-2]))/(3600.))
        except Exception,e:
            print "IFAR calculation error: %s"%e
            self.ifar1_f="none"
            self.ifar_latency_f="none"
        #if(len(self.DQ)>1):
        #    self.format_DQ()
        #    z=self.cat_passed_html.replace("cat","").replace("ifo","<a href=\""+self.DQ_http+"\">DQ</a>")
        if (self.cat_passed_html!=""):
            z=self.cat_passed_html
        else:
            z="n/a"
        #<td align=right>%s<br>%s</td comment="ced_link/qscans">
        msg = "<tr>"
        if (id!=-1):
          msg += "%s"%ID
          msg += """<td align=right>%s</td comment="ced_link">"""%(self.ced_link)
        msg += """<td align=right>%s</td comment="central time">"""%("<br>".join(self.time_f) + "<br><a href=\""+ self.job.link +"/job.html\">"+self.utc+"</a>")
        msg += """<td align=right>%s</td comment="rho">"""%(self.rho_f)
        msg += """<td align=right>%s</td comment="netCC">"""%(self.netcc_f[cWB_conf.id_cc])
        msg += """<td align=right>%s<br>%s</td comment="ifar">"""%(self.ifar1_f, self.ifar_latency_f)
        msg += """<td align=right>%s</td comment="frequency">"""%(self.frequency_f[0])
        msg += """<td align=right>%s</td comment="snr">"""%("<br>".join(self.snr_f)+"<br>")
        msg += """<td align=right>%s<br>%s</td comment="netED/penalty">"""%(self.netED, self.penalty_f)
        msg += """<td align=right>%s</td comment="hrss">"""%("<br>".join(self.hrss_f)+"<br>")
        msg += """<td align=right>%s<br>%s</td comment="duration/bandwidth">"""%(self.duration_f[0],self.bandwidth_f[0])
        msg += """<td align=right>%s</td comment="DEC/RA/phi">"""%("<br>".join([self.theta_f[2],self.phi_f[2]])+"<br>"+self.phi_f[0])
        #if (id==-1):
          #msg += """<td align=right>%s<br>%s<br>%s<br>%s</td comment="error_region_link">"""%(self.erA_link,self.error_region_link_txt,self.error_region_link_html,self.error_region_link_gif)
        msg += "%s"%(class_table)
        #msg += "%s"%(qveto_table)
        #msg += """
        #<td align=right>%s<br>%s</td comment="left">
        #<td align=right>%s<br>%s</td comment="right">
        #"""%(self.left_f[0], "%.1f"%float(self.left_science), self.right_f[0], "%.1f"%float(self.right_science))
        msg += """<td align=right>%s</td comment="DQ">"""%(z)
        if (id!=-1):
          msg += """<td align=right>%d</td comment="selected">"""%(self.selected)
        msg += "%s"%(pe_table)
        if (id!=-1):
          msg += """<td align=left>%s</td comment="comments">"""%(self.comments)
        msg += "</tr>"
       
        return msg
    
    def erA_table(self):
        self.erA_http="%s/OUTPUT.merged/TRIGGERS/erA_%s.html"%(self.job.link.replace("job.html",""),self.time[0])
        self.erA_link="<a href=\"%s\">erA</a>"%(self.erA_http)
        self.erA_fn="%s/erA_%s.html"%(self.dir,self.time[0])
        f=open(self.erA_fn,"w")
        print >>f,"<table border=1>"
        print >>f,"<tr><th>%</th><th>erA</th></tr>"
        p=range(10,100,10)
        for (pc,v) in zip(p,self.erA_f[1:-1]):
            row="<tr><td>%s</td><td>%s</td></tr>"%(pc,v)
            print >>f,row
        print >>f,"</table>"
        f.close()
            
    #def error_region_gif(self):
    #    self.error_region_gif_http="%s/OUTPUT.merged/TRIGGERS/error_region_%s.gif"%(self.job.link.replace("job.html",""),self.time[0])
    #    self.error_region_link_gif="<a href=\"%s\">gif</a>"%(self.error_region_gif_http)
    #    self.error_region_gif_fn="%s/OUTPUT.merged/TRIGGERS/error_region_%s.gif"%(self.job.dir,self.time[0])        

    def format_DQ(self):
        self.DQ_http="%s/OUTPUT.merged/TRIGGERS/DQ_%s.html"%(self.job.link.replace("job.html",""),self.time[0])
        self.DQ_link="<a href=\"%s\">cat</a>"%(self.DQ_http)
        self.DQ_fn="%s/OUTPUT.merged/TRIGGERS/DQ_%s.html"%(self.job.dir,self.time[0])
        self.DQ_txt_fn="%s/OUTPUT.merged/TRIGGERS/DQ_%s.txt"%(self.job.dir,self.time[0])

        ts=[]

        for dq in self.DQ:
            tokens=dq.split()
            ifo="?"
            type="?"
            try:
                ifo=tokens[0]
                type=tokens[1]
                description=" ".join(tokens[2:])
            except:
                print "Short DQ line: %s"%(dq)
            category="dummy"
            ts.append([category,ifo,type,description])

        ts.sort(cmp=comp_DQ)

        f=open(self.DQ_fn,"w")
        print >>f,"<table border=1><tr><th>category</th><th>ifo</th><th>type</th><th>description</th></tr>"
        for dq in ts:
            print >>f,"<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"%tuple(dq)
        print >>f,"</table>"
        f.close()
    def significance(self):
        ##  temporary fix
        ## self.sgnf[0]=-1
        ## self.sgnf[1]=-1
        ## self.sgnf[2]=-1
        ## self.sgnf[3]=-1
        ## return
        ##
        qvetodir="plot%s"%self.cutclass

        olddir=os.getcwd()
        bkg_postprod_dir="%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.postprod_dir)
        os.chdir("%s/python"%cWB_conf.online_dir)
        #aa=commands.getstatusoutput("./last_N.py -n 1 -m 1 -r 0 -o %s -R %f -d %s"%(bkg_postprod_dir,float(self.rho[cWB_conf.id_rho]),qvetodir))
        #a=commands.getstatusoutput("./last_N.py -n 1 -r 0 -o %s -R %f -d %s"%(bkg_postprod_dir,float(self.rho[cWB_conf.id_rho]),qvetodir))
        b=commands.getstatusoutput("./last_N.py -n 24 -r 0 -o %s -R %f -d %s"%(bkg_postprod_dir,float(self.rho[cWB_conf.id_rho]),qvetodir))
        #c=commands.getstatusoutput("./last_N.py -n 10000 -r 0 -o %s -R %f -d %s"%(bkg_postprod_dir,float(self.rho[cWB_conf.id_rho]),qvetodir))
        c=commands.getstatusoutput("./last_N.py -n 168 -r 0 -o %s -R %f -d %s"%(bkg_postprod_dir,float(self.rho[cWB_conf.id_rho]),qvetodir))
        f=open(self.path,"a")
        #f=open("%s/OUTPUT.merged/TRIGGERS/trigger_%s.txt"%(self.job.link.replace(cWB_conf.web_link,cWB_conf.online_dir),self.time[0]),"a")
        print >>f,"## number of events with rho>= rho of the event | corresponding rate | start time of the analyzed segment list | end time of the analyzed segment list | duration of the analyzed segment list"
        try:
           trials=cWB_conf.Cuts_trials
        except:
          try:
            trials=len(cWB_conf.Cuts_list)
          except:
            trials=1
        print >>f,"Applied trials factor: %i"%(trials)
        #print >>f,"#significance based on the last processed job (about %i seconds)"%(cWB_conf.seg_duration)
        #print >>f,aa[1]        
        #print >>f,"#significance based on the last hour"
        #print >>f,a[1]
        print >>f,"#significance based on the last day"
        print >>f,b[1]
        #print >>f,"#significance based on all the processed livetime in the run"
        print >>f,"#significance based on the last week"
        print >>f,c[1]
        try:
          num = int(b[1].split(" ")[0])
          if (num==0):
           print >>f,"far_is_upper_limit: True"
        except Exception,e:
          print "Error, significance not found"
          print e
        print >>f,"\n%s: %s"%(os.environ['SITE_CLUSTER'],self.path)
        lines=open("%s/input/burst.in"%(self.job.dir)).readlines()
        print >>f,"CAT1 %s"%("".join(lines))
        for ifo in cWB_conf.ifos:
          lines=open("%s/input/%s_cat2.in"%(self.job.dir,ifo)).readlines()
          print >>f,"CAT2 %s %s"%(ifo,"".join(lines))
        f.flush()
        f.close()
        print "After adding significance to the trigger"
        os.chdir(olddir)
        #self.sgnf[0]=aa[1]        
        #self.sgnf[1]=a[1]
        self.sgnf[2]=b[1]
        self.sgnf[3]=c[1]
    def record_passed(self):
        print "In Trigger.record_passed"
        f=open(self.path,"a")
        print >>f,str(self.cat_passed_txt)
        f.flush()
        f.close()
    def record_dq(self):
        print "In Trigger.record_dq"
        a=commands.getstatusoutput("ssh ldas-pcdev1i %s %s"%(cWB_conf.ligolw_dq_active_cats,self.time[0]))
        f=open(self.path,"a")
        lines=filter(lambda x: x.find("Welcome")==-1 and x.find("Agent")==-1, a[1].split("\n"))
        print >>f,"#DQ at the time of the online job"
        print >>f,"\n".join(lines)+"\n"
        f.flush()
        f.close()
        self.ligolw_dq_active_cats=a[1]

def ascii_dump_2_triggers(fn,segs,job):
    triggers=[]
    trigger=[]
    print "in ascii_dump_2_triggers: fn=%s"%fn
    lines=filter(lambda y: len(y)>0, map(lambda x: x.strip(), open(fn).readlines()))
    print "in ascii_dump_2_triggers: len(lines)=%d"%len(lines)
    for line in lines[1:]:
        if(line.find("# trigger")!=-1):
            if(len(trigger)!=0):
                triggers.append(Trigger(trigger,segs,job))
                print "---"
                print "in ascii_dump_2_triggers: len(trigger)=%d"%len(trigger)
#                print trigger
                print "---"
            trigger=[]
        else:
            trigger.append(line)
    print "in ascii_dump_2_triggers: End of loop"
    try:
        if(len(trigger)>0):
            triggers.append(Trigger(trigger,segs,job))
            print "---"
            print "in ascii_dump_2_triggers: len(trigger)=%d"%len(trigger)
#            print trigger
            print "---"
    except:
        pass
    return triggers


def comp_DQ(dq1,dq2):
    if(dq1[0]!=dq2[0]):
        return dq1[0]<dq2[0]
    if(dq1[1]!=dq2[1]):
        return dq1[1]<dq2[1]
    if(dq1[2]!=dq2[2]):
        return dq1[2]<dq2[2]
    return 0


if(__name__=="__main__"):
    segs=segmentlist([segment(700000000,999999999)])
    fn="/archive/home/igor/OUTPUT_JW1_L1H1V1_run65/wave_870497460_600_JW1_L1H1V1_run65_id551.txt"
    fn="../JOBS/9246/924697976-924698036/OUTPUT.merged/TRIGGERS/triggers.txt"
    triggers=ascii_dump_2_triggers(fn,segs,None)
    print len(triggers)
    for t in triggers:
        print t
        

