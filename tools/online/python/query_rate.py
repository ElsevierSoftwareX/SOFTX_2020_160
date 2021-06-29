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
from last_N import *

def intersects(path,interval):
	s=segmentlist([segment(map(int,path.split("/")[-2].split("-")))])
	intersection=s & interval; intersection.coalesce()
	if(abs(intersection)>0):
		return True
	else:
		return False

def print_histogram(H,rho):
	print "histogram segment list=%s"%(repr(H.segs))
	print "livetime in histogram=%g"%(abs(H.segs))
	print "query segment list=%s"%(repr(interval))
	print "livetime in query=%g"%(abs(interval))
	z=H.significance(options.rho)
	print """rate=%.3g events/day, rate_error=%.3g events/day, n=%d, non-zero lag livetime=%d"""%(z[1]*3600*24,math.sqrt(z[0])/H.nonzero_LT*3600*24, z[0], abs(H.nonzero_LT))

if(__name__=="__main__"):
        parser = OptionParser()
	parser.add_option("-s","--start",dest="start",type=int, help="start of the interval on which FAR is computed")
	parser.add_option("-e","--end",dest="end",type=int, help="end of the interval on which FAR is computed")
	parser.add_option("-r","--rho",dest="rho",type=float,help="for which rho to estimate FAR")
	#parser.add_option("-o","--output",dest="output",type="string",help="if specified, a histogram will be dumped into this file")
	parser.add_option("-d","--dir",dest="dir",type="string",help="directory where save information")	
        (options, args) = parser.parse_args()
	if(not options.start):
		parser.error("start must be specified")
	if(not options.end):
		parser.error("end must be specified")
	if(not options.rho):
		parser.error("rho must be specified")
	if(not options.dir):
		parser.error("dir must be specified")

	interval=segmentlist([segment(options.start, options.end)])
	
        job_dir="%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.jobs_dir)
	#a=map(lambda x: "/".join(x.split("/")[:-1]), glob.glob("%s/?????/*/output/finished"%job_dir))
	a=map(lambda x: "%s/output"%("/".join(x.split("/")[:-1])), glob.glob("%s/?????/*/finished"%job_dir))
	#print a
	s=map(lambda x: "/".join(x.split("/")), glob.glob("%s/?????/*/start"%job_dir))
	#print s

	considered_jobs=filter(lambda x: intersects(x,interval), a)
	#print type(considered_jobs)
	#print "considered jobs=%s"%(repr(considered_jobs))
	started_jobs=filter(lambda x: intersects(x,interval), s)
	#print "started jobs=%s"%(repr(started_jobs))

	if(len(considered_jobs)==0):
		print "No completed jobs found in the specified time interval"
		sys.exit(0)

        considered_segments=segmentlist(map(lambda x: segment(map(int,x.split("/")[-2].split("-"))), considered_jobs));considered_segments.coalesce()
        #print considered_segments
        started_segments=segmentlist(map(lambda x: segment(map(int,x.split("/")[-2].split("-"))), started_jobs));started_segments.coalesce()
        #print started_segments
        missing_segments=started_segments - considered_segments
        #print "missing segments=%s"%(repr(len(missing_segments)))
        rootfiles=[]
        for b in considered_jobs:
                rfiles=glob.glob("%s/w*.root"%b)
                rootfiles+=rfiles
        print options.dir
        segments_file="%s/segments.txt"%(options.dir)
        dodir=1
        if(os.path.exists(segments_file)):
             #done = open("%s/%s/seg_global.txt"%(cWB_conf.online_dir,cWB_conf.seg_dir))
             done = open(segments_file)
             done_segments = fromsegwizard(done); done_segments.coalesce()
             done.close()
             notdone=considered_segments & ~done_segments; notdone.coalesce()
             if (len(notdone)==0):
                 dodir=0
        print "dodir: %i"%dodir
        if(dodir==1):    
            #commands.getstatusoutput("rm -rf %s;mkdir %s"%(options.dir,options.dir))
            commands.getstatusoutput("mkdir %s"%(options.dir))
            f=open("%s/segments.txt"%(options.dir),"w")
            for seg in considered_segments:
                print >>f,"%s %s"%(seg[0],seg[1])
            f.close()
            merge(options.dir,rootfiles)
        if(len(missing_segments)>0):
            f=open("%s/missing.txt"%(options.dir),"w")
            for seg in missing_segments:
                print >>f,"%s %s"%(seg[0],seg[1])
            f.close()
        else:
            commands.getstatusoutput("rm %s/missing.txt"%(options.dir))

#		c="%s/OUTPUT.merged/histogram.pickle"%("/".join(b.split("/")[:-1]))
#		print "c=%s"%(repr(c))
#		try:
#			f=open(c,"rb")
#			h=pickle.load(f)
#			f.close()
#			print "="*15+" job " +"="*15
#			print_histogram(h,options.rho)
#		except:
#			continue
#		try:
#			histograms.add(h)
#		except Exception,e:
#			histograms=h
#
#	print "="*40 + " merged " + "="*40
#	try:
#		print_histogram(histograms,options.rho)
#	except:
#		print "Exiting..."
#		sys.exit(1)
#	if(options.output):
#		f=open(options.output,"w")
#		print >>f,histograms
#		f.flush()
#		f.close()
#	if(options.dump):
#		histograms.dump(options.dump)
#		
