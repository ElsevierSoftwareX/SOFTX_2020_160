#!/usr/bin/python

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

from glue.segments import *
from glue.segmentsUtils import *
import commands, sys, os, glob

sys.path.insert(0, os.getcwd())
import cWB_conf

segments_dir="%s/%s"%(cWB_conf.bkg_run_dir,cWB_conf.seg_dir)
run_segments_file="%s/%s"%(segments_dir,cWB_conf.run_segments_file)
f=open(run_segments_file)
lines = f.readlines()
f.close()
 
for line in lines:
    #print line,
    if(line.find("#")==-1):
        seg=line.split("\t")
        start=int(seg[1])
        stop=int(seg[2])
        n = int(stop-start)/cWB_conf.bkg_job_duration
        #print "%d %d %d"%(start,stop,n)
        s=start
        e=stop
        if (n>1):
            for i in range(0,n-1):
                s = start+ i*cWB_conf.bkg_job_duration
                e = start+(i+1)*cWB_conf.bkg_job_duration
                print "%d %d"%(s,e)
            s = e
            e = stop
        n = int(e-s)/cWB_conf.bkg_job_duration
        if (n>1):
            print "error"
        if (n==0):
            if (e-s>cWB_conf.bkg_job_minimum):
                print "%d %d"%(s,e)
        if (n==1):
            if (cWB_conf.bkg_job_minimum==cWB_conf.bkg_job_duration):
                e = s+ cWB_conf.bkg_job_duration
                print "%d %d"%(s,e)
            else:
                half = (e-s)/2
                if (half>cWB_conf.bkg_job_minimum):
                     e = s+ half
                     print "%d %d"%(s,e)
                     s = e
                     e = s+ half
                     print "%d %d"%(s,e)
                else:
                     e = s+cWB_conf.bkg_job_duration 
                     print "%d %d"%(s,e)
