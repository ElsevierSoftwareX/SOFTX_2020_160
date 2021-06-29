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

#import conf
if (os.environ['SITE_CLUSTER']=="CASCINA"):
  from glue.segments import *
  from glue.segmentsUtils import *
else:
  from ligo.segments import *
  from ligo.segments.utils import *
import math


class RHO_HISTOGRAM:
	def __init__(self,dir):
		print "In RHO_HISTOGRAM %s"%(dir)
		self.dir=dir
		self.lt={}
		self.nonzero_LT=0
		self.zero_LT=0
		self.lt={}
		self.triggers={}
		self.triggers0={}
		self.ctriggers={}
		self._livetime()
		self._triggers()
	def _livetime(self):
		f=open("%s/livetime.txt"%self.dir)
		lines=map(lambda x: filter(lambda y: len(y)>0, x.split(" ")), f.readlines())
	        f.close()
		for line in lines:
			(lag0,lag1,lag2,livetime)=line
#			[lag0,lag1,lag2]=self._ifoorder([lag0,lag1,lag2])
			self.lt[(lag0,lag1,lag2)]=float(livetime)
			if(not (float(lag0)==0.0 and float(lag1)==0.0 and float(lag2)==0.0)):
				self.nonzero_LT+=float(livetime)
			else:
				self.zero_LT+=float(livetime)

	def __repr__(self):
		llags=self.lt.keys()
		llags.sort()
		tlags=self.triggers.keys()
		tlags.sort(cmp=lambda a,b: -cmp(a,b))
		msg="""#livetime in different lags
%s
# lag0 lag1 lag2 livetime
nonzero_livetime = %.4f
zero_livetime = %.4f
# number of triggers in different rho bins (step = 0.1)
%s
# number of triggers with rho greater or equal than current
%s
# segments
%s
"""%("\n".join(map(lambda x: "%.4f %.4f %.4f %.4f"%(float(x[0]),float(x[1]),float(x[2]),float(self.lt[x])), llags)), \
     self.nonzero_LT, self.zero_LT,\
     "\n".join(map(lambda x: "%.1f %d"%(float(x),self.triggers[x]), tlags)), \
     "\n".join(map(lambda x: "%.1f %d"%(float(x),self.ctriggers[x]), tlags)), repr(self.segs))
		return msg
	
	def _triggers(self):
		f=open("%s/triggers.txt"%self.dir)
		lines=filter(lambda x:x.find("#")==-1,f.readlines())
		f.close()
		for l in lines:
			line=filter(lambda x:len(x)>0, l.strip().split(" "))
			(lag0,lag1,lag2,rho)=(line[5],line[6],line[7],line[0])
			if(float(rho)>100):
				print line
				print "%s %s %s %s"%(lag0,lag1,lag2,rho)
			if(not (float(lag0)==0.0 and float(lag1)==0.0 and float(lag2)==0.0)):
				try:
					self.triggers["%.1f"%float(rho)]+=1
				except:
					self.triggers["%.1f"%float(rho)]=1
			else:
				try:
					self.triggers0["%.1f"%float(rho)]+=1
				except:
					self.triggers0["%.1f"%float(rho)]=1
		rhos=self.triggers.keys()
		rhos.sort(cmp=lambda a,b: -cmp(float(a),float(b)))
		for r,i in zip(rhos,range(len(rhos))):
			if(i==0):
				self.ctriggers[r]=self.triggers[r]
			else:
				self.ctriggers[r]=self.triggers[r]+self.ctriggers[rhos[i-1]]
		self.segs=segmentlist([segment(map(lambda x: int(float(x)), self.dir.split("/")[-2].split("-")))])
		print "--------->segs=%s"%(repr(self.segs))
		print "self.dir=%s"%(self.dir)
	def add(self,other):
		for k in self.lt.keys():
			if(other.lt.has_key(k)):
				self.lt[k]+=other.lt[k]
		for k in other.lt.keys():
			if(not self.lt.has_key(k)):
				self.lt[k]=other.lt[k]
		self.nonzero_LT+=other.nonzero_LT
		self.zero_LT+=other.zero_LT
		for k in self.triggers.keys():
			if(other.triggers.has_key(k)):
				self.triggers[k]+=other.triggers[k]
		for k in other.triggers.keys():
			if(not self.triggers.has_key(k)):
				self.triggers[k]=other.triggers[k]
		for k in self.triggers0.keys():
			if(other.triggers0.has_key(k)):
				self.triggers0[k]+=other.triggers0[k]
		for k in other.triggers0.keys():
			if(not self.triggers0.has_key(k)):
				self.triggers0[k]=other.triggers0[k]
		self.ctriggers={}
		rhos=self.triggers.keys()
		rhos.sort(cmp=lambda a,b: -cmp(float(a),float(b)))
		for r,i in zip(rhos,range(len(rhos))):
			if(i==0):
				self.ctriggers[r]=self.triggers[r]
			else:
				self.ctriggers[r]=self.triggers[r]+self.ctriggers[rhos[i-1]]
		self.segs=self.segs | other.segs; self.segs.coalesce()
		return self
	def significance(self, rho):
		bin="%.1f"%(float(rho))
		try:
			n=self.ctriggers[bin]
		except:
			rhos=self.triggers.keys()
			rhos.sort(cmp=lambda x,y: cmp(float(x),float(y)))
			if(rho<float(rhos[0])):
				n=self.ctriggers[rhos[0]]
			elif(rho>float(rhos[-1])):
				n=0
			else:
				for i in range(1,len(rhos)):
					if(float(rhos[i-1])<=rho and rho<float(rhos[i])):
						n=self.ctriggers[rhos[i-1]]
						break
		return [n, n/float(self.nonzero_LT)]

	def dump(self,fn):
		f=open(fn,"w")
		rhos=self.triggers.keys()
		rhos.sort(cmp=lambda x,y: cmp(float(x),float(y)))
		for rho in rhos:
			print >>f,"%s %d %g %g"%(rho,self.ctriggers[rho],self.ctriggers[rho]/float(self.nonzero_LT),math.sqrt(self.ctriggers[rho])/float(self.nonzero_LT))
		f.close()
															
			
		
if(__name__=="__main__"):
	hs=["histogram_10000_0_6.pickle", "histogram_1_0_1.pickle", "histogram_1_0_6.pickle", "histogram_24_0_6.pickle"]
	for h in hs:
		try:
			f=open(h,"rb")
		except:
			continue
		ht=pickle.load(f)
		f.close()
		fn=h.replace("pickle","txt")
		ht.dump(fn)
		fn=h.replace("pickle","segs")
		f=open(fn,"w")
		tosegwizard(f,ht.segs)
		f.close()

