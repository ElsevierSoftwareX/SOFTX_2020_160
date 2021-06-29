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

import commands, sys

log=sys.argv[1]
executable=sys.argv[2]
dir=sys.argv[3]
opt=sys.argv[4]

com="fuser %s"%(log)
a=commands.getstatusoutput(com)
# print a
processes=[]
try:
    processes=filter(lambda x: len(x)>0, map(lambda x: x.strip(), a[-1].split(":")[-1].split(" ")))
#    print processes
except Exception,e:
    pass
#    print e
if(len(processes)==0):
    t=commands.getstatusoutput("date")[-1].strip()
    f=open(log,"a")
    print >>f,"=|"*30
    print >>f, "Restarting at %s"%(t)
    print >>f,"=|"*30    
    f.close()
    com="%s %s %s %s"%(executable,log,dir,opt)
    print com
    a=commands.getstatusoutput(com)
    print a[1]
