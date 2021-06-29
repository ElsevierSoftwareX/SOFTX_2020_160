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

import commands, cWB_conf

def get_time():
    com="lalapps_tconvert now"
    a=commands.getstatusoutput(com)
    return int(a[1].strip()) - cWB_conf.run_offset

def get_time_nooff():
    com="lalapps_tconvert now"
    a=commands.getstatusoutput(com)
    return int(a[1].strip()) 

def tconvert(gps,gtype=int):
    com="lalapps_tconvert %s"%(gtype(float(gps)))
    a=commands.getstatusoutput(com)
    return a[1].strip()

def t_cmp_threshold(a,b):
    return -cmp(float(a.rho[cWB_conf.id_rho]),float(b.rho[cWB_conf.id_rho]))

def t_cmp_time(a,b):
    return cmp(float(a.time[1]),float(b.time[1]))

def totable(header,rows):
    t="<table border=1>"
    h="<tr><th>"+"</th><th>".join(header)+"</th></tr>"
    b=map(lambda x:"<tr><td>"+"</td><td>".join(x)+"</td></tr>",rows)
    f="</table>"
    return t+"\n"+h+"\n".join(b)+"\n"+f

def add_counter_column(rows):
    return map(lambda x: [str(x[0])]+x[1], zip(range(1,len(rows)+1),rows))

          
