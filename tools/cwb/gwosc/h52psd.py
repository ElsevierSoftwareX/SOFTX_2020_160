# Copyright (C) 2020 Gabriele Vedovato
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

# extract psd from h5 file and save to dat file
# uses pesummary package

from pesummary.io import write
from pesummary.gw.file.read import read
import sys
import os

ifile=sys.argv[1]
model=sys.argv[2]
ofile_tag=sys.argv[3]

def main(ifile,model,ofile_tag):

    cwd = os.getcwd()

    print('\nPlease wait while reading file '+cwd+'/'+ifile+' ...')

    f = read(ifile)
    print(f.labels)

    print('\nPlease wait while reading psd '+model+' from file '+cwd+'/'+ifile+'...')

    psd_dict = f.psd[model]
    psd_ifos = sorted(list(psd_dict.keys()))
    print(psd_ifos)


    for ifo in psd_ifos:
      print('\nPlease wait while writing psd '+model+' to file '+cwd+'/'+ofile_tag+'_'+ifo+'.dat'+' ...')
      psd_dict[ifo].save_to_file(ofile_tag+"_"+ifo+".dat", delimiter="\t")

main(ifile,model,ofile_tag)
