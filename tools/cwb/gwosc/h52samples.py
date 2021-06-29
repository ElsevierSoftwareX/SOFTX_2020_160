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

# extract posterior samples (model) from h5 file and save to dat file
# uses pesummary package

from pesummary.io import write
from pesummary.gw.file.read import read
import numpy as np
import sys
import os

ifile=sys.argv[1]
model=sys.argv[2]
ofile=sys.argv[3]

def main(ifile,model,ofile):

    cwd = os.getcwd()

    print('\nPlease wait while reading file '+cwd+'/'+ifile+' ...')

    f = read(ifile)
    print(f.labels)

    print('\nPlease wait while reading samples '+model+' from file '+cwd+'/'+ifile+' ...')

    samples_dict = f.samples_dict
    posterior_samples = samples_dict[model]
    parameters = sorted(list(posterior_samples.keys()))
    print(parameters)

    print('\nPlease wait while writing samples '+model+' to file '+cwd+'/'+ofile+' ...')

    parameters = list(posterior_samples.keys())
    samples_array = np.array([posterior_samples[param] for param in parameters]).T
    write(parameters, samples_array, package="core", file_format="dat", filename=ofile)

main(ifile,model,ofile)
