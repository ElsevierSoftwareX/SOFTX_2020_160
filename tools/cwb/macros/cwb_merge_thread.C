/*
# Copyright (C) 2019 Gabriele Vedovato
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
*/


// thread macro used to merge output root job files in multi thread mode 
// used by CWB::Toolbox::mergeCWBTrees

void cwb_merge_thread(TString flistName, bool simulation, TString odir, 
                      TString label, bool brms, bool bvar, bool bpsm) {

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  TB.checkFile(flistName);

  vector<TString> fileList = TB.readFileList(flistName);

  TB.mergeCWBTrees(fileList, simulation, odir, label, brms, bvar, bpsm);

  exit(0);
}
