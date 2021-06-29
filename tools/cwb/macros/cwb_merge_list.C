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


//  list the labels of the already merged files : used by the cwb_merge command

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  int iversion;
  TString cwb_merge_label;
  if(gSystem->Getenv("CWB_MERGE_LABEL")==NULL) {
    cout << "Error : environment CWB_MERGE_LABEL is not defined!!!" << endl;exit(1);
  } else {
    cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  }
  // check if label has the correct format (M#)
  if(cwb_merge_label[0]!='M') {
    cout << "Error : label " << cwb_merge_label.Data() << " has bad format (M#)" << endl;exit(1);
  } else {
    TString lcheck=cwb_merge_label;
    lcheck.Remove(0,1);
    if(!lcheck.IsDigit()) {
      cout << "Error : label " << cwb_merge_label.Data() << " has bad format (M#)" << endl;exit(1);
    } else {
      iversion=lcheck.Atoi();
    }
  }

  vector<int> jobList = TB.getMergeJobList(merge_dir,data_label,iversion);

  for(int i=0;i<jobList.size();i++) {
    cout << i << " jobId : " << jobList[i] << endl;
  }
  cout << endl;
  cout << "List Size : " << jobList.size() << endl;

  exit(0);
}
