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


{
  CWB::Toolbox TB;
  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  TString cwb_merge_label;
  if(gSystem->Getenv("CWB_MERGE_LABEL")==NULL) {
    cout << "Error : environment CWB_MERGE_LABEL is not defined!!!" << endl;exit(1);
  } else {
    cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  }

  char rep_dir[1024];
  sprintf(rep_dir,"%s/%s.%s",report_dir,cwb_merge_label.Data()); 
  cout << rep_dir << endl;
  if(TB.rmDir(rep_dir,true)==true) {
    cout << "clean publish dir" << endl;
    char cmd[1024];
    sprintf(cmd,"rm %s/%s/%s",www_dir,data_label,cwb_merge_label.Data());
    cout << endl;
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  exit(0);
}
