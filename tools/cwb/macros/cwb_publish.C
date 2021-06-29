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


// obsolete 

{
  CWB::Toolbox TB;
  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  TString cwb_pp_data_dir=TString(gSystem->Getenv("CWB_PP_DATA_DIR"));
  if(cwb_pp_data_dir.Sizeof()>1) {
    sprintf(pp_dir,"%s/%s",pp_dir,cwb_pp_data_dir.Data());
  }


  TString spp_dir = pp_dir;
  spp_dir.ReplaceAll(report_dir+TString("/"),"");

  char cmd[512];
  if(TString(gSystem->Getenv("CWB_PUBLISH_OPTION")).CompareTo("clean")==0) {
    cout << "clean publish dir" << endl;
    sprintf(cmd,"rm %s/%s/%s",www_dir,data_label,spp_dir.Data());
    cout << endl;
    cout << cmd << endl;
    gSystem->Exec(cmd);
  } else {
    TB.checkFile(pp_dir);

    sprintf(cmd,"mkdir %s/%s",www_dir,data_label);
    cout << endl;
    cout << cmd << endl;
    gSystem->Exec(cmd);

    sprintf(cmd,"ln -sf %s/%s %s/%s/%s",work_dir,pp_dir,www_dir,data_label,spp_dir.Data());
    cout << cmd << endl;
    cout << endl;
    gSystem->Exec(cmd);
  }

  exit(0);
}
