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


// remove broken symbolic links in condor log dir (avoid init condor failure)
// to be used when jobs are in held status

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  // condor log dirs
  char full_condor_log_dir[1024];
  sprintf(full_condor_log_dir,"%s/%s",work_dir,log_dir);

  vector<TString> fileList = TB.getFileListFromDir(full_condor_log_dir, "", "","",true);

  int cnt=0;
  char cmd[1024];
  for(int i=0;i<fileList.size();i++) {
    TString path;
    Long_t id,size,flags,mt;
    path = CWB::Toolbox::getFileName((char*)fileList[i].Data());
    if(path!="") {
      int estat = gSystem->GetPathInfo(path.Data(),&id,&size,&flags,&mt); 
      if(estat!=0) {	// condor log out,err symbolic link is broken
        sprintf(cmd,"rm -f %s",fileList[i].Data());
        gSystem->Exec(cmd);
        cnt++;
      }
    }
  } 
  cout << endl;
  cout << "Cleaned condor log files : " << cnt << endl;
  cout << endl;

  gSystem->Exit(0);
}
