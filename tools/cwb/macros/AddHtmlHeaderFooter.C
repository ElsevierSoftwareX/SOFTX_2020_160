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



void AddHtmlHeaderFooter(TString idir="cwb", bool use_doc_header=true) {

  CWB::Toolbox TB;
  TString comment;
  Bool_t compile;
  TString fName;

  if(!idir.EndsWith("/")) idir+="/";

  vector<TString> fileList = TB.getFileListFromDir(idir, ".html");

  char cmd[1024];
  for(int i=0;i<fileList.size();i++) {
    TString tmpFile = fileList[i];
    tmpFile.ReplaceAll(".html",".html.tmp");
    //cout << fileList[i].Data() << endl;  
    //cout << tmpFile.Data() << endl;  
    sprintf(cmd,"mv %s %s",fileList[i].Data(),tmpFile.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"cat $HOME_WAT/html/etc/html/header.html %s $HOME_WAT/html/etc/html/footer.html > %s",
            tmpFile.Data(),fileList[i].Data());
    gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",tmpFile.Data());
    gSystem->Exec(cmd);
  }

  exit(0);
}

