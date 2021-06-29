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


// dump the list of labels of the report files in the report dir : used by the commands cwb_report
{
  CWB::Toolbox TB;

  vector<TString> repList = TB.getFileListFromDir(pp_dir, "M");

  TString www_rep = TString(www_dir)+TString("/")+TString(data_label)+TString("/")+TString(pp_dir);
  www_rep.ReplaceAll(report_dir+TString("/"),"");
  TString home="";
  if(gSystem->Getenv("HOME")==NULL) {
    cout << "Error : environment HOME is not defined!!!" << endl;exit(1);
  } else {
    home=TString(gSystem->Getenv("HOME"));
  }
  www_rep.ReplaceAll("~",home);
  cout << endl;
  cout << www_rep << endl;
  vector<TString> wwwList = TB.getFileListFromDir(www_rep, "M");
  cout << endl;
  cout << " -----------------------------" << endl;
  cout << " List of report directories" << endl;
  cout << " -----------------------------" << endl;
  cout << endl;

  for(int i=0;i<repList.size();i++) {
    repList[i].ReplaceAll(pp_dir,"");
    repList[i].ReplaceAll("/","");
    bool published=false;
    for(int j=0;j<wwwList.size();j++) {
      wwwList[j].ReplaceAll(www_rep,"");
      wwwList[j].ReplaceAll("/","");
      if(wwwList[j].CompareTo(repList[i])==0) published=true;;
    }
    if(published) {
      cout << " x : " << repList[i].Data() << endl;
    } else {
      cout << " - : " << repList[i].Data() << endl;
    }
  }
  cout << endl;
  cout << " -----------------------------" << endl;
  cout << " x/- : published/not-published" << endl;
  cout << endl;

  exit(0);
}
