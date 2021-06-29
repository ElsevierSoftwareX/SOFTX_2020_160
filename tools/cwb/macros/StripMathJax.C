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


// Strip MathJax Tag from file
// MathJax javascript is used to convert tex formulae inside an html block
// to a figure to be displayed from www browser 
// Since the brakets { } are special characters for texinfo, all
// brackets used in the tex formulae must be defined with the escape 
// character @. To avoid to use the escape we use a workaround : 
// each formula must be commented with "@c MJX " or "@c mjx "
// The command "makeinfo --html convert the info comment into a html comment"
// This macro strip the comment "<!-- <MJ> " and "-->"  or  "<!-- <mj> " and "-->"

void StripMathJaxTag(TString ifile, TString ofile);

void StripMathJax(TString idir="cwb") {

  CWB::Toolbox TB;

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
    StripMathJaxTag(tmpFile,fileList[i]);
    sprintf(cmd,"rm %s",tmpFile.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  exit(0);
}

void StripMathJaxTag(TString ifile, TString ofile) {

  ifstream in;
  in.open(ifile.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << ifile.Data() << endl;exit(1);}

  ofstream out;
  out.open(ofile.Data(),ios::out);

  char str[2048];
  while(true) {
    in.getline(str,2048);
    if (!in.good()) break;
    //cout << str << endl;
    TString line = str;
    // remove MathJax html comment
    if(line.BeginsWith("<!-- <MJ> ")) {
      line.ReplaceAll("<!-- <MJ> ","");
      line.ReplaceAll("-->","");
    }
    if(line.BeginsWith("<!-- <mj> ")) {
      line.ReplaceAll("<!-- <mj> ","");
      line.ReplaceAll("-->","");
    }
    out << line.Data() << endl;
  }

  out.close();

  return;
}

