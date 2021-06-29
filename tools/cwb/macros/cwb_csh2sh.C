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


vector<TString> getFileListFromDir(TString dir_name, TString endString="", TString beginString="");
TString RemoveMultipleBlank(TString string, TString left, TString right);

void cwb_csh2sh(TString script = "") {

  vector<TString> cshList = getFileListFromDir(".", ".csh");

  for(int i=0;i<cshList.size();i++) {

    TString ifile=cshList[i];
    cout << "ifile : " << ifile.Data() << endl;
    if(script!="") if(!ifile.Contains(script)) continue;

    ifstream in;
    in.open(ifile.Data(),ios::in);
    if (!in.good()) {cout << "Error Opening File : " << ifile.Data() << endl;exit(1);}

    TString ofile=ifile;
    ofile.ReplaceAll(".csh",".sh");
    cout << "ofile : " << ofile.Data() << endl;

    ofstream out;
    out.open(ofile.Data(),ios::out);
    if (!out.good()) {cout << "Error Opening File : " << ofile.Data() << endl;exit(1);}

    char istring[1024];
    while (1) {
      in.getline(istring,1024);
      if (!in.good()) break;
      TString ostring(istring);
      if(ostring.Contains("_SKIP_CSH2SH_")) {
        out << ostring.Data() << endl;
        continue;
      }
      if(!ostring.Contains("_SKIP_csh2sh_")) ostring.ReplaceAll(".csh",".sh");
      ostring.ReplaceAll("endif","fi");
      ostring.ReplaceAll("unsetenv","unset");
      ostring.ReplaceAll("#!/bin/tcsh","#!/bin/bash");

      bool tag=true;
      if(ostring.Contains("if")&&ostring.Contains("(")&&ostring.Contains(")")) {
        tag=true;
        for(int j=0;j<ostring.Sizeof();j++) if(ostring[j]=='(') {if(tag) tag=false; else ostring[j]=' ';}
        tag=true; 
        for(int j=ostring.Sizeof()-1;j>=0;j--) if(ostring[j]==')') {if(tag) tag=false; else ostring[j]=' ';}
   
        // remove multiple blank characters within "! $?"  
        ostring = RemoveMultipleBlank(ostring,"!","$?");
 
        ostring.ReplaceAll("! $?","-z $");

        ostring.ReplaceAll("$1 == \'\'","-z $1");
        ostring.ReplaceAll("$2 == \'\'","-z $2");
        ostring.ReplaceAll("$2 != \'\'","-n $2");
        ostring.ReplaceAll("$3 == \'\'","-z $3");
        ostring.ReplaceAll("$3 != \'\'","-n $3");
        ostring.ReplaceAll("(","[ ");
        ostring.ReplaceAll("==","=");
        ostring.ReplaceAll("\'","\"");
        ostring.ReplaceAll(")"," ]");
        ostring.ReplaceAll(" then","; then");
        if(ostring.ReplaceAll("\"$PATH\" !~","! $PATH =~"))
        {ostring.ReplaceAll(" [ "," [[ ");ostring.ReplaceAll(" ];"," ]];");}
        if(ostring.ReplaceAll("\"$LD_LIBRARY_PATH\" !~","! $LD_LIBRARY_PATH =~"))
        {ostring.ReplaceAll(" [ "," [[ ");ostring.ReplaceAll(" ];"," ]];");}
        if(ostring.ReplaceAll("\"$MANPATH\" !~","! $MANPATH =~"))
        {ostring.ReplaceAll(" [ "," [[ ");ostring.ReplaceAll(" ];"," ]];");}
        if(ostring.ReplaceAll("\"$PYTHONPATH\" !~","! $PYTHONPATH =~"))
        {ostring.ReplaceAll(" [ "," [[ ");ostring.ReplaceAll(" ];"," ]];");}
        if(ostring.ReplaceAll("\"$INFOPATH\" !~","! $INFOPATH =~"))
        {ostring.ReplaceAll(" [ "," [[ ");ostring.ReplaceAll(" ];"," ]];");}
        ostring.ReplaceAll("*","");
        ostring.ReplaceAll(" !~ "," != ");
      }

      TRegexp re4("set.*=.*\".*setenv.*$");
      if(ostring.Contains(re4)) {
        ostring.ReplaceAll("setenv ","export ");
        ostring.ReplaceAll("set ","");
        ostring.ReplaceAll("$","= $");
        ostring.ReplaceAll(" ","");
        ostring.ReplaceAll("\"export","\"export ");
        out << ostring.Data() << endl;
        continue;
      }

      ostring.ReplaceAll("exit","return 0");
      ostring.ReplaceAll("setenv","export");
      ostring.ReplaceAll("\t"," ");
      ostring.ReplaceAll("else if","elif");

      int iblanks=0;
      for(int j=0;j<ostring.Sizeof();j++) if(ostring[j]==' ') iblanks++; else break;
      if(!ostring.Contains("$#")) if(ostring.First("#")>0) ostring.Remove(ostring.First("#"));
      if((ostring.Contains("export")&&!ostring.Contains("unexport"))||ostring.Contains("alias")){
        TObjArray* token = TString(ostring).Tokenize(TString(' '));
        TObjString* tok[20];
        for(int j=0;j<token->GetEntries();j++){ tok[j]=(TObjString*)token->At(j);}
        ostring = tok[0]->GetString()+" "+tok[1]->GetString()+"="+tok[2]->GetString();

        if(token->GetEntries()>3){for(int j=3;j<token->GetEntries();j++)
          {ostring.Append(" ");ostring.Append(tok[j]->GetString());}}
      } else iblanks=0;

      TRegexp re1("set.*=.*basename");
      if(ostring.Contains(re1)) {
        ostring.ReplaceAll("set","");
        ostring.ReplaceAll(" ","");
        ostring.ReplaceAll("basename","basename ");
      }

      TRegexp re2("set.*=.*$");
      if(ostring.Contains(re2)&&!ostring.Contains("alias")) {
        ostring.ReplaceAll("set","");
        //ostring.ReplaceAll(" ","");
      }

      TRegexp re3("@.*[-=+*:].*[0-9]");
      if(ostring.Contains(re3)) {
        ostring.ReplaceAll("@","(( ");
        ostring = ostring+" ))";
      }

      for(int j=0;j<iblanks;j++) out << " ";    // restore initials blanks
      out << ostring.Data() << endl;
    }
    out.close();
    in.close();
  }
  exit(0);
}

TString RemoveMultipleBlank(TString string, TString left, TString right) {

  TString reg0=left+" .*"+right;
  TString pat0=left+" "+right;

  // remove multiple blank characters within "left right"
  TRegexp t(reg0);
  int index = string.Index(t);
  if(index>=0) {
     TString pat=left;
     for(int j=index+1; j<string.Sizeof();j++) {if(string[j]=='$') break; pat+=" ";}
     pat+=right;
     string.ReplaceAll(pat,pat0);
  }
  return string;
}  

vector<TString> getFileListFromDir(TString dir_name, TString endString, TString beginString) {

  TString wdir = gSystem->WorkingDirectory();
  TSystemDirectory gdir("", dir_name);
  TList *dfiles = gdir.GetListOfFiles();
  TIter dnext(dfiles);
  TSystemFile *dfile;
  TString fname;
  vector<TString> fileList;
  char path[1024];

  while ((dfile = (TSystemFile*)dnext())) {
    fname = dfile->GetName();
    sprintf(path,"%s/%s",dir_name.Data(),fname.Data());
    bool fsave=true;
    if ((endString!="")&&!fname.EndsWith(endString)) fsave=false;
    if ((beginString!="")&&!fname.BeginsWith(beginString)) fsave=false;
    if(fsave) fileList.push_back(path);
  }

  gSystem->ChangeDirectory(wdir);  // restore original dir

  return fileList;
}

