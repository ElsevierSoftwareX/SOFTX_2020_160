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


// convert file into html 
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

void cwb_mkhtml_file(TString file_name, TString options="") {

  // get title
  TString option_title = CWB::Toolbox::getParameter(options,"--title");
  option_title.ReplaceAll("#"," ");

  // get subtitle
  TString option_subtitle = CWB::Toolbox::getParameter(options,"--subtitle");
  option_subtitle.ReplaceAll("#"," ");

  // get multi, if multi=='true' the plots are presented with 2 columns format
  TString option_multi = CWB::Toolbox::getParameter(options,"--multi");
  option_multi.ToUpper(); 

  if(file_name=="") {
    cout << "cwb_mkhtml_file.C - Error : empty file name !!!" << endl;exit(1);
  }
  if(file_name.Contains("..")) {
    cout << "cwb_mkhtml_file.C - Error : file_name must not contains '..' !!!" << endl;exit(1);
  }
  if(file_name==".")  file_name=gSystem->WorkingDirectory();
  // check if file_name is a directory
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(file_name.Data(),&id,&size,&flags,&mt);
  // create an html file with the list of png plots contained in the file_name dir 
  if(flags&2 || file_name==".") { // is a directory
    TString odir = file_name;
    // get the list of png files
    vector<TString> pngList = CWB::Toolbox::getFileListFromDir(odir, ".png");
    if(!pngList.size()) {
      cout << "cwb_mkhtml_file.C - Error : no *.png files found !!!" << endl;exit(1);
    }
    // sort list
    pngList = CWB::Toolbox::sortStrings(pngList);

    TString texiName=odir+"/png_html_index.texi";
    bool overwrite=CWB::Toolbox::checkFile(texiName,true);
    if(!overwrite) gSystem->Exit(0);
    ofstream out;
    out.open(texiName,ios::out);
    out << "@c include texi macros" << endl;
    out << "@include macros.texi" << endl;
    out << "@include mathjax.texi" << endl << endl;
    out << "@*" << endl;
    if(option_title!="") {
      out << "@center @txtfont{"<<"<br>"<<option_title<<"<br>"<<", blue, h1}" << endl;
      if(option_subtitle!="") out << "@center @txtcolor{"<<option_subtitle.Data()<<",black}" << endl;
      out << "@*" << endl;
      out << "@drawline" << endl;
    }
    out << "@*" << endl;
    if(option_multi=="TRUE") {
      for(int i=0;i<pngList.size();i++) {
        TString name=pngList[i];name.ReplaceAll(".png","");
        name=gSystem->BaseName(name);
        if(i%2==0) {
          out << "@multitable @columnfractions .5 .5" << endl;
          out << "@item @center @txtfont{"<<i+1<<", red, h2}" << endl;
          out << "@center @displayimage{../.,"<<name<<",470}"<<endl;
        } else {
          out << "@tab @center @txtfont{"<<i+1<<", red, h2}" << endl;
          out << "@center @displayimage{../.,"<<name<<",470}"<<endl;
          out << "@end multitable" << endl;
        }
      }
      if(pngList.size()%2==1) out << "@end multitable" << endl;
    } else {
      for(int i=0;i<pngList.size();i++) {
        TString name=pngList[i];name.ReplaceAll(".png","");
        name=gSystem->BaseName(name);
        out << "@center @txtfont{"<<i+1<<", red, h2}" << endl;
        out << "@center @displayimage{../.,"<<name<<",1000}"<<endl;
      }
    }
    out.close();
    // convert texi into html
    TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
    TString exec_cmd = TString::Format("%s/cwb_mkhtml.csh %s wheader;rm %s",
                       cwb_scripts.Data(),texiName.Data(),texiName.Data());
    int ret=gSystem->Exec(exec_cmd);
    if(ret) {
      cout << "cwb_mkhtml_file.C : Error while executing cwb_mkhtml png_html_index.texi !!!" << endl;
      exit(1);
    }

    // Change default html title: CWB Display -> CWB Report 
    TString indexName=texiName;
    indexName.ReplaceAll(".texi","/index.html");
    char cmd[1024];
    sprintf(cmd,"sed -i 's/<title>CWB Display/<title>CWB Report/g' %s",indexName.Data());
    gSystem->Exec(cmd);

    gSystem->Exit(0);
  }

  if(!file_name.EndsWith(".C")) {
  } else if(!file_name.EndsWith(".cc")) {
  } else if(!file_name.EndsWith(".hh")) {
  } else if(!file_name.EndsWith(".c")) {
  } else if(!file_name.EndsWith(".h")) {
  } else {
    cout << "cwb_mkhtml_file.C - Error : file name must ends with .C/.c/.h/.cc/.hh !!!" << endl;exit(1);
  }
  // if file_name is not an absolute path add working dir path 
  if(!file_name.BeginsWith("/")) {
    TString tmp = file_name;
    file_name = TString(gSystem->WorkingDirectory())+"/"+tmp;
  }

  // create output html dir
  TString html_dir = file_name;
  html_dir.ReplaceAll(".C","");
  html_dir.ReplaceAll(".cc","");
  html_dir.ReplaceAll(".hh","");
  html_dir.ReplaceAll(".c","");
  html_dir.ReplaceAll(".h","");
  CWB::Toolbox::mkDir(html_dir,true);

  THtml html;
  html.SetEtcDir(gSystem->ExpandPathName("$HOME_WAT/html/etc/html"));
  html.SetProductName("CWB");
  TString html_input_dir=html_dir;
  html.SetInputDir(html_input_dir.Data());

  char cmd[1024];
  sprintf(cmd,"cp %s/html/etc/html/ROOT.css %s/",gSystem->ExpandPathName("$HOME_WAT"),html_dir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp %s/html/etc/html/ROOT.js %s/",gSystem->ExpandPathName("$HOME_WAT"),html_dir.Data());
  gSystem->Exec(cmd);

  // redirect stderr to /dev/null to getrid of messages produced by html.Convert
  fpos_t poserr; fflush(stderr); fgetpos(stderr, &poserr);
  int fderr = dup(fileno(stderr)); freopen("/dev/null", "w", stderr);
  // redirect stdout to /dev/null to getrid of messages produced by html.Convert
  fpos_t posout; fflush(stdout); fgetpos(stdout, &posout);
  int fdout = dup(fileno(stdout)); freopen("/dev/null", "w", stdout);

  // convert file to html
  char title[1024];
  sprintf(title,"<h2 class=\"convert\" align=\"center\"> %s </h2>",file_name.Data());
  html.Convert(file_name.Data(),"",html_dir,"", 0, title);

  // restore the stderr output
  fflush(stderr); dup2(fderr, fileno(stderr)); close(fderr);
  clearerr(stderr); fsetpos(stderr, &poserr);
  // restore the stdout output
  fflush(stdout); dup2(fdout, fileno(stdout)); close(fdout);
  clearerr(stdout); fsetpos(stdout, &posout);

  TString html_file = html_dir+"/"+gSystem->BaseName(file_name)+".html";
  cout << endl << "created html file : " << html_file << endl << endl;

  gSystem->Exit(0);
}
