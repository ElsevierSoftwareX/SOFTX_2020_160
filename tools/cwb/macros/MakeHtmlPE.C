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


#define WWW_PUBLIC		"https://ldas-jobs.ligo.caltech.edu/~waveburst/reports/"
#define WWW_LAG_MANUAL          "https://gwburst.gitlab.io/documentation/latest/html/faq.html?highlight=lag#what-are-lags-and-how-to-use-them"
#define WWW_SLAG_MANUAL         "https://gwburst.gitlab.io/documentation/latest/html/faq.html?highlight=lag#what-are-super-lags-and-how-to-use-them"


void mkhtml_index(vector<TString> gw_event_report, TString odir, TString label, TString title);
void ModifyFontSizeCSS(TString tabber);
int  ReadConfig(TString ifconfig, vector<TString> &gwname, vector<TString> &wdir, vector<TString> &rdir, vector<TString> &sdir);
void MakeHtmlStat(TString type, TString odir, TString label, TString ifconfig);
void MakeHtmlRep(TString type, TString odir, TString label, TString ifconfig);

void MakeHtmlPE(TString type, TString odir, TString label, TString ifconfig) {

  bool btype=false;

  if(type=="ced" || type=="postprod")                               {MakeHtmlRep(type, odir, label, ifconfig); btype=true;}
  if(type=="full" || type=="left" || type=="right" || type=="wave") {MakeHtmlStat(type, odir, label, ifconfig);btype=true;}

  if(!btype) cout << endl << "MakeHtmlPE.C - Error: type parameter ( " << type << " ) not a valid value" << endl << endl;

  exit(0);
}

void MakeHtmlRep(TString type, TString odir, TString label, TString ifconfig) {

  cout<<"MakeHtmlPE.C init ..."<<endl;

  CWB::Toolbox::checkFile(odir);	// check if odir exist

  // get CWB_USER_URL
  TString cwb_user_url;
  if(gSystem->Getenv("CWB_USER_URL")!=NULL) {
    cwb_user_url=TString(gSystem->Getenv("CWB_USER_URL"));
  }

  gSystem->Exec("date");

  vector<TString> gwname; 
  vector<TString> wdir;
  vector<TString> rdir;
  vector<TString> sdir;
  int gwsize = ReadConfig(ifconfig, gwname, wdir, rdir, sdir);

  vector<TString> link(gwsize);

  for(int i=0;i<gwsize;i++) link[i] = type.Contains("ced") ? cwb_user_url+"/"+sdir[i]+"/" : cwb_user_url+"/"+wdir[i]+"/";

  // create html index
  vector<TString> gw_event_report;

  TString title="";
  TString dir="";
  if(type.Contains("ced")) {
    title = "cWB - O3 BBH CED Reports";
    dir = "ced"; 
  } else {
    title = "cWB - O3 BBH PostProduction Reports";
    dir = "postprod"; 
  }

  char options[1024];
  int high = type.Contains("ced") ? 6700 : 6200;
  if(gwsize<10) {
    for(int n=0;n<gwsize;n++) {
      sprintf(options,"--link %s --label %s --name %s --high %d",link[n].Data(),gwname[n].Data(),dir.Data(),high);
      gw_event_report.push_back(options);
    }
  } else {
    int set=0;
    for(int n=0;n<gwsize;n++) {
      if(n%10==0) {
        if(set!=0) gw_event_report.push_back("--link </tab>/ --label dummy");
        gw_event_report.push_back(TString::Format("--link <tab>/ --label SET%d",++set));
      }
      sprintf(options,"--link %s --label %s --name %s --high %d",link[n].Data(),gwname[n].Data(),dir.Data(),high);
      gw_event_report.push_back(options);
    }
    gw_event_report.push_back("--link </tab>/ --label dummy");
  }

  mkhtml_index(gw_event_report, odir, label, title);

  exit(0);
}

void MakeHtmlStat(TString type, TString odir, TString label, TString ifconfig) {

  cout<<"MakeHtmlPE.C init ..."<<endl;

  CWB::Toolbox::checkFile(odir);	// check if odir exist

  // get CWB_USER_URL
  TString cwb_user_url;
  if(gSystem->Getenv("CWB_USER_URL")!=NULL) {
    cwb_user_url=TString(gSystem->Getenv("CWB_USER_URL"));
  }

  gSystem->Exec("date");

  vector<TString> gwname; 
  vector<TString> wdir;
  vector<TString> rdir;
  vector<TString> sdir;
  int gwsize = ReadConfig(ifconfig, gwname, wdir, rdir, sdir);

  vector<TString> link(gwsize);

  if(type.Contains("wave")) {
    if(odir.Contains("wave/time")) {
      for(int i=0;i<gwsize;i++) link[i] = cwb_user_url+"/"+wdir[i]+"/dump/"+rdir[i]+"/"+type+"/time/png_html_index/";
    } else if(odir.Contains("wave/spectrum")) {
      for(int i=0;i<gwsize;i++) link[i] = cwb_user_url+"/"+wdir[i]+"/dump/"+rdir[i]+"/"+type+"/spectrum/png_html_index/";
    } else if(odir.Contains("wave/envelope")) {
      for(int i=0;i<gwsize;i++) link[i] = cwb_user_url+"/"+wdir[i]+"/dump/"+rdir[i]+"/"+type+"/envelope/png_html_index/";
    } 
  } else {
    for(int i=0;i<gwsize;i++) link[i] = cwb_user_url+"/"+wdir[i]+"/dump/"+rdir[i]+"/"+type+"/distributions/png_html_index/";
  }

  // create html index
  vector<TString> gw_event_report;

  char options[1024];
  int high = type.Contains("wave") ? 4990 : 7900;
  if(gwsize<10) {
    for(int n=0;n<gwsize;n++) {
      sprintf(options,"--link %s --label %s --name index.html --high %d",link[n].Data(),gwname[n].Data(),high);
      gw_event_report.push_back(options);
    }
  } else {
    int set=0;
    for(int n=0;n<gwsize;n++) {
      if(n%10==0) {
        if(set!=0) gw_event_report.push_back("--link </tab>/ --label dummy");
        gw_event_report.push_back(TString::Format("--link <tab>/ --label SET%d",++set));
      }
      sprintf(options,"--link %s --label %s --name index.html --high %d",link[n].Data(),gwname[n].Data(),high);
      gw_event_report.push_back(options);
    }
    gw_event_report.push_back("--link </tab>/ --label dummy");
  }

  TString title="";
  if(type.Contains("wave")) {
    title = "cWB - O3 BBH Waveforms (" + label + ")";
  } else {
    title = "cWB - O3 BBH Matching Factors / Residuals (" + label + ")";
  }

  mkhtml_index(gw_event_report, odir, label, title);

  exit(0);
}

void mkhtml_index(vector<TString> gw_event_report, TString odir, TString label, TString title) {

  CWB::Toolbox TB;

  ofstream out;
  char ofile[1024];
  if(odir=="") {
    sprintf(ofile,"index.html");
  } else {
    sprintf(ofile,"%s/index.html", odir.Data());
  }
  cout << endl << "make index html file : " << ofile << endl << endl;
  out.open(ofile,ios::out);
  if (!out.good()) {cout << "MakeHtmlPE::mkhtml_index : Error Opening File : " << ofile << endl;exit(1);}

  // open input index template file
  char  html_index_template[1024]="";

  if(gSystem->Getenv("CWB_HTML_INDEX")==NULL) {
    cout << "MakeHtmlPE::mkhtml_index: Error - environment CWB_HTML_INDEX is not defined!!!" << endl;exit(1);
  } else {
    strcpy(html_index_template,gSystem->Getenv("CWB_HTML_INDEX"));
  }
  TB.checkFile(html_index_template);


  ifstream in;
  in.open(html_index_template,ios::in);
  if (!in.good()) {
    cout << "MakeHtmlPE::mkhtml_index : Error Opening File : " << html_index_template << endl;
    exit(1);
  }

  int ndiv_in=0;
  int ndiv_out=0;
  char istring[1024];
  while (1) {
    in.getline(istring,1024);
    if (!in.good()) break;
    TString ostring(istring);
    if(ostring.Contains("<div"))   ndiv_in++;
    if(ostring.Contains("</div")) ndiv_out++;
    if(ndiv_in>2 && ndiv_out<4) continue;
    out << ostring.Data() << endl;
  }


  out << "<html>" << endl;
  out << "<br>" << endl;
  out << "<div align=\"center\"><font color=\"red\"><h1>" << title << "</h1></font>" << endl;
  out << "<br>" << endl;


  // make tabber
  char sbody_height[256];
  sprintf(sbody_height,"%d",1900);
  out << "<div class=\"tabber\">" << endl;
  for(int i=0;i<gw_event_report.size();i++) if(gw_event_report[i]!="") {

    TString gw_event_report_link = CWB::Toolbox::getParameter(gw_event_report[i],"--link");
    if(gw_event_report_link=="" && i!=0) {
      cout<<"MakeHtmlPE::mkhtml_index : Error : gw_event_report --link not defined"<<endl;exit(1);}

    TString gw_event_report_label = CWB::Toolbox::getParameter(gw_event_report[i],"--label");
    if((gw_event_report_link!="</tab>/")&&(gw_event_report_label=="")) {
      cout<<"MakeHtmlPE::mkhtml_index : Error : gw_event_report --label not defined"<<endl;exit(1);}

    TString gw_event_report_high = CWB::Toolbox::getParameter(gw_event_report[i],"--high");
    if(gw_event_report_high=="") gw_event_report_high=sbody_height;
    int igw_event_report_high = gw_event_report_high.Atoi();

    TString gw_event_report_name = CWB::Toolbox::getParameter(gw_event_report[i],"--name");

    if(gw_event_report_link=="<tab>/") {           // open a sub tab
      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>" << gw_event_report_label << "</h2>" << endl;
      out << "<div class=\"tabber\">" << endl;
    } else if(gw_event_report_link=="</tab>/") {   // close sub tab
      out << "</div>" << endl;
      out << "</div>" << endl;
    } else {                                  // add a tab
      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>" << gw_event_report_label << "</h2>" << endl;

      if(gw_event_report_name=="") {
        out << "  <iframe src=\"" << gw_event_report_link << "header.html\" width=\"100%\" height=\"900px\" "
            << "marginwidth=\"15\" marginheight=\"15\" frameborder=\"0\"></iframe>" << endl;

        out << "  <iframe src=\"" << gw_event_report_link << "body.html\" width=\"100%\" "
            << " height=\"" << igw_event_report_high << "px\" frameborder=\"0\"></iframe>" << endl;
      } else {
        out << "  <iframe src=\"" << gw_event_report_link << gw_event_report_name << "\" width=\"100%\" "
            << " height=\"" << igw_event_report_high << "px\" frameborder=\"0\"></iframe>" << endl;
      }

      out << "</div>" << endl;
    }
  }
  out << "</div>" << endl;

  out << "</html>" << endl;

  in.close();
  out.close();

  // copy javascripts & Cascading Style Sheets to report the directory
  char cmd[1024];
  sprintf(cmd,"cp %s/html/etc/html/ROOT.css %s/",gSystem->ExpandPathName("$HOME_WAT"),odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp %s/html/etc/html/ROOT.js %s/",gSystem->ExpandPathName("$HOME_WAT"),odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp %s/html/etc/html/tabber.css %s/",gSystem->ExpandPathName("$HOME_WAT"),odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp %s/html/etc/html/tabber.js %s/",gSystem->ExpandPathName("$HOME_WAT"),odir.Data());
  gSystem->Exec(cmd);

  ModifyFontSizeCSS(odir+"/tabber.css");
}

void ModifyFontSizeCSS(TString tabber) {

  ifstream in;
  in.open(tabber.Data(),ios::in);
  if (!in.good()) {cout << "MakeHtmlPE::ModifyFontSizeCSS - Error Opening File : " << tabber.Data() << endl;exit(1);}

  ofstream out;
  TString tabber_tmp = tabber+".tmp";
  out.open(tabber_tmp,ios::out);
  if (!out.good()) {cout << "MakeHtmlPE::ModifyFontSizeCSS - Error Opening File : " << tabber_tmp << endl;exit(1);}

  char str[1024];
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    TString ostr = str;
    ostr.ReplaceAll("0.8em","0.75em");
    out << ostr.Data() << endl; 
  }
  out.close();
  in.close();

  char cmd[1024];
  sprintf(cmd,"mv %s %s",tabber_tmp.Data(),tabber.Data());
  //cout << cmd << endl;
  gSystem->Exec(cmd);
}

int ReadConfig(TString ifconfig, vector<TString> &gwname, vector<TString> &wdir, vector<TString> &rdir, vector<TString> &sdir) {

  ifstream in;
  in.open(ifconfig.Data(),ios::in);
  if (!in.good()) {cout << "MakeHtmlPE::ReadConfig - Error Opening File : " << ifconfig.Data() << endl;exit(1);}

  bool start=false;
  char log[1024];
  while(true) {
    in.getline(log,1024);
    if (!in.good()) break;
    if(log[0]=='#') continue;
    TString line = log;
    if(line.Contains("ifeq") && line.Contains("$(GW_NAME)")) {
      start=true;
      TString gw = line(line.Index(',')+1, line.Last(')')-line.Index(',')-1);
      gwname.push_back(gw);
      //cout << gw << endl;
    }
    if(line.Contains("endif") && start==true) start=false;
    if(start) {
      if(line.Contains("PE_OFFPATH")) {
        TString dir = line(line.Index('/')+1, line.Sizeof()-line.Index(',')-1);
        dir = gSystem->BaseName(dir);
        wdir.push_back(dir);
      }
      if(line.Contains("PE_RDIR")) {
        TString dir = line(line.Index('=')+1, line.Sizeof()-line.Index(',')-2);
        dir.ReplaceAll(" ","");
        rdir.push_back(dir);
        //cout << dir << endl;
      }
      if(line.Contains("PE_ONPATH")) {
        TString dir = line(line.Index('/')+1, line.Sizeof()-line.Index(',')-1);
        dir = gSystem->BaseName(dir);
        sdir.push_back(dir);
      }
    }
  }

  in.close();

  if((wdir.size() != gwname.size())) {
    cout << "MakeHtmlPE::ReadConfig - missing PE_OFFPATH config declarations: MISSED = " 
         << gwname.size()-wdir.size() << " in file " << ifconfig << endl;exit(1);
  }
  if((rdir.size() != gwname.size())) {
    cout << "MakeHtmlPE::ReadConfig - missing PE_RDIR config declarations: MISSED = " 
         << gwname.size()-rdir.size() << " in file " << ifconfig << endl;exit(1);
  }
  if((sdir.size() != gwname.size())) {
    cout << "MakeHtmlPE::ReadConfig - missing PE_ONPATH config declarations: MISSED = " 
         << gwname.size()-sdir.size() << " in file " << ifconfig << endl;exit(1);
  }

  for(int i=0;i<gwname.size();i++) cout << i << "\t" << gwname[i] << endl;
  for(int i=0;i<gwname.size();i++) cout << i << "\t" << wdir[i] << endl;
  for(int i=0;i<gwname.size();i++) cout << i << "\t" << rdir[i] << endl;
  for(int i=0;i<gwname.size();i++) cout << i << "\t" << sdir[i] << endl;

  return gwname.size();
}

