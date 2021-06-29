/*
# Copyright (C) 2019 Gabriele Vedovato, Francesco Salemi
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

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  // get user defined html index file name (used by cwb_mkrep command)
  TString cwb_mkrep_index_dir="";
  TString cwb_mkrep_index_file="";
  if(gSystem->Getenv("CWB_MKREP_INDEX_FILE")!=NULL) {
    cwb_mkrep_index_file=TString(gSystem->Getenv("CWB_MKREP_INDEX_FILE"));
    if(!cwb_mkrep_index_file.EndsWith(".html")) {
      cout<<"cwb_mkhtml_index.C : Error : the index html must have the .html extension"<<endl;
      exit(1);
    }
    // get and create dir name contained in cwb_mkrep_index_file path
    cwb_mkrep_index_dir = cwb_mkrep_index_file;
    cwb_mkrep_index_dir = cwb_mkrep_index_dir.ReplaceAll(gSystem->BaseName(cwb_mkrep_index_dir),"");
    char cmd[1024];
    sprintf(cmd,"mkdir -p %s",cwb_mkrep_index_dir.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
  } else {
    cwb_mkrep_index_file="";
  }

  // If fad is enable & simulation==4 then execute cwb_mkfad.C macro
  if(pp_fad!="" && simulation==4) {
    gROOT->LoadMacro("$HOME_WAT/tools/cwb/macros/cwb_mkfad.C");
    // get loud parameter from pp_fad 
    TString pp_fad_nzbins = CWB::Toolbox::getParameter(pp_fad,"--nzbins");
    if(pp_fad_nzbins=="") pp_fad_nzbins="0";
    int ipp_fad_nzbins = pp_fad_nzbins.Atoi();
    // execute cwb_mkfad.C macro
    cwb_mkfad("fad",ipp_fad_nzbins,0,false); 
    // add FAD report to the main report
    for(int i=0;i<pp_sreport.size();i++) {
      if(pp_sreport[i]=="") {pp_sreport[i]="--link fad --label FAD --name fad.html --high 4800";break;}
    }
  } 

  // check if there are auxiliary reports to be added to the main report
  int pp_sreport_size=0;
  for(int i=0;i<pp_sreport.size();i++) if(pp_sreport[i]!="") pp_sreport_size++;
  if((pp_sreport_size>0)&&(simulation>=0)) {
    vector<TString> pp_sreport_tmp = pp_sreport;
    for(int i=0;i<pp_sreport.size();i++) pp_sreport[i]="";
    pp_sreport[0]=TString("--label ")+data_label;
    pp_sreport_size=1;
    for(int i=0;i<pp_sreport_tmp.size();i++) if(pp_sreport_tmp[i]!="") {
      pp_sreport[pp_sreport_size]=pp_sreport_tmp[i];
      pp_sreport_size++;
    }
  }

  // open output index file
  ofstream out;
  char ofile[256];
  if(cwb_mkrep_index_file!="") {
    sprintf(ofile,"%s", cwb_mkrep_index_file.Data());
  } else {
    sprintf(ofile,"%s/index.html", pp_dir);
  }
  cout << "make index html file : " << ofile << endl;
  out.open(ofile,ios::out);
  if (!out.good()) {cout << "cwb_mkhtml_index.C : Error Opening File : " << ofile << endl;exit(1);}

  // open input index template file
  ifstream in;
  in.open(html_index_template,ios::in);
  if (!in.good()) {
    cout << "cwb_mkhtml_index.C : Error Opening File : " << html_index_template << endl;
    exit(1);
  }

  char sbody_height[256]; 
  if(simulation==0) {
    sprintf(sbody_height,"%d",2350+26*pp_max_nloudest_list);
  } else {
    if(TString(gSystem->Getenv("CWB_REPORT_PE"))!="") {
      sprintf(sbody_height,"%d",5000);
    } else {
      sprintf(sbody_height,"%d",5800);
    }
  }

  char istring[1024];
  while (1) {
    in.getline(istring,1024);
    if (!in.good()) break;
    TString ostring(istring);
    out << ostring.Data() << endl;
  }

  out << "<html>" << endl;
  out << "<br>" << endl;

  if(pp_sreport_size==0) {

    out << "  <iframe src=\"header.html\" width=\"100%\" height=\"970px\" " 
        << "marginwidth=\"15\" marginheight=\"15\" frameborder=\"0\"></iframe>" << endl;

    if(TString(gSystem->Getenv("CWB_REPORT_PE"))!="") {

      out << "<div class=\"tabber\">" << endl;

      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>Sky Localization</h2>" << endl;
      out << "  <iframe src=\"prc_body.html\" width=\"100%\" " 
          << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
      out << "</div>" << endl;

      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>Waveform</h2>" << endl;
      out << "  <iframe src=\"wrc_body.html\" width=\"100%\" " 
          << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
      out << "</div>" << endl;

      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>Chirp Mass</h2>" << endl;
      out << "  <iframe src=\"mchirp_body.html\" width=\"100%\" " 
          << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
      out << "</div>" << endl;

      out << "</div>" << endl;
    } else if(TString(gSystem->Getenv("CWB_REPORT_CBC"))!="") {

	// if CWB_DOC_URL is define then man infos are added to web pages
	  TString cwb_doc_url="";
	  if(gSystem->Getenv("CWB_DOC_URL")!=NULL) {
	    cwb_doc_url=TString(gSystem->Getenv("CWB_DOC_URL"));
	  }

	  char plot_list[1024];
	  sprintf(plot_list,"<a target=\"_parent\" name=\"Full Plot List\" href=\"%s\">Full Plot List</a>",pp_data_dir);
	  if(cwb_doc_url!="") out<<"<table align=\"center\" width=97%> <tr> <td align=\"left\">"<<endl;
	  out << plot_list << endl;
	  if(cwb_doc_url!="") {
	    out<<"</td>"<<endl;
	    out<<"<td align=\"right\">"<<endl;
	    out<<"<a target=\"_parent\" href=\""<<cwb_doc_url.Data()
	       <<"/cwb/man/Simulation-part.html#Simulation-part\">infos</a>"<<endl;
	    out<<"</td> </tr> </table>"<<endl;
	  }

      out << "<div class=\"tabber\">" << endl;

      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>Main plots</h2>" << endl;
      out << "  <iframe src=\"main_body.html\" width=\"100%\" " 
          << " height=\"" << 1900 << "px\" frameborder=\"0\"></iframe>" << endl;
      out << "</div>" << endl;

      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>ROC curves</h2>" << endl;
      out << "  <iframe src=\"ROC_body.html\" width=\"100%\" " 
          << " height=\"" << 3500 << "px\" frameborder=\"0\"></iframe>" << endl;
      out << "</div>" << endl;

      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>Sampled parameter space</h2>" << endl;
      out << "  <iframe src=\"parspace_body.html\" width=\"100%\" " 
          << " height=\"" << 1800 << "px\" frameborder=\"0\"></iframe>" << endl;
      out << "</div>" << endl;

      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>Distance plots</h2>" << endl;
      out << "  <iframe src=\"distance_body.html\" width=\"100%\" " 
          << " height=\"" << 2180 << "px\" frameborder=\"0\"></iframe>" << endl;
      out << "</div>" << endl;

      out << "<div class=\"tabbertab\">" << endl;
      out << "  <h2>SNR plots</h2>" << endl;
      out << "  <iframe src=\"snr_body.html\" width=\"100%\" " 
          << " height=\"" << 1300 << "px\" frameborder=\"0\"></iframe>" << endl;
      out << "</div>" << endl;

      out << "</div>" << endl;


    } else {
      out << "  <iframe src=\"body.html\" width=\"100%\" " 
          << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    }
  } else {

    sprintf(sbody_height,"%d",5800);
    out << "<div class=\"tabber\">" << endl;
    for(int i=0;i<pp_sreport.size();i++) if(pp_sreport[i]!="") {

      TString pp_sreport_link = CWB::Toolbox::getParameter(pp_sreport[i],"--link");
      if(pp_sreport_link=="" && i!=0) {
        cout<<"cwb_mkhtml_index.C : Error : pp_sreport --link not defined"<<endl;exit(1);}
      if((i!=0)||(simulation<0)) {
        if(!pp_sreport_link.EndsWith("/")) pp_sreport_link=pp_sreport_link+"/";
      }

      TString pp_sreport_label = CWB::Toolbox::getParameter(pp_sreport[i],"--label");
      if((pp_sreport_link!="</tab>/")&&(pp_sreport_label=="")) {
        cout<<"cwb_mkhtml_index.C : Error : pp_sreport --label not defined"<<endl;exit(1);}

      TString pp_sreport_high = CWB::Toolbox::getParameter(pp_sreport[i],"--high");
      if(pp_sreport_high=="") pp_sreport_high=sbody_height;
      int ipp_sreport_high = pp_sreport_high.Atoi();

      TString pp_sreport_name = CWB::Toolbox::getParameter(pp_sreport[i],"--name");

      if(pp_sreport_link=="<tab>/") {		// open a sub tab
        out << "<div class=\"tabbertab\">" << endl;
        out << "  <h2>" << pp_sreport_label << "</h2>" << endl;
        out << "<div class=\"tabber\">" << endl;
      } else if(pp_sreport_link=="</tab>/") {	// close sub tab
        out << "</div>" << endl;
        out << "</div>" << endl;
      } else {					// add a tab
        out << "<div class=\"tabbertab\">" << endl;
        out << "  <h2>" << pp_sreport_label << "</h2>" << endl;

        if(pp_sreport_name=="") {
          out << "  <iframe src=\"" << pp_sreport_link << "header.html\" width=\"100%\" height=\"900px\" " 
              << "marginwidth=\"15\" marginheight=\"15\" frameborder=\"0\"></iframe>" << endl;

          out << "  <iframe src=\"" << pp_sreport_link << "body.html\" width=\"100%\" " 
              << " height=\"" << ipp_sreport_high << "px\" frameborder=\"0\"></iframe>" << endl;
        } else {
          out << "  <iframe src=\"" << pp_sreport_link << pp_sreport_name << "\" width=\"100%\" " 
              << " height=\"" << ipp_sreport_high << "px\" frameborder=\"0\"></iframe>" << endl;
        }

        out << "</div>" << endl;
      }
    }
    out << "</div>" << endl;
  }

  out << "</html>" << endl;

  in.close();
  out.close();

  // copy javascripts & Cascading Style Sheets to report the directory
  TString odir = (cwb_mkrep_index_file!="") ? cwb_mkrep_index_dir : pp_dir;
  char cmd[1024];
  sprintf(cmd,"cp %s/html/etc/html/ROOT.css %s/",gSystem->ExpandPathName("$HOME_WAT"),odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp %s/html/etc/html/ROOT.js %s/",gSystem->ExpandPathName("$HOME_WAT"),odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp %s/html/etc/html/tabber.css %s/",gSystem->ExpandPathName("$HOME_WAT"),odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp %s/html/etc/html/tabber.js %s/",gSystem->ExpandPathName("$HOME_WAT"),odir.Data());
  gSystem->Exec(cmd);

  exit(0);
}
