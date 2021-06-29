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


// produce skymap html report : use -> https://github.com/reedessick/skymap_statistics

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  TString home_skymap_lib;
  if(gSystem->Getenv("HOME_SKYMAP_LIB")==NULL) {
    cout << "cwb_report_skymap.C - Error : environment HOME_SKYMAP_LIB is not defined!!!" << endl;exit(1);
  } else {
    home_skymap_lib=TString(gSystem->Getenv("HOME_SKYMAP_LIB"));
  }

  TString cwb_user_url;
  if(gSystem->Getenv("CWB_USER_URL")==NULL) {
    cout << "cwb_report_skymap.C - Error : environment CWB_USER_URL is not defined!!!" << endl;exit(1);
  } else {
    cwb_user_url=TString(gSystem->Getenv("CWB_USER_URL"));
  }

  TString cwb_rep_url="";
  if(gSystem->Getenv("CWB_REP_URL")!=NULL) {
    cwb_rep_url=TString(gSystem->Getenv("CWB_REP_URL"));
  }

  TString cwb_skymap_file;
  if(gSystem->Getenv("CWB_SKYMAP_FILE")==NULL) {
    cout << "cwb_report_skymap.C - Error : environment CWB_SKYMAP_FILE is not defined!!!" << endl;exit(1);
  } else {
    cwb_skymap_file=TString(gSystem->Getenv("CWB_SKYMAP_FILE"));
  }

  // check skymap extension
  if(cwb_skymap_file.EndsWith(".fits")) {
  } else if(cwb_skymap_file.EndsWith(".fits.gz")) {
  } else {
    cout << "cwb_report_skymap.C - Error : file name must ends with .fits/.fits.gz !!!" << endl;exit(1);
  }

  TB.checkFile(cwb_skymap_file);

  TString cwb_skymap_dir = cwb_skymap_file;	// default output dir

  TString cwb_report_options = TString(gSystem->Getenv("CWB_REPORT_OPTIONS"));

  // check if options contains the output dir & if contains multiple skymaps
  bool singleSkymap = true;
  if(cwb_report_options!="") {
    cout << cwb_report_options << endl;

    TObjArray* token = TString(cwb_report_options).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++){

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();

      if(stok=="-o") {
        if(j<token->GetEntries()-1) {
          TObjString* otoken = (TObjString*)token->At(j+1);
          cwb_skymap_dir = otoken->GetString();
        }
      }
      if(stok=="--output-dir") {
        if(j<token->GetEntries()-1) {
          TObjString* otoken = (TObjString*)token->At(j+1);
          cwb_skymap_dir = otoken->GetString();
        }
      }

      if(stok=="-O") {
        if(j<token->GetEntries()-1) {
          TObjString* otoken = (TObjString*)token->At(j+1);
          cwb_rep_url = otoken->GetString();
        }
      }
      if(stok=="--output-url") {
        if(j<token->GetEntries()-1) {
          TObjString* otoken = (TObjString*)token->At(j+1);
          cwb_rep_url = otoken->GetString();
        }
      }
      if(stok=="--fits") {
        if(j<token->GetEntries()-1) {
          TObjString* otoken = (TObjString*)token->At(j+1);
          TString cwb_skymap_comp = otoken->GetString();
          // check skymap extension
          if(cwb_skymap_comp.EndsWith(".fits")) {
          } else if(cwb_skymap_comp.EndsWith(".fits.gz")) {
          } else {
            cout << "cwb_report_skymap.C - Error : file name of comparison skymap must ends with .fits/.fits.gz !!!" << endl;exit(1);
          }
          // remove --fits string
          cwb_report_options.ReplaceAll("--fits",""); 
          singleSkymap = false;
        }
      }
    }
  }

  if(cwb_rep_url=="") {
    cout << "cwb_report_skymap.C - Error : --output-url not defined !!!" << endl;exit(1);
  }

  TString cwb_skymap_name = gSystem->BaseName(cwb_skymap_file);
  // for multiple skymaps add 'c' in front of of fits file directory
  if(!singleSkymap) cwb_skymap_dir.ReplaceAll(cwb_skymap_name,TString("c")+cwb_skymap_name);
  cwb_skymap_name.ReplaceAll(".gz","");
  cwb_skymap_name.ReplaceAll(".fits","");

  if(cwb_skymap_dir.EndsWith(".fits.gz")) cwb_skymap_dir.ReplaceAll(".fits.gz","");
  if(cwb_skymap_dir.EndsWith(".fits")) cwb_skymap_dir.ReplaceAll(".fits","");

  // creates dir for skymap report
  bool overwrite = TB.checkFile(cwb_skymap_dir,true,"skymap statistic directory already exist");
  if(!overwrite) {cout << "cwb_report_skymap.C terminated !!!" << endl<<endl;gSystem->Exit(1);}

  TB.mkDir(cwb_skymap_dir,false,true);	// remove old directory

  char ifostr[64]="";
  for(int n=0; n<nIFO; n++) sprintf(ifostr,"%s -i %c",ifostr,ifo[n][0]);
  cout << "Network : " << ifostr << endl;

  unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files

  // create temporary skymap script file
  char skymap_script_file[1024];
  sprintf(skymap_script_file,"%s/skymap_%d.csh",tmp_dir,Pid);

  TString rel_skymap_dir = cwb_skymap_dir;
  // strip report_dir from rel_skymap_dir
  if(rel_skymap_dir.BeginsWith(report_dir)) if(rel_skymap_dir.First("/")>0) rel_skymap_dir.Remove(0,rel_skymap_dir.First("/")+1);
  char ourl[1024];
  sprintf(ourl,"%s/%s/%s",cwb_user_url.Data(),data_label,rel_skymap_dir.Data());
  TString odir = cwb_skymap_dir;

  TString pyCmd = singleSkymap ? "snglFITShtml.py" : "multFITShtml.py";

  // write output skymap script file
  ofstream out;
  out.open(skymap_script_file,ios::out);
  char ostring[1024];

  sprintf(ostring,"bash -c \\");
  out << ostring << endl;
  sprintf(ostring,"' \\");
  out << ostring << endl;
  sprintf(ostring,"USER_DIR=${PWD}; \\");
  out << ostring << endl;
  sprintf(ostring,"cd %s; \\",home_skymap_lib.Data());
  out << ostring << endl;
  sprintf(ostring,". setup.sh; \\");
  out << ostring << endl;
  sprintf(ostring,"cd $USER_DIR; \\");
  out << ostring << endl;
  sprintf(ostring,"fits=%s; \\",cwb_skymap_file.Data());
  out << ostring << endl;
  sprintf(ostring,"time %s -v %s --dT-nside 128 --dT-Nsamp 500 -o %s -t cwb -O %s --no-margticks $fits %s \\",
                  pyCmd.Data(),ifostr,odir.Data(),ourl,cwb_report_options.Data());
  out << ostring << endl;
  sprintf(ostring,"'");
  out << ostring << endl;

  out.close();

  // execute skymap script file -> generate skymap statistics + skyprobcc/skyprobcc-skymapSummary_cwb.html
  char cmd[1024];
  sprintf(cmd,". %s",skymap_script_file);
  cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"rm %s",skymap_script_file);
  gSystem->Exec(cmd);

  TString cwb_skymap_parent = gSystem->DirName(cwb_skymap_dir); 

  // if symbolic link exist -> copy link to index.html
  TString symLink;
  char index_html_file[1024];
  Long_t id,size,flags,mt;
  sprintf(index_html_file,"%s/index.html",cwb_skymap_parent.Data());
  symLink = CWB::Toolbox::getFileName(index_html_file);
  if(symLink!="") {
    int estat = gSystem->GetPathInfo(symLink.Data(),&id,&size,&flags,&mt);
    if(estat==0) {  
      sprintf(cmd,"rm %s;cp %s %s",index_html_file,symLink.Data(),index_html_file);
      cout << cmd << endl;
      gSystem->Exec(cmd);
    }
  } 

  int estat = gSystem->GetPathInfo(index_html_file,&id,&size,&flags,&mt);
  if(estat!=0) gSystem->Exit(0);

  // add link "skyprobcc/skyprobcc-skymapSummary_cwb.html" to index.html file

  TString cwb_skymap_html = cwb_skymap_name+"/"+cwb_skymap_name+"-skymapSummary_cwb.html";
  if(!singleSkymap) cwb_skymap_html = TString("c")+cwb_skymap_name+"/"+"multFITS-skymapComparison_cwb.html";

  ifstream in;
  in.open(index_html_file,ios::in);
  if(!in.good()) {cout << "cwb_report_skymap.C - Error Opening File : " << index_html_file << endl;exit(1);}

  ofstream out2;
  char index_html_file_tmp[1024];
  sprintf(index_html_file_tmp,"%s.tmp", index_html_file);
  cout << index_html_file_tmp << endl;
  out2.open(index_html_file_tmp,ios::out);
  if (!out2.good()) {cout << "cwb_report_skymap.C - Error Opening File : " << index_html_file_tmp << endl;exit(1);}

  bool found=false;
  char istr[1024];
  while(1) {
    in.getline(istr,1024);
    if (!in.good()) break;
    TString ostr(istr);
    if(ostr.Contains(cwb_skymap_html)) found=true;
    if(!found && ostr.Contains("cwb_parameters.C.html")) {
      out2 << "<li> <a href=\"" << cwb_skymap_html << "\" target=\"_blank\">Skymap Statistics</a>";
           if(cwb_skymap_html=="skyprobcc/skyprobcc-skymapSummary_cwb.html")    out2 << " ( point estimate )" << endl; 
      else if(cwb_skymap_html=="mskyprobcc/mskyprobcc-skymapSummary_cwb.html")  out2 << " ( median )" << endl; 
      else if(cwb_skymap_html=="cskyprobcc/multFITS-skymapComparison_cwb.html") out2 << " ( comparison )" << endl; 
      else out2 << endl; 
      found=false;
    }
    out2 << ostr.Data() << endl;
  }
  in.close();
  out2.close();

  // copy index_html_file_tmp -> index_html_file
  sprintf(cmd,"rm %s;mv %s %s",index_html_file,index_html_file_tmp,index_html_file);
  cout << cmd << endl;
  gSystem->Exec(cmd);

  gSystem->Exit(0);
}
