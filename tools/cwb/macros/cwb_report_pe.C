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


// post-production macro for the simulation report : used by the cwb_report command
{

  #define NMDC_MAX 64

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

#ifndef _USE_ROOT6
  // the CWB_CAT_NAME declared in CWB_EPPARAMETERS_FILE is not visible. why?
  // the include is defined again  
  #undef GTOOLBOX_HH
#else
  #include "GToolbox.hh"
#endif

  if(nIFO==1) { // Single Detector Mode
    CWB::config config;
    config.Import();
    config.SetSingleDetectorMode();
    config.Export();
  }

  // get merge label
  TString mlabel=TString(gSystem->Getenv("CWB_MERGE_LABEL"));

  // initialize the ifo name (takes into account the non built-in detectors)
  vector<TString> IFO;
  for(int n=0;n<nIFO;n++) IFO.push_back(ifo[n]);
 
  // create output report dir
  CWB::Toolbox::mkDir(netdir,!pp_batch);

  TString pp_pe_mdc = CWB::Toolbox::getParameter(pp_pe,"--mdc");
  // check if format is correct 
  pp_pe_mdc.ReplaceAll("/",""); 
  if(!pp_pe_mdc.IsDigit() && pp_pe_mdc!="-1") {
    cout << "cwb_report_pe.C : Error - bad option '--mdc type' in pp_pe parameter" << endl;
    cout << "                   mdc type must be a integer number [-1:N] or i/j/.../k=[1:N]" << endl;
    cout << endl;
    exit(1);
  }
  // extract list of mdc types
  vector<int> mtype;
  TObjArray* token = TString(pp_pe_mdc).Tokenize(TString("/"));
  if(token) {
    for(int i=0;i<token->GetEntries();i++) {
      int itype = ((TObjString*)token->At(i))->GetString().Atoi();
      mtype.push_back(itype);
    }
    delete token;
  }
  if(mtype.size()>1) {	// multiple mdc types
    for(int i=0;i<mtype.size();i++) {
      if(mtype[i]==-1 || mtype[i]==0) {
        cout << "cwb_report_pe.C : Error - bad option '--mdc type' in pp_pe parameter" << endl;
        cout << "                   multiple mdc types must be in the range [1:N]" << endl;
        cout << endl;
        exit(1);
      }
    }
  } else {
    if(mtype[0]<-1) {
      cout << "cwb_report_pe.C : Error - bad option '--mdc type' in pp_pe parameter" << endl;
      cout << "                   mdc type must be in the range [-1:N]        " << endl;
      cout << endl;
      exit(1);
    }
  }

  vector<int> vtype;
  vector<TString> vname;
  if(mtype[0]!=0) { 	// read injection file types
    char    imdc_set[NMDC_MAX][128];    // injection set
    size_t  imdc_type[NMDC_MAX];        // injection type
    char    imdc_name[NMDC_MAX][128];   // injection name
    double  imdc_fcentral[NMDC_MAX];    // injection central frequencies
    double  imdc_fbandwidth[NMDC_MAX];  // injection bandwidth frequencies
    size_t  imdc_index[NMDC_MAX];       // type reference array
    size_t  imdc_iset[NMDC_MAX];        // injection set index

    int ninj=ReadInjType(mdc_inj_file,NMDC_MAX,imdc_set,imdc_type,
                         imdc_name,imdc_fcentral,imdc_fbandwidth);
    if(ninj==0) {
      cout << "cwb_report_pe.C : Error - no injection - terminated" << endl;
      exit(1);
    }
    if(mtype[0]>0) { 	// mdc_type>0  : select types from injection file types
      for(int i=0;i<mtype.size();i++) {
        for(int j=0;j<ninj;j++) {
          if(mtype[i]==(imdc_type[j]+1)) {
            vtype.push_back(imdc_type[j]+1);
            vname.push_back(imdc_name[j]);
          }
        } 
      }
    } else {		// mdc_type=-1 : use full list of types from injection file types
      for(int j=0;j<ninj;j++) {
        vtype.push_back(imdc_type[j]+1);
        vname.push_back(imdc_name[j]);
      }
    }
  } else {		// mdc_type=0  : all types are displayed together
    vtype.push_back(0);
    vname.push_back("ALL");
  }

  Color_t colors[16] = {2, 3, 6, 4, 8, 43, 7, 8, 4, 5, 2, 43, 1, 3, 2, 1};

  // create setup file for the DrawMedianPRCvsSNR.C macro
  ofstream out;
  char listFile[256];
  sprintf(listFile,"%s/DrawMedianPRCvsSNR.lst", netdir);
  cout << listFile << endl;
  out.open(listFile,ios::out);
  if(!out.good()) {cout << "cwb_report_pe.C : Error Opening File : " << listFile << endl;exit(1);}
  for(int i=0;i<vtype.size();i++) {
    char wave_file_name[1024];
    sprintf(wave_file_name,"%s/wave_%s.%s.root",merge_dir,data_label,mlabel.Data());
    out << wave_file_name << "\t" << vname[i] << "\t" << vtype[i] 
        << "\t" << (int)colors[i%16] << "\t" << 0 << endl;
  }
  out.close();

  // load PRC macros
  gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawWRC.C"));
  gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawRECvsINJ.C"));
  gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawSkyDistributionPRC.C"));
  gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawSearchAreaPRC.C"));
  gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawCosOmegaPRC.C"));
  gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawCoverageVsPercentagePRC.C"));
  gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawMedianPRCvsSNR.C"));

  // get plugin options (polarization=SCALAR/TENSOR)

  TString polarization="TENSOR";
  if(TString(parPlugin)!="") {
    cout << "pe polarization options : " << parPlugin << endl;
    TObjArray* token = TString(parPlugin).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++) {

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();

      if(stok.Contains("--pe_polarization=")) {
        polarization=stok;
        polarization.Remove(0,polarization.Last('=')+1);
        polarization.ToUpper();;
      }
    }
  }

  // execute PRC macros

  DrawMedianPRCvsSNR(listFile,"","MEDIAN50",TString(netdir)+"/median50_vs_snr.png", 
                     T_win, pp_inetcc, T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);
  DrawMedianPRCvsSNR(listFile,"","MEDIAN90",TString(netdir)+"/median90_vs_snr.png", 
                     T_win, pp_inetcc, T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawSkyDistributionPRC(data_label, netdir, mlabel, IFO, detParms, T_win, pp_inetcc, 
                         T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar, true, polarization, true);
  DrawSkyDistributionPRC(data_label, netdir, mlabel, IFO, detParms, T_win, pp_inetcc, 
                         T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar, false, polarization);

  DrawWRC("nre",data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
          T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawWRC("ff",data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
          T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawWRC("of",data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
          T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawWRC("ofnre",data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
          T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawWRC("offf",data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
          T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawWRC("snr",data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
          T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  for(int n=1;n<13;n++) {
    char gtype[16];sprintf(gtype,"mch%d",n);
    DrawWRC(gtype,data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
            T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);
  }

  DrawRECvsINJ("snr", data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
               T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawRECvsINJ("distance", data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
               T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawCosOmegaPRC(data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
                  T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawSearchAreaPRC(data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
                    T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawCoverageVsPercentagePRC("prc", data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
                              T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawCoverageVsPercentagePRC("nre", data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
                              T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  DrawCoverageVsPercentagePRC("fre", data_label, netdir, mlabel, nIFO, T_win, pp_inetcc, 
                              T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  exit(0);
}
