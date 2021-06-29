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

//
// This macro implements the combine statistic for the IMBHB and BBH searches
// See CWB::Toolbox::CombineCBC for a description of the algorithm
//

// WARNING: the combine procedure is based:
//          MDC        -> on comparison of the injected times time[nIFO]
//                        the time[nIFO] could be different for BBH and IMBHB 
//                        because are the weighted of whitened waveforms.
//                        The estimation could be different for BBH and IMBHB
//          WAVE       -> on comparison of time[0]
//                        The estimation could be different for BBH and IMBHB
//          The time uncertainty is defined in CWB::mdc::CombineCBC by the parameter
//                        #define TIME_UNCERTAINTY    0.1 // sec

{
  CWB::Toolbox TB;

  	  cwb_merge_label     = TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  TString cwb_combine_options = TString(gSystem->Getenv("CWB_COMBINE_OPTIONS"));

  TString mdir = merge_dir;

  // check if cwb_combine_options contains '--' than cwb_combine_options is used to extract all combine parameters
  TString cwb_combine_sfthr     = "";
  TString cwb_combine_sifarthr  = "";
  TString cwb_combine_search    = "";
  TString cwb_combine_wdir      = "";
  TString cwb_combine_mlabel    = "";
  TString cwb_combine_ulabel    = "";
  TString cwb_combine_check     = "";
  int cwb_combine_lag  = 0;
  int cwb_combine_slag = 0;
  if(cwb_combine_options.Contains("--")) {
    TString option="";
    // get the frequency threshold used to combine IMBHB & BBH searches
    cwb_combine_sfthr  = TB.getParameter(cwb_combine_options,"--fthr");
    // get the ifar threshold used to select zero lag events
    cwb_combine_sifarthr  = TB.getParameter(cwb_combine_options,"--ifarthr");
    // get the search type (IMBHB or BBH) to be combined 
    cwb_combine_search = TB.getParameter(cwb_combine_options,"--search");
    // get the working dir of the search to be combined
    cwb_combine_wdir = TB.getParameter(cwb_combine_options,"--wdir");
    // get the merge label of the search to be combined
    cwb_combine_mlabel = TB.getParameter(cwb_combine_options,"--mlabel");
    // get the user label to be used to tag the combined output wave/mdc files
    cwb_combine_ulabel = TB.getParameter(cwb_combine_options,"--ulabel");
    // get the lag used to select background data to be combined
    TString lag_str = TB.getParameter(cwb_combine_options,"--lag");
    if(lag_str.IsDigit()) cwb_combine_lag=lag_str.Atoi();
    // get the slag used to select background data to be combined
    TString slag_str = TB.getParameter(cwb_combine_options,"--slag");
    if(slag_str.IsDigit()) cwb_combine_slag=slag_str.Atoi();
    // get check option (background -> printout zero lag livetimes)
    cwb_combine_check = TB.getParameter(cwb_combine_options,"--check");
    cwb_combine_check.ToUpper();
    if(cwb_combine_check=="") cwb_combine_check  = "FALSE";
  }

  cout << endl;
  // check if the check option is correct
  if(cwb_combine_check!="FALSE" && cwb_combine_check!="TRUE") {
    cout << endl << "cwb_combine_cbc : Error - check available values are : true/false!!!" << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_combine_lag and cwb_combine_slag are>0
  if(cwb_combine_lag<0 || cwb_combine_slag<0) {
    cout << "cwb_combine_cbc.C : Error - lag,slag must be >= 0 : --lag " 
         << cwb_combine_lag << " --slag " << cwb_combine_slag << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_combine_fthr is an integer number
  int cwb_combine_fthr = 0;
  if(cwb_combine_sfthr.IsDigit()) {
    cwb_combine_fthr = cwb_combine_sfthr.Atoi();
  } else {
    cout << "cwb_combine_cbc.C : Error - fthr is not an interger number: --fthr " << cwb_combine_sfthr << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_combine_ifarthr is an integer number
  float cwb_combine_ifarthr = 0.;
  if(cwb_combine_sifarthr.IsFloat()) {
    cwb_combine_ifarthr = cwb_combine_sifarthr.Atof();
    if(cwb_combine_ifarthr<0) {
      cout << "cwb_combine_cbc.C : Error - ifarthr must be>=0: --ifarthr " << cwb_combine_sifarthr << endl << endl;
      gSystem->Exit(1);
    } 
  } else if(cwb_combine_sifarthr!="") {
    cout << "cwb_combine_cbc.C : Error - ifarthr is not a float number: --ifarthr " << cwb_combine_sifarthr << endl << endl;
    gSystem->Exit(1);
  }
  // simulation -> check if cwb_merge_label contains the unique events tag 
  if(simulation && !cwb_merge_label.Contains(".C_U")) {
    cout << "cwb_combine_cbc.C : Error - _merge label not contains the unique events tag '.C_U': merge label = " << cwb_merge_label << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_combine_search is defined
  if(cwb_combine_search!="IMBHB" && cwb_combine_search!="BBH") {
    cout << "cwb_combine_cbc.C : Error - search not defined, valid values are: IMBHB,BBH" << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_combine_wdir is defined
  if(cwb_combine_wdir=="") {
    cout << "cwb_combine_cbc.C : Error - wdir not defined" << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_combine_mlabel is defined
  if(cwb_combine_mlabel=="") {
    cout << "cwb_combine_cbc.C : Error - mlabel not defined" << endl << endl;
    gSystem->Exit(1);
  }
  // simulation -> check if cwb_combine_mlabel contains the unique events tag 
  if(simulation && !cwb_combine_mlabel.Contains(".C_U")) {
    cout << "cwb_combine_cbc.C : Error - mlabel not contains the unique events tag '.C_U': --mlabel " << cwb_combine_mlabel << endl << endl;
    gSystem->Exit(1);
  }
  // check if all characters in cwb_combine_ulabel are alphanumeric
  if(cwb_combine_ulabel!="" && !cwb_combine_ulabel.IsAlnum()) {
    cout << "cwb_combine_cbc.C : Error - ulabel is not alphanumeric: --ulabel " << cwb_combine_ulabel << endl << endl;
    gSystem->Exit(1);
  }

  vector<TString> ifwave(2);
  // create input wave root combine file name
  TString ifile1 = TString::Format("%s/wave_%s.%s.root",merge_dir,data_label,cwb_merge_label.Data());
  TString ifile2 = TString::Format("%s/%s/wave_%s.%s.root",cwb_combine_wdir.Data(),merge_dir,
                                   gSystem->BaseName(cwb_combine_wdir),cwb_combine_mlabel.Data());
  ifwave[0] = (cwb_combine_search=="BBH") ? ifile1 : ifile2;
  ifwave[1] = (cwb_combine_search=="BBH") ? ifile2 : ifile1;
  // check if input wave files exist
  for(int n=0;n<2;n++) CWB::Toolbox::checkFile(ifwave[n]);

  vector<TString> ifmdc;
  vector<TString> iflive;
  if(simulation) {	// simulation
    ifmdc.resize(2);
    ifmdc[0] = ifwave[0]; ifmdc[0].ReplaceAll("wave_","mdc_");
    ifmdc[1] = ifwave[1]; ifmdc[1].ReplaceAll("wave_","mdc_");
    // check if input mdc files exist
    for(int n=0;n<2;n++) CWB::Toolbox::checkFile(ifmdc[n]);
  } 
  if(!simulation && cwb_combine_check=="TRUE") {		// background
    // check livetimes zero lag IMBHB vs BBH
    cout << "Check livetimes zero lag IMBHB vs BBH ..." << endl << endl;
    iflive.resize(2);
    iflive[0] = ifwave[0]; iflive[0].ReplaceAll("wave_","live_");
    iflive[1] = ifwave[1]; iflive[1].ReplaceAll("wave_","live_");
    // check if input live files exist
    for(int n=0;n<2;n++) CWB::Toolbox::checkFile(iflive[n]);
    // check if lag,slag livetime is the same for IMBHB and BBH
    TChain live1("liveTime");
    live1.Add(iflive[0]);
    TChain live2("liveTime");
    live2.Add(iflive[1]);
    double livetime1,livetime2;
    if(cwb_combine_lag==0 && cwb_combine_slag==0) {	// zero lag 
      livetime1=CWB::Toolbox::getZeroLiveTime(nIFO,live1);
      livetime2=CWB::Toolbox::getZeroLiveTime(nIFO,live2);
    } else {						// fake zero lag
      wavearray<double> Trun(500000); Trun = 0.;
      wavearray<double> Wlag[NIFO_MAX+1];
      wavearray<double> Wslag[NIFO_MAX+1];
      wavearray<double> Tlag;
      wavearray<double> Tdlag;
      livetime1=CWB::Toolbox::getLiveTime(nIFO,live1,Trun,Wlag,Wslag,Tlag,Tdlag,cwb_combine_lag,cwb_combine_slag);
      livetime2=CWB::Toolbox::getLiveTime(nIFO,live2,Trun,Wlag,Wslag,Tlag,Tdlag,cwb_combine_lag,cwb_combine_slag);
    }
    cout << "-----------------------------------------------------------------------------------" << endl;
    cout << "cwb_combine_cbc.C : zero lag livetime of IMBHB " 
         << livetime1 << " (sec) " << livetime1/(24*3600) << " (days) " << endl;
    cout << "cwb_combine_cbc.C : zero lag livetime of BBH   " 
         << livetime2 << " (sec) " << livetime2/(24*3600) << " (days) " << endl;
    cout << "-----------------------------------------------------------------------------------" << endl;
    cout << endl;
    exit(0);
  }

  // input files 
  cout << endl;
  cout << "Input IMBHB files : " << endl;
  cout << ifwave[0] << endl;
  if(simulation) cout << ifmdc[0] << endl;
  cout << endl;
  cout << "Input BBH files : " << endl;
  cout << ifwave[1] << endl;
  if(simulation) cout << ifmdc[1] << endl;
  cout << endl;

  // create output wave root combine file name
  TString ofwave = ifile1;
  TString cwb_combine_tag = TString::Format("J_%s_f%dHz",cwb_combine_search.Data(),cwb_combine_fthr);
  if(!simulation) cwb_combine_tag = cwb_combine_tag+TString::Format("_lag%d_slag%d",cwb_combine_lag,cwb_combine_slag);
  if(cwb_combine_ulabel!="") cwb_combine_tag = cwb_combine_tag+"_"+cwb_combine_ulabel;
  ofwave.ReplaceAll(".root",TString::Format(".%s.root",cwb_combine_tag.Data()));
  // create output mdc root combine file name
  TString ofmdc = ofwave;
  ofmdc.ReplaceAll("wave_","mdc_");

  // combine searches
  TString msearch = (cwb_combine_search=="BBH") ? "IMBHB" : "BBH";             // extract master search (used in history)
  TString infos = "("+cwb_combine_tag+") - ( "+cwb_combine_options+" )";       // infos to be stored into history
  int err = TB.CombineCBC(ifwave, ifmdc, ofwave, ofmdc, nIFO, cwb_combine_fthr, 
                          msearch, infos, cwb_combine_lag, cwb_combine_slag, cwb_combine_ifarthr);
  if(err) {
    cout << "cwb_combine_cbc.C : Warning - no events selected, create a symbolic link to wave file" << endl << endl;
    Long_t id,size,flags,mt;
    int estat = gSystem->GetPathInfo(ifile1,&id,&size,&flags,&mt);
    if(estat==0) {
      char cmd[1024]; sprintf(cmd,"ln -sf ../%s %s",ifile1.Data(),ofwave.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
    }
  }

  // create output merge list combine file name
  TString ofmerge = ofwave;
  ofmerge.ReplaceAll("wave_","merge_");
  ofmerge.ReplaceAll(".root",".lst");
  char cmd[1024]; sprintf(cmd,"touch %s",ofmerge.Data());
  gSystem->Exec(cmd);

  TString oflive;
  if(simulation==0) {	// background
    // WARNING: this is correct if the zeo lag livetime IMBHB = BBH
    // create symbolic link to live root file
    TString iflive = ifile1;
    iflive.ReplaceAll("wave_","live_");
    iflive.Remove(0,iflive.Last('/')+1);    // strip path
    oflive = ofwave;
    oflive.ReplaceAll("wave_","live_");
    oflive.Remove(0,oflive.Last('/')+1);    // strip path
    //cout << oflive << endl;
    Long_t id,size,flags,mt;
    int estat = gSystem->GetPathInfo(mdir+"/"+iflive,&id,&size,&flags,&mt);
    if(estat==0) {
      sprintf(cmd,"cd %s;ln -sf %s %s",mdir.Data(),iflive.Data(),oflive.Data());
      //cout << cmd << endl;
      gSystem->Exec(cmd);
    }
  }

  // output files 
  cout << endl;
  cout << "Output Combined files : " << endl;
  cout << ofwave << endl;
  if(simulation) cout << ofmdc << endl;
  else           cout << mdir+"/"+oflive << endl;
  cout << ofmerge << endl;
  cout << endl;

  exit(0);
}
