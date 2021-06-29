/*
# Copyright (C) 2019 Gabriele Vedovato, Marco Drago
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


// post-production macro for the background report : used by the cwb_report command

{
  cout<<"cwb_report_prod_1.C starts..."<<endl;

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
#endif
  #include "GToolbox.hh"

  CWB::Toolbox::mkDir(netdir,!pp_batch);

  gStyle->SetTitleOffset(1.4,"X");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelOffset(0.014,"X");
  gStyle->SetLabelOffset(0.010,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(kFALSE);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  cout<<"netplot starts..."<<endl;

  if(nIFO==1) { // Single Detector Mode
    CWB::config config;
    config.Import();
    config.SetSingleDetectorMode();
    config.Export();
  }

  // ----------------------------------------------------------
  // canvas
  // ----------------------------------------------------------
  TCanvas *c1 = new TCanvas("c","C",0,0,600,600);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetLogx(kFALSE);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);
  c1->SetBottomMargin(0.103939);
  //gStyle->SetPalette(1,0);
  //Re-implementing red-yellow palette
  TColor::InitializeColors();
  const Int_t nRGBs = 3;
  Double_t stops[nRGBs] = { 0.00, 0.50, 1.00};
  Double_t red[nRGBs]   = { 0.559, 1.00, 1.00};
  Double_t green[nRGBs] = { 0.00, 0.00, 1.00};
  Double_t blue[nRGBs]  = { 0.10, 0.00, 0.00};
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255);

  TChain wave("waveburst");
  TChain live("liveTime");
  wave.Add(net_file_name);
  live.Add(liv_file_name);
  netevent  W(&wave,nIFO);
  wave.SetEstimate(wave.GetEntries());

  // check vetoes
  bool bcat2  = false;
  bool bhveto = false;
  bool bcat3  = false;

  for(int i=0;i<nVDQF;i++) if(VDQF[i].cat==CWB_CAT2)  bcat2=true;
  for(int i=0;i<nVDQF;i++) if(VDQF[i].cat==CWB_HVETO) bhveto=true;
  for(int i=0;i<nVDQF;i++) if(VDQF[i].cat==CWB_CAT3)  bcat3=true;

  for(int i=0;i<nVDQF;i++) {
    if((VDQF[i].cat!=CWB_CAT2)&&(VDQF[i].cat!=CWB_HVETO)&&(VDQF[i].cat!=CWB_CAT3)) continue;
    char veto_name[256];sprintf(veto_name,"%s_%s",CWB_CAT_NAME[VDQF[i].cat].Data(),VDQF[i].ifo);
    if(!CWB::Toolbox::isLeafInTree(W.fChain,veto_name)) {
       cout << "cwb_report_prod_1.C Error : Leaf " << veto_name
            << " is not present in tree" << endl;
       gSystem->Exit(1);
    }
  }

  gSystem->Exec("date");

  char fname[1024];
  char ch1[256];
  char ch3[256];

  int n;
  double sTARt, sTOp;

  W.GetEntry(0); sTARt = W.gps[0];
  sTOp = W.gps[0];
  for(int k=1;k<wave.GetEntries();k++){
    W.GetEntry(k);
    if(int(k/100000)*100000 == k) cout << "processed: " << k <<endl;
    if(W.gps[0]>sTOp)  {sTOp  = W.gps[0];}
    if(W.gps[0]<sTARt) {sTARt = W.gps[0];}
  }
  double T_bgn = sTARt;
  double T_end = sTOp;
  T_bgn-=(T_end-T_bgn)/100.;
  T_end+=(T_end-T_bgn)/100.;
  int nBins = 100; 
  sprintf(fname,"%s/time.txt",netdir);
  FILE* ftrig = fopen(fname,"w");
  if(ftrig==NULL) {
    cout << "cwb_report_prod_1.C - File open error : " << fname << endl;
    gSystem->Exit(1);
  }
  fprintf(ftrig,"#Start\tStop\n");
  fprintf(ftrig,"%d\t%d\n",(int)T_bgn,(int)T_end);
  fclose(ftrig);
  sprintf(fname,"Red dots vetoed CAT3 or hveto");

  TH2F* rho_T = new TH2F("rho_T","",nBins,T_bgn,T_end,pp_rho_bin,pp_rho_min,pp_rho_max);
  rho_T->SetTitleOffset(1.3,"Y");
  rho_T->GetXaxis()->SetTitle(fname);
  rho_T->GetXaxis()->SetTimeDisplay(1);
  rho_T->GetYaxis()->SetTitle("#rho");
  rho_T->SetMarkerStyle(20);
  rho_T->SetMarkerColor(1);
  rho_T->SetMarkerSize(0.8);
  rho_T->SetStats(kFALSE);

  TH2F* rho_T2 = new TH2F("rho_T2","",nBins,T_bgn,T_end,pp_rho_bin,pp_rho_min,pp_rho_max);
  rho_T2->SetTitle("Red dots vetoed CAT3 or hveto");
  rho_T2->GetXaxis()->SetTitle(fname);
  rho_T2->GetXaxis()->SetTimeDisplay(1);
  rho_T2->GetYaxis()->SetTitle("#rho");
  rho_T2->SetMarkerStyle(20);
  rho_T2->SetMarkerColor(2);
  rho_T2->SetMarkerSize(0.8);
  rho_T2->SetStats(kFALSE);

  if(c1) delete c1; 
  c1 = new TCanvas("c","C",0,0,800,600);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetBottomMargin(0.143939);
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);

  int nfreq = (freqHigh-freqLow)/1.;
  TH2F* ff=new TH2F("ff","ff",nfreq,freqLow,freqHigh,pp_rho_bin,pp_rho_min,pp_rho_max);
  ff->SetTitleOffset(1.3,"Y");
  ff->SetTitle("");
  ff->GetXaxis()->SetTitle("frequency, Hz");
  ff->SetStats(kFALSE);

  c1->Clear();
  sprintf(ch1,"rho[%d]:frequency[0]>>ff",pp_irho);
  c1->SetLogz(kTRUE);
  ff->GetYaxis()->SetTitle("#rho(ECOR)");
  Draw(wave,ch1,ch2,"colz");
  sprintf(fname,"%s/ECOR_frequency.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  c1->Clear();
  sprintf(ch1,"rho[0]:frequency[0]>>ff");
  c1->SetLogz(kTRUE);
  ff->GetYaxis()->SetTitle("#rho(ecor)");
  Draw(wave,ch1,ch2,"colz");
  sprintf(fname,"%s/ecor_frequency.eps",netdir);
  cout << fname << endl;
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  if(pp_rho_log) c1->SetLogy(kTRUE);
  sprintf(ch1,"rho[%d]:frequency[0]>>ff",pp_irho);
  c1->SetLogz(kTRUE);
  if(fHigh-fLow>256) c1->SetLogx(kTRUE);
  ff->SetTitle("after pp-cuts (loudest not vetoed : black)");
  ff->GetYaxis()->SetTitle("#rho");
  ff->GetXaxis()->SetMoreLogLabels(kTRUE);
  ff->GetXaxis()->SetNoExponent(kTRUE);
  Draw(wave,ch1,ch2,"colz");

  // set marker on the first not vetoed pp_max_nloudest_list loudest events
  sprintf(ch1,"rho[%d]:frequency[0]:Entry$",pp_irho);
  Draw(wave,ch1,ch2,"goff");
  int sel_size = wave.GetSelectedRows();
  double* sel_rho   = wave.GetV1();
  double* sel_freq  = wave.GetV2();
  double* sel_entry = wave.GetV3();
  Int_t *_index = new Int_t[sel_size];
  TMath::Sort(sel_size,sel_rho,_index,true);
  if(sel_size>pp_max_nloudest_list) sel_size=pp_max_nloudest_list;
  // get veto branches
  UChar_t VETO[100];
  TString VETO_NAME[100];
  int nveto=0;
  size_t nbranch = wave.GetListOfBranches()->GetEntries();
  for(int n=0;n<nbranch;++n) {
    TBranch *br =dynamic_cast<TBranch*>(wave.GetListOfBranches()->At(n));
    if(TString(br->GetName()).Contains("veto_")) {
      // cat2 is a special case, the events tagged with cat2 are not vetoed
      if(!TString(br->GetName()).Contains("veto_cat2")) {
        VETO_NAME[nveto] = br->GetName();
        wave.SetBranchAddress(br->GetName(),&VETO[nveto++]);
      }
    }
  }
  for(int i=0;i<sel_size;i++) {
    int l=_index[i];
    int j=int(sel_entry[l]);
    wave.GetEntry(j);
    bool veto=false;
    for(int k=0;k<nveto;k++) {
#ifdef CAT3_VETO
      if(VETO[k] && VETO_NAME[k].Contains("veto_cat3"))  veto=true;
#endif
#ifdef HVETO_VETO
      if(VETO[k] && VETO_NAME[k].Contains("veto_hveto")) veto=true;
#endif
#ifdef USER_VETO
      if(VETO[k] && VETO_NAME[k].Contains("veto_user"))  veto=true;
#endif
    }
    if(!veto) {
      // set black marker for not vetoed loudest events
      TMarker *mP = new TMarker(sel_freq[l], sel_rho[l], 20);
      mP->SetMarkerSize(0.9);
      mP->SetMarkerColor(kBlack);
      mP->Draw();
    }
  }
  sprintf(fname,"%s/rho_frequency.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
 
  c1->Clear();
  sprintf(ch1,"rho[%d]:time[0]>>rho_T",pp_irho);
  c1->SetLogz(kFALSE);
  c1->SetLogx(kFALSE);
  if(pp_rho_log) c1->SetLogy(kTRUE); else c1->SetLogy(kFALSE);
  if(bhveto||bcat3) { 
    sprintf(ch1,"rho[%d]:time[0]>>rho_T",pp_irho);
    Draw(wave,ch1,veto_not_vetoed,"");
    sprintf(ch1,"rho[%d]:time[0]>>rho_T2",pp_irho);
    Draw(wave,ch1,veto_vetoed,"P SAME");
    rho_T->SetTitle("after pp-cuts (black), vetoed by cat3 or hveto (red)");
  } else {
    Draw(wave,ch1,ch2,"");
    rho_T->SetTitle("after pp-cuts");
  }
  sprintf(fname,"%s/rho_time.eps",netdir);
  c1->Update(); c1->Print(fname);

 /* if (online)
  {
  //ONLINE
  if(c1) delete c1;
  c1 = new TCanvas("c","C",0,0,1200,800);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetBottomMargin(0.143939);
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);

  char s_gracedb[2048];
  sprintf(s_gracedb,"%s && rho[%i]>%f && rho[%i]<%f",ch2,pp_irho,t_rho_lum,pp_irho,t_rho_off);
  char s_offline[2048];
  sprintf(s_offline,"%s && rho[%i]>%f",ch2,pp_irho,t_rho_off);
  char notch2[2048];
  sprintf(notch2,"!(%s && rho[%i]>%f)",ch2,pp_irho,t_rho_lum);
  c1->Clear();
  c1->SetLogz(kFALSE);
  c1->SetLogx(kFALSE);

  TH2F* rho_on = new TH2F("rho_on","",nBins,T_bgn,T_end,pp_rho_bin,pp_rho_min,pp_rho_max);
  rho_on->SetTitleOffset(1.3,"Y");
  rho_on->GetXaxis()->SetTitle(fname);
  rho_on->GetXaxis()->SetTimeDisplay(1);
  rho_on->GetYaxis()->SetTitle("#rho");
  rho_on->SetMarkerStyle(20);
  rho_on->SetMarkerColor(1);
  rho_on->SetMarkerSize(0.8);
  rho_on->SetStats(kFALSE);

  TH2F* rho_on2 = new TH2F("rho_on2","",nBins,T_bgn,T_end,pp_rho_bin,pp_rho_min,pp_rho_max);
  rho_on2->GetXaxis()->SetTitle(fname);
  rho_on2->GetXaxis()->SetTimeDisplay(1);
  rho_on2->GetYaxis()->SetTitle("#rho");
  rho_on2->SetMarkerStyle(20);
  rho_on2->SetMarkerColor(2);
  rho_on2->SetMarkerSize(0.8);
  rho_on2->SetStats(kFALSE);

  TH2F* rho_on3 = new TH2F("rho_on3","",nBins,T_bgn,T_end,pp_rho_bin,pp_rho_min,pp_rho_max);
  rho_on3->GetXaxis()->SetTitle(fname);
  rho_on3->GetXaxis()->SetTimeDisplay(1);
  rho_on3->GetYaxis()->SetTitle("#rho");
  rho_on3->SetMarkerStyle(20);
  rho_on3->SetMarkerColor(3);
  rho_on3->SetMarkerSize(0.8);
  rho_on3->SetStats(kFALSE);

  sprintf(ch1,"rho[%d]:time[0]>>rho_on",pp_irho);
  Draw(wave,ch1,notch2,"P");
  sprintf(ch1,"rho[%d]:time[0]>>rho_on2",pp_irho);
  Draw(wave,ch1,s_gracedb,"P SAME");
  sprintf(ch1,"rho[%d]:time[0]>>rho_on3",pp_irho);
  Draw(wave,ch1,s_offline,"P SAME");
  rho_on->GetXaxis()->SetTitle("time, gps");
  rho_on->SetTitle("all (black), gracedb (red), offline (green)");
  sprintf(fname,"%s/rho_time_online.eps",netdir);
  c1->Update(); c1->Print(fname);

  TH2F* ff_on=new TH2F("ff_on","ff_on",nfreq,freqLow,freqHigh,pp_rho_bin,pp_rho_min,pp_rho_max);
  ff_on->SetTitleOffset(1.3,"Y");
  ff_on->SetTitle("");
  ff_on->GetXaxis()->SetTitle("frequency, Hz");
  ff_on->SetStats(kFALSE);
  ff_on->SetMarkerStyle(20);
  ff_on->SetMarkerColor(1);
  ff_on->SetMarkerSize(0.8);

  TH2F* ff_on2=new TH2F("ff_on2","ff_on2",nfreq,freqLow,freqHigh,pp_rho_bin,pp_rho_min,pp_rho_max);
  ff_on2->SetTitleOffset(1.3,"Y");
  ff_on2->SetTitle("");
  ff_on2->GetXaxis()->SetTitle("frequency, Hz");
  ff_on2->SetStats(kFALSE);
  ff_on2->SetMarkerStyle(20);
  ff_on2->SetMarkerColor(2);
  ff_on2->SetMarkerSize(0.8);

  TH2F* ff_on3=new TH2F("ff_on3","ff_on3",nfreq,freqLow,freqHigh,pp_rho_bin,pp_rho_min,pp_rho_max);
  ff_on3->SetTitleOffset(1.3,"Y");
  ff_on3->SetTitle("");
  ff_on3->GetXaxis()->SetTitle("frequency, Hz");
  ff_on3->SetStats(kFALSE);
  ff_on3->SetMarkerStyle(20);
  ff_on3->SetMarkerColor(3);
  ff_on3->SetMarkerSize(0.8);

  c1->Clear();
  ff_on->SetTitle("all (black), lumin (red), offline (green)");
  ff_on->GetYaxis()->SetTitle("#rho");
  ff_on->GetXaxis()->SetMoreLogLabels(kTRUE);
  ff_on->GetXaxis()->SetNoExponent(kTRUE);
  sprintf(ch1,"rho[%d]:frequency[0]>>ff_on",pp_irho);
  Draw(wave,ch1,notch2,"P");
  sprintf(ch1,"rho[%d]:frequency[0]>>ff_on2",pp_irho);
  Draw(wave,ch1,s_gracedb,"P SAME");
  sprintf(ch1,"rho[%d]:frequency[0]>>ff_on3",pp_irho);
  Draw(wave,ch1,s_offline,"P SAME");
  sprintf(fname,"%s/rho_frequency_online.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  TH2F* tf_on=new TH2F("tf_on","tf_on",nBins,T_bgn,T_end,nfreq,freqLow,freqHigh);
  tf_on->SetTitleOffset(1.3,"Y");
  tf_on->SetTitle("");
  tf_on->GetXaxis()->SetTitle("frequency, Hz");
  tf_on->SetStats(kFALSE);
  tf_on->SetMarkerStyle(20);
  tf_on->SetMarkerColor(1);
  tf_on->SetMarkerSize(0.8);

  TH2F* tf_on2=new TH2F("tf_on2","tf_on2",nBins,T_bgn,T_end,nfreq,freqLow,freqHigh);
  tf_on2->SetTitleOffset(1.3,"Y");
  tf_on2->SetTitle("");
  tf_on2->GetXaxis()->SetTitle("frequency, Hz");
  tf_on2->SetStats(kFALSE);
  tf_on2->SetMarkerStyle(20);
  tf_on2->SetMarkerColor(2);
  tf_on2->SetMarkerSize(0.8);

  TH2F* tf_on3=new TH2F("tf_on3","tf_on3",nBins,T_bgn,T_end,nfreq,freqLow,freqHigh);
  tf_on3->SetTitleOffset(1.3,"Y");
  tf_on3->SetTitle("");
  tf_on3->GetXaxis()->SetTitle("frequency, Hz");
  tf_on3->SetStats(kFALSE);
  tf_on3->SetMarkerStyle(20);
  tf_on3->SetMarkerColor(3);
  tf_on3->SetMarkerSize(0.8);

  c1->Clear();
  tf_on->SetTitle("all (black), lumin (red), offline (green)");
  tf_on->GetYaxis()->SetTitle("frequency");
  tf_on->GetXaxis()->SetTitle("time");
  tf_on->GetXaxis()->SetTimeDisplay(1);
  tf_on->GetXaxis()->SetMoreLogLabels(kTRUE);
  tf_on->GetXaxis()->SetNoExponent(kTRUE);
  sprintf(ch1,"frequency[0]:time[0]>>tf_on");
  Draw(wave,ch1,notch2,"P");
  sprintf(ch1,"frequency[0]:time[0]>>tf_on2");
  Draw(wave,ch1,s_gracedb,"P SAME");
  sprintf(ch1,"frequency[0]:time[0]>>tf_on3");
  Draw(wave,ch1,s_offline,"P SAME");
  sprintf(fname,"%s/time_frequency_online.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  }
*/

  gSystem->Exit(0);
}
