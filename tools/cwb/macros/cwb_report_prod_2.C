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
  cout<<"cwb_report_prod_2.C starts..."<<endl;

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
  #include <map>
  
  gStyle->SetTitleOffset(1.4,"X");
  //gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleOffset(2.0,"Y");
  gStyle->SetLabelOffset(0.010,"X");
  gStyle->SetLabelOffset(0.010,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetLabelSize(0.06,"X");
  gStyle->SetLabelSize(0.06,"Y");

  gStyle->SetTitleH(0.060);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.99);
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

  if(nIFO==1) { // Single Detector Mode
    CWB::config config;
    config.Import();
    config.SetSingleDetectorMode();
    config.Export();
  }

  // initialize the ifo name (takes into account the non built-in detectors)
  char IFO[NIFO_MAX][32];
  for(int n=0;n<nIFO;n++) {
    if(strlen(ifo[n])!=0) strcpy(IFO[n],ifo[n]); 
    else                  strcpy(IFO[n],detParms[n].name);
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
  gStyle->SetPalette(1,0);

  TChain wave("waveburst");
  TChain live("liveTime");
  wave.Add(net_file_name);
  live.Add(liv_file_name);
  netevent  W(&wave,nIFO);

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
       cout << "cwb_report_prod_2.C Error : Leaf " << veto_name
            << " is not present in tree" << endl;
       gSystem->Exit(1);
    }
  }

  // check if slag is presents in liv tree
  TBranch* branch;
  bool slagFound=false;
  TIter next(wave.GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if (TString("slag").CompareTo(branch->GetName())==0) slagFound=true;
  }
  next.Reset();
  if(!slagFound) {
    cout << "cwb_report_prod_2.C : Error - wave tree do not contains slag leaf : " 
         << net_file_name << endl;
    gSystem->Exit(1);
  }

  bool    saveVETO[100];      
  UChar_t     VETO[100];      
  int    countVETO[100];
  TString nameVETO[100];      // name of vetoes
 
  // build XDQF structure used for internal's computations
  int nXDQF=0;
  dqfile* XDQF = new dqfile[nVDQF];
  for(int j=0;j<nVDQF;j++) {
    if(VDQF[j].cat==CWB_CAT0) continue;
    if(VDQF[j].cat==CWB_CAT1) continue;
    bool bifo=false;
    for(int i=0;i<nIFO;i++) if(TString(IFO[i]).CompareTo(VDQF[j].ifo)==0) bifo=true;
    if(bifo) XDQF[nXDQF++]=VDQF[j];
  }
  for(int j=0;j<nXDQF;j++) {
    cout << "XDQF[" << j << "] " << CWB_CAT_NAME[XDQF[j].cat] << "_" << XDQF[j].ifo << endl;
  }

  for(int j=0;j<nXDQF;j++) {
    char veto_name[256];
    sprintf(veto_name,"%s_%s",CWB_CAT_NAME[XDQF[j].cat].Data(),XDQF[j].ifo);
    nameVETO[j]=veto_name;
    countVETO[j]=0;
    W.fChain->SetBranchAddress(nameVETO[j].Data(),&VETO[j]);
  }

  gSystem->Exec("date");

  double xcor, acor, rho, asnr, xCOR, xcor2, xcor3, vED;
  int i,j,n,m,k;
  bool save;
  bool save_veto;

  size_t ntype = 0;
  char fname[1024];
  char _title[1024];
  char ch1[256];
  char ch3[256];

  wavearray<double> Xcut(pp_rho_bin);
  wavearray<double> Ycut(pp_rho_bin); Ycut = 0.;
  wavearray<double> yCUT(pp_rho_bin); yCUT = 0.;
  wavearray<double> Xerr(pp_rho_bin); Xerr = 0.;
  wavearray<double> Yerr(pp_rho_bin);
  wavearray<double> yERR(pp_rho_bin);
  wavearray<double> Ycut_veto(pp_rho_bin); Ycut_veto = 0.;
  wavearray<double> Yerr_veto(pp_rho_bin);
  wavearray<double> Ycut_multi(pp_rho_bin); Ycut_multi = 0.;  //MULTI
  wavearray<double> Yerr_multi(pp_rho_bin);  //MULTI

  // create/init structures for FAR output statistic

  // find the min and max in the tree and set rho limits for far_rho.txt file
  int nsize = wave.GetEntries();
  wave.SetEstimate(nsize);
  wave.Draw(TString::Format("rho[%d]",pp_irho).Data(),"","goff");
  double rho_min = TMath::MinElement(nsize, wave.GetV1());
  double rho_max = TMath::MaxElement(nsize, wave.GetV1());

  float far_rho_min = TMath::Min(pp_rho_min,(double)TMath::Nint(rho_min-0.5));
  float far_rho_max = TMath::Max(pp_rho_max,(double)TMath::Nint(rho_max+0.5));
  float far_rho_wbin = 0.01;
  int far_rho_bin = float(far_rho_max-far_rho_min)/far_rho_wbin;
  double far_drho = double(far_rho_max-far_rho_min)/double(far_rho_bin);

  wavearray<double> Xfar(far_rho_bin);
  wavearray<double> Yfar(far_rho_bin); Yfar = 0.;


// calculate live time

  float xlag;
  Int_t xrun, rUn;
  double xlive, sTARt, sTOp, gps;
  double liveTot = 0.;
  double liveZero = 0.;
  double Live;

  wavearray<double> Trun(500000); Trun = 0.;
  wavearray<double> Wlag[NIFO_MAX+1];
  wavearray<double> Wslag[NIFO_MAX+1];
  wavearray<double> Tdlag;
  wavearray<double> Tlag;
  wavearray<double> Rlag;
  wavearray<double> Elag;
  wavearray<double> Olag;

  cout << "Start CWB::Toolbox::getLiveTime : be patient, it takes a while ..." << endl; 
  liveTot=CWB::Toolbox::getLiveTime(nIFO,live,Trun,Wlag,Wslag,Tlag,Tdlag,cwb_lag_number,cwb_slag_number);

  // use std::map to make faster the search 
  std::map<std::pair<int,int>, int> lagslag;
  for(j=0; j<Wlag[nIFO].size(); j++) {
    lagslag[std::make_pair(int(Wslag[nIFO].data[j]),int(Wlag[nIFO].data[j]))] = j;
  }

  Rlag = Tlag; Rlag = 0.;
  Elag = Tlag; Elag = 0.;
  Olag = Tlag; Olag = 0.;

  sprintf(fname,"%s/live.txt",netdir);
  FILE* flive = fopen(fname,"w");
  if(flive==NULL) {
    cout << "cwb_report_prod_2.C Error : Error opening " << fname << endl;
    gSystem->Exit(1);
  }

  // get live time zero (when is not included into Tlag array)
  if((cwb_lag_number!=0)&&(cwb_slag_number!=0)) { 
#ifdef _USE_ROOT6
    // with ROOT6 live crash, it is necessary to open live again, why?   
    TChain live1("liveTime");
    live1.Add(liv_file_name);
    liveZero=CWB::Toolbox::getZeroLiveTime(nIFO,live1);
#else
    liveZero=CWB::Toolbox::getZeroLiveTime(nIFO,live);
#endif
    fprintf(flive,"slag=%i\t lag=%i\t ",0,0);
    for (int jj=0; jj<nIFO; jj++) fprintf(flive,"%s:slag=%5.2f\tlag=%5.2f\t",IFO[jj],0.,0.);
    fprintf(flive,"live=%9.2f\n",liveZero);
    printf("live time: zero lags = %10.1f \n",liveZero); 
  }
  for(j=0; j<Wlag[nIFO].size(); j++){
    fprintf(flive,"slag=%i\t lag=%i\t ",(int)Wslag[nIFO].data[j],(int)Wlag[nIFO].data[j]);
    for (int jj=0; jj<nIFO; jj++) {
      fprintf(flive,"%s:slag=%5.2f\tlag=%5.2f\t",IFO[jj],Wslag[jj].data[j],Wlag[jj].data[j]);
    }
    fprintf(flive,"live=%9.2f\n",Tlag.data[j]);
  }
  fprintf(flive,"nonzero lags live=%10.1f \n",liveTot); 
  fclose(flive);
  printf("total live time: non-zero lags = %10.1f \n",liveTot); 
  cout << "End   CWB::Toolbox::getLiveTime : ";cout.flush();gSystem->Exec("/bin/date");

  //int const nlag=Wlag[nIFO].size()+1;
  int const nlag=Wlag[nIFO].size();  // FIX
  int min_lag=kMaxInt;
  int max_lag=0;
  for(int i=0;i<nlag;i++) if(Wlag[nIFO][i]<min_lag) min_lag=Wlag[nIFO][i];
  for(int i=0;i<nlag;i++) if(Wlag[nIFO][i]>max_lag) max_lag=Wlag[nIFO][i];
  int* iy_lag = new int[max_lag+1]; // this array is used because the lag numbers in general are not contiguous
  for(int i=0;i<=max_lag;i++) iy_lag[i]=-1;
  for(int i=0;i<nlag;i++) iy_lag[int(Wlag[nIFO][i])-min_lag]=i;
  double** y_lag = new double*[nlag]; //Ycut = 0.;
  double** y_lag_veto = new double*[nlag]; //Ycut = 0.;
  for(int i=0;i<nlag;i++) {
    y_lag[i]      = new double[pp_rho_bin]; 
    y_lag_veto[i] = new double[pp_rho_bin]; 
  }
//  double y_lag[nlag][pp_rho_bin]; //Ycut = 0.;
//  double y_lag_veto[nlag][pp_rho_bin]; //Ycut = 0.;
  for (int i=0;i<nlag;i++) for (int j=0; j<pp_rho_bin; j++) {y_lag[i][j] = 0.;y_lag_veto[i][j] = 0.;}

// histograms

  TH1F* xCor = new TH1F("xCor","",100,0.4,1.);
  xCor->SetTitleOffset(1.3,"Y");
  xCor->SetTitle("");
  xCor->GetXaxis()->SetTitle("network correlation");
  xCor->GetYaxis()->SetTitle("events");

  TH1F* dens = new TH1F("dens","",100,0.,3);
  dens->SetTitleOffset(1.3,"Y");
  dens->SetTitle("");
  dens->GetXaxis()->SetTitle("-log10(cluster density)");
  dens->GetYaxis()->SetTitle("events");

  TH1F* hpf = new TH1F("hpf","",100,0.,1.);
  hpf->SetTitleOffset(1.3,"Y");
  hpf->SetTitle("");
  hpf->GetXaxis()->SetTitle("penalty factor");
  hpf->GetYaxis()->SetTitle("events");

  TH1F* hvED = new TH1F("hvED","",200,0,2.);
  hvED->SetTitleOffset(1.3,"Y");
  hvED->SetTitle("");
  hvED->GetXaxis()->SetTitle("hvED consistency");
  hvED->GetYaxis()->SetTitle("events");

  TH1F* ecor = new TH1F("ecor","",200,0,10.);  
  ecor->SetTitleOffset(1.3,"Y");
  ecor->SetTitle("");
  ecor->GetXaxis()->SetTitle("correlated amplitude");
  ecor->GetYaxis()->SetTitle("events");

  TH1F* ECOR = new TH1F("ECOR","",200,0,10.);  
  ECOR->SetTitleOffset(1.3,"Y");
  ECOR->SetTitle("");
  ECOR->GetXaxis()->SetTitle("#rho(ECOR))");
  ECOR->GetYaxis()->SetTitle("events");

  TH1F* Like = new TH1F("Like","",100,1,4.);
  Like->SetTitleOffset(1.3,"Y");
  Like->SetTitle("");
  Like->GetXaxis()->SetTitle("log10(likelihood)");
  Like->GetYaxis()->SetTitle("events");
  
  TH1F* Mchirp = new TH1F("Mchirp","",100,0.,40.);
  Mchirp->SetTitleOffset(1.3,"Y");
  Mchirp->SetTitle("After pp cuts");
  Mchirp->GetXaxis()->SetTitle("reconstructed chirp mass(M_{#odot})");
  Mchirp->GetYaxis()->SetTitle("events");

  TH1F* McErr = new TH1F("McErr","",100,0.,1.);
  McErr->SetTitleOffset(1.3,"Y");
  McErr->SetTitle("After pp cuts");
  McErr->GetXaxis()->SetTitle("pixels used in reconstruction/total");
  McErr->GetYaxis()->SetTitle("events");

  TH1F* cutF = new TH1F("cutF","",24,60,2048.);//,300,2,6);
  cutF->SetTitleOffset(1.2,"Y");
  cutF->GetXaxis()->SetTitle("frequency, Hz");
  //cutF->GetYaxis()->SetTitle("#rho");
  cutF->SetStats(kFALSE);
  cutF->SetFillColor(2);
  //cutF->SetMarkerStyle(20);
  //cutF->SetMarkerSize(0.6);
  //cutF->SetMarkerColor(2);
 
  TH2F* rhocc = new TH2F("rhocc","",100,0.,1.,100,pp_rho_min,pp_rho_max);
  rhocc->SetTitle("0 < cc < 1");
  rhocc->SetTitleOffset(1.3,"Y");
  rhocc->GetXaxis()->SetTitle("network correlation");
  rhocc->GetYaxis()->SetTitle("#rho");
  rhocc->SetStats(kFALSE);
  rhocc->SetMarkerStyle(20);
  rhocc->SetMarkerSize(0.5);
  rhocc->SetMarkerColor(1);

  TH2F* ECOR_cc = new TH2F("ECOR_cc","",100,0.,1.,500,0,15.);
  ECOR_cc->SetTitleOffset(1.3,"Y");
  ECOR_cc->GetXaxis()->SetTitle("network correlation");
  ECOR_cc->GetYaxis()->SetTitle("#rho(ECOR)");
  ECOR_cc->SetStats(kFALSE);
  ECOR_cc->SetMarkerStyle(20);
  ECOR_cc->SetMarkerSize(0.5);
  ECOR_cc->SetMarkerColor(1);

  TH2F* ecor_cc = new TH2F("ecor_cc","",100,0.,1.,500,0,15.);
  ecor_cc->SetTitleOffset(1.3,"Y");
  ecor_cc->GetXaxis()->SetTitle("network correlation");
  ecor_cc->GetYaxis()->SetTitle("#rho(ecor)");
  ecor_cc->SetStats(kFALSE);
  ecor_cc->SetMarkerStyle(20);
  ecor_cc->SetMarkerSize(0.5);
  ecor_cc->SetMarkerColor(1);

  TH2F* rho_subnet = new TH2F("rho_subnet","",100,0.,1.,100,pp_rho_min,pp_rho_max);
  rho_subnet->SetTitle("0 < subnet < 1");
  rho_subnet->SetTitleOffset(1.3,"Y");
  rho_subnet->GetXaxis()->SetTitle("subnet");
  rho_subnet->GetYaxis()->SetTitle("#rho");
  rho_subnet->SetMarkerStyle(20);
  rho_subnet->SetMarkerColor(1);
  rho_subnet->SetMarkerSize(0.5);
  rho_subnet->SetStats(kFALSE);
 
  TH2F* rho_vED = new TH2F("rho_vED","",100,0.,1.,100,pp_rho_min,pp_rho_max);
  rho_vED->SetTitle("0 < vED < 1");
  rho_vED->SetTitleOffset(1.3,"Y");
  rho_vED->GetXaxis()->SetTitle("vED");
  rho_vED->GetYaxis()->SetTitle("#rho");
  rho_vED->SetMarkerStyle(20);
  rho_vED->SetMarkerColor(1);
  rho_vED->SetMarkerSize(0.5);
  rho_vED->SetStats(kFALSE);

  TH2F* rho_pf;  
  if(TString(analysis)=="1G") {
    rho_pf = new TH2F("rho_pf","",100,0.,1.,100,pp_rho_min,pp_rho_max);
    rho_pf->SetTitle("0 < penalty < 1");
    rho_pf->GetXaxis()->SetTitle("penalty");
  } else {                       
    rho_pf = new TH2F("rho_pf","",100,-1.,2.,100,pp_rho_min,pp_rho_max);
    rho_pf->SetTitle("chi2");
    rho_pf->GetXaxis()->SetTitle("log10(chi2)");
  }
  rho_pf->SetTitleOffset(1.3,"Y");
  rho_pf->GetYaxis()->SetTitle("#rho");
  rho_pf->SetMarkerStyle(20);
  rho_pf->SetMarkerColor(1);
  rho_pf->SetMarkerSize(0.5);
  rho_pf->SetStats(kFALSE);
 
  TH2F* pf_cc = new TH2F("pf_cc","",100,0.,1.,100,0,1.);
  pf_cc->SetTitleOffset(1.3,"Y");
  pf_cc->GetXaxis()->SetTitle("network correlation");
  pf_cc->GetYaxis()->SetTitle("penalty");
  pf_cc->SetMarkerStyle(20);
  pf_cc->SetMarkerColor(1);
  pf_cc->SetMarkerSize(0.8);
  pf_cc->SetStats(kFALSE);
 
  TH2F* pf_vED = new TH2F("pf_vED","",100,0.,1.,100,0,1.);
  pf_vED->SetTitleOffset(1.3,"Y");
  pf_vED->GetXaxis()->SetTitle("vED consistency");
  pf_vED->GetYaxis()->SetTitle("penalty");
  pf_vED->SetMarkerStyle(20);
  pf_vED->SetMarkerColor(1);
  pf_vED->SetMarkerSize(0.8);
  pf_vED->SetStats(kFALSE);
 
  TH2F* pf_ed = new TH2F("pf_ed","",100,0.,1.,100,0.,1.);
  pf_ed->SetTitleOffset(1.3,"Y");
  pf_ed->GetXaxis()->SetTitle("energy disbalance");
  pf_ed->GetYaxis()->SetTitle("penalty");
  pf_ed->SetMarkerStyle(20);
  pf_ed->SetMarkerColor(1);
  pf_ed->SetMarkerSize(0.8);
  pf_ed->SetStats(kFALSE);
  
  TH2F* rho_L = new TH2F("rho_L","",500,0.,5.,500,0,5.);
  rho_L->SetTitleOffset(1.3,"Y");
  rho_L->GetYaxis()->SetTitle("log10(average SNR)");
  rho_L->GetXaxis()->SetTitle("log10(2*#rho^2)");
  rho_L->SetStats(kFALSE);
  rho_L->SetMarkerStyle(20);
  rho_L->SetMarkerSize(0.5);
  rho_L->SetMarkerColor(1);
  
  TH2F* rho_mchirp = new TH2F("rho_mchirp","",200,0.,30.,100,pp_rho_min,pp_rho_max);
  rho_mchirp->SetTitle("rho vs Mchirp(after pp cuts)");
  rho_mchirp->SetTitleOffset(1.3,"Y");
  rho_mchirp->GetXaxis()->SetTitle("Mchirp");
  rho_mchirp->GetYaxis()->SetTitle("#rho");
  rho_mchirp->SetMarkerStyle(20);
  rho_mchirp->SetMarkerColor(1);
  rho_mchirp->SetMarkerSize(0.5);
  rho_mchirp->SetStats(kFALSE);  

  //ONLINE
  skymap sm_theta_phi(3);
  for (int qq=0;qq<sm_theta_phi.size(); qq++) sm_theta_phi.set(qq,0);
  skymap sm_ra_dec(3);
  for (int qq=0;qq<sm_ra_dec.size(); qq++) sm_ra_dec.set(qq,0);

  sprintf(fname,"%s/time.txt",netdir);
  ifstream in;
  in.open(fname);
  char lab[256];
  in.getline(lab,256);
  in >> sTARt >> sTOp; 
  in.close();
  
  double T_bgn = sTARt;
  double T_end = sTOp;
  T_bgn-=hours*3600;
  T_end+=hours*3600;
  int nBins = int((T_end-T_bgn)/hours/3600.);
  sprintf(fname,"time (%d hours per bin)",int(hours));
  
  TH1F* hrho[2];
  int hrho_col[2] = {1,3};
  for(int j=0; j<2; j++) {
    sprintf(ch1,"hrho%d",j);
    hrho[j] = new TH1F(ch1,ch1,Xcut.size(),pp_rho_min,pp_rho_max);
    hrho[j]->SetTitleOffset(1.3,"Y");
    sprintf(_title,"#rho (bin width = %2.1f)",(pp_rho_max-pp_rho_min)/pp_rho_bin);
    hrho[j]->GetXaxis()->SetTitle(_title);
    hrho[j]->GetYaxis()->SetTitle("events");
    hrho[j]->GetYaxis()->SetTitleOffset(2.0);
    hrho[j]->GetYaxis()->CenterTitle(true);
    hrho[j]->SetLineColor(hrho_col[j]);
    hrho[j]->SetFillColor(hrho_col[j]);
    hrho[j]->SetStats(kFALSE);
    hrho[j]->SetMinimum(0.5);

  }
  hrho[0]->SetTitle("full events (black), after pp cuts (green)");

  TH1F* hF[NIFO_MAX+1];
  TH1F* hF2[NIFO_MAX+1];
  for(int j=0; j<nIFO+1; j++){
    sprintf(ch1,"hF%d",j);
    hF[j] = new TH1F(ch1,ch1,nBins,T_bgn,T_end);
    hF[j]->GetXaxis()->SetTitle(fname);
    hF[j]->GetXaxis()->SetTitleSize(0.1);
    hF[j]->GetXaxis()->SetLabelSize(0.1);
    hF[j]->GetXaxis()->SetTitleOffset(0.9);
    sprintf(ch1,"hF2%d",j);
    hF2[j] = new TH1F(ch1,ch1,24,32,2048.);
    hF2[j]->GetXaxis()->SetTitle("frequency, (80 Hz per bin)");
    hF2[j]->GetXaxis()->SetTitleSize(0.1);
    hF2[j]->GetXaxis()->SetLabelSize(0.1);
    hF2[j]->GetXaxis()->SetTitleOffset(0.9);
    if(j) {
      sprintf(_title,"%s : fraction (%%)",IFO[j-1]);
      hF[j]->GetYaxis()->SetTitle(_title);
      hF2[j]->GetYaxis()->SetTitle(_title);
    }
    hF[j]->GetYaxis()->SetTitleSize(0.12);
    hF[j]->GetYaxis()->SetLabelSize(0.1);
    hF[j]->GetYaxis()->SetTitleOffset(0.4);
    hF[j]->SetStats(kFALSE);
    hF[j]->SetLineColor(1);
    hF[j]->SetFillColor(1);
    hF[j]->SetTitle("full events");
    hF2[j]->GetYaxis()->SetTitleSize(0.12);
    hF2[j]->GetYaxis()->SetLabelSize(0.1);
    hF2[j]->GetYaxis()->SetTitleOffset(0.4);
    hF2[j]->SetStats(kFALSE);
    hF2[j]->SetLineColor(1);
    hF2[j]->SetFillColor(1);
    hF2[j]->SetTitle("full events");
  }

  TH1F* hR[3];
  int hR_col[3] = {2,1,3};
  for(int j=0; j<3; j++){
    sprintf(fname,"hR%d",j);
    hR[j] = new TH1F(fname,"",nBins,T_bgn,T_end);
    hR[j]->SetTitleOffset(1.3,"Y");
    sprintf(fname,"time (%d hours per bin)",int(hours));
    hR[j]->GetXaxis()->SetTitle(fname);
    hR[j]->GetYaxis()->SetTitle("rate, Hz");
    hR[j]->GetYaxis()->SetTitleOffset(2.0);
    hR[j]->SetTitle("full events (black), after pp cuts (green)");
    hR[j]->SetStats(kFALSE);
    hR[j]->SetLineColor(hR_col[j]);
    hR[j]->SetFillColor(hR_col[j]);
  }
 
  // read events
 
  sprintf(fname,"%s/events.txt",netdir);
  FILE* ftrig = fopen(fname,"w");
  fprintf(ftrig,"# correlation threshold = %f \n",T_cor);
  fprintf(ftrig,"# network rho threshold = %f\n",T_cut);
  fprintf(ftrig,"# -/+ - not passed/passed final selection cuts\n");
  fprintf(ftrig,"#  1 - effective correlated amplitude rho \n");
  fprintf(ftrig,"#  2 - correlation coefficient 0/1 (1G/2G)\n");
  fprintf(ftrig,"#  3 - correlation coefficient 2\n");
  fprintf(ftrig,"#  4 - correlation coefficient 3\n");
  fprintf(ftrig,"#  5 - correlated amplitude \n");
  fprintf(ftrig,"#  6 - time shift \n");
  fprintf(ftrig,"#  7 - time super shift \n");  //SLAG
  fprintf(ftrig,"#  8 - likelihood \n");
  fprintf(ftrig,"#  9 - penalty factor \n");
  fprintf(ftrig,"# 10 - energy disbalance \n");
  fprintf(ftrig,"# 11 - central frequency \n");
  fprintf(ftrig,"# 12 - bandwidth \n");
  fprintf(ftrig,"# 13 - duration \n");
  fprintf(ftrig,"# 14 - number of pixels \n");
  fprintf(ftrig,"# 15 - frequency resolution \n");
  fprintf(ftrig,"# 16 - cwb run number \n");
  i = 16;
  for(int j=0; j<nIFO; j++) fprintf(ftrig,"# %2d - time for %s detector \n",++i,IFO[j]);
  for(int j=0; j<nIFO; j++) fprintf(ftrig,"# %2d - sSNR for %s detector \n",++i,IFO[j]);
  for(int j=0; j<nIFO; j++) fprintf(ftrig,"# %2d - hrss for %s detector \n",++i,IFO[j]);
  fprintf(ftrig,"# %2d - PHI \n",++i);
  fprintf(ftrig,"# %2d - THETA \n",++i);
  fprintf(ftrig,"# %2d - PSI \n",++i);
  
  int ntrg = wave.GetEntries();
  cout<<"total events: "<<ntrg<<endl;
  printf("data start GPS: %15.5f \n",sTARt);
  printf("data stop  GPS: %15.5f \n",sTOp);
  rUn = -1;
 
  if(bhveto||bcat3) {
//  int vgtrg = W.fChain->GetEntries();
//  if (vgtrg!=ntrg){cout<<"Warning: vgtree has "<<vgtrg<<" events!"<<endl;}
//  cout << "events passing the " << vetocut << " veto: " << W.fChain->Draw("",vetocut,"goff")<<endl;
  } 

  wavearray<int> LAG_PASS(ntrg);LAG_PASS=0;  	//LAG_PASS ARRAY 
  wavearray<int> SAVE(ntrg);SAVE=0;  		//SAVE ARRAY 
  wavearray<int> Wsel(ntrg);  			//MULTI 
  for(int i=0;i<ntrg;i++) Wsel[i]=0;  		//MULTI 

  cout << "cwb_report_prod_2.C : Start event loop ..." << endl;
  for(int i=0; i<ntrg; i++) {

    W.GetEntry(i);
    if(i%10000==0) cout << "cwb_report_prod_2.C : loop 1 - processed: " << i << "/" << ntrg << endl;

    save_veto=1;

    // apply veto cuts
    bool bveto=false; 
    for(int j=0;j<nXDQF;j++) {
      saveVETO[j]=1; 
      if((XDQF[j].cat==CWB_CAT2)&&(!VETO[j])) {countVETO[j]++;bveto=true;}
      if((XDQF[j].cat==CWB_CAT3)&&(VETO[j]))  {countVETO[j]++;saveVETO[j]=0;save_veto=0;}
      if((XDQF[j].cat==CWB_HVETO)&&(VETO[j])) {countVETO[j]++;saveVETO[j]=0;save_veto=0;}
    }
    if(bveto) continue;

    double WLag = W.lag[nIFO]; 
    double WSLag = W.slag[nIFO]; 

    if((cwb_lag_number==-1)&&(cwb_slag_number==-1)) {
      if((WLag==0)&&(WSLag==0)) continue;  // LAG=0 && SLAG=0
    }
    if((cwb_lag_number>=0)&&(cwb_slag_number==-1)) {
      if(WLag!=cwb_lag_number) continue;  // SINGLE LAG &  ANY SLAG
      if((WLag==0)&&(WSLag==0)) continue;  // LAG=0 && SLAG=0
    }
    if((cwb_lag_number==-1)&&(cwb_slag_number>=0)) {
      if(WSLag!=cwb_slag_number) continue;  // ANY LAG & SINGLE SLAG
      if((WLag==0)&&(WSLag==0)) continue;  // LAG=0 && SLAG=0
    }
    if((cwb_lag_number>=0)&&(cwb_slag_number>=0)) {
      if((WLag!=cwb_lag_number)||(WSLag!=cwb_slag_number)) continue;  // SINGLE LAG & SINGLE SLAG
    }
    LAG_PASS[i]=1;
 
    hrho[0]->Fill(W.rho[pp_irho]);		// rho distribution (no cuts)

    if(W.run != rUn) {
      rUn = int(W.run+0.1);
      hR[0]->Fill(float(W.gps[0]),Trun.data[rUn]);
    }

    hR[1]->Fill(float(W.gps[0]));

    // frequencies user band cuts
    bool bcut=false;
    for(int j=0;j<nFCUT;j++) {
      if((W.frequency[0]>lowFCUT[j])&&(W.frequency[0]<=highFCUT[j])) bcut=true;
    }
    if(bcut) continue;

    vED   = W.neted[0]/W.ecor;               // network disbalance
    xcor  = W.netcc[pp_inetcc];              // network correlation 0/1
    xcor2 = W.netcc[2];                      // network correlation 2
    xcor3 = W.netcc[3];                      // network correlation 3
    rho   = W.rho[pp_irho];                  // effective SNR per detector
    acor  = sqrt(W.ecor/nIFO);

    save = true; 
    if (T_scc>0)    {if(W.netcc[2] < T_scc)     save = false;}
    if (T_pen>0)    {if(W.penalty < T_pen)      save = false;}
    if (T_vED>0)    {if(vED > T_vED)            save = false;} 
    if (T_cor>0)    {if(xcor < T_cor)           save = false;}
    if (T_cut>0)    {if(rho < T_cut)            save = false;}
    if (T_acor>0)   {if(acor < T_acor)          save = false;}
    if (T_mchirp>0) {if(W.chirp[1] < T_mchirp)  save = false;}
    if (T_sphr>0)   {if(W.chirp[3] < T_sphr)    save = false;}
    if (T_efrac>0)  {if(W.chirp[5]*W.chirp[3] + log10(W.size[1])/10< T_efrac)   save = false;}

    SAVE[i]=save;

    if(!save) continue;

    // write to file all events with (rho > T_out) 
    if(rho>T_out) {

      pf_cc->Fill(xcor,W.penalty);
      if(W.frequency[0]>0 && W.duration[0]>0) dens->Fill(-log10(W.size[0]*0.5/W.duration[0]/W.frequency[0]));
      pf_vED->Fill(vED,W.penalty);
      pf_ed->Fill(vED,W.penalty);
      hvED->Fill(fabs(vED));
      hpf->Fill(W.penalty);

      if(save) fprintf(ftrig,"+");
      else     fprintf(ftrig,"-");
      for(int j=0;j<nXDQF;j++) {
        if((XDQF[j].cat==CWB_HVETO)&&(!saveVETO[j])) fprintf(ftrig,"%c",XDQF[j].ifo[0]);
      }
      fprintf(ftrig," -"); 
      for(int j=0;j<nXDQF;j++) {
        if((XDQF[j].cat==CWB_CAT3)&&(!saveVETO[j])) fprintf(ftrig,"%c",XDQF[j].ifo[0]);
      }
      fprintf(ftrig," "); 

      int wslag=0; wslag=W.slag[nIFO]; 
      fprintf(ftrig,"%5.4f %4.3f %4.3f %4.3f %5.2f %7.3f %7d %5.1e %4.3f %4.3f ",
	      rho,xcor,xcor2,xcor3,acor,W.lag[nIFO],wslag,W.likelihood,W.penalty,vED);  
      fprintf(ftrig,"%4.0f %4.0f %4.3f %3d %3.1d %5d ",
	      W.frequency[0],W.bandwidth[0],
	      W.duration[0],W.size[0],W.rate[0],int(W.run));

      for(j=0; j<nIFO; j++) fprintf(ftrig,"%14.4f ",W.time[j]); 
      for(j=0; j<nIFO; j++) fprintf(ftrig,"%5.1e ",W.sSNR[j]); 
      for(j=0; j<nIFO; j++) fprintf(ftrig,"%5.1e ",W.hrss[j]);
      if(TMath::IsNaN(W.psi[0])) W.psi[0]=-1000;	// fix psi value when is nan
      fprintf(ftrig,"%4.2f %4.2f %4.2f ",W.phi[0], W.theta[0], W.psi[0]);
      fprintf(ftrig,"\n");
    }

    if(rho>T_out) {
      int indx = lagslag[std::make_pair(int(W.slag[nIFO]),int(W.lag[nIFO]))];
      Rlag.data[indx] += 1.;
    }

    int rho_size=0;    
    double drho = (pp_rho_max-pp_rho_min)/pp_rho_bin;
    for(j=0; j<Xcut.size(); j++) {
      Xcut.data[j] = pp_rho_min+j*drho;
      if(Xcut.data[j]>pp_rho_max) continue;
      if(rho > Xcut.data[j]) {Ycut.data[j] += 1.;}  
      if(rho > Xcut.data[j] && save_veto==1) {Ycut_veto.data[j] += 1.;}
      if(rho > Xcut.data[j]) {y_lag[iy_lag[int(W.lag[nIFO])-min_lag]][j] += 1.;}
      if(rho > Xcut.data[j] && save_veto==1) {y_lag_veto[iy_lag[int(W.lag[nIFO])-min_lag]][j] += 1.;}
      rho_size++;
    }
    Xcut.resize(rho_size);
 
    // fill high resolution array (drho=0.01) 'dump to far_rho.txt file' with not vetoed bkg events
    if(save_veto) AddRho2FAR(rho, &Xfar, &Yfar, far_rho_min, far_drho);	

    Wsel[i]=1;  //MULTI

  } // End event loop

  delete [] iy_lag;
  for(int i=0;i<nlag;i++) {
    delete [] y_lag[i]; 
    delete [] y_lag_veto[i]; 
  }

  // WARNING !!! the main loop has been splitted to make faster the loop. I don't know why!!!
  for(int i=0; i<ntrg; i++) {

    W.GetEntry(i);
    if(i%10000==0) cout << "cwb_report_prod_2.C : loop 2 - processed: " << i << "/" << ntrg << endl;

    if (!LAG_PASS[i]) continue;			 // skip events if not pass lag/slag cuts

    //ONLINE
    sm_theta_phi.add(sm_theta_phi.getSkyIndex(W.theta[0],W.phi[0]),1);
    sm_ra_dec.add(sm_ra_dec.getSkyIndex(-(W.theta[2]-90.),W.phi[2]),1);

    xcor = W.netcc[pp_inetcc];                   // network correlation
    rho  = W.rho[pp_irho];                       // effective SNR per detector
    vED  = W.neted[0]/W.ecor;   		 // energy disbalance

    // get dectector with max SNR
    double mSNR = 0.;
    int detId = 0;
    for(j=0; j<nIFO; j++) {
      if(W.snr[j]>mSNR) {detId=j+1; mSNR=W.snr[j];}
    }

    // fill histograms with events which pass pp cuts

    rho_L->Fill(log10(2*rho*rho),log10(W.likelihood/nIFO));
    xCor->Fill(xcor);
    rhocc->Fill(xcor,rho);
    ecor_cc->Fill(xcor,W.rho[0]);
    ECOR_cc->Fill(xcor,W.rho[1]);
    rho_subnet->Fill(W.netcc[2],rho);
    rho_vED->Fill(vED,rho);  
    float chi2 = W.penalty>0 ? log10(W.penalty) : 0; 
    if(TString(analysis)=="1G") rho_pf->Fill(W.penalty,rho);  
    else			rho_pf->Fill(chi2,rho);

    float Wgps=float(W.gps[0]);
    float Wfreq=float(W.frequency[0]);

    hF[0]->Fill(Wgps); 
    hF[detId]->Fill(Wgps,100.);  
    hF2[0]->Fill(Wfreq); 
    hF2[detId]->Fill(Wfreq,100.);  
    cutF->Fill(Wfreq); 
   
    ecor->Fill(sqrt(W.ecor/nIFO));
    ECOR->Fill(sqrt(W.ECOR/nIFO));
    Like->Fill(log10(W.likelihood));
    
    Mchirp->Fill(W.chirp[1]);
    McErr->Fill(W.chirp[5]);
    rho_mchirp->Fill(W.chirp[1],rho);

    if (!SAVE[i]) continue;			 // skip events if not pass pp cuts

    hR[2]->Fill(Wgps);
    hrho[1]->Fill(rho);			

  } // End event loop

  fclose(ftrig);

  for(int j=0;j<nXDQF;j++) {
    cout << "Events vetoed by " << nameVETO[j].Data() << " - " << countVETO[j] << endl;
  }

  sprintf(fname,"%s/lags.txt",netdir);
  FILE* filelags = fopen(fname,"w");
  for(j=0; j<Wlag[nIFO].size(); j++){
    Elag.data[j] = sqrt(Rlag.data[j])/Tlag.data[j];
    if(Wlag[nIFO].data[j]>0.) {
      fprintf(filelags,"%6.3f %6.3f %d %d\n",
              Wslag[nIFO].data[j],Wlag[nIFO].data[j],(int)Tlag.data[j],(int)Rlag.data[j]);
    } 
    Rlag.data[j]/= Tlag.data[j];
  }
  fclose(filelags);

  if(pp_rho_vs_rate_no_multiplicity) cout << "start rate vs mutiplicity ..." << endl;
  for(j=0; j<Xcut.size(); j++){
    if(pp_rho_vs_rate_no_multiplicity) {
      // rate correct taking into account the multiplicity
      wavecomplex rate_multi=CWB::Toolbox::getRate(Xcut.data[j],Tgap,nIFO,wave,Wsel,Wlag,Wslag,Tlag);  //MULTI 
      Ycut_multi.data[j] = rate_multi.real();
      Yerr_multi.data[j] = rate_multi.imag();
    }
    Yerr.data[j] = sqrt(Ycut.data[j])/liveTot;
    Ycut.data[j] = Ycut.data[j]/liveTot;
    yERR.data[j] = sqrt(yCUT.data[j])/liveTot;
    yCUT.data[j] = yCUT.data[j]/liveTot;
    //veto
    Yerr_veto.data[j] = sqrt(Ycut_veto.data[j])/liveTot;
    Ycut_veto.data[j] = Ycut_veto.data[j]/liveTot;
#ifdef PERC_CAT3
    Ycut_veto.data[j] /= PERC_CAT3;
#endif
  }

  // write rho statistic distribution
  ofstream out_far;
  sprintf(fname,"%s/far_rho.txt",netdir);
  out_far.open(fname,ios::out);
  if(!out_far.good()) {cout << "Error Opening File " << fname << endl;exit(1);}
  ofstream out_efar;
  sprintf(fname,"%s/efar_rho.txt",netdir);
  out_efar.open(fname,ios::out);
  if(!out_efar.good()) {cout << "Error Opening File " << fname << endl;exit(1);}
  for(int i=0;i<Xfar.size(); i++) if (Yfar[i]>0) {
    double x  = Xfar[i];
    double ex = 0;
    double y  = Yfar[i]/liveTot;
    double ey = sqrt(Yfar[i])/liveTot;
    if(y>=1) out_far << x << "\t\t" << y << endl;
    if(y<1)  out_far << x << "\t\t" << y << endl;
    //cout << i << " " << x << " " << y << " " << ex << " " << ey << endl;
    if(y>=1) out_efar << x << "\t\t" << y << "\t\t" << ex << "\t" << ey << endl;
    if(y<1)  out_efar << x << "\t\t" << y << "\t" << ex << "\t" << ey << endl;
  }
  out_far.close();
  out_efar.close();

  sprintf(fname,"%s/rate_threshold.txt",netdir);
  FILE* frate = fopen(fname,"w");
  for(j=0; j<Xcut.size(); j++){
    fprintf(frate,"%3.2f %4.3e\n",Xcut.data[j],Ycut.data[j]);
  }
  fclose(frate);

  sprintf(fname,"%s/rate_threshold_veto.txt",netdir);
  FILE* fratev = fopen(fname,"w");
  for(j=0; j<Xcut.size(); j++){
    fprintf(fratev,"%3.2f %4.3e\n",Xcut.data[j],Ycut_veto.data[j]);
  }
  fclose(fratev);

  //MULTI
  sprintf(fname,"%s/rate_threshold_multi.txt",netdir);
  FILE* fratem= fopen(fname,"w");
  for(j=0; j<Xcut.size(); j++){
    fprintf(fratem,"%3.2f %4.3e\n",Xcut.data[j],Ycut_multi.data[j]);
  }
  fclose(fratem);

  char output[256];
  sprintf(fname,"%s/sigma.txt",netdir);
  ofstream out_lf(fname,ios::out);
  for(j=0; j<pp_rho_bin; j++){
    int elag=0;
    double mean=0;
    double sigma=0;
    double sum=0;
    for (int k=1; k<nlag; k++) {   // FIX 
      float rate=0.0;
      if(Tlag[k]>0) rate=y_lag[k][j]/Tlag[k];   // FIX Tlag[j] -> Tlag[k]
      mean+=rate;
      sum+=rate*rate;
      elag++;
    }
    if(elag>0) {
      mean=mean/elag; 
      sigma= sqrt((sum-elag*mean*mean)/elag);
    }
    sprintf(output,"%3.2f %4.3e %4.3e\n",Xcut[j],Ycut[j],sigma);
    out_lf << output;
  }
  cout << fname << endl;
  out_lf.close();

  sprintf(fname,"%s/sigma_veto.txt",netdir);
  ofstream out_veto(fname,ios::out);
  for(j=0; j<pp_rho_bin; j++){
    int elag=0;
    double mean=0;
    double sigma=0;
    double sum=0;
    for (int k=1; k<nlag-1; k++) {   // FIX
      float rate=0.0;
      if(Tlag[k]>0) rate=y_lag_veto[k][j]/Tlag[k];    // FIX : Tlag[j] -> Tlag[k]
      mean+=rate;
      sum+=rate*rate;
      elag++;
    }
    if(elag>0) {
      mean=mean/elag;
      sigma= sqrt((sum-elag*mean*mean)/elag);
    }
    sprintf(output,"%3.2f %4.3e %4.3e\n",Xcut[j],Ycut_veto[j],sigma);
    out_veto << output;
  }
  cout << fname << endl;
  out_veto.close();

  CWB::Toolbox::doPoissonPlot(nIFO, Wlag, Tlag, Rlag, netdir);
  
  c1->Clear();
  pf_vED->Draw("colz");
  sprintf(fname,"%s/penalty_vED.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  c1->Clear();
  rho_L->Draw("colz");
  sprintf(fname,"%s/rho_L.eps",netdir);
  c1->Update(); c1->SaveAs(fname); 

  c1->Clear();
  pf_ed->Draw("colz");
  sprintf(fname,"%s/penalty_ed.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
 
  c1->Clear();
  pf_cc->Draw("colz");
  sprintf(fname,"%s/penalty_cc.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  if(c1) delete c1;
  c1 = new TCanvas("c","C",0,0,1400,400*nIFO);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetBottomMargin(0.143939);
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);

  c1->Clear();
  c1->Divide(1,nIFO);
  for(i=1; i<nIFO+1; i++) { 
    c1->cd(i);
    c1->SetLogy(kFALSE);
    hF[i]->Divide(hF[0]); hF[i]->SetMaximum(100.);
    gStyle->SetOptStat(kFALSE);
    gPad->SetBottomMargin(0.2);
    hF[i]->Draw("HIST"); 
    hF[i]->GetYaxis()->SetLabelSize(0.06);
    hF[i]->GetXaxis()->SetTimeDisplay(1);
  }
  c1->Update();
  sprintf(fname,"%s/fraction_time.eps",netdir);
  c1->SaveAs(fname);

  c1->Clear();
  c1->Divide(1,nIFO);
  for(i=1; i<nIFO+1; i++) { 
    c1->cd(i);
    c1->SetLogy(kFALSE);
    hF2[i]->Divide(hF2[0]); hF2[i]->SetMaximum(100.);
    gStyle->SetOptStat(kFALSE);
    gPad->SetBottomMargin(0.2);
    hF2[i]->Draw("HIST"); 
    hF2[i]->GetYaxis()->SetLabelSize(0.06);
  }
  c1->Update();
  sprintf(fname,"%s/fraction_frequency.eps",netdir);
  c1->SaveAs(fname);

  c1->Clear();
  cutF->Draw("HIST");
  c1->Update();
  sprintf(fname,"%s/frequency.eps",netdir);
  c1->SaveAs(fname);

  if(c1) delete c1; 
  c1 = new TCanvas("c","C",0,0,1400,1050);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetBottomMargin(0.143939);
  //c1->SetRightMargin(0.1517039);
  c1->SetRightMargin(0.10);
  c1->SetLeftMargin(0.15);
  c1->SetTopMargin(0.0772727);

  c1->Clear();
  c1->SetLogy(kTRUE);
  hR[1]->Divide(hR[0]);
  hR[1]->Draw("HIST");
  hR[1]->SetMinimum(1.e-5);
  gStyle->SetOptStat(kFALSE);
  hR[1]->Draw("HIST");
  hR[1]->GetXaxis()->SetTimeDisplay(1);
  hR[1]->GetXaxis()->SetLabelSize(0.04);
  for(i=2; i<3; i++) { hR[i]->Divide(hR[0]); hR[i]->Draw("same"); }
  sprintf(fname,"%s/rate_time.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  c1->Clear();
  c1->SetLogy(kFALSE);
  TGraphErrors* gr2;
  gr2=new TGraphErrors(Wlag[nIFO].size(),Tdlag.data,Rlag.data,Olag.data,Elag.data);
  sprintf(_title,"after pp-cuts & rho > %3.2f",T_out);
  gr2->SetTitle(_title);
  gr2->SetLineColor(1);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(1);
  gr2->SetMarkerSize(1);
  gr2->Draw("AP");
  gr2->GetHistogram()->SetXTitle("lag");
  gr2->GetHistogram()->SetYTitle("rate, Hz");
  gr2->GetHistogram()->GetXaxis()->SetNdivisions(508);

  sprintf(fname,"%s/rate_lag.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  c1->Clear();
  c1->SetLogy(kFALSE);
  for(int l=0;l<Wlag[nIFO].size();l++) Tlag.data[l]/=86400.;
  if(gr2) delete gr2;
  gr2=new TGraphErrors(Wlag[nIFO].size(),Tdlag.data,Tlag.data,Olag.data,Elag.data);
  gr2->SetTitle("live vs lag");
  gr2->SetLineColor(1);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(1);
  gr2->SetMarkerSize(1);
  gr2->Draw("AP");
  gr2->GetHistogram()->GetXaxis()->SetNdivisions(508);
  gr2->GetHistogram()->SetXTitle("lag");
  gr2->GetHistogram()->SetYTitle("live, days");


  sprintf(fname,"%s/live_lag.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  #define YEAR (24.*3600.*365.)

  c1->Clear();
  c1->SetLogy(kTRUE);
  if(pp_rho_log) c1->SetLogx(kTRUE); else c1->SetLogx(kFALSE);
  wavearray<double> Ycut_year = Ycut;
  wavearray<double> Yerr_year = Yerr;
  Ycut_year *= YEAR;
  Yerr_year *= YEAR;
  TGraphErrors* gr1;
  gr1=new TGraphErrors(Xcut.size(),Xcut.data,Ycut_year.data,Xerr.data,Yerr_year.data);
  gr1->SetTitle("after pp-cuts");
  gr1->SetLineColor(1);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(1);
  gr1->SetMarkerSize(1);
  gr1->SetMinimum(0.5/liveTot*YEAR);
  c1->SetLogy(kTRUE);
  gr1->Draw("AP");
  gr1->GetHistogram()->SetXTitle("#rho");
  gr1->GetHistogram()->SetYTitle("rate, 1/year");

  sprintf(fname,"%s/rate_threshold.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  if(bhveto||bcat3) {
    c1->Clear();
    c1->SetLogy(kTRUE);
    if(pp_rho_log) c1->SetLogx(kTRUE); else c1->SetLogx(kFALSE);
    wavearray<double> Ycut_veto_year = Ycut_veto;
    wavearray<double> Yerr_veto_year = Yerr_veto;
    Ycut_veto_year *= YEAR;
    Yerr_veto_year *= YEAR;
    TGraphErrors* gr1_veto;
    gr1_veto=new TGraphErrors(Xcut.size(),Xcut.data,Ycut_veto_year.data,Xerr.data,Yerr_veto_year.data);
    gr1_veto->SetTitle("after pp-cuts & vetoes");
    gr1_veto->SetLineColor(1);
    gr1_veto->SetMarkerStyle(20);
    gr1_veto->SetMarkerColor(1);
    gr1_veto->SetMarkerSize(1);
    gr1_veto->SetMinimum(0.5/liveTot*YEAR);
    c1->SetLogy(kTRUE);
    gr1_veto->Draw("AP");
    gr1_veto->GetHistogram()->SetXTitle("#rho");
    gr1_veto->GetHistogram()->SetYTitle("rate, 1/year");
  
    sprintf(fname,"%s/rate_threshold_veto.eps",netdir);
    c1->Update(); c1->SaveAs(fname);
  
    c1->Clear();
    if(pp_rho_log) c1->SetLogx(kTRUE); else c1->SetLogx(kFALSE);
    gr1->SetTitle("after pp-cuts (black), after pp-cuts & vetoes (red) ");
    gr1->Draw("AP");
    gr1_veto->SetLineColor(2);
    gr1_veto->SetMarkerColor(2);
    gr1_veto->Draw("P,same");
    sprintf(fname,"%s/rate_threshold_veto_compare.eps",netdir);
    c1->Update(); c1->SaveAs(fname);
  }

  //MULTI
  if(pp_rho_vs_rate_no_multiplicity) {
    c1->Clear();
    c1->SetLogy(kTRUE);
    if(pp_rho_log) c1->SetLogx(kTRUE); else c1->SetLogx(kFALSE);
    wavearray<double> Ycut_multi_year = Ycut_multi;
    wavearray<double> Yerr_multi_year = Yerr_multi;
    Ycut_multi_year *= YEAR;
    Yerr_multi_year *= YEAR;
    TGraphErrors* gr1_multi;
    gr1_multi=new TGraphErrors(Xcut.size(),Xcut.data,Ycut_multi_year.data,Xerr.data,Yerr_multi_year.data);
    gr1_multi->SetLineColor(1);
    gr1_multi->SetMarkerStyle(20);
    gr1_multi->SetMarkerColor(1);
    gr1_multi->SetMarkerSize(1);
    gr1_multi->SetMinimum(0.5/liveTot);
    c1->SetLogy(kTRUE);
    gr1_multi->Draw("AP");
    gr1_multi->GetHistogram()->SetXTitle("#rho");
    gr1_multi->GetHistogram()->SetYTitle("rate, Hz");

    sprintf(fname,"%s/rate_threshold_multi.eps",netdir);
    c1->Update(); c1->SaveAs(fname);
  
    c1->Clear();
    if(pp_rho_log) c1->SetLogx(kTRUE); else c1->SetLogx(kFALSE);
    gr1->Draw("AP");
    gr1->SetTitle("after pp-cuts (black), no-multiplicity (green) ");
    gr1_multi->SetLineColor(kGreen+2);
    gr1_multi->SetMarkerColor(kGreen+2);
    gr1_multi->Draw("P,same");
    sprintf(fname,"%s/rate_threshold_multi_compare.eps",netdir);
    c1->Update(); c1->SaveAs(fname);
  }

  c1->Clear();
  c1->SetLogy(kFALSE);
  dens->Draw("HIST");
  sprintf(fname,"%s/density.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  c1->Clear();
  c1->SetLogy(kTRUE);
  hpf->Draw("HIST");
  sprintf(fname,"%s/penalty.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  c1->Clear();
  c1->SetLogy(kTRUE);
  hvED->Draw("HIST");
  sprintf(fname,"%s/vED.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  gStyle->SetOptStat(kTRUE);
  c1->Clear();
  c1->SetLogy(kTRUE);
  xCor->Draw("HIST");
  sprintf(fname,"%s/netcor.eps",netdir);
  c1->Update(); c1->SaveAs(fname);

  c1->Clear();
  if(pp_rho_log) c1->SetLogx(kTRUE); else c1->SetLogx(kFALSE);
  ecor->Draw("HIST");
  sprintf(fname,"%s/ecor.eps",netdir);
  c1->Update();c1->SaveAs(fname);
  
  c1->Clear();
  hrho[0]->Draw("HIST");
  hrho[1]->Draw("same");
  sprintf(fname,"%s/rho.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
#ifdef WRITE_ASCII
  gROOT->LoadMacro(wat_dir+"/tools/cwb/postproduction/burst/h12ascii.C");
  sprintf(fname,"%s/rho.txt",netdir);
  h12ascii(hrho[0],fname);
  sprintf(fname,"%s/rho_veto.txt",netdir);
  h12ascii(hrho[1],fname);
#endif 
  
  c1->Clear();
  ECOR->Draw("HIST");
  sprintf(fname,"%s/ECOR.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  Like->Draw("HIST");
  sprintf(fname,"%s/likelihood.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  Mchirp->Draw("HIST");
  sprintf(fname,"%s/Mchirp.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  c1->SetLogy(kFALSE);
  McErr->Draw("HIST");
  sprintf(fname,"%s/McErr.eps",netdir);
  c1->Update(); c1->SaveAs(fname);    
  
  c1->Clear();
  c1->SetLogx(kFALSE);
  if(pp_rho_log) c1->SetLogy(kTRUE); else c1->SetLogy(kFALSE);
  rho_subnet->Draw("colz");
  sprintf(fname,"%s/rho_subnet.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  c1->SetLogx(kFALSE);
  if(pp_rho_log) c1->SetLogy(kTRUE); else c1->SetLogy(kFALSE);
  rho_vED->Draw("colz");
  sprintf(fname,"%s/rho_vED.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  c1->SetLogx(kFALSE);
  if(pp_rho_log) c1->SetLogy(kTRUE); else c1->SetLogy(kFALSE);
  rho_pf->Draw("colz");
  sprintf(fname,"%s/rho_pf.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  c1->SetLogx(kFALSE);
  if(pp_rho_log) c1->SetLogy(kTRUE); else c1->SetLogy(kFALSE);
  rhocc->Draw("colz");
  sprintf(fname,"%s/rho_cc.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  ECOR_cc->Draw("HIST");
  sprintf(fname,"%s/ECOR_cc.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  ecor_cc->Draw("HIST");
  sprintf(fname,"%s/ecor_cc.eps",netdir);
  c1->Update(); c1->SaveAs(fname);
  
  c1->Clear();
  c1->SetLogx(kTRUE);
  if(pp_rho_log) c1->SetLogy(kTRUE); else c1->SetLogy(kFALSE);
  rho_mchirp->Draw("colz");
  sprintf(fname,"%s/rho_mchirp.eps",netdir);
  c1->Update(); c1->SaveAs(fname);  
  c1->SetLogx(kFALSE);

  // WARNING: The following declarations must be creted with new otherwise in ROOT6 macro crash when exit !!!

  //ONLINE
  gskymap* gsm_theta_phi = new gskymap(sm_theta_phi);
  gsm_theta_phi->SetTitle("Latitude - Longitude");
  gsm_theta_phi->SetWorldMap();
  gsm_theta_phi->Draw();
  sprintf(fname,"%s/phi_vs_theta_online.eps",netdir);
  gsm_theta_phi->GetCanvas()->SaveAs(fname);
  delete gsm_theta_phi;

  gskymap* gsm_ra_dec = new gskymap(sm_ra_dec);
  gsm_ra_dec->SetTitle("right ascension - declination");
  //gsm_ra_dec->SetOptions("hammer","celestial");
  gsm_ra_dec->SetGalacticDisk();
  gsm_ra_dec->Draw();
  gsm_ra_dec->GetHistogram()->GetXaxis()->SetTitle("right ascension, degree");
  gsm_ra_dec->GetHistogram()->GetYaxis()->SetTitle("declination, degree");
  sprintf(fname,"%s/ra_vs_dec_online.eps",netdir);
  gsm_ra_dec->GetCanvas()->SaveAs(fname);
  delete gsm_ra_dec;

  gSystem->Exit(0);  
}
