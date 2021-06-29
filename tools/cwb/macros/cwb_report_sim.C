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
  #define EPSILON  0.001
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

  // initialize the ifo name (takes into account the non built-in detectors)
  char IFO[NIFO_MAX][32];
  for(int n=0;n<nIFO;n++) {
    if(strlen(ifo[n])!=0) strcpy(IFO[n],ifo[n]); else strcpy(IFO[n],detParms[n].name);
  }

  CWB::Toolbox::mkDir(netdir,!pp_batch);

  gStyle->SetTitleOffset(1.4,"X");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelOffset(0.014,"X");
  gStyle->SetLabelOffset(0.010,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetLabelSize(0.03,"X");
  gStyle->SetLabelSize(0.03,"Y");

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  //gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);
//  gStyle->SetOptStat(kFALSE);

  gStyle->SetStatX(0.878);
  gStyle->SetStatY(0.918);
  gStyle->SetStatW(0.200);
  gStyle->SetStatH(0.160);
//  gStyle->SetStatBorderSize(1);
  gStyle->SetStatColor(0);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  cout<<"cwb_report_sim.C starts..."<<endl;


  Color_t colors[NMDC_MAX] = { 6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7,
                               6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7,
                               6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7,
                               6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7 };
  Style_t markers[NMDC_MAX]= {20,21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,
                              21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,20, 
                              21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,20, 
                              21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,20 };
  // histograms
  TH2F* HH[NIFO_MAX][NMDC_MAX];
  TH2F* NRE[NIFO_MAX][NMDC_MAX];
  TH1F* hT[NMDC_MAX];
  TH1F* hF[NMDC_MAX];
  TGraphErrors* EFF[NMDC_MAX];

  for(int i=0;i<NMDC_MAX;i++) {
    for(int j=0;j<nIFO;j++) HH[j][i]=NULL;
    for(int j=0;j<nIFO;j++) NRE[j][i]=NULL;
    hT[i]=NULL;
    hF[i]=NULL;
    EFF[i]=NULL;
  }

  // transform T_ifar from years into sec
  T_ifar*=(365.*24.*3600.); 

  // XDQF structure used for internal's computations
  dqfile* XDQF = new dqfile[nVDQF];

  // read injection file types
  char    imdc_set[NMDC_MAX][128];    // injection set
  size_t  imdc_type[NMDC_MAX];        // injection type
  char    imdc_name[NMDC_MAX][128];   // injection name
  double  imdc_flow[NMDC_MAX];        // injection low frequencies
  double  imdc_fhigh[NMDC_MAX];       // injection high frequencies
  double  imdc_fcentral[NMDC_MAX];    // injection low frequencies
  double  imdc_fbandwidth[NMDC_MAX];  // injection high frequencies
  size_t  imdc_index[NMDC_MAX];       // type reference array
  size_t  imdc_iset[NMDC_MAX];        // injection set index
  size_t  imdc_tset[NMDC_MAX];        // injection set type index

  int ninj=ReadInjType(mdc_inj_file,NMDC_MAX,imdc_set,imdc_type,
                       imdc_name,imdc_fcentral,imdc_fbandwidth);
  if(ninj==0) {
    cout << endl;
    cout << "cwb_report_sim.C : Error - no injection types or bad format" << endl;
    cout << "mdc injection file list : " << mdc_inj_file << endl; 
    cout << "format must be : mdc_set mdc_type mdc_name mdc_fcentral mdc_fbandwidth" << endl;  
    cout << "process terminated !!!" << endl;
    exit(1);
  } 

  for(int i=0; i<NMDC_MAX; i++) {
    imdc_flow[i]  = imdc_fcentral[i]-imdc_fbandwidth[i];
    if(imdc_flow[i]<fLow)  imdc_flow[i]=fLow;
    if(imdc_flow[i]>fHigh) imdc_flow[i]=fLow;
    imdc_fhigh[i] = imdc_fcentral[i]+imdc_fbandwidth[i];
    if(imdc_fhigh[i]<fLow)  imdc_fhigh[i]=fHigh;
    if(imdc_fhigh[i]>fHigh) imdc_fhigh[i]=fHigh;
  }

  for(int i=0; i<NMDC_MAX; i++) imdc_index[i] = NMDC_MAX;
  for(int j=0;j<ninj;j++) imdc_index[imdc_type[j]]=j;
  for(int j=0;j<ninj;j++) {
    cout << j << " " << imdc_index[j] << endl;
    if(imdc_index[j]>=ninj) {
      cout << endl;
      cout << "cwb_report_sim.C : Error - mdc type must be < num max type inj = " << ninj << endl;
      cout << "check mdc_type injection file list : " << mdc_inj_file << endl; 
      cout << "format : mdc_set mdc_type mdc_name mdc_fcentral mdc_fbandwidth" << endl;  
      cout << "process terminated !!!" << endl;
      exit(1);
    }
  }

  TString* imdc_set_name = new TString[ninj];
  int nset=0;
  for(int i=0;i<ninj;i++) {
    bool bnew=true;
    for(int j=0;j<nset;j++) if(imdc_set[i]==imdc_set_name[j]) bnew=false;
    if(bnew) imdc_set_name[nset++]=imdc_set[i];
  }
  cout << "nset : " << nset << endl;

  for(int i=0;i<nset;i++) {
    for(int j=0;j<ninj;j++) if(imdc_set[j]==imdc_set_name[i]) imdc_iset[j]=i;
  }

  for(int i=0;i<ninj;i++) {
    for(int j=0;j<nset;j++) if(imdc_set[i]==imdc_set_name[j]) imdc_tset[imdc_type[i]]=j;
  }

  for(int i=0;i<ninj;i++) {
    cout << i << " " << imdc_name[i] << " " << imdc_set[i] << " " << imdc_iset[i] 
         << " " << imdc_tset[i] << " " << imdc_type[i] << " " << imdc_index[i] << endl;
  }

  // compute efficiency vs threshold for each mdc type (used for ROC)
  TString pp_eff_vs_thr_mode = CWB::Toolbox::getParameter(pp_eff_vs_thr,"--mode");
  if((pp_eff_vs_thr_mode!="DISABLED")&&(pp_eff_vs_thr_mode!="disabled")) {

    TString pp_eff_vs_thr_eff = CWB::Toolbox::getParameter(pp_eff_vs_thr,"--eff");
    int pp_eff = 50;
    if(pp_eff_vs_thr_eff.IsFloat()) pp_eff = pp_eff_vs_thr_eff.Atoi();
    if(pp_eff<0) pp_eff=0;
    if(pp_eff>100) pp_eff=100;

    TChain sim("waveburst");
    TChain mdc("mdc");
    sim.Add(sim_file_name);
    mdc.Add(mdc_file_name);

    char fname[1024];
    wavearray<double> factorXX(ninj); factorXX=0;

    // if mode=FILE_NAME -> read factorXX from custom input file
    if((pp_eff_vs_thr_mode!="AUTO")&&(pp_eff_vs_thr_mode!="auto")) {
      sprintf(fname,"%s",pp_eff_vs_thr_mode.Data());
      //cout << fname << endl;
      FILE* fp = fopen(fname,"r");
      if(fp==NULL) {
        cout << "cwb_report_sim.C - Error opening file : " << fname << endl;
        gSystem->Exit(1);
      }
      char xmdc_name[256];
      int  xmdc_type;
      float xfactors;  
      for(int i=0; i<ninj; i++) {
        fscanf(fp,"%s %d %f",xmdc_name, &xmdc_type, &xfactors);
        cout << i << " " << xmdc_name <<" "<<  xmdc_type <<" "<< xfactors << endl;
        for(int j=0;j<ninj;j++) {
          if(TString(xmdc_name).Contains(imdc_name[j])) {factorXX[j]=xfactors;strcpy(imdc_name[j],xmdc_name);}
          if(TString(imdc_name[j]).Contains(xmdc_name)) {factorXX[j]=xfactors;strcpy(imdc_name[j],xmdc_name);}
        }
      }
      fclose(fp);
      // check if input parameters are consistent to the current ones
      for(int i=0; i<ninj; i++) {
        if(factorXX[i]==0) {
          cout << "cwb_report_sim.C - missing factor for : " << imdc_name[i] << endl << endl;
          gSystem->Exit(1);
        } 
      }
    }

    for(int i=0; i<ninj; i++) {

      int  nmdc,nsim;
      char sim_cuts[1024];
      char mdc_cuts[1024];
      char fac_cuts[1024];
      char rho_cuts[1024];

      // if mode=AUTO -> select the nearest factor to %XX efficiency
      if((pp_eff_vs_thr_mode=="AUTO")||(pp_eff_vs_thr_mode=="auto")) {

        double effXX=0;
        for(int k=0;k<nfactor;k++) {
          sprintf(rho_cuts,"rho[%d]>%f",pp_irho,T_cut);
          sprintf(fac_cuts,"abs(factor-%f)<EPSILON*%f",factors[k],factors[k]);
          if(pp_merge_types) {
            sprintf(sim_cuts,"%s && %s && %s",ch2,fac_cuts,rho_cuts);
          } else {
            sprintf(mdc_cuts,"type[0]==%d && %s",(int)(imdc_type[i]+1),fac_cuts);
            sprintf(sim_cuts,"%s && type[1]==%lu && %s && %s",ch2,imdc_type[i]+1,fac_cuts,rho_cuts);
          }
          //cout << "mdc_cuts " << mdc_cuts << endl;
          //cout << "sim_cuts " << sim_cuts << endl;
          mdc.Draw("Entry$",mdc_cuts,"goff");
          nmdc = mdc.GetSelectedRows();
          sim.Draw("Entry$",sim_cuts,"goff");
          nsim = sim.GetSelectedRows();
          double xeff = double(nsim)/double(nmdc);
          if(fabs(100*xeff-pp_eff)<fabs(100*effXX-pp_eff)) {effXX=xeff;factorXX[i]=factors[k];}
          //cout << "factor : " << factors[k] << " eff : " << nsim << "/" << nmdc << " % " << xeff << endl;
        }
        //cout << "factorXX : " << factorXX[i] << endl;
      } 

      cout << i << "\tCompute efficiency vs threshold @ factor : " 
           << factorXX[i] << "\tmdc : " <<  imdc_name[i] << endl;

      // compute mdc injections @ factorXX
      sprintf(fac_cuts,"abs(factor-%f)<EPSILON*%f",factorXX[i],factorXX[i]);
      if(pp_merge_types) {
      } else {
        sprintf(mdc_cuts,"type[0]==%d && %s",(int)(imdc_type[i]+1),fac_cuts);
      }
      //cout << "mdc_cuts " << mdc_cuts << endl;
      mdc.Draw("Entry$",mdc_cuts,"goff");
      nmdc = mdc.GetSelectedRows();
      //cout << "mdc size : " << nmdc << endl;

      // get the minimum rho for the efficiency computation
      float rho_min=0;
      TString pp_eff_vs_thr_rho_min = CWB::Toolbox::getParameter(pp_eff_vs_thr,"--rho_min");
      if(pp_eff_vs_thr_rho_min.IsFloat()) rho_min = pp_eff_vs_thr_rho_min.Atof();

      // compute efficiency vs threshold @ factorXX
      sprintf(fname,"%s/eff_%d_threshold_%s.txt",netdir,pp_eff,imdc_name[i]);
      cout << fname << endl;
      ofstream fout;
      fout.open(fname, ios::out);
      if (!fout.good()) {cout << "Error Opening File : " << fname << endl;exit(1);}
      int rho_size=0;
      double drho = (pp_rho_max-pp_rho_min)/pp_rho_bin;
      for(int j=0; j<pp_rho_bin; j++) {
        double rho = pp_rho_min+j*drho;
        double xeff = 0;
        if(rho>rho_min) {
          if(pp_merge_types) {
            sprintf(sim_cuts,"%s && %s && rho[%d]>%f", ch2,fac_cuts,pp_irho,rho);
          } else {
            sprintf(sim_cuts,"%s && type[1]==%lu && %s && rho[%d]>%f", ch2,imdc_type[i]+1,fac_cuts,pp_irho,rho);
          }
          //cout << "sim_cuts " << sim_cuts << endl;
          sim.Draw("Entry$",sim_cuts,"goff");
          nsim = sim.GetSelectedRows();
          double xeff = double(nsim)/double(nmdc);
        } 
        //cout << "rho : " << rho << " eff : " << nsim << "/" << nmdc << " % " << xeff << endl;
        fout << rho << " " << xeff << endl;
      }
      fout.close();
    }

    sprintf(fname,"%s/eff_%d_threshold_factors.txt",netdir,pp_eff);
    cout << fname << endl;
    FILE* fp = fopen(fname,"w");
    if(fp==NULL) {
      cout << "cwb_report_sim.C - Error opening file : " << fname << endl;
      gSystem->Exit(1);
    }
    for(int i=0; i<ninj; i++) {
      fprintf(fp,"%s %d %g\n",imdc_name[i], (int)imdc_type[i], factorXX[i]);
    }
    fclose(fp);

    TString pp_eff_vs_thr_quit = CWB::Toolbox::getParameter(pp_eff_vs_thr,"--quit");
    pp_eff_vs_thr_quit.ToUpper();
    if(pp_eff_vs_thr_quit=="TRUE") {
      cout << "efficiency " << pp_eff << " vs threshold terminated, exit" << endl; 
      gSystem->Exit(1);
    }
  }

  char fname[4096];
  sprintf(fname,"%s/fit_parameters_ALL.txt", netdir);
  FILE* in_all = fopen(fname,"w");
  sprintf(fname,"%s/evt_parameters_ALL.txt", netdir);
  FILE* evt_all = fopen(fname,"w");
 
  for (int iset=0;iset<nset;iset++) {

    TChain sim("waveburst");
    TChain mdc("mdc");
    sim.Add(sim_file_name);
    mdc.Add(mdc_file_name);
    netevent  ww(&sim,nIFO);
    injection mm(&mdc,nIFO);

    // check vetoes
    UChar_t VETO[100];
    int nXDQF=0;
    for(int i=0;i<nVDQF;i++) {
      if((VDQF[i].cat!=CWB_CAT2)&&(VDQF[i].cat!=CWB_HVETO)&&(VDQF[i].cat!=CWB_CAT3)) continue;
      char veto_name[256];sprintf(veto_name,"%s_%s",CWB_CAT_NAME[VDQF[i].cat].Data(),VDQF[i].ifo);
      if(!CWB::Toolbox::isLeafInTree(ww.fChain,veto_name)) {
         cout << "cwb_report_sim.C Error : Leaf " << veto_name
              << " is not present in tree" << endl;
         exit(1);
      }
    }

    // build XDQF structure used for internal's computations
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

    // set veto branches 
    for(int j=0;j<nXDQF;j++) {
      char veto_name[256];	
      sprintf(veto_name,"%s_%s",CWB_CAT_NAME[XDQF[j].cat].Data(),XDQF[j].ifo);
      ww.fChain->SetBranchAddress(veto_name,&VETO[j]);
    }

    // sigmoid fit function
    double inf = simulation==2 ? log10(factors[0]) : -25;
    double sup = simulation==2 ? log10(factors[nfactor-1]) : -19;

    if(simulation==1 && pp_factor2distance) {
      inf = log10(pp_factor2distance/factors[nfactor-1]);
      sup = log10(pp_factor2distance/factors[0]);
    }
   
    TF1 *gfit=new TF1("logNfit",logNfit,pow(10.0,inf),pow(10.0,sup),5);
    gfit->SetParameters((inf+sup)/2.,0.3,1.,1.);
    gfit->SetParNames("hrss50","sigma","betam","betap");
    gfit->FixParameter(4,pp_factor2distance);
    gfit->SetParLimits(0,inf,sup);
    if(simulation==1 && pp_factor2distance) {
//      gfit->SetParLimits(1,0.05,10.);
//      gfit->SetParLimits(2,0.05,3.);
    } else if(simulation==2) {
      gfit->SetParLimits(1,0.05,10.);
      gfit->SetParLimits(2,0.05,3.);
    } else {
      gfit->SetParLimits(1,0.08,10.);
      gfit->SetParLimits(2,0.00,3.);
    }
    gfit->SetParLimits(3,0.0,2.50);

    // legend pannel
    TLegend *legend = new TLegend(0.53,0.17,0.99,0.7,"","brNDC");
    legend->SetLineColor(1);
    legend->SetTextSize(0.033);
    legend->SetBorderSize(1);
    legend->SetLineStyle(1);
    legend->SetLineWidth(1);
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);

    //double xcor,esnr,null,etot,ecor,asnr,acor,a2or, nill, rho;
    double xcor,esnr,null,etot,asnr,acor,a2or, nill, rho;
    double a,b,vED,enrg[5];
    int j,n,m,k;
    size_t indx;

    wavearray<double>  sMDC[NMDC_MAX];
    wavearray<double>   eff[NMDC_MAX];
    wavearray<double>   err[NMDC_MAX];
    wavearray<double>   e00[NMDC_MAX];
    wavearray<double>  nMDC[NMDC_MAX];
    wavearray<double>  nTRG[NMDC_MAX];

    bool save;
  
    int ntrg = sim.GetEntries();
    int nmdc = mdc.GetEntries();

    cout<<"  total injections: " << nmdc <<endl;
  
    for(int i=0; i<ninj+1; i++) {
      if((i<ninj)&&(imdc_iset[i]!=iset)) continue;   // skip unwanted injection types

      sprintf(fname,"hT%d",i);
      hT[i] = new TH1F(fname,fname,250,-T_win,T_win);
      hT[i]->SetTitleOffset(1.3,"Y");
      hT[i]->GetXaxis()->SetTitle("detected-injected time, sec");
      hT[i]->GetYaxis()->SetTitle("events");
      sprintf(fname,"%s",imdc_name[i]); 
      hT[i]->SetTitle(fname); 
   
      sprintf(fname,"hF%d",i);
      hF[i] = new TH1F(fname,fname,512,imdc_flow[i],imdc_fhigh[i]);
      hF[i]->SetTitleOffset(1.3,"Y");
      if(i<ninj) hF[i]->GetXaxis()->SetTitle("frequency, Hz");
      else       hF[i]->GetXaxis()->SetTitle("detected-injected, Hz");
      hF[i]->GetYaxis()->SetTitle("events");
      sprintf(fname,"%s",imdc_name[i]);
      hF[i]->SetTitle(fname);

      for(int j=0; j<nIFO; j++) {
        sprintf(fname,"%s_hrss_%d",IFO[j],i);
        HH[j][i] = new TH2F(fname,"",500,-23.,-19.,500,-23.,-19.);
        HH[j][i]->SetTitleOffset(1.7,"Y");
        HH[j][i]->SetTitleOffset(1.7,"Y");
        HH[j][i]->GetXaxis()->SetTitle("injected log10(hrss)");
        HH[j][i]->GetYaxis()->SetTitle("detected log10(hrss)");
        HH[j][i]->SetStats(kFALSE);
        sprintf(fname,"%s: %s",IFO[j],imdc_name[i]);
        HH[j][i]->SetTitle(fname);
        HH[j][i]->SetMarkerColor(2);
      }
      for(int j=0; j<nIFO; j++) {
        sprintf(fname,"%s_nre_%d",IFO[j],i);
        NRE[j][i] = new TH2F(fname,"",500,0.,7.,500,-3,1.);
        NRE[j][i]->SetTitleOffset(1.7,"Y");
        NRE[j][i]->SetTitleOffset(1.7,"Y");
        NRE[j][i]->GetXaxis()->SetTitle("injected log10(snr^{2})");
        NRE[j][i]->GetYaxis()->SetTitle("log10(normalized residual energy)");
        NRE[j][i]->SetStats(kFALSE);
        sprintf(fname,"%s: %s",IFO[j],imdc_name[i]);
        NRE[j][i]->SetTitle(fname);
        NRE[j][i]->SetMarkerColor(2);
      }
    }

    TH1F* hcor = new TH1F("hcor","",100,0,1.);
    hcor->SetTitleOffset(1.3,"Y");
    hcor->GetXaxis()->SetTitle("network correlation (ecor)");
    hcor->GetYaxis()->SetTitle("events");
    hcor->SetStats(kFALSE);
  
    TH1F* hcc = new TH1F("hcc","",100,0,1.);
    hcc->SetTitleOffset(1.3,"Y");
    hcc->GetXaxis()->SetTitle("network correlation");
    hcc->GetYaxis()->SetTitle("events");
    hcc->SetStats(kFALSE);
  
    TH1F* hrho = new TH1F("hrho","",500,0,10.);
    hrho->SetTitleOffset(1.3,"Y");
    hrho->GetXaxis()->SetTitle("log10(#rho^2)");
    hrho->GetYaxis()->SetTitle("events");

    TH1F* eSNR = new TH1F("eSNR","",500,0,10.);
    eSNR->SetTitleOffset(1.3,"Y");
    eSNR->GetXaxis()->SetTitle("log10(SNR/det)");
    eSNR->GetYaxis()->SetTitle("events");
    eSNR->SetStats(kFALSE);
  
    TH2F* rhocc = new TH2F("rhocc","",100,0.,1.,500,0,10.);
    rhocc->SetTitleOffset(1.3,"Y");
    rhocc->GetXaxis()->SetTitle("network correlation");
    rhocc->GetYaxis()->SetTitle("log10(rho^2)");
    rhocc->SetStats(kFALSE);
  
    TH2F* skyc = new TH2F("skyc","",180,-180,180.,90,-90.,90.);
    skyc->SetTitleOffset(1.3,"Y");
    skyc->GetXaxis()->SetTitle("#Delta#phi, deg.");
    skyc->GetYaxis()->SetTitle("#Delta#theta, deg.");
    skyc->SetStats(kFALSE);
  
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
    pf_vED->GetXaxis()->SetTitle("H1H2 consistency");
    pf_vED->GetYaxis()->SetTitle("penalty");
    pf_vED->SetMarkerStyle(20);
    pf_vED->SetMarkerColor(1);
    pf_vED->SetMarkerSize(0.8);
    pf_vED->SetStats(kFALSE);

    TH2F* hrs1[nIFO];
    TH1F* ifoT[nIFO];
    TH2F* hrs2[nIFO];
 
    for(int i=0; i<nIFO; i++) {

      sprintf(fname,"ifoT%s",IFO[i]);
      ifoT[i] = new TH1F(fname,"",250,-0.05,0.05);
      ifoT[i]->SetTitleOffset(1.3,"Y");
      ifoT[i]->GetXaxis()->SetTitle("detected-injected time, sec");
      ifoT[i]->GetYaxis()->SetTitle("events");
      ifoT[i]->SetTitle(IFO[i]);
    
      sprintf(fname,"hrs1%s",IFO[i]);
      hrs1[i] = new TH2F(fname,"",500,-23.,-19.,500,-23.,-19.);
      hrs1[i]->SetTitleOffset(1.3,"Y");
      hrs1[i]->GetXaxis()->SetTitle("detected");
      hrs1[i]->GetYaxis()->SetTitle("injected");
      hrs1[i]->SetStats(kFALSE);
      hrs1[i]->SetTitle(IFO[i]);
      hrs1[i]->SetMarkerColor(2);
    
      sprintf(fname,"hrs2%s",IFO[i]);
      hrs2[i] = new TH2F(fname,"",500,-23.,-19.,500,-23.,-19.);
      hrs2[i]->SetTitleOffset(1.3,"Y");
      hrs2[i]->SetMarkerColor(2);
      hrs2[i]->SetStats(kFALSE);
    }

    // init simulations=4 stuff
    int nhrss=10;
    vector<double> vhrss[NMDC_MAX]; for(int i=0;i<NMDC_MAX;i++) vhrss[i].resize(nhrss+1);
    wavearray<double>* mhrss[NMDC_MAX]; for(int i=0;i<NMDC_MAX;i++) mhrss[i] = new wavearray<double>[nhrss];
    if(simulation==4) {

      // loop over injection tree : compute min/max injected hrss
      wavearray<double> hrss_min(NMDC_MAX); for(int i=0;i<NMDC_MAX;i++) hrss_min[i]=1e20;
      wavearray<double> hrss_max(NMDC_MAX); for(int i=0;i<NMDC_MAX;i++) hrss_max[i]=0.;
      for(int i=0; i<nmdc; i++) {
        mm.GetEntry(i);
        if(pp_merge_types) {
          indx=0;
        } else {
          if(imdc_tset[mm.type-1] != iset) continue;   		// skip unwanted injections
          indx = imdc_index[mm.type-1];
        }
        if(mm.strain<hrss_min[indx]) if(mm.strain>0) hrss_min[indx]=mm.strain;
        if(mm.strain>hrss_max[indx]) hrss_max[indx]=mm.strain;
      }
      for(int i=0;i<NMDC_MAX;i++) {
        sMDC[i].resize(nhrss); 
        sMDC[i]=0.;
        if(hrss_min[i]==1e20) continue;
        if(hrss_max[i]==0) continue;
        vhrss[i][0]=hrss_min[i];
        vhrss[i][nhrss]=hrss_max[i];
        double fhrss=exp(log(vhrss[i][nhrss]/vhrss[i][0])/nhrss);
        for(int j=1;j<=nhrss;j++) vhrss[i][j]=fhrss*vhrss[i][j-1];
        for(int j=0;j<nhrss;j++) sMDC[i][j]=(vhrss[i][j+1]+vhrss[i][j])/2.; 
      }
      for(int j=0;j<nhrss;j++) for(int k=0; k<NMDC_MAX; k++) {nMDC[k].append(1.);}
    }

    // loop over injection tree
    for(int i=0; i<nmdc; i++) {

      mm.GetEntry(i);
      if(pp_merge_types) {
        indx=0;	
      } else {
        if(imdc_tset[mm.type-1] != iset) continue;   		// skip unwanted injections
        indx = imdc_index[mm.type-1];
      }  
      save = true;

      n = sMDC[0].size();
      for(int j=0; j<n; j++) {
        switch(simulation) {
        case 1:							// factors are hrss multipliers
          if(pp_factor2distance) {
            if(fabs(sMDC[0].data[j]-mm.factor)<EPSILON*mm.factor) { 
  	      save = false; nMDC[indx].data[j]+=1.; j=n; 
            }
          } else {
            if(fabs(sMDC[0].data[j]-mm.strain)<0.1*mm.strain) { 	// WARNING DOUBLE HRSS
  	      save = false; nMDC[indx].data[j]+=1.; j=n; 
            }
          }
          break;
        case 2:							// factors are snr
          if(fabs(sMDC[0].data[j]-mm.factor)<EPSILON*mm.factor) { 
  	    save = false; nMDC[indx].data[j]+=1.; j=n; 
          }
          break;
        case 3:							// factors are time shift
        case 4:							// nfactor is number of retries
          if(mm.strain>vhrss[indx][j] && mm.strain<=vhrss[indx][j+1]) { 
  	    nMDC[indx].data[j]+=1.; 
            mhrss[indx][j].append(mm.strain); 
            j=n; 
          }
          break;
        } 
      }
      // simulation=4 is already performed within the previous for loop
      if(save && (simulation!=4)) { 
        switch(simulation) {
        case 1:
          if(pp_factor2distance) {
            for(int k=0; k<NMDC_MAX; k++) sMDC[k].append(mm.factor); 
          } else {
            for(int k=0; k<NMDC_MAX; k++) sMDC[k].append(mm.strain); 
          } 
          break;
        case 2:
          for(int k=0; k<NMDC_MAX; k++) sMDC[k].append(mm.factor); 
          break;
        case 3:
        case 4:
          break;
        } 
        for(int k=0; k<NMDC_MAX; k++) {nMDC[k].append(0.);}
        nMDC[indx].data[nMDC[k].size()-1]+=1.;
      }
    }
    if(simulation==1 && sMDC[0].size()>nfactor) {
       cout << "cwb_report_sim.C Error : number of detected strains " << sMDC[0].size()
            << " is greater of nfactor declared in the user_parameters.C : " << nfactor << endl;
       bool answer = CWB::Toolbox::question("do you want to ignore this warning ? ");
       if(!answer) exit(1);
    }

    // for simulation==4 -> compute sMDC 
    if(simulation==4) {
      for(int k=0; k<NMDC_MAX; k++) {
        if(sMDC[k].size()==0) {delete [] mhrss[k];continue;}
        // compute sMDC using median50
        for(int i=0;i<nhrss;i++) {
          int K = mhrss[k][i].size();
          int* hrss_index = new int[K];
          if(K>0) TMath::Sort(K,mhrss[k][i].data,hrss_index,kFALSE);
          //if(K>0) cout << "XDEB " << i << " " << mhrss[k][i].data[K/2] << endl;
          if(K>0) sMDC[k][i] = mhrss[k][i].data[hrss_index[K/2]];
          delete [] hrss_index;
          mhrss[k][i].resize(0);
        }
/*
        // compute sMDC using weighing hrss with inj distribution
        for(int i=0;i<nhrss;i++) {
          int K = mhrss[k][i].size();
          double A=0; 
          double B=0; 
          for(int j=0;j<K;j++) {
            double X=mhrss[k][i].data[j];
            A+=X;
            B+=1;
          }
          if(K>0) sMDC[k][i] = A/B;
          mhrss[k][i].resize(0);
        }
*/
        delete [] mhrss[k];
      }
    }

    int* mdc_index[NMDC_MAX];
    for(int k=0; k<NMDC_MAX; k++) {
      if(sMDC[k].size()==0) continue;
      mdc_index[k] = new int[sMDC[k].size()];
      TMath::Sort((int)sMDC[k].size(),sMDC[k].data,mdc_index[k],kFALSE);
    }

    n = sMDC[0].size();
    for(int i=0; i<NMDC_MAX; i++){
      nTRG[i].resize(n); nTRG[i]=0.;
      err[i].resize(n);  err[i]=0.;
      e00[i].resize(n);  e00[i]=0.;
      eff[i].resize(n);  eff[i]=0.;
    }
  
    cout<<"total events: " << ntrg << endl;

    // check if ifar exist in the wave tree
    TBranch* branch;
    bool ifar_exists=false;
    TIter next(ww.fChain->GetListOfBranches());
    while ((branch=(TBranch*)next())) {
      if(TString("ifar").CompareTo(branch->GetName())==0) ifar_exists=true;
    }
    float ifar=0;
    if (ifar_exists) ww.fChain->SetBranchAddress("ifar",&ifar);

    // loop over trigger tree
    for(int i=0; i<ntrg; i++) {

      ww.GetEntry(i);
      if(ww.type[1]<1) continue;
      if(pp_merge_types) {
      } else {
        if(imdc_tset[ww.type[1]-1] != iset) continue;   // skip unwanted injections
      }

      // check veto cuts
      bool bveto=false;
      for(int j=0;j<nXDQF;j++) {
        if((XDQF[j].cat==CWB_CAT2)&&(!VETO[j])) bveto=true;
        if((XDQF[j].cat==CWB_CAT3)&&(VETO[j]))  bveto=true;
        if((XDQF[j].cat==CWB_HVETO)&&(VETO[j])) bveto=true;
      }
      if(bveto) continue;

      vED = ww.neted[0]/ww.ecor;

      if (T_vED>0) { if(vED > T_vED) continue;}

      if (T_ifar>0) if(ifar_exists) {if(ifar < T_ifar) continue;}

      // frequencies user bud cut
      bool bcut=false;
      for(int j=0;j<nFCUT;j++) {
        if((ww.frequency[0]>lowFCUT[j])&&(ww.frequency[0]<=highFCUT[j])) bcut=true;
      }
      if(bcut) continue;

      if(ww.ecor<T_acor) continue;
      //if(ww.ecor<0) continue;
    
      if(fabs(ww.time[0]-ww.time[nIFO])>T_win) continue;    // skip random coincidence
      if(pp_merge_types) {
        indx=0;	
      } else {
        indx = imdc_index[ww.type[1]-1];   
      }

      null = nill = etot = 0.;
      for(int j=0; j<nIFO; j++) {
        null += ww.null[j];
        nill += ww.nill[j];
        enrg[j]= ww.snr[j]-ww.null[j];
        etot += enrg[j];
      }

      a2or = 1.e50;
      for(int j=0; j<nIFO; j++) {
        if(a2or > etot-enrg[j]) a2or = etot-enrg[j]; 
      }
      //    if(a2or<0.) continue;

      acor = sqrt(ww.ECOR/nIFO);
      //rho = sqrt(ww.ECOR*fabs(ww.netcc[0])/nIFO);	// = rho[1]
      rho = ww.rho[pp_irho];

      xcor = ww.netcc[pp_inetcc];

      hcor->Fill(ww.netcc[pp_inetcc]);
      hcc->Fill(xcor);
      rhocc->Fill(xcor,log10(rho*rho));
      pf_cc->Fill(xcor,ww.penalty);
    
      if(xcor<T_cor) continue;

      pf_vED->Fill(vED,ww.penalty);
      hvED->Fill(fabs(vED));
      hpf->Fill(ww.penalty);

      if(ww.netcc[3] < T_scc) continue;
      if(ww.penalty < T_pen) continue;
      if (T_hrss>0) if (TMath::Abs((TMath::Log10(ww.hrss[i_hrss1])-TMath::Log10(ww.hrss[i_hrss2])))>T_hrss) continue;

      bool save=false;
      if(rho > T_cut) save=true;
      if (!save) continue;

      hrho->Fill(log10(rho*rho));
      eSNR->Fill(log10(ww.likelihood/nIFO));
        
      a = ww.phi[0]-ww.phi[1];
      if(a>  180) a =  360-a;
      if(a< -180) a = -360-a;
      b = ww.theta[0]-ww.theta[1];
      if(b>  90) b =  180-b;
      if(b< -90) b = -180-b;
      skyc->Fill(a,b); 

      for(int j=0; j<nIFO; j++) {
        HH[j][ninj]->Fill(log10(ww.hrss[j+nIFO]),log10(ratio_hrss*ww.hrss[j]));
        HH[j][indx]->Fill(log10(ww.hrss[j+nIFO]),log10(ratio_hrss*ww.hrss[j]));

        double nre = ww.iSNR[j] ? (ww.iSNR[j]+ww.oSNR[j]-2*ww.ioSNR[j])/ww.iSNR[j] : 10;
        if(nre<=0) nre=10; 
        NRE[j][ninj]->Fill(log10(ww.iSNR[j]),log10(nre));  
        NRE[j][indx]->Fill(log10(ww.iSNR[j]),log10(nre));  

        ifoT[j]->Fill(ww.time[j]-ww.time[j+nIFO]);
      }
      
      hT[indx]->Fill(ww.time[0]-ww.time[nIFO]);
      hF[indx]->Fill(float(ww.frequency[0]));
      hT[ninj]->Fill(ww.time[0]-ww.time[nIFO]);
      hF[ninj]->Fill(float(ww.frequency[0]-(imdc_flow[indx]+imdc_fhigh[indx])/2));
    
      n = sMDC[0].size();

      for(int j=0; j<n; j++) {
        switch(simulation) {
        case 1:
          if(pp_factor2distance) {
            if(fabs(sMDC[0].data[j]-ww.factor)<EPSILON*ww.factor) { 
  	      nTRG[indx].data[j]+=1.; j=n; 
            }
          } else {
            if(fabs(sMDC[0].data[j]-ww.strain[1])<0.1*ww.strain[1]) { //WARNING DOUBLE HRSS
  	      nTRG[indx].data[j]+=1.; j=n; 
            }
          }
          break;
        case 2:
          if(fabs(sMDC[0].data[j]-ww.factor)<EPSILON*ww.factor) { 
  	    nTRG[indx].data[j]+=1.; j=n; 
          }
          break;
        case 3:
        case 4:
          if(ww.strain[1]>vhrss[indx][j] && ww.strain[1]<=vhrss[indx][j+1]) { 
  	    nTRG[indx].data[j]+=1.; j=n; 
          }
          break;
        } 
      }
    }

    if(simulation==1 && pp_factor2distance) {  
      for(int k=0; k<NMDC_MAX; k++) for(int i=0;i<sMDC[k].size();i++) sMDC[k].data[i] = pp_factor2distance/sMDC[k].data[i]; 
    }

    n = sMDC[0].size();
    cout << "INJ HRSS NUMBER: " << n << endl;

    FILE* in;
    char f_file[256];
    for(int i=0; i<ninj; i++) {
      if(imdc_iset[i]!=iset) continue;   // skip unwanted injection types
      sprintf(f_file,"%s/eff_%s.txt",netdir, imdc_name[i]);
      in = fopen(f_file,"w");
      if(in==NULL) {cout << "cwb_report_sim : Error opening file " << f_file << endl;exit(1);}
      for(int k=0; k<n; k++){
        int j=mdc_index[i][k];
        if(nMDC[i].data[j] <= 0.) continue;
        eff[i].data[j] = nTRG[i].data[j]/nMDC[i].data[j];
        err[i].data[j] = eff[i].data[j]>1 ? 0.01 : 
        sqrt(eff[i].data[j]*(1-eff[i].data[j])/nMDC[i].data[j]);
        fprintf(in,"%e %i %i %.3f\n", sMDC[i].data[j], (int)nTRG[i].data[j], (int)nMDC[i].data[j], eff[i].data[j]);
      }
      fclose(in);
    }
 
    TCanvas *c1 = new TCanvas("c","C",0,0,800,800);
    c1->SetBorderMode(0);
    c1->SetFillColor(0);
    c1->SetBorderSize(2);
    c1->SetLogx(kFALSE);
    c1->SetGridx();
    c1->SetGridy();
    c1->SetLeftMargin(0.137039);
    c1->SetRightMargin(0.117039);
    c1->SetTopMargin(0.0772727);
    c1->SetBottomMargin(0.143939);
    c1->ToggleEventStatus();

    for(int i=0; i<ninj; i++) {
      if(imdc_iset[i]!=iset) continue;   // skip unwanted injection types

      fprintf(evt_all,"%2d %s ", i,imdc_name[i]);

      c1->SetLogx(kFALSE);
      c1->SetLogy();
   
      c1->Clear();
      hF[i]->Draw("HIST");
      sprintf(fname,"%s/f_%s.gif",netdir,imdc_name[i]);
      fprintf(evt_all,"%3.2f %3.2f ", hF[i]->GetMean(), hF[i]->GetRMS());
      cout<<fname<<endl;
      c1->Update(); c1->SaveAs(fname);
   
      c1->Clear();
      hT[i]->Draw("HIST");
      sprintf(fname,"%s/t_%s.gif",netdir,imdc_name[i]);
      fprintf(evt_all,"%3.6f %3.6f ", hT[i]->GetMean(), hT[i]->GetRMS());
      c1->Update(); c1->SaveAs(fname);
   
      c1->SetLogx(kFALSE);
      c1->SetLogy(kFALSE);
      for(int j=0;j<nIFO;j++){
        c1->Clear();
        HH[j][i]->Draw("HIST");
        sprintf(fname,"%s/hrss_%s_%s.gif",netdir,IFO[j],imdc_name[i]);
        c1->Update(); c1->SaveAs(fname); 

        c1->Clear();
        NRE[j][i]->Draw();
        sprintf(fname,"%s/nre_%s_%s.gif",netdir,IFO[j],imdc_name[i]);
        c1->Update(); c1->SaveAs(fname); 
      }

      fprintf(evt_all,"\n");
    }

    c1->Clear();
    c1->SetLogy(kTRUE);
    c1->SetLogx(kFALSE);
    hcor->Draw("HIST");
    sprintf(fname,"%s/netcor_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);
    c1->SetLogy(kFALSE);

    c1->Clear();
    c1->SetLogy(kTRUE);
    c1->SetLogx(kFALSE);
    hcc->Draw("HIST");
    sprintf(fname,"%s/cc_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);
    c1->SetLogy(kFALSE);

    c1->Clear();
    c1->SetLogy(kTRUE);
    hrho->Draw("HIST");
    sprintf(fname,"%s/rho_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);

    c1->Clear();
    hpf->Draw("HIST");
    sprintf(fname,"%s/penalty_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);

    c1->Clear();
    eSNR->Draw("HIST");
    sprintf(fname,"%s/average_snr_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);
    c1->SetLogy(kFALSE);

    c1->Clear();
    c1->SetLogz(kTRUE);
    rhocc->Draw("colz");
    sprintf(fname,"%s/rho_cc_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);

    c1->Clear();
    c1->SetLogz(kTRUE);
    pf_cc->Draw("colz");
    sprintf(fname,"%s/penalty_cc_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);

    c1->Clear();
    c1->SetLogz(kTRUE);
    skyc->Draw("colz");
    sprintf(fname,"%s/coordinate_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);
    c1->SetLogz(kFALSE);

    c1->Clear();
    c1->SetLogz(kTRUE);
    pf_vED->Draw("colz");
    sprintf(fname,"%s/penalty_vED_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);
    c1->SetLogz(kFALSE);

    c1->Clear();
    c1->SetLogy(kTRUE);
    hvED->Draw("HIST");
    sprintf(fname,"%s/vED_%s.gif",netdir,imdc_set_name[iset].Data());
    c1->Update(); c1->SaveAs(fname);
    c1->SetLogy(kFALSE);
 
    // Efficiency
    if(c1) delete c1;
    c1 = new TCanvas("c","C",0,0,1000,800);
    c1->SetBorderMode(0);
    c1->SetFillColor(0);
    c1->SetBorderSize(2);
    c1->SetLogx(kFALSE);
    c1->SetGridx();
    c1->SetGridy();
    c1->SetRightMargin(0.0517039);
    c1->SetTopMargin(0.0772727);
    c1->SetBottomMargin(0.143939);
    c1->ToggleEventStatus();

    double chi2, hrss50, sigma, betam, betap, hrssEr;

    n = sMDC[0].size();

    sprintf(fname,"%s/fit_parameters_%s.txt", netdir,imdc_set_name[iset].Data());
    in = fopen(fname,"w");

    int iset_size=0;	// iset size
    for(int i=0; i<ninj; i++) if(imdc_iset[i]==iset) iset_size++;  
    double hleg = 0.18+iset_size*0.06;
    delete legend;
    if(pp_show_eff_fit_curve) legend = new TLegend(0.45,0.18,0.99,hleg,NULL,"brNDC");
    else                      legend = new TLegend(0.55,0.18,0.99,hleg,NULL,"brNDC");

    for(int i=0; i<ninj; i++) {
      if(imdc_iset[i]!=iset) continue;   // skip unwanted injection types

      c1->Clear();
      c1->SetLogx();
      EFF[i]=new TGraphErrors(n,sMDC[i].data,eff[i].data,e00[i].data,err[i].data);
      if(simulation==1 && pp_factor2distance) EFF[i]->GetXaxis()->SetTitle("distance (Kpc)");
      else if(simulation==2)                  EFF[i]->GetXaxis()->SetTitle("snr");
      else                                    EFF[i]->GetXaxis()->SetTitle("hrss, #frac{strain}{#sqrt{Hz}}");
      EFF[i]->GetYaxis()->SetTitle("Efficiency    ");
      EFF[i]->GetXaxis()->SetTitleOffset(1.7);
      EFF[i]->GetXaxis()->SetLabelSize(0.04);
      EFF[i]->GetYaxis()->SetTitleOffset(1.3);
      EFF[i]->GetYaxis()->SetRangeUser(0,1);
      EFF[i]->GetYaxis()->SetLabelSize(0.04);

      //gfit->SetParameters(-21.,1.,1.,1.);
      //gfit->SetParNames("hrss50","sigma","betam","betap");
      gfit->SetLineColor(colors[i]);
      //gfit->SetParLimits(0,-22.,-19.);
      //gfit->SetParLimits(1,0.01,50.);
      //gfit->SetParLimits(2,0.01,50.);
      //gfit->SetParLimits(3,0.01,50.);
      gfit->SetNpx(100000);
      if(pp_show_eff_fit_curve) EFF[i]->Fit("logNfit"); else EFF[i]->Fit("logNfit","N");  
      chi2=gfit->GetChisquare();
      hrss50=pow(10.,gfit->GetParameter(0));
      hrssEr=(pow(10.,gfit->GetParError(0))-1.)*hrss50;
      sigma=gfit->GetParameter(1);
      betam=gfit->GetParameter(2);
      betap=gfit->GetParameter(3);

      printf("%d %s %E %E +- %E %E %E %E\n",
    	     i,imdc_name[i],chi2,hrss50,hrssEr,sigma,betam,betap);
      fprintf(in,"%2d %9.3E %9.3E +- %9.3E %9.3E %9.3E %9.3E %s\n",
              i,chi2,hrss50,hrssEr,sigma,betam,betap,imdc_name[i]);
      fprintf(in_all,"%2d %9.3E %9.3E +- %9.3E %9.3E %9.3E %9.3E %s\n",
              i,chi2,hrss50,hrssEr,sigma,betam,betap,imdc_name[i]);
  
      sprintf(fname,"%s, hrss50=%5.2E", imdc_name[i],hrss50);
      EFF[i]->SetTitle(fname);
      //gStyle->SetOptStat(01);
      gfit->SetLineColor(colors[i]);
      EFF[i]->SetLineColor(colors[i]);
      EFF[i]->SetMarkerStyle(markers[i]);
      EFF[i]->SetMarkerColor(colors[i]);
      EFF[i]->SetMarkerSize(2);

      legend->SetBorderSize(1);
      legend->SetTextAlign(22);
      legend->SetTextFont(12);
      legend->SetLineColor(1);
      legend->SetLineStyle(1);
      legend->SetLineWidth(1);
      legend->SetFillColor(0);
      legend->SetFillStyle(1001);
      legend->SetTextSize(0.04);
      legend->SetLineColor(kBlack);
      legend->SetFillColor(kWhite);

      if(pp_show_eff_fit_curve) sprintf(fname,"%s, %5.2E",imdc_name[i], hrss50);
      else                      sprintf(fname,"%s",imdc_name[i]);
      legend->AddEntry(EFF[i], fname, "lp");
      if(pp_show_eff_fit_curve) EFF[i]->Draw("AP"); else EFF[i]->Draw("APL");  
      sprintf(fname,"%s/%s.gif",netdir,imdc_name[i]);
      c1->SaveAs(fname);
      c1->SetLogx(kFALSE);
      c1->Clear();
    }
    fclose(in);

    char gname[512];
    sprintf(fname,"%s/eff_%s.eps", netdir,imdc_set_name[iset].Data());
    sprintf(gname,"%s/eff_%s.gif", netdir,imdc_set_name[iset].Data());
    TPostScript epsfile(fname,113); 

    c1->Clear();
    if(simulation==2) c1->SetLogx(kFALSE); else c1->SetLogx();
    c1->SetLogy(kFALSE);

    char ptitle[256];
    bool first=true;
    TMultiGraph* mg=new TMultiGraph();
    sprintf(ptitle,"%s   rho=%3.2f, cc=%3.2f",imdc_set_name[iset].Data(),T_cut,T_cor);
    mg->SetTitle(ptitle);
    double xmin=1.79769e+308;
    double xmax=0;
    for(int i=0;i<ninj;i++) {
      if(imdc_iset[i]!=iset) continue;   // skip unwanted injection types
      mg->Add(EFF[i]);
      if(sMDC[i].data[mdc_index[i][0]]<xmin) xmin=sMDC[i].data[mdc_index[i][0]];
      if(sMDC[i].data[mdc_index[i][sMDC[i].size()-1]]>xmax) xmax=sMDC[i].data[mdc_index[i][sMDC[i].size()-1]];
    }
    if(pp_show_eff_fit_curve) mg->Draw("AP"); else mg->Draw("APL"); 
    if(simulation==1 && pp_factor2distance) mg->GetHistogram()->GetXaxis()->SetTitle("distance (Kpc)");
    else if(simulation==2)                  mg->GetHistogram()->GetXaxis()->SetTitle("snr");
    else                                    mg->GetHistogram()->GetXaxis()->SetTitle("hrss, #frac{strain}{#sqrt{Hz}}");
    mg->GetHistogram()->GetYaxis()->SetTitle("Efficiency    ");
    mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.7);
    mg->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
    mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
    mg->GetHistogram()->GetXaxis()->SetRangeUser(xmin,xmax);
    mg->GetHistogram()->GetYaxis()->SetRangeUser(0,1);
    mg->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
    legend->Draw("SAME");
    c1->Update(); epsfile.Close();
    char cmd[1024];
    sprintf(cmd,"convert %s %s",fname,gname);cout<<cmd<<endl;
    gSystem->Exec(cmd);

    // delete histograms
    for (int i=0;i<ninj+1;i++) {
      if((i<ninj)&&(imdc_iset[i]!=iset)) continue;   // skip unwanted injection types
      for(int j=0;j<nIFO;j++) if(HH[j][i]) delete HH[j][i];
      for(int j=0;j<nIFO;j++) if(NRE[j][i]) delete NRE[j][i];
      if(hT[i]) delete hT[i];
      if(hF[i]) delete hF[i];
    }

    delete hcor;
    delete hcc;
    delete hrho;
    delete eSNR;
    delete rhocc;
    delete skyc;
    delete hpf;
    delete hvED;
    delete pf_cc;
    delete pf_vED;

    for (int j=0;j<nIFO;j++) {
      delete hrs1[j];
      delete hrs2[j];
      delete ifoT[j];
    }

    if(simulation==4) {
      for(int k=0; k<NMDC_MAX; k++) {
        if(sMDC[k].size()==0) continue;
        delete [] mdc_index[k];
      }
    }
    for(int i=0;i<NMDC_MAX;i++) if(EFF[i]) {delete EFF[i];EFF[i]=NULL;}
    delete mg;
  }

  fclose(in_all);
  fclose(evt_all);

  exit(0);
}
