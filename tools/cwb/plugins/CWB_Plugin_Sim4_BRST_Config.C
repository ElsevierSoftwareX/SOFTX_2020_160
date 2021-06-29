{
  //!NOISE_MDC_SIMULATION
  // Config Plugin to generate simulated MDC injected 'on the fly' using simulation=4


  cout << "-----> CWB_Plugin_Sim4_BRST_Config.C -> " << " run : " << net->nRun << " gIFACTOR " << gIFACTOR << endl;

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int seed = 100*(*net)->nRun+gIFACTOR; 

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  par.resize(1);

  double F_GA[4]   = {0.0001,0.001,0.0025,0.004};
  for(int n=0;n<4;n++) {
    par[0].name="duration"; par[0].value=F_GA[n];
    MDC->AddWaveform(MDC_GA, par);
  }

  par.resize(2);

  double F_Q3[4] = {70,235,849,1615};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q3[n];
    par[1].name="Q"; par[1].value=3.;
    MDC->AddWaveform(MDC_SG, par);
  }

 double F_Q8d9[7] = {70,100,153,235,361,554,849};
  for(int n=0;n<7;n++) {
    par[0].name="frequency"; par[0].value=F_Q8d9[n];
    par[1].name="Q"; par[1].value=8.9;
    MDC->AddWaveform(MDC_SG, par);
  }

  double F_Q9[4]   = {1053,1304,1615,2000};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q9[n];
    par[1].name="Q"; par[1].value=9.;
    MDC->AddWaveform(MDC_SG, par);
  }

  double F_Q100[4] = {70,235,849,1615};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q100[n];
    par[1].name="Q"; par[1].value=100.;
    MDC->AddWaveform(MDC_SG, par);
  }
  par.resize(6);

  double F_WNB[5] = {100,250,1000,1000,1000};
  double B_WNB[5] = {100,100,10,1000,1000};
  double D_WNB[5] = {0.1,0.1,0.1,0.1,0.01};
  for(int n=0;n<5;n++) {
    par[0].name="frequency"; par[0].value=F_WNB[n];
    par[1].name="bandwidth"; par[1].value=B_WNB[n];
    par[2].name="duration";  par[2].value=D_WNB[n];
    for(int m=0;m<30;m++) {
      par[3].name="pseed";     par[3].value=seed+100000+n*100+m;
      par[4].name="xseed";     par[4].value=seed+100001+n*100+m;
      par[5].name="mode";      par[5].value=1;  // simmetric
      MDC->AddWaveform(MDC_WNB, par);
    }
  }

  MDC->Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjHrss(2.5e-21);
  MDC->SetInjRate(0.01);
  MDC->SetInjJitter(10.0);
  MDC->SetSourceListSeed(gIFACTOR);

  // --------------------------------------------------------
  // define the MDC hrss 
  // --------------------------------------------------------
  int K=MDC.wfList.size();
  std::vector<double> mdcHRSS(K);
  for(int i=0;i<K;i++) {
    if(TString(MDC->wfList[i].name)=="GA0d1") 			mdcHRSS[i]=8.980e-22;
    if(TString(MDC->wfList[i].name)=="GA1d0") 			mdcHRSS[i]=5.097e-22;
    if(TString(MDC->wfList[i].name)=="GA2d5") 			mdcHRSS[i]=9.252e-22;
    if(TString(MDC->wfList[i].name)=="GA4d0") 			mdcHRSS[i]=1.440e-21;
    if(TString(MDC->wfList[i].name)=="SG70Q3") 			mdcHRSS[i]=1.041e-21;
    if(TString(MDC->wfList[i].name)=="SG235Q3") 		mdcHRSS[i]=3.180e-22;
    if(TString(MDC->wfList[i].name)=="SG849Q3") 		mdcHRSS[i]=9.420e-22;
    if(TString(MDC->wfList[i].name)=="SG1615Q3") 		mdcHRSS[i]=1.747e-21;
    if(TString(MDC->wfList[i].name)=="SG70Q8d9") 		mdcHRSS[i]=1.219e-21;
    if(TString(MDC->wfList[i].name)=="SG100Q8d9") 		mdcHRSS[i]=5.879e-22;
    if(TString(MDC->wfList[i].name)=="SG153Q8d9") 		mdcHRSS[i]=3.170e-22;
    if(TString(MDC->wfList[i].name)=="SG235Q8d9") 		mdcHRSS[i]=3.148e-22;
    if(TString(MDC->wfList[i].name)=="SG361Q8d9") 		mdcHRSS[i]=5.268e-22;
    if(TString(MDC->wfList[i].name)=="SG554Q8d9") 		mdcHRSS[i]=5.772e-22;
    if(TString(MDC->wfList[i].name)=="SG849Q8d9") 		mdcHRSS[i]=9.610e-22;
    if(TString(MDC->wfList[i].name)=="SG1053Q9") 		mdcHRSS[i]=1.025e-21;
    if(TString(MDC->wfList[i].name)=="SG1304Q9") 		mdcHRSS[i]=1.230e-21;
    if(TString(MDC->wfList[i].name)=="SG1615Q9") 		mdcHRSS[i]=1.681e-21;
    if(TString(MDC->wfList[i].name)=="SG2000Q9") 		mdcHRSS[i]=2.930e-21;
    if(TString(MDC->wfList[i].name)=="SG70Q100") 		mdcHRSS[i]=1.650e-21;
    if(TString(MDC->wfList[i].name)=="SG235Q100") 		mdcHRSS[i]=3.557e-22;
    if(TString(MDC->wfList[i].name)=="SG849Q100") 		mdcHRSS[i]=8.918e-22;
    if(TString(MDC->wfList[i].name)=="SG1615Q100") 		mdcHRSS[i]=1.555e-21;
    if(TString(MDC->wfList[i].name)=="WNB100_100_0d100") 	mdcHRSS[i]=4.378e-22;
    if(TString(MDC->wfList[i].name)=="WNB250_100_0d100") 	mdcHRSS[i]=4.405e-22;
    if(TString(MDC->wfList[i].name)=="WNB1000_10_0d100") 	mdcHRSS[i]=9.035e-22;
    if(TString(MDC->wfList[i].name)=="WNB1000_1000_0d100") 	mdcHRSS[i]=3.373e-21;
    if(TString(MDC->wfList[i].name)=="WNB1000_1000_0d010") 	mdcHRSS[i]=1.841e-21;
  }


  // --------------------------------------------------------
  // define event sky distribution
  // --------------------------------------------------------
  par.resize(3);
  par[0].name="entries";par[0].value=100000;    // pool of events
  par[1].name="rho_min";par[1].value=.2;        // min rho  Kpc
  par[2].name="rho_max";par[2].value=40.;       // max rho  Kpc
  MDC->SetSkyDistribution(MDC_RANDOM,"x+50./x",par,seed);

}
