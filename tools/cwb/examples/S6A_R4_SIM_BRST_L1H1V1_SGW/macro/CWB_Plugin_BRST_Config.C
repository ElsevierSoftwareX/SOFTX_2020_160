{
  cout << "-----> macro/CWB_Plugin_BRST_Config.C" << endl;

  int seed = net->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc MDC(net);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  par.resize(1);

  double F_GA[4]   = {0.0001,0.001,0.0025,0.004};
  for(int n=0;n<4;n++) {
    par[0].name="duration"; par[0].value=F_GA[n];
    MDC.AddWaveform(MDC_GA, par);
  }

  par.resize(2);

  double F_Q3[4] = {70,235,849,1615};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q3[n];
    par[1].name="Q"; par[1].value=3.;
    MDC.AddWaveform(MDC_SG, par);
  }

  double F_Q8d9[7] = {70,100,153,235,361,554,849};
  for(int n=0;n<7;n++) {
    par[0].name="frequency"; par[0].value=F_Q8d9[n];
    par[1].name="Q"; par[1].value=8.9;
    MDC.AddWaveform(MDC_SG, par);
  }

  double F_Q9[4]   = {1053,1304,1615,2000};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q9[n];
    par[1].name="Q"; par[1].value=9.;
    MDC.AddWaveform(MDC_SG, par);
  }

  double F_Q100[4] = {70,235,849,1615};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q100[n];
    par[1].name="Q"; par[1].value=100.;
    MDC.AddWaveform(MDC_SG, par);
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
      MDC.AddWaveform(MDC_WNB, par);
    }
  }

  MDC.Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC.SetInjHrss(2.5e-21);
  MDC.SetInjRate(0.01);
  MDC.SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  vector<mdcpar> par(3);
  par[0].name="entries"; par[0].value=100000;
  par[1].name="rho_min"; par[1].value=0;
  par[2].name="rho_max"; par[2].value=0;
  MDC.SetSkyDistribution(MDC_RANDOM,par,seed);
}
