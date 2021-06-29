{
  //!NOISE_MDC_SIMULATION
  // Config Plugin to generate simulated gaussian noise and injected 'on the fly' Burst MDC

  cout << "-----> plugins/CWB_Plugin_TestClassMDC_Config.C" << endl;

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int seed = (*net)->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  par.resize(5);
  par[0].name="frequency"; par[0].value=250.;
  par[1].name="bandwidth"; par[1].value=100.;
  par[2].name="duration";  par[2].value=0.1;
  for(int n=0;n<30;n++) {
    par[3].name="pseed";     par[3].value=seed+100000+n;
    par[4].name="xseed";     par[4].value=seed+100001+n;
    MDC->AddWaveform(MDC_WNB, par);
  }

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=3.;
  MDC->AddWaveform(MDC_SG, par);

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=8.9;
  MDC->AddWaveform(MDC_SG, par);

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=9.;
  MDC->AddWaveform(MDC_SGC, par);

  MDC->Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjHrss(2.5e-21);
  MDC->SetInjRate(0.0333333);
  MDC->SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  par.resize(3);
  par[0].name="entries";par[0].value=100000;    // pool of events
  par[1].name="rho_min";par[1].value=1;         // min rho // Mpc
  par[2].name="rho_max";par[2].value=100;       // max rho // Mpc
  MDC->SetSkyDistribution(MDC_RANDOM,par,seed);
}
