{

  // Config Plugin to generate simulated gaussian noise and injected 'on the fly' Burst MDC

  cout << "-----> plugins/CWB_Plugin_BRST2_Config.C" << endl;

  int seed = net->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc MDC(net);

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
    MDC.AddWaveform(MDC_WNB, par);
  }
/*
  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=3.;
  MDC.AddWaveform(MDC_SG, par);

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=8.9;
  MDC.AddWaveform(MDC_SG, par);

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=9.;
  MDC.AddWaveform(MDC_SGC, par);
*/
  MDC.Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC.SetInjHrss(2.5e-21);
  MDC.SetInjRate(0.0333333);
  MDC.SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
/*
  vector<mdcpar> par(5);
  par[0].name="theta"; par[0].value=40;
  par[1].name="phi";   par[1].value=70;
  par[2].name="psi";   par[2].value=0;
  par[3].name="rho";   par[3].value=0;
  par[4].name="gps";   par[4].value=931158395;
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par);
*/

  par.resize(0);
  MDC.SetSkyDistribution(MDC_CUSTOM,"macro/mdc.lst",par);
}
