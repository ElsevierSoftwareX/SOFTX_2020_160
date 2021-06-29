{
  cout << "-----> macro/CWB_Plugin_TestClassMDC_Config.C" << endl;

  int seed = net->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc MDC(net);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=9.;
  MDC.AddWaveform(MDC_SGC, par);

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
  vector<mdcpar> par(4);
  par[0].name="theta"; par[0].value=40;
  par[1].name="phi"; par[1].value=70;
  par[2].name="psi"; par[2].value=0;
  par[3].name="rho"; par[3].value=1;
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par);
}
