{
  cout << "-----> macro/CWB_Plugin_TestClassMDC_Config.C" << endl;

  int seed = net->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc MDC(net);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  par.resize(2);
  par[0].name="frequency"; par[0].value=1053.;
  //par[1].name="Q"; par[1].value=3.;
  par[1].name="Q"; par[1].value=9.;
  //par[1].name="Q"; par[1].value=100.;
  MDC.AddWaveform(MDC_SG, par);
  //MDC.AddWaveform(MDC_SGC, par);

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

  vector<mdcpar> par; par.resize(4);
//  par[0].name="theta";  par[0].value=20;     // theta
//  par[1].name="phi";    par[1].value=150;        // phi

//  par[0].name="theta";  par[0].value=60;     // theta
//  par[1].name="phi";    par[1].value=200;        // phi

  par[0].name="theta";  par[0].value=90-125.96426;     // theta
  par[1].name="phi";    par[1].value=122.93704;        // phi
  par[2].name="psi";    par[2].value=38.769279;         // psi

//  par[0].name="theta";  par[0].value=90-119.912102;     // theta
//  par[1].name="phi";    par[1].value=249.969513;        // phi
//  par[2].name="psi";    par[2].value=38.783863;         // psi
  par[3].name="gps";    par[3].value=931158395;         // gps
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par);
}
