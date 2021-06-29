{
  cout << "-----> macro/CWB_Plugin_EBBH_Config.C" << endl;

  int seed = net->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc MDC(net);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  par.resize(1);
  par[0].name="macro/eBBH_3.5e06_10_25_2.lst"; 
  MDC.AddWaveform(MDC_EBBH, par);
/*
  par.resize(2);
  par[0].name="macro/eBBH_3.5e06_10_25_2.root"; 
  par[1].name="rp0>26"; 
  MDC.AddWaveform(MDC_EBBH, par);
*/
  MDC.Print(1);

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC.SetInjHrss(2.5e-21);
  MDC.SetInjRate(0.005);
  MDC.SetInjJitter(10.0);
  MDC.SetInjLength(20.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  vector<mdcpar> par(6);
  par[0].name="theta"; par[0].value=40;
  par[1].name="phi";   par[1].value=70;
  par[2].name="psi";   par[2].value=0;
  par[3].name="rho";   par[3].value=1;
  par[4].name="gps";   par[4].value=931158300;
  par[5].name="iota";  par[5].value=0;          // the angle iota (deg) is the inclination of the system
                                                // which originates the burst with respect to the line of sight
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par);
}
