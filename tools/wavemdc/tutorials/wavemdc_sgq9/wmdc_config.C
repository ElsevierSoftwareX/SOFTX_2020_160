{
  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "TestWaveMDC";

  // define frame segment list
  TString segmentList = "Segments/S6a_LS-segs.txt";

  // select segments
  int jobmin  = 210;    // start segment
  int jobmax  = 312;    // end segment
  int jobstep =  10;    // frames per job 

  // -------------------------------------------------------- 
  // define network
  // -------------------------------------------------------- 
  int nIFO=3;
  TString ifo[nIFO]={"L1","H1","V1"};
  CWB::mdc MDC(nIFO,ifo); 

  // --------------------------------------------------------
  // define channel names
  // --------------------------------------------------------
  vector<TString> chName;

  // -------------------------------------------------------- 
  // define waveforms
  // -------------------------------------------------------- 

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par;

  par.resize(2);

  double F_Q8d9[7] = {70,100,153,235,361,554,849};
  for(int n=0;n<7;n++) {
    par[0].name="frequency"; par[0].value=F_Q8d9[n];
    par[1].name="Q"; par[1].value=8.9;
    MDC.AddWaveform(MDC_SGE, par);
  }

  MDC.Print();

  // define injection parameters
  // --------------------------------------------------------
  MDC.SetInjHrss(2.5e-21);
  MDC.SetInjRate(0.01);
  MDC.SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  vector<mdcpar> par(4);
  par[0].name="entries"; par[0].value=100000;
  par[1].name="rho_min"; par[1].value=0;
  par[2].name="rho_max"; par[2].value=0;
  par[3].name="iota";    par[3].value=-1.;       // random iota
  MDC.SetSkyDistribution(MDC_RANDOM,par);

}
