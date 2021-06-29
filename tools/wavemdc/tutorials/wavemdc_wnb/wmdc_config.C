{
  #include <vector>

  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "TestWNB";

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
  int nWF=8;
  double F[nWF] = {100. ,100. ,100  ,1000. ,1000. ,1000. ,1000. ,1000.};
  double B[nWF] = {10.  ,100. ,100  ,100.  ,100.  ,1000. ,1000. ,1000.};
  double D[nWF] = {0.1  ,0.1  ,0.01 ,0.1   ,0.01  ,0.1   ,0.01  ,0.001};

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; par.resize(3);
  for (int n=0;n<nWF;n++) {
    par[0].name="frequency"; par[0].value=F[n];
    par[1].name="bandwidth"; par[1].value=B[n];
    par[2].name="duration";  par[2].value=D[n];
    for(int m=0;m<15;m++) MDC.AddWaveform(MDC_WNB, par);
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
  vector<mdcpar> par; par.resize(3); 
  par[0].name="entries";par[0].value=100000;  // pool of events
  par[1].name="rho_min";par[1].value=1;         // min rho // Mpc
  par[2].name="rho_max";par[2].value=100;       // max rho // Mpc
  MDC.SetSkyDistribution(MDC_RANDOM,par);

  //vector<mdcpar> par; par.resize(1); 
  //par[0].name="entries";par[0].value=100000;  // number of entries
  //MDC.SetSkyDistribution(MDC_MNGD,par);

  //vector<mdcpar> par; par.resize(1); 
  //par[0].name="distance_thr";par[0].value=100;  // distance max = 100 Mpc
  //MDC.SetSkyDistribution(MDC_GWGC,"../data/GWGCCatalog_Rev1d2.txt",par);

  //vector<mdcpar> par; par.resize(3); 
  //par[0].name="theta";par[0].value=30;     // theta
  //par[1].name="phi";  par[1].value=60;     // phi
  //par[2].name="psi";  par[2].value=90;     // psi
  //MDC.SetSkyDistribution(MDC_EARTH_FIX,par);
}
