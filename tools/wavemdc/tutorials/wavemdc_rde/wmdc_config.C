{
  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "RDE_TST1";

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
  int nWF=3;
  double F[nWF]   = {1590,2090,2590};
  double TAU[nWF] = {0.2,0.2,0.2};

  MDC.SetInjLength(2.0);

  char wf_name[256];
  waveform wf;
  for (int n=0;n<nWF;n++) {
    sprintf(wf_name,"RDE_%d_%1.1f",F[n],TAU[n]);
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    // add 20 different iota for each waveform
    for(int m=0;m<20;m++) {
      double cos_iota=gRandom->Uniform(-1,1);
      double iota=acos(cos_iota)*180./TMath::Pi();
      wf.hp = MDC.GetRD(F[n],TAU[n],iota,0);
      wf.hx = MDC.GetRD(F[n],TAU[n],iota,1);
      wf.hpPath = wf.name;
      wf.hxPath = wf.name;
      wf.par.resize(3);
      wf.par[0].name="freq";
      wf.par[0].value=F[n];
      wf.par[1].name="tau";
      wf.par[1].value=TAU[n];
      wf.par[2].name="iota";
      wf.par[2].value=iota;
      MDC.AddWaveform(wf);
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
