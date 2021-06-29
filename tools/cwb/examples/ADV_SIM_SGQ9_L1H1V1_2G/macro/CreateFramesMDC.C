{
  #define  N_IFO	3

  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "SGQ9";

  // -------------------------------------------------------- 
  // define network
  // -------------------------------------------------------- 
  TString ifo[N_IFO]={"L1","H1","V1"};
  CWB::mdc MDC(N_IFO,ifo); 

  // -------------------------------------------------------- 
  // define waveforms
  // -------------------------------------------------------- 
  #define nWF	15
  double F[nWF] = {70,100,153,235,361,554,849,945,1053,1172,1304,1451,1615,1797,2000};
  double Q[nWF] = {8.9,8.9,8.9,8.9,8.9,8.9,8.9,9,9,9,9,9,9,9,9};

  char wf_name[256];
  waveform wf;
  wf.par.resize(2);
  for (int n=0;n<nWF;n++) {
    wf.type = MDC_SG;
    sprintf(wf_name,"SG%dQ%1.1f",(int)F[n],Q[n]);
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    wf.name.ReplaceAll("d0","");
    wf.hp = MDC.GetSGQ(F[n],Q[n]);
    wf.hpPath = wf.name;
    wf.hx = wf.hp;
    wf.hx = 0;
    wf.hxPath = "";
    wf.par[0].name="freq"; wf.par[0].value=F[n];
    wf.par[1].name="Q";    wf.par[1].value=Q[n];
    MDC.AddWaveform(wf);
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
/*
  vector<mdcpar> par2; par2.resize(3);
  par2[0].name="entries";par2[0].value=100000;  // pool of events
  par2[1].name="rho_min";par2[1].value=1;         // min rho // Mpc
  par2[2].name="rho_max";par2[2].value=100;       // max rho // Mpc
  MDC.SetSkyDistribution(MDC_RANDOM,par2);
*/
/*
  vector<mdcpar> par3; par3.resize(1);
  par3[0].name="entries";par3[0].value=100000;  // number of entries
  MDC.SetSkyDistribution(MDC_MNGD,par3);
*/
/*
  vector<mdcpar> par4; par4.resize(1);
  par4[0].name="distance_thr";par4[0].value=100;  // distance max = 100 Mpc
  MDC.SetSkyDistribution(MDC_GWGC,"../data/GWGCCatalog_Rev1d2.txt",par4);
*/

  vector<mdcpar> par5; par5.resize(3);
  par5[0].name="theta";par5[0].value=30;     // theta
  par5[1].name="phi";  par5[1].value=60;     // phi
  par5[2].name="psi";  par5[2].value=90;     // psi
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par5);

  bool   blog     = true;
  size_t gps      = 931158100;
  size_t length   = 600;
  MDC.WriteFrameFile(frDir, frLabel, gps, length, blog);

  exit(0);
}
