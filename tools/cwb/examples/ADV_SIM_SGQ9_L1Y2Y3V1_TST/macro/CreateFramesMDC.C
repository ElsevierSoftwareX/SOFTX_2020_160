{
  #include "wat.hh"

  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "SGQ9";

  // -------------------------------------------------------- 
  // define network
  // -------------------------------------------------------- 
  int nIFO=4;
  TString ifo[nIFO]={"L1","Y2","Y3","V1"};

  detectorParams detParms[4] = {
                          {"L1", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, (    -197.716)},

                          {"Y2", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, ( +45-125.999)},  // H2 LCI
                          {"Y3", 46.4551, -119.408, 0.0, 0, ( +45-125.999), 0, (    -125.999)},  // H3 LCI

                          {"V1", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},
                         };


  detector* pD[NIFO_MAX];
  for(int i=0; i<nIFO; i++) pD[i] = new detector(detParms[i]);   // user define detector
  CWB::mdc MDC(nIFO,pD); 

  // -------------------------------------------------------- 
  // define waveforms
  // -------------------------------------------------------- 
  int nWF=15;
  double F[nWF] = {70,100,153,235,361,554,849,945,1053,1172,1304,1451,1615,1797,2000};
  double Q[nWF] = {8.9,8.9,8.9,8.9,8.9,8.9,8.9,9,9,9,9,9,9,9,9};

  char wf_name[256];
  waveform wf;
  wf.par.resize(2);
  for (int n=0;n<nWF;n++) {
    wf.type = MDC_SG;
    sprintf(wf_name,"SG%dQ%1.1f",F[n],Q[n]);
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
  vector<mdcpar> par; par.resize(3);
  par[0].name="entries";par[0].value=100000;  // pool of events
  par[1].name="rho_min";par[1].value=1;         // min rho // Mpc
  par[2].name="rho_max";par[2].value=100;       // max rho // Mpc
  MDC.SetSkyDistribution(MDC_RANDOM,par);
*/
/*
  vector<mdcpar> par; par.resize(1);
  par[0].name="entries";par[0].value=100000;  // number of entries
  MDC.SetSkyDistribution(MDC_MNGD,par);
*/
/*
  vector<mdcpar> par; par.resize(1);
  par[0].name="distance_thr";par[0].value=100;  // distance max = 100 Mpc
  MDC.SetSkyDistribution(MDC_GWGC,"../data/GWGCCatalog_Rev1d2.txt",par);
*/

  vector<mdcpar> par; par.resize(3);
  par[0].name="theta";par[0].value=30;     // theta
  par[1].name="phi";  par[1].value=60;     // phi
  par[2].name="psi";  par[2].value=90;     // psi
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par);

  bool   log      = true;
  size_t gps      = 931158100;
  size_t length   = 600;
  MDC.WriteFrameFile(frDir, frLabel, gps, length, log);

  exit(0);
}
