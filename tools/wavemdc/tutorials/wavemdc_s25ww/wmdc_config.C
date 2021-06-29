{
  #include <vector>

  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "TestS25WW";

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
  TString ifo[nIFO]={"L1","H1","G1"};
  CWB::mdc MDC(nIFO,ifo); 

  // --------------------------------------------------------
  // define channel names
  // --------------------------------------------------------
  vector<TString> chName;

  // --------------------------------------------------------
  // Add s25WW
  // --------------------------------------------------------

  TString fName = "input/s25WW.h.dat";
  wavearray<double> x;
  MDC.SetInjLength(1.5);
  MDC.ReadWaveform(x, fName);

  waveform wf;
  wf.name = "s25WW";
  wf.hp = x;
  wf.hx = wf.hp;
  wf.hx = 0;
  MDC.AddWaveform(wf);

  MDC.Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC.SetInjHrss(9.22693e-21);
  MDC.SetInjRate(0.0333333);
  MDC.SetInjJitter(10.0);

  // --------------------------------------------------------
  // define MNGD sky distribution
  // --------------------------------------------------------
  vector<mdcpar> par(1);
  par[0].name="entries"; par[0].value=100000;
  MDC.SetSkyDistribution(MDC_MNGD,par);
}
