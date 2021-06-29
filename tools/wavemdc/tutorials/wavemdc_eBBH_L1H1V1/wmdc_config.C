{
  #include <vector>

  // --------------------------------------------------------
  // define job parameters
  // --------------------------------------------------------
  TString frDir       = "frames";
  TString frLabel     = "wavemdc_eBBH";

  // define frame segment list
  TString segmentList = "Segments/eBBH_SegmentsList.txt";

  // select segments
  int jobmin  =   1;    // start segment
  int jobmax  =   1;    // end segment
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
  vector<TString> chName(nIFO);
  chName[0]="L1:EBBH-STRAIN_INJ";
  chName[1]="H1:EBBH-STRAIN_INJ";
  chName[2]="V1:EBBH-STRAIN_INJ";

  // ---------------------------------
  // Set eBBH parms
  // ---------------------------------

  vector<mdcpar> par;

  par.resize(1);
  par[0].name="macro/eBBH_3.5e06_10_25_2.lst";
  MDC.AddWaveform(MDC_EBBH, par);

/*
  par.resize(2);
  par[0].name="macro/eBBH_3.5e06_10_25_2.root";
//  par[1].name="rp0>26";
  MDC.AddWaveform(MDC_EBBH, par);
*/

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC.SetInjHrss(2.5e-21);
  MDC.SetInjRate(0.01);
  MDC.SetInjJitter(10.0);
  MDC.SetInjLength(20.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  vector<mdcpar> par(4);
  par[0].name="entries"; par[0].value=100000;
  par[1].name="rho_min"; par[1].value=0;
  par[2].name="rho_max"; par[2].value=0;
  par[3].name="iota";    par[3].value=-1;          // random
  MDC.SetSkyDistribution(MDC_RANDOM,par);

}
