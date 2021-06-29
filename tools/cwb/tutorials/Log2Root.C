{
  //
  // Convert text log mdc file to root file
  // Author : Gabriele Vedovato

  #define N_IFO 3
 
  TString ifo[N_IFO]={"L1","H1","VIRGO"};

  CWB::mdc MDC(N_IFO,ifo);
  vector<mdcpar> par;

  MDC.SetSkyDistribution(MDC_LOGFILE,"GHLTV-GWGC_v1-968600000-300000.log",par);
  MDC.DumpLog("GHLTV-GWGC_v1-968600000-300000.root","CIAO");


  gSystem->Exit(0);
}
