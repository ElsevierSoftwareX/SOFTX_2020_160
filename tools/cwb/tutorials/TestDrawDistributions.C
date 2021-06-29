// Draw sky distributions
// Author : Gabriele Vedovato
//
#include <vector>

void TestDrawDistributions() {

  #define N_IFO 3

  //TString ifo[N_IFO]={"L1","H1","V1"};
  TString ifo[N_IFO]={"L1","H1","VIRGO"};
  CWB::mdc MDC(N_IFO,ifo); 

  vector<mdcpar> par(3);
/*
  par.resize(3);
  par[0].name="entries"; par[0].value=100000;  
  par[1].name="rho_min"; par[1].value=5;  
  par[2].name="rho_max"; par[2].value=10;  
  MDC.SetSkyDistribution(MDC_RANDOM,par);
*/

  //par.resize(1);
  //par[0].name="entries"; par[0].value=100000;  
  //MDC.SetSkyDistribution(MDC_MNGD,par);
/*
  par.resize(1);
  par[0].name="distance_thr"; par[0].value=100;  
  MDC.SetSkyDistribution(MDC_GWGC,"../data/GWGCCatalog_Rev1d7.txt",par);
*/
/*
  par.resize(3);
  par[0].name="theta"; par[0].value=30;  
  par[1].name="phi"; par[1].value=60;  
  par[2].name="rho"; par[2].value=1;  
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par);
*/

  MDC.SetSkyDistribution(MDC_LOGFILE,"GHLTV-GWGC_v1-968600000-300000.log",par);


  //MDC.DrawSkyDistribution("skyplot","","geographic",2,true);
  MDC.DrawSkyDistribution("skyplot","","celestial",2,true);
  //MDC.DrawSkyDistribution("skyplot","","celestial",2,false);
  //MDC.DrawSkyDistribution("skyplot","","geographic",2,false);

  //exit(0);
}
