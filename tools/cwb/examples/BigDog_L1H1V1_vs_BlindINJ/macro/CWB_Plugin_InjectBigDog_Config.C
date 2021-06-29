#include "CWB_Plugin.h"

void CWB_PluginConfig() {

  //!NOISE_MDC_SIMULATION
  // Config Plugin to injected 'on the fly' BigDog MDC

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int* gIFACTOR;
  CWB_PLUGIN_IMPORT(int*,gIFACTOR);

  int seed = 100*((*net)->nRun)+*gIFACTOR; 
  char SEED[64]; sprintf(SEED,"%d",seed);

  cout << "gIFACTOR : " << *gIFACTOR << " net->nRun : " << (*net)->nRun << " SEED " << SEED << endl; 

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  CWB::config** cfg;
  CWB_PLUGIN_IMPORT(CWB::config**,cfg);

  // ---------------------------------
  // set BigDog inspiral parameters
  // ---------------------------------
  TString inspOptions="--xml config/CBC_BLINDINJ_968654558_adj.xml ";
          inspOptions+= "--dir "+TString((*cfg)->nodedir);

  MDC->SetInspiral("BigDog",inspOptions);
}
