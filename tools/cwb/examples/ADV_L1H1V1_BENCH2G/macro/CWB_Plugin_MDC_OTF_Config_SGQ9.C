// THIS IS THE NEW VERSION OF CONFIG PLUGIN !!!
// THE NEW VERSION IS REQUIRED BY ROOT6
// IN ROOT5 IS STILL POSSIBLE TO USE OLD & NEW VERSION 
// AS COMPARISON THE OLD VERSION IS REPORTED BELOW

#include "CWB_Plugin.h"

void CWB_PluginConfig() {

  cout << "-----> macro/CWB_Plugin_MDC_OTF_Config_SGQ9.C" << endl;

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int seed = (*net)->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  par.resize(2);

  double F_Q9[1] = {150};
  for(int n=0;n<1;n++) {
    par[0].name="frequency"; par[0].value=F_Q9[n];
    par[1].name="Q"; par[1].value=9.;
    MDC->AddWaveform(MDC_SG, par);
  }

  MDC->Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjHrss(2.5e-21);
  MDC->SetInjRate(0.01);
  MDC->SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  par.resize(3);
  par[0].name="entries"; par[0].value=100000;
  par[1].name="rho_min"; par[1].value=0;
  par[2].name="rho_max"; par[2].value=0;
  MDC->SetSkyDistribution(MDC_RANDOM,par,seed);

}

// THIS IS THE OLD VERSION OF CONFIG PLUGIN !!!
/*
{                                                                                                                                        
  cout << "-----> macro/CWB_Plugin_MDC_OTF_Config_SGQ9.C" << endl;                                                                       

  int seed = net->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc MDC(net);

  char wf_name[256];
  waveform wf;      
  vector<mdcpar> par; 

  par.resize(2);

  double F_Q9[1] = {150};
  for(int n=0;n<1;n++) {
    par[0].name="frequency"; par[0].value=F_Q9[n];
    par[1].name="Q"; par[1].value=9.;
    MDC.AddWaveform(MDC_SG, par);
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
  par.resize(3);
  par[0].name="entries"; par[0].value=100000;
  par[1].name="rho_min"; par[1].value=0;
  par[2].name="rho_max"; par[2].value=0;
  MDC.SetSkyDistribution(MDC_RANDOM,par,seed);
}
*/
