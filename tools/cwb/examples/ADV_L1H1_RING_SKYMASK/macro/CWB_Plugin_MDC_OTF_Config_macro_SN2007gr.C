#include "CWB_Plugin.h"

#define SN2007gr_GPS            870775372
#define SN2007gr_RA             40.8666
#define SN2007gr_DEC            37.3458

void CWB_PluginConfig() {

  cout << "-----> macro/CWB_Plugin_MDC_OTF_Config_macro_SN2007gr.C" << endl;

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
  par[0].name="theta"; par[0].value=SN2007gr_DEC;
  par[1].name="phi";   par[1].value=SN2007gr_RA;
  par[2].name="gps";   par[2].value=SN2007gr_GPS;
  MDC.SetSkyDistribution(MDC_CELESTIAL_FIX,par,seed);    // theta,phi are in celestial coordinates

}

