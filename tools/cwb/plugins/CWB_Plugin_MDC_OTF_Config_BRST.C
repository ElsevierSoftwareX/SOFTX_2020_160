//!NOISE_MDC_SIMULATION
// Config Plugin to generate injected 'on the fly' Burst MDC

{
  #define GWGCCatalog "$CWB_GWAT/data/GWGCCatalog_Rev1d8.txt"

  #define SKY_RANDOM
  //#define SKY_FIX
  //#define SKY_GWGC

  cout << "Execute CWB_Plugin_MDC_OTF_Config_BRST.C ..." << endl;

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int seed = (*net)->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  // declaration of the GA waveforms
  par.resize(1);
  double F_GA[4]   = {0.0001,0.001,0.0025,0.004};
  for(int n=0;n<4;n++) {
    par[0].name="duration"; par[0].value=F_GA[n];
    MDC->AddWaveform(MDC_GA, par);
  }

  // declaration of the SGQ3 waveforms
  par.resize(2);
  double F_Q3[4] = {70,235,849,1615};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q3[n];
    par[1].name="Q"; par[1].value=3.;
    MDC->AddWaveform(MDC_SG, par);
  }

  // declaration of the SGQ9 waveforms
  par.resize(2);
  double F_Q8d9[7] = {70,100,153,235,361,554,849};
  for(int n=0;n<7;n++) {
    par[0].name="frequency"; par[0].value=F_Q8d9[n];
    par[1].name="Q"; par[1].value=8.9;
    MDC->AddWaveform(MDC_SG, par);
  }
  double F_Q9[4]   = {1053,1304,1615,2000};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q9[n];
    par[1].name="Q"; par[1].value=9.;
    MDC->AddWaveform(MDC_SG, par);
  }

  // declaration of the SGQ100 waveforms
  par.resize(2);
  double F_Q100[4] = {70,235,849,1615};
  for(int n=0;n<4;n++) {
    par[0].name="frequency"; par[0].value=F_Q100[n];
    par[1].name="Q"; par[1].value=100.;
    MDC->AddWaveform(MDC_SG, par);
  }

  // declaration of the WNB waveforms
  par.resize(6);
  double F_WNB[5] = {100,250,1000,1000,1000};
  double B_WNB[5] = {100,100,10,1000,1000};
  double D_WNB[5] = {0.1,0.1,0.1,0.1,0.01};
  for(int n=0;n<5;n++) {
    par[0].name="frequency"; par[0].value=F_WNB[n];
    par[1].name="bandwidth"; par[1].value=B_WNB[n];
    par[2].name="duration";  par[2].value=D_WNB[n];
    for(int m=0;m<30;m++) {
      par[3].name="pseed";     par[3].value=seed+100000+n*100+m;
      par[4].name="xseed";     par[4].value=seed+100001+n*100+m;
      par[5].name="mode";      par[5].value=1;  // simmetric
      MDC->AddWaveform(MDC_WNB, par);
    }
  }

  MDC->Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjHrss(2.5e-21);	// injection hrss (the factors array in user_parameter.C file is used to rescale it)
  MDC->SetInjRate(0.01);		// injection rate (injections per sec)
  MDC->SetInjJitter(10.0);	// time jitter in sec

  // ---------------------------------------------------------------------------------------
  // define sky distribution
  // for a full description of the parameters see the CWB::mdc::SetSkyDistribution in mdc.cc
  // ---------------------------------------------------------------------------------------

  // this example show how to set a random sky distribution
#ifdef SKY_RANDOM
  par.resize(3);
  par[0].name="entries"; par[0].value=100000;	// number of entries randomly generated
  par[1].name="rho_min"; par[1].value=0;	// the source hrss is not rescaled
  par[2].name="rho_max"; par[2].value=0;	// the source hrss is not rescaled
  MDC->SetSkyDistribution(MDC_RANDOM,par,seed);	// random sky distribution (theta,phi,psi)
#endif

  // this example show how to set fix the injection parameters 
#ifdef SKY_FIX
  par.resize(6);
  par[0].name="theta"; par[0].value=40;
  par[1].name="phi";   par[1].value=70;
  par[2].name="psi";   par[2].value=0;
  par[3].name="rho";   par[3].value=1;
  par[4].name="gps";   par[4].value=931158300;
  par[5].name="iota";  par[5].value=0;          // the angle iota (deg) is the inclination of the system
                                                // which originates the burst with respect to the line of sight
  // theta,phi are in earth greographical coordinates (latitude, longitude)
  MDC->SetSkyDistribution(MDC_EARTH_FIX,par);	
// theta,phi are in celestial coordinates (DEC, RA)
//  MDC->SetSkyDistribution(MDC_CELESTIAL_FIX,par);
#endif

  // this example show how to use the Gravitational Wave Galaxy Catalog distribution 
#ifdef SKY_GWGC
  par.resize(1);
  par[0].name="distance_thr"; par[0].value=50000;  // Max source distance -> 50000 kPc
  MDC->SetSkyDistribution(MDC_GWGC,gSystem->ExpandPathName(GWGCCatalog),par,seed);
#endif

}
