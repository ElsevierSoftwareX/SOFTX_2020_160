//!NOISE_MDC_SIMULATION
// Config Plugin to generate injected 'on the fly' eBBH MDC

{                                                                                                             
  #define EBBH_LIST

  cout << "Execute CWB_Plugin_MDC_OTF_Config_eBBH.C ..." << endl;                                                    

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int seed = (*net)->nRun;  // WARNING : seed must be the same for each detector within the same job

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  char wf_name[256];
  waveform wf;      
  vector<mdcpar> par; 

#ifdef EBBH_LIST
  // list of eBBH parameters are read from ascii file 
  // format : one event for each line
  // -> m1 m2 rp0 e0 (hp,hx are generated 'On The Fly')
  // see macro trunk/tools/eBBH/tutorials/CreateListEBBH.C
  par.resize(1);
  par[0].name="macro/eBBH_3.5e+06_5_25_2.lst";  // ascii eBBH list file
  MDC->AddWaveform(MDC_EBBH, par);              
#else
  // list of eBBH parameters are read from root tree 
  // format : one event for each entry 
  // -> id m1 m2 rp0 e0 hp hx
  // see macro trunk/tools/eBBH/tutorials/CreateListEBBH.C
  // see macro trunk/tools/eBBH/tutorials/List2RootEBBH.C
  par.resize(2);                               
  par[0].name="macro/eBBH.root";        // root eBBH list file       
  par[1].name="rp0>26";			// selection cut (optional)
  MDC->AddWaveform(MDC_EBBH, par);
#endif
  MDC->Print(0);

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjRate(0.005);   // No of injections per second
  MDC->SetInjJitter(10.0);  // +- injection time

  // --------------------------------------------------------
  // define sky distribution
  // source are uniformely randomly generated between 1Mpc and 15Mpc
  // --------------------------------------------------------
  par.resize(4);
  par[0].name="entries"; par[0].value=200000;
  par[1].name="rho_min"; par[1].value=1000;      // min rho // Kpc
  par[2].name="rho_max"; par[2].value=15000;     // max rho // Kpc
  par[3].name="iota"   ; par[3].value=-1;        // the angle iota (deg) is the inclination of the system
                                                 // which originates the burst with respect to the line of sight
  MDC->SetSkyDistribution(MDC_RANDOM,par,seed);
}

