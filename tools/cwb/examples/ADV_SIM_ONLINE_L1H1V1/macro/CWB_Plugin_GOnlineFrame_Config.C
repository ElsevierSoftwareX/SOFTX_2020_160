{
  cout << "-----> macro/CWB_Plugin_GOnlineFrame_Config.C" << endl;

  int seed = net->nRun;  // WARNING : seed must be the same for each detector within the same job

  int  FRLEN = 4;        // output frame len (sec) - must be a submultiple of segLen

  CWB::mdc MDC(net);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  if(mdc_type=="signal") {

    // --------------------------------------------------------------------
    // configuration for injected SIGNALS
    // --------------------------------------------------------------------

    par.resize(2);

    double SF_Q3[2] = {235,1053};
    for(int n=0;n<2;n++) {
      par[0].name="frequency"; par[0].value=SF_Q3[n];
      par[1].name="Q"; par[1].value=3.;
      MDC.AddWaveform(MDC_SG, par);
    }

    double SF_Q9[2] = {235,1053};
    for(int n=0;n<2;n++) {
      par[0].name="frequency"; par[0].value=SF_Q9[n];
      par[1].name="Q"; par[1].value=9.;
      MDC.AddWaveform(MDC_SG, par);
    }

    double SF_Q100[2] = {235,1053};
    for(int n=0;n<2;n++) {
      par[0].name="frequency"; par[0].value=SF_Q100[n];
      par[1].name="Q"; par[1].value=100.;
      MDC.AddWaveform(MDC_SG, par);
    }

    par.resize(6);

    double SF_WNB[2] = {250,1000};
    double SB_WNB[2] = {100,1000};
    double SD_WNB[2] = {0.1,0.1};
    for(int n=0;n<2;n++) {
      par[0].name="frequency"; par[0].value=SF_WNB[n];
      par[1].name="bandwidth"; par[1].value=SB_WNB[n];
      par[2].name="duration";  par[2].value=SD_WNB[n];
      for(int m=0;m<30;m++) {
        par[3].name="pseed";     par[3].value=seed+100000+n*100+m;
        par[4].name="xseed";     par[4].value=seed+100001+n*100+m;
        par[5].name="mode";      par[5].value=1;  // simmetric
        MDC.AddWaveform(MDC_WNB, par);
      }
    }

    MDC.Print();

    // --------------------------------------------------------
    // define injection parameters
    // --------------------------------------------------------
    MDC.SetInjHrss(2.5e-21);	// injection hrss
    MDC.SetInjRate(0.01);		// number of injections per sec
    MDC.SetInjJitter(10.0);	// injection time jitter [+/- Inj_Jitter_Time]

    // --------------------------------------------------------
    // define sky distribution
    // --------------------------------------------------------
    vector<mdcpar> par(3);
    par[0].name="entries"; par[0].value=100000;	// number of random events to choose from
    // set distribution in rho (Kpc) from rho_min to rho_max according to 1/(rho*rho) power low 
    // the hrss is rescaled with the inverse of distance (standard candle @ 10Kpc)
    // if rho_min=rho_max=0 no rescaled is performed
    par[1].name="rho_min"; par[1].value=0;
    par[2].name="rho_max"; par[2].value=0;
    // theta,phi,psi are random
    // for more options see CWB::mdc::SetSkyDistribution infos in msd.cc file
    MDC.SetSkyDistribution(MDC_RANDOM,par,seed);

  } else if(mdc_type=="glitch") { 

    // --------------------------------------------------------------------
    // configuration for injected GLITCHES
    // --------------------------------------------------------------------

    par.resize(2);

    double GF_Q3[1] = {100};
    for(int n=0;n<1;n++) {
      par[0].name="frequency"; par[0].value=GF_Q3[n];
      par[1].name="Q"; par[1].value=3.;
      MDC.AddWaveform(MDC_SG, par);
    }

    double GF_Q9[1] = {100};
    for(int n=0;n<1;n++) {
      par[0].name="frequency"; par[0].value=GF_Q9[n];
      par[1].name="Q"; par[1].value=9.;
      MDC.AddWaveform(MDC_SG, par);
    }

    double GF_Q100[1] = {100};
    for(int n=0;n<1;n++) {
      par[0].name="frequency"; par[0].value=GF_Q100[n];
      par[1].name="Q"; par[1].value=100.;
      MDC.AddWaveform(MDC_SG, par);
    }

    par.resize(6);

    double GF_WNB[1] = {100};
    double GB_WNB[1] = {100};
    double GD_WNB[1] = {0.1};
    for(int n=0;n<1;n++) {
      par[0].name="frequency"; par[0].value=GF_WNB[n];
      par[1].name="bandwidth"; par[1].value=GB_WNB[n];
      par[2].name="duration";  par[2].value=GD_WNB[n];
      for(int m=0;m<30;m++) {
        par[3].name="pseed";     par[3].value=seed+100000+n*100+m;
        par[4].name="xseed";     par[4].value=seed+100001+n*100+m;
        par[5].name="mode";      par[5].value=1;  // simmetric
        MDC.AddWaveform(MDC_WNB, par);
      }
    }

    MDC.Print();

    CWB::Toolbox TB;

    // read distribution configuration from user custom file
    if(IFO=="L1") {
      TB.checkFile("input/L1_glitches.lst");
      MDC.SetSkyDistribution(MDC_CUSTOM,"input/L1_glitches.lst",seed+1);
    }
    if(IFO=="H1") {
      TB.checkFile("input/H1_glitches.lst");
      MDC.SetSkyDistribution(MDC_CUSTOM,"input/H1_glitches.lst",seed+2);
    }
    if(IFO=="V1") {
      TB.checkFile("input/V1_glitches.lst");
      MDC.SetSkyDistribution(MDC_CUSTOM,"input/V1_glitches.lst",seed+3);
    }
  } 
}
