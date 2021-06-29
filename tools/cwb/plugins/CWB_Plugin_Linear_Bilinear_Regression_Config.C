{
  //!REGRESSION
  // Config Plugin for linear/bilinear regression analysis for one detector

  //---------------------------------------------------------------------
  // REGRESSION MAIN
  //---------------------------------------------------------------------

  bool EXIT_AFTER_REGRESSION            = true;         // if true cwb terminate after the regression analysis
  bool APPLY_LINEAR_REGRESSION          = true;         // if true then linear regression is applied
  bool APPLY_BILINEAR_REGRESSION        = false;         // if true then blinear regression is applied;
  bool SAVE_INFRAME                     = false;         // if true the input data are saved to frame files
  bool SAVE_OUTFRAME                    = false;         // if true the cleaned data are saved to frame files
  bool SAVE_PSD_PLOT                    = true;         // if true psd cleaned vs noisy data are saved to png
  bool SAVE_EIGEN_PLOT                  = true;         // if true the regression filter eigenvaules are saved to png
  int  CUT_LOW_FREQ                     = 32;           // if > 0 the data in [0:CUT_LOW_FREQ] are set to 0

  TString OFRDIR     			= "oframes";
  TString FRNAME     			= "LHO_4k";
  TString FRLABEL     			= "H-H1_LDAS_C02_L2";

  bool DISABLE_MDC_FROM_FRAMES          = true;         // if true then disable read mdc from frames
  bool DISABLE_MDC_FROM_PLUGIN          = true;         // if true then disable read mdc from plugin

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS FOR H1
  //---------------------------------------------------------------------

  TString FRLIST_WITNESS             	= "input/H1_RDS_R_L1_937473000_947017000.frl";
  TString CHLIST_LINEAR     		= "input/S6B_H1_Linear_Regression_Channel_List.txt";
  TString CHLIST_BILINEAR       	= "input/S6B_H1_Bilinear_Regression_Channel_List.txt";

  TString CHNAME_LINEAR           	= "H0:PEM-COIL_MAGX";

  int     RESAMPLING_INDEX              = 11;           // resample target/witness channels to pow(2,RESAMPLING_INDEX)

  double  WFLOW                         = 0;            // lower frequency used to made the bilinear channels
  double  WFHIGH                        = 10;           // high  frequency used to made the bilinear channels

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS
  //---------------------------------------------------------------------

  double LAYER_WIDTH                    = 1.0;          // frequency layer resolution used in the regression analysis
  double fPOWERLINE                     = 60.0;         // powerline frequency (Hz)
  int    lPOWERLINE                     = 1;            // low power line harmonic (lPOWERLINE*fPOWERLINE)
  int    hPOWERLINE                     = 5;            // high power line harmonic (hPOWERLINE*fPOWERLINE)
  double FWIDTH                         = 5;            // frequency width used by regression

  // PARAMETERS FOR LINEAR REGRESSION

  int    L_NFILTER                      = 5;            // half-size length of a unit filter (setFilter)
  double L_APPLY_THRESHOLD              = 0.2;          // threshold used in apply
  double L_SOLVE_THRESHOLD              = 0.0;          // eigenvalue threshold (solve)
  double L_SOLVE_NEIGEN_PER_LAYER       = 0;            // number of selected eigenvalues (solve)
  char   L_SOLVE_REGULATOR              = 'h';          // regulator (solve)

  // PARAMETERS FOR BILINEAR REGRESSION

  int    B_NFILTER                      = 5;            // half-size length of a unit filter (setFilter)
  double B_APPLY_THRESHOLD              = 0.2;          // threshold used in apply
  double B_SOLVE_THRESHOLD              = 0.0;          // eigenvalue threshold (solve)
  double B_SOLVE_NEIGEN_PER_LAYER       = 0;            // number of selected eigenvalues (solve)
  char   B_SOLVE_REGULATOR              = 'h';          // regulator (solve)

  //---------------------------------------------------------------------
  // MDC 
  //---------------------------------------------------------------------

  if(type==CWB_PLUGIN_MDC && !DISABLE_MDC_FROM_PLUGIN) {

    network** net;
    CWB_PLUGIN_IMPORT(network**,net);

    int seed = (*net)->nRun;  // WARNING : seed must be the same for each detector within the same job

    CWB::mdc* MDC;
    CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

    char wf_name[256];
    waveform wf;
    vector<mdcpar> par;

    for(int n=1;n<=5;n++) {

      par.resize(2);
      par[0].name="frequency"; par[0].value=n*60.-2;
      par[1].name="Q"; par[1].value=100.;
      MDC->AddWaveform(MDC_SGC, par);

      par.resize(2);
      par[0].name="frequency"; par[0].value=n*60.;
      par[1].name="Q"; par[1].value=100.;
      MDC->AddWaveform(MDC_SGC, par);

      par.resize(2);
      par[0].name="frequency"; par[0].value=n*60.+2;
      par[1].name="Q"; par[1].value=100.;
      MDC->AddWaveform(MDC_SGC, par);

    }

    //MDC->Print();

    // --------------------------------------------------------
    // define injection parameters
    // --------------------------------------------------------
    MDC->SetInjHrss(2.5e-21);
    MDC->SetInjRate(0.0333333);
    MDC->SetInjJitter(1.0);

    // --------------------------------------------------------
    // define sky distribution
    // --------------------------------------------------------
    par.resize(3);
    par[0].name="entries";par[0].value=100000;    // pool of events
    par[1].name="rho_min";par[1].value=1;         // min rho // Mpc
    par[2].name="rho_max";par[2].value=1;         // max rho // Mpc
    MDC->SetSkyDistribution(MDC_RANDOM,par,seed);
  }
}
