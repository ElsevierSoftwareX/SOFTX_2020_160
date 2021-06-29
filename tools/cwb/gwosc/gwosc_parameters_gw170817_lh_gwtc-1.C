{
  strcpy(analysis,"2G");

  nIFO = 2;
  cfg_search = 'r';	
  optim=false;

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(refIFO,"L1");

  //lags
  lagSize = 1;
  lagStep = 1.;
  lagOff  = 0;
  lagMax  = 0;

  //superlags
  slagSize = 0;		     // number of super lags (simulation=1) - if slagSize=0 -> Standard Segments
  slagMin  = 0;
  slagMax  = 0;
  slagOff  = 0;

  //jobs
  segLen  = 1200;
  segMLS  = 600;
  segTHR  = 200;
  segEdge = 10;

  //frequency
  fLow  = 16.;        // low frequency of the search
  fHigh = 1024.;      // high frequency of the search

  levelR   = 3;
  l_low    = 4;       // low frequency resolution level
  l_high   = 10;      // high frequency resolution level

  strcpy(wdmXTalk,"wdmXTalk/OverlapCatalog16-1024.bin");

  healpix=7;

  bpp    = 0.001;
  subnet = 0.5;
//  subnet = 0.0;		// Used with Tukey Window
  subcut = 0.0;
  netRHO = 5.5;
  netCC  = 0.5;
  Acore  = 1.7;
  Tgap   = 0.2;
  Fgap   = 128.0;
  delta  = 0.5;
  cfg_gamma  = -1.0;	
  LOUD   = 300;

  pattern = 5;

  iwindow = 30.;

  //precision=GetPrecision(100,5);

//  nSky = -99;      // save pixels skymap probability pixels up to cumulative prob < 0.99
  nSky=196608;     // save all pixels

  //simulation
  nfactor = 1;
  simulation = 0;

  strcpy(channelNamesRaw[0],"L1:DCH-CLEAN_STRAIN_C02_T1700406_v3");	// https://dcc.ligo.org/LIGO-T1700406/public - L1 Glitch Removed
//  strcpy(channelNamesRaw[0],"L1:GWOSC-4KHZ_R1_STRAIN");	// Used with Tukey Window
  strcpy(channelNamesRaw[1],"H1:GWOSC-4KHZ_R1_STRAIN");

  strcpy(frFiles[0],"input/L1_frames.in");
  strcpy(frFiles[1],"input/H1_frames.in");

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
//  #define NDFQ 8	// Used with Tukey Window	
  #define NDFQ 7
  nDQF=NDFQ;
  dqfile dqf[NDFQ]={
                    {"L1" ,"input/cwb_period.txt",   CWB_CAT0, 0., false, false},
                    {"H1" ,"input/cwb_period.txt",   CWB_CAT0, 0., false, false},
                    {"L1" ,"../DQ/L1_cat0.txt",      CWB_CAT0, 0., false, false},
                    {"H1" ,"../DQ/H1_cat0.txt",      CWB_CAT0, 0., false, false},
                    {"L1" ,"../DQ/L1_cat1.txt",      CWB_CAT1, 0., false, false},
                    {"H1" ,"../DQ/H1_cat1.txt",      CWB_CAT1, 0., false, false},
                    //{"L1" ,"../DQ/L1_cat2.txt",      CWB_CAT2, 0., false, false},
                    //{"L1" ,"input/gw170817_twwin_cat2.txt",          CWB_CAT2, 0., false, false},	// Used with Tukey Window
                    {"H1" ,"../DQ/H1_cat2.txt",      CWB_CAT2, 0., false, false}
                   };
  for(int i=0;i<NDFQ;i++) DQF[i]=dqf[i];	

  plugin = TMacro("macro/CWB_Plugin_WF.C");      
  configPlugin = TMacro("macro/CWB_Plugin_Config.C");   // Macro config

/*
  plugin = TMacro("macro/CWB_Plugin_TukeyWindow_WF.C");        // Used with Tukey Window

  TString opttw = "";            		// NOTE : add space at the end of each line

  opttw += "tw_gps_time=1187008881.38 ";        // TW central GPS time
  opttw += "tw_width=0.34 ";                    // TW width
  opttw += "tw_range=0.14 ";                    // TW range
  opttw += "tw_ced_spectrogram_zmax=1e-3 ";     // TW set zmax in CED spectrograms
  opttw += "tw_dump=false ";                    // TW dump ascii file

  strcpy(parPlugin,opttw.Data());               // set TW plugin parameters
*/

  TString optwf = "";                  		// NOTE : add space at the end of each line

  //optwf += "wf_output_disable=root ";       	// disable output root file (to be used when QLveto is enabled)
  optwf += "wf_output_enable=root ";        	// enable output root file (to be used when QLveto is enabled)
  optwf += "wf_output_disable=inj ";        	// disable save injection to the output root file
  optwf += "wf_output_enable=rec ";         	// enable  save reconstructed waveform to the output root file
  //optwf += "wf_output_disable=wht ";        	// disable save whitened data to the output root file
  //optwf += "wf_output_disable=dat ";        	// disable save rec+null data to the output root file
  //optwf += "wf_output_disable=nul ";        	// disable save null data to the output root file
  optwf += "wf_output_enable=wht ";         	// enable save whitened data to the output root file
  optwf += "wf_output_enable=dat ";         	// enable save rec+null data to the output root file
  optwf += "wf_output_enable=nul ";         	// enable save null data to the output root file

  optwf += "wf_inj_tstep=150.000 ";             // is the injection step time (used only for PE simulations, must be >0)
 
  optwf += "sn_strain_file_name=input/L1_ASD.txt ";    // add L1 strain file name to first detector
  optwf += "sn_strain_file_name=input/H1_ASD.txt ";    // add H1 strain file name to second detector
  optwf += "sn_seed=1 ";                    // add seed for random noise generation to L1 detector
  optwf += "sn_seed=2 ";                    // add seed for random noise generation to H1 detector

  strcat(parPlugin,optwf.Data());           	// set WF plugin parameters
  strcpy(comment,"GWOSC");
}
