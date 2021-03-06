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
  fHigh = 512.;       // high frequency of the search

  levelR   = 4;
  l_low    = 3;       // low frequency resolution level
  l_high   = 8;       // high frequency resolution level

  strcpy(wdmXTalk,"wdmXTalk/OverlapCatalog-ilLev3-hLev9-iNu6-P10.xbin");  // 1KHz

  healpix=7;

  bpp    = 0.001;
  subnet = 0.5;
  subcut = 0.0;
  netRHO = 5.0;
  netCC  = 0.5;
  Acore  = 1.7;
  Tgap   = 0.2;
  Fgap   = 128.0;
  delta  = 0.5;
  cfg_gamma  = -1.0;	
  LOUD   = 300;

  pattern = 10;

  iwindow = 30.;

  //precision=GetPrecision(100,5);

//  nSky = -99;      // save pixels skymap probability pixels up to cumulative prob < 0.99
  nSky=196608;	   // save all pixels

  //simulation
  nfactor = 1;
  simulation = 0;

  strcpy(channelNamesRaw[0],"L1:GWOSC-4KHZ_R1_STRAIN");
  strcpy(channelNamesRaw[1],"H1:GWOSC-4KHZ_R1_STRAIN");

  strcpy(frFiles[0],"input/L1_frames.in");
  strcpy(frFiles[1],"input/H1_frames.in");

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=8;
  dqfile dqf[8]={
                    {"L1" ,"input/cwb_period.txt",   CWB_CAT0, 0., false, false},
                    {"H1" ,"input/cwb_period.txt",   CWB_CAT0, 0., false, false},
                    {"L1" ,"../DQ/L1_cat0.txt",      CWB_CAT0, 0., false, false},
                    {"H1" ,"../DQ/H1_cat0.txt",      CWB_CAT0, 0., false, false},
                    {"L1" ,"../DQ/L1_cat1.txt",      CWB_CAT1, 0., false, false},
                    {"H1" ,"../DQ/H1_cat1.txt",      CWB_CAT1, 0., false, false},
                    {"L1" ,"../DQ/L1_cat2.txt",      CWB_CAT2, 0., false, false},
                    {"H1" ,"../DQ/H1_cat2.txt",      CWB_CAT2, 0., false, false}
                   };
  for(int i=0;i<8;i++) DQF[i]=dqf[i];	


  plugin = TMacro("macro/CWB_Plugin_O3aConditioning_WF.C");        	// Macro source
  //plugin.SetTitle("macro/CWB_Plugin_WF_C.so");
  configPlugin = TMacro("macro/CWB_Plugin_Config.C");   // Macro config

  TString optwf = "";                  // NOTE : add space at the end of each line

  //optwf += "wf_output_disable=root ";       // disable output root file (to be used when QLveto is enabled)
  optwf += "wf_output_enable=root ";        // enable output root file (to be used when QLveto is enabled)
  optwf += "wf_output_disable=inj ";        // disable save injection to the output root file
  optwf += "wf_output_enable=rec ";         // enable  save reconstructed waveform to the output root file
  //optwf += "wf_output_disable=wht ";        // disable save whitened data to the output root file
  //optwf += "wf_output_disable=dat ";        // disable save rec+null data to the output root file
  //optwf += "wf_output_disable=nul ";        // disable save null data to the output root file
  optwf += "wf_output_enable=wht ";         // enable save whitened data to the output root file
  optwf += "wf_output_enable=dat ";         // enable save rec+null data to the output root file
  optwf += "wf_output_enable=nul ";         // enable save null data to the output root file

  optwf += "wf_inj_tstep=150.000 ";         // is the injection step time (used only for PE simulations, must be >0)

  optwf += "sn_strain_file_name=input/L1_ASD.txt ";    // add L1 strain file name to first detector
  optwf += "sn_strain_file_name=input/H1_ASD.txt ";    // add H1 strain file name to second detector
  optwf += "sn_seed=1 ";                    // add seed for random noise generation to L1 detector
  optwf += "sn_seed=2 ";                    // add seed for random noise generation to H1 detector

  strcpy(parPlugin,optwf.Data());           // set WF plugin parameters
  strcpy(comment,"GWOSC");
}
