//#define TEST_LAGS

{
  strcpy(analysis,"2G");

  nIFO = 3;
  cfg_search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");
  strcpy(refIFO,"L1");

  strcpy(channelNamesRaw[0],"L1:SIM-STRAIN");
  strcpy(channelNamesRaw[1],"H1:SIM-STRAIN");
  strcpy(channelNamesRaw[2],"V1:SIM-STRAIN");

  strcpy(channelNamesMDC[0],"L1:GW-H");
  strcpy(channelNamesMDC[1],"H1:GW-H");
  strcpy(channelNamesMDC[2],"V1:GW-H");

  strcpy(frFiles[0],"input/L1.frames");
  strcpy(frFiles[1],"input/H1.frames");
  strcpy(frFiles[2],"input/V1.frames");

  strcpy(frFiles[3],"input/L1H1V1-NSNS-Log.frl");
  strcpy(injectionList,"input/L1H1V1-NSNS-Log.txt");

  nDQF=6;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/inspiral.in",           CWB_CAT1, 0., false, false},
                    {"L1" ,"input/inspiral.in",           CWB_CAT2, 0., false, false},
                    {"H1" ,"input/inspiral.in",           CWB_CAT1, 0., false, false},
                    {"H1" ,"input/inspiral.in",           CWB_CAT2, 0., false, false},
                    {"V1" ,"input/inspiral.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/inspiral.in",           CWB_CAT2, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  segLen = 200;
  segMLS = 30;
  segTHR = 0;

  lagSize = 1;
  lagOff = 0;

//  netCC = 0.2; 		//WAVEGRAPH
//  netRHO = 1.0; 	//WAVEGRAPH

  fLow = 32;
  Tgap = 0.5;
  iwindow = 100;	// MDC signal must be inside Tinj +/- gap/2
//  TFgap =  6.;       	// threshold on the time-frequency separation between two pixels
  TFgap =  16.;        	// threshold on the time-frequency separation between two pixels 	//WAVEGRAPH

  simulation = 2;  
  nfactor = 1;
  factors[0]=60.;  

  // regulators
  delta=0.05;
  cfg_gamma = 0.5;
  Acore=1.5;
  subnet = 0.6;

  levelR = 4;        // resampling level
  fHigh = 512.;        // high frequency of the search

  healpix=7;

  cedDump   = true;   // dump ced plots with rho>cedRHO
  cedRHO    = 4.0;

  plugin = TMacro("macro/CWB_Plugin_wavegraph.C");          	// Macro source
  configPlugin = TMacro("macro/CWB_Plugin_NSNS_Config.C");  	// NSN Macro config
//  configPlugin = TMacro("macro/CWB_Plugin_2xNSNS_Config.C"); 	// 2 inj NSNS Macro config
//  configPlugin = TMacro("macro/CWB_Plugin_CBC_Config.C");  	// CBC Macro config
//  configPlugin = TMacro("macro/CWB_Plugin_Burst_Config.C");  	// BurstMacro config

#ifdef TEST_LAGS
  simulation = 0;  
  nfactor = 1;
  factors[0]=12; 	// -> SNR=60   

  // time shift analysis lags
  //lagSize    = 11;     // number of lags (simulation=1)
  lagSize    = 301;     // number of lags (simulation=1)
  lagStep    = 1.;    // [sec] time interval between lags
  lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  lagMax     = 100;   // 0/>0 -  standard/extended lags
//  lagFile    = "config/lag_list.txt";  // lag file list
//  lagMode[0] = 'r';  // w/r  -  write/read lag list
#endif
}
