{
  strcpy(analysis,"2G");

  nIFO = 3;
  cfg_search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");
  strcpy(refIFO,"L1");

  lagSize = 1;
  lagOff = 0;
  segTHR = 0;

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
  dqfile dqf[6]={
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

  fLow = 32;
  Tgap = 0.5;
  gap = 100;		// MDC signal must be inside Tinj +/- gap/2

  simulation = 2;  
  nfactor = 1;
  factors[0]=120.; 	// network snr 

  netRHO=5;

  // regulators
  delta=0.05;
  cfg_gamma=0.5;
  Acore=1.5;
  subnet = 0.6;

  levelR = 4;        // resampling level
  fHigh = 512.;        // high frequency of the search

  healpix=7;

  plugin = TMacro("plugins/CWB_Plugin_MDC_OTF.C");               // Macro source
  configPlugin = TMacro("plugins/CWB_Plugin_MDC_OTF_Config_NSNS.C");  // Macro config
}
