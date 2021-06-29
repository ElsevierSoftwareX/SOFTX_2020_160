{
  nIFO = 4;
  cfg_search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"");    // if empty then select user detectors
  strcpy(ifo[2],"");    // if empty then select user detectors
  strcpy(ifo[3],"V1");
  strcpy(refIFO,"L1");

  detectorParams detParms[4] = {
                          {"L1", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, (    -197.716)},

                          {"Y2", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, ( +45-125.999)},  // H2 LCI
                          {"Y3", 46.4551, -119.408, 0.0, 0, ( +45-125.999), 0, (    -125.999)},  // H3 LCI

                          {"V1", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},  
                         };

//  simulation = 1;
  lagSize = 1;
  lagOff = 0;
  segTHR = 0;
  strcpy(channelNamesRaw[0],"L1:LDAS-STRAIN");
  strcpy(channelNamesRaw[1],"Y2:LDAS-STRAIN");
  strcpy(channelNamesRaw[2],"Y3:LDAS-STRAIN");
  strcpy(channelNamesRaw[3],"V1:h_16384Hz");
  strcpy(channelNamesMDC[0],"L1:GW-H");
  strcpy(channelNamesMDC[1],"Y2:GW-H");
  strcpy(channelNamesMDC[2],"Y3:GW-H");
  strcpy(channelNamesMDC[3],"V1:GW-H");
  strcpy(frFiles[0],"input/L1.frames");
  strcpy(frFiles[1],"input/Y2.frames");
  strcpy(frFiles[2],"input/Y3.frames");
  strcpy(frFiles[3],"input/V1.frames");
  strcpy(frFiles[4],"input/MDC.frames");
  strcpy(injectionList,"/home/waveburst/soft/eced/eced-1.3/tst/eced_ga1d0/BurstMDC-BRST50MPC_S6-Log.txt");
  nfactor = 1;
//  factors[0] = 1.000000;
  nDQF=8;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/inspiral.in",           CWB_CAT1, 0., false, false},
                    {"L1" ,"input/inspiral.in",           CWB_CAT2, 0., false, false},
                    {"Y2" ,"input/inspiral.in",           CWB_CAT1, 0., false, false},
                    {"Y2" ,"input/inspiral.in",           CWB_CAT2, 0., false, false},
                    {"Y3" ,"input/inspiral.in",           CWB_CAT1, 0., false, false},
                    {"Y3" ,"input/inspiral.in",           CWB_CAT2, 0., false, false},
                    {"V1" ,"input/inspiral.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/inspiral.in",           CWB_CAT2, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];
  cedRHO = 1.000000;
  fLow=32;
//  cedDump=true;
//  factors[0]=0.064;
//  segLen = 60;
//  segMLS = 30;

//  segEdge=4;
  simulation = 2;  // snr mode
  factors[0]=35;  // snr/detector = 17

  plugin = TMacro("plugins/CWB_Plugin_TestClassCBC.C");        // Macro source
  configPlugin = TMacro("plugins/CWB_Plugin_TestClassCBC_Config.C");  // Macro config
  
}
