{
  nIFO = 3;
  cfg_search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");
  strcpy(refIFO,"L1");
//  simulation = 1;
  lagSize = 1;
  lagOff = 0;
  segTHR = 0;
  strcpy(channelNamesRaw[0],"L1:LDAS-STRAIN");
  strcpy(channelNamesRaw[1],"H1:LDAS-STRAIN");
  strcpy(channelNamesRaw[2],"V1:h_16384Hz");
  strcpy(channelNamesMDC[0],"L1:GW-H");
  strcpy(channelNamesMDC[1],"H1:GW-H");
  strcpy(channelNamesMDC[2],"V1:GW-H");
  strcpy(frFiles[0],"input/L1.frames");
  strcpy(frFiles[1],"input/H1.frames");
  strcpy(frFiles[2],"input/V1.frames");
  strcpy(frFiles[3],"input/MDC.frames");
  strcpy(injectionList,"/home/waveburst/soft/eced/eced-1.3/tst/eced_ga1d0/BurstMDC-BRST50MPC_S6-Log.txt");
  nfactor = 1;
//  factors[0] = 1.000000;
  nDQF=6;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"L1" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                    {"H1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"H1" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                    {"V1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/burst.in",           CWB_CAT2, 0., false, false},
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
  factors[0]=30;

  sprintf(www_dir,"%s/public_html/reports",TString(gSystem->Getenv("HOME")).Data());

  plugin = TMacro("plugins/CWB_Plugin_TestClassMDC.C");        // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_TestClassMDC_Config.C");  // Macro config
}
