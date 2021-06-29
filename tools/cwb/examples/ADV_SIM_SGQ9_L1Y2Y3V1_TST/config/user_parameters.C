{
  strcpy(analysis,"1G"); 
  //strcpy(analysis,"2G"); 

  nIFO = 4;
  cfg_search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"");    // if empty then select user detectors
  strcpy(ifo[2],"");    // if empty then select user detectors
  strcpy(ifo[3],"V1");
  strcpy(refIFO,"L1");

  detectorParams detPARAMS[4] = 
                {
                 {"L1", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, (    -197.716)},

                 {"Y2", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, ( +45-125.999)},  // H2 LCI
                 {"Y3", 46.4551, -119.408, 0.0, 0, ( +45-125.999), 0, (    -125.999)},  // H3 LCI

                 {"V1", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},
                };
  for(int i=0;i<4;i++) detParms[i] = detPARAMS[i];


  lagSize = 1;
  lagOff = 0;
  segTHR = 0;

  strcpy(channelNamesRaw[0],"L1:SIM-STRAIN");
  strcpy(channelNamesRaw[1],"Y2:SIM-STRAIN");
  strcpy(channelNamesRaw[2],"Y3:SIM-STRAIN");
  strcpy(channelNamesRaw[3],"V1:SIM-STRAIN");

  strcpy(channelNamesMDC[0],"L1:GW-H");
  strcpy(channelNamesMDC[1],"Y2:GW-H");
  strcpy(channelNamesMDC[2],"Y3:GW-H");
  strcpy(channelNamesMDC[3],"V1:GW-H");

  strcpy(frFiles[0],"input/L1.frames");
  strcpy(frFiles[1],"input/Y2.frames");
  strcpy(frFiles[2],"input/Y3.frames");
  strcpy(frFiles[3],"input/V1.frames");
  strcpy(frFiles[4],"input/MDC.frames");
  strcpy(injectionList,"frames/L1Y2Y3V1-SGQ9-9311/L1Y2Y3V1-SGQ9-931158100-600-Log.txt");

  nDQF=8;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"L1" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                    {"Y2" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"Y2" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                    {"Y3" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"Y3" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                    {"V1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  cedRHO = 1.000000;
  fLow=32;
  segLen = 60;
  segMLS = 30;

  nfactor = 1;
  simulation = 2;  // snr mode
  factors[0]=30;

//  www_dir is the directory where the symbolic links to report/* are created  
//  sprintf(www_dir,"%s/Sites/waveburst/reports",TString(gSystem->Getenv("HOME")).Data());

//  plugin = TMacro("macro/CWB_Plugin_TestClassMDC.C");        // Macro source
//  configPlugin = TMacro("macro/CWB_Plugin_TestClassMDC_Config.C");  // Macro config
}
