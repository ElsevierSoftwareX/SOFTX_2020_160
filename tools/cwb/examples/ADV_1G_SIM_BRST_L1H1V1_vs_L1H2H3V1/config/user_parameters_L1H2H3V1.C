{
  strcpy(analysis,"1G");

  nIFO = 4;
  search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"");
  strcpy(ifo[2],"");
  strcpy(ifo[3],"V1");
  strcpy(refIFO,"L1");

  detectorParams detPARAMS[nIFO] = {
                          {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, (    -197.716 )},  // L1
                          {"H2", 46.4551, -119.408,  0.0, 0, ( +90-125.998),    0, (+45-125.998 )},   // H1
                          {"H3", 46.4551, -119.408,  0.0, 0, ( +45-125.998),    0, (   -125.998 )},   // H1
                          {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675),    0, (    -70.5675)},   // V1
                         };
  for(int i=0;i<nIFO;i++) detParms[i] = detPARAMS[i];

  lagSize = 1;
  lagOff = 0;
  segTHR = 0;

  strcpy(channelNamesRaw[0],"L1:SIM-STRAIN");
  strcpy(channelNamesRaw[1],"H2:SIM-STRAIN");
  strcpy(channelNamesRaw[2],"H3:SIM-STRAIN");
  strcpy(channelNamesRaw[3],"V1:SIM-STRAIN");

  strcpy(frFiles[0],"input/L1.frames");
  strcpy(frFiles[1],"input/H2.frames");
  strcpy(frFiles[2],"input/H3.frames");
  strcpy(frFiles[3],"input/V1.frames");

  nDQF=4;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"H2" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"H3" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];
  cedRHO = 3.000000;
  fLow=64;

  segLen = 60;
  segMLS = 30;

  gap=5;

  nfactor = 1;
//  simulation = 2;  // snr mode
//  factors[0]=10*sqrt(3);  
  simulation = 1;  // hrss mode
  factors[0]=1;  

  healpix = 7; 

  plugin = TMacro("macro/CWB_Plugin_BRST2.C");               // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_BRST2_Config.C");  // Macro config
}
