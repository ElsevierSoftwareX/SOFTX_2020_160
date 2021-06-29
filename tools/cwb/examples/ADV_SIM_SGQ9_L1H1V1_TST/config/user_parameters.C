{

  strcpy(analysis,"2G"); 

  nIFO = 3;
  SEARCH() = 'r';
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
  strcpy(frFiles[3],"input/MDC.frames");
  strcpy(injectionList,"frames/L1H1V1-SGQ9-9311/L1H1V1-SGQ9-931158100-600-Log.txt");

  const int NDQF = 6;
  nDQF=NDQF;
  dqfile dqf[NDQF]={
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
  segLen = 60;
  segMLS = 60;

  dump = true;

  nfactor = 1;
  simulation = 2;  // snr mode
  factors[0]=30;
}
