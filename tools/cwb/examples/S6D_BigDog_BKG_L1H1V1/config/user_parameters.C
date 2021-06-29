{
  strcpy(analysis,"1G");
//  healpix=7;

  nIFO = 3;
  cfg_search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");
  strcpy(refIFO,"L1");
  simulation = 0;
  lagSize = 1;
  lagOff = 0;
  segLen = 60;
  segMLS = 30;
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
  nDQF=6;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/ws_cat1.in",           CWB_CAT1, 0., false, false},
                    {"L1" ,"input/ws_cat2.in",           CWB_CAT2, 0., false, false},
                    {"H1" ,"input/ws_cat1.in",           CWB_CAT1, 0., false, false},
                    {"H1" ,"input/ws_cat2.in",           CWB_CAT2, 0., false, false},
                    {"V1" ,"input/ws_cat1.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/ws_cat2.in",           CWB_CAT2, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];
  sprintf(refIFO,"L1");
  bpp = 0.000100;
  cedRHO = 3.000000;
  nSky = 2000;
  Psize = 4000;
  healpix=7;
}
