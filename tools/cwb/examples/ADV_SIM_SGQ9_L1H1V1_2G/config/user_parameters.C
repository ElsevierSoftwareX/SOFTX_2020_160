{

  strcpy(analysis,"2G");

  nIFO = 3;
  cfg_search = 'r';
  optim=false;

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");
  strcpy(refIFO,"L1");

  //lags
  lagSize = 1;
  lagStep = 1.;
  lagOff  = 0;
  lagMax  = 0;

  //superlags
  slagSize = 0;        // number of super lags (simulation=1) - if slagSize=0 -> Standard Segments
  slagMin  = 0;
  slagMax  = 0;
  slagOff  = 0;


  //jobs
  segLen  = 60;
  segMLS  = 60;
  segTHR  = 0;
  segEdge = 10;

  //frequency
  fLow  = 16.;        // low frequency of the search
  fHigh = 2048.;      // high frequency of the search

  levelR   = 2;
  l_low    = 4;       // low frequency resolution level
  l_high   = 10;      // high frequency resolution level

  strcpy(wdmXTalk,"wdmXTalk/OverlapCatalog16-1024.bin");

  healpix=7;

  bpp    = 0.001;
  subnet = 0.5;
  subcut = 0.0;
  netRHO = 5.0;
//  netCC  = 0.5;
netCC  = 0.0;		// the threshold has been lowered to permit the production of complete skymaps in the CED
  Acore  = 1.7;
  Tgap   = 0.2;
  Fgap   = 128.0;
  delta  = 0.5;
  cfg_amma  = -1.0;
  LOUD   = 300;

  pattern = 10;

  //precision=GetPrecision(100,5);

  // simulation 
  nfactor = 1;
  simulation = 2;   
  factors[0]=30;   // Network SNR

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

  nDQF=6;
  dqfile dqf[6]={
                    {"L1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"L1" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                    {"H1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"H1" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                    {"V1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/burst.in",           CWB_CAT2, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  strcpy(condor_tag,"ligo.dev.o2.burst.allsky.cwboffline");
}
