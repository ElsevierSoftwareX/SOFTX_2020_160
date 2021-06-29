{
  strcpy(analysis,"2G");

  nIFO = 2;
  search = 'r';
  optim=false;

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(refIFO,"L1");

  //lags
  lagSize = 590;
  lagStep = 1.;
  lagOff = 0;
  lagMax = 0;

  //superlags
  slagSize   = 1;     // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  slagMin    = 0;      // if slagMax=0 ->  slagMin must be < slagMax
  slagMax    = 1;    // if slagMax=0 -> slagMax=slagSegs
  slagOff    = 0;

  //jobs
  segLen = 600;
  segMLS = 600;
  segTHR = 0;
  segEdge = 10;

  //frequency
  fLow  = 16.;       // low frequency of the search
  fHigh = 1024.;      // high frequency of the search

  levelR = 3;
  l_low    = 4;       // low frequency resolution level
  l_high   = 10;      // high frequency resolution level

  strcpy(wdmXTalk,"wdmXTalk/OverlapCatalog16-1024.bin");

  healpix=7;

  bpp   = 0.002;
  subnet = 0.6;
  netRHO = 5.0;
  netCC = 0.5;
  Acore = 1.7;
  Tgap  = 0.0;
  Fgap  = 0.0;
  delta = 0.5;
  gamma = -0.5;
  LOUD = 300;

  //simulation
  nfactor = 1;
  simulation = 0;
  strcpy(channelNamesRaw[0],"L1:DCS-CALIB_STRAIN_C01");
  strcpy(channelNamesRaw[1],"H1:DCS-CALIB_STRAIN_C01");
  strcpy(frFiles[0],"input/L1_C01_ATLAS.frames");
  strcpy(frFiles[1],"input/H1_C01_ATLAS.frames");

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=6;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/period_first_100jobs.txt",   CWB_CAT0, 0., false, false},
                    {"H1" ,"input/period_first_100jobs.txt",   CWB_CAT0, 0., false, false},
                    {"L1" ,"input/L1_CAT0_C01.txt",            CWB_CAT0, 0., false, false},
                    {"H1" ,"input/H1_CAT0_C01.txt",            CWB_CAT0, 0., false, false},
                    {"L1" ,"input/L1_CAT1_C01.txt",            CWB_CAT1, 0., true,  false},
                    {"H1" ,"input/H1_CAT1_C01.txt",            CWB_CAT1, 0., true,  false}
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  plugin = TMacro("macro/CWB_Plugin_QLWveto_Gating.C");        // Macro source
  plugin.SetTitle("macro/CWB_Plugin_QLWveto_Gating_C.so");

  strcpy(condor_tag,"ligo.prod.o1.burst.allsky.cwboffline");

  jobfOptions |= CWB_JOBF_SAVE_TRGFILE;

}
