{
  strcpy(analysis,"2G");

  // NET setting
  nIFO = 3;           // number of detectors from the array ifo[]

  cfg_search = 'i';

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");

  strcpy(refIFO,"L1"); // reference IFO

  // cWB settings
  bpp   = 0.0001;      // probability for pixel selection
  fLow  = 64.;       // low frequency of the search
  fHigh = 2048.;      // high frequency of the search

  netCC = 0.5;        // threshold on network correlation

  // lags
  lagSize    = 1;     // number of lags (simulation=1)
  lagStep    = 1.;    // time interval between lags [sec]
  lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  lagMax     = 300;   // 0/>0 -  standard/extended lags
  lagFile    = NULL;  // lag file list

  // super lags
  slagSize   = 0;    // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  slagMin    = 0;     // if slagMax=0 ->  slagMin must be < slagMax
  slagMax    = 0;    // if slagMax=0 -> slagMax=slagSegs

  // simulation parameters
  simulation = 2;     // 2 for simulation snr const
  nfactor = 1;        // number of strain factors
  factors[0]=10.*sqrt(nIFO);

  healpix = 7;

  bpp   = 0.001;      // probability for black pixel selection (netpixel)
  subnet = 0.6;       // [<0.7] subnetwork threshold (supercluster)
  netRHO = 5.0;       // [>4.0] coherent network SNR (supercluster, likelihood)
  netCC = 0.5;        // network correlation (supercluster, likelihood)
  Acore = 1.5;        // threshold of core pixels (supercluster, likelihood)
  Tgap  = 3.0;        // defragmentation time gap between clusters (sec)
  Fgap  = 130.;       // defragmentation frequency gap between clusters (Hz)
  delta = 0.05;       // [0-1] regularize sky locations with low |fx|^2,
  cfg_gamma = 0.5;        // [0-1] defines threshold on coherent energy: cfg_gamma = 0 - Ec>0
  LOUD = 300;         // number of pixel per cluster to load TD amplitudes


  segLen = 60;
  segMLS = 60;
  segEdge = 8;
  segOverlap = 15;

  // file dump mode
  dump        = true;    // dump triggers into ascii file
  cedDump     = false;   // dump ced plots with rho>cedRHO
  cedRHO      = 3.0;
 

  strcpy(channelNamesRaw[0],"L1:SIM-STRAIN");
  strcpy(channelNamesRaw[1],"H1:SIM-STRAIN");
  strcpy(channelNamesRaw[2],"V1:SIM-STRAIN");

  strcpy(channelNamesMDC[0],"L1:GW-H");
  strcpy(channelNamesMDC[1],"H1:GW-H");
  strcpy(channelNamesMDC[2],"V1:GW-H");

  strcpy(frFiles[0],"input/L1.frames");
  strcpy(frFiles[1],"input/H1.frames");
  strcpy(frFiles[2],"input/V1.frames");

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


  plugin = TMacro("macro/CWB_Plugin_MDC_OTF.C");                   // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_EOBNRv2_Config.C");      // Macro config
}
