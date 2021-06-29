{
  strcpy(analysis,"2G");

  // NET setting
  nIFO = 3;           // number of detectors from the array ifo[]
  SEARCH() = 'g';     // see description in parameters.C
  optim=false;

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");

  strcpy(refIFO,"L1"); // reference IFO

  // cWB settings
  fLow  = 32.;        // low frequency of the search
  fHigh = 1024.;      // high frequency of the search

  iwindow = 10;

  levelR = 3;

  // job segments
  segLen     = 1400.; // Segment length [sec]
  segMLS     = 700.;  // Minimum Segment Length after DQ_CAT1 [sec]
  segEdge    = 12.;   // wavelet boundary offset [sec]

  // lags
  lagSize    = 1;     // number of lags (simulation=1)
  lagStep    = 1;     // time interval between lags [sec]
  lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  lagMax     = 0;     // 0/>0 -  standard/extended lags
  lagFile    = NULL;  // lag file list

  // super lags
  slagSize   = 0;    // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  slagMin    = 0;     // if slagMax=0 ->  slagMin must be < slagMax
  slagMax    = 0;    // if slagMax=0 -> slagMax=slagSegs

  // simulation parameters
  simulation = 2;     // 1 for simulation, 0 for production
  nfactor = 2;        // number of strain factors
  double FACTORS[] ={20,30}; // array of strain factors
  for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];

  l_low    = 3;       // low frequency resolution level
  l_high   = 10;      // high frequency resolution level

  strcpy(wdmXTalk,"wdmXTalk/OverlapCatalog8-1024.bin");

  healpix = 7;

  bpp   = 0.001;      // probability for black pixel selection (netpixel)
  subnet = 0.6;       // [<0.7] subnetwork threshold (supercluster)
  netRHO = 5.0;       // [>4.0] coherent network SNR (supercluster, likelihood)
  netCC = 0.5;        // network correlation (supercluster, likelihood)
  Acore = 1.7;
  Tgap  = 0.0;        // defragmentation time gap between clusters (sec)
  Fgap  = 0.0;        // defragmentation frequency gap between clusters (Hz)
  delta = 0.5;       // [0-1] regularize sky locations with low |fx|^2,
  GAMMA() = 0.5;      // [0-1] defines threshold on coherent energy: gamma=0 - Ec>0
  LOUD = 300;         // number of pixel per cluster to load TD amplitudes

  // file dump mode
  dump        = true;    // dump triggers into ascii file
  cedDump     = false;   // dump ced plots with rho>cedRHO
  cedRHO      = 3.0;

  // directories, file names
  strcpy(injectionList,"input/BBH.inj");

  strcpy(channelNamesRaw[0],"L1:SIM-STRAIN");
  strcpy(channelNamesRaw[1],"H1:SIM-STRAIN");
  strcpy(channelNamesRaw[2],"V1:SIM-STRAIN");

  strcpy(channelNamesMDC[0],"L1:GW-H");
  strcpy(channelNamesMDC[1],"H1:GW-H");
  strcpy(channelNamesMDC[2],"V1:GW-H");

  // dq file list 
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  const int NDQF=6;
  dqfile dqf[NDQF]={
                    {"L1" ,"input/period.in",           CWB_CAT1, 0., false, false},
                    {"L1" ,"input/period.in",           CWB_CAT2, 0., false, false},
                    {"H1" ,"input/period.in",           CWB_CAT1, 0., false, false},
                    {"H1" ,"input/period.in",           CWB_CAT2, 0., false, false},
                    {"V1" ,"input/period.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/period.in",           CWB_CAT2, 0., false, false},
                   };
  nDQF=NDQF;for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  plugin = TMacro("macro/CWB_Plugin_SimNoise_MDC_OTF.C");        // Macro source
  plugin.SetTitle("macro/CWB_Plugin_SimNoise_MDC_OTF_C.so");
  configPlugin = TMacro("macro/CWB_Plugin_MDC_OTF_Config_BBH.C");  // Macro config

}
