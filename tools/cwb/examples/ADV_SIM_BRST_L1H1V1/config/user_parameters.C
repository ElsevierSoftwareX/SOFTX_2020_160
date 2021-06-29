{
  strcpy(analysis,"1G");

  // NET setting
  nIFO = 3;           // number of detectors from the array ifo[]
  cfg_search = 'r';       // see description in parameters.C

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");

  strcpy(refIFO,"L1"); // reference IFO

  // cWB settings
  bpp   = 0.0001;      // probability for pixel selection
  fLow  = 64.;       // low frequency of the search
  fHigh = 2048.;      // high frequency of the search

  netRHO= 3.5;        // threshold on rho
  netCC = 0.5;        // threshold on network correlation

  // wavelet transformation settings

  levelR   = 2;        // resampling level
  levelF   = 6;        // level where second LPR filter is applied
  levelD   = 8;        // decomposition level
  l_low    = 3;        // low frequency resolution level
  l_high   = 8;        // high frequency resolution level

  // lags
  lagSize    = 1;     // number of lags (simulation=1)
  lagStep    = 1.;    // time interval between lags [sec]
  lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  lagMax     = 0;     // 0/>0 -  standard/extended lags
  lagFile    = NULL;  // lag file list

  // super lags
  slagSize   = 0;     // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  slagMin    = 0;     // if slagMax=0 -> slagMin must be < slagMax
  slagMax    = 0;     // if slagMax=0 -> slagMax=slagSegs

  // simulation parameters
  simulation = 2;     // factor are snr
  nfactor = 1;       // number of strain factors
  double snr = 10.*sqrt(3.);
  double FACTORS[] ={snr};

  for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];

  // regulator
  delta = 1.0;       
  cfg_gamma = 0.2;  

  healpix = 7;

  // file dump mode
  dump        = true;    // dump triggers into ascii file
  cedDump     = false;   // dump ced plots with rho>cedRHO
  cedRHO      = 4.0;

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
 
  // directories, file names
  strcpy(injectionList,"input/TST_SIM_BRST.inj");

  // dq file list 
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
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

  plugin = TMacro("macro/CWB_Plugin_BRST.C");        // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_BRST_Config.C");  // Macro config
}
