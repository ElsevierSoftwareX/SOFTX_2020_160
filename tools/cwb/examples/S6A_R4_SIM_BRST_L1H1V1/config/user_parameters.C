{
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
  slagSize   = 1;     // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  slagMin    = 0;     // if slagMax=0 -> slagMin must be < slagMax
  slagMax    = 1;     // if slagMax=0 -> slagMax=slagSegs

  // simulation parameters
  simulation = 1;     // 1 for simulation, 0 for production
  nfactor = 12;       // number of strain factors
  double FACTORS[] ={0.0375,
                     0.075,
                     0.15,
                     0.3,
                     0.6,
                     1.2,
                     2.4,
                     4.8,
                     9.6,
                     19.2,
                     38.4,
                     76.8}; // array of strain factors

  for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];

  // regulator
  delta = 1.0;       
  cfg_gamma = 0.2;   

  // file dump mode
  dump        = true;    // dump triggers into ascii file
  cedDump     = false;   // dump ced plots with rho>cedRHO
  cedRHO      = 4.0;
 
  // directories, file names
  strcpy(injectionList,"input/BurstMDC-BRST_S6-Log.txt");

  strcpy(channelNamesRaw[0],"L1:LDAS-STRAIN");
  strcpy(channelNamesRaw[1],"H1:LDAS-STRAIN");
  strcpy(channelNamesRaw[2],"V1:h_16384Hz");

  strcpy(channelNamesMDC[0],"L1:GW-H");
  strcpy(channelNamesMDC[1],"H1:GW-H");
  strcpy(channelNamesMDC[2],"V1:GW-H");

  // frame files list (MDC must be declared in the nIFO+1 position)

  strcpy(frFiles[0],"input/S6A_R3_L1_LDAS_C02_L2.frames");
  strcpy(frFiles[1],"input/S6A_R3_H1_LDAS_C02_L2.frames");
  strcpy(frFiles[2],"input/S6A_R1_V1_HrecV3.frames");
  strcpy(frFiles[nIFO],"input/S6_MDC_BRST.frames");

  // dq file list 
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=12;
  dqfile dqf[nDQF]={
		    {"L1" ,"input/S6A_OFFLINE_L1SCIENCE.txt",         CWB_CAT0, 0., false, false},
                    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true, false},
                    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"H1" ,"input/S6A_OFFLINE_H1SCIENCE.txt",         CWB_CAT0, 0., false, false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true, false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"V1" ,"input/S6A_OFFLINE_V1SCIENCE.txt",         CWB_CAT0, 0., false, false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true, false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true, false} };

  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  // uncomment plugin to apply phase mis-calibration to data
  //plugin = TMacro("plugins/CWB_Plugin_PhaseMisCal.C");        // Macro source

}
