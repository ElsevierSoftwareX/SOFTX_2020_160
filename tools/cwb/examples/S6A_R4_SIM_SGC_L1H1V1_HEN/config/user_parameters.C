{
  strcpy(analysis,"1G");

  // NET setting
  nIFO = 3;           // number of detectors from the array ifo[]
  cfg_search = 'g';       // see description in parameters.C

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");

  strcpy(refIFO,"L1"); // reference IFO

  // cWB settings

  bpp   = 0.0001;   // probability for pixel selection
  Tgap  = 0.05;     // time gap between clusters (sec)
  Fgap  = 128.;     // frequency gap between clusters (Hz)
  fLow  = 64.;      // low frequency of the search
  fHigh = 2048.;    // high frequency of the search
  Acore = sqrt(2);  // threshold for selection of core pixels
  Tlpr  = 120.;     // training time for LPR filter

  x2or  = 1.5;      // 2 OR threshold
  netRHO= 2.5;      // threshold on rho
  netCC = 0.5;      // threshold on network correlation

  // wavelet transformation settings

  levelR   = 2;        // resampling level
  levelF   = 6;        // level where second LPR filter is applied
  levelD   = 8;        // decomposition level
  l_low    = 3;        // low frequency resolution level
  l_high   = 8;        // high frequency resolution level

  // segments
  segLen     = 1000.-16;  // Segment length [sec]
  segMLS     = 300;  // Minimum Segment Length after DQ_CAT1 [sec]
  segTHR     = 30.;   // Minimum Segment Length after DQ_CAT2 [sec] (to disable put segTHR=0)
  segEdge    = 8.;    // wavelet boundary offset [sec]

  // lags
  lagSize    = 1;     // number of lags (simulation=1)
  lagStep    = 1.;    // time interval between lags [sec]
  lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  lagMax     = 0;     // 0/>0 -  standard/extended lags
  lagFile    = NULL;  // lag file list

  // super lags
  slagSize   = 1;     // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  slagMin    = 0;     // if slagMax=0 -> slagMin must be < slagMax
  slagMax    = 0;     // if slagMax=0 -> slagMax=slagSegs

  // regulator must be disable when skymask is applied !!!

  delta = 0.0;       // -1/0/1/2 - soft/weak/mild/hard regulator
  cfg_gamma = 0.0;

  // simulation parameters
  simulation = 1;     // 1 for simulation, 0 for production
  nfactor = 12;       // number of strain factors
  double FACTORS[] ={0.038,
                     0.06,
                     0.10,
                     0.17,
                     0.28,
                     0.45,
                     0.7,
                     1.3,
                     2.1,
                     3.4,
                     5.6,
                     9.2}; // array of strain factors

  // file dump mode
  dump        = true;    // dump triggers into ascii file
  cedDump     = false;   // dump ced plots with rho>cedRHO
  cedRHO      = 3.0;

  // directories, file names
  strcpy(injectionList,"config/S6A_R4_SIM_SGC_L1H1V1_HEN.inj");
 
  strcpy(channelNamesRaw[0],"L1:LDAS-STRAIN");
  strcpy(channelNamesRaw[1],"H1:LDAS-STRAIN");
  strcpy(channelNamesRaw[2],"V1:h_16384Hz");

  // frame files list (MDC must be declared in the nIFO+1 position)

  strcpy(frFiles[0],"input/S6A_R3_L1_LDAS_C02_L2.frames");
  strcpy(frFiles[1],"input/S6A_R3_H1_LDAS_C02_L2.frames");
  strcpy(frFiles[2],"input/S6A_R1_V1_HrecV3.frames");

  // dq file list 
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=14;
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
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true, false},

		    {"L1" ,"input/HEN_S6A_L1H1V1_OnSource_Segments_Cat1.txt",CWB_CAT1, 0., false, false},
		    {"L1" ,"input/HEN_S6A_L1H1V1_OnSource_Segments_Cat2.txt",CWB_CAT2, 0., false, false} };

  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  plugin = TMacro("macro/CWB_Plugin_HEN_SIM.C");               // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_HEN_SIM_Config.C");  // Macro config
}
