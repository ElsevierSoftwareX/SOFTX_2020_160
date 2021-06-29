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

  // lags
  lagSize    = 201;   // number of lags (simulation=1)
  lagStep    = 1.;    // time interval between lags [sec]
  lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  lagMax     = 300;   // 0/>0 -  standard/extended lags
  lagFile    = NULL;  // lag file list

  // super lags
  slagSize   = 10;    // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  slagMin    = 0;     // if slagMax=0 ->  slagMin must be < slagMax
  slagMax    = 10;    // if slagMax=0 -> slagMax=slagSegs
  slagOff    = 0;     // first slag id (slagOff=0 - include zero slag )

  // simulation parameters
  simulation = 0;     // 1 for simulation, 0 for production
  nfactor = 1;        // number of strain factors
  factors[0]=1.;

  // regulator
  delta = 1.0;       
  cfg_gamma = 0.2;   

  // file dump mode
  dump        = true;    // dump triggers into ascii file
  cedDump     = false;   // dump ced plots with rho>cedRHO
  cedRHO      = 4.0;
 
  // directories, file names
  strcpy(channelNamesRaw[0],"L1:LDAS-STRAIN");
  strcpy(channelNamesRaw[1],"H1:LDAS-STRAIN");
  strcpy(channelNamesRaw[2],"V1:h_16384Hz");

  // frame files list (MDC must be declared in the nIFO+1 position)
  strcpy(frFiles[0],"input/S6A_R3_L1_LDAS_C02_L2.frames");
  strcpy(frFiles[1],"input/S6A_R3_H1_LDAS_C02_L2.frames");
  strcpy(frFiles[2],"input/S6A_R1_V1_HrecV3.frames");

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
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true, false},
                   };

  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];
}
