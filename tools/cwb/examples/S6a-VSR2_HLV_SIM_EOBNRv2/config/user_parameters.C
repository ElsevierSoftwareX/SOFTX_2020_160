{


  // NET setting
  nIFO = 3;           // number of detectors from the array ifo[]
  cfg_search = 'I';       // see description in parameters.C


  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");

  strcpy(refIFO,"L1"); // reference IFO


  // cWB settings
  bpp   = 0.0001;      // probability for pixel selection
  fLow  = 32.;         // low frequency of the search
  fHigh = 512.;        // high frequency of the search
  netRHO= 3.0;         // threshold on rho
  netCC = 0.5;         // threshold on network correlation
  int levelR  = 4;


  // lags
  lagSize    = 1;     // number of lags (simulation=1)
  lagStep    = 1.;    // time interval between lags [sec]
  lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  lagMax     = 300;   // 0/>0 -  standard/extended lags
  lagFile    = NULL;  // lag file list


  // super lags
  slagSize   = 0;    // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  slagMin    = 0;    // if slagMax=0 ->  slagMin must be < slagMax
  slagMax    = 0;    // if slagMax=0 -> slagMax=slagSegs


  // simulation parameters
  int simulation = 1;      // 1 for simulation, 0 for production
  int   nfactor = 16;      // number of strain factors
  double FACTORS[] = {0.1317,0.1975,0.2963,0.4444,0.6667,1.0000,1.5000,2.2500,3.3750,5.0625,7.5937,11.3906,17.0859,25.6289,38.4434,57.6650}; // array of strain factors for S6 IMBH search
  for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];


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
  strcpy(frFiles[0],"input/S6A_L1.frames");
  strcpy(frFiles[1],"input/S6A_H1.frames");
  strcpy(frFiles[2],"input/VSR2_V1.frames"); 
 

  // dq file list 
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=12;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/S6A_OFFLINE_L1SCIENCE.txt",         CWB_CAT0, 0., false,  false},
		    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true,   false},
                    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true,   false},
                    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true,   false},
                    {"H1" ,"input/S6A_OFFLINE_H1SCIENCE.txt",         CWB_CAT0, 0., false,  false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true,   false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true,   false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true,   false},
                    {"V1" ,"input/S6A_OFFLINE_V1SCIENCE.txt",         CWB_CAT0, 0., false,  false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true,   false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true,   false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true,   false}
  };


  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];


  plugin = TMacro("plugins/CWB_Plugin_TestClassCBC.C");        // Macro source
  configPlugin = TMacro("plugins/CWB_Plugin_TestClassCBC_Config.C");  // Macro config


}
