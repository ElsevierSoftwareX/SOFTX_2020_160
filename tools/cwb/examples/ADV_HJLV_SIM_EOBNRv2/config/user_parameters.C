{


  // NET setting
  nIFO = 4;           // number of detectors from the array ifo[]
  cfg_search = 'I';       // see description in parameters.C


  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1"); 
  strcpy(ifo[3],"J1");

  strcpy(refIFO,"L1"); // reference IFO


  // cWB settings
  bpp   = 0.0001;      // probability for pixel selection
  fLow  = 16.;         // low frequency of the search
  fHigh = 512.;        // high frequency of the search
  netRHO= 2.5;         // threshold on rho
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
  int nfactor = 1;      // number of strain factors
  double FACTORS[] = {1.0000}; // array of strain factors for IMBH ADE search
  for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];


  // regulator
  delta = 1.0;       
  cfg_gamma = 0.2;   


  // file dump mode
  dump        = true;    // dump triggers into ascii file
  cedDump     = false;   // dump ced plots with rho>cedRHO
  cedRHO      = 4.0;


  // dq file list 
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=8;
  dqfile dqf[nDQF]={
		    {"L1" ,"input/inspiral.in", CWB_CAT1, 0., false, false},
                    {"L1" ,"input/inspiral.in", CWB_CAT2, 0., false, false},
                    {"H1" ,"input/inspiral.in", CWB_CAT1, 0., false, false},
                    {"H1" ,"input/inspiral.in", CWB_CAT2, 0., false, false},
                    {"V1" ,"input/inspiral.in", CWB_CAT1, 0., false, false},
                    {"V1" ,"input/inspiral.in", CWB_CAT2, 0., false, false},
                    {"J1" ,"input/inspiral.in", CWB_CAT1, 0., false, false},
                    {"J1" ,"input/inspiral.in", CWB_CAT2, 0., false, false}
  };


  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];


  plugin = TMacro("plugins/CWB_Plugin_TestClassCBC.C");        // Macro source
  configPlugin = TMacro("plugins/CWB_Plugin_TestClassCBC_Config.C");  // Macro config


}
