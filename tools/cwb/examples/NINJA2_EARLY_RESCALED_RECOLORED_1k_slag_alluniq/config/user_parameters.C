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
  fHigh = 1024.;        // high frequency of the search
  netRHO= 3.5;         // threshold on rho
  netCC = 0.5;         // threshold on network correlation

 // wavelet transformation settings
    levelR   = 3;        // resampling level
    levelF   = 6;        // level where second LPR filter is applied
    levelD   = 8;        // decomposition level
    l_low    = 3;        // low frequency resolution level  
    l_high   = 8;        // high frequency resolution level           

  // lags
  lagSize    = 201;     // number of lags (simulation=1)
  lagStep    = 2.985;    // time interval between lags [sec]
  lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  lagMax     = 600;   // 0/>0 -  standard/extended lags
  lagFile    = "lag_uniq200.txt";  // lag file list
  lagMode[0]  = 'r';

  // super lags
  slagFile = new char[1024];
  strcpy(slagFile,"slag_uniq200.txt"); 
  slagSize   = 10;    // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
//  slagMin    = 5;    // if slagMax=0 ->  slagMin must be < slagMax
//  slagMax    = 15;    // if slagMax=0 -> slagMax=slagSegs
  slagOff    = 0;  
 
 // simulation parameters
  int simulation = 0;      // 1 for simulation, 0 for production
 // int   nfactor = 16;      // number of strain factors
 // double FACTORS[] = {0.026,0.039,0.0585,0.0878,0.1317,0.1975,0.2963,0.4444,0.6667,1.0000,1.5000,2.2500,3.3750,5.0625,7.5937,11.3906}; // array of strain factors for ADV IMBH search


 // for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];

  // regulator
  delta = 1.0;       
  cfg_gamma = 0.2;   

  // file dump mode
  dump        = true;    // dump triggers into ascii file
  cedDump     = false;   // dump ced plots with rho>cedRHO
  cedRHO      = 4.0;

  strcpy(channelNamesRaw[0],"L1:LDAS-STRAIN");
  strcpy(channelNamesRaw[1],"H1:LDAS-STRAIN");
  strcpy(channelNamesRaw[2],"V1:h_16384Hz");

 strcpy(frFiles[0],"input/L1_NINJA2_G1000176_EARLY_RECOLORED_ATLAS.list");
 strcpy(frFiles[1],"input/H1_NINJA2_G1000176_EARLY_RECOLORED_ATLAS.list");
 strcpy(frFiles[2],"input/V1_NINJA2_G1000176_EARLY_RESCALED_RECOLORED_ATLAS.list");


  // dq file list 
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=12;
  dqfile dqf[nDQF]={
		    {"L1" ,"input/L1_NINJA2_G1000176_EARLY_RECOLORED.txt", CWB_CAT0, 0., false, false},
                    {"L1" ,"input/S6D_OFFLINE_L1_DQCAT1SEGMENTS.txt", CWB_CAT1, -66858896, true, false},
                    {"L1" ,"input/S6D_OFFLINE_L1_DQCAT2SEGMENTS.txt", CWB_CAT2, -66858896, true, false},
                    {"L1" ,"input/S6D_OFFLINE_L1_DQCAT4SEGMENTS.txt", CWB_CAT1, -66858896, true, false},
                    {"H1" ,"input/H1_NINJA2_G1000176_EARLY_RECOLORED.txt", CWB_CAT0, 0., false, false},
                    {"H1" ,"input/S6D_OFFLINE_H1_DQCAT1SEGMENTS.txt", CWB_CAT1, -61545543, true, false}
                    {"H1" ,"input/S6D_OFFLINE_H1_DQCAT2SEGMENTS.txt", CWB_CAT2, -61545543, true, false}
                    {"H1" ,"input/S6D_OFFLINE_H1_DQCAT4SEGMENTS.txt", CWB_CAT1, -61545543, true, false}
                    {"V1" ,"input/V1_NINJA2_G1000176_EARLY_RECOLORED.txt", CWB_CAT0, 0., false, false}
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT1SEGMENTS.txt", CWB_CAT1, -31035296., true, false}
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT2SEGMENTS.txt", CWB_CAT2, -31035296., true, false}
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT4SEGMENTS.txt", CWB_CAT1, -31035296., true, false}
};
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

//UserGroup_t* uinfo = gSystem->GetUserInfo();
 // TString uname = uinfo->fUser;
 // sprintf(nodedir,"/data/%s/%s",gSystem->HostName(),uname.Data());

  // condor declarations
  //  sprintf(condor_log,"/usr1/%s",uinfo->fUser.Data());
  
    //  sprintf(www_dir,"%s/public_html/reports",TString(gSystem->Getenv("HOME")).Data());
      
 plugin = TMacro("macro/CWB_Plugin_pixeLHood.C");        // Macro source

// plugin = TMacro("macro/CWB_Plugin_MakeSpectrum.C");        // Macro source
// plugin = TMacro("plugins/CWB_Plugin_TestClassCBC.C");        // Macro source
// configPlugin = TMacro("plugins/CWB_Plugin_TestClassCBC_Config.C");  // Macro config

}
