// main cwb configuration file for the 1G pipeline
// parameters can be redefined/added in user_parameters.C
// for reference contact
// S.Klimenko, University of Florida. klimenko@phys.ufl.edu
// G.Vedovato, University of Padova, vedovato@lnl.infn.it
{
#ifndef CWB_PARAMETER_FILE
#define CWB_PARAMETER_FILE

  #include "xroot.hh"      // defines macro to manage ROOT5 vs ROOT6

  char analysis[8]="1G";   // 1G pipeline
  CheckAnalysis();         // check consistency with the environment CWB_ANALYSIS

  bool online=false;       // true/false -> online/offline

  char ifo[NIFO_MAX][8];
  if(NIFO_MAX>0) strcpy(ifo[0],"L1");  		 // LIGO Livingston                   
  if(NIFO_MAX>1) strcpy(ifo[1],"H1");  		 // LIGO Hanford                   
  if(NIFO_MAX>2) strcpy(ifo[2],"V1");  		 // Virgo                   
  if(NIFO_MAX>3) strcpy(ifo[3],"I1");  		 // LIGO India                   
  if(NIFO_MAX>4) strcpy(ifo[4],"J1");  		 // KAGRA Japan                  
  if(NIFO_MAX>5) strcpy(ifo[5],"G1");  		 // GEO-600 Hannover
  for(int i=6;i<NIFO_MAX;i++) strcpy(ifo[i],""); // ifo[] can be redefined by user 

  int  nIFO = 3;           // size of network starting with first detector ifo[]
  char refIFO[4] = "L1";   // reference IFO
  SEARCH(char) = 'r';      // see description below
  // user define detectors list : is selected if detectorParams[n].name!=""
  // {name, latitude, longitude, elevation, AltX, AzX, AltY, AzY
  detectorParams detParms[NIFO_MAX];
  detectorParams _detParms = {"",0.,0.,0.,0,0.,0.,0.};
  for(int i=0;i<NIFO_MAX;i++) detParms[i] = _detParms;

  // cWB settings

  int    runID = 0;        // run number, set in the production job
  size_t inRate= 16384;    // input data rate
  double bpp   = 0.0001;   // probability for pixel selection
  double Tgap  = 0.05;     // time gap between clusters (sec)
  double Fgap  = 128.;     // frequency gap between clusters (Hz)
  double fLow  = 64.;      // low frequency of the search
  double fHigh = 2048.;    // high frequency of the search
  size_t fResample = 0;    // if>0 the inRate is resampled to fResample
  double Acore = sqrt(2);  // threshold for selection of core pixels
  double Tlpr  = 120.;     // 1G: training time for LPR filter 

  double x2or  = 1.5;      // 1G: 2 OR threshold 
  double netRHO= 3.5;      // threshold on rho
  double netCC = 0.5;      // threshold on network correlation

  // time-frequency transformation settings
  
  int levelR   = 2;        // resampling level : inRate[fResample]/(2^levelR) Hz  
  int levelF   = 6;        // level where second LPR filter is applied
  int levelD   = 8;        // decomposition level
  int l_low    = 3;        // low frequency resolution level (2^l_low Hz)
  int l_high   = 8;        // high frequency resolution level (2^l_high Hz)

  // time shift analysis settings

  // segments
  double segLen     = 600.;  // Segment length [sec]
  double segMLS     = 300.;  // Minimum Segment Length after DQ_CAT1 [sec]
  double segTHR     = 30.;   // Minimum Segment Length after DQ_CAT2 [sec] (to disable put segTHR=0)
  double segEdge    = 8.;    // wavelet boundary offset [sec]
  double segOverlap = 0.;    // overlap between job segments [sec]

  // lags
  size_t lagSize    = 1;     // number of lags (simulation=1)
  double lagStep    = 1.;    // time interval between lags [sec]
  size_t lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  size_t lagMax     = 150;   // 0/>0 -  standard/extended lags
  char*  lagFile    = NULL;  // lag file list
  char   lagMode[]  = "w";   // w/r  -  write/read lag list
  size_t* lagSite   = NULL;  // site index starting with 0
  double shift[NIFO_MAX];    // use for standard shifts
  for(int i=0;i<NIFO_MAX;i++) shift[i] = 0.;

  // multi lags
  int    mlagStep   = 0;     // if mlagStep=0 then 'standard lag mode' else cicle over lags with step mlagStep 

  // super lags
  int    slagSize   = 0;     // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  int    slagMin    = 0;     // if slagMax=0 ->  slagMin must be < slagMax
  int    slagMax    = 0;     // if slagMax=0 
  int    slagOff    = 0;     // first slag id (slagOff=0 - include zero slag )
  size_t* slagSite  = NULL;  // site index starting with 0
  char*  slagFile   = NULL;  // slag file list

  // whitening parameters
  double whiteWindow = 60.;  // time window dT. if = 0 - dT=T, where T is wavearray duration
  double whiteStride = 20.;  // noise sampling interval (window stride)

  // Skymap probability to be saved in the final output root file (saved if !=0 : see nSky)
  int Psave = 0;

  // DC corrections
  double dcCal[NIFO_MAX];
  for(int i=0;i<NIFO_MAX;i++) dcCal[i] = 1.0;
 
  // simulation parameters
  int simulation = 0;        // 1/2/3 for simulation(hrss/snr/tshift), 0 for production
  int nfactor=0;             // number of strain factors
  double factors[100];       // array of strain factors
  for(int i=0;i<100;i++) factors[i]=0.;
  double iwindow = 5.;       // injection time window (Tinj +/- iwindow/2)

  // noise shift data
  double dataShift[NIFO_MAX];
  for(int i=0;i<NIFO_MAX;i++) dataShift[i] = 0.;

  // use this parameter to shift in time the injections (sec)
  // use {0,0,0} to set mdc_shift to 0
  // if {-1,0,0} the shift is automaticaly selected
  // {startMDC, stopMDC}
  mdcshift mdc_shift = {0, 0, 0};   

  // delay filter

  char   filter[1024] = "up2"; // 1G: delay filter suffix: "", or "up1", or "up2"  

  // regulator
 
  double delta = 1.0;          // 1G: [0/1] -> [weak/soft]    
  GAMMA(double) = 0.2;         // 1G: set params in net5, [0/1]->net5=[nIFO/0], 
                               // 1G: if net5>[threshold=(nIFO-1)] weak/soft[according to delta] else hard 

  bool   eDisbalance = true;   // 1G:
 
  // sky settings

  bool   EFEC   =   true;  // Earth Fixed / Selestial coordinates
  size_t mode   =   0;     // sky search mode
  double angle  =   0.4;   // angular resolution
  double Theta1 =   0.;    // start theta
  double Theta2 = 180.;    // end theta
  double Phi1   =   0.;    // start theta
  double Phi2   = 360.;    // end theta
  double mask   = 0.00;    // sky mask fraction
  size_t healpix= 0;       // if not 0 use healpix sky map (number of sky pixels = 12*pow(4,healpix))

  // error regions settings
  
  double precision = 0.001; //  No = nIFO*(K+KZero)+precision*E
  
  // file dump mode

  CWB_JOBF_OPTIONS jobfOptions = CWB_JOBF_SAVE_DISABLE;   // job file options
  CWB_OUTF_OPTIONS outfOptions = CWB_OUTF_SAVE_DISABLE;   // output root file options

  bool dumpHistory = true;    // dump history into output root file
  bool dump        = false;   // dump triggers into ascii file
  bool savemode    = true;    // temporary save clusters on disc
  bool cedDump     = false;   // dump ced plots with rho>cedRHO
  double cedRHO    = 4.0;
  long nSky        = 0;       // if nSky>0 -> # of skymap prob pixels dumped to ascii
                              // if nSky=0 -> (#pixels==1000 || cum prob > 0.99)
                              // if nSky<0 -> nSky=-XYZ... save all pixels with prob < 0.XYZ...

  // directories, file names

  char filter_dir[1024];
  if(gSystem->Getenv("HOME_WAT_FILTERS")==NULL) {
    cout << "Error : environment HOME_WAT_FILTERS is not defined!!!" << endl;exit(1);
  } else {
    strcpy(filter_dir,TString(gSystem->Getenv("HOME_WAT_FILTERS")).Data());
  }

  char injectionList[1024]="";
  char skyMaskFile[1024]="";
  char skyMaskCCFile[1024]="";
  char channelNamesRaw[NIFO_MAX][50];
  char channelNamesMDC[NIFO_MAX][50];
  for(int i=0;i<NIFO_MAX;i++) strcpy(channelNamesRaw[i],"");
  for(int i=0;i<NIFO_MAX;i++) strcpy(channelNamesMDC[i],"");

  // working dir
  char work_dir[1024];
  sprintf(work_dir,"%s",gSystem->WorkingDirectory());

  char config_dir[1024] = "config";
  char input_dir[1024]  = "input";
  char output_dir[1024] = "output";
  char merge_dir[1024]  = "merge";
  char condor_dir[1024] = "condor";
  char report_dir[1024] = "report";
  char macro_dir[1024]  = "macro";
  char log_dir[1024]    = "log";
  char data_dir[1024]   = "data";
  char tmp_dir[1024]    = "tmp";
  char ced_dir[1024]    = "report/ced";
  char pp_dir[1024]     = "report/postprod";
  char dump_dir[1024]   = "report/dump";
  char www_dir[1024];

  UserGroup_t* uinfo = gSystem->GetUserInfo();
  TString uname = uinfo->fUser;
  TString uhome = TString(gSystem->Getenv("HOME")).Data();

  // set public dir
  if(gSystem->Getenv("WWW_PUBLIC_DIR")!=NULL) {
    // read from env
    TString www_public_dir = TString(gSystem->Getenv("WWW_PUBLIC_DIR"));
    www_public_dir.ReplaceAll("X_HOME",uhome.Data());
    www_public_dir.ReplaceAll("X_USER_NAME",uname.Data());
    www_public_dir.ReplaceAll("X_HOST_NAME",gSystem->HostName());
    sprintf(www_dir,"%s",www_public_dir.Data());
  } else {
    cout << "cwb_parameters.C : Error env WWW_PUBLIC_DIR not defined !!!" << endl;
    cout << "www_public_dir can not be defined" << endl;
    exit(1);                   
  }
  cout << "www_dir    : " << www_dir << endl;

  // data label
  char data_label[1024];
  sprintf(data_label,gSystem->BaseName(work_dir));
  cout << "data_label : " << data_label << endl;

  // condor declarations
  char condor_log[1024];
  if(gSystem->Getenv("CONDOR_LOG_DIR")!=NULL) {
    // get host name & strip string after the first '.'
    // strip is necessary to get the correct condor path
    // condor log dir is in the local header nodes
    // used in ATLAS cluster to fix issue with HSM
    TString host_name = gSystem->HostName();
    if(host_name.First(".")>=0) host_name.Resize(host_name.First("."));
    // read from env
    TString condor_log_dir = TString(gSystem->Getenv("CONDOR_LOG_DIR"));
    condor_log_dir.ReplaceAll("X_HOME",uhome.Data());
    condor_log_dir.ReplaceAll("X_USER_NAME",uname.Data());
    condor_log_dir.ReplaceAll("X_HOST_NAME",host_name);
    sprintf(condor_log,"%s",condor_log_dir.Data());
  } else {
    cout << "cwb_parameters.C : Error env CONDOR_LOG_DIR not defined !!!" << endl;
    cout << "condor_log can not be defined" << endl;
    exit(1);                   
  }
  cout << "condor_log : " << condor_log << endl;

  // Define a Unique Tag for Condor Jobs
  // See the following link:
  // https://ldas-gridmon.ligo.caltech.edu/accounting/condor_groups/determine_condor_account_group.html
  // Ex : ligo.dev.o1.burst.allsky.cwboffline, ligo.prod.o1.burst.allsky.cwboffline
  char condor_tag[1024] = "";

  // frame files list : [0:nIFO-1]/[nIFO:2*nIFO-1] contains strain/mdc file names
  // If all mdc channels are in a single frame file -> mdc must be declared in the nIFO position
  char frFiles[2*NIFO_MAX][1024];
  for(int i=0;i<2*NIFO_MAX;i++) strcpy(frFiles[i],"");
  // frame reading retry time (sec) : 0 -> disable
  // retry time = frRetryTime*(num of trials) : max trials = 3
  int  frRetryTime=60;	

  // dq file list 
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  int nDQF = 0;
  dqfile DQF[20];

  // read and dump data on local disk (nodedir)
  char nodedir[1024] = "";
  if(gSystem->Getenv("NODE_DATA_DIR")!=NULL) {
    // read from env
    TString node_data_dir = TString(gSystem->Getenv("NODE_DATA_DIR"));
    node_data_dir.ReplaceAll("X_HOME",uhome.Data());
    node_data_dir.ReplaceAll("X_USER_NAME",uname.Data());
    node_data_dir.ReplaceAll("X_HOST_NAME",gSystem->HostName());
    sprintf(nodedir,"%s",node_data_dir.Data());
  } else {
    cout << "cwb_parameters.C : Error env NODE_DATA_DIR not defined !!!" << endl;
    cout << "nodir can not be defined" << endl;
    exit(1);                   
  }
  cout << "nodename   : " << gSystem->HostName() << endl;
  cout << "nodedir    : " << nodedir << endl;

  // Plugin

  TMacro plugin;                // Macro source
  TMacro configPlugin;  	// Macro config
  plugin.SetName("");
  configPlugin.SetName("");
  bool dataPlugin = false;      // if dataPlugin=true disable read data from frames
  bool mdcPlugin = false;       // if dataPlugin=true disable read mdc from frames

  // statistics:
  // L  - likelihood
  // c  - network correlation coefficient
  // A  - energy disbalance asymmetry
  // P  - penalty factor based on correlation coefficients <x,s>/sqrt(<x,x>*<s,s>)
  // E  - total energy in the data streams

  // 1G search modes
  // 'c' - un-modeled search, fast S5 cWB version, requires constraint settings
  // 'h' - un-modeled search, S5 cWB version, requires constraint settings
  // 'B' - un-modeled search, max(P*L*c/E) 
  // 'b' - un-modeled search, max(P*L*c*A/E) 
  // 'I' - elliptical polarisation, max(P*L*c/E) 
  // 'S' - linear polarisation, max(P*L*c/E) 
  // 'G' - circular polarisation, max(P*L*c/E) 
  // 'i' - elliptical polarisation, max(P*L*c*A/E) 
  // 's' - linear polarisation, max(P*L*c*A/E) 
  // 'g' - circular polarisation, max(P*L*c*A/E) 

  // ----------------------------------------------
  // obsolete parameters
  // ----------------------------------------------

  double &gap = iwindow;        // alias of iwindow
  char   wdmXTalk[1024]="";	// 2G parameter
  double TFgap = 0.;		// 2G parameter
  size_t upTDF = 0;		// 2G parameter
  size_t TDSize = 0;		// 2G parameter
  int    BATCH = 0; 		// 2G parameter
  int    LOUD = 0; 		// 2G parameter
  bool   dcPlugin   = false; 	// 2G parameter
  bool   outPlugin  = false; 	// 2G parameter
  bool   optim = true;          // 2G parameter
  double subnet = 0.;           // 2G parameter
  bool   dataPlugin = false;    // 2G parameter
  bool   mdcPlugin  = false;    // 2G parameter
  bool   dcPlugin   = false;    // 2G parameter
  bool   cohPlugin  = false;    // 2G parameter
  bool   scPlugin   = false;    // 2G parameter
  bool   outPlugin  = false;    // 2G parameter

#endif   
}
