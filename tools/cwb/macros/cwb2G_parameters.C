// main cwb configuration file for the 2G pipeline
// parameters can be redefined/added in user_parameters.C
// for reference contact
// S.Klimenko, University of Florida. klimenko@phys.ufl.edu
// G.Vedovato, University of Padova, vedovato@lnl.infn.it
{
#ifndef CWB_PARAMETER_FILE
#define CWB_PARAMETER_FILE

  #include "xroot.hh"        // defines macro to manage ROOT5 vs ROOT6 

  /****** search ******/

  char analysis[8]="2G";     // 2G pipeline 
  CheckAnalysis();           // check consistency with the environment CWB_ANALYSIS
  bool online=false;         // true/false -> online/offline
  SEARCH(char) = 'r';        // see description below
  bool optim=false;          // true -> optimal resolution likelihood analysis 
  double fLow  = 64.;        // low frequency of the search
  double fHigh = 2048.;      // high frequency of the search

  /****** network configuration ******/

  char ifo[NIFO_MAX][8];
  if(NIFO_MAX>0) strcpy(ifo[0],"L1");  		 // LIGO Livingston                   
  if(NIFO_MAX>1) strcpy(ifo[1],"H1");  		 // LIGO Hanford                   
  if(NIFO_MAX>2) strcpy(ifo[2],"V1");  		 // Virgo                   
  if(NIFO_MAX>3) strcpy(ifo[3],"I1");  		 // LIGO India                   
  if(NIFO_MAX>4) strcpy(ifo[4],"J1");  		 // KAGRA Japan                  
  if(NIFO_MAX>5) strcpy(ifo[5],"G1");  		 // GEO-600 Hannover
  for(int i=6;i<NIFO_MAX;i++) strcpy(ifo[i],""); // ifo[] can be redefined by user 

  int  nIFO = 3;            // size of network starting with first detector ifo[]
  char refIFO[4] = "L1";    // reference IFO
  detectorParams detParms[NIFO_MAX];   // user defined detectors
  detectorParams _detParms = {"",0.,0.,0.,0,0.,0.,0.};
  for(int i=0;i<NIFO_MAX;i++) detParms[i] = _detParms;

  /****** cWB production thresholds & regulators ******/

  double bpp   = 0.001;      // probability for black pixel selection (netpixel) 
  double subnet = 0.7;       // [0,0.7] sub network threshold (supercluster)
  double subcut = 0.33;      // [0,1]   sub network threshold in the skyloop (supercluster) 
  double netRHO = 4;         // [>4.0] coherent network SNR (supercluster, likelihood)
  double netCC = 0.5;        // network correlation (supercluster, likelihood)
  double Acore = sqrt(2);    // threshold of core pixels (supercluster, likelihood)
  double Tgap  = 3.0;        // defragmentation time gap between clusters (sec)
  double Fgap  = 130.;       // defragmentation frequency gap between clusters (Hz)
  double TFgap =  6.;        // threshold on the time-frequency separation between two pixels 

  // regulators 
  double delta = 0.5;	     // [-1:1] regulate 2 Detector sky locations
                             // delta=0 : regulator is disabled, delta<0 : select Lo as skystat instead of Lr
  GAMMA(double) = 0.5;       // [-1:1] regulate |fx|<<|f+| and |f+|<<1 sky locations  
                             // gamma=0 : regulator is disabled, gamma<0 : sky prior is applied 

  /****** run configuration ******/

  int    runID = 0;          // run number, set in the production job
  size_t inRate= 16384;      // data sampling rate
  size_t fResample = 0;      // if>0 the inRate is resampled to fResample

  // time-frequency transformation settings
  int l_low    = 3;          // low  frequency resolution level (2^l_low Hz)
  int l_high   = 8;          // high frequency resolution level (2^l_high Hz)
  int levelR   = 2;          // resampling level : inRate[fResample]/(2^levelR) Hz

  // job segments
  double segLen     = 600.;  // Segment length [sec]
  double segMLS     = 300.;  // Minimum Segment Length after DQ_CAT1 [sec]
  double segTHR     = 30.;   // Minimum Segment Length after DQ_CAT2 [sec] (to disable put segTHR=0)
  double segEdge    = 8.;    // wavelet boundary offset [sec]
  double segOverlap = 0.;    // overlap between job segments [sec]

  // time shift analysis lags
  size_t lagSize    = 1;     // number of lags (simulation=1)
  double lagStep    = 1.;    // [sec] time interval between lags
  size_t lagOff     = 0;     // first lag id (lagOff=0 - include zero lag )
  size_t lagMax     = 150;   // 0/>0 -  standard/extended lags
  char*  lagFile    = NULL;  // lag file list
  char   lagMode[]  = "w";   // w/r  -  write/read lag list
  size_t* lagSite   = NULL;  // site index starting with 0

  // time shift analysis super lags
  int    slagSize   = 0;     // number of super lags (simulation=1) - if slagSize=0 -> Standard Segments
  int    slagMin    = 0;     // select the minimum available slag distance : slagMin must be <= slagMax
  int    slagMax    = 0;     // select the maximum available slag distance
  int    slagOff    = 0;     // first slag id (slagOff=0 - include zero slag )
  size_t* slagSite  = NULL;  // site index starting with 0
  char*  slagFile   = NULL;  // slag file list

  // whitening parameters
  double whiteWindow = 60.;  // [sec] time window dT. if = 0 - dT=T, where T is segment duration
  double whiteStride = 20.;  // [sec] noise sampling time stride

  // Skymap probability to be saved in the final output root file (saved if !=0 : see nSky)
  int Psave = 0;

  // DC corrections
  double dcCal[NIFO_MAX];
  for(int i=0;i<NIFO_MAX;i++) dcCal[i] = 1.0;
 
  // simulation parameters
  int simulation = 0;        // 1/2/3/4 for simulation(hrss/snr/tshift/multi), 0 for production
  int nfactor=0;             // number of simulation factors
  double factors[FACTORS_MAX]; // array of simulation factors (when sim=4 factors[0] is used as offset [must be int])
  for(int i=0;i<FACTORS_MAX;i++) factors[i]=0.;
  double iwindow = 5.;       // injection time window (Tinj +/- iwindow/2)

  // noise shift data
  double dataShift[NIFO_MAX];
  for(int i=0;i<NIFO_MAX;i++) dataShift[i] = 0.;

  // use this parameter to shift in time the injections (sec)
  // use {0,0,0} to set mdc_shift to 0
  // if {-1,0,0} the shift is automatically selected
  // {startMDC, stopMDC, offset}
  // see description in the method CWB::Toolbox::getMDCShift 
  mdcshift mdc_shift = {0, 0, 0};   

  // delay filter
                               // catalog of WDM cross-talk coefficients 
  char   wdmXTalk[1024] = "wdmXTalk/OverlapCatalog_Lev_8_16_32_64_128_256_iNu_4_Prec_10.bin";       
  size_t upTDF = 4;            // upsample factor to obtain rate of TD filter : TDRate = (inRate>>levelR)*upTDF 
  size_t TDSize = 12;          // time-delay filter size (max 20) 

  // coherence stage settings
  // select pixel pattern used to produce the energy max maps for pixel's selection 
  // patterns: "/" - ring-up, "\" - ring-down, "|" - delta, "-" line, "*" - single
  //
  // pattern =  0 - "*"   1-pixel  standard search
  // pattern =  1 - "3|"  3-pixels vertical packet (delta)
  // pattern =  2 - "3-"  3-pixels horizontal packet (line)
  // pattern =  3 - "3/"  3-pixels diagonal packet (ring-up)
  // pattern =  4 - "3\"  3-pixels anti-diagonal packet (ring-down)
  // pattern =  5 - "5/"  5-pixels diagonal packet (ring-up)
  // pattern =  6 - "5\"  5-pixels anti-diagonal packet (ring-down)
  // pattern =  7 - "3+"  5-pixels plus packet (plus)
  // pattern =  8 - "3x"  5-pixels cross packet (cross)
  // pattern =  9 - "9p"  9-pixels square packet (box)
  // pattern = else - "*" 1-pixel  packet (single)
  //
  // ------------------------------------------------------------------------------------
  // pattern==0                   Standard Search : std-pixel    selection + likelihood2G
  // pattern!=0 && pattern<0      Mixed    Search : packet-pixel selection + likelihood2G 
  // pattern!=0 && pattern>0      Packed   Search : packet-pixel selection + likelihoodWP

  int pattern = 0;

  // supercluster stage settings

  int    BATCH = 10000;        // max number of pixel to process in one loadTDamp batch
  int    LOUD = 200;           // number of pixel per cluster to load TD amplitudes 

  // sky settings

  bool   EFEC   =   true;    // Earth Fixed / Selestial coordinates
  size_t mode   =   0;       // sky search mode
  double angle  =   0.4;     // angular resolution
  double Theta1 =   0.;      // start theta
  double Theta2 = 180.;      // end theta
  double Phi1   =   0.;      // start theta
  double Phi2   = 360.;      // end theta
  size_t healpix= 7;         // if not 0 use healpix sky map (number of sky pixels = 12*pow(4,healpix))

  double precision = GetPrecision(0,0);    // set parameters for big clusters events management
                                           // par1 (csize) : cluster size threshold 
                                           // par2 (order) : order of healpix resampled skymap (<healpix) 
                                           // default (0,0) = disabled
                                           // if enabled the skyloop of the events with volume>=csize 
                                           // is downsampled to skymap(order)

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

  //****** directories, file names ******//

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
  sprintf(data_label,"%s",gSystem->BaseName(work_dir));
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
  dqfile DQF[DQF_MAX];

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
  if(gSystem->Getenv("_USE_PEGASUS")!=NULL) strcpy(nodedir,".");
  cout << "nodename   : " << gSystem->HostName() << endl;
  cout << "nodedir    : " << nodedir << endl;

  // get CWB_CONFIG
  char cwb_config_env[1024] = "";
  if(gSystem->Getenv("CWB_CONFIG")!=NULL) {
    strcpy(cwb_config_env,TString(gSystem->Getenv("CWB_CONFIG")).Data());
  }

  // get SITE_CLUSTER
  char site_cluster_env[1024] = "";
  if(gSystem->Getenv("SITE_CLUSTER")!=NULL) {
    strcpy(site_cluster_env,TString(gSystem->Getenv("SITE_CLUSTER")).Data());
  }

  // Plugin

  TMacro plugin;                // Macro source
  TMacro configPlugin;  	// Macro config
  plugin.SetName("");
  configPlugin.SetName("");
  char parPlugin[1024] = "";	// user defined parameters (used in the plugin)
  bool dataPlugin = false;      // if dataPlugin=true disable read data from frames
  bool mdcPlugin  = false;      // if dataPlugin=true disable read mdc from frames
  bool dcPlugin   = false;      // if dcPlugin=true disable built-in data conditioning (only 2G)
  bool cohPlugin  = false;      // if cohPlugin=true disable built-in coherence stage (only 2G)
  bool scPlugin   = false;      // if scPlugin=true disable built-in supercluster function (only 2G)
  bool outPlugin  = false;      // if outPlugin=true disable built-in output wave file (only 2G)

  char comment[1024] = "";	// user defined comment

  // 2G search modes
  //   r - un-modeled
  //   i - iota - wave: no,partial dispersion correction
  //   p - Psi - wave (no dispersion correction)
  // l,s - linear, loose linear
  // c,g - circular. loose circular
  // e,b - elliptical (no dispersion correction), b=p for now
  // 
  // low/upper case search (like 'i'/'I') & optim=false - standard MRA
  // low case search (like 'i')           & optim=true  - extract PCs from a single resolution
  // upper case (like 'I')                & optim=true  - standard single resolution analysis

  // ----------------------------------------------
  // obsolete parameters
  // ----------------------------------------------

  double &gap = iwindow;        // alias of iwindow
  int    levelF = 0;		// 1G parameter 
  int    levelD = 0;		// 1G parameter 
  int    mlagStep = 0;		// 1G parameter
  bool   eDisbalance = false;	// 1G parameter 
  double mask = 0.00;           // 1G parameter
  double Tlpr  = 0.;            // 1G parameter
  double x2or  = 0.;            // 1G parameter
  double shift[NIFO_MAX];       // 1G parameter
  char   filter[1024] = "";     // 1G parameter

#endif
}
