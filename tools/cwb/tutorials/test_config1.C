{
#ifndef CWB_PARAMETER_FILE_TEST
#define CWB_PARAMETER_FILE_TEST

  #include "xroot.hh"        // defines macro to manage ROOT5 vs ROOT6


  char analysis[8] 	= "1G";		// 1G or 2G analysis

  int nIFO		= 0;		// size of network starting with first detector ifo[]
  SEARCH(char) 		= 'r';		// see description below

  char ifo[NIFO_MAX][8];
  char refIFO[4];			// reference IFO
  // user define detectors list : is selected if detectorParams[n].name!=""
  // {name, latitude, longitude, elevation, AltX, AzX, AltY, AzY}
  detectorParams detParms[NIFO_MAX];

  // cWB settings

  size_t rateANA	= 4096;		// analysis data rate
  double bpp		= 0.000100;	// probability for pixel selection
  double Tgap		= 0.050000;	// time gap between clusters (sec)
  double Fgap		= 128.000000;	// frequency gap between clusters (Hz)
  double fLow		= 64.000000;	// low frequency of the search
  double fHigh		= 2048.000000;	// high frequency of the search
  double fResample	= 0.000000;	// if zero resampling is not applied
  double Acore		= 1.414214;	// threshold for selection of core pixels
  double Tlpr		= 120.000000;	// training time for LPR filter

  double x2or		= 1.500000;	// 2 OR threshold
  double netRHO		= 3.500000;	// threshold on rho
  double netCC		= 0.500000;	// threshold on network correlation

  // wavelet transformation settings

  int levelR		= 2;		// resampling level
  int levelF		= 6;		// level where second LPR filter is applied
  int levelD		= 8;		// decomposition level
  int l_low		= 3;		// low frequency resolution level
  int l_high		= 8;		// high frequency resolution level

  // time shift analysis settings

  // segments
  double segLen		= 600.000000;	// Segment length [sec]
  double segMLS		= 300.000000;	// Minimum Segment Length after DQ_CAT1 [sec]
  double segTHR		= 30.000000;	// Minimum Segment Length after DQ_CAT2 [sec]
  double segEdge	= 8.000000;	// wavelet boundary offset [sec]

  // lags
  size_t lagSize	= 1;		// number of lags (simulation=1)
  double lagStep	= 1.000000;	// time interval between lags [sec]
  size_t lagOff		= 0;		// first lag id (lagOff=0 - include zero lag )
  size_t lagMax		= 150;		// 0/>0 -  standard/extended lags
  char*  lagFile	= NULL;		// lag file list
  char   lagMode[2]	= "w";		// w/r  -  write/read lag list
  size_t* lagSite	= NULL;		// site index starting with 0
  double shift[NIFO_MAX] = {0,0,0,0,0,0,0,0 };	// use for standard shifts

  // multi lags
  int    mlagStep	= 0;		// if mlagStep=0 then 'standard lag mode' else cicle over lags with step mlagStep

  // super lags
  int    slagSize	= 0;		// number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  int    slagMin	= 0;		// if slagMax=0 ->  slagMin must be < slagMax
  int    slagMax	= 0;		// if slagMax=0
  int    slagOff	= 0;		// first slag id (slagOff=0 - include zero slag )
  size_t* slagSite	= NULL;		// site index starting with 0

 // whitening parameters
  double whiteWindow		= 60.000000;	// time window dT. if = 0 - dT=T, where T is wavearray duration
  double whiteStride		= 20.000000;	// noise sampling interval (window stride)

  int Psave		= 0;		//  Skymap probability pixels to be saved in the final output root file

  // DC corrections
  double dcCal[NIFO_MAX] = {0,0,0,0,0,0,0,0 };	// use for standard dcCals

  // simulation parameters
  int simulation 	= 0;		// 1 for simulation, 0 for production
  double gap 		= 5.000000;		// analysis time window for injections
  int nfactor 		= 0;		// number of strain factors
  double factors[100];			// array of strain factors

  // noise shift data
  double dataShift[NIFO_MAX] = {0,0,0,0,0,0,0,0 };	// use for standard dataShifts

  // use this parameter to shift in time the injections (sec)
  // use {0,0,0} to set mdc_shift to 0
  // if {-1,0,0} the shift is automaticaly selected
  // {startMDC, stopMDC}
  mdcshift mdc_shift = {0.000000, 0.000000, 0.000000};

  // delay filter

  double TDRate	= 16384.000000;		// time-delay filter rate
  char filter[1024] = "up2";		// delay filter suffix: "", or "up1", or "up2"

  // regulator

  double delta	= 1.000000;		//  [0/1] -> [weak/soft]
  GAMMA(double) = 0.200000;		//  set params in net5, [0/1]->net5=[nIFO/0],
  					//  if net5>[threshold=(nIFO-1)] weak/soft[according to delta] else hard
  bool   eDisbalance	= 1;

  // sky settings

  bool   EFEC	= 1;			// Earth Fixed / Selestial coordinates
  size_t mode		= 0;		// sky search mode
  double angle		= 0.400000;	// angular resolution
  double Theta1		= 0.000000;	// start theta
  double Theta2		= 180.000000;	// end theta
  double Phi1		= 0.000000;	// start theta
  double Phi2		= 360.000000;	// end theta
  double mask		= 0.000000;	// sky mask fraction
  size_t healpix	= 0;		// if not 0 use healpix sky map (healpix order)

  // error regions settings

  double precision	= 0.001000;	//  No = nIFO*(K+KZero)+precision*E

  // file dump mode

  bool saveTemp		= 0;		// save temporary root data file
  bool dumpHistory	= 1;		// dump history into output root file
  bool dump		= 1;		// dump triggers into ascii file
  bool savemode		= 1;		// temporary save clusters on disc
  bool cedDump		= 0;		// dump ced plots with rho>cedRHO
  double cedRHO		= 4.000000;
  double nSky		= 0.000000;	// # of skymap prob pixels dumped to ascii : (nSky=0 -> (#pixels==1000 || cum prob > 0.99))

  // directories, file names

  char filter_dir[1024]	= "/home/waveburst/soft/filters/data64_wat-4.8.2";
  char injectionList[1024]	= "";
  char skyMaskFile[1024]	= "";
  char skyMaskCCFile[1024]	= "";
  char channelNamesRaw[NIFO_MAX][50];
  char channelNamesMDC[NIFO_MAX][50];

  // working dir
  char work_dir[512]	= "/home/vedovato/Y2/coherent/SVN/watrepo/wat/trunk/tools/cwb/tutorials";
  char config_dir[512]	= "config";
  char input_dir[512]	= "input";
  char output_dir[512]	= "output";
  char merge_dir[512]	= "merge";
  char condor_dir[512]	= "condor";
  char report_dir[512]	= "report";
  char macro_dir[512]	= "macro";
  char log_dir[512]	= "log";
  char data_dir[512]	= "data";
  char tmp_dir[512]	= "tmp";
  char ced_dir[512]	= "report/ced";
  char pp_dir[512]	= "report/postprod";
  char dump_dir[512]	= "report/dump";
  char www_dir[512]	= "/home/vedovato/WWW/LSC/reports";

  // data label
  char data_label[512]	= "tutorials";

  // condor declarations
  char condor_log[512]	= "/local/user/vedovato";

  // frame files list (MDC must be declared in the nIFO+1 position)
  char frFiles[NIFO_MAX+1][256];

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  int nDQF = 12;
  dqfile DQF[12] = {
                   {"L1", "input/S6A_OFFLINE_L1SCIENCE.txt", (CWB_CAT)0, 0.000000, 0, 0},
                   {"L1", "input/S6A_OFFLINE_L1_DQCAT1SEGMENTS.txt", (CWB_CAT)1, 0.000000, 1, 0},
                   {"L1", "input/S6A_OFFLINE_L1_DQCAT2SEGMENTS.txt", (CWB_CAT)2, 0.000000, 1, 0},
                   {"L1", "input/S6A_OFFLINE_L1_DQCAT4SEGMENTS.txt", (CWB_CAT)1, 0.000000, 1, 0},
                   {"H1", "input/S6A_OFFLINE_H1SCIENCE.txt", (CWB_CAT)0, 0.000000, 0, 0},
                   {"H1", "input/S6A_OFFLINE_H1_DQCAT1SEGMENTS.txt", (CWB_CAT)1, 0.000000, 1, 0},
                   {"H1", "input/S6A_OFFLINE_H1_DQCAT2SEGMENTS.txt", (CWB_CAT)2, 0.000000, 1, 0},
                   {"H1", "input/S6A_OFFLINE_H1_DQCAT4SEGMENTS.txt", (CWB_CAT)1, 0.000000, 1, 0},
                   {"V1", "input/S6A_OFFLINE_V1SCIENCE.txt", (CWB_CAT)0, 0.000000, 0, 0},
                   {"V1", "input/S6A_OFFLINE_V1_DQCAT1SEGMENTS.txt", (CWB_CAT)1, 0.000000, 1, 0},
                   {"V1", "input/S6A_OFFLINE_V1_DQCAT2SEGMENTS.txt", (CWB_CAT)2, 0.000000, 1, 0},
                   {"V1", "input/S6A_OFFLINE_V1_DQCAT4SEGMENTS.txt", (CWB_CAT)1, 0.000000, 1, 0},
                 };

  // read and dump data on local disk (nodedir)
  char nodedir[1024]	= "/atlas/node/www.atlas.aei.uni-hannover.de/user/vedovato";

  // plugin name
  TMacro plugin;
  TMacro configPlugin;
  bool dataPlugin		= 0;		// if dataPlugin=true disable read data from frames
  bool mdcPlugin		= 0;		// if dataPlugin=true disable read data from frames

  // statistics
  // L  - likelihood
  // c  - network correlation coefficient
  // A  - energy disbalance asymmetry
  // P  - penalty factor based on correlation coefficients <x,s>/sqrt(<x,x>*<s,s>)
  // E  - total energy in the data streams

  // search modes
  // 'c' - un-modeled search, fast S5 cWB version, requires constraint settings
  // 'h' - un-modeled search, S5 cWB version, requires constraint settings
  // 'B' - un-modeled search, max(P*L*c/E)
  // 'b' - un-modeled search, max(P*L*c*A/E)
  // 'I' - elliptical plarisation, max(P*L*c/E)
  // 'S' - linear plarisation, max(P*L*c/E)
  // 'G' - circular plarisation, max(P*L*c/E)
  // 'i' - elliptical plarisation, max(P*L*c*A/E)
  // 's' - linear plarisation, max(P*L*c*A/E)
  // 'g' - circular plarisation, max(P*L*c*A/E)

#endif
}
