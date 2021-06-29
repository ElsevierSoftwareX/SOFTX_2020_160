/*
# Copyright (C) 2019 Gabriele Vedovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


/**********************************************************
 * Package:      config Class Library
 * File name:    config.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef CONFIG_HH
#define CONFIG_HH

#include "TObjString.h"
#include "TObjArray.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TNamed.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "wavecomplex.hh"
#include "wavearray.hh"
#include "detector.hh"
#include "wat.hh"

#include "Toolbox.hh"
#include "History.hh"

// job file options
enum CWB_JOBF_OPTIONS {
  CWB_JOBF_SAVE_DISABLE       = 0x00000000,
  CWB_JOBF_SAVE_CONFIG        = 0x00000001,
  CWB_JOBF_SAVE_NETWORK       = 0x00000002,
  CWB_JOBF_SAVE_HISTORY       = 0x00000004,
  CWB_JOBF_SAVE_STRAIN        = 0x00000008,
  CWB_JOBF_SAVE_MDC           = 0x00000010,
  CWB_JOBF_SAVE_CSTRAIN       = 0x00000020,
  CWB_JOBF_SAVE_COHERENCE     = 0x00000040,
  CWB_JOBF_SAVE_SUPERCLUSTER  = 0x00000080,
  CWB_JOBF_SAVE_LIKELIHOOD    = 0x00000100,
  CWB_JOBF_SAVE_CED           = 0x00000200,
  CWB_JOBF_SAVE_JNET_MACRO    = 0x00000400,
  CWB_JOBF_SAVE_CWB           = 0x00000800,
  CWB_JOBF_SAVE_SPARSE        = 0x00001000,
  CWB_JOBF_SAVE_WFINJ         = 0x00002000,
  CWB_JOBF_SAVE_WFREC         = 0x00004000,
  CWB_JOBF_SAVE_NODE          = 0x00008000,
  CWB_JOBF_SAVE_TRGFILE       = 0x00010000,
  CWB_JOBF_SAVE_CSPARSE       = 0x00020000,
  CWB_JOBF_SAVE_ALL           = 0xFFFFFFFF
};

// output root file options
enum CWB_OUTF_OPTIONS {
  CWB_OUTF_SAVE_DISABLE       = 0x00000000,
  CWB_OUTF_SAVE_VAR           = 0x00000001,
  CWB_OUTF_SAVE_NOISE         = 0x00000002,
  CWB_OUTF_SAVE_ALL           = 0xFFFFFFFF
};

const int FACTORS_MAX	= 1000;
const int DQF_MAX	= 100;

namespace CWB {

class config : public TNamed {

public:
  
  config(TString umacro="");         // umacro=="" -> Reset config parameters 
                                     // unamed!="" -> Import config parameters from umacro
  ~config();  

  void Init();                       // Reset parameters to 0/NULL/""/' '
  void SetVar(bool MODE);
  void Export(TString fname="");     // Export config values to CINT or to file fname
  void Import(TString umacro="");    // umacro=="" -> Import config parameters from CINT
                                     // unamed!="" -> Import config parameters from umacro
  void Print(Option_t* option="");
  int  Compare(CWB::config config);
  void View();			     					// *MENU*
  void DumpConfig(const char* filename = "", Option_t* option = "");  	// *MENU*
  void DumpPlugin(const char* filename = "");    			// *MENU*
  void DumpConfigPlugin(const char* filename = "");			// *MENU*
  void SetSingleDetectorMode();      // Set configuration to works as 1 detector 
  void Check();			     // check consistency of parameters

  virtual void Browse(TBrowser *b);  

// config parameters

  char analysis[8];        // 1G or 2G analysis
  bool online;		   // true/false -> online/offline

  int nIFO;                // size of network starting with first detector ifo[]
  char search;             // see description below
  bool optim;              // true -> optimal resolution likelihood analysis

  char ifo[NIFO_MAX][8];
  char refIFO[4];          // reference IFO
  // user define detectors list : is selected if detectorParams[n].name!=""
  // {name, latitude, longitude, elevation, AltX, AzX, AltY, AzY
  detectorParams detParms[NIFO_MAX]; 

  // cWB settings

  size_t inRate;           // input data rate 
  double bpp;              // probability for pixel selection
  double Tgap;             // 1G: time gap between clusters (sec)
                           // 2G: defragmentation time gap between clusters (sec)
  double Fgap;             // 1G: frequency gap between clusters (Hz)
                           // 2G: defragmentation frequency gap between clusters (Hz)
  double TFgap;            // threshold on the time-frequency separation between two pixels

  double fLow;             // low frequency of the search
  double fHigh;            // high frequency of the search
  size_t fResample;        // if zero resampling is not applied
  double Acore;            // threshold for selection of core pixels
  double Tlpr;             // training time for LPR filter

  double x2or;             // 2 OR threshold
  double netRHO;           // threshold on rho
  double netCC;            // threshold on network correlation

  // wavelet transformation settings

  int levelR;              // resampling level
  int levelF;              // level where second LPR filter is applied
  int levelD;              // decomposition level
  int l_low;               // low frequency resolution level
  int l_high;              // high frequency resolution level

  // time shift analysis settings

  // segments
  double segLen;           // Segment length [sec]
  double segMLS;           // Minimum Segment Length after DQ_CAT1 [sec]
  double segTHR;           // Minimum Segment Length after DQ_CAT2 [sec]
  double segEdge;          // wavelet boundary offset [sec]
  double segOverlap;       // overlap between segments [sec]

  // lags
  size_t lagSize;          // number of lags (simulation=1)
  double lagStep;          // time interval between lags [sec]
  size_t lagOff;           // first lag id (lagOff=0 - include zero lag )
  size_t lagMax;           // 0/>0 -  standard/extended lags
  char*  lagFile;          // lag file list
  char   lagMode[2];       // w/r  -  write/read lag list
  size_t* lagSite;         // site index starting with 0
  double shift[NIFO_MAX];  // use for standard shifts

  // multi lags
  int    mlagStep;         // if mlagStep=0 then 'standard lag mode' else cicle over lags with step mlagStep

  // super lags
  int    slagSize;         // number of super lags (simulation=1) - if slagSize=0 -> Igor Segments
  int    slagMin;          // if slagMax=0 ->  slagMin must be < slagMax
  int    slagMax;          // if slagMax=0
  int    slagOff;          // first slag id (slagOff=0 - include zero slag )
  size_t* slagSite;        // site index starting with 0
  char*  slagFile;         // slag file list

  // whitening parameters
  double whiteWindow;	   // time window dT. if = 0 - dT=T, where T is wavearray duration
  double whiteStride;      // noise sampling interval (window stride)

  // Skymap probability pixels to be saved in the final output root file
  int Psave;

  // DC corrections
  double dcCal[NIFO_MAX];

  // simulation parameters
  int simulation;          // 1/2/3 for simulation [1=strain, 2=snr, 3=tshift], 0 for production
  double iwindow;          // analysis time window for injections (Range = Tinj +/- iwindow/2)
  int nfactor;             // number of strain/snr factors
  double factors[FACTORS_MAX];    // array of strain/snr factors

  // noise shift data
  double dataShift[NIFO_MAX];

  // use this parameter to shift in time the injections (sec)
  // use {0,0,0} to set mdc_shift to 0
  // if {-1,0,0} the shift is automaticaly selected
  // {startMDC, stopMDC}
  mdcshift mdc_shift;

  // delay filter

  char   wdmXTalk[1024];   // 2G: catalog of WDM cross-talk coefficients 
  size_t upTDF;            // 2G: upsample factor to obtain rate of TD filter : TDrate = (inRate>>levelR)*upTDF 
  size_t TDSize;           // 2G: time-delay filter size (max 20) 
  char   filter[1024];     // 1G: delay filter suffix: "", or "up1", or "up2" 

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

  int pattern;

  // supercluster stage

  int    BATCH;            // 2G: max number of pixel to process in one loadTDamp batch 
  int    LOUD;             // 2G: number of pixel per cluster to load TD amplitudes 
  double subnet;           // 2G: sub network threshold (supercluster) 
  double subcut;           // 2G: sub network threshold in the skyloop (supercluster) 

  // regulator

  double delta;            // 1G: [0/1] -> [weak/soft]  

                           // 2G: (0. - 0.1) - typical value is around |fx|^2,
                           // 2G: if high (>0.1) and than force the hard regulator
                           // 2G: this regulator is unlikely to change

  double gamma;            // 1G: set params in net5, [0/1]->net5=[nIFO/0],
                           // if net5>[threshold=(nIFO-1)] weak/soft[according to delta] else hard 

                           // 2G:
                           // (0. - 0.1) - typical value is around |fx|^2,
                           // if high (>0.1) and than force the hard regulator
                           // this regulator is unlikely to change
                           // (0. - 1.) - defines threshold on coherent energy.
                           // If gamma=0, than all pixels with negative
                           // coherent energy are forced to produce zero signal response.
                           // It is not desirable to have gamma>1.

  bool   eDisbalance;      // 1G:

  // sky settings

  bool   EFEC;             // Earth Fixed / Selestial coordinates
  size_t mode;             // sky search mode
  double angle;            // angular resolution
  double Theta1;           // start theta
  double Theta2;           // end theta
  double Phi1;             // start theta
  double Phi2;             // end theta
  double mask;             // sky mask fraction
  size_t healpix;          // if not 0 use healpix sky map (healpix order)

  // error regions settings

  double precision;        // 1G : No = nIFO*(K+KZero)+precision*E 
                           // 2G : precision of energy calculation with time delay filter

  // file dump mode

  CWB_JOBF_OPTIONS jobfOptions; // job file options
  CWB_OUTF_OPTIONS outfOptions; // output file options

  bool dumpHistory;        // dump history into output root file
  bool dump;               // dump triggers into ascii file
  bool savemode;           // temporary save clusters on disc
  bool cedDump;            // dump ced plots with rho>cedRHO
  double cedRHO;
  long nSky;               // if nSky>0 -> # of skymap prob pixels dumped to ascii 
                           // if nSky=0 -> (#pixels==1000 || cum prob > 0.99)
                           // if nSky<0 -> nSky=-XYZ... save all pixels with prob < 0.XYZ...

  // directories, file names

  char filter_dir[1024];

  char injectionList[1024];
  char skyMaskFile[1024];
  char skyMaskCCFile[1024];
  char channelNamesRaw[NIFO_MAX][50];
  char channelNamesMDC[NIFO_MAX][50];

  // working dir
  char work_dir[1024];

  char config_dir[1024];
  char input_dir[1024];
  char output_dir[1024];
  char merge_dir[1024];
  char condor_dir[1024];
  char report_dir[1024];
  char macro_dir[1024];
  char log_dir[1024];
  char data_dir[1024];
  char tmp_dir[1024];
  char ced_dir[1024];
  char pp_dir[1024];
  char dump_dir[1024];
  char www_dir[1024];

  // data label
  char data_label[1024];

  // condor declarations
  char condor_log[1024];

  // Define a Unique Tag for Condor Jobs
  char condor_tag[1024];

  // frame files list : [0:nIFO-1]/[nIFO:2*nIFO-1] contains strain/mdc file names
  // If all mdc channels are in a single frame file -> mdc must be declared in the nIFO position
  char frFiles[2*NIFO_MAX][1024];
  // frame reading retry time (sec) : 0 -> disable
  int  frRetryTime;

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  int nDQF;
  dqfile DQF[DQF_MAX];

  // read and dump data on local disk (nodedir)
  char nodedir[1024];

  // cwb config path
  char cwb_config_env[1024];

  // cluster site name
  char site_cluster_env[1024];

  // Plugin 
  TMacro plugin;	// Macro source
  TMacro configPlugin;	// Macro config
  char parPlugin[1024]; // user defined parameters (used in the plugin)
  bool dataPlugin;	// if dataPlugin=true disable read data from frames 
  bool mdcPlugin;	// if mdcPlugin=true disable read mdc from frames
  bool dcPlugin;	// if dcPlugin=true disable built-in data conditioning (only 2G)
  bool cohPlugin;	// if cohPlugin=true disable built-in coherence stage (only 2G)
  bool scPlugin;	// if scPlugin=true disable built-in supercluster function (only 2G)
  bool outPlugin;	// if outPlugin=true disable built-in output wave file (only 2G)

  char comment[1024];   // user defined comment

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

  // 2G search modes 
  //   r - un-modeled
  //   i - iota - wave (no dispersion correction)
  //   p - Psi - wave
  // l,s - linear
  // c,g - circular
  // e,b - elliptical (no dispersion correction)
  
  ClassDef(config,24)
};  

} // end namespace

#endif
