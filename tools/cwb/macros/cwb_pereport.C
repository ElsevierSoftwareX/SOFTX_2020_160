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


// This macro is used to for the post-processing analysis of simulations
// done with PE posteriors samples
// Author : Gabriele Vedovato

#define XIFO 4

#pragma GCC system_header

#define _USE_LAL

#include "cwb.hh"
#include "config.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TMath.h"
#include "mdc.hh"
#include "cwb2G.hh"
#include "watplot.hh"
#include "gwavearray.hh"
#include "network.hh"
#include "gnetwork.hh"
#include "constants.hh"
#include <fstream>
#include <vector>
#include "TKDE.h"
#include "TSpline.h"

// ---------------------------------------------------------------------------------
// INFOS
// ---------------------------------------------------------------------------------

// the injected PE posteriors source parameters are defined in the xml file obtained
// from the conversion of the posterior_samples.dat to posterior_samples.xml
// The simulation with posteriors must be done using the plugin CWB_Plugin_WF.C which is
// used to save the inj/rec/nul waveforms into the output root file

// NOTE: the wave_dir (see input options) must include a symbolic link 
//       to the output root file of the onsource rec event
//       the name of the symbolic link must contains the string defined in data_label

// NOTE: to permit the selection of the events based on IFAR the root files 
//       in the output dir must be merged and the cwb_setifar must be applied

// ---------------------------------------------------------------------------------
// HOW RUN THIS MACRO (EXAMPLE)
// ---------------------------------------------------------------------------------

/*
root -b -l 'cwb_pereport.C+("--data_label=O2_K19_C02c_LH_BBH_SIM_tst1 \
                             --gps_event=1185389807.333 --inj_time_step=150 --off_event=0 --wave_dir=merge \
                             --nifo=2 --ifar_thr=50 --plot_type=time --plot_efraction=0.9")'

root -b -l 'cwb_pereport.C+("--data_label=O2_K21_C02c_LH_BBH_SIM_tst1 --gps_event=1187008879.45 \
                             --inj_time_step=150 --off_event=3 --wave_dir=output --nifo=2 --ifar_thr=0 --plot_type=phase \
                             --nevents=100 --sync_phase=true")'
*/

// ---------------------------------------------------------------------------------
// DEFINES : Default values used to initialize the user options input parameters
// ---------------------------------------------------------------------------------

#define PE_MAX_EVENTS		15000

#define PE_GW_NAME		"GWXXYYZZ"
#define PE_DATA_LABEL		"wave_"
#define PE_WAVE_DIR		"output"
#define PE_NIFO			2
#define PE_NCYCLES		800
#define PE_IFAR_THR		0
#define PE_GPS_EVENT		0
#define PE_TCOA_SHIFT		0.
#define PE_GPS_ERROR		0.01
#define PE_OFF_EVENT		0
#define PE_INJ_TIME_STEP	150
#define PE_INJ_ONSRC		false

#define PE_CCUT			""
#define PE_RCUT			""

#define PE_LOGL			""

#define PE_RSTAT		""

#define PE_DUMP_CR		true
#define PE_USE_ONSRC		""

#define PE_PLOT_ONSRC		false

#define PE_SYNC_PHASE		false
#define PE_SYNC_TIME		false

#define PE_SENSITIVITY		false
#define PE_RESIDUALS		false
#define PE_STATISTICS		false
#define PE_EFRACTION		0.01

#define PE_SIDE                 "full"
#define PE_COVERAGE             "local"

#define PE_PLOT_CR		true
#define PE_PLOT_TYPE		"time"
#define PE_PLOT_FNAME		"pewave"
#define PE_PLOT_ODIR		""
#define PE_PLOT_EFRACTION	0.999
#define PE_PLOT_ZOOM		1.0
#define PE_PLOT_99PERC		false

#define PE_MAKE_FF		true
#define PE_MAKE_OF		false
#define PE_MAKE_RE		true
#define PE_MAKE_NE		false
#define PE_MAKE_ONSRC		false

// ---------------------------------------------------------------------------------
// INPUT USER OPTIONS
// ---------------------------------------------------------------------------------

// rec event     : is the real detected event reconstructed by cWB
// rec posterior : is the injected posterior reconstructed by cWB 

struct uoptions {
  TString   gw_name;   		// GW event name
  TString   data_label;   	// data label of output root files (used to select files)
  TString   wave_dir;     	// wave root files directory : ex -> (output/merge) 

  int       nifo;		// number of detectors 
  int       nevents;		// max number of events read from root files (must be <= PE_MAX_EVENTS)
  int       ncycles;		// phase cycles display in plot_phase; xaxis [0:ncycles] 

  double    ifar_thr;		// ifar threshold (years) used to to select the detected events (0: -> select all)

  double    gps_event;		// gps time of rec event (first ifo)
  double    tcoa_shift;		// user tcoa time shift (default = 0 sec)
  double    gps_error;		// errors gps time of rec event used to identify event in root file
  int       off_event;		// time offset of rec_gps_time and merger time of posteriors (integer)
  int       inj_time_step;	// injection time step of posteriors (integer)
  bool      inj_onsrc;		// def=false. true -> injection are maded in at on-source time. It is setup to true when inj_time_step==0

  bool	    dump_cr;		// true/false: if enabled rec & confidence regions are dumped to ascii file 
  TString   use_onsrc;	        // if enabled (!="") the onsource vFF,vOF,vRE is replaced with 
                                //             the sample associated with median xFF.xOF,xRE onsource distribution, allowed values: 
                                //             mleft/mfull -> median of matching values left/full
                                //             eleft/efull -> median of residual energy values left/full

  bool	    sync_phase;		// true/false: if enabled the rec posteriors are 
                                // syncronized in phase respect to the rec event 
  bool	    sync_time;		// true/false: if enabled the rec posteriors are 
                                // syncronized in time respect to the rec event 
				// if sync_phase=true & sync_time=true the sync is made in both time and phase	

  bool	    sensitivity;	// true/false: if enabled the computation of sensitivity 
  bool	    residuals;		// true/false: if enabled the residuals (null) are produced instead of waveforms 
  bool	    statistics;		// true/false: if enabled the matching factors of reconstructed vs injected posteriors samples are computed
  double    efraction;	        // used to select the time range for FF/OF/RE computation

  TString   side;               // full/left/right/chirp, if left/right then FF/OF/RE are computed before/after the coalescence time
                                // if chirp then the ccut is applied according to the ccut parameters 

  TString   coverage;           // local/full - if 'full' then if coverage is computed over the time range of median waveform within the plot_efraction of energy (def=0.999)

  TString   plot_fname;		// label used for the output plot name
                                // ex: plot_name=pewave -> L1_pewave_ifar_0_efrac_0d999_Time_A33.png
  TString   plot_odir;		// output plot directory
  TString   plot_type;		// plot type used for plot_cr: available options are: time/phase/envelope/frequency/spectrum
  double    plot_efraction;	// used in time plot, select the time range based on the 
                                // signal energy fraction respect to the maximum amplitude
  double    plot_zoom;	        // used to zooming the time CR plots
  bool	    plot_cr;		// true(def)/false -> plot waveform cr
  bool	    plot_99perc;	// true/false(false) -> plot 99 percentage belt

  bool	    make_ff;		// true(def)/false -> make fitting factor
  bool	    make_of;		// true/false(def) -> make overlap factor
  bool	    make_re;		// true(def)/false -> make residual energy factor
  bool	    make_ne;		// true/false(def) -> make null energy factor
  bool	    make_onsrc;		// true/false(def) -> make onsource disctribution (green): cWB onsource vs inj-offsource 

  TString   ccut;		// Is used to select TF area using chirp TF cut. If="" then it is not applied.
				// The selected TF region is defined by the parameters declared in the ccut string.
				// format: wdm_fres:mc_lower_factor:mc_upper_factor:left_time:right_time
				// default: 16:1.35:1.65:0.5:0.0 

  // the following parameters are setup parsing the user ccut parameter. Are initialized with the default values
  int       ccut_wdm_fres; 
  double    ccut_bchirp; 
  double    ccut_uchirp; 
  double    ccut_ltime; 
  double    ccut_rtime; 
  double    ccut_fmax;

  TString   rcut;		// Is used to select TF area using pixels above the rcut_thr threshold. If="" then it is not applied.
				// The selected TF region is defined by the parameters declared in the rcut string.
				// format: wdm_fres:thr
				// default: 64:0.95

  // the following parameters are setup parsing the user rcut parameter. Are initialized with the default values
  int       rcut_wdm_fres; 
  double    rcut_thr; 

  TString   logl;		// used to define the parameters to compute logl	(EXPERIMENTAL OPTION)
				// format: enable:flow:fhigh:arthr:ifo_mask:resample
				// default: false:0:4096:0.001:0x11111111:1

  // the following parameters are setup parsing the user logl parameter. Are initialized with the default values
  bool      logl_enabled;
  float	    logl_flow;	
  float	    logl_fhigh;	
  float	    logl_arthr;	
  int       logl_ifo_mask;
  int       logl_resample;

  TString   rstat;		// used to define the parameters to compute rstat	(EXPERIMENTAL OPTION)

  // the following parameters are setup parsing the user rstat parameter. Are initialized with the default values

  bool      rstat_enabled;
  TString   rstat_type;
  int       rstat_rtrials;
  int       rstat_jtrials;

};

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions           gOPT;        		// global User Options
int                gEVENTS;     		// number of read events
static TString     gIFO[3] = {"L1","H1","V1"};	// network ifos. WARNING: this setup is hard coded !!!  
WDM<double>*       gWDM;			// wdm used for chirp TF cut 

// chirp tf cut default parameters

#define CCUT_WDM_FRES      	16	// the WDM frequency resolution (Hz) used to decompose the signal in TF domain
#define CCUT_BCHIRP		1.35	// mchirp factor used to select the bottom region boundary
#define CCUT_UCHIRP		1.65	// mchirp factor used to select the upper region boundary
#define CCUT_LTIME		0.5	// left time from tcoa (sec) used to select the the time region: [tcoa-left_time:tcoa-right_time]
#define CCUT_RTIME		0.0 	// right time from tcoa (sec) used to select the the time region: [tcoa-left_time:tcoa-right_time]
#define CCUT_FMAX		512.0 	// hard coded init value for freq max used by the tchirp and fchirp functions

#define RCUT_WDM_FRES      	8	// the WDM frequency resolution (Hz) used to decompose the signal in TF domain
#define RCUT_THR 		0.95	// normalized residual energy threshold [0,1]

// logl default parameters	(EXPERIMENTAL OPTION)

#define LOGL_ENABLED	false		// true/false -> enable/disable logl computation
#define LOGL_FLOW	0		// logl low frequency start
#define LOGL_FHIGH	8192		// logl high frequency stop
#define LOGL_ARTHR	0.001		// logl amplitude ratio threshold -> min_freq_amp vREC > max_freq_amp vREC
#define LOGL_IFO_MASK	0x1111111	// logl ifo mask used to select ifos for logl computation
#define LOGL_RESAMPLE	1		// logl resample, used to speedup logl computation. logl is computed every 'resample' samples 

//#define LOGL_USE_KDE			// select the probability distribution type: KDE or binning Histogram for logl computation
//#define LOGL_TEST_KDE	1000		// used to check comparison KDE vs binning Histogram

// rstat defaul parameters _. rstat is the statistic based on symmetric null hyphotesis distribution

#define RSTAT_ENABLED	false		// true/false -> enable/disable rstat computation
#define RSTAT_TYPE      "median"        // rstat1/median
#define RSTAT_RTRIALS	0		// number of entries in the null hyphotesis distribution. If 0 then rtrials=gEVENTS
#define RSTAT_JTRIALS	0		// number of injections used to build the 'green' distribution. If 0 then jtrials=gEVENTS

// ---------------------------------------------------------------------------------
// WAVEFORMS (vectors index is ifos)
// ---------------------------------------------------------------------------------

std::vector<wavearray<double> > vNUL;		// null reconstructed by wave packet 
std::vector<wavearray<double> > vREC;		// signal reconstructed by wave packet 
std::vector<wavearray<double> > vINJ;		// whitened injected waveform 
std::vector<wavearray<double> > vAUX;		// auxiliary whitened injected waveform 
std::vector<wavearray<double> > vMED;		// median reconstructed amplitude
std::vector<wavearray<double> > vL50;		// = vMED - 50% Lower Bound 
std::vector<wavearray<double> > vU50;		// = vMED - 50% Upper Bound
std::vector<wavearray<double> > vL90;		// = 90% Lower Bound - vMED
std::vector<wavearray<double> > vU90;		// = 90% Upper Bound - vMED
std::vector<wavearray<double> > vL99;		// = 99% Lower Bound - vMED
std::vector<wavearray<double> > vU99;		// = 99% Upper Bound - vMED

std::vector<wavearray<double> > pSNR;		// = on-source SNR  evolution of vREC computed with cWB::Toolbos::GetPhase
std::vector<wavearray<double> > pFRQ;		// = on-source Freq evolution of vREC computed with cWB::Toolbos::GetPhase
std::vector<wavearray<double> > pTIM;		// = on-source Time evolution of vREC computed with cWB::Toolbos::GetPhase

std::vector<double> vTREC; 			// absolute tstart of vREC. The value stored into wavearray start is module inj_time_step.

double vTCOA; 					// tcoa    	OnSource sample
double vFACTOR; 				// factor  	OnSource sample
float  vTHETA; 					// theta   	OnSource sample
float  vPHI;   					// phi     	OnSource sample
float  vM1;   					// mass1 	OnSource sample
float  vM2;   					// mass2 	OnSource sample
float  vS1;   					// spin1 	OnSource sample
float  vS2;   					// spin2 	OnSource sample
float  vDOF;   					// Degrees of Freedom 	OnSource sample
float  vSNR[NIFO_MAX];                          // reconstructed SNR	OnSource sample

double vTSYNC[NIFO_MAX];                        // time shift OnSource sample after time sync
double vPSYNC[NIFO_MAX];                        // phase shift OnSource sample after phase sync

std::vector<wavearray<double> > wREC[PE_MAX_EVENTS]; 	// reconstructed OffSource signals
std::vector<wavearray<double> > wNUL[PE_MAX_EVENTS]; 	// reconstructed null OffSource signal sdata
std::vector<wavearray<double> > wINJ[PE_MAX_EVENTS]; 	// OffSource iwhitened injected posteriors
std::vector<wavearray<double> > wAUX[PE_MAX_EVENTS]; 	// OffSource auxiliary iwhitened injected posteriors

double wTCOA[PE_MAX_EVENTS]; 			// tcoa   	OffSource injected posteriors
double wFACTOR[PE_MAX_EVENTS]; 			// factor 	OffSource injected posteriors
float  wTHETA[PE_MAX_EVENTS]; 			// theta  	OffSource injected posteriors
float  wPHI[PE_MAX_EVENTS];   			// phi    	OffSource injected posteriors
float  wM1[PE_MAX_EVENTS];   			// mass1 	OffSource injected posteriors
float  wM2[PE_MAX_EVENTS];   			// mass2 	OffSource injected posteriors
float  wS1[PE_MAX_EVENTS];   			// spin1 	OffSource injected posteriors
float  wS2[PE_MAX_EVENTS];   			// spin2 	OffSource injected posteriors
float  wDOF[PE_MAX_EVENTS];   			// Degrees of Freedom 	OffSource injected posteriors
float  wSNR[PE_MAX_EVENTS][NIFO_MAX];           // reconstructed SNR	OffSource injected posteriors

std::vector<double> vFF,vOF,vRE,vNE;	// Fitting/Overlap/ResidualEnergy/NullEnergy Factors of reconstructed vs injected Posterior Samples
					// the first element is the 
                                        // cWB point estimate Fitting/Overlap/ResidualEnergy Factor respect to Map Posterior Sample

std::vector<double> xFF,xOF,xRE,xNE;    // Matching factors onsource reconstructed vs offsource injected 


// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions(TString ofname="");

void LoadOutputRootFiles(int nIFO, TString data_label, TString output_dir, int max_events=PE_MAX_EVENTS, double ifar_thr=0);
int  SyncWaveforms(int nIFO, TString type);

TString GetParameterFromInjLog(TString log, TString param);

void GetCCutParms(TString options);
int  ComputeCCutStatistics(int nIFO);
double GetCCut(wavearray<double>* ts, double tcoa, double m1, double m2, double s1z, double s2z);

wavearray<double> GetRCut(wavearray<double>* ts, wavearray<double>* aux);
void GetRCutParms(TString options);

void GetLoglParms(TString options);

void GetRstatParms(TString options);

int  ComputeStatistics(int nIFO);
int  GetOnSourceStatistic(int nIFO, double &mFF, double &mOF, double &mRE);
void GetOffSourceStatistic(int nIFO, int recId, vector<int> injId,
                           vector<double>& oFF, vector<double>& oOF, vector<double>& oRE);
void MakeDistributions(TString type="FF");
void MakePvalueDistribution(TString type="FF");
void DumpStatistics(TString ofname, TString params, bool app); 

float TestStatistic(int nIFO, int id, float ref); // Used for tests 

int  ComputeWaveformCR(int nIFO);
void DumpWaveformCR(int nIFO);
void PlotWaveformsCR(int nIFO);

int  ComputeWaveformFCR(int nIFO);
void ComputeWaveformCR(int nIFO, int selected, double* cov50, double* cov90, double* cov99);
void GetFullCoverage(int nIFO, int selected, double* fcov50, double* fcov90, double* fcov99);

void PlotWaveforms(int nIFO, int id, bool residual=false);
void PlotWaveformAsymmErrors(TString ofname, TString title, wavearray<double> wrec,
                             wavearray<double> wmed, wavearray<double> wl50, wavearray<double> wu50,
                             wavearray<double> wl90, wavearray<double> wu90, 
                             wavearray<double> wl99, wavearray<double> wu99, wavearray<double> wref, 
                             TString pdir, double P, double Q, double tref, double tcoa, int ifoid);

double GetInjTcoa(double geocentric_tcoa, network* NET, TString ifo, double theta, double phi);
double GetDelay(network* NET, TString ifo, double theta, double phi);
double GetSyncTcoa(wavearray<double> wf, double sync_time, double sync_phase, double tcoa);

void ComputePhaseSensitivity(TString match, double lower, double upper, double& dmatch, double& dphase);
void ComputeTimeSensitivity(TString match, double lower, double upper, double& dmatch, double& dtime);
void ComputeAmpSensitivity(TString match, double lower, double upper, double& dmatch, double& damp);

void cwb_pereport(TString options) {

  // Read Input Options
  ResetUserOptions();
  ReadUserOptions(options);
  PrintUserOptions();
  PrintUserOptions("pereport_config.txt");

  if(gOPT.nifo>3) {
    cout << "cwb_pereport - Error : max nifo is 3" << endl;
    gSystem->Exit(1);
  }
  if((gOPT.plot_type!="time")&&(gOPT.plot_type!="phase")&&(gOPT.plot_type!="envelope")&&(gOPT.plot_type!="frequency")&&(gOPT.plot_type!="spectrum")) {
    cout << "cwb_pereport - Error : plot_type option (" << gOPT.plot_type << ") not allowed. Available options are: time/phase/envelope/frequency/spectrum" << endl;
    gSystem->Exit(1);
  }
  if((gOPT.coverage!="local")&&(gOPT.coverage!="full")) {
    cout << "cwb_pereport - Error : coverage option (" << gOPT.coverage << ") not allowed. Available options are: local/full" << endl;
    gSystem->Exit(1);
  }

  gEVENTS=0;	// init the read events

  // init vTSYNC, vPSYNC
  for(int n=0;n<NIFO_MAX;n++) vTSYNC[n]=0;
  for(int n=0;n<NIFO_MAX;n++) vPSYNC[n]=0;

  // read reconstructed posteriors from output root files
  LoadOutputRootFiles(gOPT.nifo, gOPT.data_label, gOPT.wave_dir, gOPT.nevents, gOPT.ifar_thr);

  // init WDM transform
  if(gOPT.side=="ccut") {	// WDM transform for ccut
    int ccut_wdm_levels = int(vREC[0].rate()/gOPT.ccut_wdm_fres/2);
    cout << endl << "RATE " << vREC[0].rate() << " WDM levels " << ccut_wdm_levels << endl << endl;
    gWDM = new WDM<double>(ccut_wdm_levels, ccut_wdm_levels, 6, 10);
  }
  if(gOPT.side=="rcut") {	// WDM transform for rcut
    int rcut_wdm_levels = int(vREC[0].rate()/gOPT.rcut_wdm_fres/2);
    cout << endl << "RATE " << vREC[0].rate() << " WDM levels " << rcut_wdm_levels << endl << endl;
    gWDM = new WDM<double>(rcut_wdm_levels, rcut_wdm_levels, 6, 10);
  }

  gStyle->SetTitleW(0.95);	// fix weird title display offset

  // compute Sample Statistics (Fitting Factors, Overlap Factors, Residual Energy)
  if(gOPT.statistics) {

    if(gOPT.side=="ccut") { 
      int selected = ComputeCCutStatistics(gOPT.nifo);
      if(selected==0) exit(0);

      gOPT.sensitivity=false;
      gOPT.make_onsrc=false;
      MakeDistributions("RE"); 	// plot Residual Energy Factor

    } else {		  
      int selected = ComputeStatistics(gOPT.nifo);
      if(selected==0) exit(0);

      if(!gOPT.rstat_enabled) PlotWaveforms(gOPT.nifo, -1, true);	// plot onsource residuals
      if(!gOPT.rstat_enabled) PlotWaveforms(gOPT.nifo, -1, false);	// plot onsource waveforms
      //PlotWaveforms(gOPT.nifo, 10, true);	// plot offsource residuals
      //PlotWaveforms(gOPT.nifo, 10, false);	// plot offsource waveforms

      // plot Fitting Factor
      if(gOPT.make_ff) MakeDistributions("FF");
      if(gOPT.make_ff && gOPT.make_onsrc) MakePvalueDistribution("FF");
      // plot Overlap Factor
      if(gOPT.make_of) MakeDistributions("OF");
      if(gOPT.make_of && gOPT.make_onsrc) MakePvalueDistribution("OF");
      // plot Residual Energy Factor
      if(gOPT.make_re) MakeDistributions("RE");
      if(gOPT.make_re && gOPT.make_onsrc) MakePvalueDistribution("RE");
      // plot Null Energy Factor
      if(gOPT.make_ne && (gOPT.residuals)) MakeDistributions("NE");
    }
  } else {	// compute waveform CR

    // compute Posterior CR
    int selected = (gOPT.coverage=="local") ? ComputeWaveformCR(gOPT.nifo) : ComputeWaveformFCR(gOPT.nifo);
    if(selected==0) exit(0);

    // plot cWB rec and posterior confidence regions
    if(gOPT.plot_cr) PlotWaveformsCR(gOPT.nifo);

    // rec & confidence regions are dumped to acii file (only for time domain plots)
    if(gOPT.dump_cr && gOPT.plot_type!="phase") DumpWaveformCR(gOPT.nifo);
  }

  delete gWDM;

  exit(0);
}

int ComputeStatistics(int nIFO) {

  while(!vFF.empty()) vFF.pop_back(); vFF.clear(); std::vector<double>().swap(vFF);	// reset off-source vector vFF
  while(!vOF.empty()) vOF.pop_back(); vOF.clear(); std::vector<double>().swap(vOF);	// reset off-source vector vOF
  while(!vRE.empty()) vRE.pop_back(); vRE.clear(); std::vector<double>().swap(vRE);	// reset off-source vector vRE

  while(!xFF.empty()) xFF.pop_back(); xFF.clear(); std::vector<double>().swap(xFF);	// reset on-source  vector xFF
  while(!xOF.empty()) xOF.pop_back(); xOF.clear(); std::vector<double>().swap(xOF);	// reset on-source  vector xOF
  while(!xRE.empty()) xRE.pop_back(); xRE.clear(); std::vector<double>().swap(xRE);	// reset on-source  vector xRE
  while(!xNE.empty()) xNE.pop_back(); xNE.clear(); std::vector<double>().swap(xNE);	// reset on-source  vector xNE

  if(gOPT.rstat_enabled==false) {	// compute standard statistic

    gnetwork* NET = new gnetwork(nIFO,gIFO);

    // The waveforms are aligned and, if requested, sync in phase and time
    // the reference waveforms are INJ, this choice is necessary to fix the original onsource waveform
    int selected=SyncWaveforms(nIFO, "statistics");	
    cout << endl << "ComputeStatistics - selected events : " << selected << endl << endl;
    if(selected==0) return 0;

    // cWB point estimate OnSource Matching Factors 

    vector<double> tstart;
    if(gOPT.side=="right") {
      tstart.resize(nIFO);
      for(int n=0;n<nIFO;n++) {
        double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
        tstart[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
      }
    }

    vector<double> tstop;
    if(gOPT.side=="left") {
      tstop.resize(nIFO);
      for(int n=0;n<nIFO;n++) {
        double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
        tstop[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
      }
    }

    double FF = CWB::mdc::GetMatchFactor("ff",vREC,vINJ,tstart,tstop);
    double OF = CWB::mdc::GetMatchFactor("of",vREC,vINJ,tstart,tstop);
    double RE = CWB::mdc::GetMatchFactor("re",vREC,vINJ,tstart,tstop);

    vFF.push_back(FF);
    vOF.push_back(OF);
    vRE.push_back(RE);

    if(gOPT.residuals) {	// injection null energy
      vector<wavearray<double> > vDAT = vNUL; for(int n=0;n<nIFO;n++) vDAT[n]+=vREC[n];
      double NE = CWB::mdc::GetMatchFactor("re",vDAT,vINJ,tstart,tstop);
      vNE.push_back(NE);
    }
 
    // cWB point estimate OffSource Matching Factors

    for(int i=0;i<gEVENTS;i++) {            

      if(wREC[i].size()==0) continue;

      vector<double> tstart;
      if(gOPT.side=="right") {
        tstart.resize(nIFO);
        for(int n=0;n<nIFO;n++) {
          double inj_tcoa = gOPT.inj_time_step+fmod(wTCOA[i],1);
          tstart[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], wTHETA[i], wPHI[i]);
        }
      }

      vector<double> tstop;
      if(gOPT.side=="left") {
        tstop.resize(nIFO);
        for(int n=0;n<nIFO;n++) {
          double inj_tcoa = gOPT.inj_time_step+fmod(wTCOA[i],1);
          tstop[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], wTHETA[i], wPHI[i]);
        }
      }

      double FF = CWB::mdc::GetMatchFactor("ff",wREC[i],wINJ[i],tstart,tstop);
      double OF = CWB::mdc::GetMatchFactor("of",wREC[i],wINJ[i],tstart,tstop);
      double RE = CWB::mdc::GetMatchFactor("re",wREC[i],wINJ[i],tstart,tstop);

      vFF.push_back(FF);
      vOF.push_back(OF);
      vRE.push_back(RE);

      //float pvalue = TestStatistic(nIFO, i, RE);
      //cout << i << "\tpvalue " << pvalue << "\tref " << RE << endl;

      if(gOPT.residuals) {	// injection null energy
        vector<wavearray<double> > vDAT = wNUL[i]; for(int n=0;n<nIFO;n++) vDAT[n]+=wREC[i][n];
        double NE = CWB::mdc::GetMatchFactor("re",vDAT,wINJ[i],tstart,tstop);
        vNE.push_back(NE);
      }
    }

    // dump match factors
    TString pdir=gOPT.plot_odir;
    char ofname[256] = "match_factors.txt";
    TString fName = pdir!="" ? TString(pdir)+TString("/")+ofname : ofname;
    cout << "write : " << fName << endl;
    ofstream out;
    out.open(fName.Data(),ios::out);
    if (!out.good()) {cout << "Error Opening Output File : " << fName << endl;exit(1);}
    out.precision(14);
    double iSNR=0; 
    double oSNR=0; 
    for(int n=0;n<vINJ.size();n++) for(int i=0;i<vINJ[n].size();i++) iSNR+=pow(vINJ[n][i],2);	// inj onsource SNR
    for(int n=0;n<vREC.size();n++) for(int i=0;i<vREC[n].size();i++) oSNR+=pow(vREC[n][i],2);	// rec onsource SNR
    if(gOPT.residuals) {
      out << "tcoa " << vTCOA << "\tiSNR " << sqrt(iSNR) <<  "\toSNR " << sqrt(oSNR) << "\tfactor " << vFACTOR
          << "\tff " << vFF[0] << "\tof " << vOF[0] << "\tre " << vRE[0] << "\tne " << vNE[0] << endl;
    } else {
      out << "tcoa " << vTCOA << "\tiSNR " << sqrt(iSNR) <<  "\toSNR " << sqrt(oSNR)
          << "\tff " << vFF[0] << "\tof " << vOF[0] << "\tre " << vRE[0] << endl;
    }
    for(int j=0;j<gEVENTS;j++) {
      iSNR=oSNR=0;
      for(int n=0;n<wINJ[j].size();n++) for(int i=0;i<wINJ[j][n].size();i++) iSNR+=pow(wINJ[j][n][i],2);	// inj offsource SNR
      for(int n=0;n<wREC[j].size();n++) for(int i=0;i<wREC[j][n].size();i++) oSNR+=pow(wREC[j][n][i],2);	// rec offsource SNR
      if(gOPT.residuals) {
        out << "tcoa " << wTCOA[j] << "\tiSNR " << sqrt(iSNR) <<  "\toSNR " << sqrt(oSNR) << "\tfactor " << wFACTOR[j]
            << "\tff " << vFF[j+1] << "\tof " << vOF[j+1] << "\tre " << vRE[j+1] << "\tne " << vNE[j+1] << endl;
      } else {
        out << "tcoa " << wTCOA[j] << "\tiSNR " << sqrt(iSNR) <<  "\toSNR " << sqrt(oSNR)
            << "\tff " << vFF[j+1] << "\tof " << vOF[j+1] << "\tre " << vRE[j+1] << endl;
      }
    }
    out.close();

    if(gOPT.make_onsrc || gOPT.use_onsrc!="") {

      double mFF,mOF,mRE;
      int imed=GetOnSourceStatistic(nIFO,mFF,mOF,mRE);	// fill xFF,xOF,xRE and compute median xFF,xOF,xRE

      // the onsource vFF,vOF,vRE,vINJ are replaced with the sample associated with median xFF onsource distribution 
      if(gOPT.use_onsrc!="") {
        cout << endl << "GetOnSourceStatistic : " << "\tmFF = " << mFF  << "\tmOF = " << mOF << "\tmRE = " << mRE << endl << endl;
        vFF[0]=mFF;
        vOF[0]=mOF;
        vRE[0]=mRE;

        for(int n=0;n<nIFO;n++) vINJ[n]=wINJ[imed][n];
      }
    }
 
    delete NET;
    return selected;

  } else { 	// compute rstatistic

    // init
    if(gOPT.rstat_rtrials==0) {gOPT.rstat_rtrials=gEVENTS;cout<<endl << "---> Set rstat_rtrials = " << gOPT.rstat_rtrials << endl<<endl;}
    if(gOPT.rstat_jtrials==0) {gOPT.rstat_jtrials=gEVENTS;cout<<endl << "---> Set rstat_jtrials = " << gOPT.rstat_jtrials << endl<<endl;}
    if(gOPT.rstat_rtrials>gEVENTS) gOPT.rstat_rtrials=gEVENTS;
    if(gOPT.rstat_jtrials>gEVENTS) gOPT.rstat_jtrials=gEVENTS;

    // fill recId with random reconstructed off-source events
    vector<int> recId;
    std::map<int,bool> runique;
    while(recId.size()<gOPT.rstat_rtrials) {
      int k = int(gRandom->Uniform(0,gEVENTS));
      if(runique[k]==false) {runique[k]=true;recId.push_back(k);}
    }

    // fill injId with random reconstructed off-source events
    vector<int> injId;
    std::map<int,bool> junique;
    while(injId.size()<gOPT.rstat_jtrials) {
      int k = int(gRandom->Uniform(0,gEVENTS));
      if(junique[k]==false) {junique[k]=true;injId.push_back(k);}
    }

    // get match reconstructed offsource vs injected offsource
    int kEVENTS = gOPT.rstat_rtrials;
    int pc = 0; int npc= 0;
    int ipc = double(kEVENTS)/10.; if(ipc==0) ipc=1;
    cout << endl << "Compute offsource distribution : be patient, it takes a while ..." << endl << endl;
    vector<vector<double> > voFF,voOF,voRE;
    for(int k=0;k<kEVENTS;k++) {
      if((npc++)%ipc==0) {if(kEVENTS>10) {cout << pc<<"%";if (pc<100) cout << " - ";pc+=10;cout.flush();}}
      vector<double> oFF,oOF,oRE;
      GetOffSourceStatistic(nIFO,recId[k],injId,oFF,oOF,oRE);	// fill oFF,oOF,oRE
      voFF.push_back(oFF);
      voOF.push_back(oOF);
      voRE.push_back(oRE);
    }
    if(pc==100) cout << pc<<"%";;
    cout << endl << endl;

    // compute null hypothesis r-distribution

    vFF.push_back(0);					// reserve place on-source FF r-statistic
    vOF.push_back(0);					// reserve place on-source OF r-statistic
    vRE.push_back(0);					// reserve place on-source RE r-statistic

    if(gOPT.rstat_type=="rstat1") {
      for(int k=0;k<kEVENTS;k++) {
        int nTOT=0; int nFF=0; int nOF=0; int nRE=0;
        for(int i=0;i<kEVENTS;i++) {
          for(int j=0;j<gOPT.rstat_jtrials;j++) {
            nTOT++;
            if(voFF[k][j]>voFF[i][j]) nFF++;
            if(voOF[k][j]>voOF[i][j]) nOF++;
            if(voRE[k][j]>voRE[i][j]) nRE++;
          }
        }
        // save off-source r-statistic
        vFF.push_back(double(nFF)/double(nTOT));
        vOF.push_back(double(nOF)/double(nTOT));
        vRE.push_back(double(nRE)/double(nTOT));
      }
    }

    if(gOPT.rstat_type=="median") {
      int nentries = gOPT.rstat_jtrials;
      int *index = new int[nentries];
      double *value = new double[nentries];
      int imed = (nentries*50.)/100.; if(imed>=nentries) imed=nentries-1;

      for(int k=0;k<kEVENTS;k++) {
        for(int i=0;i<nentries;i++) value[i]=voFF[k][i];
        TMath::Sort(nentries,value,index,false);
        vFF.push_back(value[index[imed]]);		// save off-source FF r-statistic

        for(int i=0;i<nentries;i++) value[i]=voOF[k][i];
        TMath::Sort(nentries,value,index,false);
        vOF.push_back(value[index[imed]]);		// save off-source OF r-statistic

        for(int i=0;i<nentries;i++) value[i]=voRE[k][i];
        TMath::Sort(nentries,value,index,false);
        vRE.push_back(value[index[imed]]);		// save off-source RE r-statistic
      }

      delete [] index;
      delete [] value;
    }

    // get match reconstructed onsource vs injected offsource
    GetOffSourceStatistic(nIFO,-1,injId,xFF,xOF,xRE);	// fill xFF,xOF,xRE (green distribution)

    // compute on-source r-statistic value

    if(gOPT.rstat_type=="rstat1") {
      int nTOT=0; int nFF=0; int nOF=0; int nRE=0;
      for(int i=0;i<kEVENTS;i++) {
        for(int j=0;j<gOPT.rstat_jtrials;j++) {
          nTOT++;
          if(xFF[j]>voFF[i][j]) nFF++;
          if(xOF[j]>voOF[i][j]) nOF++;
          if(xRE[j]>voRE[i][j]) nRE++;
        }
      }
      // save on-source r-statistic
      vFF[0]=double(nFF)/double(nTOT);
      vOF[0]=double(nOF)/double(nTOT);
      vRE[0]=double(nRE)/double(nTOT);

      gOPT.make_onsrc=false;				// disable plot on-source distribution
    }

    if(gOPT.rstat_type=="median") {
      int nentries = gOPT.rstat_jtrials;
      int *index = new int[nentries];
      double *value = new double[nentries];
      int imed = (nentries*50.)/100.; if(imed>=nentries) imed=nentries-1;

      for(int i=0;i<nentries;i++) value[i]=xFF[i];
      TMath::Sort(nentries,value,index,false);
      vFF[0]=value[index[imed]];		// save on-source FF r-statistic

      for(int i=0;i<nentries;i++) value[i]=xOF[i];
      TMath::Sort(nentries,value,index,false);
      vOF[0]=value[index[imed]];		// save on-source OF r-statistic

      for(int i=0;i<nentries;i++) value[i]=xRE[i];
      TMath::Sort(nentries,value,index,false);
      vRE[0]=value[index[imed]];		// save on-source RE r-statistic

      delete [] index;
      delete [] value;
    }

    return kEVENTS;
  }

  return 0;
}

int ComputeCCutStatistics(int nIFO) {

  if(gOPT.side!="ccut") {cout << "ComputeCCutStatistics Error: enabled only with gOPT.side=ccut" << endl;exit(1);}

  // The waveforms are aligned and, if requested, sync in phase and time
  // the reference waveforms are INJ, this choice is necessary to fix the original onsource waveform
  int selected=SyncWaveforms(nIFO, "statistics");	
  cout << endl << "ComputeCCutStatistics - selected events : " << selected << endl << endl;
  if(selected==0) return 0;

  cout << endl << "Apply CCut : be patient, it takes a while ..." << endl << endl;
  gnetwork NET(nIFO,gIFO);
  int pc = 0;
  int npc= 0;
  int ipc = double(nIFO*gEVENTS)/10.; if(ipc==0) ipc=1;

  double RE=0;
  for(int n=0;n<nIFO;n++) {                      
    double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1)-vINJ[n].start();
           inj_tcoa = GetInjTcoa(inj_tcoa, &NET, gIFO[n], vTHETA, vPHI);
    wavearray<double> vDif = CWB::mdc::GetDiff(&vREC[n], &vINJ[n]);
    RE += GetCCut(&vDif, inj_tcoa, vM1, vM2, vS1, vS2);
  }
  vFF.push_back(0);		// dummy value
  vOF.push_back(0);		// dummy value
  vRE.push_back(RE);

  for(int i=0;i<gEVENTS;i++) {            
    RE=0;
    for(int n=0;n<nIFO;n++) {                      
      if((npc++)%ipc==0) {if(gEVENTS>100) {cout << pc<<"%";if (pc<100) cout << " - ";pc+=10;cout.flush();}}
      double inj_tcoa = gOPT.inj_time_step+fmod(wTCOA[i],1)-wINJ[i][n].start();
             inj_tcoa = GetInjTcoa(inj_tcoa, &NET, gIFO[n], wTHETA[i], wPHI[i]);
      wavearray<double> wDif = CWB::mdc::GetDiff(&wREC[i][n], &wINJ[i][n]); 
      RE += GetCCut(&wDif, inj_tcoa, wM1[i], wM2[i], wS1[i], wS2[i]);
    } 
    vFF.push_back(0);		// dummy value
    vOF.push_back(0);		// dummy value
    vRE.push_back(RE);
  }                                                                    

  if(pc==100) cout << pc<<"%";;
  cout << endl << endl;

  return selected;
}

void
PlotWaveforms(int nIFO, int id, bool residual) {

  if(id>=gEVENTS) {
    cout << "PlotWaveforms - Error : sample id exceed max number " << gEVENTS << endl;
    gSystem->Exit(1);
  }

  wavearray<double> wrec,winj;
  gnetwork* NET = new gnetwork(nIFO,gIFO);

  // plots rec vs inj
  for(int n=0;n<nIFO;n++) {

    TString tagid=""; 

    if(id<0) {wrec=vREC[n];winj=vINJ[n];tagid="onsource";}				// onsource posterior sample
    else     {wrec=wREC[id][n];winj=wINJ[id][n];tagid=TString::Format("sample %d",id);}	// offsource posterior sample

    // restore original time/phase
    if(gOPT.sync_phase || gOPT.sync_time) {
      CWB::mdc::TimePhaseSync(wrec,-vTSYNC[n],-vPSYNC[n]);
      CWB::mdc::TimePhaseSync(winj,-vTSYNC[n],-vPSYNC[n]);
    }

    // resample data to increase the plot resolution 
    wrec.Resample(16384); 
    winj.Resample(16384); 

    gwavearray<double> gw(&wrec);
    gw.Draw(GWAT_TIME);                             
    if(residual) {
      wavearray<double> wres=wrec;
      for(int i=0;i<wres.size();i++) wres[i]-=winj[i]; 
      gw.Draw(&wres,GWAT_TIME,"SAME",kRed);      
    } else {  
      gw.Draw(&winj,GWAT_TIME,"SAME",kRed);      
    }

    watplot* plot = gw.GetWATPLOT();              // get pointer to watplot object

    // set x,y axis titles and plot title
    char gtitle[256];
    if(residual) {
      sprintf(gtitle,"%s - %s - Residual(Red) - Reconstructed(Black) - Tcoa(Green) - (%s)",gOPT.gw_name.Data(),gIFO[n].Data(),tagid.Data());
    } else {
      sprintf(gtitle,"%s - %s - Injected(Red) - Reconstructed(Black) - Tcoa(Green) - (%s)",gOPT.gw_name.Data(),gIFO[n].Data(),tagid.Data());
    }
    plot->gtitle(gtitle,"time(sec)","amplitude");

    // set the time range to be displayed
    double bT, eT;
    double P=0.999;
    double tmax; double amax=0; 			// get time of max abs(amp) waveform envelope
    wavearray<double> W=wrec; W+=winj; 			// wrec+winj is used to compute the boundaries
    wavearray<double> wenv = CWB::mdc::GetEnvelope(&W);
    for(int i=0;i<wenv.size();i++) if(fabs(wenv[i])>amax) {amax=fabs(wenv[i]);tmax=wenv.start()+i/wenv.rate();}
    CWB::mdc::GetTimeBoundaries(W, P, bT, eT, tmax, P);
    bT-=W.start(); eT-=W.start();
    plot->graph[0]->GetHistogram()->GetXaxis()->SetRangeUser(bT,eT);

    // set x title
    TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %d",int(vTREC[n]));
    plot->graph[0]->GetHistogram()->GetXaxis()->SetTitle(xtitle.Data());

    bool is_tcoa_defined=false; 
    double inj_tcoa;
    if(id<0) {							// onsource posterior sample	
      if(vTCOA>=0) is_tcoa_defined=true; 
      inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
      inj_tcoa = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
      if(gOPT.sync_phase || gOPT.sync_time) {   		// inj_tcoa of vINJ[n] after time/phase sync
        inj_tcoa = GetSyncTcoa(vINJ[n], -vTSYNC[n], -vPSYNC[n], inj_tcoa);
      }
    } else {							// offsource posterior sample
      if(wTCOA[id]>=0) is_tcoa_defined=true; 
      inj_tcoa = gOPT.inj_time_step+fmod(wTCOA[id],1);
      inj_tcoa = GetInjTcoa(inj_tcoa, NET, gIFO[n], wTHETA[id], wPHI[id]);
      if(gOPT.sync_phase || gOPT.sync_time) {   		// inj_tcoa of wINJ[id][n] after time/phase sync
        inj_tcoa = GetSyncTcoa(wINJ[id][n], -vTSYNC[n], -vPSYNC[n], inj_tcoa);
      }
    }

    // add vertical green line at tcoa time
    wavearray<double> wline = wrec;
    double dt=1./wline.rate();
    for(int i=0;i<wline.size();i++) {
      double time=i*dt+wline.start();
      if(fabs(time-inj_tcoa)<dt) wline[i]=-1.0; else wline[i]=-1000.;
    }
    if(is_tcoa_defined) gw.Draw(&wline,GWAT_TIME,"SAME",kGreen);      

    TString gfile;
    if(residual) {
      gfile="x"+gIFO[n]+"_rec_vs_res_"+tagid+".png";
    } else {
      gfile="x"+gIFO[n]+"_rec_vs_inj_"+tagid+".png";
    }

    (*plot) >> gfile;

    if(residual) {
      gfile="x"+gIFO[n]+"_rec_vs_res_"+tagid+".root";
    } else {
      gfile="x"+gIFO[n]+"_rec_vs_inj_"+tagid+".root";
    }

    (*plot) >> gfile;

    cout << "created plot file name : " << gfile << endl;
  }
}

int ComputeWaveformCR(int nIFO) {

  // The waveforms are aligned and, if requested, sync in phase and time
  // the reference waveform is REC, this choice is necessary to fix the original onsource waveform
  int selected=SyncWaveforms(nIFO, "waveform");	
  cout << endl << "ComputeWaveformCR - selected events : " << selected << endl << endl;
  if(selected==0) return 0;

  if(gOPT.logl_enabled && gOPT.plot_type=="spectrum") {
    cout << endl << "Compute log-likelihood : be patient, it takes a while ..." << endl << endl;
  }
  double fmin=1e20;
  double fmax=0;
  double logl;			// Log_e Likelihood (used to compute the Bayes Factor)
  int nlogl=0;			// number of samples used in logl
  int pc = 0;
  int ipc = double(nIFO*wREC[0][0].size())/10.; if(ipc==0) ipc=1;

  float l90_dof,u90_dof;
  if(gOPT.logl_enabled && gOPT.plot_type=="spectrum") {
    int *index = new int[gEVENTS];
    TMath::Sort(gEVENTS,wDOF,index,false);
    int il90 = (gEVENTS*5.)/100.;   if(il90>=gEVENTS) il90=gEVENTS-1;
    int iu90 = (gEVENTS*95.)/100.;  if(iu90>=gEVENTS) iu90=gEVENTS-1;
    l90_dof = wDOF[index[il90]];
    u90_dof = wDOF[index[iu90]];
  }

  // compute vMED, vL50, vU50, vL90, vU90, vL99, vU99
  wavearray<double> wmed[NIFO_MAX];                
  wavearray<double> wl50[NIFO_MAX];                
  wavearray<double> wu50[NIFO_MAX];                
  wavearray<double> wl90[NIFO_MAX];                
  wavearray<double> wu90[NIFO_MAX];                
  wavearray<double> wl99[NIFO_MAX];                
  wavearray<double> wu99[NIFO_MAX];                
  for(int n=0;n<nIFO;n++) {                      
    wmed[n] = wREC[0][n];  wmed[n] = 0;  
    wl50[n] = wREC[0][n];  wl50[n] = 0;  
    wu50[n] = wREC[0][n];  wu50[n] = 0;  
    wl90[n] = wREC[0][n];  wl90[n] = 0;  
    wu90[n] = wREC[0][n];  wu90[n] = 0;  
    wl99[n] = wREC[0][n];  wl99[n] = 0;  
    wu99[n] = wREC[0][n];  wu99[n] = 0;  

    int nentries = selected;			// number of detected events in the offsource
    int *index = new int[nentries];
    double *value = new double[nentries];
    for(int j=0;j<wREC[0][n].size();j++) {
      if(gOPT.residuals) {
        int k=0; for(int i=0;i<gEVENTS;i++) if(wNUL[i][n].size()) value[k++] = wNUL[i][n][j];  // select detected events
      } else {
        int k=0; for(int i=0;i<gEVENTS;i++) if(wREC[i][n].size()) value[k++] = wREC[i][n][j];  // select detected events
      }
      TMath::Sort(nentries,value,index,false);

      int imed = (nentries*50.)/100.; if(imed>=nentries) imed=nentries-1;
      wmed[n][j] = value[index[imed]];

      int il50 = (nentries*25.)/100.;  if(il50>=nentries) il50=nentries-1;
      int iu50 = (nentries*75.)/100.;  if(iu50>=nentries) iu50=nentries-1;
      int il90 = (nentries*5.)/100.;   if(il90>=nentries) il90=nentries-1;
      int iu90 = (nentries*95.)/100.;  if(iu90>=nentries) iu90=nentries-1;
      int il99 = (nentries*0.5)/100.;  if(il99>=nentries) il99=nentries-1;
      int iu99 = (nentries*99.5)/100.; if(iu99>=nentries) iu99=nentries-1;

      double med = wmed[n][j];
      double l50 = value[index[il50]];
      double u50 = value[index[iu50]];
      double l90 = value[index[il90]];
      double u90 = value[index[iu90]];
      double l99 = value[index[il99]];
      double u99 = value[index[iu99]];

      bool check=true;
      if(!(l50<=u50)) check=false;  
      if(!(l90<=u90)) check=false;  
      if(!(l99<=u99)) check=false;  
      if(!(med>=l50 && med<=u50)) check=false;  
      if(!(med>=l90 && med<=u90)) check=false;  
      if(!(med>=l99 && med<=u99)) check=false;  
      if(!check) {cout << "ComputeWaveformCR : standard wrong median, lower, upper bound !!! " << l50 << " " << med << " " << u50 << endl;}  
      if(!(med>=l50 && med<=u50)) med=l50;  
      wmed[n][j]=med;

      wl50[n][j] = fabs(med-l50);
      wu50[n][j] = fabs(u50-med);
      wl90[n][j] = fabs(med-l90);
      wu90[n][j] = fabs(u90-med);
      wl99[n][j] = fabs(med-l99);
      wu99[n][j] = fabs(u99-med);

      // compute log likelihood used for the the computation of bayes factor
      if(gOPT.logl_enabled && gOPT.plot_type=="spectrum") {
        int ifo_mask=1<<4*n;
        if(!(gOPT.logl_ifo_mask&ifo_mask)) continue;	// select ifos
        double freq = j/vREC[n].rate();
        if((n*wREC[0][n].size()+j)%ipc==0) {cout << pc<<"%";if (pc<100) cout << " - ";pc+=10;cout.flush();}
	if(vREC[n][j]/vREC[n].max() < gOPT.logl_arthr) continue;	// skip frequency bin if too small
        if(freq>gOPT.logl_flow && freq<gOPT.logl_fhigh && j%gOPT.logl_resample==0) {
          double xmin=value[index[0]]; 
          double xmax=value[index[nentries-1]]; 
#ifdef LOGL_TEST_KDE
          cout << j << " freq = " << freq << " vREC = " << vREC[n][j] << " vREC[n][j]/vREC[n].max() = " << vREC[n][j]/vREC[n].max() << endl;
          TCanvas* clogl = NULL;
          TH1D* hprob = NULL;
          if(j==LOGL_TEST_KDE) {
            clogl = new TCanvas();
            double dbin = (u90-l90)/(nentries/40.);
            int nbins = (xmax-xmin)/dbin; 
            hprob = new TH1D("hprob","",nbins,xmin,xmax);
            for(int i=0;i<nentries;i++) hprob->Fill(value[i]);	// fill histogram
            hprob->SetLineColor(kBlue);
            hprob->Draw();
            hprob->Scale(1./hprob->Integral(),"width" );
          }
          // smooth probability distribution with kernel density estimation
          double rho = 1.0; 					//default value
          TKDE * kde = new TKDE(nentries, value, xmin, xmax, "", rho);
          if(j==LOGL_TEST_KDE) {
            //kde->Draw("ConfidenceInterval@0.95 Same");
            kde->GetFunction()->SetLineColor(kRed);
            kde->Draw("SAME");
            clogl->SetLogx();
            clogl->SetLogy();
            clogl->Print("prob_distr_for_logl.png"); 
            exit(0);
          }
          // get the probability of the on-source reconstructed event
          double prob = kde->GetValue(vREC[n][j]); 
          delete kde;
          if(hprob!=NULL) delete hprob;
          if(clogl!=NULL) delete clogl;
#else
#ifdef LOGL_USE_KDE
          // smooth probability distribution with kernel density estimation
          double rho = 1.0; 					//default value
          TKDE * kde = new TKDE(nentries, value, xmin, xmax, "", rho);
          // get the probability of the on-source reconstructed event
/*
          double x = vREC[n][j]; 
          if(x<l90) x=l90;
          if(x>u90) x=u90;
          double prob = kde->GetValue(x); 
*/
          double prob = kde->GetValue(vREC[n][j]); 
          delete kde;
#else
          double dbin = (u90-l90)/(nentries/40.);
          int nbins = (xmax-xmin)/dbin; 
          TH1D* hprob = new TH1D("hprob","",nbins,xmin,xmax);
          for(int i=0;i<nentries;i++) hprob->Fill(value[i]);	// fill histogram
          hprob->Draw();
          hprob->Scale(1./hprob->Integral(),"width" );
          int bin = hprob->FindBin(vREC[n][j]);
          double prob = hprob->GetBinContent(bin);
          delete hprob;
#endif
#endif
	  if(prob<(1./nentries)) continue;	// skip frequency bin if prob is less than available inverse entries
          //cout << "logl -> " << j << "/" << n << " freq = " << freq 
          //     << " onsrc/l90 = " << vREC[n][j] << "/" << l90 << " logl = " << log(prob) << endl;
          nlogl++;
          logl+=log(prob);
          if(freq<fmin) fmin=freq;
          if(freq>fmax) fmax=freq;
        }
      }        
    }        
    delete [] index;
    delete [] value;

    vMED.push_back(wmed[n]);
    vL50.push_back(wl50[n]);
    vU50.push_back(wu50[n]);
    vL90.push_back(wl90[n]);
    vU90.push_back(wu90[n]);
    vL99.push_back(wl99[n]);
    vU99.push_back(wu99[n]);
  }                                                                    

  if(gOPT.logl_enabled && gOPT.plot_type=="spectrum") {
    if(pc==100) cout << pc<<"%";
    cout << endl << endl;
    // dump log likelihood
    logl = nlogl>0 ? logl/nlogl : 0;	// normalized logl
    logl*= vDOF;			// multiplied by the on-source degrees of freedom
    float l90_logl = l90_dof*logl/vDOF;
    float u90_logl = u90_dof*logl/vDOF;
    TString parm = TString::Format("log likelihood = %g [%g,%g] ( frequency bins = %d - DoF = %g [%g:%g]) fmin:fmax = %g:%g (Hz)",
                                    logl,l90_logl,u90_logl,vDOF,l90_dof,u90_dof,nlogl,fmin,fmax); 
    DumpStatistics("spectrum_logl.txt", parm, false);
    cout << parm << endl;
  }

  // Compute time boundaries which contain the P energy fraction of median

  double P = gOPT.plot_efraction;
  double estart, estop;				
  for(int n=0;n<nIFO;n++) {
    double bT, eT;
    CWB::mdc::GetTimeBoundaries(vMED[n], P, bT, eT);
    if(n==0) {
      estart=bT;estop=eT;
    } else {
      if(bT<estart) estart=bT;
      if(eT>estop)  estop=eT;
    }
  }
  cout << endl;
  cout << "estart=" << estart << "\testop=" << estop << endl;
  cout << endl;

  int nstart = (estart-vREC[0].start())*vREC[0].rate();
  int nstop  = (estop -vREC[0].start())*vREC[0].rate();

  // Compute the effective coverages

  for(int n=0;n<nIFO;n++) {
    int n50=0;
    int n90=0;
    int n99=0;
    for(int i=0;i<gEVENTS;i++) {
      if(!wREC[i][n].size()) continue;  // select detected events
      bool b50=false;
      bool b90=false;
      bool b99=false;
      for(int j=nstart;j<nstop;j++) {
        double vl50 = vMED[n][j]-fabs(vL50[n][j]);
        double vu50 = vMED[n][j]+fabs(vU50[n][j]);
        double vl90 = vMED[n][j]-fabs(vL90[n][j]);
        double vu90 = vMED[n][j]+fabs(vU90[n][j]);
        double vl99 = vMED[n][j]-fabs(vL99[n][j]);
        double vu99 = vMED[n][j]+fabs(vU99[n][j]);

        double value = wREC[i][n][j];
        if(value<vl50 || value>vu50) b50=true;
        if(value<vl90 || value>vu90) b90=true;
        if(value<vl99 || value>vu99) b99=true;
      }
      if(b50) n50++;
      if(b90) n90++;
      if(b99) n99++;
    }
    cout << endl;
    cout << "IFO " << gIFO[n] << "\tEffective Coverage at 50% = " << int(100.-100.*n50/(double)selected) << endl;
    cout << "IFO " << gIFO[n] << "\tEffective Coverage at 90% = " << int(100.-100.*n90/(double)selected) << endl;
    cout << "IFO " << gIFO[n] << "\tEffective Coverage at 99% = " << int(100.-100.*n99/(double)selected) << endl;
    cout << endl;
  }

  return selected;
}

void LoadOutputRootFiles(int nIFO, TString data_label, TString output_dir, int max_events, double ifar_thr) {

  if(ifar_thr<=0) ifar_thr=-1e-10;	// force all events which do not pass the PP cuts to be removed

  gEVENTS=0;

  char beginString[1024];
  sprintf(beginString,"wave_");	
  char endString[1024];
  sprintf(endString,".root");	
  char containString[1024];
  sprintf(containString,"%s",data_label.Data());

  vector<TString> fileList = CWB::Toolbox::getFileListFromDir(output_dir, endString, beginString, containString);

  wavearray<double>* sample_wREC[NIFO_MAX];
  for(int i=0;i<nIFO;i++) sample_wREC[i] = new wavearray<double>;
  wavearray<double>* sample_wNUL[NIFO_MAX];
  if(gOPT.residuals) for(int i=0;i<nIFO;i++) sample_wNUL[i] = new wavearray<double>;
  wavearray<double>* sample_wINJ[NIFO_MAX];
  for(int i=0;i<nIFO;i++) sample_wINJ[i] = new wavearray<double>;
  wavearray<double>* sample_wAUX[NIFO_MAX];
  if(gOPT.side=="rcut") for(int i=0;i<nIFO;i++) sample_wAUX[i] = new wavearray<double>;
  int ndim;
  double* time = new double[2*nIFO];
  float* sSNR = new float[nIFO];
  int*  size = new int[nIFO];
  float norm;
  float theta[4];
  float phi[4];
  float ifar=0;
  float likelihood;
  float factor;
  double tcoa=-1;		// geocentric tcoa
  string* log = new string;     // injection log

  TBranch* branch;
  bool ifar_exists=false;
  bool tcoa_exists=false;
  bool vnul_exists=false;
  bool vaux_exists=false;

  bool efound=false;
  char command[1024];
  int nfile = fileList.size();
  int nevt=0;
  for(int n=0;n<nfile;n++) {
    cout << n << " " << fileList[n].Data() << endl;

    TFile* froot = new TFile(fileList[n].Data(),"READ");
    if(froot==NULL) {
      cout << "LoadOutputRootFiles - Error : Failed to open file : " <<  fileList[n].Data() << endl;
      gSystem->Exit(1);
    }
    TTree* itree = (TTree*)froot->Get("waveburst");
    if(itree==NULL) {
      cout << "LoadOutputRootFiles - Error : Failed to open tree waveburst from file : " <<  fileList[n].Data() << endl;
      //continue;
      gSystem->Exit(1);
    }

    ifar_exists=false;
    TIter ifar_next(itree->GetListOfBranches());
    while ((branch=(TBranch*)ifar_next())) {
      if(TString("ifar").CompareTo(branch->GetName())==0) ifar_exists=true;
    }
    ifar_next.Reset();

    tcoa_exists=false;
    TIter tcoa_next(itree->GetListOfBranches());
    while ((branch=(TBranch*)tcoa_next())) {
      if(TString("wf_tcoa").CompareTo(branch->GetName())==0) tcoa_exists=true;
    }
    tcoa_next.Reset();

    vnul_exists=false;
    vaux_exists=false;
    TIter vnul_next(itree->GetListOfBranches());
    while ((branch=(TBranch*)vnul_next())) {
      if(TString("wNUL_0").CompareTo(branch->GetName())==0) vnul_exists=true;
      if(TString("wAUX_0").CompareTo(branch->GetName())==0) vaux_exists=true;
    }
    if(!vnul_exists) gOPT.residuals=false;	// disable residuals
    vnul_next.Reset();
    if(gOPT.side=="rcut" && !vaux_exists) {
      cout << "LoadOutputRootFiles - Error : auxiliary injection do not exist!" << endl;
      cout << "                      with option side=rcut must be provided the auxiliary injections" << endl;
      gSystem->Exit(1);
    }

    for(int k=0;k<itree->GetEntries();k++) {

      if(k%100==0) cout << "ENTRY = " << k << " / " << itree->GetEntries() << endl;
      if(!(nevt<max_events || fileList[n].Contains("job0.root"))) continue;

      // check if user input nIFO is consistent with the ndim read from the root file  
      itree->SetBranchAddress("ndim",&ndim);
      itree->GetEntry(k);
      if(ndim!=nIFO) {
        cout << "LoadOutputRootFiles - Error : number of ifos=" << ndim 
             << " declared in the output root files is different from the user declared value " << nIFO << endl;
        cout << endl << "Use NIFO=" << ndim << " or check file; " << endl << endl << fileList[n] << endl << endl; 
        gSystem->Exit(1);
      }

      for(int i=0;i<nIFO;i++) itree->SetBranchAddress(TString::Format("wREC_%d",i).Data(),&sample_wREC[i]);
      if(gOPT.residuals) {
        for(int i=0;i<nIFO;i++) itree->SetBranchAddress(TString::Format("wNUL_%d",i).Data(),&sample_wNUL[i]);
      }
      for(int i=0;i<nIFO;i++) itree->SetBranchAddress(TString::Format("wINJ_%d",i).Data(),&sample_wINJ[i]);
      if(gOPT.side=="rcut") {
        for(int i=0;i<nIFO;i++) itree->SetBranchAddress(TString::Format("wAUX_%d",i).Data(),&sample_wAUX[i]);
      }
      itree->SetBranchAddress("time",time);
      itree->SetBranchAddress("sSNR",sSNR);
      itree->SetBranchAddress("size",size);
      itree->SetBranchAddress("norm",&norm);
      itree->SetBranchAddress("theta",theta);
      itree->SetBranchAddress("phi",phi);
      itree->SetBranchAddress("likelihood",&likelihood);
      itree->SetBranchAddress("factor",&factor);
      itree->SetBranchAddress("log",&log);
      if(ifar_exists) itree->SetBranchAddress("ifar",&ifar);
      if(tcoa_exists) itree->SetBranchAddress("wf_tcoa",&tcoa);
      else if(gOPT.side!="full") {
        cout << endl << "LoadOutputRootFiles - Error : wf_tcoa not found, side!=full requires tcoa" << endl << endl; 
        gSystem->Exit(1);
      }

      itree->GetEntry(k);						// read wREC,skyprob objects

      ifar/=(24.*3600.*365.);						// sec -> year

      bool onsource=false;
      if(fabs(gOPT.gps_event-time[0])<gOPT.gps_error) onsource=true;
      if(!onsource && ifar<ifar_thr) continue;				// skip offsource events below the ifar threshold

      bool fail=false;

      // check if objects are not null
      for(int i=0;i<nIFO;i++) {
        if(sample_wREC[i]==NULL) {
          cout << "LoadOutputRootFiles - Error : Object wavearray not exist !!! " <<  endl;
          //gSystem->Exit(1);
          fail=true;
        }
      }
      if(fail) continue; 

      tcoa+=gOPT.tcoa_shift;						// apply time shift to tcoa

      // fill object vectors
      cout.precision(14);
      if(fabs(gOPT.gps_event-time[0])<gOPT.gps_error && fileList[n].Contains("_job0.root")) {	// OnSource
        vTCOA=tcoa;						// store OnSource sample geoacentric tcoa
        vFACTOR=factor;						// store OnSource factor
        vTHETA=theta[1];					// store OnSource sample theta
        vPHI=phi[1];						// store OnSource sample phi
        vDOF= norm>0?size[0]/norm:0.;				// store OnSource degrees of freedom
        TString m1 = GetParameterFromInjLog(log->c_str(), "mass1");	// get mass1 from injection parameters
        TString m2 = GetParameterFromInjLog(log->c_str(), "mass2");	// get mass2 from injection parameters
        vM1=m1.Atof();						// store OnSource sample m1
        vM2=m2.Atof();						// store OnSource sample m2
        TString s1 = GetParameterFromInjLog(log->c_str(), "spin1");	// get spin1z from injection parameters
        TString s2 = GetParameterFromInjLog(log->c_str(), "spin2");	// get spin2z from injection parameters
        vS1=s1.Atof();						// store OnSource sample s1
        vS2=s2.Atof();						// store OnSource sample s2
        for(int i=0;i<nIFO;i++) vSNR[i]=sqrt(sSNR[i]);		// store OnSource SNR 
        efound=true;
        vTREC.resize(nIFO);
        for(int i=0;i<nIFO;i++) {
          if(gOPT.side=="rcut") {
            CWB::mdc::Align(*sample_wAUX[i], *sample_wREC[i]);	// align onsource aux wrt rec
            sample_wAUX[i]->start(fmod(sample_wAUX[i]->start(),gOPT.inj_time_step));
            vAUX.push_back(*sample_wAUX[i]);  			// store reconstructed null event waveform
          }
          CWB::mdc::Align(*sample_wINJ[i], *sample_wREC[i]);	// align onsource inj wrt rec
          sample_wINJ[i]->start(fmod(sample_wINJ[i]->start(),gOPT.inj_time_step));
          vINJ.push_back(*sample_wINJ[i]);  			// store whitened injected waveform
          vTREC[i]=sample_wREC[i]->start();	// save absolute start time
          sample_wREC[i]->start(fmod(sample_wREC[i]->start(),gOPT.inj_time_step));
          vREC.push_back(*sample_wREC[i]);  			// store reconstructed event waveform
          if(gOPT.residuals) {
            sample_wNUL[i]->start(fmod(sample_wNUL[i]->start(),gOPT.inj_time_step));
            vNUL.push_back(*sample_wNUL[i]);  			// store reconstructed null event waveform
          }
        }
        double toff = gOPT.inj_time_step-int(fmod(gOPT.gps_event,gOPT.inj_time_step))-gOPT.off_event;
        cout << "EVENT FOUND @ GPS time -> " << time[0] << "\tSNR -> " << sqrt(likelihood) << "\tIFAR -> " << ifar << " yrs" << " toff -> " << toff << endl;
        cout << endl;
        cout << "                  File -> " << fileList[n] << endl << endl;
        for(int i=0;i<nIFO;i++) vREC[i].start(int(vREC[i].start()+toff)%gOPT.inj_time_step);
        for(int i=0;i<nIFO;i++) vINJ[i].start(int(vINJ[i].start()+toff)%gOPT.inj_time_step);
        if(gOPT.residuals) {
          for(int i=0;i<nIFO;i++) vNUL[i].start(int(vNUL[i].start()+toff)%gOPT.inj_time_step);
        }
        if(gOPT.side=="rcut") {
          for(int i=0;i<nIFO;i++) vAUX[i].start(int(vAUX[i].start()+toff)%gOPT.inj_time_step);
        }

        // read ifo names from the tree user list and updated gIFO global array according to the input data network setup
        TList* ifoList = itree->GetUserInfo();
        detectorParams dParams[NIFO_MAX];
        for (int i=0;i<ifoList->GetSize();i++) {
          detector* pDetector = (detector*)ifoList->At(i);
          dParams[i] = pDetector->getDetectorParams();
          gIFO[i] = dParams[i].name;
        }

      } else {							// OffSource
        wTCOA[gEVENTS]=tcoa;					// store injected geocentric tcoa
        wFACTOR[gEVENTS]=factor;				// store injected factor
        wTHETA[gEVENTS]=theta[1];				// store injected theta
        wPHI[gEVENTS]=phi[1];					// store injected phi
        wDOF[gEVENTS]= norm>0?size[0]/norm:0.;			// store injected degrees of freedom
        TString m1 = GetParameterFromInjLog(log->c_str(), "mass1");
        TString m2 = GetParameterFromInjLog(log->c_str(), "mass2");
        wM1[gEVENTS]=m1.Atof();					// store injected m1
        wM2[gEVENTS]=m2.Atof();					// store injected m2
        TString s1 = GetParameterFromInjLog(log->c_str(), "spin1");
        TString s2 = GetParameterFromInjLog(log->c_str(), "spin2");
        wS1[gEVENTS]=s1.Atof();					// store injected s1
        wS2[gEVENTS]=s2.Atof();					// store injected s2
        for(int i=0;i<nIFO;i++) wSNR[gEVENTS][i]=sqrt(sSNR[i]);	// store reconstructed SNR 
        for(int i=0;i<nIFO;i++) {
          double toff=0;
          if(gOPT.inj_onsrc) {
            toff = gOPT.inj_time_step-int(fmod(gOPT.gps_event,gOPT.inj_time_step))-gOPT.off_event;
          }
          sample_wREC[i]->start(fmod(sample_wREC[i]->start()+toff,gOPT.inj_time_step));
          wREC[gEVENTS].push_back(*sample_wREC[i]);  		// store reconstructed posterior sample waveform
          sample_wINJ[i]->start(fmod(sample_wINJ[i]->start()+toff,gOPT.inj_time_step));
          wINJ[gEVENTS].push_back(*sample_wINJ[i]);  		// store whitened posterior sample waveform
          if(gOPT.residuals) {
            sample_wNUL[i]->start(fmod(sample_wNUL[i]->start()+toff,gOPT.inj_time_step));
            wNUL[gEVENTS].push_back(*sample_wNUL[i]);  		// store reconstructed null posterior sample waveform
          }
          if(gOPT.side=="rcut") {
            sample_wAUX[i]->start(fmod(sample_wAUX[i]->start()+toff,gOPT.inj_time_step));
            wAUX[gEVENTS].push_back(*sample_wAUX[i]);  		// store reconstructed null posterior sample waveform
          }
        }
        gEVENTS++;
      }

      if(gEVENTS>=PE_MAX_EVENTS) {
        cout << "LoadOutputRootFiles - Warning : Number of events > PE_MAX_EVENTS !!! " <<  PE_MAX_EVENTS << endl;
        break;
      }

      nevt++;
    }

    froot->Close();
  }

  if(vREC.size()!=nIFO || vINJ.size()!=nIFO) {
    cout << "LoadOutputRootFiles - Error : Number of OnSource Events = " << vREC.size()/nIFO << " -> must be 1" << endl;
    cout << "                              check the onsrc links in the ofsrc merge directory " << endl;
    cout << "                              check the onsrc event time " << endl;
    exit(1);
  }

  if(gEVENTS==0) {
    cout << "LoadOutputRootFiles - Error : Found 0 OffSource Events !!! " << endl;
    exit(1);
  }


  if(!efound) {
    cout << "LoadOutputRootFiles - Error : GPS EVENT NOT FOUND !!! " <<  endl;
    cout << endl;
    cout << "Possible reasons:" << endl;
    cout << endl;
    cout << "1) The root file with EVENT has not been included" << endl; 
    cout << "2) gps_error is too rectrictive, current value is " << gOPT.gps_error << " , try a new value with option --gps_error=... " << endl;
    cout << endl;
    gSystem->Exit(1);
  }

  double toffset = wREC[0][0].start()-vREC[0].start();
  if(toffset!=0) {
    cout << endl;
    cout << "LoadOutputRootFiles - Warning  : detected event and posterios have a possible inconsistent GPS time" << endl;
    cout << "                                 check the detected event root file" << endl;
    cout << "                                 try --off_event=" << gOPT.off_event-toffset << endl;
    cout << endl;
//    gSystem->Exit(1);
  }

  for(int i=0;i<nIFO;i++) delete sample_wREC[i];
  if(gOPT.residuals) for(int i=0;i<nIFO;i++) delete sample_wNUL[i];
  if(gOPT.side=="rcut") for(int i=0;i<nIFO;i++) delete sample_wAUX[i];
  if(gOPT.statistics) for(int i=0;i<nIFO;i++) delete sample_wINJ[i];
  delete [] time;
  delete log;
}

void PlotWaveformsCR(int nIFO) {

  gnetwork* NET = new gnetwork(nIFO,gIFO);

  TString ifar_thr_label  = TString::Format("%g",gOPT.ifar_thr).ReplaceAll(".","d");
  TString ncycles_label   = TString::Format("%d",gOPT.ncycles);

  TString type_label = "";
  type_label += "ifar_"+ifar_thr_label;
  if(gOPT.plot_type=="phase")      type_label+="_Phase";
  if(gOPT.plot_type=="time")       type_label+="_Time";
  if(gOPT.plot_type=="envelope")   type_label+="_Envelope";
  if(gOPT.plot_type=="frequency")  type_label+="_Frequency";
  if(gOPT.plot_type=="spectrum")   type_label+="_Spectrum";

  TString pdir=gOPT.plot_odir;
  char title[256]; char ofname[256];
  for(int n=0; n<nIFO; n++) {
    if(gOPT.plot_99perc) {
      sprintf(title,"%s (ifar=%sy) : cWB REC (red) - Errors (entries=%d) : med (black) : 50% (dark_grey) : 90% (medium gray) : 99% (light_gray)",
              gIFO[n].Data(),ifar_thr_label.Data(),gEVENTS);
    } else {
      sprintf(title,"%s (ifar=%sy) : cWB REC (red) - Errors (entries=%d) : med (black) : 50% (dark_grey) : 90% (light gray)",
              gIFO[n].Data(),ifar_thr_label.Data(),gEVENTS);
    }

    int K = (gOPT.plot_type=="spectrum" || gOPT.plot_type=="phase") ? 1 : 3;
    for(int k=0;k<K;k++) {	// plot at 3 different time zooms
      // get tcoa of vINJ[n] after time/phase sync
      double tcoa = GetInjTcoa(vTCOA, NET, gIFO[n], vTHETA, vPHI);
             tcoa = tcoa-vTREC[n]+vREC[n].start();
      if(gOPT.sync_phase || gOPT.sync_time) {
        tcoa = GetSyncTcoa(vINJ[n], vTSYNC[n], vPSYNC[n], tcoa);
      }
      // get time of max abs(amp) waveform envelope
      double tmax=tcoa; 
      if(gOPT.plot_type=="time" || gOPT.plot_type=="envelope" || gOPT.plot_type=="frequency") {
        wavearray<double> wenv = CWB::mdc::GetEnvelope(&vMED[n]);
        double amax=0;
        for(int i=0;i<wenv.size();i++) if(fabs(wenv[i])>amax) {amax=fabs(wenv[i]);tmax=wenv.start()+i/wenv.rate();}
      }
      double tref = (tcoa>tmax) ? tcoa : tmax;	// time used as reference to zoom the waveform
      double P=1; 
      TString zoom;	// the A/B/C characters are use to sort the plots when cwb_mkhtml is used  
      if(k==0) {P=0.65/gOPT.plot_zoom; if(P>1) P=1; zoom = TString::Format("C%d",int(P*100));}
      if(k==1) {P=0.85/gOPT.plot_zoom; if(P>1) P=1; zoom = TString::Format("B%d",int(P*100));}
      if(k==2) {P=0.99/gOPT.plot_zoom; if(P>1) P=1; zoom = TString::Format("A%d",int(P*100));}
      double Q=P+0.5*(1-P);	// increases the time range on the right side
      if(gOPT.plot_type=="spectrum" || gOPT.plot_type=="phase") zoom="A";
      sprintf(ofname,"%s_%s_%s_%s.%s",gIFO[n].Data(),gOPT.plot_fname.Data(), type_label.Data(), zoom.Data(), "png");
      if(gOPT.residuals) {
        PlotWaveformAsymmErrors(ofname, title, vNUL[n], vMED[n], vL50[n], vU50[n], vL90[n], vU90[n], vL99[n], vU99[n], vNUL[n], pdir, P, Q, tref, tcoa, n);
      } else {
        PlotWaveformAsymmErrors(ofname, title, vREC[n], vMED[n], vL50[n], vU50[n], vL90[n], vU90[n], vL99[n], vU99[n], vREC[n], pdir, P, Q, tref, tcoa, n);
      }
    }
  }
}

void PlotWaveformAsymmErrors(TString ofname, TString title, wavearray<double> wrec,
                             wavearray<double> wmed, wavearray<double> wl50, wavearray<double> wu50,
                             wavearray<double> wl90, wavearray<double> wu90, 
                             wavearray<double> wl99, wavearray<double> wu99, wavearray<double> wref, 
                             TString pdir, double P, double Q, double tref, double tcoa, int ifoid) {

  // resample data to increase the plot resolution 
  if(gOPT.plot_type=="time") {
    wrec.Resample(16384);
    wmed.Resample(16384);
    wl50.Resample(16384);
    wu50.Resample(16384);
    wl90.Resample(16384);
    wu90.Resample(16384);
    wl99.Resample(16384);
    wu99.Resample(16384);
    wref.Resample(16384);
  }

  int size = wrec.size();

  wavearray<double> x(size);
  wavearray<double> ex(size); ex=0;
  if(gOPT.plot_type=="spectrum") {
    for (int i=0; i<size; i++) x[i] = i/wrec.rate();
  } else if(gOPT.plot_type=="phase") {
    for (int i=0; i<size; i++) x[i] = i;
  } else {
    for (int i=0; i<size; i++) x[i] = i/wrec.rate()+(wref.start());
  }

  TString xtitle;
  if(gOPT.plot_type=="phase")     xtitle = "Phase (cycles)";
  if(gOPT.plot_type=="envelope")  xtitle = TString::Format("Time Envelope (sec) : GPS OFFSET = %d",int(gOPT.gps_event-gOPT.inj_time_step));
  if(gOPT.plot_type=="time")      xtitle = TString::Format("Time (sec) : GPS OFFSET = %d",int(gOPT.gps_event-gOPT.inj_time_step));
  if(gOPT.plot_type=="frequency") xtitle = TString::Format("Time Frequency (sec) : GPS OFFSET = %d",int(gOPT.gps_event-gOPT.inj_time_step));
  if(gOPT.plot_type=="spectrum")  xtitle = "Frequency (Hz)";

  TGraphAsymmErrors* egr99 = new TGraphAsymmErrors(size,x.data,wmed.data,ex.data,ex.data,wl99.data,wu99.data);
  egr99->SetLineColor(18);
  egr99->SetFillStyle(1001);
  egr99->SetFillColor(18);
  egr99->GetXaxis()->SetTitle(xtitle);
  if(gOPT.plot_type=="phase")     egr99->GetYaxis()->SetTitle("degrees");
  if(gOPT.plot_type=="envelope")  egr99->GetYaxis()->SetTitle("magnitude");
  if(gOPT.plot_type=="time")      egr99->GetYaxis()->SetTitle("magnitude");
  if(gOPT.plot_type=="frequency") egr99->GetYaxis()->SetTitle("Hz");
  if(gOPT.plot_type=="spectrum")  egr99->GetYaxis()->SetTitle("magnitude");
  egr99->SetTitle(title);
  egr99->GetXaxis()->SetTitleFont(42);
  egr99->GetXaxis()->SetLabelFont(42);
  egr99->GetXaxis()->SetLabelOffset(0.012);
  egr99->GetXaxis()->SetTitleOffset(1.5);
  egr99->GetYaxis()->SetTitleFont(42);
  egr99->GetYaxis()->SetLabelFont(42);
  egr99->GetYaxis()->SetLabelOffset(0.01);
  egr99->GetYaxis()->SetTitleOffset(1.4);

  TGraphAsymmErrors* egr90 = new TGraphAsymmErrors(size,x.data,wmed.data,ex.data,ex.data,wl90.data,wu90.data);
  egr90->SetLineColor(17);
  egr90->SetFillStyle(1001);
  egr90->SetFillColor(17);
  egr90->GetXaxis()->SetTitle(xtitle);
  if(gOPT.plot_type=="phase")     egr90->GetYaxis()->SetTitle("degrees");
  if(gOPT.plot_type=="envelope")  egr90->GetYaxis()->SetTitle("magnitude");
  if(gOPT.plot_type=="time")      egr90->GetYaxis()->SetTitle("magnitude");
  if(gOPT.plot_type=="frequency") egr90->GetYaxis()->SetTitle("Hz");
  if(gOPT.plot_type=="spectrum")  egr90->GetYaxis()->SetTitle("magnitude");
  egr90->SetTitle(title);
  egr90->GetXaxis()->SetTitleFont(42);
  egr90->GetXaxis()->SetLabelFont(42);
  egr90->GetXaxis()->SetLabelOffset(0.012);
  egr90->GetXaxis()->SetTitleOffset(1.5);
  egr90->GetYaxis()->SetTitleFont(42);
  egr90->GetYaxis()->SetLabelFont(42);
  egr90->GetYaxis()->SetLabelOffset(0.01);
  egr90->GetYaxis()->SetTitleOffset(1.4);

  TGraphAsymmErrors* egr50 = new TGraphAsymmErrors(size,x.data,wmed.data,ex.data,ex.data,wl50.data,wu50.data);
  egr50->SetLineColor(15);
  egr50->SetFillStyle(1001);
  egr50->SetFillColor(15);

  TGraph* agr = new TGraph(size,x.data,wmed.data);
  agr->SetLineWidth(1);
  agr->SetLineColor(kBlack);
  agr->SetLineStyle(1);

  TGraph* gr = new TGraph(size,x.data,wrec.data);
  gr->SetLineWidth(1);
  gr->SetLineColor(2);

  TCanvas* canvas = new TCanvas("wave_pe_cwb","wave_pe_cwb",200,20,800,500);
  canvas->cd();
  canvas->SetGridx();
  canvas->SetGridy();
  if(gOPT.plot_type=="spectrum") {
    canvas->SetLogx();
    canvas->SetLogy();
  }
  
  if(gOPT.plot_type=="time" || gOPT.plot_type=="envelope" || gOPT.plot_type=="frequency") {
    double bT, eT;
    wavearray<double> W=wmed; W+=wu99; W+=wref;	// the 99% belt + wref is used to compute the boundaries
    CWB::mdc::GetTimeBoundaries(W, P, bT, eT, tref, Q);
    if(gOPT.plot_99perc) egr99->GetXaxis()->SetRangeUser(bT, eT);
    else                 egr90->GetXaxis()->SetRangeUser(bT, eT);
  }
  if(gOPT.plot_type=="spectrum") {
    if(gOPT.plot_99perc) egr99->GetYaxis()->SetRangeUser(1e-3, 1.5 * egr99->GetHistogram()->GetMaximum());
    else                 egr90->GetYaxis()->SetRangeUser(1e-3, 1.5 * egr90->GetHistogram()->GetMaximum());
    if(gOPT.plot_99perc) egr99->GetXaxis()->SetRangeUser(20, wrec.size()/wrec.rate());
    else                 egr90->GetXaxis()->SetRangeUser(20, wrec.size()/wrec.rate());
  }

  if(gOPT.plot_99perc) {
    egr99->Draw("A3");
    egr90->Draw("3same");
  } else {
    egr90->Draw("A3");
  }
  egr50->Draw("3same");
  agr->Draw("Lsame");
  gr->Draw("Lsame");

  // add vertical green line at tcoa time
  TGraph* grline = NULL;
  if(gOPT.plot_type=="time" || gOPT.plot_type=="envelope" || gOPT.plot_type=="frequency") {
    wavearray<double> wline = wrec;
    double dt=1./wline.rate();
    bool bline=false;
    for(int i=0;i<wline.size();i++) {
      double time=i*dt+wline.start();
      if(fabs(time-tcoa)<dt && !bline) {
        if(gOPT.plot_type=="time")      wline[i]=-1.0;
        if(gOPT.plot_type=="envelope")  wline[i]= 0.0;
        if(gOPT.plot_type=="frequency") wline[i]= 0.0;
        bline=true;
      } else wline[i]=-1000.;
    }
    grline = new TGraph(size,x.data,wline.data);
    grline->SetLineWidth(1);
    grline->SetLineColor(kGreen);
    grline->Draw("Lsame");
  }

  TGraph* grpfrq = NULL;
  TGraph* grpsnr = NULL;
  if(gOPT.plot_type=="phase") {

    grpsnr = new TGraph(x.size(),x.data,pSNR[ifoid].data);
    grpsnr->SetLineWidth(1);
    grpsnr->SetLineColor(kBlue);
    grpsnr->Draw("Lsame");

    grpfrq = new TGraph(x.size(),x.data,pFRQ[ifoid].data);
    grpfrq->SetLineWidth(1);
    grpfrq->SetLineColor(kGreen);
    grpfrq->Draw("Lsame");

    if(gOPT.plot_99perc) egr99->GetXaxis()->SetRangeUser(0, gOPT.ncycles);
    else                 egr90->GetXaxis()->SetRangeUser(0, gOPT.ncycles);
    if(gOPT.plot_99perc) egr99->GetYaxis()->SetRangeUser(-180, 180);
    else                 egr90->GetYaxis()->SetRangeUser(-180, 180);
  }

  if(ofname!="") {
    TString pfname = pdir!="" ? TString(pdir)+TString("/")+ofname : ofname;
    canvas->Print(pfname);
    cout << "write : " << pfname << endl;
    pfname.ReplaceAll(".png",".root");
    canvas->Print(pfname);
    cout << "write : " << pfname << endl;
  }

  delete canvas;
  delete egr50;
  delete egr90;
  if(gOPT.plot_99perc) delete egr99;
  delete agr;
  delete gr;
  if(grline!=NULL) delete grline;
  if(grpfrq!=NULL) delete grpfrq;
  if(grpsnr!=NULL) delete grpsnr;
}

void ReadUserOptions(TString options) {

  // get plugin options 

  if(TString(options)!="") {

    //cout << "WF options : " << options << endl;
    TObjArray* token = TString(options).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++) {

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();
      stok.ReplaceAll(" ",""); 	// remove white spaces
      stok.ReplaceAll("	","");	// remove tabs

      if(stok.Contains("dump_cr=")) {
        TString pe_dump_cr=stok;
        pe_dump_cr.Remove(0,pe_dump_cr.Last('=')+1);
        if(pe_dump_cr=="true")  gOPT.dump_cr=true;
        if(pe_dump_cr=="false") gOPT.dump_cr=false;
      }

      if(stok.Contains("use_onsrc=")) {
        TString pe_use_onsrc=stok;
        pe_use_onsrc.Remove(0,pe_use_onsrc.Last('=')+1);
        if(pe_use_onsrc!="mfull" && pe_use_onsrc!="mleft" && pe_use_onsrc!="efull" && pe_use_onsrc!="eleft") pe_use_onsrc="";
        gOPT.use_onsrc=pe_use_onsrc;
      }

      if(stok.Contains("sync_phase=")) {
        TString pe_sync_phase=stok;
        pe_sync_phase.Remove(0,pe_sync_phase.Last('=')+1);
        if(pe_sync_phase=="true")  gOPT.sync_phase=true;
        if(pe_sync_phase=="false") gOPT.sync_phase=false;
      }

      if(stok.Contains("sync_time=")) {
        TString pe_sync_time=stok;
        pe_sync_time.Remove(0,pe_sync_time.Last('=')+1);
        if(pe_sync_time=="true")  gOPT.sync_time=true;
        if(pe_sync_time=="false") gOPT.sync_time=false;
      }

      if(stok.Contains("sensitivity=")) {
        TString pe_sensitivity=stok;
        pe_sensitivity.Remove(0,pe_sensitivity.Last('=')+1);
        if(pe_sensitivity=="true")  gOPT.sensitivity=true;
        if(pe_sensitivity=="false") gOPT.sensitivity=false;
      }

      if(stok.Contains("residuals=")) {
        TString pe_residuals=stok;
        pe_residuals.Remove(0,pe_residuals.Last('=')+1);
        if(pe_residuals=="true")  gOPT.residuals=true;
        if(pe_residuals=="false") gOPT.residuals=false;
      }

      if(stok.Contains("statistics=")) {
        TString pe_statistics=stok;
        pe_statistics.Remove(0,pe_statistics.Last('=')+1);
        if(pe_statistics=="true")  gOPT.statistics=true;
        if(pe_statistics=="false") gOPT.statistics=false;
      }

      if(stok.Contains("plot_cr=")) {
        TString pe_plot_cr=stok;
        pe_plot_cr.Remove(0,pe_plot_cr.Last('=')+1);
        if(pe_plot_cr=="true")  gOPT.plot_cr=true;
        if(pe_plot_cr=="false") gOPT.plot_cr=false;
      }

      if(stok.Contains("plot_99perc=")) {
        TString pe_plot_99perc=stok;
        pe_plot_99perc.Remove(0,pe_plot_99perc.Last('=')+1);
        if(pe_plot_99perc=="true")  gOPT.plot_99perc=true;
        if(pe_plot_99perc=="false") gOPT.plot_99perc=false;
      }

      if(stok.Contains("make_ff=")) {
        TString pe_make_ff=stok;
        pe_make_ff.Remove(0,pe_make_ff.Last('=')+1);
        if(pe_make_ff=="true")  gOPT.make_ff=true;
        if(pe_make_ff=="false") gOPT.make_ff=false;
      }

      if(stok.Contains("efraction=")) {
        TString pe_efraction=stok;
        pe_efraction.Remove(0,pe_efraction.Last('=')+1);
        if(pe_efraction.IsFloat()) gOPT.efraction=pe_efraction.Atof();
      }

      if(stok.Contains("side=")) {
        TString pe_side=stok;
        pe_side.Remove(0,pe_side.Last('=')+1);
        gOPT.side=pe_side;
      }

      if(stok.Contains("coverage=")) {
        TString pe_coverage=stok;
        pe_coverage.Remove(0,pe_coverage.Last('=')+1);
        gOPT.coverage=pe_coverage;
      }

      if(stok.Contains("make_of=")) {
        TString pe_make_of=stok;
        pe_make_of.Remove(0,pe_make_of.Last('=')+1);
        if(pe_make_of=="true")  gOPT.make_of=true;
        if(pe_make_of=="false") gOPT.make_of=false;
      }

      if(stok.Contains("make_re=")) {
        TString pe_make_re=stok;
        pe_make_re.Remove(0,pe_make_re.Last('=')+1);
        if(pe_make_re=="true")  gOPT.make_re=true;
        if(pe_make_re=="false") gOPT.make_re=false;
      }

      if(stok.Contains("make_ne=")) {
        TString pe_make_ne=stok;
        pe_make_ne.Remove(0,pe_make_ne.Last('=')+1);
        if(pe_make_ne=="true")  gOPT.make_ne=true;
        if(pe_make_ne=="false") gOPT.make_ne=false;
      }

      if(stok.Contains("make_onsrc=")) {
        TString pe_make_onsrc=stok;
        pe_make_onsrc.Remove(0,pe_make_onsrc.Last('=')+1);
        if(pe_make_onsrc=="true")  gOPT.make_onsrc=true;
        if(pe_make_onsrc=="false") gOPT.make_onsrc=false;
      }

      if(stok.Contains("gw_name=")) {
        TString pe_gw_name=stok;
        pe_gw_name.Remove(0,pe_gw_name.Last('=')+1);
        gOPT.gw_name=pe_gw_name;
      }

      if(stok.Contains("data_label=")) {
        TString pe_data_label=stok;
        pe_data_label.Remove(0,pe_data_label.Last('=')+1);
        gOPT.data_label=pe_data_label;
      }

      if(stok.Contains("wave_dir=")) {
        TString pe_wave_dir=stok;
        pe_wave_dir.Remove(0,pe_wave_dir.Last('=')+1);
        gOPT.wave_dir=pe_wave_dir;
      }

      if(stok.Contains("nifo=")) {
        TString pe_nifo=stok;
        pe_nifo.Remove(0,pe_nifo.Last('=')+1);
        if(pe_nifo.IsDigit()) gOPT.nifo=pe_nifo.Atoi();
      }

      if(stok.Contains("nevents=")) {
        TString pe_nevents=stok;
        pe_nevents.Remove(0,pe_nevents.Last('=')+1);
        if(pe_nevents.IsDigit()) gOPT.nevents=pe_nevents.Atoi();
      }

      if(stok.Contains("ncycles=")) {
        TString pe_ncycles=stok;
        pe_ncycles.Remove(0,pe_ncycles.Last('=')+1);
        if(pe_ncycles.IsDigit()) gOPT.ncycles=pe_ncycles.Atoi();
      }

      if(stok.Contains("ifar_thr=")) {
        TString pe_ifar_thr=stok;
        pe_ifar_thr.Remove(0,pe_ifar_thr.Last('=')+1);
        if(pe_ifar_thr.IsFloat()) gOPT.ifar_thr=pe_ifar_thr.Atof();
      }

      if(stok.Contains("gps_event=")) {
        TString pe_gps_event=stok;
        pe_gps_event.Remove(0,pe_gps_event.Last('=')+1);
        if(pe_gps_event.IsFloat()) gOPT.gps_event=pe_gps_event.Atof();
      }

      if(stok.Contains("tcoa_shift=")) {
        TString pe_tcoa_shift=stok;
        pe_tcoa_shift.Remove(0,pe_tcoa_shift.Last('=')+1);
        if(pe_tcoa_shift.IsFloat()) gOPT.tcoa_shift=pe_tcoa_shift.Atof();
      }

      if(stok.Contains("gps_error=")) {
        TString pe_gps_error=stok;
        pe_gps_error.Remove(0,pe_gps_error.Last('=')+1);
        if(pe_gps_error.IsFloat()) gOPT.gps_error=pe_gps_error.Atof();
      }

      if(stok.Contains("off_event=")) {
        TString pe_off_event=stok;
        pe_off_event.Remove(0,pe_off_event.Last('=')+1);
        if(pe_off_event.IsFloat()) gOPT.off_event=pe_off_event.Atoi();
      }

      if(stok.Contains("inj_time_step=")) {
        TString pe_inj_time_step=stok;
        pe_inj_time_step.Remove(0,pe_inj_time_step.Last('=')+1);
        if(pe_inj_time_step.IsDigit()) {
          gOPT.inj_time_step=pe_inj_time_step.Atoi();
          if(gOPT.inj_time_step<=0) {
            gOPT.inj_time_step=PE_INJ_TIME_STEP;
            gOPT.inj_onsrc=true;
          }
        } else {
	  cout << endl;
          cout << "cwb_report - Error : input parameter inj_time_step=" << pe_inj_time_step << " is not an integer !!! " <<  endl;
	  cout << endl;
          exit(1);
        }
      }

      if(stok.Contains("ccut=")) {
        TString pe_ccut=stok;
        pe_ccut.Remove(0,pe_ccut.Last('=')+1);
        gOPT.ccut=pe_ccut;
        GetCCutParms(gOPT.ccut);
      }

      if(stok.Contains("rcut=")) {
        TString pe_rcut=stok;
        pe_rcut.Remove(0,pe_rcut.Last('=')+1);
        gOPT.rcut=pe_rcut;
        GetRCutParms(gOPT.rcut);
      }

      if(stok.Contains("logl=")) {
        TString pe_logl=stok;
        pe_logl.Remove(0,pe_logl.Last('=')+1);
        gOPT.logl=pe_logl;
        GetLoglParms(gOPT.logl);
      }

      if(stok.Contains("rstat=")) {
        TString pe_rstat=stok;
        pe_rstat.Remove(0,pe_rstat.Last('=')+1);
        gOPT.rstat=pe_rstat;
        GetRstatParms(gOPT.rstat);
      }

      if(stok.Contains("plot_fname=")) {
        TString pe_plot_fname=stok;
        pe_plot_fname.Remove(0,pe_plot_fname.Last('=')+1);
        gOPT.plot_fname=pe_plot_fname;
      }

      if(stok.Contains("plot_odir=")) {
        TString pe_plot_odir=stok;
        pe_plot_odir.Remove(0,pe_plot_odir.Last('=')+1);
        gOPT.plot_odir=pe_plot_odir;
      }

      if(stok.Contains("plot_type=")) {
        TString pe_plot_type=stok;
        pe_plot_type.Remove(0,pe_plot_type.Last('=')+1);
        gOPT.plot_type=pe_plot_type;
      }

      if(stok.Contains("plot_efraction=")) {
        TString pe_plot_efraction=stok;
        pe_plot_efraction.Remove(0,pe_plot_efraction.Last('=')+1);
        if(pe_plot_efraction.IsFloat()) gOPT.plot_efraction=pe_plot_efraction.Atof();
      }

      if(stok.Contains("plot_zoom=")) {
        TString pe_plot_zoom=stok;
        pe_plot_zoom.Remove(0,pe_plot_zoom.Last('=')+1);
        if(pe_plot_zoom.IsFloat()) gOPT.plot_zoom=pe_plot_zoom.Atof();
        if(gOPT.plot_zoom<=0) gOPT.plot_zoom=1.;
      }
    }
  }
}

void ResetUserOptions() {

    gOPT.gw_name	= PE_GW_NAME;  
    gOPT.data_label	= PE_DATA_LABEL;  
    gOPT.wave_dir	= PE_WAVE_DIR;  
    gOPT.nifo		= PE_NIFO;
    gOPT.nevents	= PE_MAX_EVENTS;
    gOPT.ncycles	= PE_NCYCLES;
    gOPT.ifar_thr	= PE_IFAR_THR;
    gOPT.gps_event	= PE_GPS_EVENT;
    gOPT.tcoa_shift	= PE_TCOA_SHIFT;
    gOPT.gps_error	= PE_GPS_ERROR;
    gOPT.off_event	= PE_OFF_EVENT;
    gOPT.inj_time_step	= PE_INJ_TIME_STEP;
    gOPT.inj_onsrc	= PE_INJ_ONSRC;
    gOPT.ccut		= PE_CCUT;
    gOPT.rcut     	= PE_RCUT;
    gOPT.logl     	= PE_LOGL;
    gOPT.rstat     	= PE_RSTAT;
    gOPT.dump_cr	= PE_DUMP_CR;
    gOPT.use_onsrc      = PE_USE_ONSRC;
    gOPT.sync_phase	= PE_SYNC_PHASE;
    gOPT.sync_time	= PE_SYNC_TIME;
    gOPT.sensitivity	= PE_SENSITIVITY;
    gOPT.residuals	= PE_RESIDUALS;
    gOPT.statistics	= PE_STATISTICS;
    gOPT.efraction	= PE_EFRACTION;
    gOPT.side	        = PE_SIDE;
    gOPT.coverage       = PE_COVERAGE;
    gOPT.plot_fname	= PE_PLOT_FNAME;
    gOPT.plot_odir	= PE_PLOT_ODIR;
    gOPT.plot_type	= PE_PLOT_TYPE;
    gOPT.plot_efraction	= PE_PLOT_EFRACTION;
    gOPT.plot_zoom	= PE_PLOT_ZOOM;
    gOPT.plot_cr	= PE_PLOT_CR;
    gOPT.plot_99perc	= PE_PLOT_99PERC;
    gOPT.make_ff	= PE_MAKE_FF;
    gOPT.make_of	= PE_MAKE_OF;
    gOPT.make_re	= PE_MAKE_RE;
    gOPT.make_ne	= PE_MAKE_NE;
    gOPT.make_onsrc	= PE_MAKE_ONSRC;
}

void PrintUserOptions(TString ofname) {

    // redirect cout to buffer and return contents
    std::streambuf *old = NULL;
    std::stringstream buffer;
    if(ofname!="") {
      old = std::cout.rdbuf(buffer.rdbuf());
    }

    cout.precision(6);

    TString sync_phase    = gOPT.sync_phase  ?"true":"false";
    TString sync_time     = gOPT.sync_time   ?"true":"false";
    TString sensitivity   = gOPT.sensitivity ?"true":"false";
    TString residuals     = gOPT.residuals   ?"true":"false";
    TString statistics    = gOPT.statistics  ?"true":"false";
    TString plot_cr       = gOPT.plot_cr     ?"true":"false";
    TString plot_99perc   = gOPT.plot_99perc ?"true":"false";
    TString make_ff       = gOPT.make_ff     ?"true":"false";
    TString make_of       = gOPT.make_of     ?"true":"false";
    TString make_re       = gOPT.make_re     ?"true":"false";
    TString make_ne       = gOPT.make_ne     ?"true":"false";
    TString make_onsrc    = gOPT.make_onsrc  ?"true":"false";
    TString logl_enabled  = gOPT.logl_enabled  ?"true":"false";
    TString rstat_enabled = gOPT.rstat_enabled ?"true":"false";

    cout << endl;

    cout << "-----------------------------------------"     << endl;
    cout << "cwb_pereport config options              "     << endl;
    cout << endl; 
    cout << "UTC -  "; cout.flush();
    wat::Time date("now"); cout << date.GetDateString()     << endl; 
    cout << "-----------------------------------------"     << endl << endl;

    cout << "PE_GW_NAME          " << gOPT.gw_name          << endl;
    cout << "PE_DATA_LABEL       " << gOPT.data_label       << endl;
    cout << "PE_WAVE_DIR         " << gOPT.wave_dir         << endl;
    cout << "PE_NIFO             " << gOPT.nifo             << endl;
    cout << "PE_NEVENTS          " << gOPT.nevents          << endl;
    cout << "PE_NCYCLES          " << gOPT.ncycles          << endl;
    cout << "PE_IFAR_THR         " << gOPT.ifar_thr         << endl;
    cout << "PE_GPS_EVENT        " << std::setprecision(14) << gOPT.gps_event << std::setprecision(6) << endl;
    cout << "PE_TCOA_SHIFT       " << gOPT.tcoa_shift       << endl;
    cout << "PE_GPS_ERROR        " << gOPT.gps_error        << endl;
    cout << "PE_OFF_EVENT        " << gOPT.off_event        << endl;
    cout << "PE_INJ_TIME_STEP    " << gOPT.inj_time_step    << endl;
    cout << "PE_INJ_ONSRC        " << gOPT.inj_onsrc        << endl;
    cout << "PE_CCUT             " << gOPT.ccut             << endl;
    cout << "PE_RCUT             " << gOPT.rcut             << endl;
    cout << "PE_LOGL             " << gOPT.logl             << endl;
    cout << "PE_RSTAT            " << gOPT.rstat            << endl;
    cout << "PE_DUMP_CR          " << gOPT.dump_cr          << endl;
    cout << "PE_USE_ONSRC        " << gOPT.use_onsrc        << endl;
    cout << "PE_SYNC_PHASE       " << sync_phase            << endl;
    cout << "PE_SYNC_TIME        " << sync_time             << endl;
    cout << "PE_SENSITIVITY      " << sensitivity           << endl;
    cout << "PE_RESIDUALS        " << residuals             << endl;
    cout << "PE_STATISTICS       " << statistics            << endl;
    cout << "PE_EFRACTION        " << gOPT.efraction        << endl;
    cout << "PE_SIDE             " << gOPT.side             << endl;
    cout << "PE_COVERAGE         " << gOPT.coverage         << endl;
    cout << "PE_PLOT_FNAME       " << gOPT.plot_fname       << endl;
    cout << "PE_PLOT_ODIR        " << gOPT.plot_odir        << endl;
    cout << "PE_PLOT_TYPE        " << gOPT.plot_type        << endl;
    cout << "PE_PLOT_EFRACTION   " << gOPT.plot_efraction   << endl;
    cout << "PE_PLOT_ZOOM        " << gOPT.plot_zoom        << endl;
    cout << "PE_PLOT_CR          " << plot_cr               << endl;
    cout << "PE_PLOT_99PERC      " << plot_99perc           << endl;
    cout << "PE_MAKE_FF          " << make_ff               << endl;
    cout << "PE_MAKE_OF          " << make_of               << endl;
    cout << "PE_MAKE_RE          " << make_re               << endl;
    cout << "PE_MAKE_NE          " << make_ne               << endl;
    cout << "PE_MAKE_ONSRC       " << make_onsrc            << endl;

    cout << endl;

    cout << "CCUT_WDM_FRES       " << gOPT.ccut_wdm_fres   << endl;   	
    cout << "CCUT_BCHIRP         " << gOPT.ccut_bchirp     << endl;	
    cout << "CCUT_UCHIRP         " << gOPT.ccut_uchirp     << endl;
    cout << "CCUT_LTIME          " << gOPT.ccut_ltime      << endl;
    cout << "CCUT_RTIME          " << gOPT.ccut_rtime      << endl;
    cout << "CCUT_FMAX           " << gOPT.ccut_fmax       << endl;

    cout << endl;

    cout << "RCUT_WDM_FRES       " << gOPT.rcut_wdm_fres   << endl;   	
    cout << "RCUT_THR            " << gOPT.rcut_thr        << endl;	

    cout << endl;

    cout << "LOGL_ENABLED        " << logl_enabled         << endl;
    cout << "LOGL_FLOW           " << gOPT.logl_flow       << endl;   	
    cout << "LOGL_FHIGH          " << gOPT.logl_fhigh      << endl;   	
    cout << "LOGL_ARTHR          " << gOPT.logl_arthr      << endl;   	
    cout << "LOGL_IFO_MASK       " << "0x" << std::hex << gOPT.logl_ifo_mask << std::dec << endl;   	
    cout << "LOGL_RESAMPLE       " << gOPT.logl_resample   << endl;   	

    cout << endl;

    cout << "RSTAT_ENABLED       " << rstat_enabled        << endl;
    cout << "RSTAT_TYPE          " << gOPT.rstat_type      << endl;
    cout << "RSTAT_RTRIALS       " << gOPT.rstat_rtrials   << endl;
    cout << "RSTAT_JTRIALS       " << gOPT.rstat_jtrials   << endl;

    cout << endl;

    if(ofname!="") {
      // restore cout 
      std::cout.rdbuf(old);
      // dump config to ofname file
      ofstream out;
      out.open(ofname.Data(),ios::out);
      out.precision(14);
      out << buffer.str() << endl;
      out.close();
    }
}

// Dumps reconstructed waveform: time rec posterior med and CR in ASCII & ROOT format
void DumpWaveformCR(int nIFO) {

  // save rec/med/cr to output ascii file
  char ofname[1024];
  for(int n=0; n<nIFO; n++) {

    sprintf(ofname,"%s_%s_cr.txt",gIFO[n].Data(),gOPT.plot_fname.Data());

    ofstream out;
    out.open(ofname,ios::out);
    if (!out.good()) {cout << "Error Opening Output File : " << ofname << endl;exit(1);}
    cout << "Create Output File : " << ofname << endl;
    out.precision(19);

    // write header
    if(gOPT.plot_type=="time" || gOPT.plot_type=="envelope") { 
      out << "#whitened data : time, amp_cwb_rec, " 
          << "amp_post_median, amp_post_lower_50_cr, amp_post_lower_90_cr, amp_post_upper_50_cr, amp_post_upper_90_cr" << endl;
    } else if(gOPT.plot_type=="spectrum") {
      out << "#whitened data : frequency, amp_cwb_rec, " 
          << "amp_post_median, amp_post_lower_50_cr, amp_post_lower_90_cr, amp_post_upper_50_cr, amp_post_upper_90_cr" << endl;
    } else if(gOPT.plot_type=="phase") {
      out << "#whitened data : cycles, amp_cwb_rec, " 
          << "amp_post_median, amp_post_lower_50_cr, amp_post_lower_90_cr, amp_post_upper_50_cr, amp_post_upper_90_cr" << endl;
    }

    // write data
    int size = vREC[n].size();
    double dx=1./vREC[n].rate();
    for (int i=0; i<size; i++) {
      double x = i*dx;
      if(gOPT.plot_type=="time" || gOPT.plot_type=="envelope") x+=vTREC[n];
      double vl50 = vMED[n][i]-fabs(vL50[n][i]);
      double vu50 = vMED[n][i]+fabs(vU50[n][i]);
      double vl90 = vMED[n][i]-fabs(vL90[n][i]);
      double vu90 = vMED[n][i]+fabs(vU90[n][i]);

      out << x
          << " " << vREC[n][i] << " " << vMED[n][i] 
          << " " << vl50 << " " << vl90 << " " << vu50 << " " << vu90
          << endl;
    }

    out.close();
  }

  // open output root file to store rec,med,cr
  sprintf(ofname,"%s_cr.root",gOPT.plot_fname.Data());
  TFile *froot = new TFile(ofname, "RECREATE");
  if(froot==NULL) {
    cout << "Failed to create file !!! " <<  ofname << endl;
    gSystem->Exit(1);
  }
  for(int n=0; n<nIFO; n++) {
    // save start times
    double tREC = vREC[n].start();
    double tMED = vMED[n].start();
    double tL50 = vL50[n].start();
    double tL90 = vL90[n].start();
    double tL99 = vL99[n].start();
    double tU50 = vU50[n].start();
    double tU90 = vU90[n].start();
    double tU99 = vU99[n].start();

    if(gOPT.plot_type=="time" || gOPT.plot_type=="envelope") { 
      // set the start GPS time
      vREC[n].start(vTREC[n]);
      vMED[n].start(vTREC[n]);
      vL50[n].start(vTREC[n]);
      vL90[n].start(vTREC[n]);
      vL99[n].start(vTREC[n]);
      vU50[n].start(vTREC[n]);
      vU90[n].start(vTREC[n]);
      vU99[n].start(vTREC[n]);
    }

    // write data to output root file
    vREC[n].Write(TString::Format("REC_%s",gIFO[n].Data()));
    vMED[n].Write(TString::Format("MED_%s",gIFO[n].Data()));
    vL50[n].Write(TString::Format("L50_%s",gIFO[n].Data()));
    vL90[n].Write(TString::Format("L90_%s",gIFO[n].Data()));
    vL99[n].Write(TString::Format("L99_%s",gIFO[n].Data()));
    vU50[n].Write(TString::Format("U50_%s",gIFO[n].Data()));
    vU90[n].Write(TString::Format("U90_%s",gIFO[n].Data()));
    vU99[n].Write(TString::Format("U99_%s",gIFO[n].Data()));

    // restore start times 
    vREC[n].start(tREC);
    vMED[n].start(tMED);
    vL50[n].start(tL50);
    vL90[n].start(tL90);
    vL99[n].start(tL99);
    vU50[n].start(tU50);
    vU90[n].start(tU90);
    vU99[n].start(tU99);
  }
  froot->Close();
}

void MakeDistributions(TString type) {

  if(type!="FF" && type!="OF" && type!="RE" && type!="NE") {
    cout << "MakeDistributions - Error : plot type not available, must be FF/OF/RE/NE" << endl;
    gSystem->Exit(1);
  }

  TCanvas* cmf = new TCanvas("MatchingFactor", "MatchingFactor", 200, 20, 1000, 600);
  cmf->SetLogx(false);
  cmf->SetLogy(true);
  cmf->SetGridx(true);
  cmf->SetGridy(true);

  gStyle->SetOptStat("emrou");
  gStyle->SetStatFont(63);
  gStyle->SetStatFontSize(12);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.12);

  double xmin=1e20;
  double xmax=-1e-20;
  if(type=="FF") {
    for(int i=0;i<vFF.size();i++) {
      if(vFF[i]<xmin) xmin=vFF[i];
      if(vFF[i]>xmax) xmax=vFF[i];
    }
    //xmin=-1.0;
    //xmin=0.0;
    xmax=1.0;
  }
  if(type=="OF") {
    for(int i=0;i<vOF.size();i++) {
      if(vOF[i]<xmin) xmin=vOF[i];
      if(vOF[i]>xmax) xmax=vOF[i];
    }
  }
  if(type=="RE") {
    double xmean=0.;
    for(int i=0;i<vRE.size();i++) {
      if(vRE[i]<xmin) xmin=vRE[i];
      if(vRE[i]>xmax) xmax=vRE[i]+0.1*xmax;
      xmean+=vRE[i];
    }
    // limit the xmax when big outliers are present
    xmean/=vRE.size();
    if(xmax>6*xmean) xmax=6*xmean;
    if(xmax<vRE[0])  xmax=2*vRE[0];
  }
  if(type=="NE") {
    for(int i=0;i<vNE.size();i++) {
      if(vNE[i]<xmin) xmin=vNE[i];
      if(vNE[i]>xmax) xmax=vNE[i]+0.1*xmax;
    }
  }

  // histogram of matching factor onsource reconstructed vs offsource injected
  TH1F* xhist = new TH1F("xhist","xhist",100,xmin,xmax);
  //xhist->SetStats(kFALSE);
  xhist->SetLineWidth(2);
  xhist->SetLineColor(kGreen+2);
  for(int i=1;i<xFF.size();i++) {
    if(type=="FF") xhist->Fill(xFF[i]);
    if(type=="OF") xhist->Fill(xOF[i]);
    if(type=="RE") xhist->Fill(xRE[i]);
    if(type=="NE") xhist->Fill(xNE[i]);
  }

  int nFF=0;
  int nOF=0;
  int nRE=0;
  int nNE=0;
  // histogram of matching factor onsource/offsource reconstructed vs onsource/offsource injected 
  TH1F* hist = new TH1F("hist","hist",100,xmin,xmax);
  //hist->SetStats(kFALSE);
  if(type=="FF") for(int i=1;i<vFF.size();i++) hist->Fill(vFF[i]);
  if(type=="OF") for(int i=1;i<vOF.size();i++) hist->Fill(vOF[i]);
  if(type=="RE") for(int i=1;i<vRE.size();i++) hist->Fill(vRE[i]);
  if(type=="NE") for(int i=1;i<vNE.size();i++) hist->Fill(vNE[i]);

  int nTOT=0;
  for(int i=1;i<vFF.size();i++) {
    if(type=="FF") {nTOT++;if(vFF[i]<vFF[0]) nFF++;}
    if(type=="OF") {nTOT++;if(vOF[i]<vOF[0]) nOF++;}
    if(type=="RE") {nTOT++;if(vRE[i]>vRE[0]) nRE++;}
    if(type=="NE") {nTOT++;if(vNE[i]>vNE[0]) nNE++;}
  }

  double pvalueFF = double(nFF)/double(nTOT);
  double pvalueOF = double(nOF)/double(nTOT);
  double pvalueRE = double(nRE)/double(nTOT);
  double pvalueNE = double(nNE)/double(nTOT);

  cout << endl;
  if(type=="FF") cout << "pvalue FF=" << pvalueFF << endl;
  if(type=="OF") cout << "pvalue OF=" << pvalueOF << endl;
  if(type=="RE") cout << "pvalue RE=" << pvalueRE << endl;
  if(type=="NE") cout << "pvalue NE=" << pvalueNE << endl;
  cout << endl;

  hist->SetLineWidth(2);

  hist->Draw("HIST");
  cmf->Modified(); cmf->Update();
  TPaveStats *sthist = (TPaveStats*)hist->FindObject("stats");
  sthist->SetTextColor(kBlack);
  sthist->SetLineColor(kBlack);
  if(gOPT.make_onsrc) {
    sthist->SetTextColor(kBlue);
    xhist->Draw("SAMES");
    cmf->Modified(); cmf->Update();
    TPaveStats *stxhist = (TPaveStats*)xhist->FindObject("stats");
    stxhist->SetTextColor(kGreen+2);
    stxhist->SetLineColor(kBlack);
    stxhist->SetY2NDC(.68);
    stxhist->SetOptStat(1110);
  }
  float yhmax = hist->GetMaximum()>xhist->GetMaximum() ? hist->GetMaximum() : xhist->GetMaximum();
  hist->GetYaxis()->SetRangeUser(0.9,1.1*yhmax);
  hist->GetYaxis()->SetNdivisions(3);

  double mean = hist->GetMean();
  double rms  = hist->GetRMS();

  // compute offsource median. 50%, 90%, 99%
  int nentries=vFF.size()-1;	// index=0 contains onsource
  int *index = new int[nentries];
  float *value = new float[nentries];
  for(int i=0;i<nentries;i++) {
    if(type=="FF") value[i]=vFF[i+1];
    if(type=="OF") value[i]=vOF[i+1];
    if(type=="RE") value[i]=vRE[i+1];
    if(type=="NE") value[i]=vNE[i+1];
  }
  TMath::Sort(nentries,value,index,false);

  int imed = (nentries*50.)/100.;  if(imed>=nentries) imed=nentries-1;
  int il50 = (nentries*25.)/100.;  if(il50>=nentries) il50=nentries-1;
  int iu50 = (nentries*75.)/100.;  if(iu50>=nentries) iu50=nentries-1;
  int il90 = (nentries*5.)/100.;   if(il90>=nentries) il90=nentries-1;
  int iu90 = (nentries*95.)/100.;  if(iu90>=nentries) iu90=nentries-1;
  int il99 = (nentries*0.5)/100.;  if(il99>=nentries) il99=nentries-1;
  int iu99 = (nentries*99.5)/100.; if(iu99>=nentries) iu99=nentries-1;
  int il100 = 0;
  int iu100 = nentries-1;

  double med = value[index[imed]];
  double l50 = value[index[il50]];
  double u50 = value[index[iu50]];
  double l90 = value[index[il90]];
  double u90 = value[index[iu90]];
  double l99 = value[index[il99]];
  double u99 = value[index[iu99]];
  double l100 = value[index[il100]];
  double u100 = value[index[iu100]];

  double match=-1;
  if(type=="FF") match=vFF[0];
  if(type=="OF") match=vOF[0];
  if(type=="RE") match=vRE[0];
  if(type=="NE") match=vNE[0];

  double pvalue=0;
  if(type=="FF") pvalue=pvalueFF;
  if(type=="OF") pvalue=pvalueOF;
  if(type=="RE") pvalue=pvalueRE;
  if(type=="NE") pvalue=pvalueNE;

  TString xtitle="";
  if(type=="FF") xtitle="Fitting Factor";
  if(type=="OF") xtitle="Overlap Factor";
  if(type=="RE") xtitle="Residual Energy";
  if(type=="NE") xtitle="Null Energy";

  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetLabelOffset(0.02);
  hist->GetXaxis()->SetNoExponent();

  hist->GetYaxis()->SetTitle("counts");
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->GetYaxis()->SetLabelOffset(0.02);
  hist->GetYaxis()->SetNoExponent();
  if(gOPT.rstat_enabled) hist->GetXaxis()->SetTitle("R statistic");
  char htitle[1024];
  if(gOPT.rstat_enabled) {
    sprintf(htitle,"%s %s (%s) - on-source (red - R-stat=%0.4f, pvalue=%0.4f) vs off-source (blue)",
            gOPT.gw_name.Data(),xtitle.Data(),gOPT.side.Data(),match,pvalue);
  } else {
    sprintf(htitle,"%s %s (%s) - cWB (red - %s=%0.4f, pvalue=%0.4f) vs off-source (blue)",
            gOPT.gw_name.Data(),xtitle.Data(),gOPT.side.Data(),type.Data(),match,pvalue);
  }
  hist->SetTitle(htitle);
  hist->SetName("Statistic");

  // draw vertical line at the cWB point estimate FF/OF/RE/NE value
  float ymax = hist->GetMaximum();
  float xval = 0; 
  if(type=="FF") xval=vFF[0];
  if(type=="OF") xval=vOF[0];
  if(type=="RE") xval=vRE[0];
  if(type=="NE") xval=vNE[0];
  TLine *line = new TLine(xval,0,xval,ymax);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->Draw();

  // set stats box line color black
//  cmf->Update();
//  TPaveStats *ptstats = (TPaveStats*)cmf->GetPrimitive("stats");
//  ptstats->SetLineColor(kBlack);
//  cmf->Modified();

  TString pdir=gOPT.plot_odir;
  xtitle.ReplaceAll(" ","");
  char ofname[256];
  sprintf(ofname,"%s_%s.png",gOPT.plot_fname.Data(),xtitle.Data());
  TString fName = pdir!="" ? TString(pdir)+TString("/")+ofname : ofname;
  cmf->Print(fName);
  cout << "write : " << fName << endl;
  fName.ReplaceAll(".png",".root");
  cmf->Print(fName);
  cout << "write : " << fName << endl;

  //if(type=="RE") {xtitle="ResidualSNR"; match=sqrt(match);}
  // dump summary statistic
  fName.ReplaceAll(".root",".txt");
  TString params = TString::Format("%s : %s %0.4f pvalue %0.4f median %0.4f l99 %0.4f l90 %0.4f l50 %0.4f u50 %0.4f u90 %0.4f u99 %0.4f",
                                   gOPT.gw_name.Data(),xtitle.Data(),match,pvalue,med,l99,l90,l50,u50,u90,u99);
  DumpStatistics(fName, params, false);

  if(gOPT.side=="ccut") { 	// dump ccut parameters & match & pvalue
    params = TString::Format("%s : %s %0.4f pvalue %0.4f ccut_wdm_fres %d ccut_bchirp %0.4f ccut_uchirp %0.4f ccut_ltime %0.4f ccut_rtime %0.4f",
                             gOPT.gw_name.Data(),xtitle.Data(),match,pvalue,gOPT.ccut_wdm_fres,gOPT.ccut_bchirp,gOPT.ccut_uchirp,gOPT.ccut_ltime,gOPT.ccut_rtime);
    DumpStatistics("pestat_ccut.txt", params, false);
  }

  if(gOPT.side=="ccut") { 	// dump ccut parameters & match & pvalue
    DumpStatistics("pestat_ccut.lst", "id\tpvalue\tresidual", false);
    for(int n=0;n<vFF.size();n++) {
      int nRE=0; for(int i=1;i<vFF.size();i++) if(vRE[i]>vRE[n]) nRE++;
      double pvalueRE = nRE/(double)(vRE.size()-1);
      TString parm = TString::Format("%d\t%.6f\t%.6f",n,pvalueRE,vRE[n]); 
      DumpStatistics("pestat_ccut.lst", parm, true);
    }
  }

  TString sfName = fName;	// save output file name, used as template for the following auxiliary output files

  // dump sorted statistic
  fName=sfName;fName.ReplaceAll(".txt","_statistic.txt");
  for(int i=0;i<nentries;i++) {
    TString parm = TString::Format("%f",value[index[i]]); 
    if(i==0) DumpStatistics(fName, parm, false);
    else     DumpStatistics(fName, parm, true);
  }

  // dump factors, used with ONSPE multi injections
  if(type=="NE") {
    int fmed  = wFACTOR[index[imed]];
    int fl50  = wFACTOR[index[il50]];
    int fu50  = wFACTOR[index[iu50]];
    int fl90  = wFACTOR[index[il90]];
    int fu90  = wFACTOR[index[iu90]];
    int fl99  = wFACTOR[index[il99]];
    int fu99  = wFACTOR[index[iu99]];
    int fl100 = wFACTOR[index[il100]];
    int fu100 = wFACTOR[index[iu100]];
    fName=sfName;fName.ReplaceAll(".txt","_factor.txt");
    TString params = TString::Format("%s : %s %0.4f pvalue %0.4f median %d l100 %d l99 %d l90 %d l50 %d u50 %d u90 %d u99 %d u100 %d",
                                   gOPT.gw_name.Data(),xtitle.Data(),match,pvalue,fmed,fl100,fl99,fl90,fl50,fu50,fu90,fu99,fu100);
    DumpStatistics(fName, params, false);
  }

  if(gOPT.sensitivity) {

    double dphase50,dtime50,damp50;
    double dphase90,dtime90,damp90;
    if(type=="FF") {
      double dff_phase,dff_time,dff_amp;
      ComputePhaseSensitivity("ff",l90,u90,dff_phase,dphase90);
      ComputePhaseSensitivity("ff",l50,u50,dff_phase,dphase50);
      printf("\n\n------> FF Phase Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dff_phase);
      ComputeTimeSensitivity("ff",l90,u90,dff_time,dtime90);
      ComputeTimeSensitivity("ff",l50,u50,dff_time,dtime50);
      printf("\n\n------> FF Time  Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dff_time);
      ComputeAmpSensitivity("ff",l90,u90,dff_amp,damp90);
      ComputeAmpSensitivity("ff",l50,u50,dff_amp,damp50);
      printf("\n\n------> FF Amp   Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dff_amp);
      params = TString::Format("%s : %s\tdff_phase\t%0.4f\tdff_time\t%0.4f\tdff_amp\t%0.4f\tdphase50\t%0.4f\tdtime50\t%0.5f\tdamp50\t%0.4f\tdphase90\t%0.4f\tdtime90\t%0.5f\tdamp90\t%0.4f",
                               gOPT.gw_name.Data(),xtitle.Data(),dff_phase,dff_time,dff_amp,dphase50,dtime50,damp50,dphase90,dtime90,damp90);
    }
    if(type=="OF") {
      double dof_phase,dof_time,dof_amp;
      ComputePhaseSensitivity("of",l90,u90,dof_phase,dphase90);
      ComputePhaseSensitivity("of",l50,u50,dof_phase,dphase50);
      printf("\n\n------> OF Phase Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dof_phase);
      ComputeTimeSensitivity("of",l90,u90,dof_time,dtime90);
      ComputeTimeSensitivity("of",l50,u50,dof_time,dtime50);
      printf("\n\n------> OF Time  Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dof_time);
      ComputeAmpSensitivity("of",l90,u90,dof_amp,damp90);
      ComputeAmpSensitivity("of",l50,u50,dof_amp,damp50);
      printf("\n\n------> OF Amp   Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dof_amp);
      params = TString::Format("%s : %s\tdof_phase\t%0.4f\tdof_time\t%0.4f\tdof_amp\t%0.4f\tdphase50\t%0.4f\tdtime50\t%0.5f\tdamp50\t%0.4f\tdphase90\t%0.4f\tdtime90\t%0.5f\tdamp90\t%0.4f",
                               gOPT.gw_name.Data(),xtitle.Data(),dof_phase,dof_time,dof_amp,dphase50,dtime50,damp50,dphase90,dtime90,damp90);
    }
    if(type=="RE") {
      double dre_phase,dre_time,dre_amp;
      ComputePhaseSensitivity("re",l90,u90,dre_phase,dphase90);
      ComputePhaseSensitivity("re",l50,u50,dre_phase,dphase50);
      printf("\n\n------> RE Phase Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dre_phase);
      ComputeTimeSensitivity("re",l90,u90,dre_time,dtime90);
      ComputeTimeSensitivity("re",l50,u50,dre_time,dtime50);
      printf("\n\n------> RE Time  Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dre_time);
      ComputeAmpSensitivity("re",l90,u90,dre_amp,damp90);
      ComputeAmpSensitivity("re",l50,u50,dre_amp,damp50);
      printf("\n\n------> RE Amp   Sensitivity (l50=%.2f - u50=%.2f) = %.4f\n\n",l50,u50,dre_amp);
      params = TString::Format("%s : %s\tdre_phase\t%0.4f\tdre_time\t%0.4f\tdre_amp\t%0.4f\tdphase50\t%0.4f\tdtime50\t%0.5f\tdamp50\t%0.4f\tdphase90\t%0.4f\tdtime90\t%0.5f\tdamp90\t%0.4f",
                               gOPT.gw_name.Data(),xtitle.Data(),dre_phase,dre_time,dre_amp,dphase50,dtime50,damp50,dphase90,dtime90,damp90);
    }
  
    if(type!="NE") {
      fName=sfName;fName.ReplaceAll(".txt","_sensitivity.txt");
      DumpStatistics(fName, params, false);
    }
  }

  // compute 90% percentile bounds of the on-source distribution
  if(gOPT.make_onsrc) {
    int xnentries=xFF.size();
    int* xindex = new int[xnentries];
    float* xvalue = new float[xnentries];
    for(int i=0;i<xnentries;i++) {
      if(type=="FF") xvalue[i]=xFF[i];
      if(type=="OF") xvalue[i]=xOF[i];
      if(type=="RE") xvalue[i]=xRE[i];
      if(type=="NE") xvalue[i]=xNE[i];
    }
    TMath::Sort(xnentries,xvalue,xindex,false);
    int ixl50 = (xnentries*25.)/100.;  if(ixl50>=xnentries) ixl50=xnentries-1;
    int ixu50 = (xnentries*75.)/100.;  if(ixu50>=xnentries) ixu50=xnentries-1;
    int ixl90 = (xnentries*5.)/100.;   if(ixl90>=xnentries) ixl90=xnentries-1;
    int ixu90 = (xnentries*95.)/100.;  if(ixu90>=xnentries) ixu90=xnentries-1;
    int ixl99 = (xnentries*0.5)/100.;  if(ixl99>=xnentries) ixl99=xnentries-1;
    int ixu99 = (xnentries*99.5)/100.; if(ixu99>=xnentries) ixu99=xnentries-1;
    double xl50 = xvalue[xindex[ixl50]];
    double xu50 = xvalue[xindex[ixu50]];
    double xl90 = xvalue[xindex[ixl90]];
    double xu90 = xvalue[xindex[ixu90]];
    double xl99 = xvalue[xindex[il99]];
    double xu99 = xvalue[xindex[ixu99]];
    // compute on-source distribution probability @ 50%,90%,99%
    int bxl50 = hist->FindBin(xl50);				// lower 50% xhist bound bin
    int bxu50 = hist->FindBin(xu50);				// upper 50% xhist bound bin
    int bxl90 = hist->FindBin(xl90);				// lower 90% xhist bound bin
    int bxu90 = hist->FindBin(xu90);				// upper 90% xhist bound bin
    int bxl99 = hist->FindBin(xl99);				// lower 99% xhist bound bin
    int bxu99 = hist->FindBin(xu99);				// upper 99% xhist bound bin
    hist->Scale(1./hist->Integral());                             // distribution is normalized to 1
    double prob50=0;
    for(int i=bxl50;i<=bxu50;i++) prob50+=hist->GetBinContent(i); // get probability of xhist distr @ 50%
    double prob90=0;
    for(int i=bxl90;i<=bxu90;i++) prob90+=hist->GetBinContent(i); // get probability of xhist distr @ 90%
    double prob99=0;
    for(int i=bxl99;i<=bxu99;i++) prob99+=hist->GetBinContent(i); // get probability of xhist distr @ 99%
    // dump to file on-source distribution probability @ 50%,90%,99%
    fName=sfName;fName.ReplaceAll(".txt","_onsource_prob.txt");
    params = TString::Format("%s : prob50 %0.4f prob90 %0.4f prob99 %0.4f",
                             gOPT.gw_name.Data(),prob50,prob90,prob99);
    DumpStatistics(fName, params, false);
    cout << endl << xtitle << " on-source probability -> " << params << endl << endl;

    // dump sorted statistic
    fName=sfName;fName.ReplaceAll(".txt","_onsource_statistic.txt");
    for(int i=0;i<xnentries;i++) {
      TString parm = TString::Format("%f",xvalue[xindex[i]]); 
      if(i==0) DumpStatistics(fName, parm, false);
      else     DumpStatistics(fName, parm, true);
    }

    delete [] xindex;
    delete [] xvalue;
/*
    // dump on-source marginal likelihood used for bayes factor
    double prob=0;
    for(int i=0;i<=hist->GetNbinsX();i++) {
      prob += hist->GetBinContent(i)*xhist->GetBinContent(i);
    }
    if(hist->Integral()*xhist->Integral()>0) prob /= hist->Integral()*xhist->Integral();
    else prob=0.;
    TString parm = TString::Format("on-source marginal likelihood = %g", prob);
    fName=sfName;fName.ReplaceAll(".txt","_mlike.txt");
    DumpStatistics(fName, parm, false);
    cout << endl << sfName.ReplaceAll(".txt"," -> ") << parm << endl << endl;
*/
  }

  delete hist;
  delete xhist;
  delete cmf;
  delete [] index;
  delete [] value;
}

double
GetInjTcoa(double geocentric_tcoa, network* NET, TString ifo, double theta, double phi) {
// compute coalescence time

  // Add Time Delay respect to geocenter
  CWB::mdc MDC(NET);
  double tdelay = MDC.GetDelay(ifo,"",phi,theta);
  double ifo_tcoa = geocentric_tcoa+tdelay;

  return ifo_tcoa;
}

double
GetDelay(network* NET, TString ifo, double theta, double phi) {

  // Add Time Delay respect to geocenter
  CWB::mdc MDC(NET);
  double tdelay = MDC.GetDelay(ifo,"",phi,theta);

  return tdelay;
}

void
DumpStatistics(TString ofname, TString params, bool app) {

   ofstream out;
   if(app) out.open(ofname.Data(),ios::app);
   else    out.open(ofname.Data(),ios::out);
   out.precision(14);
   out << params.Data() << endl;
   out.close();
}

void 
ComputePhaseSensitivity(TString match, double lower, double upper, double& dmatch, double& dphase) {

  if((match!="ff")&&(match!="of")&&(match!="re")) {dmatch=dphase=0.;return;}

  vector<wavearray<double> > xINJ = vINJ;		 

  bool   lfound,ufound;
  lfound=ufound=false;
  double lphase,uphase;
  lphase=uphase=0.;
  cout.precision(4); 
  for(int i=0;i<180;i++) {
    if(lfound&&ufound) continue;
    for(int n=0;n<xINJ.size();n++) CWB::mdc::PhaseShift(xINJ[n],1.);
    double value = CWB::mdc::GetMatchFactor(match,xINJ,vINJ);
    if(match=="ff" || match=="of") {
      if((!ufound)&&(value>upper)) uphase=i+1; else ufound=true;
      if((!lfound)&&(value>lower)) lphase=i+1; else lfound=true;
    } else {
      if((!ufound)&&(value<upper)) uphase=i+1; else ufound=true;
      if((!lfound)&&(value<lower)) lphase=i+1; else lfound=true;
    }
    //cout << match.ToUpper() << "\tComputePhaseSensitivity " << i+1 << "\t" << value << "\t" << lower << "\t" << upper << endl;
  }
  dmatch = (uphase!=lphase) ? (upper-lower)/(uphase-lphase) : 0;
  dmatch/=(upper+lower)/2.;
  dmatch=fabs(dmatch);
  dphase=fabs(uphase-lphase);
}

void 
ComputeTimeSensitivity(TString match, double lower, double upper, double& dmatch, double& dtime) {

  if((match!="ff")&&(match!="of")&&(match!="re")) {dmatch=dtime=0.;return;}

  vector<wavearray<double> > xINJ = vINJ;		 

  double dt=0.00001;

  bool   lfound,ufound;
  lfound=ufound=false;
  double ltime,utime;
  ltime=utime=0.;
  cout.precision(4); 
  for(int i=0;i<10000;i++) {
    if(lfound&&ufound) continue;
    double t=(i+1)*dt;
    for(int n=0;n<xINJ.size();n++) CWB::mdc::TimeShift(xINJ[n],dt);
    double value = CWB::mdc::GetMatchFactor(match,xINJ,vINJ);
    if(match=="ff" || match=="of") {
      if((!ufound)&&(value>upper)) utime=t; else ufound=true;
      if((!lfound)&&(value>lower)) ltime=t; else lfound=true;
    } else {
      if((!ufound)&&(value<upper)) utime=t; else ufound=true;
      if((!lfound)&&(value<lower)) ltime=t; else lfound=true;
    }
    //cout << match.ToUpper() << "\tComputeTimeSensitivity " << t << "\t" << value << "\t" << lower << "\t" << upper << endl;
  }
  dmatch = (utime!=ltime) ? (upper-lower)/(utime-ltime) : 0;
  dmatch/=(upper+lower)/2.;
  dmatch=fabs(dmatch);
  dtime=fabs(utime-ltime);
}

void 
ComputeAmpSensitivity(TString match, double lower, double upper, double& dmatch, double& damp) {

  if((match!="of")&&(match!="re")) {dmatch=damp=0.;return;}

  double da=0.99;

  bool lfound,ufound;
  double lamp,uamp;
  cout.precision(4); 
  if(match=="of") {
    lamp=uamp=0.;
    lfound=ufound=false;
    vector<wavearray<double> > xINJ = vINJ;		 
    for(int n=0;n<xINJ.size();n++) xINJ[n]*=pow(da,-200);
    for(int i=0;i<401;i++) {
      if(lfound&&ufound) continue;
      double a = pow(da,(i-200));
      double value = CWB::mdc::GetMatchFactor(match,xINJ,vINJ);
      if((!ufound)&&(value>upper)) uamp=a; else ufound=true;
      if((!lfound)&&(value>lower)) lamp=a; else lfound=true;
      //cout << i-200 << " " << "\tComputeAmpSensitivity " << a << "\t" << value << "\t" << lower << "\t" << upper << endl;
      for(int n=0;n<xINJ.size();n++) xINJ[n]*=da;
    }
    dmatch = (uamp!=lamp) ? (upper-lower)/(uamp-lamp) : 0;
    dmatch/=(upper+lower)/2.;
    dmatch=fabs(dmatch);
    damp=fabs(uamp-lamp);
    return;
  }
  if(match=="re") {
    // amplitude factor <=1 
    lamp=uamp=0.;
    lfound=ufound=false;
    vector<wavearray<double> > xINJ = vINJ;		 
    for(int i=0;i>-200;i--) {
      if(lfound&&ufound) continue;
      double a = pow(da,i);
      double value = CWB::mdc::GetMatchFactor(match,xINJ,vINJ);
      if((!ufound)&&(value<upper)) uamp=a; else ufound=true;
      if((!lfound)&&(value<lower)) lamp=a; else lfound=true;
      //cout << i << " " << "\tComputeAmpSensitivity " << a << "\t" << value << "\t" << lower << "\t" << upper << endl;
      for(int n=0;n<xINJ.size();n++) xINJ[n]*=1./da;
    }
    double damp1=fabs(uamp-lamp);
    double dmatch1 = (uamp!=lamp) ? (upper-lower)/(uamp-lamp) : 0;
    dmatch1/=(upper+lower)/2.;

    // amplitude factor >=1 
    lamp=uamp=0.;
    lfound=ufound=false;
    xINJ = vINJ;		 
    for(int i=0;i<200;i++) {
      if(lfound&&ufound) continue;
      double a = pow(da,i);
      double value = CWB::mdc::GetMatchFactor(match,xINJ,vINJ);
      if((!ufound)&&(value<upper)) uamp=a; else ufound=true;
      if((!lfound)&&(value<lower)) lamp=a; else lfound=true;
      //cout << i << " " << "\tComputeAmpSensitivity " << a << "\t" << value << "\t" << lower << "\t" << upper << endl;
      for(int n=0;n<xINJ.size();n++) xINJ[n]*=da;
    }
    double damp2=fabs(uamp-lamp);
    double dmatch2 = (uamp!=lamp) ? (upper-lower)/(uamp-lamp) : 0;
    dmatch2/=(upper+lower)/2.;

    dmatch=(fabs(dmatch1)+fabs(dmatch2))/2.;
    damp=(damp1+damp2)/2.;
    return;
  }
}

double 
GetSyncTcoa(wavearray<double> wf, double sync_time, double sync_phase, double tcoa) {

// this function return the time shift of the tcoa after the time/phase sync
//
// wf is the wINJ after time/phase sync
// sync_time, sync_phase are the time/phase sync
// tcoa is the coalescence time of wINJ

  CWB::mdc::TimePhaseSync(wf, -sync_time, -sync_phase);   // return to the original wf before time/phase sync

  wavearray<double> fwf = CWB::Toolbox::getHilbertIFrequency(wf);       // get frequency vs time

  double dt = 1./wf.rate();
  int n = int((tcoa-wf.start())/dt);
  float inj_fcoa = fwf[n];                      // get frequency at the coalescence time

  double tshift = inj_fcoa>0 ? sync_time+(1./inj_fcoa)*(sync_phase/360.) : sync_time; // get time shift after time/phase sync

  return tcoa+tshift;
}

float 
TestStatistic(int nIFO, int id, float ref) {
// Used for tests 

  if(wREC[id].size()==0) return -1.;

  int nTRY = 100;

  gRandom->SetSeed(150914+id);

  // select nTRY random onsource samples
  vector<int> irnd(nTRY);
  for(int i=0;i<nTRY;i++) irnd[i] = int(gRandom->Uniform(0,gEVENTS));

  gnetwork* NET = new gnetwork(nIFO,gIFO);

  int N=0;
  int M=0;
  float min=1e20;
  std::vector<wavearray<double> > wrec(nIFO);
  for(int j=0;j<nTRY;j++) {
    int i=irnd[j];
    if(wINJ[i].size()==0) continue;
    // phase time,sync the wREC sample waveforms respect to wINJ
    for(int n=0;n<nIFO;n++) {
      // align wINJ wrt wREC
      wrec[n] = CWB::mdc::GetAligned(id<0 ? &vREC[n] : &wREC[id][n], &wINJ[i][n]);
      char sout[1024];
      double sync_phase, sync_time, sync_xcor;
      if(gOPT.sync_time && !gOPT.sync_phase) {
        sync_xcor = CWB::mdc::TimeSync(wrec[n], wINJ[i][n], sync_time);
        sprintf(sout,"wINJ Time Sync -> ifo = %*s - evtId = %*d / %*d - sync_time = %*.4f sec - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 5, i, 5 , nTRY, 7, sync_time, 6, sync_xcor);
        //cout << sout << endl;
      }
      if(gOPT.sync_phase && !gOPT.sync_time) {
        sync_xcor = CWB::mdc::PhaseSync(wrec[n], wINJ[i][n], sync_phase);
        sprintf(sout,"wINJ Phase Sync -> ifo = %*s - evtId = %*d / %*d - sync_phase = %*.1f deg - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 5, i, 5 , nTRY, 6, sync_phase, 6, sync_xcor);
        //cout << sout << endl;
      }
      if(gOPT.sync_phase && gOPT.sync_time) {
        sync_xcor = CWB::mdc::TimePhaseSync(wrec[n], wINJ[i][n], sync_time, sync_phase);
        sprintf(sout,"wINJ TimePhase Sync -> ifo = %*s - evtId = %*d / %*d - sync_time = %*.4f sec - sync_phase = %*.1f deg - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 5, i, 5 , nTRY, 7, sync_time, 6, sync_phase, 6, sync_xcor);
        //cout << sout << endl;
      }
    }

    vector<double> tstart;
    vector<double> tstop(nIFO);
    // select the left side of the waveform
    for(int n=0;n<nIFO;n++) {
      double inj_tcoa = gOPT.inj_time_step+fmod(wTCOA[i],1);
      tstop[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], wTHETA[i], wPHI[i]);
    }

    //float ET=0; for(int n=0;n<nIFO;n++) for(int j=0;j<wINJ[i][n].size();j++) ET+=pow(wINJ[i][n][j],2);
    double re = CWB::mdc::GetMatchFactor("re",wrec, wINJ[i],tstart,tstop);
    if(re>ref) M++;
    if(re<min) min=re;
    N++;
  }
  delete NET;

  float pvalue=(float)M/(float)N;
  return pvalue;
}

int
GetOnSourceStatistic(int nIFO, double &mFF, double &mOF, double &mRE) {

  // compute statistics onsource/offsource reconstructed vs onsource/offsource injected 
  // results are stored in xFF,xOF,xRE,xNE
  // WARNING!!! the offsource whitened waveforms are aligned in time to the onsource reconstructed waveform
  //            -> from this point the offsource whitened waveforms can not be used for the standard comparison

  cout << endl << "Compute onsource distribution : be patient, it takes a while ..." << endl << endl;

  gnetwork* NET = new gnetwork(nIFO,gIFO);

  // align the injected sample waveforms respect to the point estimated 
  for(int n=0;n<nIFO;n++) {
    for(int i=0;i<gEVENTS;i++) {
      if(wREC[i].size()) {      // event detected
        wINJ[i][n] = CWB::mdc::GetAligned(&wINJ[i][n], &vREC[n]);
      } else {                  // event rejected (set to 0)
        wREC[i][n] = wINJ[i][n]; wREC[i][n]=0;
      }
    }
  }

  // phase time,sync the injected waveforms respect to the point estimated waveforms
  // NOTE: by default the time sync is applied, if requested also the phase sync is applied
  int pc = 0;
  int ipc = double(gEVENTS)/10.; if(ipc==0) ipc=1;
  for(int i=0;i<gEVENTS;i++) {
    if(i%ipc==0) {if(gEVENTS>100) {cout << pc<<"%";if (pc<100) cout << " - ";pc+=10;cout.flush();}}
    for(int n=0;n<nIFO;n++) {
      char sout[1024];
      double sync_phase, sync_time, sync_xcor;
      if(!gOPT.sync_phase) {
        sync_xcor = CWB::mdc::TimeSync(wINJ[i][n], vREC[n], sync_time);
        sprintf(sout,"Time Sync -> ifo = %*s - evtId = %*d / %*d - sync_time = %*.4f sec - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 5, i+1, 5 , gEVENTS, 7, sync_time, 6, sync_xcor);
        //cout << sout << endl;
      }
      if(gOPT.sync_phase && gOPT.sync_time) {
        sync_xcor = CWB::mdc::TimePhaseSync(wINJ[i][n], vREC[n], sync_time, sync_phase);
        sprintf(sout,"TimePhase Sync -> ifo = %*s - evtId = %*d / %*d - sync_time = %*.4f sec - sync_phase = %*.1f deg - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 5, i+1, 5 , gEVENTS, 7, sync_time, 6, sync_phase, 6, sync_xcor);
        //cout << sout << endl;
      }
    }
  }
  if(pc==100) cout << pc<<"%";;
  cout << endl << endl;

  // cWB point estimate Matching Factors respect to Posterior Samples Median
  std::vector<double> zFF,zRE;
  for(int i=0;i<gEVENTS;i++) {

    if(wREC[i].size()==0) continue;

    vector<double> tstart;
    if(gOPT.side=="right") {
      tstart.resize(nIFO);
      for(int n=0;n<nIFO;n++) {
        double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
        tstart[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
      }
    }

    vector<double> tstop;
    if(gOPT.side=="left") {
      tstop.resize(nIFO);
      for(int n=0;n<nIFO;n++) {
        double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
        tstop[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
      }
    }

    double FF = CWB::mdc::GetMatchFactor("ff",vREC,wINJ[i],tstart,tstop);
    double OF = CWB::mdc::GetMatchFactor("of",vREC,wINJ[i],tstart,tstop);
    double RE = CWB::mdc::GetMatchFactor("re",vREC,wINJ[i],tstart,tstop);

    xFF.push_back(FF);
    xOF.push_back(OF);
    xRE.push_back(RE);

    if(gOPT.residuals) {	// null energy
      vector<wavearray<double> > vDAT = vNUL; for(int n=0;n<nIFO;n++) vDAT[n]+=vREC[n];
      double NE = CWB::mdc::GetMatchFactor("re",vDAT,wINJ[i],tstart,tstop);
      xNE.push_back(NE);
    }

    // compute FF/RE full/left side, used to get median sample for onsource comparision
    if(gOPT.use_onsrc=="mfull" || gOPT.use_onsrc=="mleft") {
      tstart.resize(0);
      if(gOPT.use_onsrc=="mfull") {
        tstop.resize(0);
      } else {
        tstop.resize(nIFO);
        for(int n=0;n<nIFO;n++) {
          double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
          tstop[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
        }
      }
      FF = CWB::mdc::GetMatchFactor("ff",vREC,wINJ[i],tstart,tstop);
      zFF.push_back(FF);
    } else {	// efull/eleft   
      tstart.resize(0);
      if(gOPT.use_onsrc=="efull") {
        tstop.resize(0);
      } else {
        tstop.resize(nIFO);
        for(int n=0;n<nIFO;n++) {
          double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
          tstop[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
        }
      }
      RE = CWB::mdc::GetMatchFactor("re",vREC,wINJ[i],tstart,tstop);
      zRE.push_back(RE);
    }
  }

  int imed=-1;
  if(gOPT.use_onsrc!="") {
    // compute onsource median
    int nentries = (gOPT.use_onsrc=="mfull" || gOPT.use_onsrc=="mleft") ? zFF.size() : zRE.size();
    int *index = new int[nentries];
    float *value = new float[nentries];

    // compute median FF/RE full/left side onsource distribution
    for(int i=0;i<nentries;i++) value[i] = (gOPT.use_onsrc=="mfull" || gOPT.use_onsrc=="mleft") ? zFF[i] : zRE[i];
    TMath::Sort(nentries,value,index,false);

    imed = int((nentries*50.)/100.);

    // get medians 
    mFF = xFF[index[imed]];
    mOF = xOF[index[imed]];
    mRE = xRE[index[imed]];

    delete [] index;
    delete [] value;
  }

  delete NET;

  return imed;
}

void
GetOffSourceStatistic(int nIFO, int recId, vector<int> injId,
                      vector<double>& oFF, vector<double>& oOF, vector<double>& oRE) {

  // compute statistics offsource reconstructed vs offsource injected
  // results are stored in oFF,oOF,oRE
  // WARNING!!! the offsource whitened waveforms are aligned in time to the onsource reconstructed waveform
  //            -> from this point the offsource whitened waveforms can not be used for the standard comparison

  if(recId>=gEVENTS) return;

  gnetwork* NET = new gnetwork(nIFO,gIFO);

  int kEVENTS = gOPT.rstat_jtrials;

  // set reconstructed off-source event
  std::vector<wavearray<double> > oREC = recId>=0 ? wREC[recId] : vREC;

  std::vector<wavearray<double> > oINJ[PE_MAX_EVENTS]; 	// OffSource iwhitened injected posteriors

  // align the injected sample waveforms respect to the point estimated
  for(int k=0;k<kEVENTS;k++) {
    int i=injId[k];
    for(int n=0;n<nIFO;n++) oINJ[i].push_back(CWB::mdc::GetAligned(&wINJ[i][n], &oREC[n]));
  }

  // phase time,sync the injected waveforms respect to the point estimated waveforms
  // NOTE: by default the time sync is applied, if requested also the phase sync is applied
  for(int k=0;k<kEVENTS;k++) {
    int i=injId[k];

    // compute optimal phase/time shifts for each detector
    char sout[1024];
    double sync_phase[NIFO_MAX], sync_time[NIFO_MAX], sync_xcor[NIFO_MAX];
    for(int n=0;n<nIFO;n++) {
      if(gOPT.sync_phase) {	// apply time/phase sync
        sync_xcor[n] = CWB::mdc::TimePhaseSync(oINJ[i][n], oREC[n], sync_time[n], sync_phase[n]);
        sprintf(sout,"TimePhase Sync -> ifo = %*s - evtId = %*d / %*d - sync_time = %*.4f sec - sync_phase = %*.1f deg - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 5, i+1, 5 , gEVENTS, 7, sync_time[n], 6, sync_phase[n], 6, sync_xcor[n]);
      } else {			// apply time sync
        sync_xcor[n] = CWB::mdc::TimeSync(oINJ[i][n], oREC[n], sync_time[n]);
        sprintf(sout,"Time Sync -> ifo = %*s - evtId = %*d / %*d - sync_time = %*.4f sec - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 5, i+1, 5 , gEVENTS, 7, sync_time[n], 6, sync_xcor[n]);
      }
      //cout << sout << endl;
    }
/*
    // computed weighted average sync phase/time
    double wsnr2=0.;
    double wsync_phase=0.;
    double wsync_time=0.;
    for(int n=0;n<nIFO;n++) {
        wsync_phase+=pow(wSNR[k][n],2)*sync_phase[n];
        wsync_time+=pow(wSNR[k][n],2)*sync_time[n];
        wsnr2+=pow(wSNR[k][n],2);
    }
    wsync_phase/=wsnr2;
    wsync_time /=wsnr2;

    // restore original time/phase and apply global weighted average sync phese/time
    for(int n=0;n<nIFO;n++) {
      if(!gOPT.sync_phase) {
        CWB::mdc::TimeSync(oINJ[i][n], wsync_time-sync_time[n]);
        if(n==0) sprintf(sout,"Time Sync -> - evtId = %*d / %*d - sync_time = %*.4f sec",
                         5, k+1, 5 , gEVENTS, 7, wsync_time);
      }
      if(gOPT.sync_phase && gOPT.sync_time) {
        CWB::mdc::TimePhaseSync(oINJ[i][n], wsync_time-sync_time[n], wsync_phase-sync_phase[n]);
        if(n==0) sprintf(sout,"TimePhase Sync -> - evtId = %*d / %*d - sync_time = %*.4f sec - sync_phase = %*.1f deg",
                         5, k+1, 5 , gEVENTS, 7, wsync_time, 6, wsync_phase);
      }
    }
*/
    //cout << sout << endl;
  }

  // cWB point estimate Matching Factors respect to Posterior Samples Median
  for(int k=0;k<kEVENTS;k++) {
    int i=injId[k];

    if(wREC[i].size()==0) continue;

    vector<double> tstart;
    if(gOPT.side=="right") {
      tstart.resize(nIFO);
      for(int n=0;n<nIFO;n++) {
        double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
        tstart[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
      }
    }

    vector<double> tstop;
    if(gOPT.side=="left") {
      tstop.resize(nIFO);
      for(int n=0;n<nIFO;n++) {
        double inj_tcoa = gOPT.inj_time_step+fmod(vTCOA,1);
        tstop[n] = GetInjTcoa(inj_tcoa, NET, gIFO[n], vTHETA, vPHI);
      }
    }

    if(oINJ[i].size()!=oREC.size()) {
      cout << "GetOffSourceStatistic - Error: the elements of input injId vector are not unique !!!" << endl; 
      gSystem->Exit(1);
    }

    double FF = CWB::mdc::GetMatchFactor("ff",oREC,oINJ[i],tstart,tstop);
    double OF = CWB::mdc::GetMatchFactor("of",oREC,oINJ[i],tstart,tstop);
    double RE = CWB::mdc::GetMatchFactor("re",oREC,oINJ[i],tstart,tstop);

    oFF.push_back(FF);
    oOF.push_back(OF);
    oRE.push_back(RE);
  }

  delete NET;

  return;
}

void MakePvalueDistribution(TString type) {

  if(type!="FF" && type!="OF" && type!="RE") {
    cout << "MakePvalueDistribution - Error : plot type not available, must be FF/OF/REE" << endl;
    gSystem->Exit(1);
  }

  TCanvas* cpd = new TCanvas("PvalueDistribution", "PvalueDistribution", 200, 20, 1000, 600);
  cpd->SetLogx(false);
  cpd->SetLogy(true);
  cpd->SetGridx(true);
  cpd->SetGridy(true);

  gStyle->SetOptStat("emr");

  // histogram of matching factor onsource reconstructed vs offsource injected
  TH1F* hist = new TH1F("hpvalue","hpvalue",100,0,1);
  //hist->SetStats(kFALSE);
  hist->SetLineWidth(2);
  hist->SetLineColor(kGreen+2);

  // sort null hyphotesis values
  int nentries=vFF.size()-1;	// index=0 contains onsource
  int *index = new int[nentries];
  float *value = new float[nentries];
  for(int i=0;i<nentries;i++) {
    if(type=="FF") value[i]=vFF[i+1];
    if(type=="OF") value[i]=vOF[i+1];
    if(type=="RE") value[i]=vRE[i+1];
  }
  TMath::Sort(nentries,value,index,false);

  int onsize=xFF.size();
  float *pvalue = new float[onsize];
  for(int i=0;i<onsize;i++) {
    double x; 
    if(type=="FF") x=xFF[i];
    if(type=="OF") x=xOF[i];
    if(type=="RE") x=xRE[i];
    bool loop=true;
    // find pvalue of x using dichotomic search
    int imed=0;
    int imin=0;
    int imax=nentries-1;
    while(loop) {
      imed = imin+int((imax-imin)*50./100.);  
      if(imed>=imax) imed=imax;
      double xmed=value[index[imed]];
      imin = (x<xmed) ? imin : imed;
      imax = (x<xmed) ? imed : imax;
      if(imax==imin+1) loop=false;
    }
    if(type=="FF") pvalue[i]=double(imed)/double(nentries);
    if(type=="OF") pvalue[i]=double(imed)/double(nentries);
    if(type=="RE") pvalue[i]=double(nentries-imed)/double(nentries);
    hist->Fill(pvalue[i]);
  }
  delete [] index;
  delete [] value;

  // compute pvalue median and 90% lower/upper limits
  int *pindex = new int[onsize];
  TMath::Sort(onsize,pvalue,pindex,false);
  float pvalue_med = pvalue[pindex[int(onsize*50./100.)]];
  float pvalue_l90 = pvalue[pindex[int(onsize* 5./100.)]];
  float pvalue_u90 = pvalue[pindex[int(onsize*95./100.)]];

  hist->GetXaxis()->SetTitle("pvalue");
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetLabelOffset(0.02);
  hist->GetXaxis()->SetNoExponent();

  hist->GetYaxis()->SetTitle("counts");
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->GetYaxis()->SetLabelOffset(0.02);
  hist->GetYaxis()->SetNoExponent();
  hist->GetYaxis()->SetNdivisions(3);

  TString xtitle="";
  if(type=="FF") xtitle="Fitting Factor";
  if(type=="OF") xtitle="Overlap Factor";
  if(type=="RE") xtitle="Residual Energy";

  char htitle[256];
  sprintf(htitle,"%s %s (%s) - %s on-source p-value distribution (l90=%f, med=%f, u90=%f)",
          gOPT.gw_name.Data(),xtitle.Data(),gOPT.side.Data(),type.Data(),pvalue_l90,pvalue_med,pvalue_u90);
  hist->SetTitle(htitle);
  hist->SetName("Statistic");
  hist->Draw("HIST");

  // draw vertical lines at the pvalue med/l90/u90
  float ymax = hist->GetMaximum();
  TLine *line_med = new TLine(pvalue_med,0,pvalue_med,ymax);
  line_med->SetLineWidth(2);
  line_med->SetLineColor(kRed);
  line_med->Draw();
  TLine *line_l90 = new TLine(pvalue_l90,0,pvalue_l90,ymax);
  line_l90->SetLineWidth(2);
  line_l90->SetLineColor(kBlue);
  line_l90->Draw();
  TLine *line_u90 = new TLine(pvalue_u90,0,pvalue_u90,ymax);
  line_u90->SetLineWidth(2);
  line_u90->SetLineColor(kBlue);
  line_u90->Draw();

  // set stats box line color black
  cpd->Update();
  TPaveStats *ptstats = (TPaveStats*)cpd->GetPrimitive("stats");
  ptstats->SetLineColor(kBlack);
  cpd->Modified();

  TString pdir=gOPT.plot_odir;
  xtitle.ReplaceAll(" ","");
  char ofname[256];
  sprintf(ofname,"%s_%s_%s.png",gOPT.plot_fname.Data(),"onsource_pvalue_distribution",xtitle.Data());
  TString fName = pdir!="" ? TString(pdir)+TString("/")+ofname : ofname;
  cpd->Print(fName);
  cout << "write : " << fName << endl;
  fName.ReplaceAll(".png",".root");
  cpd->Print(fName);
  cout << "write : " << fName << endl;

  // dump sorted pvalue onsource statistic
  sprintf(ofname,"%s_%s_%s.txt",gOPT.plot_fname.Data(),xtitle.Data(),"onsource_pvalue_statistic");
  fName = pdir!="" ? TString(pdir)+TString("/")+ofname : ofname;
  for(int i=0;i<onsize;i++) {
    TString parm = TString::Format("%f",pvalue[pindex[i]]);
    if(i==0) DumpStatistics(fName, parm, false);
    else     DumpStatistics(fName, parm, true);
  }

  delete [] pindex;
  delete [] pvalue;

  delete line_med;
  delete line_l90;
  delete line_u90;
  delete hist;
  delete cpd;
}

int SyncWaveforms(int nIFO, TString type) {

// for type=statistics the reference is INJ, this choice permits to fix the tcoa 
// for type=waveform the reference is REC, this choice is necessary to fix the original onsource waveform

  if(type!="statistics" && type!="waveform") return 0;

  int selected=0;
  if(type=="statistics" || gOPT.plot_type=="phase") {
    // align onsource/offsource sample waveforms respect to the reconstructed 
    for(int n=0;n<nIFO;n++) {                      
      selected=0;
      CWB::mdc::Align(vINJ[n], vREC[n]); 							// onsource waveform
      if(gOPT.residuals) vNUL[n]=CWB::mdc::GetAligned(&vNUL[n], &vREC[n]);			// onsource null
      for(int i=0;i<gEVENTS;i++) {            
        if(wREC[i].size()) {									// event detected
          CWB::mdc::Align(wINJ[i][n], wREC[i][n]); 						// offsource waveform
          if(gOPT.residuals) wNUL[i][n]=CWB::mdc::GetAligned(&wNUL[i][n], &wREC[i][n]);		// offsource null
          if(gOPT.side=="rcut") wAUX[i][n]=CWB::mdc::GetAligned(&wAUX[i][n], &wREC[i][n]);	// offsource auxiliary
          selected++;
        }
      } 
      if(selected==0) return 0;	// return if no events are detects
    }                                                                    
  } else {
    // align the reconstructed waveforms respect to the point estimated 
    for(int n=0;n<nIFO;n++) {                      
      selected=0;
      for(int i=0;i<gEVENTS;i++) {            
        if(wREC[i].size()) {	// event detected
          wREC[i][n] = CWB::mdc::GetAligned(&wREC[i][n], &vREC[n]); 
          if(gOPT.residuals) wNUL[i][n] = CWB::mdc::GetAligned(&wNUL[i][n], &vNUL[n]); 
          if(gOPT.side=="rcut") wAUX[i][n] = CWB::mdc::GetAligned(&wAUX[i][n], &vAUX[n]); 
          selected++;
        } else {			// event rejected (set to 0)
          wREC[i][n] = vREC[n]; wREC[i][n]=0;
          if(gOPT.residuals) wNUL[i][n] = vNUL[n]; wNUL[i][n]=0;
          if(gOPT.side=="rcut") wAUX[i][n] = vAUX[n]; wAUX[i][n]=0;
        }
      } 
      if(selected==0) return 0;	// return if no events are detects
    }                                                                    
  }

  if(gOPT.sync_phase || gOPT.sync_time) {

    // OnSource - phase time,sync the reconstructed waveforms respect to the injected waveforms
    for(int n=0;n<nIFO;n++) {

      wavearray<double>* wref = type=="statistics" || gOPT.plot_type=="phase" ? &vINJ[n] : &vREC[n]; 
      wavearray<double>* wsht = type=="statistics" || gOPT.plot_type=="phase" ? &vREC[n] : &vREC[n]; 

      char sout[1024];
      double sync_phase, sync_time, sync_xcor;
      if(gOPT.sync_time && !gOPT.sync_phase) {
        sync_xcor = CWB::mdc::TimeSync(*wsht, *wref, sync_time);
        vTSYNC[n]+=sync_time;
        sprintf(sout,"Time Sync -> ifo = %*s - sync_time = %*.4f sec - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 7, sync_time, 6, sync_xcor);
        cout << sout << endl;
        if(gOPT.side=="rcut") CWB::mdc::TimeSync(vAUX[n], sync_time);
      }
      if(gOPT.sync_phase && !gOPT.sync_time) {
        sync_xcor = CWB::mdc::PhaseSync(*wsht, *wref, sync_phase);
        vPSYNC[n]+=sync_phase;
        sprintf(sout,"Phase Sync -> ifo = %*s - sync_phase = %*.1f deg - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 6, sync_phase, 6, sync_xcor);
        cout << sout << endl;
        if(gOPT.side=="rcut") CWB::mdc::PhaseSync(vAUX[n], sync_phase);
      }
      if(gOPT.sync_phase && gOPT.sync_time) {
        sync_xcor = CWB::mdc::TimePhaseSync(*wsht, *wref, sync_time, sync_phase);
        vTSYNC[n]+=sync_time;
        vPSYNC[n]+=sync_phase;
        sprintf(sout,"TimePhase Sync -> ifo = %*s - sync_time = %*.4f sec - sync_phase = %*.1f deg - sync_xcor = %*.3f",
                     2, gIFO[n].Data(), 7, sync_time, 6, sync_phase, 6, sync_xcor);
        cout << sout << endl;
        if(gOPT.side=="rcut") CWB::mdc::TimePhaseSync(vAUX[n], sync_time, sync_phase);
      }
    }

    // OffSource - phase time,sync the reconstructed waveforms respect to the injected waveforms
    for(int i=0;i<gEVENTS;i++) {
      for(int n=0;n<nIFO;n++) {
  
        wavearray<double>* wref = type=="statistics" || gOPT.plot_type=="phase" ? &wINJ[i][n] : &vREC[n]; 
        wavearray<double>* wsht = type=="statistics" || gOPT.plot_type=="phase" ? &wREC[i][n] : &wREC[i][n]; 

        char sout[1024];
        double sync_phase, sync_time, sync_xcor;
        if(gOPT.sync_time && !gOPT.sync_phase) {
          sync_xcor = CWB::mdc::TimeSync(*wsht, *wref, sync_time);
          sprintf(sout,"Time Sync -> ifo = %*s - evtId = %*d / %*d - sync_time = %*.4f sec - sync_xcor = %*.3f",
                       2, gIFO[n].Data(), 5, i+1, 5 , gEVENTS, 7, sync_time, 6, sync_xcor);
          cout << sout << endl;
          if(gOPT.side=="rcut") CWB::mdc::TimeSync(wAUX[i][n], sync_time);
        }
        if(gOPT.sync_phase && !gOPT.sync_time) {
          sync_xcor = CWB::mdc::PhaseSync(*wsht, *wref, sync_phase);
          sprintf(sout,"Phase Sync -> ifo = %*s - evtId = %*d / %*d - sync_phase = %*.1f deg - sync_xcor = %*.3f",
                       2, gIFO[n].Data(), 5, i+1, 5 , gEVENTS, 6, sync_phase, 6, sync_xcor);
          cout << sout << endl;
          if(gOPT.side=="rcut") CWB::mdc::PhaseSync(wAUX[i][n], sync_phase);
        }
        if(gOPT.sync_phase && gOPT.sync_time) {
          sync_xcor = CWB::mdc::TimePhaseSync(*wsht, *wref, sync_time, sync_phase);
          sprintf(sout,"TimePhase Sync -> ifo = %*s - evtId = %*d / %*d - sync_time = %*.4f sec - sync_phase = %*.1f deg - sync_xcor = %*.3f",
                       2, gIFO[n].Data(), 5, i+1, 5 , gEVENTS, 7, sync_time, 6, sync_phase, 6, sync_xcor);
          cout << sout << endl;
          if(gOPT.side=="rcut") CWB::mdc::TimePhaseSync(wAUX[i][n], sync_time, sync_phase);
        }
      }
    }
  } 

  if(gOPT.side=="rcut") {
    cout << endl << "Apply RCut : be patient, it takes a while ..." << endl << endl;
    int pc = 0;
    int ipc = double(gEVENTS)/10.; if(ipc==0) ipc=1;
    for(int n=0;n<nIFO;n++) {                      
      for(int i=0;i<gEVENTS;i++) {            
        if(n==nIFO-1 && i%ipc==0) {if(gEVENTS>100) {cout << pc<<"%";if (pc<100) cout << " - ";pc+=10;cout.flush();}}
        wINJ[i][n] = GetRCut(&wINJ[i][n],&wAUX[i][n]);
        wREC[i][n] = GetRCut(&wREC[i][n],&wAUX[i][n]);
      } 
      vINJ[n] = GetRCut(&vINJ[n],&vAUX[n]);
      vREC[n] = GetRCut(&vREC[n],&vAUX[n]);
    }                                                                    
    if(pc==100) cout << pc<<"%";;
    cout << endl << endl;
  } 

  // convert time wavearray into envelope/frequency/spectrum 
  if(type=="waveform") {
    for(int n=0;n<nIFO;n++) {                      
      for(int i=0;i<gEVENTS;i++) {            
        // convert wREC[i][n] time into envelope
        if(gOPT.plot_type=="envelope")  wREC[i][n] = CWB::Toolbox::getHilbertEnvelope(wREC[i][n]);
        // convert wREC[i][n] time into frequency
        if(gOPT.plot_type=="frequency") wREC[i][n] = CWB::Toolbox::getHilbertIFrequency(wREC[i][n]);
        // convert wREC[i][n] time into spectrum oneside
        if(gOPT.plot_type=="spectrum")  wREC[i][n] = CWB::mdc::GetSpectrum(&wREC[i][n],true);
      } 
      // convert vREC[n] time into envelope
      if(gOPT.plot_type=="envelope")  vREC[n] = CWB::Toolbox::getHilbertEnvelope(vREC[n]);
      // convert vREC[n] time into frequency
      if(gOPT.plot_type=="frequency") vREC[n] = CWB::Toolbox::getHilbertIFrequency(vREC[n]);
      // convert vREC[n] time into spectrum oneside
      if(gOPT.plot_type=="spectrum")  vREC[n] = CWB::mdc::GetSpectrum(&vREC[n],true);
    }                                                                    

    // convert vREC[n] time into phase
    wavearray<double> fi,fr,s,t;
    if(gOPT.plot_type=="phase") {  
      for(int n=0;n<nIFO;n++) {                      
        for(int i=0;i<gEVENTS;i++) {            
          wREC[i][n] = CWB::Toolbox::GetPhase(wINJ[i][n], wREC[i][n], fi, fr, s, t);
          wREC[i][n].resize(gOPT.ncycles);
          wREC[i][n]*=180.;
        }
      }
      for(int n=0;n<nIFO;n++) {
        pSNR.push_back(s);	// cumulative SNR of vINJ[n]+vREC[n]
        pFRQ.push_back(fi); 	// vINJ[n] frequency
        pTIM.push_back(t); 	// vINJ[n] time
        vREC[n] = CWB::Toolbox::GetPhase(vINJ[n], vREC[n], pFRQ[n], fr, pSNR[n], pTIM[n]);
        vREC[n].resize(gOPT.ncycles);
        vREC[n]*=180.;
        // set max pFRQ and max pSNR to 180
        pFRQ[n]*=180./pFRQ[n].max();
        pSNR[n]*=180./pSNR[n].max();
      }
    }
  }

  return selected;
}

TString GetParameterFromInjLog(TString log, TString param) {
// get params from log string                                                     

  TObjArray* token = log.Tokenize(TString(" "));
  for(int i=0;i<token->GetEntries();i++) {
    TObjString* otoken = (TObjString*)token->At(i);
    TString stoken = otoken->GetString();
    // exctract value
    if(stoken==param) {
      if(i<token->GetEntries()-1) {
        if(stoken=="spin1" || stoken=="spin2") {
          otoken = (TObjString*)token->At(i+3);		// get spinz
          return otoken->GetString();
        } else {
          otoken = (TObjString*)token->At(i+1);
          return otoken->GetString();
        }
      }
    }
  }
  if(token) delete token;

  return "";
}

double GetCCut(wavearray<double>* ts, double tcoa, double m1, double m2, double s1z, double s2z) {
// apply chirp TF cuts and compute total energy inside the chirp cut region

  // NOTE: the current inmplementation works only with SEOBNRv4ROM model

  const auto& fchirp = static_cast<double(*)(double,double,double,double,double)>(CWB::mdc::SimIMRSEOBNRv4ROMFrequencyOfTime);
  const auto& tchirp = static_cast<double(*)(double,double,double,double,double)>(CWB::mdc::SimIMRSEOBNRv4ROMTimeOfFrequency);

  // produce the TF map:
  WSeries<double> W;                            // TF map container
  W.Forward(*ts, *gWDM);                        // apply the WDM to the time series 

  double fmin=16;
  // set the initial value of fmax using the global variable gOPT.ccut_fmax
  // the LAL functions XLALSimIMRSEOBNRv4ROMFrequencyOfTime/XLALSimIMRSEOBNRv4ROMTimeOfFrequency use spline and works only
  // for a given frequency interval which depends on masses and spins of the binary system
  double fmax= ts->rate()/2.>gOPT.ccut_fmax ? gOPT.ccut_fmax : ts->rate()/2.;

  int nPx=30;					// number of points (nPx) used to compute the up/bottom chirp lines
  double T[nPx],Fb[nPx],Fu[nPx];
  double dF = (fmax-fmin)/nPx;
  for(int i=0;i<nPx;i++) {
    double F = fmin+i*dF;
    Fb[i] = gOPT.ccut_bchirp*F;			// alpha bottom
    Fu[i] = gOPT.ccut_uchirp*F;			// alpha up
    double t=tchirp(F,m1,m2,s1z,s2z);
    if(t<0) {					// the LAL function return error because frequency is not in the available frequency range
      cout << "cwb_pereport:GetCCut - update gOPT.ccut_fmax frequency to " << F
           << " Hz in order to avoid errors in LAL functions XLALSimIMRSEOBNRv4ROMFrequencyOfTime/XLALSimIMRSEOBNRv4ROMTimeOfFrequency" << endl;
      gOPT.ccut_fmax=fmin+(i-1)*dF;		// update gOPT.ccut_fmax frequency
      nPx=i;					// update the number of points (nPx) used to compute the up/bottom chirp lines
    } else {
      T[i] = tcoa-t;
    }
  }
  TSpline3 tbchirp("tbchirp",Fb,T,nPx);		// bottom chirp line 	(time vs freq)
  TSpline3 tuchirp("tuchirp",Fu,T,nPx);		// up chirp line	(time vs freq)
  TSpline3 fbchirp("fbchirp",T,Fb,nPx);		// bottom chirp line	(freq vs time)
  TSpline3 fuchirp("fuchirp",T,Fu,nPx);		// up chirp line	(freq vs time)

  double df = W.resolution();                   // frequency pixel resolution (hz)
  double dt = 1./(2*df);                        // time pixel resolution (sec)

  int layers = W.maxLayer()+1;                  // numbers of frequency bins (first & last bins have df/2)
  int slices = W.sizeZero();                    // number of time bins

  double tl = tcoa-gOPT.ccut_ltime;		// left time cut
  double tr = tcoa-gOPT.ccut_rtime;		// right time cut

  // rescale the phase 00/90 pixel amplitudes

  for(int i=1;i<slices;i++) {                   // loop over time bins, skip first slice

    double time  = i*dt;			// pixel central time

    double tm = time-dt/2;			// pixel minimum time 
    double tM = time+dt/2;			// pixel maximum time

    if(tm<tl-dt || tM>tr+dt) {                  // remove pixels outside the time interval [tl-dt:tr+dt]
      for(int j=1;j<layers;j++) {               // loop over frequency bins
        W.putSample(0,i,j+0.01);    		// zero a00 pixel amplitude
        W.putSample(0,i,-(j+0.01)); 		// zero a90 pixel amplitude
      }
      continue;
    }

    double fb = fbchirp.Eval(time);		// fbchirp frequency at pixel central time
    double fu = fuchirp.Eval(time);		// fuchirp frequency at pixel central time

    if(TMath::IsNaN(fb)) fb = ts->rate()/2.;	// check if fb is NaN
    if(TMath::IsNaN(fu)) fu = ts->rate()/2.;	// check if fu is NaN

    double f1 = fbchirp.Eval(tm);		// fbchirp frequency at pixel minimum time
    double f2 = fbchirp.Eval(tM);		// fbchirp frequency at pixel maximum time
    double f3 = fuchirp.Eval(tm);		// fuchirp frequency at pixel minimum time
    double f4 = fuchirp.Eval(tM);		// fuchirp frequency at pixel maximum time

    for(int j=1;j<layers;j++) {                 // loop over frequency bins, skip first layer

      double freq = j*df;			// pixel central frequency

      double fm = freq-df/2;			// pixel minimum frequency
      double fM = freq+df/2;			// pixel maximum frequency

      if(fm<fmin || fM>fmax) {
        W.putSample(0,i,j+0.01);    		// zero a00 pixel amplitude
        W.putSample(0,i,-(j+0.01)); 		// zero a90 pixel amplitude
        continue;
      }

      double t1 = tbchirp.Eval(fm);		// fbchirp time at pixel minimum frequency
      double t2 = tbchirp.Eval(fM);		// fbchirp time at pixel maximum frequency
      double t3 = tuchirp.Eval(fm);		// fuchirp time at pixel minimum frequency
      double t4 = tuchirp.Eval(fM);		// fuchirp time at pixel maximum frequency

      // check if pixel is inside the time cut range tl,tr
      bool ctime = (tM>tl && tm<tr) ? true : false; 

      // check if pixel is inside the frequency cut range fb,fu
      bool cfreq = (fM>fb && fm<fu) ? true : false; 

      double x1=t1,x2=t2,x3=t3,x4=t4,xm=tm,xM=tM;	// assign x* = t*
      double y1=f1,y2=f2,y3=f3,y4=f4,ym=fm,yM=fM;	// assign y* = f*

      // for pixels crossed by the cut times tl,tr assign new time pixel boundaries xm,xM and new f*
      if(tl>tm && tl<tM) {xm=tl; f1 = fbchirp.Eval(xm); f3 = fuchirp.Eval(xm);}
      if(tr>tm && tr<tM) {xM=tr; f2 = fbchirp.Eval(xM); f4 = fuchirp.Eval(xM);} 

      y1=f1,y2=f2,y3=f3,y4=f4;				// update y*

      // check if pixel is crossed by fbchirp,fuchirp lines
      // a pixel is crossed by the chirp lines if at least one of t*,f* points are inside the time/freq pixel ranges

      bool cline=false;							// false/true = not-crossed/crossed

      if(t1>xm && t1<xM) cline=true; 
      if(t2>xm && t2<xM) cline=true; 
      if(t3>xm && t3<xM) cline=true; 
      if(t4>xm && t4<xM) cline=true; 

      if(f1>ym && f1<yM) cline=true; 
      if(f2>ym && f2<yM) cline=true; 
      if(f3>ym && f3<yM) cline=true; 
      if(f4>ym && f4<yM) cline=true; 

      double Ap=(xM-xm)*(yM-ym);					// area pixel (is <=0.5 = total pixel area) 

      if(ctime&&cline) {	// pixels crossed by fbchirp,fuchirp lines and inside the time cut range tl,tr

        // for pixels with frequency outside the frequency pixel range assign ym,yM
        if(y1<ym) y1=ym; if(y1>yM) y1=yM;
        if(y2<ym) y2=ym; if(y2>yM) y2=yM;
        if(y3<ym) y3=ym; if(y3>yM) y3=yM;
        if(y4<ym) y4=ym; if(y4>yM) y4=yM;

	// compute area Ab for pixels crossed by fbchirp line 

        double Ab=0;							// area pixel above the fbchirp line

        // for pixels with times t1,t2 outside the time pixel range assign xm,xM and compute the corresponding new frequencies y1,y2
        if(t1<xm) {x1=xm; y1 = fbchirp.Eval(x1);}
        if(t1>xM) {x1=xM; y1 = fbchirp.Eval(x1);}
        if(t2<xm) {x2=xm; y2 = fbchirp.Eval(x2);}
        if(t2>xM) {x2=xM; y2 = fbchirp.Eval(x2);}

        // for pixels with frequency outside the frequency pixel range assign ym,yM
        if(y1<ym) y1=ym; if(y1>yM) y1=yM;
        if(y2<ym) y2=ym; if(y2>yM) y2=yM;

        // compute area pixel above the fbchirp line
        if(y1 >ym && y2==yM) Ab = (x2-xm)*(yM-y1)/2;
        if(y1 >ym && y2 <yM) Ab = (xM-xm)*(yM-y2)+(xM-xm)*(y2-y1)/2; 
        if(y1==ym && y2==yM) Ab = (x1-xm)*(yM-ym)+(x2-x1)*(yM-ym)/2;
        if(y1==ym && y2 <yM) Ab = Ap-(xM-x1)*(y2-ym)/2;

	// compute area Au for pixels crossed by fuchirp line 

        double Au=0;							// area pixel below the fuchirp line

        // for pixels with times t3,t4 outside the time pixel range assign xm,xM and compute the corresponding new frequencies y3,y4
        if(t3<xm) {x3=xm; y3 = fuchirp.Eval(x3);}
        if(t3>xM) {x3=xM; y3 = fuchirp.Eval(x3);}
        if(t4<xm) {x4=xm; y4 = fuchirp.Eval(x4);}
        if(t4>xM) {x4=xM; y4 = fuchirp.Eval(x4);}

        // for pixels with frequency outside the frequency pixel range assign ym,yM
        if(y3<ym) y3=ym; if(y3>yM) y3=yM;
        if(y4<ym) y4=ym; if(y4>yM) y4=yM;

        // compute area pixel below the fuchirp line
        if(y3 >ym && y4==yM) Au = Ap-(x4-xm)*(yM-y3)/2;
        if(y3 >ym && y4 <yM) Au = (xM-xm)*(y3-ym)+(xM-xm)*(y4-y3)/2; 
        if(y3==ym && y4==yM) Au = (xM-x4)*(yM-ym)+(x4-x3)*(yM-ym)/2;
        if(y3==ym && y4 <yM) Au = (xM-x3)*(y4-ym)/2;

        if(Ab<0) Ab=0; if(Ab>Ap) Ab=Ap;					// fix area inconsistency due to precision errors
        if(Au<0) Au=0; if(Au>Au) Au=Au;					// fix area inconsistency due to precision errors
        double A = Ab+Au; A-=Ap;					// combine Ab and Au when pixel is crossed by fbchirp fuchirp lines

        if(A<0)   {cout << "GetCCut Warning: pixel area < 0  " << endl; A=0;}
        if(A>Ap)  {cout << "GetCCut Warning: pixel area > Ap " << endl; A=Ap;}
        double R=sqrt(A/0.5);						// compute the rescaling amplitude factor (0.5 is the total pixel area)
/*
        TString tag = R==0 ? "*" : "";
        TString sR = TString::Format("GetCCut Rescaling Factor: Ap %0.3f \ttime %0.3f \tfreq %0.3f \tR %0.3f \t%s",Ap,time,freq,R,tag.Data());
        cout << sR << endl;
*/
        double a00 = W.getSample(i,j+0.01); W.putSample(a00*R,i,j+0.01);	// rescaling a00 pixel amplitude
        double a90 = W.getSample(i,-(j+0.01));W.putSample(a90*R,i,-(j+0.01));	// rescaling a90 pixel amplitude
      }

      if((ctime && cfreq) && !cline) { 	// pixels inside/outside the chirp cut area not belonging to the cline pixels 

        double R = (x4<=tl || x1>=tr) ? 0 : sqrt(Ap/0.5);       // compute the rescaling amplitude factor (0.5 is the total pixel area)
                                                                // remove pixels outside the time interval [tl,tr]
/*
        TString tag = "+";
        TString sR = TString::Format("GetCCut Rescaling Factor: Ap %0.3f \ttime %0.3f \tfreq %0.3f \tR %0.3f \t%s",Ap,time,freq,R,tag.Data());
        cout << sR << endl;
*/
        double a00 = W.getSample(i,j+0.01); W.putSample(a00*R,i,j+0.01);	// rescaling a00 pixel amplitude
        double a90 = W.getSample(i,-(j+0.01));W.putSample(a90*R,i,-(j+0.01));	// rescaling a90 pixel amplitude
      }

      //if(ctime && cfreq) { 
      //if(ctime && cline) { 
      //if((ctime && cfreq) && !cline) { 
      if((ctime && cline)||(ctime && cfreq)) { 
      } else {
        W.putSample(0,i,j+0.01);
        W.putSample(0,i,-(j+0.01));
      }
    }
  }

  double Et=0;                                  // Energy in time-frequency chirp cut region
  for(int i=1;i<slices;i++) {                   // loop over time bins
    for(int j=1;j<layers;j++) {                 // loop over frequency bins
      double a00 = W.getSample(i,j+0.01);       // read a00 pixel amplitude
      double a90 = W.getSample(i,-(j+0.01));    // read a90 pixel amplitude
      Et+=(a00*a00+a90*a90)/2;
    }
  }

  return Et;
}

void GetCCutParms(TString options) {

  // initialize with the default values
  gOPT.ccut_wdm_fres   = CCUT_WDM_FRES;
  gOPT.ccut_bchirp     = CCUT_BCHIRP;
  gOPT.ccut_uchirp     = CCUT_UCHIRP;
  gOPT.ccut_ltime      = CCUT_LTIME;
  gOPT.ccut_rtime      = CCUT_RTIME;
  gOPT.ccut_fmax       = CCUT_FMAX;

  // read the user parameter ccut
  TObjArray* token = TString(options).Tokenize(TString(':'));
  if(token->GetEntries()!=5) {
    cout << endl;
    cout << "cwb_pereport - Error : ccut format - must be: wdm_fres:bchirp:uchirp:ltime:rtime" <<  endl;
    cout << "               Default setup is: 16:1.35:1.65:0.5:0.0" << endl;
    cout << "               To setup a default value use *" << endl;
    cout << "               Example: *:1.25:1.75:*:*" << endl;
    cout << endl;
    exit(1);
  } 
  for(int j=0;j<token->GetEntries();j++) {
    TObjString* tok = (TObjString*)token->At(j);
    TString stok = tok->GetString();

    if(stok=="*") continue;

    if(j==0) {
      TString wdm_fres=stok;
      if(wdm_fres.IsDigit()) gOPT.ccut_wdm_fres=wdm_fres.Atoi();
      else {cout<<"cwb_pereport - Error : wdm_fres is not integer"<<endl;exit(1);}                   
    }
    if(j==1) {
      TString mc_lower_factor=stok;
      if(mc_lower_factor.IsFloat()) gOPT.ccut_bchirp=mc_lower_factor.Atof();
      else {cout<<"cwb_pereport - Error : mc_lower_factor is not float"<<endl;exit(1);}                   
    }
    if(j==2) {
      TString mc_upper_factor=stok;
      if(mc_upper_factor.IsFloat()) gOPT.ccut_uchirp=mc_upper_factor.Atof();
      else {cout<<"cwb_pereport - Error : mc_upper_factor is not float"<<endl;exit(1);}                   
    }
    if(j==3) {
      TString left_time=stok;
      if(left_time.IsFloat()) gOPT.ccut_ltime=left_time.Atof();
      else {cout<<"cwb_pereport - Error : left_time is not float"<<endl;exit(1);}                   
    }
    if(j==4) {
      TString right_time=stok;
      if(right_time.IsFloat()) gOPT.ccut_rtime=right_time.Atof();
      else {cout<<"cwb_pereport - Error : right_time is not float"<<endl;exit(1);}                   
    }
  }

  int x=gOPT.ccut_wdm_fres;        // must be power of 2
  if(!((x != 0) && ((x & (~x + 1)) == x)) || gOPT.ccut_wdm_fres<=0) {
    cout<<"cwb_pereport : upTDF  parameter non valid : must be power of 2 : "<<gOPT.ccut_wdm_fres<<endl;
    exit(1);
  }
  if(gOPT.ccut_bchirp < 0) {
    cout<<"cwb_pereport - Error :.ccut_bchirp must be >= 0"<<endl;
    exit(1);
  }
  if(gOPT.ccut_bchirp > gOPT.ccut_uchirp) {
    cout<<"cwb_pereport - Error :.ccut_bchirp must be < ccut_uchirp"<<endl;
    exit(1);
  }
  if(gOPT.ccut_rtime < 0) {
    cout<<"cwb_pereport - Error :.ccut_rtime >=0"<<endl;
    exit(1);
  }
  if(gOPT.ccut_ltime < gOPT.ccut_rtime) {
    cout<<"cwb_pereport - Error :.ccut_ltime must be >.ccut_rtime"<<endl;
    exit(1);
  }
}

wavearray<double> GetRCut(wavearray<double>* ts, wavearray<double>* aux) {  
// apply residual cut

  // produce the TF aux map:
  WSeries<double> WAUX;           		// TF map container
  WAUX.Forward(*aux, *gWDM);         		// apply the WDM to the auxiliary time series 

  int layers = WAUX.maxLayer()+1;  		// numbers of frequency bins (first & last bins have df/2)
  int slices = WAUX.sizeZero();    		// number of time bins

  int esize = slices*layers;
  float* energy = new float[esize];

  // set to 0 the phase 00/90 amplitudes outside the chirp band
  float etot=0;
  for(int i=0;i<slices;i++) {			// loop over time bins
    for(int j=0;j<layers;j++) {			// loop over frequency bins
      float A00 = WAUX.getSample(i,j+0.01);     // get phase 00 amplitude
      float A90 = WAUX.getSample(i,-(j+0.01));  // get phase 90 amplitude
      energy[j+i*layers] = A00*A00+A90*A90;
      etot+=energy[j+i*layers];
    }
  }

  // produce the TF ts map:
  WSeries<double> WTS;          		// TF map container
  WTS.Forward(*ts, *gWDM);         		// apply the WDM to the time series 

  int *index = new int[esize];
  TMath::Sort(esize,energy,index,true);
  float ecum=0;
  for(int k=0;k<esize;k++) {
    int m = index[k];
    ecum+=energy[m];
    if(ecum <= etot*gOPT.rcut_thr) {
      //cout << k << "/" << esize << " " << energy[m] << endl;
    } else {
      int j = m%layers;
      int i = (m-j)/layers;
      WTS.putSample(0,i,j+0.01);
      WTS.putSample(0,i,-(j+0.01));
    }
  }

  // return to time domain
  WSeries<double> xWTS = WTS;
  WTS.Inverse();
  xWTS.Inverse(-2);
  WTS += xWTS;
  WTS *= 0.5;

  return WTS;
}

void GetRCutParms(TString options) {

  // initialize with the default values
  gOPT.rcut_wdm_fres   = RCUT_WDM_FRES;
  gOPT.rcut_thr        = RCUT_THR;

  // read the user parameter rcut
  TObjArray* token = TString(options).Tokenize(TString(':'));
  if(token->GetEntries()!=2) {
    cout << endl;
    cout << "cwb_pereport - Error : rcut format - must be: wdm_fres:thr" <<  endl;
    cout << "               Default setup is: 64:1.0" << endl;
    cout << "               To setup a default value use *" << endl;
    cout << "               Example: *:0.9" << endl;
    cout << endl;
    exit(1);
  } 
  for(int j=0;j<token->GetEntries();j++) {
    TObjString* tok = (TObjString*)token->At(j);
    TString stok = tok->GetString();

    if(stok=="*") continue;

    if(j==0) {
      TString wdm_fres=stok;
      if(wdm_fres.IsDigit()) gOPT.rcut_wdm_fres=wdm_fres.Atoi();
      else {cout<<"cwb_pereport - Error : wdm_fres is not integer"<<endl;exit(1);}                   
    }
    if(j==1) {
      TString thr=stok;
      if(thr.IsFloat()) gOPT.rcut_thr=thr.Atof();
      else {cout<<"cwb_pereport - Error : thr is not float"<<endl;exit(1);}                   
    }
  }

  int x=gOPT.rcut_wdm_fres;        // must be power of 2
  if(!((x != 0) && ((x & (~x + 1)) == x)) || gOPT.rcut_wdm_fres<=0) {
    cout<<"cwb_pereport : upTDF  parameter non valid : must be power of 2 : "<<gOPT.rcut_wdm_fres<<endl;
    exit(1);
  }
  if(gOPT.rcut_thr < 0 || gOPT.rcut_thr > 1) {
    cout<<"cwb_pereport - Error :rcut_thr must be >= 0 && <=1"<<endl;
    exit(1);
  }
}

void GetLoglParms(TString options) {

  // initialize with the default values
  gOPT.logl_enabled    = LOGL_ENABLED;
  gOPT.logl_flow       = LOGL_FLOW;
  gOPT.logl_fhigh      = LOGL_FHIGH;
  gOPT.logl_arthr      = LOGL_ARTHR;
  gOPT.logl_ifo_mask   = LOGL_IFO_MASK;
  gOPT.logl_resample   = LOGL_RESAMPLE;

  // read the user parameter logl
  TObjArray* token = TString(options).Tokenize(TString(':'));
  if(token->GetEntries()!=6) {
    cout << endl;
    cout << "cwb_pereport - Error : logl format - must be: enabled:flow:fhigh:arthr:ifo_mask:resample" <<  endl;
    cout << "               Default setup is: false:0:8192:0.001:0x11111111:1" << endl;
    cout << "               Example: true:20:100:*:0x100:*  (select only the third ifo)" << endl;
    cout << endl;
    exit(1);
  } 
  for(int j=0;j<token->GetEntries();j++) {
    TObjString* tok = (TObjString*)token->At(j);
    TString stok = tok->GetString();

    if(stok=="*") continue;

    if(j==0) {
      TString enabled=stok;
      enabled.Remove(0,enabled.Last('=')+1);
      if(enabled=="true")  gOPT.logl_enabled=true;
      if(enabled=="false") gOPT.logl_enabled=false;
    }
    if(j==1) {
      TString flow=stok;
      if(flow.IsFloat()) gOPT.logl_flow=flow.Atof();
      else {cout<<"cwb_pereport - Error : logl flow is not float"<<endl;exit(1);}                   
    }
    if(j==2) {
      TString fhigh=stok;
      if(fhigh.IsFloat()) gOPT.logl_fhigh=fhigh.Atof();
      else {cout<<"cwb_pereport - Error : logl fhigh is not float"<<endl;exit(1);}                   
    }
    if(j==3) {
      TString arthr=stok;
      if(arthr.IsFloat()) gOPT.logl_arthr=arthr.Atof();
      else {cout<<"cwb_pereport - Error : logl arthr is not float"<<endl;exit(1);}                   
    }
    if(j==4) {
      TString ifo_mask=stok;
      if(ifo_mask.IsHex()) {
        std::string sifo_mask = ifo_mask.Data();
        gOPT.logl_ifo_mask=std::stoi(sifo_mask, 0, 16);      
      } else {cout<<"cwb_pereport - Error : logl ifo_mask is not hex"<<endl;exit(1);}                   
    }
    if(j==5) {
      TString resample=stok;
      if(resample.IsDigit()) gOPT.logl_resample=resample.Atoi();
      else {cout<<"cwb_pereport - Error : logl resample is not integer"<<endl;exit(1);}                   
    }
  }

  if(gOPT.logl_flow > gOPT.logl_fhigh) {
    cout<<"cwb_pereport - Error : logl_flow must be < logl_fhigh "<<endl;
    exit(1);
  }
  if(gOPT.logl_arthr>=1) {
    cout<<"cwb_pereport - Error : logl_arthr must be < 1"<<endl;
    exit(1);
  }
}

void GetRstatParms(TString options) {

  // initialize with the default values
  gOPT.rstat_enabled    = RSTAT_ENABLED;
  gOPT.rstat_type       = RSTAT_TYPE;
  gOPT.rstat_rtrials    = RSTAT_RTRIALS;
  gOPT.rstat_jtrials    = RSTAT_JTRIALS;

  // read the user parameter rstat
  TObjArray* token = TString(options).Tokenize(TString(':'));
  if(token->GetEntries()!=4) {
    cout << endl;
    cout << "cwb_pereport - Error : rstat format - must be: true:median:rtrials:jtrials" <<  endl;
    cout << "               Default setup is: false:median:0:0" << endl;
    cout << "               Example: true:median:1000:100 " << endl;
    cout << endl;
    exit(1);
  }
  for(int j=0;j<token->GetEntries();j++) {
    TObjString* tok = (TObjString*)token->At(j);
    TString stok = tok->GetString();

    if(stok=="*") continue;

    if(j==0) {
      TString enabled=stok;
      enabled.Remove(0,enabled.Last('=')+1);
      if(enabled=="true")  gOPT.rstat_enabled=true;
      if(enabled=="false") gOPT.rstat_enabled=false;
    }
    if(j==1) {
      TString type=stok;
      type.Remove(0,type.Last('=')+1);
      if((type!="rstat1")&&(type!="median")) {
        cout<<"cwb_pereport - Error : available options for rstat_type are: rstat1/median "<<endl;
        exit(1);
      } else gOPT.rstat_type=type;
    }
    if(j==2) {
      TString rtrials=stok;
      if(rtrials.IsDigit()) gOPT.rstat_rtrials=rtrials.Atoi();
      else {cout<<"cwb_pereport - Error : rstat rtrials is not integer"<<endl;exit(1);}
    }
    if(j==3) {
      TString jtrials=stok;
      if(jtrials.IsDigit()) gOPT.rstat_jtrials=jtrials.Atoi();
      else {cout<<"cwb_pereport - Error : rstat jtrials is not integer"<<endl;exit(1);}
    }
  }

  if(gOPT.rstat_enabled==false) {	// when rstat is disabled set default
    gOPT.rstat_rtrials=RSTAT_RTRIALS;
    gOPT.rstat_jtrials=RSTAT_JTRIALS;
  } else {
    if(gOPT.rstat_rtrials < 0) {
      cout<<"cwb_pereport - Error : rstat_rtrials must be >=0 "<<endl;
      exit(1);
    }
    if(gOPT.rstat_jtrials < 0) {
      cout<<"cwb_pereport - Error : rstat_jtrials must be >=0 "<<endl;
      exit(1);
    }
  }
}

int ComputeWaveformFCR(int nIFO) {

#define MAX_TRIALS	30
#define MAX_ECOV	1.0

  // The waveforms are aligned and, if requested, sync in phase and time
  // the reference waveform is REC, this choice is necessary to fix the original onsource waveform
  int selected=SyncWaveforms(nIFO, "waveform");	
  cout << endl << "ComputeWaveformFCR - selected events : " << selected << endl << endl;
  if(selected==0) return 0;

  double cov50[NIFO_MAX],cov90[NIFO_MAX],cov99[NIFO_MAX];		// local coverage

  double fcov50[NIFO_MAX],fcov90[NIFO_MAX],fcov99[NIFO_MAX];		// full coverage
  double bcov50[NIFO_MAX],bcov90[NIFO_MAX],bcov99[NIFO_MAX];		// bottom coverage
  double ucov50[NIFO_MAX],ucov90[NIFO_MAX],ucov99[NIFO_MAX];		// up coverage

  for(int n=0;n<nIFO;n++) {bcov50[n]=0.;   bcov90[n]=0.;   bcov99[n]=0.;  }
  for(int n=0;n<nIFO;n++) { cov50[n]=50.;   cov90[n]=90.;   cov99[n]=99.; }
  for(int n=0;n<nIFO;n++) {ucov50[n]=100.; ucov90[n]=100.; ucov99[n]=100.;}

  int  trials=0;
  bool done=false;
  while(!done && trials<MAX_TRIALS) { 		// dichotomic search

    trials++;

    ComputeWaveformCR(nIFO, selected, cov50, cov90, cov99);

    GetFullCoverage(nIFO, selected, fcov50, fcov90, fcov99);
    for(int n=0;n<nIFO;n++) {
      cout.precision(6);
      cout << endl;
      cout << "Trials = " << trials << "\tIFO " << gIFO[n] << "\tFull/Local Coverage at 50% = " << (int)fcov50[n] << " / " << cov50[n] << endl;
      cout << "Trials = " << trials << "\tIFO " << gIFO[n] << "\tFull/Local Coverage at 90% = " << (int)fcov90[n] << " / " << cov90[n] << endl;
      cout << "Trials = " << trials << "\tIFO " << gIFO[n] << "\tFull/Local Coverage at 99% = " << (int)fcov99[n] << " / " << cov99[n] << endl;
      cout << endl;
    }

    done=true;
    for(int n=0;n<nIFO;n++) {                      
      if(fabs(fcov50[n]-50.)>MAX_ECOV) {
        if(fcov50[n]<50.) {bcov50[n]=cov50[n]; cov50[n] += (ucov50[n]- cov50[n])/2.;}
        else              {ucov50[n]=cov50[n]; cov50[n] -= ( cov50[n]-bcov50[n])/2.;}
        done=false;
      }
      if(fabs(fcov90[n]-90.)>MAX_ECOV) {
        if(fcov90[n]<90.) {bcov90[n]=cov90[n]; cov90[n] += (ucov90[n]- cov90[n])/2.;}
        else              {ucov90[n]=cov90[n]; cov90[n] -= ( cov90[n]-bcov90[n])/2.;}
        done=false;
      }
      if(fabs(fcov99[n]-99.)>MAX_ECOV) {
        if(fcov99[n]<99.) {bcov99[n]=cov99[n]; cov99[n] += (ucov99[n]- cov99[n])/2.;}
        else              {ucov99[n]=cov99[n]; cov99[n] -= ( cov99[n]-bcov99[n])/2.;}
        done=false;
      }
    }
  }

  // dump coverage
  TString parm="";
  DumpStatistics("full_coverage.txt", "", false);
  for(int n=0;n<nIFO;n++) {                      
    parm = TString::Format("%s : Trials = %d/%d \tIFO %s \tFull/Local Coverage at 50% = %.4f / %.4f", 
           gOPT.gw_name.Data(), trials, MAX_TRIALS, gIFO[n].Data(), fcov50[n], cov50[n]);
    DumpStatistics("full_coverage.txt", parm, true);
    parm = TString::Format("%s : Trials = %d/%d \tIFO %s \tFull/Local Coverage at 90% = %.4f / %.4f", 
           gOPT.gw_name.Data(), trials, MAX_TRIALS, gIFO[n].Data(), fcov90[n], cov90[n]);
    DumpStatistics("full_coverage.txt", parm, true);
    parm = TString::Format("%s : Trials = %d/%d \tIFO %s \tFull/Local Coverage at 99% = %.4f / %.4f", 
           gOPT.gw_name.Data(), trials, MAX_TRIALS, gIFO[n].Data(), fcov99[n], cov99[n]);
    DumpStatistics("full_coverage.txt", parm, true);
    DumpStatistics("full_coverage.txt", "", true);
  }

  return selected;
}

void ComputeWaveformCR(int nIFO, int selected, double* cov50, double* cov90, double* cov99) {

  // clean vectors
  while(!vMED.empty()) vMED.pop_back();
  vMED.clear(); std::vector<wavearray<double> >().swap(vMED);
  while(!vL50.empty()) vL50.pop_back();
  vL50.clear(); std::vector<wavearray<double> >().swap(vL50);
  while(!vU50.empty()) vU50.pop_back();
  vU50.clear(); std::vector<wavearray<double> >().swap(vU50);
  while(!vL90.empty()) vL90.pop_back();
  vL90.clear(); std::vector<wavearray<double> >().swap(vL90);
  while(!vU90.empty()) vU90.pop_back();
  vU90.clear(); std::vector<wavearray<double> >().swap(vU90);
  while(!vL99.empty()) vL99.pop_back();
  vL99.clear(); std::vector<wavearray<double> >().swap(vL99);
  while(!vU99.empty()) vU99.pop_back();
  vU99.clear(); std::vector<wavearray<double> >().swap(vU99);

  // compute vMED, vL50, vU50, vL90, vU90, vL99, vU99
  wavearray<double> wmed[NIFO_MAX];                
  wavearray<double> wl50[NIFO_MAX];                
  wavearray<double> wu50[NIFO_MAX];                
  wavearray<double> wl90[NIFO_MAX];                
  wavearray<double> wu90[NIFO_MAX];                
  wavearray<double> wl99[NIFO_MAX];                
  wavearray<double> wu99[NIFO_MAX];                
  for(int n=0;n<nIFO;n++) {                      
    wmed[n] = wREC[0][n];  wmed[n] = 0;  
    wl50[n] = wREC[0][n];  wl50[n] = 0;  
    wu50[n] = wREC[0][n];  wu50[n] = 0;  
    wl90[n] = wREC[0][n];  wl90[n] = 0;  
    wu90[n] = wREC[0][n];  wu90[n] = 0;  
    wl99[n] = wREC[0][n];  wl99[n] = 0;  
    wu99[n] = wREC[0][n];  wu99[n] = 0;  

    int nentries = selected;			// number of detected events in the offsource
    int *index = new int[nentries];
    double *value = new double[nentries];
    for(int j=0;j<wREC[0][n].size();j++) {
      if(gOPT.residuals) {
        int k=0; for(int i=0;i<gEVENTS;i++) if(wNUL[i][n].size()) value[k++] = wNUL[i][n][j];  // select detected events
      } else {
        int k=0; for(int i=0;i<gEVENTS;i++) if(wREC[i][n].size()) value[k++] = wREC[i][n][j];  // select detected events
      }
      TMath::Sort(nentries,value,index,false);

      int imed = (nentries*50.)/100.; if(imed>=nentries) imed=nentries-1;
      wmed[n][j] = value[index[imed]];

      int il50 = (nentries*(50.-cov50[n]/2.))/100.; if(il50>=nentries) il50=nentries-1;
      int iu50 = (nentries*(50.+cov50[n]/2.))/100.; if(iu50>=nentries) iu50=nentries-1;
      int il90 = (nentries*(50.-cov90[n]/2.))/100.; if(il90>=nentries) il90=nentries-1;
      int iu90 = (nentries*(50.+cov90[n]/2.))/100.; if(iu90>=nentries) iu90=nentries-1;
      int il99 = (nentries*(50.-cov99[n]/2.))/100.; if(il99>=nentries) il99=nentries-1;
      int iu99 = (nentries*(50.+cov99[n]/2.))/100.; if(iu99>=nentries) iu99=nentries-1;

      double med = wmed[n][j];
      double l50 = value[index[il50]];
      double u50 = value[index[iu50]];
      double l90 = value[index[il90]];
      double u90 = value[index[iu90]];
      double l99 = value[index[il99]];
      double u99 = value[index[iu99]];

      bool check=true;
      if(!(l50<=u50)) check=false;  
      if(!(l90<=u90)) check=false;  
      if(!(l99<=u99)) check=false;  
      if(!(med>=l50 && med<=u50)) check=false;  
      if(!(med>=l90 && med<=u90)) check=false;  
      if(!(med>=l99 && med<=u99)) check=false;  
      if(!check) {cout << "ComputeWaveformFCR : standard wrong median, lower, upper bound !!! " << l50 << " " << med << " " << u50 << endl;}  
      if(!(med>=l50 && med<=u50)) med=l50;  
      wmed[n][j]=med;

      wl50[n][j] = fabs(med-l50);
      wu50[n][j] = fabs(u50-med);
      wl90[n][j] = fabs(med-l90);
      wu90[n][j] = fabs(u90-med);
      wl99[n][j] = fabs(med-l99);
      wu99[n][j] = fabs(u99-med);
    }        
    delete [] index;
    delete [] value;

    vMED.push_back(wmed[n]);
    vL50.push_back(wl50[n]);
    vU50.push_back(wu50[n]);
    vL90.push_back(wl90[n]);
    vU90.push_back(wu90[n]);
    vL99.push_back(wl99[n]);
    vU99.push_back(wu99[n]);
  }                                                                    
}

void GetFullCoverage(int nIFO, int selected, double* fcov50, double* fcov90, double* fcov99) {

  // Compute time boundaries which contain the P energy fraction of median

  double P = gOPT.plot_efraction;
  double estart, estop;				
  for(int n=0;n<nIFO;n++) {
    double bT, eT;
    CWB::mdc::GetTimeBoundaries(vMED[n], P, bT, eT);
    if(n==0) {
      estart=bT;estop=eT;
    } else {
      if(bT<estart) estart=bT;
      if(eT>estop)  estop=eT;
    }
  }
  cout << endl;
  cout << "estart=" << estart << "\testop=" << estop << endl;
  cout << endl;

  int nstart = (estart-vREC[0].start())*vREC[0].rate();
  int nstop  = (estop -vREC[0].start())*vREC[0].rate();

  // Compute the full coverages

  for(int n=0;n<nIFO;n++) {
    int n50=0;
    int n90=0;
    int n99=0;
    for(int i=0;i<gEVENTS;i++) {
      if(!wREC[i][n].size()) continue;  // select detected events
      bool b50=false;
      bool b90=false;
      bool b99=false;
      for(int j=nstart;j<nstop;j++) {
        double vl50 = vMED[n][j]-fabs(vL50[n][j]);
        double vu50 = vMED[n][j]+fabs(vU50[n][j]);
        double vl90 = vMED[n][j]-fabs(vL90[n][j]);
        double vu90 = vMED[n][j]+fabs(vU90[n][j]);
        double vl99 = vMED[n][j]-fabs(vL99[n][j]);
        double vu99 = vMED[n][j]+fabs(vU99[n][j]);

        double value = wREC[i][n][j];
        if(value<vl50 || value>vu50) b50=true;
        if(value<vl90 || value>vu90) b90=true;
        if(value<vl99 || value>vu99) b99=true;
      }
      if(b50) n50++;
      if(b90) n90++;
      if(b99) n99++;
    }
    fcov50[n] = (100.-100.*n50/(double)selected);
    fcov90[n] = (100.-100.*n90/(double)selected);
    fcov99[n] = (100.-100.*n99/(double)selected);
  }
}

