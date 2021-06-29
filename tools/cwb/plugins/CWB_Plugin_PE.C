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


#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "mdc.hh"
#include "cwb2G.hh"
#include "watplot.hh"
#include "gwavearray.hh"
#include "network.hh"
#include <fstream>
#include <vector>

// ---------------------------------------------------------------------------------
// TO BE DONE
// ---------------------------------------------------------------------------------

// MEDIAN SKYMAP, reformat of xCED SkyMap TAB
// FIX WAVEFORM INJECTION in RedoAnalysis, must be excluded CAT2 , we can use TH1D for random selection
// ADD SVN Version to DumpUserOptions 							---> (DONE)
// FIX output format of DumpUserOptions 						---> (DONE)
// FIX xCED Spectrograms when pe_noise=0, only detected waveform is displayed 		---> (DONE)
// ADD GPS time offset to the time plots						---> (DONE)
// ADD WDM TF maps ?
// ADD descriptions for each TAB in the xCED
// SAVE injected whitened signal 							---> (DONE)
// SAVE all PE/waveform in the output ROOT file						---> (DONE)
// For simulation=0 set window +/- iwindow/2 around the event				---> (DONE)
// SAVE pe_config to root output file						
// DUMP to CED dir all PE parameters							---> (DONE)	
// ADD Checks of the input PE config parameters						---> (DONE)
// FIX pe_noise=2 when cedDump=false -> the TF is cleaned and whitened noise is deleted. 
// USE pe_ced_options to select tabs to be included into xCED				---> (DONE)
// IMPLEMENT command : cwb_condor mtpe jobID						---> (DONE)
// FIX iSNR,oSNR for non simulation PE : the output SNRnet plot is zero
// add check on ISTAGE, PE can only be executed with FULL stage 			---> (DONE)
// add TF PCA likelihood/null								---> (DONE)

// ---------------------------------------------------------------------------------
// HOW TO CONFIGURE PLUGIN
// the following is an example : must be included in the config/user_parameters.C file
// see the DumpUserOptions function for the full description of parameters
// ---------------------------------------------------------------------------------
/*
  plugin = TMacro(gSystem->ExpandPathName("$HOME_CWB/plugins/CWB_Plugin_PE.C"));        // Macro source

  TString optpe = "";                  		// NOTE : add space at the end of each line
  optpe += "pe_retry=100 ";			// number of trials
  optpe += "pe_gps=1167559936 ";		// only gps +/- iwindow is analyzed
  optpe += "pe_noise=0 ";              		// signal used for trials : 0/1/2 reconstructed/injected/(original=reconstructed+noise)
  optpe += "pe_signal=0 ";             		// waveform is injected in the whitened HoT and apply Coherent+SuperCluster+Likelihood stage
  optpe += "pe_amp_cal_err=0.1 ";		// max percentage of amplitude miscalibration  (uniform in -0.9,+1.1) -> det 1
  optpe += "pe_phs_cal_err=10 ";		// max phase (degrees) of phase miscalibration (uniform in -10,+10)   -> det 1
  optpe += "pe_amp_cal_err=0.1 ";		// ...                                                                -> det 2
  optpe += "pe_phs_cal_err=10 ";		// ...                                                                -> det 2     
  optpe += "pe_ced_dump=-1 ";			// dump CED at gLRETRY=PE_CED_DUMP
  optpe += "pe_skymask=0 ";			// disable -> search over all sky

  //optpe += "pe_multitask=true ";		// enable/disable multitask (only for 1 job and simulation=0/4, dag is generated with cwb_condor mtpe jobID)
						// for simulation=0 user must define : optpe += "pe_retry=100 ";
						// for simulation=4 user must define : nfactor=xx, factor[0]=yy : where xx>0 and yy>0          	

						// CED options (only if cedDump=true)
  optpe += "pe_ced_enable=tfmap ";
  optpe += "pe_ced_enable=rdr ";
  optpe += "pe_ced_enable=skymap ";
  optpe += "pe_ced_enable=rec ";
  optpe += "pe_ced_enable=inj ";
  optpe += "pe_ced_disable=rinj ";
  optpe += "pe_ced_enable=cm ";
  optpe += "pe_ced_enable=distr ";
  optpe += "pe_ced_enable=null ";
  optpe += "pe_ced_disable=pca ";

						// OUTPUT options (only if cedDump=false)
  optpe += "pe_output_enable=inj ";         	// save injection to the output root file
  optpe += "pe_output_enable=rec ";         	// save reconstructed waveform to the output root file
  optpe += "pe_output_enable=med ";         	// save median to the output root file
  optpe += "pe_output_enable=p50 ";         	// save percentile 50 to the output root file
  optpe += "pe_output_enable=p90 ";         	// save percentile 90 to the output root file
  optpe += "pe_output_enable=avr ";         	// save averaged waveform to the output root file
  optpe += "pe_output_enable=rms ";         	// save RMS to the output root file

  strcpy(parPlugin,optpe.Data());      		// set PE plugin parameters
  strcpy(comment,"pe configuration example");

*/

// ---------------------------------------------------------------------------------
// DEFINES
// ---------------------------------------------------------------------------------

#define nTRIALS 		100		// number of trials for noise statistic (GetNoiseStatistic)

#define PE_INDEX               "$HOME_CWB/www/pe_index_modern.html"	// html template used for pe index report page

//#define PLOT_TYPE		"root"
#define PLOT_TYPE		"png"
//#define SAVE_WAVE2SPARSE_PLOT
//#define SAVE_REDUCED_INJ	

#define PLOT_MEDIAN				// plot median,lower50 boud,upper50 bound,lower90 boud,upper90 bound
						// if commented we plot average,rms,2*rms 

//#define SAVE_SKYPROB				// save gSKYPROB
#define SKYMASK_RADIUS		0.1		// used to define skymask

#define MAX_TRIALS		1000

#define EFRACTION_CUT		0.98		// energy threshold for time and frequency cuts
//#define SET_WAVEFORM_CUTS			// enable time and frequency cuts

// default pe user options
#define PE_TRIALS 		100		// number of trials
#define PE_SIGNAL		0		// signal used for trials : 0/1/2 reconstructed/injected/(original=reconstructed+noise)
                                                // note : option 1 can be used only in simulation mode
#define PE_NOISE		1		// noise used for trials :  0/1/2 
                                                // 0 - waveform is injected in the whitened HoT and apply Coherent+SuperCluster+Likelihood stages
                                                // 1 - add gaussian noise to likelihood sparse map and apply Likelihood stage
                                                // 2 - add whitened HoT to likelihood sparse map and apply Likelihood stage
#define PE_AMP_CAL_ERR 		0		// max percentage of amplitude miscalibration  if>0 (uniform in -max,+max) else (gaus(max))
#define PE_PHS_CAL_ERR 		0		// max phase (degrees) of phase miscalibration if>0 (uniform in -max,+max) else (gaus(max))
#define PE_SKYMASK 		0		// skymask used for trials : 0(def)/1/2 
                                                // 0 : disable -> search over all sky
                                                // 1 : skymask with SKYMASK_RADIUS and source selected according the reconstructed skymap probability 
                                                // 2 : skymask select only the sky position where the waveform is reconstructed
#define PE_SEED 		0		// 0 : seed used  by PE for ramdom generation
#define PE_GPS 			0		// default 0 - if >0 only gps +/- iwindow is analyzed

#define PE_MULTITASK		false		// if enabled, PE trials are executed by different jobs in multitask : only simulation=0 and single event (gps>0) 

#define NOISE_SIGMA		1.0		// noise sigma used in AddNoiseAndCalErrToSparse 

#define PE_CED_DUMP 		-1    		// dump CED at gLRETRY=PE_CED_DUMP

#define PE_CED_TFMAP 		true		// Shows the Time-Frequency Maps -> Spectrograms, Scalograms, TF Like,Null  
#define PE_CED_RDR 		true		// Shows Reconstructed Detector Response tab 
#define PE_CED_PSD 		true		// Shows Power Spectral Density tab 
#define PE_CED_SKYMAP 		true		// Shows SkyMaps
#define PE_CED_REC 		true		// Shows Reconstructed Waveforms/Instantaneous-Frequency with errors
#define PE_CED_INJ 		true		// Shows Injected vs Reconstructed Waveforms/Instantaneous-Frequency  with errors
#define PE_CED_rINJ 		false		// Shows Reduced Injected vs Reconstructed Waveforms with errors 
#define PE_CED_CM		true		// Shows Chirp Mass Value/Error Distributions
#define PE_CED_DISTR		false		// Shows Residuals Distributions
#define PE_CED_NULL 		true		// Shows Null Pixels Distributions
#define PE_CED_PCA 		false		// Shows PCA likelihood TF plot

#define PE_OUTPUT_INJ		false		// save injected      into the output root file
#define PE_OUTPUT_REC		false		// save reconstructed into the output root file
#define PE_OUTPUT_WHT		false		// save whitened data into the output root file
#define PE_OUTPUT_DAT		false		// save rec+nois data into the output root file
#define PE_OUTPUT_MED		false		// save median        into the output root file
#define PE_OUTPUT_P50		false		// save percentile 50 into the output root file
#define PE_OUTPUT_P90		false		// save percentile 90 into the output root file
#define PE_OUTPUT_AVR		false		// save averaged      into the output root file
#define PE_OUTPUT_RMS		false		// save rms           into the output root file

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void ClearVectors();
void ClearWaveforms(detector* ifo);

double GetSparseMapData(SSeries<double>* SS, bool phase, int index);

void SetOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, bool dump_pe);
void DumpOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor);
TString DumpCED(network* NET, netevent* &EVT, CWB::config* cfg, double factor);
void CreateIndexHTML(TString dirCED, int nIFO, TString* ifo, bool sim=false);
void WriteBodyHTML(TString html_template, TString html_tag_beg, TString html_tag_end, ofstream* out, int nIFO=1, TString* ifo=NULL);

void Wave2Sparse(network* NET, CWB::config* cfg, char type);
void AddNoiseAndCalErrToSparse(network* NET, CWB::config* cfg, char type);
int  RedoAnalysis(TFile* jfile, CWB::config* cfg, network* NET);
void ReplaceSuperclusterData(TFile*& jfile, CWB::config* cfg, network* NET, double gps=0);

void GetChirpMassStatistic(std::vector<double>*  vCHIRP);
void GetFactorsStatistic(int nIFO);
void GetSNRnetStatistic(int nIFO);
void GetNullPixels(std::vector<double>* aNUL, std::vector<double>* ANUL, 
                      network* NET, CWB::config* cfg, int lag, int id);
void GetNoisePixels(std::vector<double>* aNSE, std::vector<double>* ANSE, 
                      network* NET, CWB::config* cfg, int lag, int id);
std::vector<int> ComputeStatisticalErrors(network* NET, CWB::config* cfg, int ID);
void SaveSkyProb(network* NET, CWB::config* cfg, int id);
void SetEventWindow(CWB::config* cfg, double gps);
skymap GetSkyProb(network* NET, int id);
skymap GetMedianSkyProb(network* NET);
void DumpSkyProb(skymap* skyprob, network* NET, netevent* &EVT, TString odir);

void SetWaveformCuts(wavearray<double>* x, double bT, double eT, double bF, double eF);
wavearray<double> GetCutWaveform(wavearray<double> x, double bT, double eT, double bF, double eF);

double GetCentralTime(wavearray<double> x);
double GetTimeBoundaries(wavearray<double> x, double P, double& bT, double& eT);
double GetFrequencyBoundaries(wavearray<double> x, double P, double& bF, double& eF);
double GetCentralFrequency(wavearray<double> x);

std::vector<netpixel> DoPCA(network* NET, CWB::config* cfg, int lag, int id);
void ResetPCA(network* NET, CWB::config* cfg, netcluster* pwc, std::vector<netpixel>* vPIX, int ID);

std::vector<wavearray<double> > GetPCAWaveform(network* NET, CWB::config* cfg, int lag, int id);
std::vector<wavearray<double> > GetInjWaveform(network* NET, netevent* EVT, int id, double factor);
std::vector<wavearray<double> > GetRecWaveform(network* NET, netevent* EVT, int id);
std::vector<wavearray<double> > GetWhitenedData(network* NET, CWB::config* cfg);
std::vector<wavearray<double> > GetWaveform(network* NET, int lag, int id, char type, bool shift=true);
std::vector<wavearray<double> > GetSigWaveform(network* NET, CWB::config* cfg, int lag, int id);

void SetSkyMask(network* NET, CWB::config* cfg, double theta, double phi, double radius);

wavearray<double> GetDifWaveform(wavearray<double>* wf1, wavearray<double>* wf2);
wavearray<double> GetAlignedWaveform(wavearray<double>* wf1, wavearray<double>* wf2);
wavearray<double> GetWaveformEnvelope(wavearray<double>* wf);
wavearray<double> AddWaveforms(wavearray<double>* wf1, wavearray<double>* wf2);

double FittingFactor(wavearray<double>* wf1, wavearray<double>* wf2);

void PlotSparse(int ifoID, network* NET, CWB::config* cfg, int ID, wavearray<double>* wf);

void PlotWaveform(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wf1, 
                  wavearray<double>* wf2, wavearray<double>* wf3, wavearray<double>* wref, 
                  bool fft=false, TString pdir="", double P=0.99, EColor col1=kRed, EColor col2=kBlue, EColor col3=(EColor)418);
void PlotWaveformErrors(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wrec,
                  wavearray<double>* wavr, wavearray<double>* werr, wavearray<double>* wref, TString pdir="", double P=0.99);
void PlotWaveformAsymmErrors(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wrec,
                  wavearray<double>* wmed, wavearray<double>* wl50, wavearray<double>* wu50, 
                  wavearray<double>* wl90, wavearray<double>* wu90, wavearray<double>* wref, TString pdir, double P, bool freq=false);
void PlotFrequencyErrors(TString ofname, TString title, CWB::config* cfg, wavearray<double>* frec,
                  wavearray<double>* favr, wavearray<double>* ferr, wavearray<double>* wref, TString pdir, double P);
void PlotWaveforms(network* NET, CWB::config* cfg, int ID, TString pdir="");
void PlotResiduals(int ID, TString pdir="", int sim=true, char type='w');
void PlotSNRnet(int nIFO, TString pdir, int sim=true);
void PlotChirpMass(int gtype, TString pdir, int sim=true);
void PlotFactors(int gtype, int nIFO, TString pdir);
void GetNullStatistic(std::vector<double>* vNUL, std::vector<double>* vNSE, int ifoId, TString ifoName, TString gtype, TString pdir="");
void PlotTimeFrequencyPCA(network* NET, netevent* EVT, CWB::config* cfg, int id, int lag, TString pdir);
void PlotSpectrogram(TString type, network* NET, netevent* &EVT, CWB::config* cfg, TString pdir);

void DumpRecWavePE(network* NET, TString pdir=""); 
void DumpInjWavePE(network* NET, TString pdir=""); 

void LoadFromMultiTaskJobOutput(int runID, CWB::config* cfg);

// ---------------------------------------------------------------------------------
// USER CONFIG OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  int    trials;
  int    signal;
  int    noise;
  float  amp_cal_err[NIFO_MAX];
  float  phs_cal_err[NIFO_MAX];
  int    id;
  int    skymask;
  int    seed;		// default 0
  double gps;		// default 0 - if >0 only gps +/- iwindow is analyzed

  bool   multitask;	// if enabled, PE trials are executed by different jobs in multitask : only simulation=0 and single event (gps>0) 

  int    ced_dump;    	// dump CED at gLRETRY=ced_dump (-1 = disable)

  bool   ced_tfmap;	// Shows the Time-Frequency Maps -> Spectrograms, Scalograms, TF Like,Null
  bool   ced_rdr;	// Shows Reconstructed Detector Response
  bool   ced_psd;	// Show  Power Spectral Density
  bool   ced_skymap;	// Shows SkyMaps
  bool   ced_rec;	// Shows Reconstructed Waveforms/Instantaneous-Frequency with errors
  bool   ced_inj;	// Shows Injected vs Reconstructed Waveforms/Instantaneous-Frequency  with errors
  bool   ced_rinj;	// Shows Reduced Injected vs Reconstructed Waveforms with errors
  bool   ced_cm;	// Shows Chirp Mass Value/Error Distributions
  bool   ced_distr;	// Shows Residuals Distributions
  bool   ced_null;	// Shows Null Pixels Distributions
  bool   ced_pca;	// Shows PCA likelihood TF plot 

  bool   output_inj;	// save injected      into the output root file
  bool   output_rec;	// save reconstructed into the output root file
  bool   output_wht;	// save whitened data into the output root file
  bool   output_dat;	// save rec+nois data into the output root file
  bool   output_med;	// save median        into the output root file
  bool   output_p50;	// save percentile 50 into the output root file
  bool   output_p90;	// save percentile 90 into the output root file
  bool   output_avr;	// save averaged      into the output root file
  bool   output_rms;	// save rms           into the output root file
};

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);
void DumpUserOptions(TString odir, CWB::config* cfg);

// ----------------------------------------------------
// ROOT Output PE Parameters
// ----------------------------------------------------

struct PE {				// structure for output for estimated parameters

  int   trials;				// number of effective trials
  float erR[11];                  	// probability distribution of residuals
  float erF[11];                  	// probability distribution of frequency residuals
  float nstat[2*7*NIFO_MAX];      	// null pixel statistic, for each detector 
					// 0->7  : Null 00
					// 7->13 : Null 90
     					// gPE.stat[0] = Null  Pixels Mean
     					// gPE.stat[1] = Null  Pixels RMS 
     					// gPE.stat[2] = Noise Pixels Mean
     					// gPE.stat[3] = Noise Pixels RMS 
     					// gPE.stat[4] = Noise Pixels Chi2
     					// gPE.stat[5] = KolmogorovTest   
     					// gPE.stat[6] = AndersonDarlingTest 
  float snet[2];			// SNRnet statistic, 0 -> avr, 1 -> rms  
  float ff[2];				// Fitting Factor statistic, 0 -> avr, 1 -> rms  
  float of[2];				// Overlap Factor statistic, 0 -> avr, 1 -> rms  
  float mch[2];				// chirp mass statistic, 0 -> avr, 1 -> rms  

  wavearray<double>* wINJ[NIFO_MAX];
  wavearray<double>* wREC[NIFO_MAX];
  wavearray<double>* wWHT[NIFO_MAX];
  wavearray<double>* wDAT[NIFO_MAX];
  wavearray<double>* wMED[NIFO_MAX];
  wavearray<double>* wL50[NIFO_MAX]; 
  wavearray<double>* wU50[NIFO_MAX];
  wavearray<double>* wL90[NIFO_MAX];
  wavearray<double>* wU90[NIFO_MAX];
  wavearray<double>* wAVR[NIFO_MAX];
  wavearray<double>* wRMS[NIFO_MAX];
};

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions 	gOPT;				// global User Options
PE		gPE;				// global output for estimated parameters
TTree* 		gTREE;				// output tree file name
TString 	gOUTPUT;			// output root file name
cwb 		gCWB;
bool		gCEDDUMP;			// save user_parameters cedDump

double 		gITHETA;			// injected theta sky position (EVT->theta[1])
double 		gIPHI;				// injected phy sky position (EVT->phi[1])
double 		gOTHETA;			// reconstructed theta sky position (EVT->theta[3]) : the maximum of detection statistic
double 		gOPHI;				// reconstructed phy sky position (EVT->phi[3]) : the maximum of detection statistic

skymap 		gSKYPROB;              		// saved reconstructed probability skymap
TH1D   		gHSKYPROB;              	// used to extract random skyloc for skymask

double 		gBT, gET;             		// time ranges used for cutted injections (EFRACTION_CUT)
double 		gBF, gEF;             		// freq ranges used for cutted injections (EFRACTION_CUT)  

wavearray<double> gHOT[NIFO_MAX];               // HoT time series used to save whitened data

double		gINRATE;			// input data sampling rate
int 		gRATEANA;          		// analysis data rate

int		gMTRIAL;			// trial ID, used to process a single trial when PE is executed in multitask mode (only with simulation=0)
int		gID;				// ID of the detected event, used in multitask mode to exclude from the processing the ID!=gID

double		gSEGGPS;			// segment start time 

// ---------------------------------------------------------------------------------
// waveforms (vectors index is ifos)
// ---------------------------------------------------------------------------------

std::vector<wavearray<double> > vINJ;		// injected
std::vector<wavearray<double> > vREC;		// signal reconstructed by wave packet 
std::vector<wavearray<double> > vWHT;		// whitened data (in the vREC time range) 
std::vector<wavearray<double> > vPCA;		// reconstructed by PCA
std::vector<wavearray<double> > vDAT;		// noise+signal reconstructed by wave packet
std::vector<wavearray<double> > vNUL;		// vDAT-vREC
std::vector<wavearray<double> > cINJ;		// injected (cutted) on TF sub-space of the reconstructed waveform vREC

std::vector<wavearray<double> > vAVR;		// average reconstructed amplitude
std::vector<wavearray<double> > vRMS;		// rms 

std::vector<wavearray<double> > vMED;		// median reconstructed amplitude
std::vector<wavearray<double> > vL50;		// = vMED - 50% Lower Bound 
std::vector<wavearray<double> > vU50;		// = vMED - 50% Upper Bound
std::vector<wavearray<double> > vL90;		// = 90% Lower Bound - vMED
std::vector<wavearray<double> > vU90;		// = 90% Upper Bound - vMED

std::vector<wavearray<double> > fREC;		// instantaneous reconstructed frequency 
std::vector<wavearray<double> > fINJ;		// instantaneous cutted injected frequency
std::vector<wavearray<double> > fAVR;		// instantaneous averaged reconstructed frequency
std::vector<wavearray<double> > fRMS;		// rms of instantaneous averaged frequency 

std::vector<wavearray<double> > fMED;		// median instantaneous reconstructed frequency
std::vector<wavearray<double> > fL50;		// 50% Lower Bound 
std::vector<wavearray<double> > fU50;		// 50% Upper Bound
std::vector<wavearray<double> > fL90;		// 90% Lower Bound 
std::vector<wavearray<double> > fU90;		// 90% Upper Bound

std::vector<double > vRES;			// normalized residual energy
std::vector<double > fRES;			// normalized frequency weighted residuals

std::vector<wavearray<double> > wREC[MAX_TRIALS]; // reconstructed signal in the trials

std::vector<double>  iSNR[NIFO_MAX]; 		// input SNR^2
std::vector<double>  oSNR[NIFO_MAX]; 		// reconstructed SNR^2
std::vector<double> ioSNR[NIFO_MAX]; 		// iSNR*oSNR
std::vector<double>     vLIKELIHOOD; 		// likelihood

std::vector<double>  aNUL[NIFO_MAX]; 		// null 00 phase amplitudes
std::vector<double>  ANUL[NIFO_MAX]; 		// null 90 phase amplitudes

std::vector<double>  aNSE[NIFO_MAX]; 		// noise 00 phase amplitudes
std::vector<double>  ANSE[NIFO_MAX]; 		// noise 90 phase amplitudes

std::vector<double>  vCHIRP[6]; 		// 0 -> injected chirp mass
   				                // 1 -> reconstructed chirp mass
         				        // 2 -> reconstructed chirp mass error
            				        // 3 -> ellipticity parameter
                   				// 4 -> pixel fraction
                   				// 5 -> energy fraction

skymap 	             wSKYPROB[MAX_TRIALS]; 	// sky probability trials
skymap		     mSKYPROB;			// sky probability median

// ---------------------------------------------------------------------------------
// sparse maps
// ---------------------------------------------------------------------------------

std::vector<SSeries<double> > vSS[NIFO_MAX];    // original sparse maps
std::vector<SSeries<double> > jSS[NIFO_MAX];    // injected 
std::vector<SSeries<double> > rSS[NIFO_MAX];    // reconstructed
std::vector<SSeries<double> > dSS[NIFO_MAX];    // noise+signal reconstructed by wave packet

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Plugin used for Parameters Estimation

  if(type==CWB_PLUGIN_CONFIG) {  
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)

    if(gCWB2G->istage!=CWB_STAGE_FULL) {cout<< "CWB_Plugin_PE - Error : PE can be executed only in CWB_STAGE_FULL mode" << endl;exit(1);}

    if(cfg->pattern<=0) {
      cout << "CWB_Plugin_PE Error : PE enable only with pattern>0 !!!" << endl;
      exit(1);
    }
    cfg->outPlugin=true;  				// disable built-in output root file
    gINRATE=cfg->inRate;				// input data sampling rate
    gRATEANA=gCWB2G->rateANA;          			// analysis data rate
    gID=0;						// used in multitask mode

    ResetUserOptions(); 				// set default config options
    ReadUserOptions(cfg->parPlugin); 			// user config options : read from parPlugin
    TString cwb_inet_options=TString(gSystem->Getenv("CWB_INET_OPTIONS"));	// overwrite parPlugin
    ReadUserOptions(cwb_inet_options); 			// user config options : read from command line

    if(gOPT.ced_inj && gOPT.ced_rinj) gOPT.ced_rinj=false;

    if(gOPT.multitask) {
      if(cfg->simulation==4 && cfg->nfactor!=1) {
        cout<< "CWB_Plugin_PE - Error : when simulation=4, PE multitask can be executed only with nfactors=1" << endl;exit(1);
      }
    } 

    if(cfg->cedDump) { 					// if at least one of the output options is enabled 
							// then we enable the output root file even when cedDump=true
      if(gOPT.output_inj || gOPT.output_rec || 
         gOPT.output_med || gOPT.output_p50 || 
         gOPT.output_p90 || gOPT.output_avr ||
         gOPT.output_rms || gOPT.output_wht || gOPT.output_dat)   cfg->online=true;  		
    }

    // get gMTRIAL
    // NOTE : trial ID, used to process a single trial when PE is executed in multitask mode
    //        can be used only for simulation=0 and gOPT.gps>0 (single event) !!!
    gMTRIAL=gOPT.trials;
    if(gOPT.multitask) {
      TString gtrial=TString(gSystem->Getenv("CWB_MDC_FACTOR"));
      if(gtrial.CompareTo("")!=0) {
        if(!gtrial.IsDigit()) {cout<< "CWB_Plugin_PE - Error : when simulation=0, CWB_MDC_FACTOR must be an interger > 0" << endl;exit(1);}
        gMTRIAL=gtrial.Atoi();
        if(gMTRIAL==0) {cout<< "CWB_Plugin_PE - Error : when simulation=0, CWB_MDC_FACTOR must be an interger > 0" << endl;exit(1);}
        if(gMTRIAL!=gOPT.trials) {			// multitask jobs
          // add trial to data_label -> different output files for each trial
          sprintf(cfg->data_label,"%s_trial%d",cfg->data_label,gMTRIAL);	
          cout << "CWB_Plugin_PE - MultiTask Mode : new data_label = " << cfg->data_label << endl;
          // disable ooptions not used in multitask mode
          cfg->cedDump=false;
          gOPT.output_inj=false;
          gOPT.output_rec=false;
          gOPT.output_wht=false;
          gOPT.output_dat=false;
          gOPT.output_med=false;
          gOPT.output_p50=false;
          gOPT.output_p90=false;
          gOPT.output_avr=false;
          gOPT.output_rms=false;
        } else {						// last multitask job : collect all output from multitask jobs
          LoadFromMultiTaskJobOutput(gCWB2G->runID, cfg);	// read wREC,gSKYPROB from multitalk output root files
        }
        bool last_trial = (gMTRIAL==gOPT.trials) ? true : false;
        if(!last_trial) {gOPT.trials=1;gMTRIAL*=-1;}
      }
    }

    gCEDDUMP=cfg->cedDump;				// save user_parameter cedDump 

    gRandom->SetSeed(gOPT.seed);			// set seed for PE random generation

    SetEventWindow(cfg, gOPT.gps);			// if gps>0 only gps +/- iwindow is analyzed

    // check input PE config parameters

    if(gOPT.signal!=0 && gOPT.signal!=1 && gOPT.signal!=2) {
      cout << endl;
      cout << "CWB_Plugin_PE Error : option pe_signal not allowed : must be (0/1/2) !!!" << endl;
      cout << endl;
      exit(1);
    }

    if(cfg->simulation==0 && gOPT.signal==1) {
      cout << endl;
      cout << "CWB_Plugin_PE Error : option pe_signal=1 is allowed only in simulation mode !!!" << endl;
      cout << endl;
      exit(1);
    }

    if(cfg->simulation==0 && gOPT.skymask==3) {
      cout << endl;
      cout << "CWB_Plugin_PE Error : option pe_skymask=3 is allowed only in simulation mode !!!" << endl;
      cout << endl;
      exit(1);
    }
  }

  if(type==CWB_PLUGIN_INIT_JOB) {  
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)

    if(gOPT.gps<gCWB2G->GetSegBegin() || gOPT.gps>gCWB2G->GetSegEnd()) {
      cout.precision(14);
      cout << "CWB_Plugin_PE - Error : gps time must be within the segment interval: " 
           << gCWB2G->GetSegBegin() << "\t" << gCWB2G->GetSegEnd() << endl;exit(1);
      cout.precision(6);
    } 
  }

  if(type==CWB_PLUGIN_NETWORK) {
    PrintUserOptions(cfg);      // print config options
  }

  if(type==CWB_PLUGIN_IDATA_CONDITIONING) {
    if(NET->mdcTime.size()>1) {
      cout << endl;
      cout << "CWB_Plugin_PE Error : detected " << NET->mdcTime.size() << " injections." 
           << " Must be only one injection per segment !!!" << endl;
      cout << endl;
      exit(1);
    }
  }

  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {
    // save whitened HoT
    int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==NET->getifo(n)->Name) {ifoID=n;break;}
    gHOT[ifoID] = *x;
  }

  if(type==CWB_PLUGIN_ILIKELIHOOD) {			// INIT CWB LIKELIHOOD STAGE (once for each lag)
    NET->wfsave=true;                                   // enable waveform save

    // export gLRETRY
    if(gOPT.multitask) EXPORT(int,gLRETRY,"gLRETRY = 1;")
    else               EXPORT(int,gLRETRY,TString::Format("gLRETRY = %d;",gOPT.trials).Data())

    ReplaceSuperclusterData(jfile, cfg, NET, gOPT.gps);	// save in jfile max size supercluster
  }

  if(type==CWB_PLUGIN_XLIKELIHOOD) {			// BEFORE EVENT RECONSTRUCTION (for each event repeated gOPT.trials)

    wavearray<double> cid = NET->getwc(0)->get((char*)"ID",  0,'S',0); // get cluster ID
    if(cid.size()==0) return;
    int ID = size_t(cid.data[cid.size()-1]+0.1);

    if(gOPT.trials==0) return;      

    // import gLRETRY
    int gLRETRY=-1; IMPORT(int,gLRETRY)
    cout << endl;
    cout << "-------------------------------------------------------------------" << endl;
    cout << "-----> CWB_Plugin_PE.C -> ID : " << ID 
         << " -> gLRETRY : " << gLRETRY << "/" << gOPT.trials << endl;
    cout << "-------------------------------------------------------------------" << endl;

    // -------------------------------------------------
    // FIRST TRIAL
    // reconstruct signal with original noise (full sky)
    // -------------------------------------------------
    if(gLRETRY==gOPT.trials) {      
      ClearVectors();
      cfg->cedDump=false;
    }

    // -----------------------------------------------------------------------------------
    // AFTER FIRST TRIAL 
    // reconstruct signal with reconstructed signal + gaussian noise & (calib. errors)
    // -----------------------------------------------------------------------------------
    if(gLRETRY!=gOPT.trials) {      

      // SKYMASK
      if(gOPT.skymask==1) {
	// select random sky position from sky probability gHSKYPROB
        int isky = (int)gHSKYPROB.GetRandom();
        double th = 90-gSKYPROB.getTheta(isky);
        double ph = gSKYPROB.getPhi(isky);

        SetSkyMask(NET, cfg, th, ph, SKYMASK_RADIUS);
      }
      if(gOPT.skymask==2) {
        // select reconstructed sky position (the maximum of detection statistic : waveform is reconstructed in such position)
        SetSkyMask(NET, cfg, gOTHETA, gOPHI, SKYMASK_RADIUS);
      }
      if(gOPT.skymask==3) {
        // select injected sky position : only in simulation mode
        SetSkyMask(NET, cfg, gITHETA, gIPHI, SKYMASK_RADIUS);
      }

      // ADD NOISE (& ADD CALIBRATION ERRORS)
      // add noise to the sparse map of the reconstructed event
      // Note : calibration errors are applied only if gOPT.noise is enabled !!!
      char jtype='r';
      if(gOPT.signal==0) jtype='r';	// use reconstructed waveform for trials
      if(gOPT.signal==1) jtype='j';	// use injected waveform for trials
      if(gOPT.signal==2) jtype='v';	// use data for trials
      if(gOPT.noise) AddNoiseAndCalErrToSparse(NET, cfg, jtype);	

      // DUMP CED
      cfg->cedDump=false;
      if((gOPT.ced_dump>=0) && (gLRETRY==gOPT.ced_dump)) cfg->cedDump=true;
    }

    // -------------------------------------------------
    // LAST TRIAL
    // reconstruct signal with original noise (full sky)
    // -------------------------------------------------
    if((gLRETRY==0)&&(gMTRIAL==gOPT.trials)) {	// use original injection setup                      		

      // SKYMASK
      // select detected sky position (waveform is reconstructed in such poosition)
      //SetSkyMask(NET, cfg, gOTHETA, gOPHI, SKYMASK_RADIUS);
      // restore full sky search for the last trial
      SetSkyMask(NET, cfg, 0, 0, 360);	// full sky search

      if(gOPT.noise!=0) {
        // restore original sparse maps
        int nIFO = NET->ifoListSize(); 	// number of detectors
        for(int n=0;n<nIFO;n++) {
          detector* pD = NET->getifo(n);
          pD->vSS=vSS[n];
        }
      }

      // DUMP CED
      cfg->cedDump=gCEDDUMP;		// restore user_parameter cedDump (executed in the last trial)
    }
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {			// AFTER EVENT RECONSTRUCTION (for each event repeated gOPT.trials)

    // import trials
    int gLRETRY=-1; IMPORT(int,gLRETRY)

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];
    double ofactor=0;
    if(cfg->simulation==4)      ofactor=-gIFACTOR;
    else if(cfg->simulation==3) ofactor=-gIFACTOR;
    else                        ofactor=factor;

    int nIFO = NET->ifoListSize();                      // number of detectors
    int K = NET->nLag;                                  // number of time lag
    int rate = 0;                                       // select all resolutions
    netevent* EVT;
    wavearray<double> id;

    SetOutputFile(NET, EVT, cfg, false);                // set output root file

    for(int k=0; k<K; k++) {  				// loop over the lags

      id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);

      if(id.size()==0) cfg->cedDump=gCEDDUMP;           // rejected -> restore user_parameter cedDump

      for(int j=0; j<(int)id.size(); j++) {  		// loop over cluster index

        int ID = size_t(id.data[j]+0.5);

        if(NET->getwc(k)->sCuts[ID-1]!=-1)  continue; 	// skip rejected/processed clusters

        EVT->output2G(NULL,NET,ID,k,ofactor);		// get reconstructed parameters

        for(int n=0; n<nIFO; n++) {			// save SNR
           iSNR[n].push_back(EVT->iSNR[n]);
           oSNR[n].push_back(EVT->oSNR[n]);
          ioSNR[n].push_back(EVT->ioSNR[n]);
        }
        for(int n=0; n<6; n++) {			// save chirp mass parameters
          vCHIRP[n].push_back(EVT->chirp[n]); 
        }
        vLIKELIHOOD.push_back(EVT->likelihood);		// save likelihood

        // print event parameters
        cout << endl;
        cout << "event parameters : ID -> " << ID << endl;
        for(int n=0;n<nIFO;n++) printf("rec_time %s : %.4f\n",NET->ifoName[n], EVT->time[n]);
        cout << "rec_theta : " << EVT->theta[0] << " rec_phi : " << EVT->phi[0] << endl;
        cout << "SNRnet : " << sqrt(EVT->likelihood) << " netcc[0] : " << EVT->netcc[0] 
             << " rho[0] : " << EVT->rho[0] << " size : " << EVT->size[0] << endl;

        if(gOPT.trials==0) {
          DumpOutputFile(NET, EVT, cfg, ID, k, ofactor); // dump event to output root file
	  continue;
        }

        // -------------------------------------------------
        // FIRST TRIAL
        // reconstruct signal with original noise (full sky)
        // -------------------------------------------------
        if((gLRETRY==gOPT.trials) || (gOPT.multitask && gLRETRY==1)) {          
          // get waveforms
          if(cfg->simulation) {
            vINJ = GetInjWaveform(NET, EVT, ID, factor);  // get injected waveform
            if(vINJ.size()!=nIFO) {
              cout << "CWB_Plugin_PE Error : Injection Waveform Not Found !!!" << endl;
              exit(1);
            }
            cINJ = vINJ;				  // cINJ is overwritten by GetSigWaveform
          }
          vREC = GetRecWaveform(NET, EVT, ID);            // get reconstructed waveform
          vWHT = GetWhitenedData(NET, cfg);            	  // get whitened data (in the vREC time range)
          vDAT = GetWaveform(NET, k, ID, 'W');    	  // get reconstructed+noise waveform
          vNUL = GetWaveform(NET, k, ID, 'N');    	  // get noise-reconstructed waveform
          for(int n=0;n<nIFO;n++) {			  // compute the istantaneous frequency of reconstructed waveform
            fREC.push_back(CWB::Toolbox::getHilbertIFrequency(vREC[n])); 
          }

          Wave2Sparse(NET,cfg,'v');                       // save original sparse map to vSS
          Wave2Sparse(NET,cfg,'r');                       // rec to rSS
          Wave2Sparse(NET,cfg,'d');                       // dat to dSS
          if(cfg->simulation) Wave2Sparse(NET,cfg,'j');   // inj to jSS
          // get injected waveform on TF sub-space of the reconstructed waveform vREC
          if(cfg->simulation) {
            cINJ = GetSigWaveform(NET, cfg, k, ID);
            for(int n=0;n<nIFO;n++) {			  // compute the istantaneous frequency of cutted injected waveform
              fINJ.push_back(CWB::Toolbox::getHilbertIFrequency(cINJ[n])); 
            }
          }

          if(gOPT.skymask==1) SaveSkyProb(NET,cfg,ID);	  // save sky probability, used for trials

          GetNullPixels(aNUL, ANUL, NET, cfg, k, ID);  	  // get null pixels in TF domain
          GetNoisePixels(aNSE, ANSE, NET, cfg, k, ID);    // get noise pixels in TF domain 
          for(int n=0; n<nIFO; n++) {			  // fill gPE.nstat output 	
            GetNullStatistic(&aNUL[n], &aNSE[n], n, NET->ifoName[n], "null_00");	
            GetNullStatistic(&ANUL[n], &ANSE[n], n, NET->ifoName[n], "null_90");
          }

          // save event parameters
          gSEGGPS = EVT->gps[0];			  // save segment start time (sec) 
          gITHETA = EVT->theta[1]; gIPHI   = EVT->phi[1]; // save injected sky position 
          gOTHETA = EVT->theta[3]; gOPHI   = EVT->phi[3]; // save reconstructed sky poosition (the maximum of detection statistic)

          //vPCA = GetPCAWaveform(NET, cfg, k, ID);	  // get pca waveform
        }

        // -----------------------------------------------------------------------------------
        // AFTER FIRST TRIAL
        // reconstruct signal with reconstructed signal + gaussian noise & (calib. errors)
        // -----------------------------------------------------------------------------------
        if(gLRETRY!=gOPT.trials) {                       
          wREC[gLRETRY] = GetRecWaveform(NET, EVT, ID);   	// get reconstructed waveform
          wSKYPROB[gLRETRY] = GetSkyProb(NET, ID);   		// get sky probability
          gID=ID;						// used in multitask mode
        }

        // -------------------------------------------------
        // LAST TRIAL
        // reconstruct signal with original noise (full sky)
        // -------------------------------------------------
        if(gLRETRY==0) { 
          std::vector<int> nrec = ComputeStatisticalErrors(NET, cfg, ID);
	  if(cfg->simulation) GetFactorsStatistic(nIFO);  	// fill gPE.ff, gPE.of
	  GetSNRnetStatistic(nIFO);			  	// fill gPE.snet
          GetChirpMassStatistic(vCHIRP);			// fill gPE.mch
          mSKYPROB = GetMedianSkyProb(NET);			// get median sky probability

          SetOutputFile(NET, EVT, cfg, true);                  	// set output root file
          DumpOutputFile(NET, EVT, cfg, ID, k, ofactor); 	// dump event to output root file
          if(cfg->cedDump) {
            TString dirCED = DumpCED(NET, EVT, cfg, ofactor);	// dump CED 
            cout << "dirCED : " << dirCED << endl;
            DumpSkyProb(&mSKYPROB, NET, EVT, dirCED); 		// dump median mSKYPROB to fits file
            TString ifo[NIFO_MAX];
            for(int n=0; n<nIFO; n++) ifo[n]=NET->ifoName[n];
            CreateIndexHTML(dirCED, nIFO,ifo,cfg->simulation);	// create html file
            if(nrec.size()) {
              PlotWaveforms(NET,cfg,ID,dirCED);			// plot waveforms
              if(gOPT.ced_distr) {
                PlotResiduals(ID,dirCED,cfg->simulation,'w');	// plot residual enegy
                PlotResiduals(ID,dirCED,cfg->simulation,'f');	// plot frequency residual enegy
                gStyle->SetOptFit(1111);
                PlotSNRnet(nIFO,dirCED,cfg->simulation);	// plot SNRnet
              }
              if(gOPT.ced_cm) {
                for(int n=1;n<6;n++) PlotChirpMass(n,dirCED,cfg->simulation); // plot mchirp
              }
              PlotFactors(0,nIFO,dirCED);			// plot fitting factor
              PlotFactors(1,nIFO,dirCED);			// plot overlap factor;
              if(gOPT.ced_null) {
                for(int n=0; n<nIFO; n++) {			// plot null statistic	
                  GetNullStatistic(&aNUL[n], &aNSE[n], n, ifo[n], "null_00", dirCED);	
                  GetNullStatistic(&ANUL[n], &ANSE[n], n, ifo[n], "null_90", dirCED);
                }
              }
              gStyle->SetOptFit(0000);
              DumpRecWavePE(NET,dirCED); 			// dumps rec waveform/time/freq/errors array in ASCII format.
              if(cfg->simulation) DumpInjWavePE(NET,dirCED); 	// dumps inj waveform/time in ASCII format.
            }
            DumpUserOptions(dirCED, cfg);			// dump PE user options

            if(gOPT.ced_pca)    PlotTimeFrequencyPCA(NET, EVT, cfg, ID, k, dirCED);	// plot tf likelihood pca
            if(gOPT.ced_null)   PlotSpectrogram("null", NET, EVT, cfg, dirCED);		// plot null spectrograms
            if(gOPT.ced_null)   PlotSpectrogram("signal", NET, EVT, cfg, dirCED);	// plot reconstructed signal spectrograms
            if(cfg->simulation) PlotSpectrogram("injection", NET, EVT, cfg, dirCED);	// plot injected signal spectrograms
	  }

          for(int n=0;n<nIFO;n++) {
            detector* pD = NET->getifo(n);
            pD->vSS=vSS[n]; 					// restore original sparse maps
            if(!cfg->simulation) ClearWaveforms(pD);     	// release waveform memory
          }
        }
      }
    }

    jfile->cd();
    if(EVT) delete EVT;
  }

  // if gOPT.noise=0 the reconstructed is injected in the whitened data and Coherence+SuperCluster is done
  if(type==CWB_PLUGIN_OLIKELIHOOD) {				
    if(gOPT.noise==0) {
      int gLRETRY=-1; IMPORT(int,gLRETRY)
      int csize=0;
      do {
        cout << endl;
        cout << "-------------------------------------------------------------------" << endl;
        cout << "-----> CWB_Plugin_PE.C -> RedoAnalysis : " 
             << " -> gLRETRY : " << gLRETRY << "/" << gOPT.trials << endl;
        cout << "-------------------------------------------------------------------" << endl;
        csize = RedoAnalysis(jfile, cfg, NET);
        EXPORT(int,gLRETRY,TString::Format("gLRETRY = %d;",--gLRETRY).Data())
      } while(csize==0 && gLRETRY>0);
      EXPORT(int,gLRETRY,TString::Format("gLRETRY = %d;",++gLRETRY).Data())
    }
  }

  return;
}

std::vector<netpixel> DoPCA(network* NET, CWB::config* cfg, int lag, int id) {

  double ee;

  size_t nIFO = NET->ifoList.size();

  float  En = 2*NET->acor*NET->acor*nIFO;     // network energy threshold in the sky loop

  int size = NET->a_00.size();
  int f_ = NIFO/4;

  netpixel* pix;
  netcluster* pwc = NET->getwc(lag);
  std::vector<netpixel> vPIX;

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, id, false);   // buffer for pixel IDs
  int V = pI.size();
  if(V>cfg->BATCH) return vPIX;                                 // attach TD amp to pixels < V

  wavearray<float>  xi(size); xi=0;           // PC 00 array
  wavearray<float>  XI(size); XI=0;           // PC 90 array

  __m128* _xi  = (__m128*) xi.data;           // set pointer to PC 00 array
  __m128* _XI  = (__m128*) XI.data;           // set pointer to PC 90 array

  __m128* _aa  = (__m128*) NET->a_00.data;    // set pointer to 00 array
  __m128* _AA  = (__m128*) NET->a_90.data;    // set pointer to 90 array

  int nPC = 0;
  for(int j=0; j<V; j++) {
    int jf = j*f_;                            // source sse pointer increment
    _sse_zero_ps(_xi+jf);                     // zero MRA amplitudes
    _sse_zero_ps(_XI+jf);                     // zero MRA amplitudes
    ee = _sse_abs_ps(_aa+jf,_AA+jf);          // total pixel energy / quadrature
    if(ee>En) nPC++; else ee=0.;              // count core pixels
    NET->rNRG.data[j] = ee;                   // init residual energy array
    NET->pNRG.data[j] = NET->rNRG.data[j];    // init residual energy array
  }

  nPC = NET->_sse_mra_ps(xi.data,XI.data,En,nPC);  // get principal components

  for(int j=0; j<V; j++) {                    // loop over principal components
    pix = pwc->getPixel(id,pI[j]);
    vPIX.push_back(*pix);                     // save original pixels
    pix->core = false;
    ee = NET->pNRG.data[j];                   // total pixel energy
    if(ee<En) continue;
    pix->core = true;
    for(int i=0; i<nIFO; i++) {
       pix->setdata(double(xi.data[j*NIFO+i]),'S',i);    // store 00 whitened response PC
       pix->setdata(double(XI.data[j*NIFO+i]),'P',i);    // store 90 whitened response PC
    }
  }

  return vPIX;
}

void ResetPCA(network* NET, CWB::config* cfg, netcluster* pwc, std::vector<netpixel>* vPIX, int ID) {

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, ID, false);   // buffer for pixel IDs
  int V = pI.size();
  if(V>cfg->BATCH) return;                                      // attach TD amp to pixels < V
  for(int j=0; j<V; j++) {
    netpixel* pix = pwc->getPixel(ID,pI[j]);
    *pix = (*vPIX)[j];
  }

  while(!vPIX->empty()) vPIX->pop_back();
  vPIX->clear(); std::vector<netpixel>().swap(*vPIX);
}

std::vector<wavearray<double> > GetPCAWaveform(network* NET, CWB::config* cfg, int lag, int id) {

  std::vector<wavearray<double> > xPCA;       // reconstructed

  std::vector<netpixel> vPIX;
  vPIX = DoPCA(NET, cfg, lag, id);            // do PCA analysis

  int nIFO = NET->ifoListSize();              // number of detectors
  netcluster* pwc = NET->getwc(lag);

  if(NET->getMRAwave(id,lag,'S',0,true)) {    // reconstruct whitened shifted pd->waveForm
    detector* pd;
    for(int i=0; i<nIFO; i++) {               // loop over detectors
       pd = NET->getifo(i);
       wavearray<double>* wf = &pd->waveForm;
       wf->start(pwc->start+pd->waveForm.start());
       xPCA.push_back(*wf);
    }
  }

  ResetPCA(NET, cfg, pwc, &vPIX, id);         // restore WP pwc
 
  return xPCA;
}

void PlotTimeFrequencyPCA(network* NET, netevent* EVT, CWB::config* cfg, int id, int lag, TString pdir) {

  std::vector<netpixel> vPIX;
  vPIX = DoPCA(NET, cfg, lag, id);            // do PCA analysis

  int nIFO = NET->ifoListSize();              // number of detectors
  netcluster* pwc = NET->getwc(lag);

  TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",EVT->gps[0]);

  watplot WTS(const_cast<char*>("wts"));
  WTS.canvas->cd();
  char fname[1024];
  sprintf(fname, "%s/pca_l_tfmap_scalogram.png",pdir.Data());
  cout << "write " << fname << endl;
  WTS.plot(pwc, 1, nIFO, 'L', 0, const_cast<char*>("COLZ"));
  WTS.hist2D->GetYaxis()->SetRangeUser(EVT->low[0],EVT->high[0]);
  WTS.hist2D->GetXaxis()->SetTitle(xtitle);
  WTS.canvas->Print(fname);
  WTS.clear();
/*
  sprintf(fname, "%s/pca_n_tfmap_scalogram.png",pdir.Data());
  cout << "write " << fname << endl;
  WTS.plot(pwc, 1, nIFO, 'N', 0, const_cast<char*>("COLZ"));
  WTS.canvas->Print(fname);
  WTS.clear();
*/

  ResetPCA(NET, cfg, pwc, &vPIX, id);         // restore WP pwc

}

std::vector<wavearray<double> > GetInjWaveform(network* NET, netevent* EVT, int ID, double factor) {

  std::vector<wavearray<double> > xINJ;        // injected

  int nIFO = NET->ifoListSize();               // number of detectors

  double recTime = EVT->time[0];               // rec event time det=0

  injection INJ(nIFO);
  // get inj ID
  int M = NET->mdc__IDSize();
  double injTime = 1.e12;
  int injID   = -1;
  for(int m=0; m<M; m++) {
     int mdcID = NET->getmdc__ID(m);
     double T = fabs(recTime - NET->getmdcTime(mdcID));
     if(T<injTime && INJ.fill_in(NET,mdcID)) {
        injTime = T;
        injID = mdcID;
     }
  }

  if(INJ.fill_in(NET,injID)) {                 // get inj parameters

     wavearray<double>** pwfINJ = new wavearray<double>*[nIFO];

     // extract whitened injected & reconstructed waveforms
     for(int n=0; n<nIFO; n++) {

        detector* pd = NET->getifo(n);

        pwfINJ[n] = INJ.pwf[n];
        if (pwfINJ[n]==NULL) {
           cout << "CWB_Plugin_Boostrap.C : Error : Injected waveform not saved !!! : detector "
                << NET->ifoName[n] << endl;
           continue;
        }

        double rFactor = 1.;
        rFactor *= factor>0 ? factor : 1;
        wavearray<double>* wf = pwfINJ[n];
        *wf*=rFactor;                       // injected wf is multiplied by the custom factor
        xINJ.push_back(*wf);
        *wf*=1./rFactor;                    // restore injected amplitude
     }
     delete [] pwfINJ;
  }

  return xINJ;
}

std::vector<wavearray<double> > GetRecWaveform(network* NET, netevent* EVT, int ID) {

  std::vector<wavearray<double> > xREC;        // reconstructed

  int nIFO = NET->ifoListSize();               // number of detectors

  wavearray<double>** pwfREC = new wavearray<double>*[nIFO];
  for(int n=0; n<nIFO; n++) {

     detector* pd = NET->getifo(n);
     int idSize = pd->RWFID.size();
     int wfIndex=-1;
     for(int mm=0; mm<idSize; mm++) if (pd->RWFID[mm]==ID) wfIndex=mm;

     if(wfIndex<0) {
        cout << "CWB_Plugin_Boostrap.C : Error : Reconstructed waveform not saved !!! : ID -> "
             << ID << " : detector " << NET->ifoName[n] << endl;
        continue;
     }
     if(wfIndex>=0) pwfREC[n] = pd->RWFP[wfIndex];

     wavearray<double>* wf = pwfREC[n];
     xREC.push_back(*wf);
  }
  delete [] pwfREC;

  return xREC;
}

std::vector<wavearray<double> > GetWaveform(network* NET, int lag, int id, char type, bool shift) {

  if(type!='S' && type!='s' && type!='W' && type!='w' && type!='N' && type!='n') {
    cout << "GetWaveform : not valid type : Abort" << endl; exit(1);
  } 

  std::vector<wavearray<double> > vWAV;       // reconstructed

  int nIFO = NET->ifoListSize();              // number of detectors
  netcluster* pwc = NET->getwc(lag);

  if(type=='S' || type=='s') {
    NET->getMRAwave(id,lag,'S',0,shift);       // reconstruct whitened shifted pd->waveForm
    detector* pd;
    for(int i=0; i<nIFO; i++) {               // loop over detectors
       pd = NET->getifo(i);
       wavearray<double>* wf = &pd->waveForm;
       wf->start(pwc->start+pd->waveForm.start());
       vWAV.push_back(*wf);
    }
  }

  if(type=='W' || type=='w') {
    NET->getMRAwave(id,lag,type,0,shift);       // reconstruct+noise whitened shifted pd->waveBand
    detector* pd;
    for(int i=0; i<nIFO; i++) {               // loop over detectors
       pd = NET->getifo(i);
       wavearray<double>* wf = &pd->waveBand;
       wf->start(pwc->start+pd->waveBand.start());
       vWAV.push_back(*wf);
    }
  }

  if(type=='N' || type=='n') {
    NET->getMRAwave(id,lag,'W',0,shift);       // reconstruct+noise whitened shifted pd->waveBand
    NET->getMRAwave(id,lag,'S',0,shift);       // reconstruct whitened shifted pd->waveForm
    detector* pd;
    for(int i=0; i<nIFO; i++) {               // loop over detectors
       pd = NET->getifo(i);
       pd->waveNull = pd->waveBand;
       pd->waveNull-= pd->waveForm;
       wavearray<double>* wf = &pd->waveNull;
       wf->start(pwc->start+pd->waveNull.start());
       vWAV.push_back(*wf);
    }
  }

  return vWAV;
}

std::vector<wavearray<double> > GetSigWaveform(network* NET, CWB::config* cfg, int lag, int id) {
// get injected waveform on TF sub-space of the reconstructed waveform vREC

  std::vector<wavearray<double> > vSIG;       	// reconstructed

  int nIFO = NET->ifoListSize();              	// number of detectors
  netcluster* pwc = NET->getwc(lag);

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, id, false); // buffer for pixel IDs

  netpixel* pix = pwc->getPixel(id,pI[0]);

  int tsize=1;

  int V = pI.size();
  if(!V) return vSIG;
  int V4  = V + (V%4 ? 4 - V%4 : 0);
  int V44 = V4 + 4;

  float En = 2*NET->acor*NET->acor*nIFO;        // network energy threshold in the sky loop

  float  *pd[NIFO], *pD[NIFO];
  float  *v00[NIFO],*v90[NIFO];
  float  *pa[NIFO], *pA[NIFO];

  float* ptmp;
  int i,j;

  std::vector<float*> _DAT;              	// vectors for packet amplitudes
  std::vector<float*> _vtd;              	// vectors of TD amplitudes
  std::vector<float*> _vTD;              	// vectors of TD amplitudes
 
  if(_DAT.size()) _avx_free_ps(_DAT);           // container for data packet amplitudes
  if(_vtd.size()) _avx_free_ps(_vtd);           // array for 00 amplitudes
  if(_vTD.size()) _avx_free_ps(_vTD);           // array for 90 amplitudes
  for(i=0; i<NIFO; i++) {
    ptmp = (float*)_mm_malloc((V4*3+8)*sizeof(float),32);
    for(j=0; j<(V4*3+8); j++) ptmp[j]=0; _DAT.push_back(ptmp);   // concatenated arrays {amp}{AMP}{norm}{n,N,c,s}
    ptmp = (float*)_mm_malloc(tsize*V4*sizeof(float),32);
    for(j=0; j<tsize*V4; j++) ptmp[j]=0; _vtd.push_back(ptmp);   // array of aligned vectors
    ptmp = (float*)_mm_malloc(tsize*V4*sizeof(float),32);
    for(j=0; j<tsize*V4; j++) ptmp[j]=0; _vTD.push_back(ptmp);   // array of aligned vectors
  }

  std::vector<float*> _AVX;              	// vectors for network pixel statistics
  if(_AVX.size()) _avx_free_ps(_AVX);
  float* p_et = (float*)_mm_malloc(V4*sizeof(float),32);      // 0
  for(j=0; j<V4; j++) p_et[j]=0; _AVX.push_back(p_et);
  float* pMSK = (float*)_mm_malloc(V44*sizeof(float),32);     // 1  - pixel mask
  for(j=0; j<V44; j++) pMSK[j]=0; _AVX.push_back(pMSK);   pMSK[V4]=nIFO;
  float* p_fp = (float*)_mm_malloc(V44*sizeof(float),32);     // 2- |f+|^2 (0:V4), +norm (V4:V4+4)
  for(j=0; j<V44; j++) p_fp[j]=0; _AVX.push_back(p_fp);
  float* p_fx = (float*)_mm_malloc(V44*sizeof(float),32);     // 3- |fx|^2 (0:V4), xnorm (V4:V4+4)
  for(j=0; j<V44; j++) p_fx[j]=0; _AVX.push_back(p_fx);
  float* p_si = (float*)_mm_malloc(V4*sizeof(float),32);      // 4
  for(j=0; j<V4; j++) p_si[j]=0; _AVX.push_back(p_si);
  float* p_co = (float*)_mm_malloc(V4*sizeof(float),32);      // 5
  for(j=0; j<V4; j++) p_co[j]=0; _AVX.push_back(p_co);
  float* p_uu = (float*)_mm_malloc((V4+16)*sizeof(float),32);  // 6 - 00+ unit vector(0:V4), norm(V4), cos(V4+4)
  for(j=0; j<V4+16; j++) p_uu[j]=0; _AVX.push_back(p_uu);
  float* p_UU = (float*)_mm_malloc((V4+16)*sizeof(float),32);  // 7 - 90+ unit vector(0:V4), norm(V4), sin(V4+4)
  for(j=0; j<V4+16; j++) p_UU[j]=0; _AVX.push_back(p_UU);
  float* p_vv = (float*)_mm_malloc((V4+16)*sizeof(float),32);  // 8- 00x unit vector(0:V4), norm(V4), cos(V4+4)
  for(j=0; j<V4+16; j++) p_vv[j]=0; _AVX.push_back(p_vv);
  float* p_VV = (float*)_mm_malloc((V4+16)*sizeof(float),32);  // 9- 90x unit vector(0:V4), norm(V4), sin(V4+4)
  for(j=0; j<V4+16; j++) p_VV[j]=0; _AVX.push_back(p_VV);
  float* p_au = (float*)_mm_malloc(V4*sizeof(float),32);      // 10
  for(j=0; j<V4; j++) p_au[j]=0; _AVX.push_back(p_au);
  float* p_AU = (float*)_mm_malloc(V4*sizeof(float),32);      // 11
  for(j=0; j<V4; j++) p_AU[j]=0; _AVX.push_back(p_AU);
  float* p_av = (float*)_mm_malloc(V4*sizeof(float),32);      // 12
  for(j=0; j<V4; j++) p_av[j]=0; _AVX.push_back(p_av);
  float* p_AV = (float*)_mm_malloc(V4*sizeof(float),32);      // 13
  for(j=0; j<V4; j++) p_AV[j]=0; _AVX.push_back(p_AV);
  float* p_uv = (float*)_mm_malloc(V4*4*sizeof(float),32);    // 14 special array for GW norm calculation
  for(j=0; j<V4*4; j++) p_uv[j]=0; _AVX.push_back(p_uv);
  float* p_ee = (float*)_mm_malloc(V4*sizeof(float),32);      // 15 + energy array
  for(j=0; j<V4; j++) p_ee[j]=0; _AVX.push_back(p_ee);
  float* p_EE = (float*)_mm_malloc(V4*sizeof(float),32);      // 16 x energy array
  for(j=0; j<V4; j++) p_EE[j]=0; _AVX.push_back(p_EE);
  float* pTMP=(float*)_mm_malloc(V4*4*NIFO*sizeof(float),32); // 17 temporary array for _avx_norm_ps()
  for(j=0; j<V4*4*NIFO; j++) pTMP[j]=0; _AVX.push_back(pTMP);
  float* p_ni = (float*)_mm_malloc(V4*sizeof(float),32);      // 18 + network index 
  for(j=0; j<V4; j++) p_ni[j]=0; _AVX.push_back(p_ni);
  float* p_ec = (float*)_mm_malloc(V4*sizeof(float),32);      // 19 + coherent energy 
  for(j=0; j<V4; j++) p_ec[j]=0; _AVX.push_back(p_ec);
  float* p_gn = (float*)_mm_malloc(V4*sizeof(float),32);      // 20 + Gaussian noise correction 
  for(j=0; j<V4; j++) p_gn[j]=0; _AVX.push_back(p_gn);
  float* p_ed = (float*)_mm_malloc(V4*sizeof(float),32);      // 21 + energy disbalance
  for(j=0; j<V4; j++) p_ed[j]=0; _AVX.push_back(p_ed);
  float* p_rn = (float*)_mm_malloc(V4*sizeof(float),32);      // 22 + sattelite noise in TF domain
  for(j=0; j<V4; j++) p_rn[j]=0; _AVX.push_back(p_rn);

  double R = NET->getifo(0)->TFmap.rate();
  for(int j=0; j<V; j++) {
    netpixel* pix = pwc->getPixel(id,pI[j]);
    if(!pix->core) continue;             // select only core pixels
    for(int n=0; n<nIFO; n++) {
      int index = pix->data[n].index;
      int ires = int(TMath::Log2(R/pix->rate))-cfg->l_low;
      double aa = GetSparseMapData(&jSS[n][ires], true, index);
      double AA = GetSparseMapData(&jSS[n][ires], false , index);
      //aa = GetSparseMapData(&rSS[n][ires], true, index);
      //AA = GetSparseMapData(&rSS[n][ires], false , index);
      _vtd[n][j] = aa;                     // copy 00 data
      _vTD[n][j] = AA;                     // copy 90 data
    }
  }

  for(int i=0; i<NIFO; i++) {            // set up zero delay and packet pointers
    pd[i] = _DAT[i]; 
    pD[i] = _DAT[i]+V4;
    pa[i] = _vtd[i];
    pA[i] = _vTD[i];
    v00[i] = pa[i];
    v90[i] = pA[i];
  }

  En=0.0001;
  int M;
  double Eo=0;
  wavearray<float> D_snr;
  Eo    = _avx_loadata_ps(v00,v90,pd,pD,En,_AVX,V4);  // calculate data stats and store in _AVX
  Eo    = _avx_packet_ps(pd,pD,_AVX,V4);              // get data packet
  D_snr = NET->_avx_norm_ps(pd,pD,_AVX,V4);           // data packet energy snr
  M     = _avx_setAMP_ps(pd,pD,_AVX,V4);              // set data packet amplitudes for data

  for(int j=0; j<V; j++) {                           // loop over principle components
     pix = pwc->getPixel(id,pI[j]);
     pix->likelihood = (pMSK[j]>0) ? (p_ee[j]+p_EE[j])/2 : 0; // total pixel energy
     if(pMSK[j]>0) pix->core = true;
     //pix->core = true;
     for(int i=0; i<nIFO; i++) {
        pix->setdata(double(pd[i][j]),'S',i);        // 00 whitened
        pix->setdata(double(pD[i][j]),'P',i);        // 90 whitened
     }
  }
  vSIG = GetWaveform(NET, lag, id, 'S', false);      // get reconstructed+noise waveform

#ifdef SAVE_REDUCED_INJ	
  for(int i=0; i<nIFO; i++) {
    char title[256];
    char ofname[256];
    sprintf(title,"%s (TIME) : vREC(red) - vSIG(blue)",NET->ifoName[i]);
    sprintf(ofname,"%s_vREC_vSIG_time_id_%d.%s",NET->ifoName[i],id,"root");
    PlotWaveform(ofname, title, cfg,&vREC[i], &vSIG[i], NULL, NULL, false);   // time
    sprintf(title,"%s (TIME) : vINJ(red) - vSIG(blue)",NET->ifoName[i]);
    sprintf(ofname,"%s_vINJ_vSIG_time_id_%d.%s",NET->ifoName[i],id,"root");
    PlotWaveform(ofname, title, cfg,&vINJ[i], &vSIG[i], NULL, NULL, false);   // time
  }
#endif

  if(_AVX.size()) _avx_free_ps(_AVX);
  if(_DAT.size()) _avx_free_ps(_DAT);           // container for data packet amplitudes
  if(_vtd.size()) _avx_free_ps(_vtd);           // array for 00 amplitudes
  if(_vTD.size()) _avx_free_ps(_vTD);           // array for 90 amplitudes

  return vSIG;
}
  
void GetFactorsStatistic(int nIFO) {

  gPE.ff[0]=0; gPE.ff[1]=0;
  gPE.of[0]=0; gPE.of[1]=0;

  int size = iSNR[0].size();
  if(size==0) return;

  for(int i=0;i<size;i++)  {
    double  isnr=0; for(int n=0;n<nIFO;n++)  isnr+= iSNR[n][i]; 
    double  osnr=0; for(int n=0;n<nIFO;n++)  osnr+= oSNR[n][i]; 
    double iosnr=0; for(int n=0;n<nIFO;n++) iosnr+=ioSNR[n][i]; 
    double ff = iosnr/sqrt(isnr*isnr);		// fitting factor
    double of = iosnr/sqrt(isnr*osnr);		// overlap factor
    
    gPE.ff[0] += ff;
    gPE.ff[1] += pow(ff,2);
    
    gPE.of[0] += of;
    gPE.of[1] += pow(of,2);
  }

  gPE.ff[0] /= size;
  gPE.ff[1] = sqrt(gPE.ff[1]/size-gPE.ff[0]*gPE.ff[0]);

  gPE.of[0] /= size;
  gPE.of[1] = sqrt(gPE.of[1]/size-gPE.of[0]*gPE.of[0]);
}

void GetChirpMassStatistic(std::vector<double>*  vCHIRP) {

  gPE.mch[0]=0;
  gPE.mch[1]=0;

  int size = vCHIRP[1].size();
  if(size==0) return;

  for(int i=0;i<size;i++)  {
    gPE.mch[0] += vCHIRP[1][i];
    gPE.mch[1] += pow(vCHIRP[1][i],2);
  }

  gPE.mch[0] /= size;
  gPE.mch[1]  = sqrt(gPE.mch[1]/size-gPE.mch[0]*gPE.mch[0]);
}

void GetSNRnetStatistic(int nIFO) {

  gPE.snet[0]=0;
  gPE.snet[1]=0;

  int size = iSNR[0].size();
  if(size==0) return;

  for(int i=0;i<size;i++)  {
    double osnr=0;
    for(int n=0;n<nIFO;n++) osnr+=oSNR[n][i];
    gPE.snet[0] += sqrt(osnr);
    gPE.snet[1] += osnr;
  }

  gPE.snet[0] /= size;
  gPE.snet[1]  = sqrt(gPE.snet[1]/size-gPE.snet[0]*gPE.snet[0]);
}

void GetNullPixels(std::vector<double>* aNUL, std::vector<double>* ANUL, 
                   network* NET, CWB::config* cfg, int lag, int id) {
// get null pixels in TF domain
// we use delayed vSS amplitudes and the reconstructed TF pixels save in the NET->a_00,NET->a_90 arrays

  size_t nIFO = NET->ifoList.size();

  netcluster* pwc = NET->getwc(lag);

  // the TD amplitudes are cleaned by the likelihoodWP function, must be computed again
  int sCuts = pwc->sCuts[id-1];			// save cluster status
  pwc->sCuts[id-1] = 0;				// cluster must be enabled (mark as done by likelihoodWP)	
  pwc->setcore(false,id);
  pwc->loadTDampSSE(*NET, 'a', cfg->BATCH, cfg->BATCH); // attach TD amp to pixels

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, id, true);   // buffer for pixel IDs
  int V = pI.size();

  netpixel* pix = pwc->getPixel(id,pI[0]);
  int tsize = pix->tdAmp[0].size();
  if(!tsize || tsize&1) {                       // tsize%1 = 1/0 = power/amplitude
     cout<<"GetNullPixels error: wrong pixel TD data\n";
     exit(1);
  }
  tsize /= 2;
  cout << "tsize :" << tsize << endl;
  int V4  = V + (V%4 ? 4 - V%4 : 0);
  std::vector<int>* vtof = &(pwc->nTofF[id-1]);

  for(int n=0; n<2*7*nIFO; n++) gPE.nstat[n]=0;

  int cnt[NIFO_MAX];
  for(int m=0; m<nIFO; m++) cnt[m]=0;

  float* a_00 = NET->a_00.data;    // set pointer to 00 array of reconstructed signal
  float* a_90 = NET->a_90.data;    // set pointer to 90 array of reconstructed signal

  double R = NET->getifo(0)->TFmap.rate();
  for(int j=0; j<V; j++) {
    netpixel* pix = pwc->getPixel(id,pI[j]);
    //if(!pix->core) continue;       // select only core pixels
    float ee=0;
    for(int m=0; m<nIFO; m++) {
      double ss = a_00[j*NIFO+m];
      double SS = a_90[j*NIFO+m];
      ee += pow(ss,2)+pow(SS,2);
    }
    if(ee==0) continue;		     // only core pixels are selected
    for(int m=0; m<nIFO; m++) {
      int index = pix->data[m].index;
      int ires = int(TMath::Log2(R/pix->rate))-cfg->l_low;
      // get original sparse map amplitudes
      //double dd = GetSparseMapData(&vSS[m][ires], true, index);
      //double DD = GetSparseMapData(&vSS[m][ires], false , index);
      //double dd = GetSparseMapData(&dSS[m][ires], true, index);
      //double DD = GetSparseMapData(&dSS[m][ires], false , index);

      // get shifted data (signal+noise) amplitudes
      int l = (*vtof)[m]+tsize/2; 
      double dd = pix->tdAmp[m].data[l];             // copy TD 00 data
      double DD = pix->tdAmp[m].data[l+tsize];       // copy TD 90 data

      //double ss = GetSparseMapData(&rSS[m][ires], true, index);
      //double SS = GetSparseMapData(&rSS[m][ires], false , index);
      // get reconstructed amplitudes (shifted)
      double ss = a_00[j*NIFO+m];
      double SS = a_90[j*NIFO+m];

      // 00 phase	
      aNUL[m].push_back(dd-ss); 

      // 90 phase	
      ANUL[m].push_back(DD-SS);

      cnt[m]++;
    }
  }
 
  pwc->sCuts[id-1] = sCuts; 		// restore cluster status 
  pwc->clean(id);			// cean TD ammplitudes

  return;
}

double GetSparseMapData(SSeries<double>* SS, bool phase, int index) {

  int layer = SS->GetLayer(index);
  int start = SS->sparseLookup[layer];                // sparse table layer offset
  int end   = SS->sparseLookup[layer+1]-1;            // sparse table layer+1 offset

  int i = SS->binarySearch(SS->sparseIndex.data,start,end,index);

  if(i<0) return 0;

  return phase ? SS->sparseMap00[i] : SS->sparseMap90[i];
}

void PlotWaveform(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wf1, 
                  wavearray<double>* wf2, wavearray<double>* wf3, wavearray<double>* wref, 
                  bool fft, TString pdir, double P, EColor col1, EColor col2, EColor col3) {

  watplot PTS(const_cast<char*>("ptspe"),200,20,800,500);

  //cout << "Print " << ofname << endl;
  double tmin=1.e20;
  if(wref==NULL) return;
  if(wf1==NULL) return;
  else          tmin=wf1->start();
  if(wf2!=NULL) if(wf2->start()<tmin) tmin=wf2->start();
  if(wf3!=NULL) if(wf3->start()<tmin) tmin=wf3->start();

  tmin=gSEGGPS;

                                           wf1->start(wf1->start()-tmin);
  if(wf2!=NULL) if(wf2!=wf1)               wf2->start(wf2->start()-tmin);
  if(wf3!=NULL) if(wf3!=wf1 && wf3!=wf2)   wf3->start(wf3->start()-tmin);
  if(wref!=wf1 && wref!=wf2  && wref!=wf3) wref->start(wref->start()-tmin);

  double bT, eT;
  GetTimeBoundaries(*wref, P, bT, eT);

  if(fft) {
    wavearray<double> wf=*wf1; int k=0;
    for(int i=0;i<wf1->size();i++) {
      double time=i/wf1->rate()+wf1->start();
      if(time>bT && time<eT) wf[k++]=wf1->data[i];
    }
    wf.resize(k);
    PTS.plot(&wf, const_cast<char*>("AL"), col1, 0, 0, true, cfg->fLow, cfg->fHigh);
  } else {
    PTS.plot(wf1, const_cast<char*>("AL"), col1, bT, eT);
    PTS.graph[0]->GetXaxis()->SetRangeUser(bT, eT);
  }
  PTS.graph[0]->SetLineWidth(1);
  PTS.graph[0]->SetTitle(title);

  TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",tmin);
  PTS.graph[0]->GetXaxis()->SetTitle(xtitle);

  if(wf2!=NULL) {
    if(fft) {
      if(wf2!=NULL) {
        wavearray<double> wf=*wf2; int k=0;
        for(int i=0;i<wf2->size();i++) {
          double time=i/wf2->rate()+wf2->start();
          if(time>bT && time<eT) wf[k++]=wf2->data[i];
        }
        wf.resize(k);
        PTS.plot(&wf, const_cast<char*>("SAME"), col2, 0, 0, true, cfg->fLow, cfg->fHigh);
      }
    } else {
      if(wf2!=NULL) PTS.plot(wf2, const_cast<char*>("SAME"), col2, 0, 0);
    }
    if(wf2!=NULL) PTS.graph[1]->SetLineWidth(2);
  }

  if(wf3!=NULL) {
    if(fft) {
      if(wf3!=NULL) {
        wavearray<double> wf=*wf3; int k=0;
        for(int i=0;i<wf3->size();i++) {
          double time=i/wf3->rate()+wf3->start();
          if(time>bT && time<eT) wf[k++]=wf3->data[i];
        }
        wf.resize(k);
        PTS.plot(&wf, const_cast<char*>("SAME"), col3, 0, 0, true, cfg->fLow, cfg->fHigh);
      }
    } else {
      if(wf3!=NULL) PTS.plot(wf3, const_cast<char*>("SAME"), col3, 0, 0);
    }
    if(wf3!=NULL) PTS.graph[2]->SetLineWidth(1);
  }

                                           wf1->start(wf1->start()+tmin);
  if(wf2!=NULL) if(wf2!=wf1)               wf2->start(wf2->start()+tmin);
  if(wf3!=NULL) if(wf3!=wf1 && wf3!=wf2)   wf3->start(wf3->start()+tmin);
  if(wref!=wf1 && wref!=wf2  && wref!=wf3) wref->start(wref->start()+tmin);

  if(fft) {
    PTS.canvas->SetLogx();
    PTS.canvas->SetLogy();
  }

  if(ofname!="") {
    ofname = TString(pdir)+TString("/")+ofname;
    PTS.canvas->Print(ofname);
    cout << "write : " << ofname << endl;
  }
}

double FittingFactor(wavearray<double>* wf1, wavearray<double>* wf2) {

   double R=wf1->rate();

   double b_wf1 = wf1->start();
   double e_wf1 = wf1->start()+wf1->size()/R;
   double b_wf2 = wf2->start();
   double e_wf2 = wf2->start()+wf2->size()/R;

   int o_wf1 = b_wf1>b_wf2 ? 0 : int((b_wf2-b_wf1)*R+0.5);
   int o_wf2 = b_wf1<b_wf2 ? 0 : int((b_wf1-b_wf2)*R+0.5);

   double startXCOR = b_wf1>b_wf2 ? b_wf1 : b_wf2;
   double endXCOR   = e_wf1<e_wf2 ? e_wf1 : e_wf2;
   int sizeXCOR  = int((endXCOR-startXCOR)*R+0.5);

   double wf11=0.;
   double wf22=0.;
   double wf12=0.;
   for(int i=0;i<sizeXCOR;i++) wf12 += wf1->data[i+o_wf1] * wf2->data[i+o_wf2];
   for(int i=0;i<wf1->size();i++) wf11 += pow(wf1->data[i],2);
   for(int i=0;i<wf2->size();i++) wf22 += pow(wf2->data[i],2);

   double FF = wf12/sqrt(wf11*wf22);

   return FF;
}

wavearray<double> GetDifWaveform(wavearray<double>* wf1, wavearray<double>* wf2) {

   double R=wf1->rate();

   double b_wf1 = wf1->start();
   double e_wf1 = wf1->start()+wf1->size()/R;
   double b_wf2 = wf2->start();
   double e_wf2 = wf2->start()+wf2->size()/R;

   int o_wf1 = b_wf1>b_wf2 ? 0 : int((b_wf2-b_wf1)*R+0.5);
   int o_wf2 = b_wf1<b_wf2 ? 0 : int((b_wf1-b_wf2)*R+0.5);

   double startXCOR = b_wf1>b_wf2 ? b_wf1 : b_wf2;
   double endXCOR   = e_wf1<e_wf2 ? e_wf1 : e_wf2;
   int sizeXCOR  = int((endXCOR-startXCOR)*R+0.5);

   wavearray<double> wfdif(sizeXCOR);
   wfdif=0.;
   wfdif.rate(R);
   wfdif.start(b_wf1+double(o_wf1)/R);

   for(int i=0;i<sizeXCOR;i++) wfdif[i] = wf1->data[i+o_wf1] - wf2->data[i+o_wf2];

   return wfdif;
}

void PlotSparse(int ifoID, network* NET, CWB::config* cfg, int ID, wavearray<double>* wf) {

   detector* pD = NET->getifo(ifoID);

   int nRES = cfg->l_high-cfg->l_low+1;     // number of frequency resolution levels
   for(int i=0;i<nRES;i++) {
      // rebuild wseries from sparse table only with core pixels
      bool core = true;                                         
      SSeries<double> ss = pD->vSS[i];                         
      ss.Expand(core);                                        
      ss.Inverse();                                           
                                                              
      char title[256];
      char ofname[256];
      sprintf(title,"%s (TIME) : INJ(red) - SM(blue)",NET->ifoName[ifoID]);
      sprintf(ofname,"%s_INJ_SM_time_id_%d_res_%d.%s",NET->ifoName[ifoID],ID,i,PLOT_TYPE);
      PlotWaveform(ofname, title, cfg, wf, &ss, NULL, NULL, false);   // time
   }                                                                                   
}                                                                              

void ReadUserOptions(TString options) {

  int n_amp_cal_err=0;
  int n_phs_cal_err=0;
  if(options.CompareTo("")!=0) {                               
    cout << options << endl;                                       
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet
                                                                               
      TObjArray* token = TString(options).Tokenize(TString(' '));         
        for(int j=0;j<token->GetEntries();j++){                                

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();            

        if(stok.Contains("pe_trials=") || stok.Contains("pe_retry=")) {
          TString pe_trials=stok;        
          pe_trials.Remove(0,pe_trials.Last('=')+1);
          if(pe_trials.IsDigit()) gOPT.trials=pe_trials.Atoi();
          if(gOPT.trials>MAX_TRIALS) {
            cout << "ReadUserOptions Error : trials must be <= MAX_TRIALS : " << MAX_TRIALS << endl;exit(1);
          }
        }                                                    

        if(stok.Contains("pe_ced_dump=")) {
          TString pe_ced_dump=stok;        
          pe_ced_dump.Remove(0,pe_ced_dump.Last('=')+1);
          if(pe_ced_dump.IsDigit()) gOPT.ced_dump=pe_ced_dump.Atoi();
        }                                                                

        if(stok.Contains("pe_ced_enable=")) {
          TString pe_ced_enable=stok;        
          pe_ced_enable.Remove(0,pe_ced_enable.Last('=')+1);
          if(pe_ced_enable=="tfmap")  gOPT.ced_tfmap  = true;   
          if(pe_ced_enable=="rdr")    gOPT.ced_rdr    = true;   
          if(pe_ced_enable=="psd")    gOPT.ced_psd    = true;   
          if(pe_ced_enable=="skymap") gOPT.ced_skymap = true;   
          if(pe_ced_enable=="rec")    gOPT.ced_rec    = true;   
          if(pe_ced_enable=="inj")    gOPT.ced_inj    = true;   
          if(pe_ced_enable=="rinj")   gOPT.ced_rinj   = true;   
          if(pe_ced_enable=="cm")     gOPT.ced_cm     = true;   
          if(pe_ced_enable=="distr")  gOPT.ced_distr  = true;   
          if(pe_ced_enable=="null")   gOPT.ced_null   = true;   
          if(pe_ced_enable=="pca")    gOPT.ced_pca    = true;   
        }                                                                

        if(stok.Contains("pe_ced_disable=")) {
          TString pe_ced_disable=stok;        
          pe_ced_disable.Remove(0,pe_ced_disable.Last('=')+1);
          if(pe_ced_disable=="tfmap")  gOPT.ced_tfmap  = false;   
          if(pe_ced_disable=="rdr")    gOPT.ced_rdr    = false;   
          if(pe_ced_disable=="psd")    gOPT.ced_psd    = false;   
          if(pe_ced_disable=="skymap") gOPT.ced_skymap = false;   
          if(pe_ced_disable=="rec")    gOPT.ced_rec    = false;   
          if(pe_ced_disable=="inj")    gOPT.ced_inj    = false;   
          if(pe_ced_disable=="rinj")   gOPT.ced_rinj   = false;   
          if(pe_ced_disable=="cm")     gOPT.ced_cm     = false;   
          if(pe_ced_disable=="distr")  gOPT.ced_distr  = false;   
          if(pe_ced_disable=="null")   gOPT.ced_null   = false;   
          if(pe_ced_disable=="pca")    gOPT.ced_pca    = false;   
        }                                                                

        if(stok.Contains("pe_seed=")) {
          TString pe_seed=stok;        
          pe_seed.Remove(0,pe_seed.Last('=')+1);
          if(pe_seed.IsDigit()) gOPT.seed=pe_seed.Atoi();
        }                                                                

        if(stok.Contains("pe_gps=")) {
          TString pe_gps=stok;        
          pe_gps.Remove(0,pe_gps.Last('=')+1);
          if(pe_gps.IsFloat()) gOPT.gps=pe_gps.Atof();
        }                                                                

        if(stok.Contains("pe_amp_cal_err=")) {
          TString pe_amp_cal_err=stok;        
          pe_amp_cal_err.Remove(0,pe_amp_cal_err.Last('=')+1);
          if(pe_amp_cal_err.IsFloat()) gOPT.amp_cal_err[n_amp_cal_err]=pe_amp_cal_err.Atof();
          if(n_amp_cal_err<(NIFO_MAX-1)) n_amp_cal_err++;
        }                                                                

        if(stok.Contains("pe_phs_cal_err=")) {
          TString pe_phs_cal_err=stok;        
          pe_phs_cal_err.Remove(0,pe_phs_cal_err.Last('=')+1);
          if(pe_phs_cal_err.IsFloat()) gOPT.phs_cal_err[n_phs_cal_err]=pe_phs_cal_err.Atof();
          if(n_phs_cal_err<(NIFO_MAX-1)) n_phs_cal_err++;
        }                                                                

        if(stok.Contains("pe_skymask=")) {
          TString pe_skymask=stok;        
          pe_skymask.Remove(0,pe_skymask.Last('=')+1);
          if(pe_skymask=="false")  gOPT.skymask=0;   
          if(pe_skymask=="true")   gOPT.skymask=1;    
          if(pe_skymask.IsDigit()) gOPT.skymask=pe_skymask.Atoi();
        }                                               

        if(stok.Contains("pe_signal=")) {
          TString pe_signal=stok;
          pe_signal.Remove(0,pe_signal.Last('=')+1);
          if(pe_signal.IsDigit()) gOPT.signal=pe_signal.Atoi();
        }

        if(stok.Contains("pe_noise=")) {
          TString pe_noise=stok;
          pe_noise.Remove(0,pe_noise.Last('=')+1);
          if(pe_noise=="false")  gOPT.noise=0;   
          if(pe_noise=="true")   gOPT.noise=1;    
          if(pe_noise.IsDigit()) gOPT.noise=pe_noise.Atoi();
        }

        if(stok.Contains("pe_multitask=")) {
          TString pe_multitask=stok;
          pe_multitask.Remove(0,pe_multitask.Last('=')+1);
          if(pe_multitask=="true")  gOPT.multitask=true;   
          if(pe_multitask=="false") gOPT.multitask=false;   
        }

        if(stok.Contains("pe_output_enable=")) {
          TString pe_output_enable=stok;
          pe_output_enable.Remove(0,pe_output_enable.Last('=')+1);
          if(pe_output_enable=="inj")  gOPT.output_inj=true;   
          if(pe_output_enable=="rec")  gOPT.output_rec=true;   
          if(pe_output_enable=="wht")  gOPT.output_wht=true;   
          if(pe_output_enable=="dat")  gOPT.output_dat=true;   
          if(pe_output_enable=="med")  gOPT.output_med=true;   
          if(pe_output_enable=="p50")  gOPT.output_p50=true;   
          if(pe_output_enable=="p90")  gOPT.output_p90=true;   
          if(pe_output_enable=="avr")  gOPT.output_avr=true;   
          if(pe_output_enable=="rms")  gOPT.output_rms=true;   
        }

        if(stok.Contains("pe_output_disable=")) {
          TString pe_output_disable=stok;
          pe_output_disable.Remove(0,pe_output_disable.Last('=')+1);
          if(pe_output_disable=="inj")  gOPT.output_inj=false;   
          if(pe_output_disable=="rec")  gOPT.output_rec=false;   
          if(pe_output_disable=="wht")  gOPT.output_wht=false;   
          if(pe_output_disable=="dat")  gOPT.output_dat=false;   
          if(pe_output_disable=="med")  gOPT.output_med=false;   
          if(pe_output_disable=="p50")  gOPT.output_p50=false;   
          if(pe_output_disable=="p90")  gOPT.output_p90=false;   
          if(pe_output_disable=="avr")  gOPT.output_avr=false;   
          if(pe_output_disable=="rms")  gOPT.output_rms=false;   
        }

      }
    }
  }
}

void PrintUserOptions(CWB::config* cfg) {

    cout << "-----------------------------------------"     << endl;
    cout << "PE config options                       "      << endl;
    cout << "-----------------------------------------"     << endl << endl;
    cout << "PE_TRIALS            " << gOPT.trials          << endl;
    cout << "PE_SIGNAL            " << gOPT.signal          << endl;
    cout << "PE_NOISE             " << gOPT.noise           << endl;
    for(int n=0;n<cfg->nIFO;n++) {
      cout << "PE_AMP_CAL_ERR  " << cfg->ifo[n] << "   " << gOPT.amp_cal_err[n] << endl;
      cout << "PE_PHS_CAL_ERR  " << cfg->ifo[n] << "   " << gOPT.phs_cal_err[n] << endl;
    }
    cout << "PE_SKYMASK           " << gOPT.skymask         << endl;
    cout << "PE_SEED              " << gOPT.seed            << endl;
    cout.precision(14);
    cout << "PE_GPS               " << gOPT.gps             << endl;
    cout.precision(6);

    cout << "PE_MULTITASK         " << gOPT.multitask       << endl;

    cout << "PE_CED_DUMP          " << gOPT.ced_dump        << endl;

    cout << "PE_CED_TFMAP         " << gOPT.ced_tfmap       << endl;
    cout << "PE_CED_RDR           " << gOPT.ced_rdr         << endl;
    cout << "PE_CED_PSD           " << gOPT.ced_psd         << endl;
    cout << "PE_CED_SKYMAP        " << gOPT.ced_skymap      << endl;
    cout << "PE_CED_REC           " << gOPT.ced_rec         << endl;
    cout << "PE_CED_INJ           " << gOPT.ced_inj         << endl;
    cout << "PE_CED_rINJ          " << gOPT.ced_rinj        << endl;
    cout << "PE_CED_CM            " << gOPT.ced_cm          << endl;
    cout << "PE_CED_DISTR         " << gOPT.ced_distr       << endl;
    cout << "PE_CED_NULL          " << gOPT.ced_null        << endl;
    cout << "PE_CED_PCA           " << gOPT.ced_pca         << endl;

    cout << "PE_OUTPUT_INJ        " << gOPT.output_inj      << endl;
    cout << "PE_OUTPUT_REC        " << gOPT.output_rec      << endl;
    cout << "PE_OUTPUT_WHT        " << gOPT.output_wht      << endl;
    cout << "PE_OUTPUT_DAT        " << gOPT.output_dat      << endl;
    cout << "PE_OUTPUT_MED        " << gOPT.output_med      << endl;
    cout << "PE_OUTPUT_P50        " << gOPT.output_p50      << endl;
    cout << "PE_OUTPUT_P90        " << gOPT.output_p90      << endl;
    cout << "PE_OUTPUT_AVR        " << gOPT.output_avr      << endl;
    cout << "PE_OUTPUT_RMS        " << gOPT.output_rms      << endl;

    cout << endl;
}

void DumpUserOptions(TString odir, CWB::config* cfg) {

    TString ofName = odir+"/pe_config.txt";

    ofstream out;
    out.open(ofName,ios::out);
    if(!out.good()) {cout << "DumpUserOptions : Error Opening File : " << ofName << endl;exit(1);}

    TString info="";
    TString tabs="\t\t\t\t";

    char version[128];
    sprintf(version,"WAT Version : %s - SVN Revision : %s - Tag/Branch : %s",watversion('f'),watversion('r'),watversion('b'));

    out << endl;
    out << "--------------------------------"     << endl;
    out << "PE version                      "     << endl;
    out << "--------------------------------"     << endl;
    out << endl;
    out << version << endl;
    out << endl;

    out << "--------------------------------"     << endl;
    out << "PE config options               "     << endl;
    out << "--------------------------------"     << endl;
    out << endl;

    out << "pe_trials            " << gOPT.trials << endl;
    info = "// number of trials (def=100)";
    out << tabs << info << endl;

    out << "pe_signal            " << gOPT.signal << endl;
    info = "// signal used for trials : 0(def)/1/2";
    out << tabs << info << endl;
    info = "// 0 - reconstructed waveform";
    out << tabs << info << endl;
    info = "// 1 - injected waveform (available only in simulation mode)";
    out << tabs << info << endl;
    info = "// 2 - original waveform = (reconstructed+null)";
    out << tabs << info << endl;

    out << "pe_noise             " << gOPT.noise << endl;
    info = "// noise used for trials :  0/1(def)/2 ";
    out << tabs << info << endl;
    info = "// 0 - waveform is injected in the whitened HoT and apply Coherent+SuperCluster+Likelihood stages";
    out << tabs << info << endl;
    info = "// 1 - add gaussian noise to likelihood sparse maps and apply Likelihood stage";
    out << tabs << info << endl;
    info = "// 2 - add whitened HoT to likelihood sparse maps and apply Likelihood stage";
    out << tabs << info << endl;

    for(int n=0;n<cfg->nIFO;n++) {
      out << "pe_amp_cal_err  " << cfg->ifo[n] << "   " << gOPT.amp_cal_err[n] << endl;
    }
    info = "// max percentage of amplitude miscalibration : def(0) -> disabled";
    out << tabs << info << endl;
    info = "// if>0 -> uniform in [1-pe_amp_cal_err,1+pe_amp_cal_err] else gaus(1,pe_amp_cal_err)";
    out << tabs << info << endl;

    for(int n=0;n<cfg->nIFO;n++) {
      out << "pe_phs_cal_err  " << cfg->ifo[n] << "   " << gOPT.phs_cal_err[n] << endl;
    }
    info = "// max phase (degrees) miscalibration : def(0) -> disabled";
    out << tabs << info << endl;
    info = "// if>0 -> uniform in [-pe_phs_cal_err,+pe_phs_cal_err] else gaus(pe_phs_cal_err)";
    out << tabs << info << endl;

    out << "pe_multitask         " << (gOPT.multitask ? "enable" : "disable") << endl;
    info = "// if enabled, PE trials are executed by different jobs in multitask : only simulation=0 and single event (gps>0)";
    out << tabs << info << endl;
    
    out << "pe_ced_dump          " << gOPT.ced_dump << endl;
    info = "// dump CED at gLRETRY=PE_CED_DUMP (def(-1) -> disabled)";
    out << tabs << info << endl;

    out << "pe_skymask           " << gOPT.skymask << endl;
    info = "// skymask used for trials : 0(def)/1/2/3 ";
    out << tabs << info << endl;
    info = "// 0 - disable -> search over all sky";
    out << tabs << info << endl;
    out << tabs << "// 1 - skymask with radius " << SKYMASK_RADIUS << "(deg) and source selected according the reconstructed skymap probability " << endl;
    info = "// 2 - skymask select only the sky position where the waveform is reconstructed (the maximum of detection statistic)";
    out << tabs << info << endl;
    info = "// 3 - skymask select only the sky position where the waveform has been injected";
    out << tabs << info << endl;

    out << "pe_seed              " << gOPT.seed << endl;
    info = "// seed used by PE for random generation - 0(def) -> random seed";
    out << tabs << info << endl;

    out.precision(14);
    out << "pe_gps               " << gOPT.gps << endl;
    info = "// if >0 only gps +/- iwindow is analyzed - 0(def) -> disabled";
    out << tabs << info << endl;
    out.precision(6);

    out << endl;
    out << "pe_ced tfmap         " << (gOPT.ced_tfmap ? "enable" : "disable") << endl;
    info = "// tfmap   - pe_ced_enable(def)/disable : Shows/Hides the Time-Frequency Maps Tab -> Spectrograms, Scalograms, TF Likelihood,Null";
    out << tabs << info << endl;
    out << "pe_ced rdr           " << (gOPT.ced_rdr ? "enable" : "disable") << endl;
    info = "// rdr     - pe_ced_enable(def)/disable : Shows/Hides Reconstructed Detector Response Tab";
    out << tabs << info << endl;
    out << "pe_ced psd           " << (gOPT.ced_psd ? "enable" : "disable") << endl;
    info = "// psd     - pe_ced_enable(def)/disable : Shows/Hides Power Spectral Density Tab";
    out << tabs << info << endl;
    out << "pe_ced skymap        " << (gOPT.ced_skymap ? "enable" : "disable") << endl;
    info = "// skymap  - pe_ced_enable(def)/disable : Shows/Hides SkyMaps Tab";
    out << tabs << info << endl;
    out << "pe_ced rec           " << (gOPT.ced_rec ? "enable" : "disable") << endl;
    info = "// rec     - pe_ced_enable(def)/disable : Shows/Hides Reconstructed Waveforms/Instantaneous-Frequency with errors Tab";
    out << tabs << info << endl;
    out << "pe_ced inj           " << (gOPT.ced_inj ? "enable" : "disable") << endl;
    info = "// inj     - pe_ced_enable(def)/disable : Showsi/Hides Injection' tab reports the comparison with the injected whitened waveform";
    out << tabs << info << endl;
    out << "pe_ced rinj          " << (gOPT.ced_rinj ? "enable" : "disable") << endl;
    info = "// rinj    - pe_ced_enable/disable(def) : Shows/Hides Injection' tab reports the comparison with the injected whitened waveform";
    out << tabs << info << endl;
    info = "                                          in time-frequency subspace of the PointEstimated waveform";
    out << tabs << info << endl;
    out << "pe_ced cm            " << (gOPT.ced_cm ? "enable" : "disable") << endl;
    info = "// cm      - pe_ced_enable(def)/disable : Shows/Hides Chirp Mass Value/Error Distributions Tab";
    out << tabs << info << endl;
    out << "pe_ced distr         " << (gOPT.ced_distr ? "enable" : "disable") << endl;
    info = "// distr   - pe_ced_enable/disable(def) : Shows/Hides Residuals Distributions Tab";
    out << tabs << info << endl;
    out << "pe_ced null          " << (gOPT.ced_null ? "enable" : "disable") << endl;
    info = "// null    - pe_ced_enable/disable(def) : Shows/Hides Null Pixels Distributions Tab";
    out << tabs << info << endl;
    out << "pe_ced pca           " << (gOPT.ced_pca ? "enable" : "disable") << endl;
    info = "// pca     - pe_ced_enable/disable(def) : Shows/Hides PCA likelihood TF plot in the Time-Frequency Maps Tab";
    out << tabs << info << endl;

    out.close();
}

void ResetUserOptions() {

    gOPT.trials          = PE_TRIALS;
    gOPT.signal          = PE_SIGNAL;
    gOPT.noise           = PE_NOISE;
    for(int n=0;n<NIFO_MAX;n++) gOPT.amp_cal_err[n] = PE_AMP_CAL_ERR;
    for(int n=0;n<NIFO_MAX;n++) gOPT.phs_cal_err[n] = PE_PHS_CAL_ERR;
    gOPT.skymask         = PE_SKYMASK;
    gOPT.seed            = PE_SEED;
    gOPT.gps             = PE_GPS;

    gOPT.multitask       = PE_MULTITASK;

    gOPT.ced_dump        = PE_CED_DUMP;

    gOPT.ced_tfmap       = PE_CED_TFMAP;
    gOPT.ced_rdr         = PE_CED_RDR;
    gOPT.ced_psd         = PE_CED_PSD;
    gOPT.ced_skymap      = PE_CED_SKYMAP;
    gOPT.ced_rec         = PE_CED_REC;
    gOPT.ced_inj         = PE_CED_INJ;
    gOPT.ced_rinj        = PE_CED_rINJ;
    gOPT.ced_cm          = PE_CED_CM;
    gOPT.ced_distr       = PE_CED_DISTR;
    gOPT.ced_null        = PE_CED_NULL;
    gOPT.ced_pca         = PE_CED_PCA;

    gOPT.output_inj	 = PE_OUTPUT_INJ;
    gOPT.output_rec	 = PE_OUTPUT_REC;
    gOPT.output_wht	 = PE_OUTPUT_WHT;
    gOPT.output_dat	 = PE_OUTPUT_DAT;
    gOPT.output_med	 = PE_OUTPUT_MED;
    gOPT.output_p50	 = PE_OUTPUT_P50;
    gOPT.output_p90	 = PE_OUTPUT_P90;
    gOPT.output_avr	 = PE_OUTPUT_AVR;
    gOPT.output_rms	 = PE_OUTPUT_RMS;
}

void SetOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, bool dump_pe) {

   // import slagShift
   float* gSLAGSHIFT=NULL; IMPORT(float*,gSLAGSHIFT)

   int nIFO = NET->ifoListSize();                       // number of detectors

   // search output root file in the system list
   TFile* froot = NULL;                         
   TList *files = (TList*)gROOT->GetListOfFiles();
   gOUTPUT="";                                    
   if (files) {                                   
     TIter next(files);                           
     TSystemFile *file;                           
     TString fname;                               
     bool check=false;                            
     while ((file=(TSystemFile*)next())) {        
        fname = file->GetName();                  
        // set output root file as the current file
        if(fname.Contains("wave_")) {              
          froot=(TFile*)file;froot->cd();          
          gOUTPUT=fname;                           
          gOUTPUT.ReplaceAll(".root.tmp",".txt");  
          //cout << "output file name : " << fname << endl;
        }                                                  
     }                                                               
     if(!froot) {                                                    
       cout << "CWB_Plugin_Boostrap.C : Error - output root file not found" << endl;
       gSystem->Exit(1);                                                                             
     }                                                                                               
   } else {                                                                                          
     cout << "CWB_Plugin_Boostrap.C : Error - output root file not found" << endl;                   
     gSystem->Exit(1);                                                                               
   }                                                                                                 

   if(dump_pe) { 
     if(gOPT.output_inj) if(cfg->simulation) for(int n=0;n<nIFO;n++) gPE.wINJ[n] = &vINJ[n];
     if(gOPT.output_rec) for(int n=0;n<nIFO;n++) gPE.wREC[n] = &vREC[n];
     if(gOPT.output_wht) for(int n=0;n<nIFO;n++) gPE.wWHT[n] = &vWHT[n];
     if(gOPT.output_dat) for(int n=0;n<nIFO;n++) gPE.wDAT[n] = &vDAT[n];
     if(gOPT.output_med) for(int n=0;n<nIFO;n++) gPE.wMED[n] = &vMED[n];
     if(gOPT.output_p50) for(int n=0;n<nIFO;n++) gPE.wL50[n] = &vL50[n];
     if(gOPT.output_p50) for(int n=0;n<nIFO;n++) gPE.wU50[n] = &vU50[n];
     if(gOPT.output_p90) for(int n=0;n<nIFO;n++) gPE.wL90[n] = &vL90[n];
     if(gOPT.output_p90) for(int n=0;n<nIFO;n++) gPE.wU90[n] = &vU90[n];
     if(gOPT.output_avr) for(int n=0;n<nIFO;n++) gPE.wAVR[n] = &vAVR[n];
     if(gOPT.output_rms) for(int n=0;n<nIFO;n++) gPE.wRMS[n] = &vRMS[n];
   }

   gTREE = (TTree *) froot->Get("waveburst");
   if(gTREE!=NULL) {                         
     EVT = new netevent(gTREE,nIFO);         
     if(dump_pe) { 
       gTREE->Branch("pe_trials",&gPE.trials,"pe_trials/I");
       gTREE->Branch("pe_erR",gPE.erR,TString::Format("pe_erR[%i]/F",11));                           
       gTREE->Branch("pe_erF",gPE.erF,TString::Format("pe_erF[%i]/F",11));                           
       gTREE->Branch("pe_nstat",gPE.nstat,TString::Format("pe_nstat[%i]/F",2*7*cfg->nIFO));
       gTREE->Branch("pe_snet",gPE.snet,TString::Format("pe_snet[%i]/F",2));
       gTREE->Branch("pe_ff",gPE.ff,TString::Format("pe_ff[%i]/F",2));
       gTREE->Branch("pe_of",gPE.of,TString::Format("pe_of[%i]/F",2));
       gTREE->Branch("pe_mch",gPE.mch,TString::Format("pe_mch[%i]/F",2));

       for(int n=0;n<nIFO;n++) {
         if(gOPT.output_inj) if(cfg->simulation) gTREE->Branch(TString::Format("pe_wINJ_%d",n).Data(),"wavearray<double>",&gPE.wINJ[n],32000,0);
         if(gOPT.output_rec) gTREE->Branch(TString::Format("pe_wREC_%d",n).Data(),"wavearray<double>",&gPE.wREC[n],32000,0);
         if(gOPT.output_wht) gTREE->Branch(TString::Format("pe_wWHT_%d",n).Data(),"wavearray<double>",&gPE.wWHT[n],32000,0);
         if(gOPT.output_dat) gTREE->Branch(TString::Format("pe_wDAT_%d",n).Data(),"wavearray<double>",&gPE.wDAT[n],32000,0);
         if(gOPT.output_med) gTREE->Branch(TString::Format("pe_wMED_%d",n).Data(),"wavearray<double>",&gPE.wMED[n],32000,0);
         if(gOPT.output_p50) gTREE->Branch(TString::Format("pe_wL50_%d",n).Data(),"wavearray<double>",&gPE.wL50[n],32000,0);
         if(gOPT.output_p50) gTREE->Branch(TString::Format("pe_wU50_%d",n).Data(),"wavearray<double>",&gPE.wU50[n],32000,0);
         if(gOPT.output_p90) gTREE->Branch(TString::Format("pe_wL90_%d",n).Data(),"wavearray<double>",&gPE.wL90[n],32000,0);
         if(gOPT.output_p90) gTREE->Branch(TString::Format("pe_wU90_%d",n).Data(),"wavearray<double>",&gPE.wU90[n],32000,0);
         if(gOPT.output_avr) gTREE->Branch(TString::Format("pe_wAVR_%d",n).Data(),"wavearray<double>",&gPE.wAVR[n],32000,0);
         if(gOPT.output_rms) gTREE->Branch(TString::Format("pe_wRMS_%d",n).Data(),"wavearray<double>",&gPE.wRMS[n],32000,0);
       }

       if(gOPT.multitask&&(gMTRIAL!=gOPT.trials)) {	// save reconstricted waveforms and skyprob when multitask mode
         for(int n=0;n<nIFO;n++) gTREE->Branch(TString::Format("mtpe_wREC_%d",n).Data(),"wavearray<double>",&wREC[0][n],32000,0);
         gTREE->Branch("mtpe_skyprob","skymap",&wSKYPROB[0],  32000,0);
       }
     }
   } else {                                                                      
     EVT = new netevent(nIFO);                                                   
     gTREE = EVT->setTree();                                                  
   }                                                                                        
   EVT->setSLags(gSLAGSHIFT);                          // set slags into netevent           
   EVT->Psave=cfg->Psave;                                                                   
}                                                                                           

void DumpOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor) {

   if(cfg->dump) EVT->dopen(gOUTPUT.Data(),const_cast<char*>("a"),false);        
   EVT->output2G(gTREE,NET,ID,k,factor);             // get reconstructed parameters
   if(cfg->dump) EVT->dclose();
   if(!cfg->cedDump) NET->getwc(k)->sCuts[ID-1]=1;   // mark as processed
}

void ClearVectors() {

   while(!vINJ.empty()) vINJ.pop_back();
   vINJ.clear(); std::vector<wavearray<double> >().swap(vINJ);
   while(!vREC.empty()) vREC.pop_back();
   vREC.clear(); std::vector<wavearray<double> >().swap(vREC);
   while(!vWHT.empty()) vWHT.pop_back();
   vWHT.clear(); std::vector<wavearray<double> >().swap(vWHT);
   while(!vPCA.empty()) vPCA.pop_back();
   vPCA.clear(); std::vector<wavearray<double> >().swap(vPCA);
   while(!vDAT.empty()) vDAT.pop_back();
   vDAT.clear(); std::vector<wavearray<double> >().swap(vDAT);
   while(!vNUL.empty()) vNUL.pop_back();
   vNUL.clear(); std::vector<wavearray<double> >().swap(vNUL);
   while(!cINJ.empty()) cINJ.pop_back();
   cINJ.clear(); std::vector<wavearray<double> >().swap(cINJ);

   while(!vAVR.empty()) vAVR.pop_back();
   vAVR.clear(); std::vector<wavearray<double> >().swap(vAVR);
   while(!vRMS.empty()) vRMS.pop_back();
   vRMS.clear(); std::vector<wavearray<double> >().swap(vRMS);

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

   while(!fREC.empty()) fREC.pop_back();
   fREC.clear(); std::vector<wavearray<double> >().swap(fREC);
   while(!fINJ.empty()) fINJ.pop_back();
   fINJ.clear(); std::vector<wavearray<double> >().swap(fINJ);
   while(!fAVR.empty()) fAVR.pop_back();
   fAVR.clear(); std::vector<wavearray<double> >().swap(fAVR);
   while(!fRMS.empty()) fRMS.pop_back();
   fRMS.clear(); std::vector<wavearray<double> >().swap(fRMS);

   while(!fMED.empty()) fMED.pop_back();
   fMED.clear(); std::vector<wavearray<double> >().swap(fMED);
   while(!fL50.empty()) fL50.pop_back();
   fL50.clear(); std::vector<wavearray<double> >().swap(fL50);
   while(!fU50.empty()) fU50.pop_back();
   fU50.clear(); std::vector<wavearray<double> >().swap(fU50);
   while(!fL90.empty()) fL90.pop_back();
   fL90.clear(); std::vector<wavearray<double> >().swap(fL90);
   while(!fU90.empty()) fU90.pop_back();
   fU90.clear(); std::vector<wavearray<double> >().swap(fU90);

   while(!vRES.empty()) vRES.pop_back();
   vRES.clear(); std::vector<double >().swap(vRES);
   while(!fRES.empty()) fRES.pop_back();
   fRES.clear(); std::vector<double >().swap(fRES);

   for(int n=0;n<MAX_TRIALS;n++) {
     while(!wREC[n].empty()) wREC[n].pop_back();
     wREC[n].clear(); std::vector<wavearray<double> >().swap(wREC[n]);
   }

   for(int n=0;n<NIFO_MAX;n++) {
     while(!iSNR[n].empty()) iSNR[n].pop_back();
     iSNR[n].clear(); std::vector<double>().swap(iSNR[n]);
     while(!oSNR[n].empty()) oSNR[n].pop_back();
     oSNR[n].clear(); std::vector<double>().swap(oSNR[n]);
     while(!ioSNR[n].empty()) ioSNR[n].pop_back();
     ioSNR[n].clear(); std::vector<double>().swap(ioSNR[n]);
   }

   for(int n=0;n<6;n++) {
     while(!vCHIRP[n].empty()) vCHIRP[n].pop_back();
     vCHIRP[n].clear(); std::vector<double>().swap(vCHIRP[n]);
   }

   while(!vLIKELIHOOD.empty()) vLIKELIHOOD.pop_back();
   vLIKELIHOOD.clear(); std::vector<double >().swap(vLIKELIHOOD);

   for(int n=0;n<NIFO_MAX;n++) {
     while(!aNUL[n].empty()) aNUL[n].pop_back();
     aNUL[n].clear(); std::vector<double>().swap(aNUL[n]);
     while(!ANUL[n].empty()) ANUL[n].pop_back();
     ANUL[n].clear(); std::vector<double>().swap(ANUL[n]);
   }

   for(int n=0;n<NIFO_MAX;n++) {
     while(!aNSE[n].empty()) aNSE[n].pop_back();
     aNSE[n].clear(); std::vector<double>().swap(aNSE[n]);
     while(!ANSE[n].empty()) ANSE[n].pop_back();
     ANSE[n].clear(); std::vector<double>().swap(ANSE[n]);
   }

   for(int n=0;n<NIFO_MAX;n++) {
     while(!vSS[n].empty()) vSS[n].pop_back();
     vSS[n].clear(); std::vector<SSeries<double> >().swap(vSS[n]);
     while(!rSS[n].empty()) rSS[n].pop_back();
     rSS[n].clear(); std::vector<SSeries<double> >().swap(rSS[n]);
     while(!jSS[n].empty()) jSS[n].pop_back();
     jSS[n].clear(); std::vector<SSeries<double> >().swap(jSS[n]);
     while(!dSS[n].empty()) dSS[n].pop_back();
     dSS[n].clear(); std::vector<SSeries<double> >().swap(dSS[n]);
  }
}

void PlotWaveforms(network* NET, CWB::config* cfg, int ID, TString pdir) {

   int nIFO = NET->ifoListSize();                      // number of detectors

   wavearray<double> wDIF,wALI,wENV;

   char title[256]; char ofname[256];
   for(int n=0; n<nIFO; n++) {

     // COMPUTE SN, TString pdirR
/*
     double SNR=0;
     for(int j=0;j<vREC[n].size();j++) SNR+=pow(vREC[n].data[j],2);
     cout << "SNR " << NET->ifoName[n] << " " << sqrt(SNR) << endl;
     double NUL=0;
     for(int j=0;j<vNUL[n].size();j++) NUL+=pow(vNUL[n].data[j],2);
     cout << "NUL " << NET->ifoName[n] << " " << sqrt(NUL) << endl;
*/
/*
     // PLOT -> SPARSE MAP
     PlotSparse(n, NET, cfg, ID, &vINJ[n]);
*/
/*
     // PLOT -> INJ : PCA : REC
     sprintf(title,"%s (TIME) : vINJ(red) - vPCA(blue) - vREC(green)",NET->ifoName[n]);
     sprintf(ofname,"%s_vINJ_vPCA_vREC_time_id_%d.%s",NET->ifoName[n],ID,PLOT_TYPE);
     PlotWaveform(ofname, title, cfg, &vINJ[n], &vPCA[n], &vREC[n],NULL, false,pdir);   // time
     sprintf(title,"%s (FFT) : vINJ(red) - vPCA(blue) - vREC(green)",NET->ifoName[n]);
     //sprintf(ofname,"%s_vINJ_vPCA_vREC_fft_id_%d.%s",NET->ifoName[n],ID,PLOT_TYPE);
     //PlotWaveform(ofname, title, cfg, &vINJ[n], &vPCA[n], &vREC[n],NULL, true,pdir);    // fft
*/

#ifdef SET_WAVEFORM_CUTS
     // PLOT -> cINJ : cINJ-vARV : vRMS
     wDIF = GetDifWaveform(&cINJ[n], &vAVR[n]);
     sprintf(title,"%s (TIME) : cINJ(red) - cINJ-vAVR(blue) - vRMS(green)",NET->ifoName[n]);
     sprintf(ofname,"%s_cINJ_cINJvAVR_vRMS_time.%s",NET->ifoName[n],PLOT_TYPE);
     PlotWaveform(ofname, title, cfg, &cINJ[n], &wDIF, &vRMS[n], NULL, false, pdir);
#endif
/*
     // PLOT -> cINJ : wREC 
     sprintf(title,"%s (TIME) : cINJ(red) - wREC(blue)",NET->ifoName[n]);
     sprintf(ofname,"%s_cINJ_wREC_time.%s",NET->ifoName[n],PLOT_TYPE);
     PlotWaveform(ofname, title, cfg, &cINJ[n], &wREC[0][n], NULL, NULL, false, pdir);
*/

#ifdef SET_WAVEFORM_CUTS
     // PLOT -> cINJ : vREC : vAVR
     sprintf(title,"%s (TIME) : cINJ(red) - vREC(blue) - vAVR(green)",NET->ifoName[n]);
     sprintf(ofname,"%s_cINJ_vREC_vAVR_time.%s",NET->ifoName[n],PLOT_TYPE);
     PlotWaveform(ofname, title, cfg, &cINJ[n], &vREC[n], &vAVR[n], NULL, false, pdir);
#endif

#ifdef SET_WAVEFORM_CUTS
     // PLOT -> cINJ : cINJ
     sprintf(title,"%s (TIME) : cINJ(red) - cINJ(blue)",NET->ifoName[n]);
     sprintf(ofname,"%s_cINJ_cINJ_time.%s",NET->ifoName[n],PLOT_TYPE);
     PlotWaveform(ofname, title, cfg, &cINJ[n], &cINJ[n], NULL, NULL, false, pdir);
#endif

     if(gOPT.ced_rec) {
       // PLOT -> vWHT vs vREC 
       // NOTE : this plot overwrite the standard IFO_wf_signal.png/IFO_wf_signal_fft.png/ plots
       //        the band signal is the now the whitend data 
       double scale = sqrt(1<<cfg->levelR);
       vWHT[n]*=1/scale;
       vREC[n]*=1/scale;
       sprintf(ofname,"%s_wf_signal.%s",NET->ifoName[n],PLOT_TYPE);
       PlotWaveform(ofname, "", cfg, &vWHT[n], &vREC[n], NULL, &vREC[n], false, pdir,0.999,(EColor)16,kRed);
       sprintf(ofname,"%s_wf_signal_fft.%s",NET->ifoName[n],PLOT_TYPE);
       PlotWaveform(ofname, "", cfg, &vWHT[n], &vREC[n], NULL, &vREC[n], true, pdir,0.999,(EColor)16,kRed);
       vWHT[n]*=scale;
       vREC[n]*=scale;
     }

     for(int j=0;j<2;j++) {

       TString tag = j==0 ? "time" : "time_zoom";	// file tag
       double    P = j==0 ? 0.995 : 0.9;		// time zoom

       // RECONSTRUCTED 

       if(gOPT.ced_rec) {

         // PLOT -> vAVR(vMED) : vRMS(PERC)
         //wENV = GetWaveformEnvelope(&vREC[n]);
         sprintf(title,"%s (TIME) : vREC (red) - vRMS:2vRMS (grey:light_gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_vREC_vRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
         PlotWaveformAsymmErrors(ofname, "", cfg, &vREC[n], &vMED[n], &vL50[n], &vU50[n], &vL90[n], &vU90[n], &vREC[n], pdir, P);
#else
         PlotWaveformErrors(ofname, "", cfg, &vREC[n], &vAVR[n], &vRMS[n], &vREC[n], pdir, P);
#endif

         // PLOT -> vREC-vAVR(vMED) : vRMS(PERC)
#ifdef PLOT_MEDIAN
         wDIF = GetDifWaveform(&vREC[n], &vMED[n]);
#else
         wDIF = GetDifWaveform(&vREC[n], &vAVR[n]);
#endif
         sprintf(title,"%s (TIME) : vREC-vAVR(red) - 2*vRMS(gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_vRECvAVR_vRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
         PlotWaveformAsymmErrors(ofname, "", cfg, &wDIF, &wDIF, &vL50[n], &vU50[n], &vL90[n], &vU90[n], &vREC[n], pdir, P);
#else
         PlotWaveformErrors(ofname, "", cfg, &wDIF, &wDIF, &vRMS[n], &vREC[n], pdir, P);
#endif

         // PLOT -> fAVR : fRMS(fMED)
         sprintf(title,"%s (FREQ) : fREC (red) - fRMS:2fRMS (grey:light_gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_fREC_fRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
         PlotWaveformAsymmErrors(ofname, "", cfg, &fREC[n], &fMED[n], &fL50[n], &fU50[n], &fL90[n], &fU90[n], &vREC[n], pdir, P, true);
#else
         PlotFrequencyErrors(ofname, "", cfg, &fREC[n], &fAVR[n], &fRMS[n], &vREC[n], pdir, P);
#endif

         // PLOT -> vREC : vDAT-vREC : vRMS
         wDIF = GetDifWaveform(&vDAT[n], &vREC[n]);
         sprintf(title,"%s (TIME) : vREC(red) - vDAT-vREC(blue) - vRMS(green)",NET->ifoName[n]);
         sprintf(ofname,"%s_vREC_vDATvREC_vRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
         PlotWaveform(ofname, "", cfg, &vREC[n], &wDIF, &vRMS[n], &vREC[n], false, pdir, P);
       }

       if(!cfg->simulation) continue; 

       // RECONSTRUCTED vs INJECTED

       if(gOPT.ced_inj) {
         // PLOT -> vINJ : vAVR(vMED) : vRMS(PERC)
         wALI = GetAlignedWaveform(&vINJ[n], &vREC[n]); 
         sprintf(title,"%s (TIME) : vINJ (red) - vRMS:2vRMS (grey:light_gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_vINJ_vRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
	 PlotWaveformAsymmErrors(ofname, "", cfg, &wALI, &vMED[n], &vL50[n], &vU50[n], &vL90[n], &vU90[n], &vREC[n], pdir, P);
#else
         PlotWaveformErrors(ofname, "", cfg, &wALI, &vAVR[n], &vRMS[n], &vREC[n], pdir, P);
#endif

         // PLOT -> fINJ : fAVE(fMED) : fRMS(PERC)
         sprintf(title,"%s (FREQ) : fINJ (red) - fRMS:2fRMS (grey:light_gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_fINJ_fRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
         PlotWaveformAsymmErrors(ofname, "", cfg, &fINJ[n], &fMED[n], &fL50[n], &fU50[n], &fL90[n], &fU90[n], &vREC[n], pdir, P, true);
#else
         PlotFrequencyErrors(ofname, "", cfg, &fINJ[n], &fAVR[n], &fRMS[n], &vREC[n], pdir, P);
#endif

         // PLOT -> vAVR(vMED)-vINJ : vRMS(PERC)
#ifdef PLOT_MEDIAN
         wDIF = GetDifWaveform(&vMED[n], &vINJ[n]);
#else
         wDIF = GetDifWaveform(&vAVR[n], &vINJ[n]);
#endif
         wDIF = GetAlignedWaveform(&wDIF, &vREC[n]); 
         sprintf(title,"%s (TIME) : vAVR-vINJ(red) - 2*vRMS(gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_vAVRvINJ_vRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
         PlotWaveformAsymmErrors(ofname, "", cfg, &wDIF, &wDIF, &vL50[n], &vU50[n], &vL90[n], &vU90[n], &vREC[n], pdir, P);
#else
         PlotWaveformErrors(ofname, "", cfg, &wDIF, &wDIF, &vRMS[n], &vREC[n], pdir, P);
#endif
       }

       // RECONSTRUCTED vs reduced INJECTED

       if(gOPT.ced_rinj) {

         // PLOT -> vAVR(vMED)-cINJ : vRMS(PERC)
#ifdef PLOT_MEDIAN
         wDIF = GetDifWaveform(&vMED[n], &cINJ[n]);
#else
         wDIF = GetDifWaveform(&vAVR[n], &cINJ[n]);
#endif
         wDIF = GetAlignedWaveform(&wDIF, &vREC[n]); 
         sprintf(title,"%s (TIME) : vAVR-cINJ(red) - 2*vRMS(gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_vAVRcINJ_vRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
	 PlotWaveformAsymmErrors(ofname, "", cfg, &wDIF, &wDIF, &vL50[n], &vU50[n], &vL90[n], &vU90[n], &vREC[n], pdir, P);
#else
         PlotWaveformErrors(ofname, "", cfg, &wDIF, &wDIF, &vRMS[n], &vREC[n], pdir, P);
#endif

         // PLOT -> cINJ : vAVR(vMED) : vRMS(PERC)
         sprintf(title,"%s (TIME) : cINJ (red) - vRMS:2vRMS (grey:light_gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_cINJ_vRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
	 PlotWaveformAsymmErrors(ofname, "", cfg, &cINJ[n], &vMED[n], &vL50[n], &vU50[n], &vL90[n], &vU90[n], &vREC[n], pdir, P);
#else
         PlotWaveformErrors(ofname, "", cfg, &cINJ[n], &vAVR[n], &vRMS[n], &vREC[n], pdir, P);
#endif

         // PLOT -> fINJ : fAVR(PERC) : fRMS
         sprintf(title,"%s (FREQ) : fINJ (red) - fRMS:2fRMS (grey:light_gray)",NET->ifoName[n]);
         sprintf(ofname,"%s_fINJ_fRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
#ifdef PLOT_MEDIAN
         PlotWaveformAsymmErrors(ofname, "", cfg, &fINJ[n], &fMED[n], &fL50[n], &fU50[n], &fL90[n], &fU90[n], &vREC[n], pdir, P, true);
#else
         PlotFrequencyErrors(ofname, "", cfg, &fINJ[n], &fAVR[n], &fRMS[n], &vREC[n], pdir, P);
#endif

         // PLOT -> cINJ : vREC : vRMS
         sprintf(title,"%s (TIME) : cINJ(red) - vREC(blue) - vRMS(green)",NET->ifoName[n]);
         sprintf(ofname,"%s_cINJ_vREC_vRMS_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
         PlotWaveform(ofname, "", cfg, &cINJ[n], &vREC[n], &vRMS[n], &vREC[n], false, pdir, P);

         // PLOT -> vINJ : vINJ-cINJ : vINJcINJ
         wDIF = GetDifWaveform(&vINJ[n], &cINJ[n]);
         wDIF = GetAlignedWaveform(&wDIF, &vREC[n]); 
         wALI = GetAlignedWaveform(&vINJ[n], &vREC[n]); 
         sprintf(title,"%s (TIME) : vINJ(red) - cINJ(blue) - vINJ-cINJ(green)",NET->ifoName[n]);
         sprintf(ofname,"%s_vINJ_cINJ_vINJcINJ_%s.%s",NET->ifoName[n],tag.Data(),PLOT_TYPE);
         PlotWaveform(ofname, "", cfg, &wALI, &cINJ[n], &wDIF, &vREC[n], false, pdir, P);
       }
     }
/*
     // FITTING FACTOR
     double FF1 = FittingFactor(&cINJ[n], &vREC[n]);
     double FF2 = FittingFactor(&cINJ[n], &vPCA[n]);
     cout << NET->ifoName[n] << " -> Fitting Factor - inj/rec = " << FF1 << " inj/pca = " << FF2 << endl;
*/
  }
}

void SetSkyMask(network* NET, CWB::config* cfg, double theta, double phi, double radius) {

  // --------------------------------------------------------
  // define SetSkyMask
  // --------------------------------------------------------
  sprintf(cfg->skyMaskFile,"--theta %f --phi %f --radius %f",90-theta,phi,radius);
  cout << endl << "SetSkyMask : " << cfg->skyMaskFile << endl << endl;
  gCWB.SetSkyMask(NET,cfg,cfg->skyMaskFile,'e');

}

void AddNoiseAndCalErrToSparse(network* NET, CWB::config* cfg, char type) {

  if(type!='j' && type!='r' && type!='v') {
    cout << "AddNoiseAndCalErrToSparse - Error : wrong wave type" << endl; exit(1);                                                             
  }                                                                                                                               

  double d2r = TMath::Pi()/180.;

  int m;
  int nIFO = NET->ifoListSize();               // number of detectors
  int nRES = NET->wdmListSize();               // number of frequency resolution levels

  // update vSS maps
  wavearray<double> WF[NIFO_MAX];
  WSeries<double>  w;
  SSeries<double>* p;
  WDM<double>* pwdm = NULL;

  for(int n=0; n<nIFO; n++) {                   // create random time series
     p = NET->getifo(n)->getSTFmap(0);
     WF[n].resize(p->wdm_nSTS);
     WF[n].rate(p->wdm_rate);
     WF[n].start(p->wdm_start);

     int    size = WF[n].size();
     double rate = WF[n].rate();
     int jS,jE,jW,k;
     double dt,df,ts;
     TComplex C;
     switch(gOPT.noise) {
     case 1:					// gaussian noise
       for(int i=0; i<WF[n].size(); i++) {
         WF[n].data[i] = gRandom->Gaus(0,NOISE_SIGMA);
       }
       cout << "AddNoiseAndCalErrToSparse - Info : add gaussian noise to waveform type : " << type << endl;                                                             
       break; 
     case 2:					// data noise
       w = *(NET->getifo(n)->getTFmap());	// get whitened data 
       w.Inverse();				// get level 0
       // apply time shift to input whitened data (integer number of sammples)
       jS = cfg->segEdge*rate;
       jE = size-jS;
       jW = rate*cfg->iwindow/2;
       k  = gRandom->Uniform(jS+jW,jE-jW); 	// select random offset (excludes +/- iwindow/2 around the event)
       WF[n] = w;
       for(int i=jS; i<jE; i++)  { 
         WF[n].data[k++] = w.data[i];
         if(k==jE) k=jS;
       }
       // apply time shift to input whitened data (fraction of dt)
       WF[n].FFTW(1);
       dt = 1./rate;
       df = rate/size;
       ts = gRandom->Uniform(0.,dt);     	// select a random time shift within the sample time range
       for (int j=0;j<size/2;j++) {		
          TComplex X(WF[n].data[2*j],WF[n].data[2*j+1]);
          X=X*C.Exp(TComplex(0.,-2*PI*j*df*ts));   
          WF[n].data[2*j]=X.Re();
          WF[n].data[2*j+1]=X.Im();
       }
       WF[n].FFTW(-1);
       printf("Info : %s add data noise with time shift : %.5f sec to waveform type : %c\n",NET->getifo(n)->Name,(double)k/rate+ts,type);
       break; 
     default:
       WF[n]=0;
       break;
     }
  }

  // fill waveform sparse maps
  for(int n=0; n<nIFO; n++) {
     double A = 1;
     if(gOPT.amp_cal_err[n]) {
       double amp_cal_err = 1;
       if(gOPT.amp_cal_err[n]>0) amp_cal_err = gRandom->Uniform(1-gOPT.amp_cal_err[n],1+gOPT.amp_cal_err[n]);
       else                      amp_cal_err = gRandom->Gaus(1,fabs(gOPT.amp_cal_err[n])); 
       cout << NET->ifoName[n] << " -> amp_cal_err : " << amp_cal_err << endl;
       A = amp_cal_err;
     }
     double C=1,S=0;
     if(gOPT.phs_cal_err[n]) {
       double phs_cal_err = 0;
       if(gOPT.phs_cal_err[n]>0) phs_cal_err = gRandom->Uniform(-gOPT.phs_cal_err[n],gOPT.phs_cal_err[n]);
       else                      phs_cal_err = gRandom->Gaus(0,fabs(gOPT.phs_cal_err[n]));
       cout << NET->ifoName[n] << " -> phs_cal_err : " << phs_cal_err << endl;
       C = cos(-phs_cal_err*d2r);
       S = sin(-phs_cal_err*d2r);
     }
     for(int j=0; j<nRES; j++) {
        w.Forward(WF[n],*(NET->wdmList[j]));
        int k = NET->getifo(n)->getSTFind(w.wrate());   // pointer to sparse TF array
        p = NET->getifo(n)->getSTFmap(k);               // pointer to sparse TF array
        for(int l=0; l<p->sparseMap00.size(); l++) {
           m = p->sparseIndex[l];
           double aa,AA;
           if(type=='r') {                              // reconstructed waveform                
             aa = rSS[n][k].sparseMap00[l];
             AA = rSS[n][k].sparseMap90[l];
           }
           if(type=='j') {                              // injected waveform                
             aa = jSS[n][k].sparseMap00[l];
             AA = jSS[n][k].sparseMap90[l];
           }
           if(type=='v') {                              // original sparse map (reconstructed+noise)                
             aa = vSS[n][k].sparseMap00[l];
             AA = vSS[n][k].sparseMap90[l];
           }
           // add calibration error
           aa*=A; AA*=A;
           // add phase error
           double bb =  C*aa + S*AA;
           double BB = -S*aa + C*AA ;
           // add gaussian noise
           p->sparseMap00[l]=bb+w.data[m];
           p->sparseMap90[l]=BB+w.data[m+w.maxIndex()+1];
        }
     }
  }
  return;
}

void Wave2Sparse(network* NET, CWB::config* cfg, char type) {

  if(type!='j' && type!='r' && type!='v' && type!='d') {
    cout << "Wave2Sparse - Error : wrong wave type" << endl; exit(1);                                                             
  }                                                                                                                               

  // init waveform sparse maps
  int nIFO = NET->ifoListSize();                       // number of detectors
  for(int n=0;n<nIFO;n++) {                                                  
     detector* pD = NET->getifo(n);                                          
     if(type=='v') vSS[n]=pD->vSS;                                      
     if(type=='r') rSS[n]=pD->vSS;                                      
     if(type=='j') jSS[n]=pD->vSS;                                      
     if(type=='d') dSS[n]=pD->vSS;                                      
  }                                                                          
  if(type=='v') return;                                                 

  int nRES = cfg->l_high-cfg->l_low+1;     // number of frequency resolution levels

  // build waveform array
  wavearray<double> WF[NIFO_MAX];
  for(int n=0; n<nIFO; n++) {    
    WF[n].resize(vSS[n][0].wdm_nSTS);                               
    WF[n].rate(vSS[n][0].wdm_rate);                                 
    WF[n].start(vSS[n][0].wdm_start);                               
    WF[n]=0;                                                        
    if(type=='j') {                                              
      int wos = double(vINJ[n].start()-WF[n].start())*WF[n].rate();   
      for (int i=0;i<vINJ[n].size();i++) WF[n][wos+i] = vINJ[n][i];   
    }                                                                 
    if(type=='r') {                                              
      int wos = double(vREC[n].start()-WF[n].start())*WF[n].rate();   
      for (int i=0;i<vREC[n].size();i++) WF[n][wos+i] = vREC[n][i];   
    }                                                                 
    if(type=='d') {                                              
      int wos = double(vDAT[n].start()-WF[n].start())*WF[n].rate();   
      for (int i=0;i<vDAT[n].size();i++) WF[n][wos+i] = vDAT[n][i];   
    }                                                                 

#ifdef SAVE_WAVE2SPARSE_PLOT
    gwavearray<double> gw(&WF[n]);
    gw.Draw(GWAT_TIME);
    watplot* plot = gw.GetWATPLOT();
    TString gfile;
    if(type=='r') gfile=TString::Format("%s/WAVE2SPARSE_REC_%s.root",".",NET->ifoName[n]);
    if(type=='j') gfile=TString::Format("%s/WAVE2SPARSE_INJ_%s.root",".",NET->ifoName[n]);
    if(type=='d') gfile=TString::Format("%s/WAVE2SPARSE_DAT_%s.root",".",NET->ifoName[n]);
    (*plot) >> gfile;
#endif
  }

  // fill waveform sparse maps
  WSeries<double> w;          
  WDM<double>* pwdm = NULL;   
  double r00,r90;             
  for(int n=0; n<nIFO; n++) { 
    for(int i=0; i<nRES; i++) {
       w.Forward(WF[n],*(NET->wdmList[i]));
       int k = NET->getifo(n)->getSTFind(w.wrate());  // pointer to sparse TF array

       int size;
       if(type=='r') size = rSS[n][k].sparseMap00.size();
       if(type=='j') size = jSS[n][k].sparseMap00.size();
       if(type=='d') size = dSS[n][k].sparseMap00.size();

       for(int j=0; j<size; j++) {
         int index;               
         if(type=='r') index = rSS[n][k].sparseIndex[j];
         if(type=='j') index = jSS[n][k].sparseIndex[j];
         if(type=='d') index = dSS[n][k].sparseIndex[j];
         double v00 = w.pWavelet->pWWS[index];               
         double v90 = w.pWavelet->pWWS[index+w.maxIndex()+1];
         if(type=='r') {                                
           rSS[n][k].sparseMap00[j]=v00;                     
           rSS[n][k].sparseMap90[j]=v90;                     
         }                                                   
         if(type=='j') {                                
           jSS[n][k].sparseMap00[j]=v00;                     
           jSS[n][k].sparseMap90[j]=v90;                     
         }                                                   
         if(type=='d') {                                
           dSS[n][k].sparseMap00[j]=v00;                     
           dSS[n][k].sparseMap90[j]=v90;                     
         }                                                   
       }
    }
  }
}

std::vector<int> ComputeStatisticalErrors(network* NET, CWB::config* cfg, int ID) {

  int nIFO = NET->ifoListSize();                // number of detectors
  std::vector<int>  nrec(nIFO);			// number of averaged 

#ifdef SET_WAVEFORM_CUTS
  // get time and frequency ranges
  gBT=1e20; gET=0;
  gBF=1e20; gEF=0;
  for(int n=0;n<nIFO;n++) { 
    double bT,eT,bF,eF;                      
    GetTimeBoundaries(vREC[n], EFRACTION_CUT, bT, eT);
    GetFrequencyBoundaries(vREC[n], EFRACTION_CUT, bF, eF);
    if(bT<gBT) gBT=bT; if(eT>gET) gET=eT;
    if(bF<gBF) gBF=bF; if(eF>gEF) gEF=eF;
  }
  cout << endl;
  cout << "CUT TIME RANGE : " << gBT-gBT << " " << gET-gBT << endl;
  cout << "CUT FREQ RANGE : " << gBF << " " << gEF << endl;
  cout << endl;
#endif

  // copy cINJ (only the vREC time range) into winj
  wavearray<double> winj[NIFO_MAX];                
  wavearray<double> finj[NIFO_MAX];                
  if(cfg->simulation) {                            
    for(int n=0;n<nIFO;n++) {                      
      winj[n] = GetAlignedWaveform(&cINJ[n], &vREC[n]); 
      finj[n] = CWB::Toolbox::getHilbertIFrequency(winj[n]); 
#ifdef SET_WAVEFORM_CUTS
      SetWaveformCuts(&winj[n], gBT, gET, gBF, gEF);
      cINJ.push_back(GetAlignedWaveform(&winj[n], &vREC[n])); 
#endif
    }                                                                    
  }                                                                      

  // compute average of wREC (only the vREC time range) into wavr
  wavearray<double> wrec[MAX_TRIALS][NIFO_MAX];                
  wavearray<double> frec[MAX_TRIALS][NIFO_MAX];                
  for(int n=0;n<nIFO;n++) {                      
    nrec[n]=0;
    for(int i=0;i<gOPT.trials;i++) {            
      if(wREC[i].size()) {	// event detected
        wrec[i][n] = GetAlignedWaveform(&wREC[i][n], &vREC[n]); 
        frec[i][n] = CWB::Toolbox::getHilbertIFrequency(wrec[i][n]); 
#ifdef SET_WAVEFORM_CUTS
        SetWaveformCuts(&wrec[i][n], gBT, gET, gBF, gEF);
#endif
        nrec[n]++;
      } else {			// event rejected (set to 0)
        wrec[i][n] = vREC[n]; wrec[i][n]=0;
        frec[i][n] = wrec[i][n]; 
      }
    } 
    if(nrec[n]==0) {nrec.clear();return nrec;}	// return if no events are detects
  }                                                                    

  gPE.trials = nrec[0];				// save number of effective trials
  cout << endl << "ComputeStatisticalErrors - detected/trials : " << nrec[0] << "/" << gOPT.trials << endl << endl;

  // compute vAVR & vRMS
  wavearray<double> wavr[NIFO_MAX];                
  wavearray<double> wrms[NIFO_MAX];                
  for(int n=0;n<nIFO;n++) {                      
    wrms[n] = wrec[0][n];  wrms[n] = 0;  
    wavr[n] = wrec[0][n];  wavr[n] = 0;  
    for(int i=0;i<gOPT.trials;i++) {            
      wavr[n] += wrec[i][n];  
      for(int j=0;j< wrec[i][n].size();j++) wrms[n][j] += wrec[i][n][j]*wrec[i][n][j];  
    }
    wavr[n] *= 1./nrec[n];  

    wrms[n] *= 1./nrec[n];  
    for(int j=0;j< wrms[n].size();j++) {
      wrms[n][j] -= wavr[n][j]*wavr[n][j];  
      wrms[n][j]  = sqrt(wrms[n][j]);  
    }

    vAVR.push_back(wavr[n]);
    vRMS.push_back(wrms[n]);
  }                                                                    

  // compute vMED, vL50, vU50, vL90, vU90
  wavearray<double> wmed[NIFO_MAX];                
  wavearray<double> wl50[NIFO_MAX];                
  wavearray<double> wu50[NIFO_MAX];                
  wavearray<double> wl90[NIFO_MAX];                
  wavearray<double> wu90[NIFO_MAX];                
  for(int n=0;n<nIFO;n++) {                      
    wmed[n] = wrec[0][n];  wmed[n] = 0;  
    wl50[n] = wrec[0][n];  wl50[n] = 0;  
    wu50[n] = wrec[0][n];  wu50[n] = 0;  
    wl90[n] = wrec[0][n];  wl90[n] = 0;  
    wu90[n] = wrec[0][n];  wu90[n] = 0;  

    int ntry = gPE.trials;			// number of detected events in the trials
    int *index = new int[ntry];
    float *value = new float[ntry];
    for(int j=0;j<wrec[0][n].size();j++) {

      int k=0; for(int i=0;i<gOPT.trials;i++) if(wREC[i].size()) value[k++] = wrec[i][n][j];  // select detected events
      TMath::Sort(ntry,value,index,false);

      int imed = (ntry*50.)/100.; if(imed>=ntry) imed=ntry-1;
      wmed[n][j] = value[index[imed]];

      int il50 = (ntry*25.)/100.; if(il50>=ntry) il50=ntry-1;
      int iu50 = (ntry*75.)/100.; if(iu50>=ntry) iu50=ntry-1;
      int il90 = (ntry*5.)/100.;  if(il90>=ntry) il90=ntry-1;
      int iu90 = (ntry*95.)/100.; if(iu90>=ntry) iu90=ntry-1;

      wl50[n][j] = wmed[n][j]>0 ? wmed[n][j]-value[index[il50]] : value[index[iu50]]-wmed[n][j];
      wu50[n][j] = wmed[n][j]>0 ? value[index[iu50]]-wmed[n][j] : wmed[n][j]-value[index[il50]];
      wl90[n][j] = wmed[n][j]>0 ? wmed[n][j]-value[index[il90]] : value[index[iu90]]-wmed[n][j];
      wu90[n][j] = wmed[n][j]>0 ? value[index[iu90]]-wmed[n][j] : wmed[n][j]-value[index[il90]];
    }        
    delete [] index;
    delete [] value;

    vMED.push_back(wmed[n]);
    vL50.push_back(wl50[n]);
    vU50.push_back(wu50[n]);
    vL90.push_back(wl90[n]);
    vU90.push_back(wu90[n]);
  }                                                                    

  // compute fAVR & fRMS
  wavearray<double> favr[NIFO_MAX];                
  wavearray<double> frms[NIFO_MAX];                
  for(int n=0;n<nIFO;n++) {                      
    frms[n] = frec[0][n];  frms[n] = 0;  
    favr[n] = frec[0][n];  favr[n] = 0;  
    for(int i=0;i<gOPT.trials;i++) {            
      favr[n] += frec[i][n];  
      for(int j=0;j< frec[i][n].size();j++) frms[n][j] += frec[i][n][j]*frec[i][n][j];  
    }
    favr[n] *= 1./nrec[n];  

    frms[n] *= 1./nrec[n];  
    for(int j=0;j< frms[n].size();j++) {
      frms[n][j] -= favr[n][j]*favr[n][j];  
      frms[n][j]  = sqrt(frms[n][j]);  
    }

    fAVR.push_back(favr[n]);
    fRMS.push_back(frms[n]);
  }                                                                    

  // compute fMED, fL50, fU50, fL90, fU90
  wavearray<double> fmed[NIFO_MAX];                
  wavearray<double> fl50[NIFO_MAX];                
  wavearray<double> fu50[NIFO_MAX];                
  wavearray<double> fl90[NIFO_MAX];                
  wavearray<double> fu90[NIFO_MAX];                
  for(int n=0;n<nIFO;n++) {                      
    fmed[n] = frec[0][n];  fmed[n] = 0;  
    fl50[n] = frec[0][n];  fl50[n] = 0;  
    fu50[n] = frec[0][n];  fu50[n] = 0;  
    fl90[n] = frec[0][n];  fl90[n] = 0;  
    fu90[n] = frec[0][n];  fu90[n] = 0;  

    int ntry = gPE.trials;			// number of detected events in the trials
    int *index = new int[ntry];
    float *value = new float[ntry];
    for(int j=0;j<frec[0][n].size();j++) {

      int k=0; for(int i=0;i<gOPT.trials;i++) if(wREC[i].size()) value[k++] = frec[i][n][j];  // select detected events
      TMath::Sort(ntry,value,index,false);

      int imed = (ntry*50.)/100.; if(imed>=ntry) imed=ntry-1;
      fmed[n][j] = value[index[imed]];

      int il50 = (ntry*25.)/100.; if(il50>=ntry) il50=ntry-1;
      int iu50 = (ntry*75.)/100.; if(iu50>=ntry) iu50=ntry-1;
      int il90 = (ntry*5.)/100.;  if(il90>=ntry) il90=ntry-1;
      int iu90 = (ntry*95.)/100.; if(iu90>=ntry) iu90=ntry-1;

      fl50[n][j] = fmed[n][j]>0 ? fmed[n][j]-value[index[il50]] : value[index[iu50]]-fmed[n][j];
      fu50[n][j] = fmed[n][j]>0 ? value[index[iu50]]-fmed[n][j] : fmed[n][j]-value[index[il50]];
      fl90[n][j] = fmed[n][j]>0 ? fmed[n][j]-value[index[il90]] : value[index[iu90]]-fmed[n][j];
      fu90[n][j] = fmed[n][j]>0 ? value[index[iu90]]-fmed[n][j] : fmed[n][j]-value[index[il90]];
    }        
    delete [] index;
    delete [] value;

    fMED.push_back(fmed[n]);
    fL50.push_back(fl50[n]);
    fU50.push_back(fu50[n]);
    fL90.push_back(fl90[n]);
    fU90.push_back(fu90[n]);
  }                                                                    

  // ---------------------------------------------------
  // compute residuals energy
  // ---------------------------------------------------
  double eavr=0;		// energy avr^2
  for(int n=0;n<nIFO;n++) { 
    for(int j=0;j< wavr[n].size();j++) eavr += pow(wavr[n][j],2);  
  }
  for(int i=0;i<gOPT.trials;i++) {            
    if(wREC[i].size()) {	// event detected
      double eres=0;		// residual energy (rec-avr)^2
      for(int n=0;n<nIFO;n++) { 
        for(int j=0;j< wavr[n].size();j++) {
          eres += pow(wrec[i][n][j]-wavr[n][j],2);  
        }
      }
      // add normalized residual energy of wrec to vRES
      vRES.push_back(eres/eavr);
    }
  }
  double jres=0;		// residual energy (inj-avr)^2
  if(cfg->simulation) {                            
    for(int n=0;n<nIFO;n++) { 
      for(int j=0;j< wavr[n].size();j++) {
        jres += pow(winj[n][j]-wavr[n][j],2);  
      }
    }
    jres/=eavr;
  } 
  // add normalized residual energy of winj to vRES
  vRES.push_back(jres);

  // compute probability distribution of residuals 
  int size = vRES.size()-1;
  wavearray<int> index(size);
  wavearray<double> eres(size);
  for(int i=0;i<size;i++) eres[i] = vRES[i]; 
  TMath::Sort(size,const_cast<double*>(eres.data),index.data,false);
  gPE.erR[0] = jres;
  for(int i=1;i<10;i++) gPE.erR[i]=eres[index[int(i*size/10.)]];
  gPE.erR[10]=0;
  for(int i=0;i<size;i++)  {
    if(eres[i]<=jres) gPE.erR[10]+=1.;
  }
  gPE.erR[10]/=size;
  // print probability distribution of residuals
  cout.precision(3);
  cout << endl << "gPE.erR : ";
  for(int i=0;i<10;i++) cout << gPE.erR[i] << "/";
  cout << gPE.erR[10] << endl << endl;

  // ---------------------------------------------------
  // compute frequency residuals 
  // ---------------------------------------------------
  eavr=0;			// energy wavr^4*favr^2
  for(int n=0;n<nIFO;n++) { 
    for(int j=0;j< wavr[n].size();j++) eavr += pow(wavr[n][j],4)*pow(favr[n][j],2);  
  }
  for(int i=0;i<gOPT.trials;i++) {            
    if(wREC[i].size()) {	// event detected
      double eres=0;		// weighted frequency residual energy wavr^4*(frec-favr)^2
      for(int n=0;n<nIFO;n++) { 
        for(int j=0;j<wavr[n].size();j++) {
          eres += pow(wavr[n][j],4)*pow(frec[i][n][j]-favr[n][j],2);  
        }
      }
      // add normalized frequency residual energy of frec to fRES
      fRES.push_back(eres/eavr);
    }
  }
  jres=0;			// weighted frequency residual energy wavr^4*(finj-favr)^2
  if(cfg->simulation) {                            
    for(int n=0;n<nIFO;n++) { 
      for(int j=0;j< wavr[n].size();j++) {
        jres += pow(wavr[n][j],4)*pow(finj[n][j]-favr[n][j],2);  
      }
    }
    jres/=eavr;
  } 
  // add normalized frequency residual energy of finj to fRES
  fRES.push_back(jres);

  // compute probability distribution of residuals 
  for(int i=0;i<size;i++) eres[i] = fRES[i]; 
  TMath::Sort(size,const_cast<double*>(eres.data),index.data,false);
  gPE.erF[0] = jres;
  for(int i=1;i<10;i++) gPE.erF[i]=eres[index[int(i*size/10.)]];
  gPE.erF[10]=0;
  for(int i=0;i<size;i++)  {
    if(eres[i]<=jres) gPE.erF[10]+=1.;
  }
  gPE.erF[10]/=size;
  // print probability distribution of frequency residuals
  cout.precision(3);
  cout << endl << "gPE.erF : ";
  for(int i=0;i<10;i++) cout << gPE.erF[i] << "/";
  cout << gPE.erF[10] << endl << endl;

  return nrec;
}

wavearray<double> GetAlignedWaveform(wavearray<double>* wf1, wavearray<double>* wf2) {

   wavearray<double> wf = *wf2;
   wf=0;

   if(wf1==NULL)      return wf;
   if(wf1->size()==0) return wf;

   double R=wf1->rate();

   double b_wf1 = wf1->start();
   double e_wf1 = wf1->start()+wf1->size()/R;
   double b_wf2 = wf2->start();
   double e_wf2 = wf2->start()+wf2->size()/R;

   int o_wf1 = b_wf1>b_wf2 ? 0 : int((b_wf2-b_wf1)*R+0.5);
   int o_wf2 = b_wf1<b_wf2 ? 0 : int((b_wf1-b_wf2)*R+0.5);

   double startXCOR = b_wf1>b_wf2 ? b_wf1 : b_wf2;
   double endXCOR   = e_wf1<e_wf2 ? e_wf1 : e_wf2;
   int sizeXCOR  = int((endXCOR-startXCOR)*R+0.5);

   for(int i=0;i<sizeXCOR;i++) wf[i+o_wf2] = wf1->data[i+o_wf1];

   return wf;
}

wavearray<double> GetWaveformEnvelope(wavearray<double>* wf) {

   wavearray<double> wfq;

   if(wf==NULL)      return wfq;
   if(wf->size()==0) return wfq;

   wfq = *wf;

   CWB::mdc::PhaseShift(wfq,90);
   for(int i=0;i<wf->size();i++) wfq[i] = sqrt(pow(wf->data[i],2)+pow(wfq[i],2));

   return wfq;
}

skymap GetSkyProb(network* NET, int id) {

   std::vector<float>* vP;
   std::vector<int>*   vI;

   skymap skyprob = NET->getifo(0)->tau;
   skyprob = 0;

   // save the probability skymap
   if(NET->wc_List[0].p_Map.size()) {

      vP = &(NET->wc_List[0].p_Map[id-1]);
      vI = &(NET->wc_List[0].p_Ind[id-1]);

      skyprob = NET->getifo(0)->tau;
      skyprob = 0.;                 

      for(int j=0; j<int(vP->size()); j++) {
         int i = (*vI)[j];                  
         double th = skyprob.getTheta(i);
         double ph = skyprob.getPhi(i);  
         int k=skyprob.getSkyIndex(th, ph);
         skyprob.set(k,(*vP)[j]);          
      }                                        
   }                                           

   return skyprob;
}

void SaveSkyProb(network* NET, CWB::config* cfg, int id) {

   std::vector<float>* vP;
   std::vector<int>*   vI;

   // save the probability skymap
   if(NET->wc_List[0].p_Map.size()) {

      vP = &(NET->wc_List[0].p_Map[id-1]);
      vI = &(NET->wc_List[0].p_Ind[id-1]);

      gSKYPROB = NET->getifo(0)->tau;
      gSKYPROB = 0.;                 

      for(int j=0; j<int(vP->size()); j++) {
         int i = (*vI)[j];                  
         double th = gSKYPROB.getTheta(i);
         double ph = gSKYPROB.getPhi(i);  
         int k=gSKYPROB.getSkyIndex(th, ph);
         gSKYPROB.set(k,(*vP)[j]);          
      }                                        
   }                                           

   gHSKYPROB.SetBins(gSKYPROB.size(),0,gSKYPROB.size()-1);

   double PROB=0;
   for(int l=0;l<gSKYPROB.size();l++) PROB+=gSKYPROB.get(l);
   cout << "PROB : " << PROB << endl;                             
   for(int l=0;l<gSKYPROB.size();l++) gHSKYPROB.SetBinContent(l,gSKYPROB.get(l));

/*
   wavearray<int> index(gSKYPROB.size());
   wavearray<float> skyprob(gSKYPROB.size());
   for(int l=0;l<gSKYPROB.size();l++) skyprob[l]=gSKYPROB.get(l);
   TMath::Sort(int(gSKYPROB.size()),const_cast<float*>(skyprob.data),index.data,true);
   double THR = PROB>0.998 ? 0.998 : PROB;                                               
   double prob=0;                                                                        
   int L=0;                                                                              
   skymap SkyMask = gSKYPROB;                                                                
   for(int l=0;l<gSKYPROB.size();l++) {                                               
      prob+=skyprob[index[l]];                                                           
      if(prob>THR && L>1000) SkyMask.set(index[l],0); else {SkyMask.set(index[l],1);L++;}
   }                                                                                     
   cout << "PROB(0.99) : " << prob << " " << L << endl;                                  
   //NET->setSkyMask(SkyMask,'e');                                                       
   //NET->setIndexMode(0);                                                               
*/

#ifdef SAVE_SKYPROB
   // plot gSKYPROB 
   gskymap gSM=gSKYPROB;
   //gSM.SetOptions("hammer","Celestial",2);
   gSM.SetOptions("","Geographic");
   gSM.SetZaxisTitle("probability");        
   gSM.Draw(1);                             
   //TH2D* hsm = gSM.GetHistogram();        
   //hsm->GetZaxis()->SetTitleOffset(0.85); 
   //hsm->GetZaxis()->SetTitleSize(0.03);   
   gSM.Print(TString::Format("%s/gSkyprob_%d.png",".",id)); 
#endif                                                        
}                                                             

void PlotResiduals(int ID, TString pdir, int sim, char type) {

  std::vector<double > xRES = type=='f' ? fRES : vRES;	

  int size = xRES.size()-1;	// last index of xRES contains the residual of injection
  double jres = xRES[size];	// residual of injection
  if(size==0) return;

  double hmin=1e20;
  double hmax=0;
  for(int i=0;i<size;i++)  {if(xRES[i]>hmax) hmax=xRES[i]; if(xRES[i]<hmin) hmin=xRES[i];}
  hmin-=0.4*(hmax-hmin);
  hmax+=0.4*(hmax-hmin);
  hmax*=1.2;if(hmax>1) hmax=1;	

  // plot residuals
  TCanvas* cres = new TCanvas("residuals", "residuals", 200, 20, 600, 600);
  TH1F* hres = new TH1F("residuals","residuals",100,hmin,hmax);
  for(int i=0;i<size;i++)  hres->Fill(xRES[i]);

  // make it cumulative
  double integral = hres->ComputeIntegral();
  if (integral==0) {cout << "PlotResiduals : Empty histogram !!!" << endl;return;}
  double* cumulative = hres->GetIntegral();
  for (int i=0;i<hres->GetNbinsX();i++) hres->SetBinContent(i,cumulative[i]/integral);
  //cout << "integral " << integral << endl;
  hres->SetLineWidth(2);

  hres->Draw();

  hres->GetXaxis()->SetTitle("residuals");
  hres->GetXaxis()->SetTitleOffset(1.4);
  hres->GetXaxis()->SetLabelOffset(0.02);
  hres->GetXaxis()->SetNoExponent();
  hres->GetXaxis()->SetMoreLogLabels();
  if(type=='f') {
    hres->GetXaxis()->SetRangeUser(0.001,1.0);
  } else {
    hres->GetXaxis()->SetRangeUser(0.01,1.0);
  }

  hres->GetYaxis()->SetTitle("probability");
  hres->GetYaxis()->SetTitleOffset(1.6);
  hres->GetYaxis()->SetLabelOffset(0.02);
  hres->GetYaxis()->SetNoExponent();
  hres->GetYaxis()->SetRangeUser(0.01,1.0);

  cres->SetLogx();
  cres->SetLogy();
  cres->SetGridx();
  cres->SetGridy();

  // draw vertical line at the jres value
  float ymax = hres->GetMaximum();
  TLine *line = new TLine(jres,0,jres,ymax);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  if(sim) line->Draw();

  // print probability distribution of residuals
  cout.precision(3);
  if(type=='f') {
    cout << endl << "---------> gPE.erF : ";
    for(int i=0;i<10;i++) cout << gPE.erF[i] << "/";
    cout << gPE.erF[10] << endl << endl;
  } else {
    cout << endl << "---------> gPE.erR : ";
    for(int i=0;i<10;i++) cout << gPE.erR[i] << "/";
    cout << gPE.erR[10] << endl << endl;
  }

  gStyle->SetLineColor(kBlack);
  char title[256];
  if(type=='f') {
    sprintf(title,"Cumulative distribution of freq residuals errors (entries=%d - prob inj = %2.2g)",size,gPE.erF[10]);
  } else { 
    sprintf(title,"Cumulative distribution of residuals errors (entries=%d - prob inj = %2.2g)",size,gPE.erR[10]);
  }
  hres->SetTitle(title);
  hres->SetStats(kFALSE);
  if(type=='f') {
    cres->Print(TString::Format("%s/fresiduals_errors.%s",pdir.Data(),PLOT_TYPE));
  } else {
    cres->Print(TString::Format("%s/residuals_errors.%s",pdir.Data(),PLOT_TYPE));
  }

  delete line;
  delete hres;
  delete cres;
}

void PlotSNRnet(int nIFO, TString pdir, int sim) {

  int size = sim ? iSNR[0].size() : vLIKELIHOOD.size();
  if(size==0) return;

  double hmin=1e20;
  double hmax=0;
  for(int i=0;i<size;i++)  {
    double osnr=0;
    if(sim) for(int n=0;n<nIFO;n++) osnr+=oSNR[n][i];
    else    osnr+=vLIKELIHOOD[i];
    osnr=sqrt(osnr);
    if(osnr>hmax) hmax=osnr;
    if(osnr<hmin) hmin=osnr;
  }
  hmin-=0.4*(hmax-hmin);
  hmax+=0.4*(hmax-hmin);

  double isnr=0;
  if(sim) {
    for(int n=0;n<nIFO;n++) isnr+=iSNR[n][0];
    isnr=sqrt(isnr);
    if(isnr<hmin) hmin=0.8*isnr;
    if(isnr>hmax) hmax=1.2*isnr;
  }

  // plot SNRnet
  TCanvas* csnr = new TCanvas("SNRnet", "SNRnet", 200, 20, 600, 600);
  TH1F* hsnr = new TH1F("SNRnet","SNRnet",10,hmin,hmax);
  for(int i=0;i<size;i++)  {
    double osnr=0;
    if(sim) for(int n=0;n<nIFO;n++) osnr+=oSNR[n][i];
    else    osnr+=vLIKELIHOOD[i];
    hsnr->Fill(sqrt(osnr));
  }
  hsnr->SetLineWidth(2);

  hsnr->Draw();

  hsnr->GetXaxis()->SetTitle("SNRnet");
  hsnr->GetXaxis()->SetTitleOffset(1.4);
  hsnr->GetXaxis()->SetLabelOffset(0.02);
  hsnr->GetXaxis()->SetNoExponent();

  hsnr->GetYaxis()->SetTitle("counts");
  hsnr->GetYaxis()->SetTitleOffset(1.6);
  hsnr->GetYaxis()->SetLabelOffset(0.02);
  hsnr->GetYaxis()->SetNoExponent();

  csnr->SetLogy();
  csnr->SetGridx();
  csnr->SetGridy();

  // draw vertical line at the jres value
  float ymax = hsnr->GetMaximum();
  TLine *line = new TLine(isnr,0,isnr,ymax);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  if(sim) line->Draw();

  hsnr->Fit("gaus","q0");
  TF1* gauss = (TF1*)hsnr->GetFunction("gaus");
  if(gauss) {
    gauss->Draw("SAME");
    gauss->SetLineColor(kGreen);
    gauss->SetLineWidth(2);
  }

  char title[256];
  if(sim) sprintf(title,"SNRnet Distribution (entries=%d - SNRnet inj = %2.2g)",size,isnr);
  else    sprintf(title,"SNRnet Distribution (entries=%d)",size);
  hsnr->SetTitle(title);
  hsnr->SetStats(kTRUE);
  csnr->Print(TString::Format("%s/SNRnet.%s",pdir.Data(),PLOT_TYPE));
  hsnr->SetStats(kFALSE);

  delete line;
  delete hsnr;
  delete csnr;
}

void PlotChirpMass(int gtype, TString pdir, int sim) {

  // gtype=1 -> reconstructed chirp mass
  // gtype=2 -> reconstructed chirp mass error
  // gtype=3 -> ellipticity parameter
  // gtype=4 -> pixel fraction
  // gtype=5 -> energy fraction

  if(gtype<1 || gtype>5) return;

  TString gname = "";
  if(gtype==1) gname="chirp_mass";
  if(gtype==2) gname="chirp_mass_err";
  if(gtype==3) gname="chirp_ell";
  if(gtype==4) gname="chirp_pfrac";
  if(gtype==5) gname="chirp_efrac";

  int size = vCHIRP[gtype].size();
  if(size==0) return;

  double hmin=1e20;
  double hmax=0;
  for(int i=0;i<size;i++)  {
    double omchirp = vCHIRP[gtype][i];	// reconstructed mchirp
    if(omchirp>hmax) hmax=omchirp;
    if(omchirp<hmin) hmin=omchirp;
  }
  hmin-=0.4*(hmax-hmin);
  hmax+=0.4*(hmax-hmin);

  double imchirp=vCHIRP[0][0];		// injected mchirp
  if(sim && gtype==1) {
    if(imchirp<hmin) hmin=0.8*imchirp;
    if(imchirp>hmax) hmax=1.2*imchirp;
  }

  // plot mchirp
  TCanvas* cmchirp = new TCanvas("mchirp", "mchirp", 200, 20, 600, 600);
  TH1F* hmchirp = new TH1F("mchirp","mchirp",10,hmin,hmax);
  for(int i=0;i<size;i++) hmchirp->Fill(vCHIRP[gtype][i]);
  hmchirp->SetLineWidth(2);

  hmchirp->Draw();

  hmchirp->GetXaxis()->SetTitle(gname);
  hmchirp->GetXaxis()->SetTitleOffset(1.4);
  hmchirp->GetXaxis()->SetLabelOffset(0.02);
  hmchirp->GetXaxis()->SetNoExponent();

  hmchirp->GetYaxis()->SetTitle("counts");
  hmchirp->GetYaxis()->SetTitleOffset(1.6);
  hmchirp->GetYaxis()->SetLabelOffset(0.02);
  hmchirp->GetYaxis()->SetNoExponent();

  cmchirp->SetLogy();
  cmchirp->SetGridx();
  cmchirp->SetGridy();

  // draw vertical line at the jres value
  float ymax = hmchirp->GetMaximum();
  TLine *line = new TLine(imchirp,0,imchirp,ymax);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  if(sim && gtype==1) line->Draw();
/*
  hmchirp->Fit("gaus","q0");
  TF1* gaus = (TF1*)hmchirp->GetFunction("gaus");
  if(gaus) {
    gaus->Draw("SAME");
    gaus->SetLineColor(kGreen);
    gaus->SetLineWidth(2);
  }
*/
  gStyle->SetLineColor(kBlack);
  char title[256];
  if(sim && gtype==1) sprintf(title,"%s distribution (entries=%d - inj chirp mass = %2.2g)",gname.Data(),size,imchirp);
  else                sprintf(title,"%s distribution (entries=%d)",gname.Data(),size);
  hmchirp->SetTitle(title);
  hmchirp->SetStats(kTRUE);
  cmchirp->Print(TString::Format("%s/%s.%s",pdir.Data(),gname.Data(),PLOT_TYPE));
  hmchirp->SetStats(kFALSE);

  delete line;
  delete hmchirp;
  delete cmchirp;
}

void PlotFactors(int gtype, int nIFO, TString pdir) {

  // if gtype=0 -> Fitting Factor
  // if gtype=1 -> Overlap Factor

  if(gtype!=0 && gtype!=1) return;

  int size = iSNR[0].size();
  if(size==0) return;

  double hmin=1e20;
  double hmax=0;
  for(int i=0;i<size;i++)  {
    double  isnr=0; for(int n=0;n<nIFO;n++)  isnr+= iSNR[n][i]; 
    double  osnr=0; for(int n=0;n<nIFO;n++)  osnr+= oSNR[n][i]; 
    double iosnr=0; for(int n=0;n<nIFO;n++) iosnr+=ioSNR[n][i]; 
    double ff = iosnr/sqrt(isnr*isnr);		// fitting factor
    double of = iosnr/sqrt(isnr*osnr);		// overlap factor
    double factor = gtype==0 ? ff : of; 
    if(factor>hmax) hmax=factor; if(factor<hmin) hmin=factor;
  }
  hmin-=0.4*(hmax-hmin);
  hmax+=0.4*(hmax-hmin);

  // plot Factors
  TCanvas* cf = new TCanvas("Factors", "Factors", 200, 20, 600, 600);
  TH1F* hf = new TH1F("Factors","Factors",10,hmin,hmax);
  for(int i=0;i<size;i++)  {
    double  isnr=0; for(int n=0;n<nIFO;n++)  isnr+= iSNR[n][i]; 
    double  osnr=0; for(int n=0;n<nIFO;n++)  osnr+= oSNR[n][i]; 
    double iosnr=0; for(int n=0;n<nIFO;n++) iosnr+=ioSNR[n][i]; 
    double ff = iosnr/sqrt(isnr*isnr);		// fitting factor
    double of = iosnr/sqrt(isnr*osnr);		// overlap factor
    double factor = gtype==0 ? ff : of; 
    hf->Fill(factor);
  }
  hf->SetLineWidth(2);

  hf->Draw();

  if(gtype==0) hf->GetXaxis()->SetTitle("Fitting Factor");
  if(gtype==1) hf->GetXaxis()->SetTitle("Overlap Factor");
  hf->GetXaxis()->SetTitleOffset(1.4);
  hf->GetXaxis()->SetLabelOffset(0.02);
  hf->GetXaxis()->SetNoExponent();

  hf->GetYaxis()->SetTitle("counts");
  hf->GetYaxis()->SetTitleOffset(1.6);
  hf->GetYaxis()->SetLabelOffset(0.02);
  hf->GetYaxis()->SetNoExponent();

  cf->SetLogy();
  cf->SetGridx();
  cf->SetGridy();

  char title[256];
  if(gtype==0) sprintf(title,"Fitting Factor Distribution (entries=%d)",size);
  if(gtype==1) sprintf(title,"Overlap Factor Distribution (entries=%d)",size);
  hf->SetTitle(title);
  gStyle->SetLineColor(kBlack);
  hf->SetStats(kTRUE);
  if(gtype==0) cf->Print(TString::Format("%s/FittingFactor.%s",pdir.Data(),PLOT_TYPE));
  if(gtype==1) cf->Print(TString::Format("%s/OverlapFactor.%s",pdir.Data(),PLOT_TYPE));
  hf->SetStats(kFALSE);

  delete hf;
  delete cf;
}

void GetNullStatistic(std::vector<double>* vNUL, std::vector<double>* vNSE, int ifoId, TString ifoName, TString gtype, TString pdir) {

  if(vNUL->size()==0 || vNSE->size()==0) return;

  double xmin=1e20;
  double xmax=0;
  for(int i=0;i<vNUL->size();i++)  {
    if((*vNUL)[i]>xmax) xmax=(*vNUL)[i]; if((*vNUL)[i]<xmin) xmin=(*vNUL)[i];
  }
  xmin *= xmin>0 ? 0.5 : 1.5;
  xmax *= xmax>0 ? 1.5 : 0.5;
  if(fabs(xmin)>fabs(xmax)) xmax=fabs(xmin)*xmax/fabs(xmax); else xmin=fabs(xmax)*xmin/fabs(xmin);

  TH1F* hnul = new TH1F("null","null",10,xmin,xmax);
  for(int i=0;i<vNUL->size();i++) hnul->Fill((*vNUL)[i]);
  hnul->SetLineWidth(2);
  hnul->SetLineColor(kRed);

  TH1F* hnse = new TH1F("noise","noise",10,xmin,xmax);
  for(int i=0;i<vNSE->size();i++) hnse->Fill((*vNSE)[i]);
  double scale = double(hnul->GetEntries())/hnse->GetEntries();
  hnse->Scale(scale);
  hnse->SetLineWidth(2);
  hnse->SetLineColor(kBlack);

  double ymin=0.9;
  double ymax=0;
  for(int i=0;i<=hnul->GetNbinsX();i++) {
    double binc = hnul->GetBinContent(i);
    if(binc>ymax) ymax=binc; 
  }
  for(int i=0;i<=hnse->GetNbinsX();i++) {
    double binc = hnse->GetBinContent(i);
    if(binc>ymax) ymax=binc; 
  }
  ymax *= 1.1;
  hnse->GetYaxis()->SetRangeUser(ymin,ymax);

  double KS = hnul->Integral() ? hnul->KolmogorovTest(hnse,"N") : 0.;
  double AD = hnul->Integral() ? hnul->AndersonDarlingTest(hnse) : 0;

  hnul->Fit("gaus","q0");
  TF1* gaus = (TF1*)hnul->GetFunction("gaus");
  double chi2 = gaus ? gaus->GetChisquare()/gaus->GetNDF() : 0; 

  if(pdir!="") {
    // plot null
    TCanvas* cnul = new TCanvas("cnull", "cnull", 200, 20, 600, 600);

    hnse->Draw();
    hnul->Draw("same");

    hnse->GetXaxis()->SetTitle(gtype);
    hnse->GetXaxis()->SetTitleOffset(1.4);
    hnse->GetXaxis()->SetLabelOffset(0.02);
    hnse->GetXaxis()->SetNoExponent();

    hnse->GetYaxis()->SetTitle("counts");
    hnse->GetYaxis()->SetTitleOffset(1.6);
    hnse->GetYaxis()->SetLabelOffset(0.02);
    hnse->GetYaxis()->SetNoExponent();

    cnul->SetLogy();
    cnul->SetGridx();
    cnul->SetGridy();

    hnse->Fit("gaus","q0");
    gaus = (TF1*)hnse->GetFunction("gaus");
    if(gaus) {
      gaus->Draw("SAME");
      gaus->SetLineColor(kGreen);
      gaus->SetLineWidth(2);
    }

    char title[256];
    sprintf(title,"%s : %s (entries=%d) : KS=%0.2f : AD=%0.2f : CHI2=%0.2f",ifoName.Data(),gtype.Data(),vNUL->size(),KS,AD,chi2);
    hnse->SetTitle(title);
    hnse->SetStats(kTRUE);
    cnul->Print(TString::Format("%s/%s_%s_dist.%s",pdir.Data(),ifoName.Data(),gtype.Data(),PLOT_TYPE));
    hnse->SetStats(kFALSE);

    delete hnul;
    delete hnse;
    delete cnul;

  } else { 		// fill gPE.nstat

    int I = 2*7*ifoId;	// offset
    if(gtype=="null_90") I+=7;

    gPE.nstat[I+0] = hnul->GetMean(); 
    gPE.nstat[I+1] = hnul->GetRMS(); 
    gPE.nstat[I+2] = hnse->GetMean(); 
    gPE.nstat[I+3] = hnse->GetRMS(); 
    gPE.nstat[I+4] = chi2; 
    gPE.nstat[I+5] = KS; 
    gPE.nstat[I+6] = AD; 

    cout << endl;
    cout << " Pixels Statistic    : " << ifoName << " " << gtype << endl;
    cout << " Null  Pixels Mean   : " << gPE.nstat[I+0] << endl;
    cout << " Null  Pixels RMS    : " << gPE.nstat[I+1] << endl;
    cout << " Noise Pixels Mean   : " << gPE.nstat[I+2] << endl;
    cout << " Noise Pixels RMS    : " << gPE.nstat[I+3] << endl;
    cout << " Noise Pixels Chi2   : " << gPE.nstat[I+4] << endl;
    cout << " KolmogorovTest      : " << gPE.nstat[I+5] << endl;
    cout << " AndersonDarlingTest : " << gPE.nstat[I+6] << endl;
    cout << endl;

    delete hnul;
    delete hnse;
  }
}

double GetCentralTime(wavearray<double> x) {

  double a;
  double E=0.,T=0.;
  int size=(int)x.size();
  double rate=x.rate();
  for(int j=0;j<size;j++) {
    a = x[j];
    T += a*a*j/rate;                   // central time
    E += a*a;                          // energy
  }
  T = E>0 ? T/E : 0.5*size/rate;

  return T;
}

double GetTimeBoundaries(wavearray<double> x, double P, double& bT, double& eT) {

  if(P<0) P=0;
  if(P>1) P=1;

  int N = x.size();

  double E = 0;                                                 // signal energy
  double avr = 0;                                               // average
  for(int i=0;i<N;i++) {avr+=i*x[i]*x[i]; E+=x[i]*x[i];}
  int M=int(avr/E);                                             // central index

  // search range which contains percentage P of the total energy E
  int jB=0;
  int jE=N-1;
  double a,b;
  double sum = ((M>=0)&&(M<N)) ? x[M]*x[M] : 0.;
  for(int j=1; j<N; j++) {
    a = ((M-j>=0)&&(M-j<N)) ? x[M-j] : 0.;
    b = ((M+j>=0)&&(M+j<N)) ? x[M+j] : 0.;
    if(a) jB=M-j;
    if(b) jE=M+j;
    sum += a*a+b*b;
    if(sum/E > P) break;
  }

  bT = x.start()+jB/x.rate();
  eT = x.start()+jE/x.rate();

  return eT-bT;
}

double GetCentralFrequency(wavearray<double> x) {

  double a;
  double E=0.,F=0.;
  int size=(int)x.size();
  double rate=x.rate();
  x.FFTW(1);
  double dF=(rate/(double)size)/2.;
  for(int j=0;j<size;j+=2) {
    a = x[j]*x[j]+x[j+1]*x[j+1];
    F += a*j*dF;                       // central frequency
    E += a;                            // energy
  }
  F = E>0 ? F/E : 0.5*rate;

  return F;
}

double GetFrequencyBoundaries(wavearray<double> x, double P, double& bF, double& eF) {

  if(P<0) P=0;
  if(P>1) P=1;

  int N = x.size();

  x.FFTW(1);

  double E = 0;                                                 // signal energy
  double avr = 0;                                               // average
  for(int i=0;i<N;i++) {avr+=i*x[i]*x[i]; E+=x[i]*x[i];}
  int M=int(avr/E);                                             // central index

  // search range which contains percentage P of the total energy E
  int jB=0;
  int jE=N-1;
  double a,b;
  double sum = ((M>=0)&&(M<N)) ? x[M]*x[M] : 0.;
  for(int j=1; j<N; j++) {
    a = ((M-j>=0)&&(M-j<N)) ? x[M-j] : 0.;
    b = ((M+j>=0)&&(M+j<N)) ? x[M+j] : 0.;
    if(a) jB=M-j;
    if(b) jE=M+j;
    sum += a*a+b*b;
    if(sum/E > P) break;
  }

  double dF=(x.rate()/(double)x.size())/2.;

  bF = jB*dF;
  eF = jE*dF;

  return eF-bF;
}

wavearray<double> GetCutWaveform(wavearray<double> x, double bT, double eT, double bF, double eF) {

  bT-=x.start();
  eT-=x.start();

  int size=(int)x.size();

  // cut time range bT,eT
  double T=0.;
  double dT = 1./x.rate();
  for(int j=0;j<size;j++) {
    T = j*dT; 
    if(T<bT || T>eT) x[j]=0;
  }

  // cut frequency range bF,eF
  double F=0.;
  double dF=(x.rate()/(double)x.size())/2.;
  x.FFTW(1);
  for(int j=0;j<size;j+=2) {
    F = j*dF; 
    if(F<bF || F>eF) {x[j]=0;x[j+1]=0;}
  }
  x.FFTW(-1);

  return x;
}

void SetWaveformCuts(wavearray<double>* x, double bT, double eT, double bF, double eF) {

  bT-=x->start();
  eT-=x->start();

  int size=(int)x->size();

  // cut time range bT,eT
  double T=0.;
  double dT = 1./x->rate();
  for(int j=0;j<size;j++) {
    T = j*dT; 
    if(T<bT || T>eT) x->data[j]=0;
  }

  // cut frequency range bF,eF
  double F=0.;
  double dF=(x->rate()/(double)x->size())/2.;
  x->FFTW(1);
  for(int j=0;j<size;j+=2) {
    F = j*dF; 
    if(F<bF || F>eF) {x->data[j]=0;x->data[j+1]=0;}
  }
  x->FFTW(-1);
}

TString DumpCED(network* NET, netevent* &EVT, CWB::config* cfg, double ofactor) {

  cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
  int gIFACTOR=-1; IMPORT(int,gIFACTOR)

  int nIFO = NET->ifoListSize();  	// number of detectors
  double factor = cfg->factors[gIFACTOR];

  char sfactor[32];
  if(cfg->simulation==3) {
    if(factor<0)  sprintf(sfactor,"n%g",fabs(factor));
    if(factor==0) sprintf(sfactor,"z%g",factor);
    if(factor>0)  sprintf(sfactor,"p%g",factor);
  } else          sprintf(sfactor,"%g",factor);

  char out_CED[1024];
  if(cfg->simulation) {
    char sim_label[512];
    sprintf(sim_label,"%d_%d_%s_%s_job%d",int(gCWB2G->Tb),int(gCWB2G->dT),cfg->data_label,sfactor,gCWB2G->runID);
    sprintf(out_CED,"%s/ced_%s_%d",cfg->nodedir,sim_label,gSystem->GetPid());
  } else {
    char prod_label[512];
    sprintf(prod_label,"%d_%d_%s_slag%d_lag%lu_%lu_job%d",
            int(gCWB2G->Tb),int(gCWB2G->dT),cfg->data_label,gCWB2G->slagID,cfg->lagOff,cfg->lagSize,gCWB2G->runID);
    sprintf(out_CED,"%s/ced_%s_%d",cfg->nodedir,prod_label,gSystem->GetPid());
  }
  cout << out_CED << endl;

  // restore whitened data into the detector TF map
  if(gOPT.noise==0) {	
    for(int n=0; n<nIFO; n++) {
      WSeries<double>* pTF = NET->getifo(n)->getTFmap();
      *pTF = gHOT[n];
    }
  }

  cout<<"DumpCED : dump ced to disk ..." <<endl;
  CWB::ced ced(NET, EVT, out_CED);
  cfg->cedRHO=0;
  switch(gCWB2G->istage) {
  case CWB_STAGE_FULL:
  case CWB_STAGE_INIT:
    // use TF map & inj signals
    ced.SetOptions(cfg->simulation,cfg->cedRHO,cfg->inRate);
    for(int n=0; n<nIFO; n++) {gCWB2G->pTF[n]->setlow(cfg->fLow); gCWB2G->pTF[n]->sethigh(cfg->fHigh);}
    break;
  default:
    // use sparse map & inj signals
    ced.SetOptions(cfg->simulation,cfg->cedRHO,cfg->inRate,true);
  }
  if(gCWB2G->singleDetector) ced.SetChannelName(cfg->channelNamesRaw[0]);
  bool fullCED = gCWB2G->singleDetector ? false : true;
  ced.Write(ofactor,fullCED);

  char ifostr[20]="";
  if(gCWB2G->singleDetector) {
    sprintf(ifostr,"%s%s",ifostr,NET->ifoName[0]);
  } else {
    for(int n=0; n<nIFO; n++) sprintf(ifostr,"%s%s",ifostr,NET->ifoName[n]);
  }

  char dirCED[1024];
  sprintf(dirCED, "%s/%s", out_CED, ifostr);
  if(gCWB2G->singleDetector) {
    sprintf(dirCED, "%s_%.3f",dirCED,EVT->start[0]);
  } else {
    for(int n=0; n<nIFO; n++) sprintf(dirCED, "%s_%.3f",dirCED,EVT->start[n]);
  }

  return dirCED;
}

void CreateIndexHTML(TString dirCED, int nIFO, TString* ifo, bool sim) {

  cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
  if(gCWB2G->singleDetector) nIFO=1;

  // ------------------------------------------------------------------------------------
  // expand html file paths                                                              
  // ------------------------------------------------------------------------------------

  TString pe_index = gSystem->ExpandPathName(PE_INDEX);

  // ------------------------------------------------------------------------------------
  // check if environmental variables are defined                                        
  // ------------------------------------------------------------------------------------
  if(gSystem->Getenv("HOME_CED_WWW")==NULL) {                                            
    cout << "Error : environment HOME_CED_WWW is not defined!!! " << endl;exit(1);}      

  char ofileName[256]="";
  sprintf(ofileName,"%s/pe_index.html",dirCED.Data());                                   
  cout << "ofileName : " <<  ofileName << endl;                             

  // -------------------------------------------------------------------------------
  // open output file
  // -------------------------------------------------------------------------------
  ofstream out;
  out.open(ofileName,ios::out);
  if(!out.good()) {cout << "Error Opening File : " << ofileName << endl;exit(1);}
  // HEADER
  WriteBodyHTML(pe_index, "PE_HEADER_BEG", "PE_HEADER_END", &out);
  // BODY 
  WriteBodyHTML(pe_index, "PE_BODY_BEG", "PE_BODY_END", &out);
  // TABBERS                                  
  int sbody_height = 0;
  out << "<div class=\"tabber\">" << endl;
  // -------------------------------------------------------------------------------
  // Time-Frequency Maps
  // -------------------------------------------------------------------------------
  if(gOPT.ced_tfmap) { 
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>Time-Frequency Maps</h2>" << endl;
    sbody_height = 1000+nIFO*400;
    if(gOPT.ced_pca) sbody_height += 600;
    out << "  <iframe src=\"tfmap_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;                                                          
  }
  // -------------------------------------------------------------------------------
  // Reconstructed Detector Responses
  // -------------------------------------------------------------------------------
  if(gOPT.ced_rdr) { 
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>Reconstructed Detector Responses</h2>" << endl;
    if(!sim) sbody_height = 350+nIFO*450;                        
    else     sbody_height = 350+nIFO*1340;                       
    out << "  <iframe src=\"rec_signal_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;                                                          
  }
  // -------------------------------------------------------------------------------
  // Power Spectral Density
  // -------------------------------------------------------------------------------
  if(gOPT.ced_psd) { 
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>PSD</h2>" << endl;
    sbody_height = 350+(nIFO+1)*750;                       
    out << "  <iframe src=\"psd_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;                                                          
  }
  // -------------------------------------------------------------------------------
  // SkyMaps
  // -------------------------------------------------------------------------------
  sbody_height = 2050;
  if(gOPT.ced_skymap && !gCWB2G->singleDetector) {
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>SkyMaps</h2>" << endl;       
    out << "  <iframe src=\"skymap_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;    
  }                                                      
  // -------------------------------------------------------------------------------
  // Reconstructed vs Data/Average
  // -------------------------------------------------------------------------------
  if(gOPT.ced_rec) { 
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>Waveform Errors</h2>" << endl;
    sbody_height = 370+nIFO*880;                        
    out << "  <iframe src=\"rec_avr_signal_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;                                                          
  }
  // -------------------------------------------------------------------------------
  // Reconstructed vs Injection
  // -------------------------------------------------------------------------------
  if(sim && gOPT.ced_inj) { 
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>Injection</h2>" << endl;
    sbody_height = 390+nIFO*880;
    out << "  <iframe src=\"rec_inj_signal_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;                                                          
  }
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // Reconstructed vs Reduced Injection
  // -------------------------------------------------------------------------------
  if(sim && gOPT.ced_rinj) { 
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>Injection</h2>" << endl;
    sbody_height = 390+nIFO*320;
    out << "  <iframe src=\"rec_inj_signal_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;                                                          
  }
  // -------------------------------------------------------------------------------
  // Chirp Mass Distributions
  // -------------------------------------------------------------------------------
  if(gOPT.ced_cm) { 
    sbody_height = 2100;
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>Chirp Mass</h2>" << endl;       
    out << "  <iframe src=\"mchirp_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;                                                          
  }
  // -------------------------------------------------------------------------------
  // Distributions
  // -------------------------------------------------------------------------------
  if(gOPT.ced_distr) { 
    sbody_height = 2000;
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>Distributions</h2>" << endl;       
    out << "  <iframe src=\"distribution_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;     
  }                                                     
  // -------------------------------------------------------------------------------
  // Null Pixels Distribution
  // -------------------------------------------------------------------------------
  if(gOPT.ced_null) { 
    sbody_height = 100+nIFO*1300;
    out << "<div class=\"tabbertab\">" << endl;
    out << "  <h2>Null Pixels</h2>" << endl;       
    out << "  <iframe src=\"nstat_body.html\" width=\"100%\" "
        << " height=\"" << sbody_height << "px\" frameborder=\"0\"></iframe>" << endl;
    out << "</div>" << endl;   
  }                                                       
  // -------------------------------------------------------------------------------
  // close output file
  // -------------------------------------------------------------------------------
  out << "</div>" << endl;
  out.close();


  if(gOPT.ced_tfmap) { 
    out.open(TString::Format("%s/tfmap_body.html",dirCED.Data()).Data(),ios::out);
    // TFMAPS
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    WriteBodyHTML(pe_index, "PE_TFMAP_HEADER_BEG", "PE_TFMAP_HEADER_END", &out);
    WriteBodyHTML(pe_index, "PE_TFMAP_BEG", "PE_TFMAP_END", &out, nIFO, ifo);
    // LIKELIHOOD
    out << "<p><br /></p> " << endl;
    out << "<p><br /></p> " << endl;
    WriteBodyHTML(pe_index, "PE_LIKELIHOOD_BEG", "PE_LIKELIHOOD_END", &out);
    // POLARGRAMS	: NOT IMPLEMENTED IN WP !!!
    //out << "<p><br /></p> " << endl;
    //out << "<p><br /></p> " << endl;
    //WriteBodyHTML(pe_index, "PE_POLARGRAM_BEG", "PE_POLARGRAM_END", &out);
    if(gOPT.ced_pca) WriteBodyHTML(pe_index, "PE_PCA_BEG", "PE_PCA_END", &out);
    out.close();                                                                          
  }

  // RECONSTRUCTED SIGNAL
  if(gOPT.ced_rdr) { 
    out.open(TString::Format("%s/rec_signal_body.html",dirCED.Data()).Data(),ios::out);
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    WriteBodyHTML(pe_index, "PE_REC_SIGNAL_HEADER_BEG", "PE_REC_SIGNAL_HEADER_END", &out);
    for(int n=0;n<nIFO;n++) {                                                                                                                     
      WriteBodyHTML(pe_index, "PE_IFO_SIGNAL_BEG", "PE_IFO_SIGNAL_END", &out, 1, &ifo[n]);
      if(sim) {                                                                                                                          
        WriteBodyHTML(pe_index, "PE_INJ_SIGNAL_BEG", "PE_INJ_SIGNAL_END", &out, 1, &ifo[n]);
      }
      WriteBodyHTML(pe_index, "PE_REC_SIGNAL_BEG", "PE_REC_SIGNAL_END", &out, 1, &ifo[n]);
      out << "<p><br /></p> " << endl;
      out << "<p><br /></p> " << endl;
    }
    out.close();
  }

  // PSD
  if(gOPT.ced_psd) { 
    out.open(TString::Format("%s/psd_body.html",dirCED.Data()).Data(),ios::out);
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    WriteBodyHTML(pe_index, "PE_PSD_HEADER_BEG", "PE_PSD_HEADER_END", &out);
    for(int n=0;n<nIFO;n++) {                                                                                                                     
      WriteBodyHTML(pe_index, "PE_IFO_PSD_BEG", "PE_IFO_PSD_END", &out, 1, &ifo[n]);
      WriteBodyHTML(pe_index, "PE_PSD_BEG", "PE_PSD_END", &out, 1, &ifo[n]);
      out << "<p><br /></p> " << endl;
      out << "<p><br /></p> " << endl;
    }
    out.close();
  }

  // SKYMAP
  if(gOPT.ced_skymap && !gCWB2G->singleDetector) {
    out.open(TString::Format("%s/skymap_body.html",dirCED.Data()).Data(),ios::out);
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    WriteBodyHTML(pe_index, "PE_SKYMAP_BEG", "PE_SKYMAP_END", &out);
    out.close();
  }

  // RECONSTRUCTED vs DATA/AVERAGE
  if(gOPT.ced_rec) { 
    out.open(TString::Format("%s/rec_avr_signal_body.html",dirCED.Data()).Data(),ios::out);
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    WriteBodyHTML(pe_index, "PE_REC_AVR_SIGNAL_HEADER_BEG", "PE_REC_AVR_SIGNAL_HEADER_END", &out);
    for(int n=0;n<nIFO;n++) {                                                                                                                     
      WriteBodyHTML(pe_index, "PE_REC_AVR_SIGNAL_BEG", "PE_REC_AVR_SIGNAL_END", &out, 1, &ifo[n]);
      out << "<p><br /></p> " << endl;
      out << "<p><br /></p> " << endl;
    }
    out.close();
  }

  // RECONSTRUCTED vs INJECTED
  if(sim && (gOPT.ced_inj || gOPT.ced_rinj)) {                                                                                                                          
    out.open(TString::Format("%s/rec_inj_signal_body.html",dirCED.Data()).Data(),ios::out);
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    if(gOPT.ced_inj)  WriteBodyHTML(pe_index, "PE_REC_INJ_SIGNAL_0_HEADER_BEG", "PE_REC_INJ_SIGNAL_0_HEADER_END", &out);
    if(gOPT.ced_rinj) WriteBodyHTML(pe_index, "PE_REC_INJ_SIGNAL_HEADER_BEG", "PE_REC_INJ_SIGNAL_HEADER_END", &out);
    if(sim) { 
      for(int n=0;n<nIFO;n++) {                                                                                                                     
        if(gOPT.ced_inj)  WriteBodyHTML(pe_index, "PE_REC_INJ_SIGNAL_0_BEG", "PE_REC_INJ_SIGNAL_0_END", &out, 1, &ifo[n]);
        if(gOPT.ced_rinj) WriteBodyHTML(pe_index, "PE_REC_INJ_SIGNAL_BEG", "PE_REC_INJ_SIGNAL_END", &out, 1, &ifo[n]);
        out << "<p><br /></p> " << endl;
        out << "<p><br /></p> " << endl;
      }
    }
    out.close();
  }


  // DISTRIBUTION
  if(gOPT.ced_distr) { 
    out.open(TString::Format("%s/distribution_body.html",dirCED.Data()).Data(),ios::out);
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    WriteBodyHTML(pe_index, "PE_DISTRIBUTION_BEG", "PE_DISTRIBUTION_END", &out);
    if(sim) WriteBodyHTML(pe_index, "PE_FACTORS_DISTRIBUTION_BEG", "PE_FACTORS_DISTRIBUTION_END", &out);
    out.close();
  }

  // NULL PIXELS DISTRIBUTION
  if(gOPT.ced_null) { 
    out.open(TString::Format("%s/nstat_body.html",dirCED.Data()).Data(),ios::out);
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    WriteBodyHTML(pe_index, "PE_NULPIX_DISTRIBUTION_HEADER_BEG", "PE_NULPIX_DISTRIBUTION_HEADER_END", &out);
    for(int n=0;n<nIFO;n++) {                                                                                                                     
      WriteBodyHTML(pe_index, "PE_NULPIX_DISTRIBUTION_BEG", "PE_NULPIX_DISTRIBUTION_END", &out, 1, &ifo[n]);
      out << "<p><br /></p> " << endl;
      out << "<p><br /></p> " << endl;
    }
    out.close();
  }

  // MCHIRP DISTRIBUTION
  if(gOPT.ced_cm) { 
    out.open(TString::Format("%s/mchirp_body.html",dirCED.Data()).Data(),ios::out);
    WriteBodyHTML(pe_index, "PE_JSCRIPT_BEG", "PE_JSCRIPT_END", &out);
    WriteBodyHTML(pe_index, "PE_MCHIRP_DISTRIBUTION_HEADER_BEG", "PE_MCHIRP_DISTRIBUTION_HEADER_END", &out);
    WriteBodyHTML(pe_index, "PE_MCHIRP_DISTRIBUTION_BEG", "PE_MCHIRP_DISTRIBUTION_END", &out);
    out.close();
  }

  // copy tabber javascript to CED directory
  char cmd[1024];
  sprintf(cmd,"cp ${HOME_WAT}//html/etc/html/tabber.css %s",dirCED.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp ${HOME_WAT}//html/etc/html/tabber.js %s",dirCED.Data());
  gSystem->Exec(cmd);
  // overwrite CED index.html
  sprintf(cmd,"mv %s/pe_index.html %s/index.html",dirCED.Data(),dirCED.Data());
  gSystem->Exec(cmd);

  return;
}

void WriteBodyHTML(TString html_template, TString html_tag_beg, TString html_tag_end, ofstream* out, int nIFO, TString* ifo) {

  cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
  if(gCWB2G->singleDetector) nIFO=1;

  // get SITE_CLUSTER
  TString site_cluster="VIRTUALBOX";                                            // default value
  if(gSystem->Getenv("SITE_CLUSTER")!=NULL) {
    site_cluster=TString(gSystem->Getenv("SITE_CLUSTER"));
  }
  TString cluster_site_logo  = "cluster_site_logo_modern.png";
  TString cluster_site_url1  = "http://www.ligo.caltech.edu/";
  TString cluster_site_name1 = "LIGO Homepage";
  TString cluster_site_url2  = "https://www.virgo-gw.eu/";
  TString cluster_site_name2 = "VIRGO Homepage";
  if(site_cluster=="ATLAS") {
    cluster_site_logo  = "atlas_logo_modern.png";
    cluster_site_url1  = "http://www.aei.mpg.de/14026/AEI_Hannover/";
    cluster_site_name1 = "AEI Hannover Homepage";
    cluster_site_url2  = "http://www.aei.mpg.de/14026/AEI_Hannover";
    cluster_site_name2 = "AEI Hannover Homepage";
  }
  if(site_cluster=="CIT") {
    cluster_site_logo  = "ligo_virgo_logo_modern.png";
    cluster_site_url1  = "http://www.ligo.caltech.edu/";
    cluster_site_name1 = "LIGO Homepage";
    cluster_site_url2  = "https://www.virgo-gw.eu/";
    cluster_site_name2 = "VIRGO Homepage";
  }
  if(site_cluster=="VIRTUALBOX") {
    cluster_site_logo  = "cluster_site_logo_modern.png";
    cluster_site_url1  = "https://www.gw-openscience.org/about/";
    cluster_site_name1 = "GWOSC Homepage";
    cluster_site_url2  = "https://www.gw-openscience.org/about/";
    cluster_site_name2 = "GWOSC Homepage";
  }

  // get HOME_WWW
  TString home_www="~waveburst/waveburst";                                      // default value
  if(gSystem->Getenv("HOME_WWW")!=NULL) {
    home_www=TString(gSystem->Getenv("HOME_WWW"));
  }

  // get HOME_CED_WWW
  TString home_ced_www="~waveburst/waveburst/ced-1.0-modern";                   // default value
  if(gSystem->Getenv("HOME_CED_WWW")!=NULL) {
    home_ced_www=TString(gSystem->Getenv("HOME_CED_WWW"));
  }

  // get CWB_DOC_URL
  TString cwb_doc_url="https://ldas-jobs.ligo.caltech.edu/~waveburst/doc";      // default value
  if(gSystem->Getenv("CWB_DOC_URL")!=NULL) {
    cwb_doc_url=TString(gSystem->Getenv("CWB_DOC_URL"));
  }

  // get CWB_GIT_URL
  TString cwb_git_url="https://git.ligo.org";                                   // default value
  if(gSystem->Getenv("CWB_GIT_URL")!=NULL) {
    cwb_git_url=TString(gSystem->Getenv("CWB_GIT_URL"));
  }

  // get issues URL (Ex: https://gitlab.com/groups/gwburst/-/issues)
  TString cwb_git_issues_url=cwb_git_url;
  TString partA=cwb_git_issues_url;
  TString partB=cwb_git_issues_url;
  partA.Remove(partA.Last('/'),partA.Sizeof()-1);
  partB.Remove(0,partB.Last('/')+1);
  cwb_git_issues_url=partA+"/groups/"+partB+"/-/issues";


  char line[1024];
  bool output=false;
  for(int n=0;n<nIFO;n++) {                                                                                                                     
    ifstream in;
    in.open(html_template.Data(),ios::in);
    if(!in.good()) {cout << "Error Opening File : " << html_template.Data() << endl;exit(1);}
    while(1) {
      in.getline(line,1024);
      if(!in.good()) break;
      TString sline = line;
      if(sline.Contains(html_tag_beg)&&!output) output=true;
      if(sline.Contains(html_tag_end)&& output) {output=false;break;}
      if(ifo!=NULL) sline.ReplaceAll("IFO",ifo[n]);
      sline.ReplaceAll("xNTRIALS",TString::Format("%d",(int)vCHIRP[0].size()-1).Data());
      sline.ReplaceAll("xINRATE",TString::Format("%d",(int)gINRATE));
      sline.ReplaceAll("xRATEANA",TString::Format("%d",(int)gRATEANA));
#ifdef PLOT_MEDIAN
      sline.ReplaceAll("xMEAN","Median");
      sline.ReplaceAll("xSIGMA","50% Percentile");
      sline.ReplaceAll("x2SIGMA","90% Percentile");
#else
      sline.ReplaceAll("xMEAN","Mean");
      sline.ReplaceAll("xSIGMA","RMS");
      sline.ReplaceAll("x2SIGMA","2RMS");
#endif
      sline.ReplaceAll("HOME_CED_WWW",home_ced_www.Data());
      sline.ReplaceAll("CWB_DOC_URL",cwb_doc_url.Data());                                 
      sline.ReplaceAll("HOME_WWW",home_www.Data());                                       
      if(cwb_doc_url!="") {
        sline.ReplaceAll("<!--CWB_DOC_URL","");
        sline.ReplaceAll("CWB_DOC_URL-->",""); 
        sline.ReplaceAll("XCWB_DOC_URL",cwb_doc_url.Data());
      } 
      // replace HOME_CED_WWW, HOME_CED_WWW, CWB_DOC_URL with used defined values (watenv)
      sline.ReplaceAll("HOME_WWW",home_www);
      sline.ReplaceAll("HOME_CED_WWW",home_ced_www);
      if(cwb_doc_url.Contains("gwburst.gitlab.io")) {     // public cWB documentation     
        sline.ReplaceAll("/CWB_DOC_URL/cwb/man",cwb_doc_url);
      } else {
        sline.ReplaceAll("CWB_DOC_URL",cwb_doc_url);
      }
      sline.ReplaceAll("CWB_GIT_URL",cwb_git_url);
      sline.ReplaceAll("CWB_GIT_ISSUES_URL",cwb_git_issues_url);
      sline.ReplaceAll("CLUSTER_SITE_LOGO",cluster_site_logo);
      sline.ReplaceAll("CLUSTER_SITE_URL1",cluster_site_url1);
      sline.ReplaceAll("CLUSTER_SITE_NAME1",cluster_site_name1);
      sline.ReplaceAll("CLUSTER_SITE_URL2",cluster_site_url2);
      sline.ReplaceAll("CLUSTER_SITE_NAME2",cluster_site_name2);

      sline.ReplaceAll(TString("=\"/http"), TString("=\"http"));
      sline.ReplaceAll(TString("=\"//"), TString("=\"/"));

      if(output) (*out) << sline.Data() << endl;
    }
    in.close();
  }
}

void
ClearWaveforms(detector* ifo) {

  int n;

  n = ifo->IWFP.size();
  for (int i=0;i<n;i++) {
    wavearray<double>* wf = (wavearray<double>*)ifo->IWFP[i];
    delete wf;
  }
  ifo->IWFP.clear();
  ifo->IWFID.clear();

  n = ifo->RWFP.size();
  for (int i=0;i<n;i++) {
    wavearray<double>* wf = (wavearray<double>*)ifo->RWFP[i];
    delete wf;
  }
  ifo->RWFP.clear();
  ifo->RWFID.clear();
}

void PlotFrequencyErrors(TString ofname, TString title, CWB::config* cfg, wavearray<double>* frec,
                  wavearray<double>* favr, wavearray<double>* ferr, wavearray<double>* wref, TString pdir, double P) {

  int size = frec->size();

  wavearray<double> time(size);
  wavearray<double> etime(size); etime=0;
  for (int i=0; i<size; i++) time.data[i] = i/frec->rate()+(wref->start()-gSEGGPS);

  wavearray<double> f2err=*ferr; f2err*=2.;

  TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",gSEGGPS);

  double bT, eT;
  GetTimeBoundaries(*wref, P, bT, eT);
  bT-=gSEGGPS;
  eT-=gSEGGPS;
  for(int i=0;i<frec->size();i++) {
    double xtime = i/wref->rate();
    if(xtime>bT && xtime<eT) continue;
    frec->data[i]=0;
    favr->data[i]=0;
    ferr->data[i]=0;
    f2err.data[i]=0;
  }

  TGraphErrors* e2gr = new TGraphErrors(size,time.data,favr->data,etime.data,f2err.data);
  e2gr->SetLineColor(17);
  e2gr->SetFillStyle(1001);
  e2gr->SetFillColor(17);
  e2gr->GetXaxis()->SetTitle(xtitle);
  e2gr->GetYaxis()->SetTitle("frequency (hz)");
  e2gr->SetTitle(title);
  e2gr->GetXaxis()->SetTitleFont(42);
  e2gr->GetXaxis()->SetLabelFont(42);
  e2gr->GetXaxis()->SetLabelOffset(0.012);
  e2gr->GetXaxis()->SetTitleOffset(1.5);
  e2gr->GetYaxis()->SetTitleFont(42);
  e2gr->GetYaxis()->SetLabelFont(42);
  e2gr->GetYaxis()->SetLabelOffset(0.01);
  e2gr->GetYaxis()->SetTitleOffset(1.4);


  TGraphErrors* egr = new TGraphErrors(size,time.data,favr->data,etime.data,ferr->data);
  egr->SetLineColor(15);
  egr->SetFillStyle(1001);
  egr->SetFillColor(15);

  TGraph* agr = new TGraph(size,time.data,favr->data);
  agr->SetLineWidth(1);
  agr->SetLineColor(kWhite);
  agr->SetLineStyle(1);

  TGraph* gr = new TGraph(size,time.data,frec->data);
  gr->SetLineWidth(1);
  gr->SetLineColor(2);

  TCanvas* canvas = new TCanvas("frec", "frec",200,20,800,500);
  canvas->cd();
  canvas->SetGridx();
  canvas->SetGridy();

  e2gr->GetXaxis()->SetRangeUser(bT, eT);
  e2gr->Draw("AP4");
  egr->Draw("P4same");
  agr->Draw("Lsame");
  gr->Draw("Lsame");

  if(ofname!="") {
    ofname = TString(pdir)+TString("/")+ofname;
    canvas->Print(ofname);
    cout << "write : " << ofname << endl;
  }

  delete canvas;
  delete e2gr;
  delete egr;
  delete agr;
  delete gr;
}

void PlotWaveformErrors(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wrec,
                  wavearray<double>* wavr, wavearray<double>* werr, wavearray<double>* wref, TString pdir, double P) {

  int size = wrec->size();

  wavearray<double> time(size);
  wavearray<double> etime(size); etime=0;
  for (int i=0; i<size; i++) time.data[i] = i/wrec->rate()+(wref->start()-gSEGGPS);

  wavearray<double> w2err=*werr; w2err*=2.;

  TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",gSEGGPS);

  TGraphErrors* e2gr = new TGraphErrors(size,time.data,wavr->data,etime.data,w2err.data);
  e2gr->SetLineColor(17);
  e2gr->SetFillStyle(1001);
  e2gr->SetFillColor(17);
  e2gr->GetXaxis()->SetTitle(xtitle);
  e2gr->GetYaxis()->SetTitle("magnitude");
  e2gr->SetTitle(title);
  e2gr->GetXaxis()->SetTitleFont(42);
  e2gr->GetXaxis()->SetLabelFont(42);
  e2gr->GetXaxis()->SetLabelOffset(0.012);
  e2gr->GetXaxis()->SetTitleOffset(1.5);
  e2gr->GetYaxis()->SetTitleFont(42);
  e2gr->GetYaxis()->SetLabelFont(42);
  e2gr->GetYaxis()->SetLabelOffset(0.01);
  e2gr->GetYaxis()->SetTitleOffset(1.4);

  TGraphErrors* egr = new TGraphErrors(size,time.data,wavr->data,etime.data,werr->data);
  egr->SetLineColor(15);
  egr->SetFillStyle(1001);
  egr->SetFillColor(15);

  TGraph* agr = new TGraph(size,time.data,wavr->data);
  agr->SetLineWidth(1);
  agr->SetLineColor(kWhite);
  agr->SetLineStyle(1);

  TGraph* gr = new TGraph(size,time.data,wrec->data);
  gr->SetLineWidth(1);
  gr->SetLineColor(2);

  TCanvas* canvas = new TCanvas("wrec", "wrec",200,20,800,500);
  canvas->cd();
  canvas->SetGridx();
  canvas->SetGridy();

  double bT, eT;
  GetTimeBoundaries(*wref, P, bT, eT);
  bT-=gSEGGPS;
  eT-=gSEGGPS;
  e2gr->GetXaxis()->SetRangeUser(bT, eT);
  e2gr->Draw("A3");
  egr->Draw("3same");
  agr->Draw("Lsame");
  gr->Draw("Lsame");

  if(ofname!="") {
    ofname = TString(pdir)+TString("/")+ofname;
    canvas->Print(ofname);
    cout << "write : " << ofname << endl;
  }

  delete canvas;
  delete e2gr;
  delete egr;
  delete agr;
  delete gr;
}

void PlotWaveformAsymmErrors(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wrec,
                             wavearray<double>* wmed, wavearray<double>* wl50, wavearray<double>* wu50, 
                             wavearray<double>* wl90, wavearray<double>* wu90, wavearray<double>* wref, TString pdir, double P, bool freq) {

  int size = wrec->size();

  wavearray<double> time(size);
  wavearray<double> etime(size); etime=0;
  for (int i=0; i<size; i++) time[i] = i/wrec->rate()+(wref->start()-gSEGGPS);

  double bT, eT;
  GetTimeBoundaries(*wref, P, bT, eT);
  bT-=gSEGGPS;
  eT-=gSEGGPS;

  // set to 0 the frequency values outside the time range -> fix the y scale autoscale
  // info : this procedure modify the frequency input data but it is not relevant
  if(freq) {  	
    for(int i=0;i<wrec->size();i++) {
      if(time[i]>bT && time[i]<eT) continue;
      wrec->data[i]=0; wmed->data[i]=0; wl50->data[i]=0; wu50->data[i]=0; wl90->data[i]=0; wu90->data[i]=0;
    }
  }

  TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",gSEGGPS);

  TGraphAsymmErrors* egr90 = new TGraphAsymmErrors(size,time.data,wmed->data,etime.data,etime.data,wl90->data,wu90->data);
  egr90->SetLineColor(17);
  egr90->SetFillStyle(1001);
  egr90->SetFillColor(17);
  egr90->GetXaxis()->SetTitle(xtitle);
  if(freq) egr90->GetYaxis()->SetTitle("frequency (hz)");
  else     egr90->GetYaxis()->SetTitle("magnitude");
  egr90->SetTitle(title);
  egr90->GetXaxis()->SetTitleFont(42);
  egr90->GetXaxis()->SetLabelFont(42);
  egr90->GetXaxis()->SetLabelOffset(0.012);
  egr90->GetXaxis()->SetTitleOffset(1.5);
  egr90->GetYaxis()->SetTitleFont(42);
  egr90->GetYaxis()->SetLabelFont(42);
  egr90->GetYaxis()->SetLabelOffset(0.01);
  egr90->GetYaxis()->SetTitleOffset(1.4);

  TGraphAsymmErrors* egr50 = new TGraphAsymmErrors(size,time.data,wmed->data,etime.data,etime.data,wl50->data,wu50->data);
  egr50->SetLineColor(15);
  egr50->SetFillStyle(1001);
  egr50->SetFillColor(15);

  TGraph* agr = new TGraph(size,time.data,wmed->data);
  agr->SetLineWidth(1);
  agr->SetLineColor(kWhite);
  agr->SetLineStyle(1);

  TGraph* gr = new TGraph(size,time.data,wrec->data);
  gr->SetLineWidth(1);
  gr->SetLineColor(2);

  TCanvas* canvas = new TCanvas("wrec", "wrec",200,20,800,500);
  canvas->cd();
  canvas->SetGridx();
  canvas->SetGridy();

  egr90->GetXaxis()->SetRangeUser(bT, eT);
  egr90->Draw("A3");
  egr50->Draw("3same");
  agr->Draw("Lsame");
  gr->Draw("Lsame");

  if(ofname!="") {
    ofname = TString(pdir)+TString("/")+ofname;
    canvas->Print(ofname);
    cout << "write : " << ofname << endl;
  }

  delete canvas;
  delete egr50;
  delete egr90;
  delete agr;
  delete gr;
}

// Dumps reconstructed waveform/time/freq/errors array in ASCII format.
void DumpRecWavePE(network* NET, TString pdir) {

  int nIFO = NET->ifoListSize();	// number of detectors

  char ofname[256];
  for(int n=0; n<nIFO; n++) {

    sprintf(ofname,"%s/%s_pe_wave.dat",pdir.Data(),NET->ifoName[n]);
 
    ofstream out;
    out.open(ofname,ios::out);
    if (!out.good()) {cout << "Error Opening Output File : " << ofname << endl;exit(1);}
    cout << "Create Output File : " << ofname << endl;
    out.precision(19); 

    // write header
    out << "#whitened data : time," << 
                           " amp_point, amp_mean, amp_rms," << 
                           " amp_median, amp_lower_50_perc, amp_lower_90_perc, amp_upper_50_perc, amp_upper_90_perc," << 
                           " frq_point, frq_mean, frq_rms," << 
                           " frq_median, frq_lower_50_perc, frq_lower_90_perc, frq_upper_50_perc, frq_upper_90_perc" << endl;

    // write data
    int size = vREC[n].size();
    double dt=1./vREC[n].rate();
    //netcluster* pwc = NET->getwc(0);
    //detector* pD = NET->getifo(n);
    //double toffset = pwc->start+pD->waveForm.start();
    for (int i=0; i<size; i++) {
      double time = i*dt+vREC[n].start();

      double vl50 = vMED[n][i]-fabs(vL50[n][i]);
      double vu50 = vMED[n][i]+fabs(vU50[n][i]);
      double vl90 = vMED[n][i]-fabs(vL90[n][i]);
      double vu90 = vMED[n][i]+fabs(vU90[n][i]);

      double fl50 = fMED[n][i]-fabs(fL50[n][i]);
      double fu50 = fMED[n][i]+fabs(fU50[n][i]);
      double fl90 = fMED[n][i]-fabs(fL90[n][i]);
      double fu90 = fMED[n][i]+fabs(fU90[n][i]);

      out << time 
          << " " << vREC[n][i] << " " << vAVR[n][i] << " " << vRMS[n][i] 
          << " " << vMED[n][i] << " " << vl50 << " " << vl90 << " " << vu50 << " " << vu90
          << " " << fREC[n][i] << " " << fAVR[n][i] << " " << fRMS[n][i]
          << " " << fMED[n][i] << " " << fl50 << " " << fl90 << " " << fu50 << " " << fu90
          << endl;
    }

    out.close();
  }
}

// Dumps injected waveform/time array in ASCII format.
void DumpInjWavePE(network* NET, TString pdir) {

  int nIFO = NET->ifoListSize();	// number of detectors

  char ofname[256];
  for(int n=0; n<nIFO; n++) {

    sprintf(ofname,"%s/%s_inj_pe.dat",pdir.Data(),NET->ifoName[n]);
 
    ofstream out;
    out.open(ofname,ios::out);
    if (!out.good()) {cout << "Error Opening Output File : " << ofname << endl;exit(1);}
    cout << "Create Output File : " << ofname << endl;
    out.precision(19); 

    // write header
    out << "# time white_amp" << endl;

    // write data
    int size = vINJ[n].size();
    double dt=1./vINJ[n].rate();
    for (int i=0; i<size; i++) {
      double time = i*dt+vINJ[n].start();
      out << time << " " << vINJ[n][i] << endl;
    }

    out.close();
  }
}

int RedoAnalysis(TFile* jfile, CWB::config* cfg, network* NET) {

  // import global variables
  int gLRETRY=-1; IMPORT(int,gLRETRY)
  int gIFACTOR=-1; IMPORT(int,gIFACTOR)
  cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)

  if(gLRETRY==0) return 1;

  bool use_original_data = (gLRETRY==1) ? true : false;
  // when multitask only in the last trial the original data is used
  if(gOPT.multitask) use_original_data = (gMTRIAL==gOPT.trials) ? true : false;

  if(use_original_data) cout << endl << "Last Trial : Analysis of the original data" << endl << endl;

  int nRES = cfg->l_high-cfg->l_low+1;     		// number of frequency resolution levels
  int nIFO = NET->ifoListSize();  			// number of detectors

  vector<TString> delObjList;
  // supercluster clusters and parse maps are removed 
  delObjList.push_back("supercluster");
  delObjList.push_back("sparse");
  TString jname = jfile->GetPath();
  jname.ReplaceAll(":/","");
  jfile->Close();
  gCWB2G->FileGarbageCollector(jname,"",delObjList);

  // circular shift in the range [jS,jE] randomly whitened ifos HoT & add rec/inj waveforms into @ the same rec time
  int size,rate;
  int jS,jE,jW,k;
  WSeries<double>* pTF;
  for(int n=0; n<nIFO; n++) {                         	// create random time series

    // select waveform (reconstructed/injected) to be added to the whitened HoT
    wavearray<double> wf;
    if(gOPT.signal==0) wf=vREC[n];					// reconstructed
    if(gOPT.signal==1) wf=GetAlignedWaveform(&vINJ[n], &vREC[n]);	// injected
    if(gOPT.signal==2) wf=vDAT[n];					// reconstructed+null

    // apply amplitude mis-calibration amp_cal_err
    if((!use_original_data) && gOPT.amp_cal_err[n]) {	// in the last try (gLRETRY==1) we use the original data
      double amp_cal_err = 1;
      if(gOPT.amp_cal_err[n]>0) amp_cal_err = gRandom->Uniform(1-gOPT.amp_cal_err[n],1+gOPT.amp_cal_err[n]);
      else                      amp_cal_err = gRandom->Gaus(1,fabs(gOPT.amp_cal_err[n])); 
      cout << "Apply Amplitude Mis-Calibration (" << gOPT.amp_cal_err[n] <<"%)  " 
           << NET->ifoName[n] << " -> amp_cal_err : " << amp_cal_err << endl;
      wf*=amp_cal_err;
    }
    // apply phase mis-calibration phs_cal_err
    if((!use_original_data) && gOPT.phs_cal_err[n]) {	// in the last try (gLRETRY==1) we use the original data
      double phs_cal_err = 0;
      if(gOPT.phs_cal_err[n]>0) phs_cal_err = gRandom->Uniform(-gOPT.phs_cal_err[n],gOPT.phs_cal_err[n]);
      else                      phs_cal_err = gRandom->Gaus(0,fabs(gOPT.phs_cal_err[n]));
      cout << "Apply Phase     Mis-Calibration (" << gOPT.phs_cal_err[n] <<" deg)  "
           << NET->ifoName[n] << " -> phs_cal_err : " << phs_cal_err << " deg" << endl;
      CWB::mdc::PhaseShift(wf,phs_cal_err);
    }

    pTF = NET->getifo(n)->getTFmap();
    size = gHOT[n].size();
    rate = gHOT[n].rate();
    // apply time shift to input whitened data (integer number of sammples)
    jS = cfg->segEdge*rate;
    jE = size-jS;
    jW = rate*cfg->iwindow/2;
    k  = gRandom->Uniform(jS+jW,jE-jW); 		// select random offset (excludes +/- iwindow/2 around the event)
    if(use_original_data) {				// in the last try we restore the original data
      *gCWB2G->hot[n] = gHOT[n];
    } else {
      printf("Info : %s data time shift : %.5f sec\n",NET->getifo(n)->Name,(double)(k-jS)/rate);
      for(int i=jS; i<jE; i++)  {
        gCWB2G->hot[n]->data[k++] = gHOT[n].data[i];
        if(k==jE) k=jS;
      }
      *(gCWB2G->hot[n]) = AddWaveforms(gCWB2G->hot[n],&wf); // add waveform (reconstructed/injected) to whitened HoT
    }
    pTF->Forward(*(gCWB2G->hot[n]),*(NET->wdmList[nRES-1]));
  }

  // perform coherence and supercluster stages
  gCWB2G->Coherence(gIFACTOR);
  gCWB2G->SuperCluster(gIFACTOR);

  NET->setDelayIndex(gCWB2G->TDRate);
  if(gCWB2G->singleDetector) NET->setIndexMode(cfg->mode);  // when nIFO=1 exclude duplicate delay configurations

  // set low-rate TD filters
  for(int k=0; k<nRES; k++) gCWB2G->pwdm[k]->setTDFilter(cfg->TDSize, cfg->upTDF);
  NET->setDelayIndex(gCWB2G->TDRate);

  jfile = new TFile(jname);
  ReplaceSuperclusterData(jfile, cfg, NET, gOPT.gps);	    // save in jfile only max size supercluster

  // read sparse map from job file
  cout << "Loading sparse TF map ... " << endl;
  for(int n=0; n<nIFO; n++) {
    detector* pD = NET->getifo(n);
    pD->sclear();   // clear vector with sparse maps 
    TString ifo=NET->ifoName[n];
    for(int i=0; i<nRES; i++) {
      char swname[32];
      if(cfg->simulation) sprintf(swname,"sparse/%s-level:%d:%d",ifo.Data(),gIFACTOR,i+cfg->l_low);
      else                sprintf(swname,"sparse/%s-level:%d",ifo.Data(),i+cfg->l_low);
      SSeries<double>* psw = (SSeries<double>*)jfile->Get(swname);
      if(psw==NULL) {
        cout << "CWB_Plugin_PE - Likelihood : sparse map " << swname
             << " not exist in job file" << endl;exit(1);
      }
      SSeries<double> SS = *psw;
      pD->vSS.push_back(SS);
      delete psw;
    }
    cout<<endl;
  }

  netcluster* pwc = NET->getwc(0);
  pwc->cData.clear();
  // read cluster list & metadata netcluster object  
  int cycle = cfg->simulation ? gIFACTOR : Long_t(NET->wc_List[0].shift);
  vector<int> clist = pwc->read(jfile,"supercluster","clusters",0,cycle);
  return clist.size();
}

wavearray<double> AddWaveforms(wavearray<double>* wf1, wavearray<double>* wf2) {

   wavearray<double> wf = *wf1;

   if(wf1==NULL)      return wf;
   if(wf1->size()==0) return wf;

   double R=wf1->rate();

   double b_wf1 = wf1->start();
   double e_wf1 = wf1->start()+wf1->size()/R;
   double b_wf2 = wf2->start();
   double e_wf2 = wf2->start()+wf2->size()/R;

   int o_wf1 = b_wf1>b_wf2 ? 0 : int((b_wf2-b_wf1)*R+0.5);
   int o_wf2 = b_wf1<b_wf2 ? 0 : int((b_wf1-b_wf2)*R+0.5);

   double startXCOR = b_wf1>b_wf2 ? b_wf1 : b_wf2;
   double endXCOR   = e_wf1<e_wf2 ? e_wf1 : e_wf2;
   int sizeXCOR  = int((endXCOR-startXCOR)*R+0.5);

   for(int i=0;i<sizeXCOR;i++) wf[i+o_wf1] += wf2->data[i+o_wf2];

   return wf;
}

void SetEventWindow(CWB::config* cfg, double gps) {

  if(gps<0) return;

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}

  for(int n=0; n<cfg->nIFO; n++) {

    strcpy(cfg->DQF[cfg->nDQF].ifo, cfg->ifo[n]);
    sprintf(cfg->DQF[cfg->nDQF].file, "%s/%s_%s.gps_%d",cfg->tmp_dir,cfg->ifo[n],cfg->data_label,int(gps));
    cfg->DQF[cfg->nDQF].cat    = CWB_CAT2;
    cfg->DQF[cfg->nDQF].shift  = 0.;
    cfg->DQF[cfg->nDQF].invert = false;
    cfg->DQF[cfg->nDQF].c4     = true;
    cfg->nDQF++;

    cout << cfg->DQF[cfg->nDQF-1].file << endl; 

    ofstream out;
    out.open(cfg->DQF[cfg->nDQF-1].file,ios::out);
    cout << "Write file : " << cfg->DQF[cfg->nDQF-1].file << endl;
    if (!out.good()) {cout << "Error Opening File : " << cfg->DQF[cfg->nDQF-1].file << endl;exit(1);}
    out.precision(14);
    int istart = int(gps)-cfg->iwindow;
    int istop  = int(gps)+cfg->iwindow;
    out << "1 " << istart << " " << istop << " " << 2*cfg->iwindow << endl;
    out.close();
  }
}

skymap GetMedianSkyProb(network* NET) {
// compute median sky probability

  skymap skyprob = NET->getifo(0)->tau;
  skyprob=0;

  int ntry=0; for(int i=0;i<gOPT.trials;i++) if(wREC[i].size()) ntry++; 	// number of detected events 

  int *index = new int[ntry];
  float *value = new float[ntry];

  int L = skyprob.size();
  double sum=0;
  for(int l=0;l<L;l++) {		// for each sky location compute median sky probability
    
    int k=0; for(int i=0;i<gOPT.trials;i++) if(wREC[i].size()) value[k++] = wSKYPROB[i].get(l); 
    TMath::Sort(ntry,value,index,false);

    int imed = (ntry*50.)/100.; if(imed>=ntry) imed=ntry-1;

    skyprob.set(l,value[index[imed]]);	// save median sky probability

    sum+=value[index[imed]]; 
  }

  // normalize sky probability
  if(sum>0) for(int l=0;l<L;l++) skyprob.set(l,skyprob.get(l)/sum);

  delete [] index;
  delete [] value;

  return skyprob;
}

void DumpSkyProb(skymap* skyprob, network* NET, netevent* &EVT, TString odir) {

  skymap skyprobcc = NET->getifo(0)->tau; 

  // convert skyprob in celestial coordinates
  int L = skyprobcc.size();
  for(int l=0; l<L; l++) {
    double th = skyprob->getTheta(l);
    double ph = skyprob->getPhi(l);

    double ra = skyprob->getRA(l);
    int k=skyprob->getSkyIndex(th, ra);
    skyprobcc.set(k,skyprob->get(l));
  }

  // dump skyprobcc to fits file
  char fname[1024];
  if(skyprobcc.getType()) { // check if it is healpix 
    sprintf(fname, "%s/mskyprobcc.%s", odir.Data(), "fits");

    // build fits configur info
    TString analysis = "2G";
    if(NET->MRA) analysis+=":MRA";
    if(NET->pattern==0)                     analysis+=":Packet(0)";
    if((NET->pattern!=0 && NET->pattern<0)) analysis+=TString::Format(":Packet(%d)",NET->pattern);
    if((NET->pattern!=0 && NET->pattern>0)) analysis+=TString::Format(":Packet(+%d)",NET->pattern);

    char configur[64]="";
    char search = NET->tYPe;
    if (search=='r')                 sprintf(configur,"%s un-modeled",analysis.Data());
    if (search=='i')                 sprintf(configur,"%s iota-wave",analysis.Data());
    if (search=='p')                 sprintf(configur,"%s psi-wave",analysis.Data());
    if((search=='l')||(search=='s')) sprintf(configur,"%s linear",analysis.Data());
    if((search=='c')||(search=='g')) sprintf(configur,"%s circular",analysis.Data());
    if((search=='e')||(search=='b')) sprintf(configur,"%s elliptical",analysis.Data());
    skyprobcc.Dump2fits(fname,EVT->time[0],configur,const_cast<char*>("PROB"),const_cast<char*>("pix-1"),'C');
  }

}

void LoadFromMultiTaskJobOutput(int runID, CWB::config* cfg) {

  int nIFO = cfg->nIFO;

  char beginString[1024];
  sprintf(beginString,"wave_");	
  char endString[1024];
  sprintf(endString,"_job%d.root",runID);	
  char containString[1024];
  sprintf(containString,"%s_trial",cfg->data_label);

  vector<TString> fileList = CWB::Toolbox::getFileListFromDir(cfg->output_dir,endString,beginString,containString,containString);

  wavearray<double>* mtpe_wREC[NIFO_MAX];
  for(int i=0;i<nIFO;i++) mtpe_wREC[i] = new wavearray<double>;
  skymap* mtpe_skyprob = new skymap(int(0));
  float* chirp = new float[6];
  float* isnr = new float[nIFO];
  float* osnr = new float[nIFO];
  float* iosnr = new float[nIFO];
  float likelihood;

  char command[1024];
  int nfile = fileList.size();
  for(int n=0;n<nfile;n++) {
    cout << n << " " << fileList[n].Data() << endl;

    TFile* froot = new TFile(fileList[n].Data(),"READ");
    if(froot==NULL) {
      cout << "CWB_Plugin_PE Error : Failed to open file : " <<  fileList[n].Data() << endl;
      gSystem->Exit(1);
    }
    TTree* itree = (TTree*)froot->Get("waveburst");
    if(itree==NULL) {
      cout << "CWB_Plugin_PE Error : Failed to open tree waveburst from file : " <<  fileList[n].Data() << endl;
      gSystem->Exit(1);
    }

    for(int i=0;i<nIFO;i++) itree->SetBranchAddress(TString::Format("mtpe_wREC_%d",i).Data(),&mtpe_wREC[i]);
    itree->SetBranchAddress("mtpe_skyprob",&mtpe_skyprob);
    itree->SetBranchAddress("chirp",chirp);
    itree->SetBranchAddress("likelihood",&likelihood);

    itree->GetEntry(0);						// read wREC,skyprob objects

    // check if objects are not null
    for(int i=0;i<nIFO;i++) {
      if(mtpe_wREC[i]==NULL) {
        cout << "CWB_Plugin_PE Error : Object wavearray not exist !!! " <<  endl;
        gSystem->Exit(1);
      }
    }
    if(mtpe_skyprob==NULL) {
      cout << "CWB_Plugin_PE Error : Object skymap not exist !!! " <<  endl;
      gSystem->Exit(1);
    }

    // fill object vectors
    wSKYPROB[gOPT.trials-n-1] = *mtpe_skyprob;   		// restore sky probability

    for(int i=0;i<nIFO;i++) 
      wREC[gOPT.trials-n-1].push_back(*mtpe_wREC[i]);  		// restore reconstructed waveform

    for(int i=0;i<nIFO;i++) {					// restore SNR
       iSNR[i].push_back(isnr[i]);
       oSNR[i].push_back(osnr[i]);
      ioSNR[i].push_back(iosnr[i]);
    }

    for(int j=0; j<6; j++) vCHIRP[j].push_back(chirp[j]);	// restore chirp mass parameters

    vLIKELIHOOD.push_back(likelihood);				// restore likelihood

    froot->Close();

    // remove temporary root file
    sprintf(command,"/bin/rm %s",fileList[n].Data());
    cout << command << endl;
    gSystem->Exec(command);
  }

  for(int i=0;i<nIFO;i++) delete mtpe_wREC[i];
  delete mtpe_skyprob;
  delete [] chirp;
  delete [] isnr;
  delete [] osnr;
  delete [] iosnr;
}

void GetNoisePixels(std::vector<double>* aNSE, std::vector<double>* ANSE, 
                    network* NET, CWB::config* cfg, int lag, int id) {
// get noise statistic in TF domain using random event TF patterns 
// this is used for Kolmogorov test of the residuals

  size_t nIFO = NET->ifoList.size();

  netcluster* pwc = NET->getwc(lag);
  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, id, false); // buffer for pixel IDs
  netpixel* pix = pwc->getPixel(id,pI[0]);
  int V = pI.size();
  if(!V) return;

  int nRES = NET->wdmListSize();               	// number of frequency resolution levels

  double R = NET->getifo(0)->TFmap.rate();

  double rate = gHOT[0].rate();
  int    size = gHOT[0].size();

  int layers_high = 1<<cfg->l_high;

  // compute the minimum index in the super cluster
  int index_min=2*gHOT[0].size();
  for(int k=0; k<V; k++) {
    netpixel* pix = pwc->getPixel(id,pI[k]);
    for(int n=0; n<nIFO; n++) {
      int index = pix->data[n].index;
      if(index<index_min) index_min=index;
    }
  }
  index_min = index_min-(index_min%layers_high); // set index_min a multiple of layers_high	

  // find range where to select noise data
  int jS = cfg->segEdge*rate;
  int jE = size-jS;
  int jW = rate*cfg->iwindow;

  WSeries<double>  w;
  for(int n=0; n<nIFO; n++) {
    for(int j=0;j<nRES;j++) {			 // j=0 -> l_high : j=nRES-1 -> l_low
      w.Forward(gHOT[n],*(NET->wdmList[nRES-1-j]));
      for(int m=0; m<nTRIALS; m++) {
        int index_off = gRandom->Uniform(jS,jE-jW); 	// select random offset (excludes last iwindow range)
        index_off = index_off-(index_off%layers_high); 	// set index_off a multiple of layers_high	
        for(int k=0; k<V; k++) {
          netpixel* pix = pwc->getPixel(id,pI[k]);
          int index = pix->data[n].index;
          int ires = int(TMath::Log2(R/pix->rate))-cfg->l_low;
          int ilayers = 1<<(ires+cfg->l_low);
          index = index_off+(index-index_min);	// shift pixels
          if(ires==j) {
            // get original sparse map amplitudes
            //double dd = GetSparseMapData(&vSS[n][ires], true, index);
            //double DD = GetSparseMapData(&vSS[n][ires], false , index);
            // get the shifted amplitudes
            double dd = w.data[index];
            double DD = w.data[index+w.maxIndex()+1];
            aNSE[n].push_back(dd); 
            ANSE[n].push_back(DD); 
          }
        }
      }
    }
  } 
  return;
}

void PlotSpectrogram(TString type, network* NET, netevent* &EVT, CWB::config* cfg, TString pdir) {

  if(type!="data" && type!="null" && type!="signal" && type!="injection") {
    cout << "CWB_Plugin_PE Error : wrong PlotSpectrogram type, must be data/null/signal" << endl;
    exit(1);
  }  

  int nIFO = NET->ifoListSize(); 		// number of detectors
  char fname[1024];

  for(int n=0; n<nIFO; n++) {
    WSeries<double>* pTF = NET->getifo(n)->getTFmap();
    if(type=="null") {
      wavearray<double> wf = vREC[n];
      wf*=-1; 
      *pTF = AddWaveforms(&gHOT[n],&wf);	// sub waveform (reconstructed/injected) to whitened HoT
    }
    if(type=="signal") {
      wavearray<double> wf = vREC[n];
      *pTF=0;
      *pTF = AddWaveforms(pTF,&wf); 		// reconstructed signal
    }
    if(type=="injection") {
      wavearray<double> wf = vINJ[n];
      *pTF=0;
      *pTF = AddWaveforms(pTF,&wf); 		// injected signal
    }

    *pTF*=1./sqrt(1<<cfg->levelR);

    TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",EVT->gps[n]);

    if(pTF->size()==0) continue;
    if(pTF->getLevel()>0) pTF->Inverse();

    int nfact=4;
    int nfft=nfact*512;
    int noverlap=nfft-10;
    double fparm=nfact*6;
    int ystart = int((EVT->start[n]-EVT->gps[n]-1)*pTF->rate());
    int ystop  = int((EVT->stop[n]-EVT->gps[n]+1)*pTF->rate());
    ystart-=nfft;
    ystop+=nfft;
    int ysize=ystop-ystart;
    wavearray<double> y;y.resize(ysize);y.rate(pTF->rate());y.start(ystart/pTF->rate());

    // stft use dt=y.rate() to normalize data but whitened data are already normalized by dt 
    // so before stft data must be divided by 1./sqrt(dt)
    for(int i=0;i<(int)y.size();i++) y.data[i]=pTF->data[i+ystart]/sqrt(y.rate());

    CWB::STFT stft(y,nfft,noverlap,"energy","gauss",fparm);

    TCanvas* canvas;
    double tstart = nfft/pTF->rate()+ystart/pTF->rate();
    double tstop = (ysize-nfft)/pTF->rate()+ystart/pTF->rate();

    tstart+=0.9;tstop-=0.9;
    stft.Draw(tstart,tstop,pTF->getlow(),pTF->gethigh(),0,0,1);
    stft.GetHistogram()->GetXaxis()->SetTitle(xtitle);
    sprintf(fname, "%s/%s_%s_spectrogram_0.png", pdir.Data(), NET->ifoName[n], type.Data());
    canvas = stft.GetCanvas();
    stft.Print(fname);
    canvas->SetLogy(true);
    stft.GetHistogram()->GetXaxis()->SetTitle(xtitle);
    sprintf(fname, "%s/%s_%s_spectrogram_logy_0.png", pdir.Data(), NET->ifoName[n], type.Data());
    stft.Print(fname);

    y.resize(0);
  } 
}

std::vector<wavearray<double> > GetWhitenedData(network* NET, CWB::config* cfg) {
  // get whitened data (in the vREC time range)

  std::vector<wavearray<double> > xWHT;        // temporary stuff

  int nIFO = NET->ifoListSize();               // number of detectors

  for(int n=0; n<nIFO; n++) {
    xWHT.push_back(GetAlignedWaveform(&gHOT[n], &vREC[n])); 
  }
  return xWHT;
}

// the following function replace the supercluster clusters with the max size supercluster  
void ReplaceSuperclusterData(TFile*& jfile, CWB::config* cfg, network* NET, double gps) {
  // gps=0 -> select cluster with max size
  // gps>0 -> select cluster nearest to gps time

  int gIFACTOR=-1; IMPORT(int,gIFACTOR)
  cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)

  int cycle = cfg->simulation ? gIFACTOR : Long_t(NET->wc_List[0].shift);

  netcluster* pwc = NET->getwc(0);
  pwc->clear();

  // read cluster list & metadata netcluster object  
  vector<int> clist = pwc->read(jfile,"supercluster","clusters",0,cycle);
  //pwc->print();

  if(clist.size()==0) return;

  // find max size supercluster index
  int maxk=0;
  int npixels=0;                		// total loaded pixels per lag
  int msize=0;
  for(int k=0;k<(int)clist.size();k++) {      	// loop over the cluster list
    // read pixels & tdAmp into netcluster pwc
    pwc->read(jfile,"supercluster","clusters",-2,cycle,0,clist[k]);
    int psize = pwc->size()-npixels;
    if(psize>msize) {maxk=k;msize=psize;}
    npixels+=psize;
  }

  // find min abs(gps-time) supercluster index
  int mink=0;
  double mtime=1.e20;
  wavearray<double> ctime = pwc->get((char*)"time",0,'L',0);            // get cluster time
  for(int k=0; k<ctime.size(); k++) {                                   // loop over clusters
    if(fabs(ctime[k]+pwc->start-gps)<mtime) {mink=k;mtime=fabs(ctime[k]+pwc->start-gps);}
  }

  int kindex = gps>0 ? mink : maxk;

  pwc->clear();
  // read max supercluster (pixels & tdAmp) into netcluster pwc
  pwc->read(jfile,"supercluster","clusters",-1,cycle,0,clist[kindex]);

  // supercluster are removed from jfile
  vector<TString> delObjList;
  delObjList.push_back("supercluster");
  TString jname = jfile->GetPath();
  jname.ReplaceAll(":/","");
  jfile->Close();
  gCWB2G->FileGarbageCollector(jname,"",delObjList);

  // max supercluster is stored in jfile
  jfile = new TFile(jname,"UPDATE");
  gCWB2G->jfile=jfile;
  if(jfile==NULL||!jfile->IsOpen())
      {cout << "CWB_Plugin_PE.C - Error opening : " << jname <<  endl;exit(1);}
  pwc->write(jfile,"supercluster","clusters",0,cycle);
  pwc->write(jfile,"supercluster","clusters",-1,cycle);

  jfile->Write();
}
