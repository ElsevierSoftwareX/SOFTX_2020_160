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

#define _USE_HEALPIX

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
#include "gskymap.hh"
#include <fstream>
#include <vector>

// ---------------------------------------------------------------------------------
// HOW TO SETUP PLUGIN WF IN CWB USER CONFIGURATION (EXAMPLE)
// ---------------------------------------------------------------------------------

/*
  TString optwf = "";                  	    // NOTE : add space at the end of each line

  optwf += "wf_output_disable=root ";       // disable output root file (to be used when QLveto is enabled)
  optwf += "wf_output_disable=inj ";        // disable save injection to the output root file
  optwf += "wf_output_enable=rec ";         // enable  save reconstructed waveform to the output root file
  optwf += "wf_output_disable=wht ";        // disable save whitened data to the output root file
  optwf += "wf_output_disable=dat ";        // disable save rec+null data to the output root file
  optwf += "wf_output_disable=nul ";        // disable save null data to the output root file
  optwf += "wf_inj_tstep=150.000 ";         // is the the injection step time 

  strcpy(parPlugin,optwf.Data());           // set WF plugin parameters
*/

// ---------------------------------------------------------------------------------
// DEFINES
// ---------------------------------------------------------------------------------

#define WF_VERSION		1.2

#define WF_OUTPUT_ROOT  	false 

#define WF_OUTPUT_INJ           false           // save injected      into the output root file
#define WF_OUTPUT_REC           true            // save reconstructed into the output root file
#define WF_OUTPUT_WHT           false           // save whitened data into the output root file
#define WF_OUTPUT_DAT           false           // save rec+null data into the output root file
#define WF_OUTPUT_NUL           false           // save null data into the output root file

#define WF_OUTPUT_STRAIN        false           // save full strain  data buffer into the output root file
#define WF_OUTPUT_MDC           false           // save full strain mdc injections data buffer into the output root file
#define WF_OUTPUT_CSTRAIN       false           // save full cstrain data buffer into the output root file

#define WF_OUTPUT_CLUSTER       false           // save event crate/cluster
#define WF_OUTPUT_POL           false           // save event polarization components (dual stream polarization frame)

#define WF_FITS_PE		""
#define WF_GPS_PE		0

#define WF_FITS_OS		""
#define WF_GPS_OS		0

#define WF_INJ_TSTEP		0		// is the injection step time
#define WF_INJ_TOFF		0		// is time offset (sec) of periodic time step injections

//#define PLOT_INJ    
//#define TEST1
//#define TEST2
//#define TEST3
//#define DUMP_FITS

#define FRACTION		0.01

// ---------------------------------------------------------------------------------
// USER SGW PLUGIN OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  bool   output_root;   // (def=false) set to false if this plugin is used with QLveto (QLveto plugin output par to output root file)
                        // NOTE: QLveto plugin must be declared as the last one in cwb_mplugin
  bool   output_inj;    // save injected      into the output root file
  bool   output_rec;    // (def=true) save reconstructed into the output root file
  bool   output_wht;    // save whitened data into the output root file
  bool   output_dat;    // save rec+null data into the output root file
  bool   output_nul;    // save null data into the output root file

  bool   output_strain; // save full strain  data buffer into the output root file
  bool   output_mdc;    // save full strain mdc injections data buffer into the output root file
  bool   output_cstrain;// save full cstrain data buffer into the output root file

  bool   output_cluster;// save event crate/cluster into the output root file
  bool   output_pol;	// save event polarization components into the output root file

  TString fits_pe;	// fits file name -> reference skymap: Parameter Estimation  
  double  gps_pe;	// gps  Parameter Estimation  gps time

  TString fits_os;	// fits file name -> reference skymap: cWB On Source 
  double  gps_os;	// gps On Source gps time

  double  inj_tstep; 	// is the injection time step 
  double  inj_toff; 	// is time offset (sec) of periodic time step injections
};

// ----------------------------------------------------
// ROOT Output WF Parameters
// ----------------------------------------------------

struct WF {                             // structure for output root parameters

  wavearray<double>* wINJ[NIFO_MAX];	// whiten injected waveforms
  wavearray<double>* wAUX[NIFO_MAX];	// whiten auxiliary injected waveforms
  wavearray<double>* cINJ[NIFO_MAX];	// injected (cutted) on TF sub-space of the reconstructed waveform vREC
  wavearray<double>* wREC[NIFO_MAX];	// whiten reconstructed waveforms
  wavearray<double>* wWHT[NIFO_MAX];	// whiten data (in the waveform time range)
  wavearray<double>* wDAT[NIFO_MAX];	// whiten rec+noise data (in the waveform time range)
  wavearray<double>* wNUL[NIFO_MAX];	// whiten null data (in the waveform time range)

  std::vector<netpixel>* pCLUSTER;      // event netcluster pointer

  wavearray<double>  rPOL00; 		// polarization of 00 amplitudes -> radial
  wavearray<double>  aPOL00; 		// polarization of 00 amplitudes -> angle
  wavearray<double>  rPOL90; 		// polarization of 90 amplitudes -> radial
  wavearray<double>  aPOL90; 		// polarization of 90 amplitudes -> angle
};

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions           gOPT;                // global User Options
TTree*             gTREE;               // output tree file name
TString            gOUTPUT;             // output root file name
WF                 gWF;                 // global ROOT Output WF Parameters
double             gSEGGPS;             // segment start time 
skymap             gSM_PE;		// reference skymap: Parameter Estimation
skymap             gSM_OS;		// reference skymap: cWB On Source
float		   gSM_PE_OF;		// skymap PE overlap factor
float		   gSM_OS_OF;		// skymap OS overlap factor
bool               gCEDDUMP;            // save user_parameters cedDump
float              gFF[3];		// Fitting Factor W1,W2    0/1/2 -> full/left/right
float              gOF[3];		// Overlap Factor W1,W2    0/1/2 -> full/left/right
float              gRE[3];		// residual energy W1-W2   0/1/2 -> full/left/right

wavearray<double> gSTRAIN[NIFO_MAX];    // HoT time series used to save strain data
wavearray<double> gMDC[NIFO_MAX];       // HoT time series used to save mdc strain data
wavearray<double> gCSTRAIN[NIFO_MAX];   // HoT time series used to save whitened data

wavearray<double> gHOT[NIFO_MAX];       // HoT time series used to save whitened data

double		    gTCOA;		// Coalescence Time
std::vector<double> gvTCOA;		// vector of gTCOA: contains the gTCOA of all injections in the job
string		    gLOG;		// List of input parameters used to create the injected waveform
std::vector<string> gvLOG;		// vector of gLOG: contains the gLOG of all injections in the job

float           	gCRATE;		// event netcluster::rate
std::vector<netpixel> 	gCLUSTER;	// event netcluster

// ---------------------------------------------------------------------------------
// waveforms (vectors index is ifos)
// ---------------------------------------------------------------------------------

std::vector<wavearray<double> > vINJ;           // injected (full)
std::vector<wavearray<double> > vAUX;           // auxiliary injected (full)
std::vector<wavearray<double> > cINJ;           // injected (cutted) on TF sub-space of the reconstructed waveform vREC
std::vector<wavearray<double> > vREC;           // signal reconstructed by wave packet 
std::vector<wavearray<double> > vWHT;           // whitened data (in the vREC time range) 
std::vector<wavearray<double> > vDAT;           // noise+signal reconstructed by wave packet
std::vector<wavearray<double> > vNUL;           // vDAT-vREC
std::vector<wavearray<double> > vRES;           // vREC-cINJ

// ---------------------------------------------------------------------------------
// sparse maps
// ---------------------------------------------------------------------------------

std::vector<SSeries<double> > vSS[NIFO_MAX];    // original sparse maps
std::vector<SSeries<double> > jSS[NIFO_MAX];    // injected (only TF sub-space of the reconstructed waveform vREC)
std::vector<SSeries<double> > rSS[NIFO_MAX];    // reconstructed
std::vector<SSeries<double> > dSS[NIFO_MAX];    // noise+signal reconstructed by wave packet
std::vector<SSeries<double> > xSS[NIFO_MAX];    // sparse map -> vRES=vREC-cINJ
 
// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);

void ClearWaveforms(detector* ifo);
void ClearVectors();

void SetOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, bool dump_wf);
void DumpOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor);

std::vector<netpixel> GetCluster(network* NET, int lag, int id);

std::vector<wavearray<double> > GetWaveform(network* NET, int lag, int id, char type, bool shift=true);
std::vector<wavearray<double> > GetWhitenedData(network* NET, CWB::config* cfg);
std::vector<wavearray<double> > GetInjWaveform(network* NET, netevent* EVT, int id, double factor);
std::vector<wavearray<double> > GetAuxWaveform(network* NET, netevent* EVT, int id, double factor);
std::vector<wavearray<double> > GetRecWaveform(network* NET, netevent* EVT, int id);

std::vector<wavearray<double> > GetSigWaveform(network* NET, CWB::config* cfg, int lag, int id);
void Wave2Sparse(network* NET, CWB::config* cfg, char type);
double GetSparseMapData(SSeries<double>* SS, bool phase, int index);

void DumpSkymap(network* NET, int lag, netevent* EVT, int ID);

void PlotWaveform(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wf1,
                  wavearray<double>* wf2, wavearray<double>* wf3, wavearray<double>* wref,
                  bool fft=false, TString pdir="", double P=0.99, EColor col1=kRed, EColor col2=kBlue, EColor col3=(EColor)418);

skymap ReadSkyMap(TString fitsFile, CWB::config* cfg, double gps);
double ComputeSkyMapOF(skymap* sm1, skymap* sm2);
double ComputeMaxSkyMapProduct(skymap* sm1, skymap* sm2);
skymap GetSkyMap(network* NET, int lag, netevent* EVT, int ID);
skymap EarthCoordinatesSM(skymap* sm, double gps);

std::vector<netpixel> SavePixels(network* NET, CWB::config* cfg, int lag, int id);
void RestorePixels(network* NET, CWB::config* cfg, netcluster* pwc, std::vector<netpixel>* vPIX, int id);

double GetInjTcoa(double geocentric_tcoa, network* NET, TString ifo, double theta, double phi);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Plugin used to save inj/rec/wht waveforms to output ROOT file

  if(type==CWB_PLUGIN_CONFIG) {  

    ResetUserOptions();                                 // set default config options
    ReadUserOptions(cfg->parPlugin);                    // user config options : read from parPlugin

    cfg->outPlugin=true;  				// disable built-in output root file

    gCEDDUMP=cfg->cedDump;                              // save user_parameter cedDump 

    if(gOPT.fits_pe!="" && gOPT.gps_pe<=0) {
        cout<< "CWB_Plugin_xWF - Error : plugin parameter wf_gps_pe not defined" << endl;exit(1);
    }
    if(gOPT.fits_os!="" && gOPT.gps_os<=0) {
        cout<< "CWB_Plugin_xWF - Error : plugin parameter wf_gps_os not defined" << endl;exit(1);
    }
    if(gOPT.inj_tstep<=0) {
        cout<< "CWB_Plugin_xWF - Error : plugin parameter inj_tstep not defined" << endl;exit(1);
    }

    if(gOPT.fits_pe!="") {CWB::Toolbox::checkFile(gOPT.fits_pe); gSM_PE = ReadSkyMap(gOPT.fits_pe, cfg, gOPT.gps_pe);}
    if(gOPT.fits_os!="") {CWB::Toolbox::checkFile(gOPT.fits_os); gSM_OS = ReadSkyMap(gOPT.fits_os, cfg, gOPT.gps_os);}
  }

  if(type==CWB_PLUGIN_NETWORK) {
    PrintUserOptions(cfg);      // print config options
  }

  if(type==CWB_PLUGIN_STRAIN) {
    // save strain
    int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==NET->getifo(n)->Name) {ifoID=n;break;}
    if(gOPT.output_strain) {
      gSTRAIN[ifoID] = *x;
      if(cfg->fResample>0) gSTRAIN[ifoID].Resample(cfg->fResample);      // RESAMPLING
      gSTRAIN[ifoID].Resample(gSTRAIN[ifoID].rate()/(1<<cfg->levelR));   // resampling
    } 
  }

  if(type==CWB_PLUGIN_MDC) {
    // save mdc strain data
    int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==NET->getifo(n)->Name) {ifoID=n;break;}
    if(gOPT.output_mdc) {
      gMDC[ifoID] = *x;
      if(cfg->fResample>0) gMDC[ifoID].Resample(cfg->fResample);      	// RESAMPLING
      gMDC[ifoID].Resample(gMDC[ifoID].rate()/(1<<cfg->levelR));   	// resampling
    }
  }

  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {
    // save whitened HoT
    int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==NET->getifo(n)->Name) {ifoID=n;break;}
    if(gOPT.output_wht) gHOT[ifoID] = *x;
    if(gOPT.output_cstrain) gCSTRAIN[ifoID] = *x;
  }

  if(type==CWB_PLUGIN_MDC && cfg->simulation) {
    // save coalescence times
    for(int k=0;k<(int)NET->mdcTime.size();k++) gvTCOA.push_back(NET->mdcTime[k]);
    // save logs
    for(int k=0;k<(int)NET->mdcList.size();k++) {
      gvLOG.push_back(NET->mdcList[k]);
      gvLOG[k].erase(std::remove(gvLOG[k].begin(), gvLOG[k].end(), '\n'), gvLOG[k].end()); // remove new line
    }
  }

  if(type==CWB_PLUGIN_XLIKELIHOOD) {			// BEFORE EVENT RECONSTRUCTION 
    cfg->cedDump = true;	// enable skymaps filling
  }

  if(type==CWB_PLUGIN_ILIKELIHOOD) {                    // DEFINE WAVE TREE - FIX CASE WHEN NO EVENTS HAVE BEEN FOUND !!!
    netevent* EVT;
    SetOutputFile(NET, EVT, cfg, false);                // set default output waveburst tree in root file
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {			// AFTER EVENT RECONSTRUCTION 

    cfg->cedDump=gCEDDUMP;      // restore user_parameter cedDump 

    cout << "CWB_Plugin_WF.C " << endl;

    ClearVectors();

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

    SetOutputFile(NET, EVT, cfg, false);               	// set output root file

    injection INJ(nIFO);

    for(int k=0; k<K; k++) {  				// loop over the lags

      id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);

      for(int j=0; j<(int)id.size(); j++) {  		// loop over cluster index

        int ID = size_t(id.data[j]+0.5);

        if(NET->getwc(k)->sCuts[ID-1]!=-1)  continue; 	// skip rejected/processed clusters

        EVT->output2G(NULL,NET,ID,k,ofactor);		// get reconstructed parameters

	// get coalescence time
        if(cfg->simulation) {
          double minDT = 1.e12;
          int injID = -1;
          int M = NET->mdc__IDSize();
          for(int m=0; m<M; m++) {
             int mdcID = NET->getmdc__ID(m);
             double DT = fabs(EVT->time[0] - NET->getmdcTime(mdcID));
             if(DT<minDT && INJ.fill_in(NET,mdcID)) {
                minDT = DT;
                injID = mdcID;
             }
          }
    	  gTCOA = gvTCOA[injID];
	  gLOG  = gvLOG[injID];
        }

        // print event parameters
        cout << endl;
        cout << "CWB_Plugin_WF - event parameters : ID -> " << ID << endl;
        for(int n=0;n<nIFO;n++) printf("rec_time %s : %.4f\n",NET->ifoName[n], EVT->time[n]);
        cout << "rec_theta : " << EVT->theta[0] << " rec_phi : " << EVT->phi[0] << endl;
        cout << "SNRnet : " << sqrt(EVT->likelihood) << " netcc[0] : " << EVT->netcc[0] 
             << " rho[0] : " << EVT->rho[0] << " size : " << EVT->size[0] << endl;

        // save event parameters
        gSEGGPS = EVT->gps[0];                          // save segment start time (sec) 

        skymap skyprob = GetSkyMap(NET, k, EVT, ID);
        skyprob = EarthCoordinatesSM(&skyprob, EVT->time[0]);
        if(gOPT.fits_pe!="") {
          gSM_PE_OF = ComputeSkyMapOF(&skyprob, &gSM_PE);
          cout << "SKYMAP PE Overlap Factor = " << gSM_PE_OF << endl;
        }
        if(gOPT.fits_os!="") {
          gSM_OS_OF = ComputeSkyMapOF(&skyprob, &gSM_OS);
          cout << "SKYMAP OS Overlap Factor = " << gSM_OS_OF << endl;
        }

        cout << "MAX SKYMAP Correlation " << NET->nCorrelation.max() << endl;
        double maxcc = ComputeMaxSkyMapProduct(&NET->nCorrelation, &NET->nProbability);
        cout << "MAX SKYMAP Correlation*CumProbability " << maxcc << endl;
        if(gOPT.fits_pe!="") {
          double maxcc_pe = ComputeMaxSkyMapProduct(&NET->nCorrelation, &gSM_PE);
          cout << "MAX SKYMAP Correlation*CumProbability PE " << maxcc_pe << endl;
        }
        if(gOPT.fits_os!="") {
          double maxcc_os = ComputeMaxSkyMapProduct(&NET->nCorrelation, &gSM_OS);
          cout << "MAX SKYMAP Correlation*CumProbability OS " << maxcc_os << endl;
        }

#ifdef DUMP_FITS
        DumpSkymap(NET, k, EVT, ID);			// save skymap to fits file
#endif

        // get waveforms
        if(cfg->simulation) {
          vINJ = GetInjWaveform(NET, EVT, ID, factor);  // get injected waveform
          if(vINJ.size()!=nIFO) {
            cout << "CWB_Plugin_WF Error : Injection Waveform Not Found !!!" << endl;
            exit(1);
          }
          vAUX = GetAuxWaveform(NET, EVT, ID, factor);  // get injected waveform
        }
        vREC = GetWaveform(NET, k, ID, 'S');  // get reconstructed waveform
        vWHT = GetWhitenedData(NET, cfg);     // get whitened data (in the vREC time range)
        vDAT = GetWaveform(NET, k, ID, 'W');  // get reconstructed+noise waveform
        vNUL = GetWaveform(NET, k, ID, 'N');  // get noise-reconstructed waveform
        for(int i=0; i<nIFO; i++) vDAT[i] = CWB::mdc::GetAligned(&vDAT[i], &vREC[i]);

#ifdef TEST1
        double sync_time;
        vREC[1]*=-1.0;
        double sync_xcor = CWB::mdc::TimeSync(vREC[0], vREC[1], sync_time);
        cout << ID << " vREC[0]xvREC[1] --------------------------------> sync_time " << sync_time << " sync_xcor : " << sync_xcor << endl;
#endif

        Wave2Sparse(NET,cfg,'v');                       // save original sparse map to vSS
        Wave2Sparse(NET,cfg,'r');                       // rec to rSS
        Wave2Sparse(NET,cfg,'d');                       // dat to dSS
        if(cfg->simulation) Wave2Sparse(NET,cfg,'j');   // inj to jSS

        // get injected waveform on TF sub-space of the reconstructed waveform vREC
        if(cfg->simulation) {
          cINJ = GetSigWaveform(NET, cfg, k, ID);
          for(int i=0; i<nIFO; i++) vINJ[i] = CWB::mdc::GetAligned(&vINJ[i], &vREC[i]);
          for(int i=0; i<nIFO; i++) cINJ[i] = CWB::mdc::GetAligned(&cINJ[i], &vREC[i]);
          for(int i=0; i<nIFO; i++) vRES.push_back(CWB::mdc::GetDiff(&vREC[i], &cINJ[i]));  // vRES=vREC-cINJ
          Wave2Sparse(NET,cfg,'x');                     // res to xSS
        }

        if(cfg->simulation) {
          double FF,OF,RE;
          vector<double> tstart;
          vector<double> tstop;

          // get FF/OF/RE for the full time range
          gFF[0] = CWB::mdc::GetMatchFactor("ff",vREC,vINJ,tstart,tstop);
          gOF[0] = CWB::mdc::GetMatchFactor("of",vREC,vINJ,tstart,tstop);
          gRE[0] = CWB::mdc::GetMatchFactor("re",vREC,vINJ,tstart,tstop);

          // get FF/OF/RE for the left time range (before the coalescence time)
          tstart.resize(0);
          tstop.resize(nIFO);
          for(int n=0;n<nIFO;n++) {
            tstop[n] = GetInjTcoa(gTCOA, NET, NET->ifoName[n], EVT->theta[1], EVT->phi[1]);
          }
          gFF[1] = CWB::mdc::GetMatchFactor("ff",vREC,vINJ,tstart,tstop);
          gOF[1] = CWB::mdc::GetMatchFactor("of",vREC,vINJ,tstart,tstop);
          gRE[1] = CWB::mdc::GetMatchFactor("re",vREC,vINJ,tstart,tstop);

          // get FF/OF/RE for the right time range (after the coalescence time)
          tstart.resize(nIFO);
          tstop.resize(0);
          for(int n=0;n<nIFO;n++) {
            tstart[n] = GetInjTcoa(gTCOA, NET, NET->ifoName[n], EVT->theta[1], EVT->phi[1]);
          }
          gFF[2] = CWB::mdc::GetMatchFactor("ff",vREC,vINJ,tstart,tstop);
          gOF[2] = CWB::mdc::GetMatchFactor("of",vREC,vINJ,tstart,tstop);
          gRE[2] = CWB::mdc::GetMatchFactor("re",vREC,vINJ,tstart,tstop);
        }

        gCRATE   = NET->getwc(k)->rate;			// get event cluster rate
        gCLUSTER = GetCluster(NET, k, ID);		// get event netcluster

        if(gOPT.output_pol) {
          gWF.rPOL00 = NET->r00_POL[0];			// polarization of 00 amplitudes -> radial
          gWF.aPOL00 = NET->r00_POL[1];			// polarization of 00 amplitudes -> angle
          gWF.rPOL90 = NET->r90_POL[0];			// polarization of 90 amplitudes -> radial
          gWF.aPOL90 = NET->r90_POL[1];			// polarization of 90 amplitudes -> angle
        }

#ifdef TEST2
        wavearray<double> vDIF[NIFO_MAX];
        for(int i=0; i<nIFO; i++) vDIF[i] = CWB::mdc::GetDiff(&vREC[i], &cINJ[i]);
        double sync_time;
        vDIF[1]*=-1.0;
        double sync_xcor = CWB::mdc::TimeSync(vDIF[0], vDIF[1], sync_time);
        cout << ID << " vDIF[0]xvDIF[1] --------------------------------> sync_time " << sync_time << " sync_xcor : " << sync_xcor << endl;
#endif

#ifdef TEST3
        char title[256];
        char ofname[256];
        for(int i=0; i<nIFO; i++) {
          // PLOT -> vRES 
          sprintf(title,"%s (TIME) : vREC(red) - cINJ(blue)",NET->ifoName[i]);
          sprintf(ofname,"%s_vREC_cINJ%d.%s",NET->ifoName[i],ID,"root");
          PlotWaveform(ofname, title, cfg, &vREC[i], &cINJ[i], NULL, &vREC[i], false);
        }
#endif

        SetOutputFile(NET, EVT, cfg, true);            	// set output root file

        if(gOPT.output_root) {
          DumpOutputFile(NET, EVT, cfg, ID, k, ofactor); 	// dump event to output root file

          for(int n=0;n<nIFO;n++) {
            detector* pD = NET->getifo(n);
            if(!cfg->simulation) ClearWaveforms(pD);     	// release waveform memory
          }
        }
      }
    }

    jfile->cd();
    if(EVT) delete EVT;
  }

  return;
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

void SetOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, bool dump_wf) {

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
       cout << "CWB_Plugin_WF.C : Error - output root file not found" << endl;
       gSystem->Exit(1);                                                                             
     }                                                                                               
   } else {                                                                                          
     cout << "CWB_Plugin_WF.C : Error - output root file not found" << endl;                   
     gSystem->Exit(1);                                                                               
   }                                                                                                 

   if(dump_wf) {
     if(gOPT.output_inj) if(cfg->simulation) for(int n=0;n<nIFO;n++) gWF.wINJ[n] = &vINJ[n];
     if(gOPT.output_inj) if(cfg->simulation && vAUX.size()==nIFO) for(int n=0;n<nIFO;n++) gWF.wAUX[n] = &vAUX[n];
     if(gOPT.output_inj) if(cfg->simulation) for(int n=0;n<nIFO;n++) gWF.cINJ[n] = &cINJ[n];
     if(gOPT.output_rec) for(int n=0;n<nIFO;n++) gWF.wREC[n] = &vREC[n];
     if(gOPT.output_wht) for(int n=0;n<nIFO;n++) gWF.wWHT[n] = &vWHT[n];
     if(gOPT.output_dat) for(int n=0;n<nIFO;n++) gWF.wDAT[n] = &vDAT[n];
     if(gOPT.output_nul) for(int n=0;n<nIFO;n++) gWF.wNUL[n] = &vNUL[n];
     if(gOPT.output_cluster) gWF.pCLUSTER = &gCLUSTER;
   }

   gTREE = (TTree *) froot->Get("waveburst");
   if(gTREE!=NULL) {                         
     EVT = new netevent(gTREE,nIFO);         
     if(dump_wf) {
       if(gOPT.fits_pe!="") gTREE->Branch("sm_pe_of",&gSM_PE_OF,"sm_pe_of/F");
       if(gOPT.fits_os!="") gTREE->Branch("sm_os_of",&gSM_OS_OF,"sm_os_of/F");
       if(cfg->simulation)  gTREE->Branch("wf_ff", gFF, TString::Format("wf_ff[%i]/F",3));
       if(cfg->simulation)  gTREE->Branch("wf_of", gOF, TString::Format("wf_of[%i]/F",3));
       if(cfg->simulation)  gTREE->Branch("wf_re", gRE, TString::Format("wf_re[%i]/F",3));
       if(cfg->simulation)  gTREE->Branch("wf_tcoa", &gTCOA, "wf_tcoa/D");
       for(int n=0;n<nIFO;n++) {
         if(gOPT.output_inj) if(cfg->simulation) gTREE->Branch(TString::Format("wINJ_%d",n).Data(),"wavearray<double>",&gWF.wINJ[n],32000,0);
         if(gOPT.output_inj) if(cfg->simulation && vAUX.size()==nIFO) gTREE->Branch(TString::Format("wAUX_%d",n).Data(),"wavearray<double>",&gWF.wAUX[n],32000,0);
         if(gOPT.output_inj) if(cfg->simulation) gTREE->Branch(TString::Format("cINJ_%d",n).Data(),"wavearray<double>",&gWF.cINJ[n],32000,0);
         if(gOPT.output_rec) gTREE->Branch(TString::Format("wREC_%d",n).Data(),"wavearray<double>",&gWF.wREC[n],32000,0);
         if(gOPT.output_wht) gTREE->Branch(TString::Format("wWHT_%d",n).Data(),"wavearray<double>",&gWF.wWHT[n],32000,0);
         if(gOPT.output_dat) gTREE->Branch(TString::Format("wDAT_%d",n).Data(),"wavearray<double>",&gWF.wDAT[n],32000,0);
         if(gOPT.output_nul) gTREE->Branch(TString::Format("wNUL_%d",n).Data(),"wavearray<double>",&gWF.wNUL[n],32000,0);
         if(gOPT.output_strain)  gTREE->Branch(TString::Format("wSTRAIN_%d",n).Data(),"wavearray<double>",&gSTRAIN[n],32000,0);
         if(gOPT.output_mdc)     gTREE->Branch(TString::Format("wMDC_%d",n).Data(),"wavearray<double>",&gMDC[n],32000,0);
         if(gOPT.output_cstrain) gTREE->Branch(TString::Format("wCSTRAIN_%d",n).Data(),"wavearray<double>",&gCSTRAIN[n],32000,0);
       }
       if(gOPT.output_cluster) {
         gTREE->Branch("seggps",&gSEGGPS,"seggps/D");
         gTREE->Branch("crate",&gCRATE,"crate/F");
         gTREE->Branch("cluster","vector<netpixel>",&gWF.pCLUSTER,32000,0);
       }
       if(gOPT.output_pol) {
         gTREE->Branch("rPOL00","wavearray<double>",&gWF.rPOL00,32000,0);
         gTREE->Branch("aPOL00","wavearray<double>",&gWF.aPOL00,32000,0);
         gTREE->Branch("rPOL90","wavearray<double>",&gWF.rPOL90,32000,0);
         gTREE->Branch("aPOL90","wavearray<double>",&gWF.aPOL90,32000,0);
       }
     }
   } else {                                                                      
     EVT = new netevent(nIFO);                                                   
     gTREE = EVT->setTree();                                                  
     if(gOPT.fits_pe!="") gTREE->SetBranchAddress("sm_pe_of",&gSM_PE_OF);
     if(gOPT.fits_os!="") gTREE->SetBranchAddress("sm_os_of",&gSM_OS_OF);
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

void ReadUserOptions(TString options) {

  // get plugin options 

  if(TString(options)!="") {

    //cout << "WF options : " << options << endl;
    TObjArray* token = TString(options).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++) {

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();

      if(stok.Contains("wf_output_enable=")) {
        TString wf_output_enable=stok;
        wf_output_enable.Remove(0,wf_output_enable.Last('=')+1);
        if(wf_output_enable=="root")    gOPT.output_root=true;
        if(wf_output_enable=="inj")     gOPT.output_inj=true;
        if(wf_output_enable=="rec")     gOPT.output_rec=true;
        if(wf_output_enable=="wht")     gOPT.output_wht=true;
        if(wf_output_enable=="dat")     gOPT.output_dat=true;
        if(wf_output_enable=="nul")     gOPT.output_nul=true;
        if(wf_output_enable=="strain")  gOPT.output_strain=true;
        if(wf_output_enable=="mdc")     gOPT.output_mdc=true;
        if(wf_output_enable=="cstrain") gOPT.output_cstrain=true;
        if(wf_output_enable=="cluster") gOPT.output_cluster=true;
        if(wf_output_enable=="pol")     gOPT.output_pol=true;
      }

      if(stok.Contains("wf_output_disable=")) {
        TString wf_output_disable=stok;
        wf_output_disable.Remove(0,wf_output_disable.Last('=')+1);
        if(wf_output_disable=="root")    gOPT.output_root=false;
        if(wf_output_disable=="inj")     gOPT.output_inj=false;
        if(wf_output_disable=="rec")     gOPT.output_rec=false;
        if(wf_output_disable=="wht")     gOPT.output_wht=false;
        if(wf_output_disable=="dat")     gOPT.output_dat=false;
        if(wf_output_disable=="nul")     gOPT.output_nul=false;
        if(wf_output_disable=="strain")  gOPT.output_strain=false;
        if(wf_output_disable=="mdc")     gOPT.output_mdc=false;
        if(wf_output_disable=="cstrain") gOPT.output_cstrain=false;
        if(wf_output_disable=="cluster") gOPT.output_cluster=false;
        if(wf_output_disable=="pol")     gOPT.output_pol=false;
      }

      if(stok.Contains("wf_fits_pe=")) {
        TString wf_fits_pe=stok;
        wf_fits_pe.Remove(0,wf_fits_pe.Last('=')+1);
        gOPT.fits_pe=wf_fits_pe;
      }
      if(stok.Contains("wf_gps_pe=")) {
        TString wf_gps_pe=stok;
        wf_gps_pe.Remove(0,wf_gps_pe.Last('=')+1);
        if(wf_gps_pe.IsFloat()) gOPT.gps_pe=wf_gps_pe.Atof();
      }

      if(stok.Contains("wf_fits_os=")) {
        TString wf_fits_os=stok;
        wf_fits_os.Remove(0,wf_fits_os.Last('=')+1);
        gOPT.fits_os=wf_fits_os;
      }
      if(stok.Contains("wf_gps_os=")) {
        TString wf_gps_os=stok;
        wf_gps_os.Remove(0,wf_gps_os.Last('=')+1);
        if(wf_gps_os.IsFloat()) gOPT.gps_os=wf_gps_os.Atof();
      }
      if(stok.Contains("wf_inj_tstep=")) {
        TString wf_inj_tstep=stok;
        wf_inj_tstep.Remove(0,wf_inj_tstep.Last('=')+1);
        if(wf_inj_tstep.IsFloat()) gOPT.inj_tstep=wf_inj_tstep.Atof();
      }
      if(stok.Contains("wf_inj_toff=")) {
        TString wf_inj_toff=stok;
        wf_inj_toff.Remove(0,wf_inj_toff.Last('=')+1);
        if(wf_inj_toff.IsFloat()) gOPT.inj_toff=wf_inj_toff.Atof();
      }
    }
  }
}

void ResetUserOptions() {

    gOPT.output_root     = WF_OUTPUT_ROOT;
    gOPT.output_inj      = WF_OUTPUT_INJ;
    gOPT.output_rec      = WF_OUTPUT_REC;
    gOPT.output_wht      = WF_OUTPUT_WHT;
    gOPT.output_dat      = WF_OUTPUT_DAT;
    gOPT.output_nul      = WF_OUTPUT_NUL;
    gOPT.output_strain   = WF_OUTPUT_STRAIN;
    gOPT.output_mdc      = WF_OUTPUT_MDC;
    gOPT.output_cstrain  = WF_OUTPUT_CSTRAIN;
    gOPT.output_cluster  = WF_OUTPUT_CLUSTER;
    gOPT.output_pol      = WF_OUTPUT_POL;
    gOPT.fits_pe         = WF_FITS_PE;
    gOPT.gps_pe          = WF_GPS_PE;
    gOPT.fits_os         = WF_FITS_OS;
    gOPT.gps_os          = WF_GPS_OS;
    gOPT.inj_tstep       = WF_INJ_TSTEP;
    gOPT.inj_toff        = WF_INJ_TOFF;
}

void PrintUserOptions(CWB::config* cfg) {

    cout << "-----------------------------------------"     << endl;
    cout << "WF config options       "     << endl;
    cout << "-----------------------------------------"     << endl << endl;

    cout << "WF_OUTPUT_ROOT       " << gOPT.output_root     << endl;

    cout << "WF_OUTPUT_INJ        " << gOPT.output_inj      << endl;
    cout << "WF_OUTPUT_REC        " << gOPT.output_rec      << endl;
    cout << "WF_OUTPUT_WHT        " << gOPT.output_wht      << endl;
    cout << "WF_OUTPUT_DAT        " << gOPT.output_dat      << endl;
    cout << "WF_OUTPUT_NUL        " << gOPT.output_nul      << endl;

    cout << "WF_OUTPUT_STRAIN     " << gOPT.output_strain   << endl;
    cout << "WF_OUTPUT_MDC        " << gOPT.output_mdc      << endl;
    cout << "WF_OUTPUT_CSTRAIN    " << gOPT.output_cstrain  << endl;

    cout << "WF_OUTPUT_CLUSTER    " << gOPT.output_cluster  << endl;
    cout << "WF_OUTPUT_POL        " << gOPT.output_pol      << endl;

    cout << "WF_FITS_PE           " << gOPT.fits_pe         << endl;
    cout << "WF_GPS_PE            " << gOPT.gps_pe          << endl;
    cout << "WF_FITS_OS           " << gOPT.fits_os         << endl;
    cout << "WF_GPS_OS            " << gOPT.gps_os          << endl;

    cout << "WF_INJ_TSTEP         " << gOPT.inj_tstep       << endl;
    cout << "WF_INJ_TOFF          " << gOPT.inj_toff        << endl;

    cout << endl;
}

std::vector<wavearray<double> > GetWhitenedData(network* NET, CWB::config* cfg) {
  // get whitened data (in the vREC time range)

  std::vector<wavearray<double> > xWHT;        // temporary stuff

  int nIFO = NET->ifoListSize();               // number of detectors

  for(int n=0; n<nIFO; n++) {
    xWHT.push_back(CWB::mdc::GetAligned(&gHOT[n], &vREC[n]));
  }
  return xWHT;
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
           cout << "CWB_Plugin_WF.C : Error : Injected waveform not saved !!! : detector "
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

std::vector<wavearray<double> > GetAuxWaveform(network* NET, netevent* EVT, int ID, double factor) {

  std::vector<wavearray<double> > xAUX;        // auxiliary injected

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

  // search the auxiliary injection, same injection time
  detector* pD = NET->getifo(0);
  int J = pD->IWFP.size();
  cout.precision(14);
  int auxID = -1;
  for(int j=0;j<J;j++) {
    int mdcID = pD->IWFID[j];
    if(mdcID<0 || mdcID==injID) continue;	// skip primary and strain injections
    wavearray<double>* wf = (wavearray<double>*)pD->IWFP[j];
    double start = wf->start();
    double stop = start+wf->size()/wf->rate();
    if(recTime>start && recTime<stop) auxID = j;
  }

  if(auxID>0) {
    for(int n=0; n<nIFO; n++) {
      pD = NET->getifo(n);
      cout.precision(14);
      wavearray<double>* wf = (wavearray<double>*)pD->IWFP[auxID];
      xAUX.push_back(*wf);
    } 
  }

  return xAUX;
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
        cout << "CWB_Plugin_WF.C : Error : Reconstructed waveform not saved !!! : ID -> "
             << ID << " : detector " << NET->ifoName[n] << endl;
        wavearray<double> wf;
        xREC.push_back(wf);
        continue;
     }
     if(wfIndex>=0) pwfREC[n] = pd->RWFP[wfIndex];

     wavearray<double>* wf = pwfREC[n];
     xREC.push_back(*wf);
  }
  delete [] pwfREC;

  return xREC;
}

void ClearVectors() {

   while(!vINJ.empty()) vINJ.pop_back();
   vINJ.clear(); std::vector<wavearray<double> >().swap(vINJ);
   while(!cINJ.empty()) cINJ.pop_back();
   cINJ.clear(); std::vector<wavearray<double> >().swap(cINJ);
   while(!vREC.empty()) vREC.pop_back();
   vREC.clear(); std::vector<wavearray<double> >().swap(vREC);
   while(!vWHT.empty()) vWHT.pop_back();
   vWHT.clear(); std::vector<wavearray<double> >().swap(vWHT);
   while(!vDAT.empty()) vDAT.pop_back();
   vDAT.clear(); std::vector<wavearray<double> >().swap(vDAT);
   while(!vNUL.empty()) vNUL.pop_back();
   vNUL.clear(); std::vector<wavearray<double> >().swap(vNUL);
   while(!vRES.empty()) vRES.pop_back();
   vRES.clear(); std::vector<wavearray<double> >().swap(vRES);
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

  std::vector<netpixel> vPIX = SavePixels(NET, cfg, lag, id);	// save original pixels

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

  RestorePixels(NET, cfg, pwc, &vPIX, id);	     // restore original pixels	

#ifdef PLOT_INJ	
  wavearray<double> wDIF[NIFO_MAX];
  for(int i=0; i<nIFO; i++) {
    char title[256];
    char ofname[256];

    sprintf(title,"%s (TIME) : vREC(red) - vSIG(blue)",NET->ifoName[i]);
    sprintf(ofname,"%s_vREC_vSIG_time_id_%d.%s",NET->ifoName[i],id,"root");
    PlotWaveform(ofname, title, cfg, &vREC[i], &vSIG[i], NULL, &vREC[i], false);   // time

    sprintf(title,"%s (TIME) : vINJ(red) - vSIG(blue)",NET->ifoName[i]);
    sprintf(ofname,"%s_vINJ_vSIG_time_id_%d.%s",NET->ifoName[i],id,"root");
    PlotWaveform(ofname, title, cfg, &vINJ[i], &vSIG[i], NULL, &vREC[i], false);   // time

    // PLOT -> vREC : vREC-vSIG 
    wDIF[i] = CWB::mdc::GetDiff(&vREC[i], &vSIG[i]);
    sprintf(title,"%s (TIME) : vREC(red) - cREC-vSIG(blue)",NET->ifoName[i]);
    sprintf(ofname,"%s_cINJ_vREvSIG_id_%d.%s",NET->ifoName[i],id,"root");
    PlotWaveform(ofname, title, cfg, &vREC[i], &wDIF[i], NULL, &vREC[i], false);
  }
#endif

  if(_AVX.size()) _avx_free_ps(_AVX);
  if(_DAT.size()) _avx_free_ps(_DAT);           // container for data packet amplitudes
  if(_vtd.size()) _avx_free_ps(_vtd);           // array for 00 amplitudes
  if(_vTD.size()) _avx_free_ps(_vTD);           // array for 90 amplitudes

  return vSIG;
}
 
void Wave2Sparse(network* NET, CWB::config* cfg, char type) {

  if(type!='j' && type!='r' && type!='v' && type!='d' && type!='x') {
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
     if(type=='x') xSS[n]=pD->vSS;                                      
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
    if(type=='x') {                                              
      int wos = double(vDAT[n].start()-WF[n].start())*WF[n].rate();   
      for (int i=0;i<vRES[n].size();i++) WF[n][wos+i] = vRES[n][i];   
    }                                                                 

#ifdef SAVE_WAVE2SPARSE_PLOT
    gwavearray<double> gw(&WF[n]);
    gw.Draw(GWAT_TIME);
    watplot* plot = gw.GetWATPLOT();
    TString gfile;
    if(type=='r') gfile=TString::Format("%s/WAVE2SPARSE_REC_%s.root",".",NET->ifoName[n]);
    if(type=='j') gfile=TString::Format("%s/WAVE2SPARSE_INJ_%s.root",".",NET->ifoName[n]);
    if(type=='d') gfile=TString::Format("%s/WAVE2SPARSE_DAT_%s.root",".",NET->ifoName[n]);
    if(type=='x') gfile=TString::Format("%s/WAVE2SPARSE_RES_%s.root",".",NET->ifoName[n]);
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
       if(type=='x') size = xSS[n][k].sparseMap00.size();

       for(int j=0; j<size; j++) {
         int index;               
         if(type=='r') index = rSS[n][k].sparseIndex[j];
         if(type=='j') index = jSS[n][k].sparseIndex[j];
         if(type=='d') index = dSS[n][k].sparseIndex[j];
         if(type=='x') index = xSS[n][k].sparseIndex[j];
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
         if(type=='x') {                                
           xSS[n][k].sparseMap00[j]=v00;                     
           xSS[n][k].sparseMap90[j]=v90;                     
         }                                                   
       }
    }
  }
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
  CWB::mdc::GetTimeBoundaries(*wref, P, bT, eT);

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
    if(pdir!="") ofname = TString(pdir)+TString("/")+ofname;
    PTS.canvas->Print(ofname);
    cout << "write : " << ofname << endl;
  }
}

void 
DumpSkymap(network* NET, int lag, netevent* EVT, int ID) {

  // Dump2fits probability skymap  (healpix)
  skymap skyprobcc = NET->getifo(0)->tau;
  skyprobcc=0.;
  skymap skyprob = skyprobcc;
  skyprob=1.e-12;

  std::vector<float>* vP;
  std::vector<int>*   vI;

  vP = &(NET->wc_List[lag].p_Map[ID-1]);
  vI = &(NET->wc_List[lag].p_Ind[ID-1]);
  double th,ph,ra;
  int k;
  for(int j=0; j<int(vP->size()); j++) {
    int i = (*vI)[j];
    th = skyprob.getTheta(i);
    ph = skyprob.getPhi(i);

    k=skyprob.getSkyIndex(th, ph);
    skyprob.set(k,(*vP)[j]);

    ra = skyprob.getRA(i);
    k=skyprob.getSkyIndex(th, ra);
    skyprobcc.set(k,(*vP)[j]);
  }

  char fname[1024];
  sprintf(fname, "skyprobcc_%d.%s", ID, "fits");
  skyprobcc.Dump2fits(const_cast<char*>(fname),EVT->time[0],const_cast<char*>(""),const_cast<char*>("PROB"),const_cast<char*>("pix-1"),'C');

  return;
} 

skymap 
GetSkyMap(network* NET, int lag, netevent* EVT, int ID) {

  skymap skyprobcc = NET->getifo(0)->tau;
  skyprobcc=0.;
  skymap skyprob = skyprobcc;
  skyprob=1.e-12;

  std::vector<float>* vP;
  std::vector<int>*   vI;

  vP = &(NET->wc_List[lag].p_Map[ID-1]);
  vI = &(NET->wc_List[lag].p_Ind[ID-1]);
  double th,ph,ra;
  int k;
  for(int j=0; j<int(vP->size()); j++) {
    int i = (*vI)[j];
    th = skyprob.getTheta(i);
    ph = skyprob.getPhi(i);

    k=skyprob.getSkyIndex(th, ph);
    skyprob.set(k,(*vP)[j]);

    ra = skyprob.getRA(i);
    k=skyprob.getSkyIndex(th, ra);
    skyprobcc.set(k,(*vP)[j]);
  }

  return skyprobcc;
} 

double
ComputeMaxSkyMapProduct(skymap* sm1, skymap* sm2) {

  skymap csm = *sm2;			// cumulative probability

  double integral=0;
  int L=csm.size();
  int* index = new int[L];              // sorted index
  double* prob = new double[L];         // probability array
  for(int l=0;l<L;l++) {                // loop over the sky grid
    prob[l] = csm.get(l);               // get skymap prob
    integral += prob[l];
  }

  TMath::Sort(L,prob,index,false);      // sort prob
  double cumul=0.0;
  for(int l=0;l<L;l++) {                // loop over the sky grid
    int m=index[l];                     // sorted index
    cumul+=prob[m];
    csm.set(m,cumul/integral);          // set norlmalized cumulative probabilty
  }

  double max=-1.e+20;
  for(int l=0;l<L;l++) {       		// loop over the sky grid
    double p1=sm1->get(l);
    double p2=csm.get(l);
    if(p1*p2>max) max=p1*p2;
  }

  return max;
}

double
ComputeSkyMapOF(skymap* sm1, skymap* sm2) {

  int L = sm1->size();         // number of pixels in the sphere

  double pp12=0.;
  double pp11=0.;
  double pp22=0.;
  for(int l=0;l<L;l++) {       // loop over the sky grid
    double p1=sm1->get(l);
    double p2=sm2->get(l);
    pp12 += p1*p2;
    pp11 += p1*p1;
    pp22 += p2*p2;
  }

  double OF = pp12/sqrt(pp11*pp22);

  return OF;
}

skymap
ReadSkyMap(TString fitsFile, CWB::config* cfg, double gps) {

  gskymap sm((char*)fitsFile.Data());
  sm.resample(cfg->healpix);

  return EarthCoordinatesSM(&sm, gps);
}

skymap 
EarthCoordinatesSM(skymap* sm, double gps) {

  gskymap xsm = *sm;

  int L = sm->size();                	// number of pixels in the sphere
  for(int l=0;l<L;l++) {                // loop over the sky grid
    double p = sm->get(l);          
    double phi = sm->getPhi(l);
    double theta = sm->getTheta(l);
    double xphi = sm->RA2phi(phi, gps);
    int k = sm->getSkyIndex(theta, xphi);
    xsm.set(k,p);            
  }

  return xsm;
}

std::vector<netpixel> 
SavePixels(network* NET, CWB::config* cfg, int lag, int id) {

  netpixel* pix;
  netcluster* pwc = NET->getwc(lag);
  std::vector<netpixel> vPIX;

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, id, false);   // buffer for pixel IDs
  int V = pI.size();
  if(V>cfg->BATCH) return vPIX;                                 // attach TD amp to pixels < V

  for(int j=0; j<V; j++) {                    // loop over principal components
    pix = pwc->getPixel(id,pI[j]);
    vPIX.push_back(*pix);                     // save original pixels
  }

  return vPIX;
}

void 
RestorePixels(network* NET, CWB::config* cfg, netcluster* pwc, std::vector<netpixel>* vPIX, int id) {

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, id, false);   // buffer for pixel IDs
  int V = pI.size();
  if(V>cfg->BATCH) return;                                      // attach TD amp to pixels < V
  for(int j=0; j<V; j++) {
    netpixel* pix = pwc->getPixel(id,pI[j]);
    *pix = (*vPIX)[j];
  }

  while(!vPIX->empty()) vPIX->pop_back();
  vPIX->clear(); std::vector<netpixel>().swap(*vPIX);
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

std::vector<netpixel> GetCluster(network* NET, int lag, int id) {

  netpixel* pix;
  netcluster* pwc = NET->getwc(lag);
  std::vector<netpixel> vPIX;

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, id, false);   // buffer for pixel IDs

  int V = pI.size();
  for(int j=0; j<V; j++) {                    // loop over principal components
    pix = pwc->getPixel(id,pI[j]);
    vPIX.push_back(*pix);                     // save original pixels
  }

  return vPIX;
}

