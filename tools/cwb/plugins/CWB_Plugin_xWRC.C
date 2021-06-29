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
#include "TF2.h"
#include "TMath.h"
#include "mdc.hh"
#include "watplot.hh"
#include "gwavearray.hh"
#include "gskymap.hh"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

//#define SAVE_WHT_PLOT			// enable event WHITE plots 
//#define SAVE_WF_PLOT
//#define SAVE_SKYPROB
//#define SAVE_PLOT_WERR

//#define SAVE_PLOT_WERR2
//#define SAVE_PLOT_RESIDUALS
//#define SAVE_MRA_PLOT			// enable event MRA plots 
//#define SAVE_WF_COMPARISON

#define PLOT_DIR	"plot"

// ----------------------------------------------------
// WRC plugin functions
// ----------------------------------------------------

void ComputeResidualEnergy(wavearray<double>* wfINJ, wavearray<double>* wfREC, 
                           double& enINJ, double& enREC, double& xcorINJ_REC);

void SetOutput(network* NET, netevent* &EVT, CWB::config* cfg);
void DumpOutput(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor);
void GetRecWaveform(network* NET, netevent* EVT, int ID, int k, int factor);
void GetInjWaveform(network* NET, netevent* EVT, int ID, int k, int factor);
void GetRInjWaveform(network* NET, netevent* EVT, CWB::config* cfg, int ID, int k);
void ReadUserOptions();
void Clear();
void SetSkyMask(network* NET, CWB::config* cfg, int seed=0);
void AddNoise2Sparse(network* NET, CWB::config* cfg, int seed=0);
void Wave2Sparse(network* NET, CWB::config* cfg, char wave_type);
void SaveSkyProb(network* NET, CWB::config* cfg, int id);
void ComputeErrorWF(wavearray<double>* wfINJ, wavearray<double>* wfREC, int ifoID);
void ComputeErrorStatistic(network* NET, CWB::config* cfg, int ID);
double GetSparseMap(SSeries<double>* SS, bool phase, int index);
int _sse_mra_ps(network* NET, float* amp, float* AMP, float Eo, int K);

void PlotFinal(int n); 
void PlotFinal2(network* NET, int ID);
void PlotFinal3(network* NET, int ID, int ifoID, wavearray<double>* w1, 
                wavearray<double>* w2, TString tag, TString gtype="png");
void PlotWaveform(TString ifo, wavearray<double>* wfINJ, wavearray<double>* wfREC,
                  CWB::config* cfg, bool fft=false, bool strain=false);
void PlotMRA(network* NET, int ID, int k, TString tag);

// ----------------------------------------------------
// Temporary WRC structures
// ----------------------------------------------------

wavearray<double> wINJ[NIFO_MAX];		// injected
wavearray<double> sINJ[NIFO_MAX];		// injected aligned to the same time
wavearray<double> rINJ[NIFO_MAX];		// injected waveform in wREC TF sub-domain
wavearray<double> wREC[NIFO_MAX];		// reconstructed
double wDINJ[NIFO_MAX];				// inj delays
double wDREC[NIFO_MAX];				// inj delays

std::vector<SSeries<double> > vSS[NIFO_MAX];	// original sparse maps
std::vector<SSeries<double> > sSS[NIFO_MAX];	// injected aligned to the same time
std::vector<SSeries<double> > rSS[NIFO_MAX];	// reconstructed
std::vector<SSeries<double> > jSS[NIFO_MAX];	// injected
std::vector<WSeries<double> > vN[NIFO_MAX];	// sim noise

std::vector<int>    wTRY[NIFO_MAX];
std::vector<double> wAVR[NIFO_MAX];
std::vector<double> wRMS[NIFO_MAX];

std::vector<wavearray<double> > vREC[NIFO_MAX];	// reconstructed
std::vector<double> vDREC[NIFO_MAX];		// detector rec waveform delays

wavearray<int> sI;				// shuffle index array

cwb     xCWB;
TTree*  net_tree = NULL;
TString outDump;
bool    detected=false;

skymap rec_skyprob;      		// save reconstructed probability skymap
TH1D   hist_skyprob;			// used to extract random skyloc for skymask

double inj_phi;				// save injected phi
double inj_theta;			// save injected theta

double rec_phi;				// save reconstructed phi
double rec_theta;			// save reconstructed theta

float rec_erA[11];			// save reconstructed error regions

// ----------------------------------------------------
// Config Parameters
// ----------------------------------------------------

//#define USE_wINJ	// use wINJ (default wREC) (are pixels used as inj+noise in the retry)
//#define USE_wWHT	// use wWHT (default wREC)

#define SKYMASK_SEED	0

#define APPLY_INJ_MRA


// user config options : default values or read from command line
bool WRC_SKYMASK_CONFIG=true;
bool WRC_NOISE_CONFIG=true;
int  WRC_RETRY_CONFIG=100;
bool WRC_AMP_CAL_ERR_CONFIG=false;
bool WRC_PHS_CAL_ERR_CONFIG=false;
int  WRC_CED_ID_CONFIG=0;
int  WRC_CED_RETRY_CONFIG=0;	// dump CED at gLRETRY=WRC_CED_RETRY
int  WRC_ID_CONFIG=0;
bool WRC_SKYMASK_REC_SKY_CONFIG=true;

// config options used in the plugin
bool WRC_SKYMASK=false;
bool WRC_NOISE=false;
int  WRC_RETRY=0;
bool WRC_AMP_CAL_ERR=false;
bool WRC_PHS_CAL_ERR=false;
int  WRC_CED_ID=0;
int  WRC_CED_RETRY=0;		// dump CED at gLRETRY=WRC_CED_RETRY
int  WRC_ID=0;
bool WRC_SKYMASK_REC_SKY=false;

#define AMP_CAL_ERR	0.1
#define PHS_CAL_ERR	10

//#define COMPUTE_RESIDUAL_ENERGY

// ----------------------------------------------------
// Output Parameters
// ----------------------------------------------------

float erR[11];			// probability distribution of residuals

#define nPAR    7

TString parName[nPAR] = { "wrc_try",
                          "wrc_wf_noise",
                          "wrc_sm_phi",
                          "wrc_sm_theta",
                          "wrc_sm_radius",
                          "wrc_cal_amp",
                          "wrc_cal_phs"
                        };
float   parValue[nPAR];

float  WRC_WF_NOISE=0;
float  WRC_SM_THETA;
float  WRC_SM_PHI;
float  WRC_SM_RADIUS=0.1;
float  WRC_CAL_AMP=0;
float  WRC_CAL_PHS=0;

#define WRC_PLUGIN_VERSION	"v1.3"

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!NOISE_MDC_SIMULATION
// This plugin is for Waveform Reconstruction Challenge (WRC) studies 
// using the last Likelihood Stage Event Loop functionality 

  cout << endl;
  cout << "-----> CWB_Plugin_xWRC.C - " << WRC_PLUGIN_VERSION << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  

    cfg->outPlugin=true;  				// disable built-in output root file

    // set default config options
    WRC_RETRY           = WRC_RETRY_CONFIG;
    WRC_SKYMASK         = WRC_SKYMASK_CONFIG;
    WRC_NOISE           = WRC_NOISE_CONFIG;
    WRC_AMP_CAL_ERR     = WRC_AMP_CAL_ERR_CONFIG;
    WRC_PHS_CAL_ERR     = WRC_PHS_CAL_ERR_CONFIG;
    WRC_CED_ID          = WRC_CED_ID_CONFIG;
    WRC_CED_RETRY       = WRC_CED_RETRY_CONFIG;
    WRC_SKYMASK_REC_SKY = WRC_SKYMASK_REC_SKY_CONFIG;

    // user config options : read from command line
    ReadUserOptions(); 

    // print config options
    cout << "-----------------------------------------" << endl; 
    cout << "WRC config options                       " << endl << endl;
    cout << "WRC_RETRY             " << WRC_RETRY	<< endl;
    cout << "WRC_SKYMASK           " << WRC_SKYMASK	<< endl;
    cout << "WRC_NOISE             " << WRC_NOISE	<< endl;
    cout << "WRC_AMP_CAL_ERR       " << WRC_AMP_CAL_ERR	<< endl;
    cout << "WRC_PHS_CAL_ERR       " << WRC_PHS_CAL_ERR	<< endl;
    cout << "WRC_CED_ID            " << WRC_CED_ID	<< endl;
    cout << "WRC_CED_RETRY         " << WRC_CED_RETRY	<< endl;
    cout << "WRC_SKYMASK_REC_SKY   " << WRC_SKYMASK_REC_SKY << endl;
    cout << endl;
  }

  if(type==CWB_PLUGIN_NETWORK) {
  }

  if(type==CWB_PLUGIN_ILIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_xWRC.C -> CWB_PLUGIN_ILIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    if(NET->nLag!=1) {
      cout << "CWB_Plugin_xWRC.C -> implemented only for zero lag" << endl;
      gSystem->Exit(1);
    }

    // export gLRETRY
    EXPORT(int,gLRETRY,TString::Format("gLRETRY = %d;",WRC_RETRY).Data())
  }

  if(type==CWB_PLUGIN_XLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_xWRC.C -> CWB_PLUGIN_XLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // import gLRETRY
    int gLRETRY=-1; IMPORT(int,gLRETRY)
    cout << "-----> CWB_Plugin_xWRC.C -> " << " gLRETRY : " << gLRETRY << "/" << WRC_RETRY << " - WRC : ";

    if(WRC_NOISE)    cout << "N ";
    if(WRC_AMP_CAL_ERR) cout << "A ";
    if(WRC_PHS_CAL_ERR) cout << "P ";
    if(WRC_SKYMASK)     cout << "S ";
    cout << " - ";
#ifdef USE_wINJ	// use wINJ (default wREC) (are pixels used as inj+noise in the retry)
    cout << "wJ ";
#endif
#ifdef USE_wWHT	// use wWHT (default wREC)
    cout << "wW ";
#endif
#if !defined (USE_wINJ) && !defined (USE_wWHT)
    cout << "wR ";
#endif
    cout << endl;
    if((WRC_NOISE||WRC_AMP_CAL_ERR||WRC_PHS_CAL_ERR) && gLRETRY<WRC_RETRY) {
       if(rSS[0].size()>0) AddNoise2Sparse(NET,cfg, 0);
    }
    if(WRC_SKYMASK && gLRETRY<WRC_RETRY) SetSkyMask(NET, cfg, SKYMASK_SEED);	// set new skymask

    if(gLRETRY==WRC_RETRY) {	// first time
       cfg->nSky=-999; NET->nSky=cfg->nSky;
    } else {
       cfg->nSky=1000; NET->nSky=cfg->nSky;
    }

    if(WRC_CED_RETRY) {
      if(gLRETRY==WRC_CED_RETRY) {
        cfg->cedDump=true;
      } else {
        cfg->cedDump=false;
      }
    }

    // for each neew event the wrc config is initialized with the user input wrc values 
    if(gLRETRY==WRC_RETRY) {	// first time
       WRC_SKYMASK         = WRC_SKYMASK_CONFIG;
       WRC_NOISE           = WRC_NOISE_CONFIG;
       WRC_AMP_CAL_ERR     = WRC_AMP_CAL_ERR_CONFIG;
       WRC_PHS_CAL_ERR     = WRC_PHS_CAL_ERR_CONFIG;
       WRC_CED_ID          = WRC_CED_ID_CONFIG;
       WRC_CED_RETRY       = WRC_CED_RETRY_CONFIG;
       WRC_SKYMASK_REC_SKY = WRC_SKYMASK_REC_SKY_CONFIG;

       // restore full skymask
       sprintf(cfg->skyMaskFile,"--theta %f --phi %f --radius %f",0.,0.,360.);
       //cout << cfg->skyMaskFile << endl;
       xCWB.SetSkyMask(NET,cfg,cfg->skyMaskFile,'e');
    }
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_xWRC.C -> CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // import retry
    int gLRETRY=-1; IMPORT(int,gLRETRY)

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];                 
    double ofactor=0;
    if(cfg->simulation==4)      ofactor=-factor;
    else if(cfg->simulation==3) ofactor=-gIFACTOR;
    else                        ofactor=factor;

    int nIFO = NET->ifoListSize();			// number of detectors
    int rate = 0;					// select all resolutions
    int zlag  = 0;					// zelect only zero lag
    netevent* EVT;

    SetOutput(NET, EVT, cfg);				// set output root file

    wavearray<double>id = NET->getwc(zlag)->get(const_cast<char*>("ID"), 0, 'L', rate);

    for(int j=0; j<(int)id.size(); j++) {  		// loop over cluster index

      int ID = size_t(id.data[j]+0.5);

      if(NET->getwc(zlag)->sCuts[ID-1]!=-1) continue;   // skip rejected/processed clusters
 
      if(WRC_ID>0 && ID!=WRC_ID) {			// skip wrc analysis if ID!=WRC_ID
        EXPORT(int,gLRETRY,"gLRETRY = 0;")
        NET->getwc(zlag)->sCuts[ID-1]=1; 		// mark as processed                              
        continue;
      }

      if(WRC_CED_ID) {					// enable/disable CED
        if(ID==WRC_CED_ID) {
          cfg->cedDump=true;
        } else {
          cfg->cedDump=false;
        }
      }

      if(gLRETRY==WRC_RETRY) {				// mark event as detected
         detected=true;
         cout << "-----> Event Detected ID = " << ID << endl;
      }

      GetRecWaveform(NET, EVT, ID, zlag, ofactor);	// get & save reconstructed waveform

#if defined (SAVE_PLOT_WERR) || defined (SAVE_PLOT_WERR2) || (SAVE_WHT_PLOT)
      if(gLRETRY<WRC_RETRY) {	
         for(int n=0; n<nIFO; n++) {
           wavearray<double>* wfINJ = &wINJ[n];			 	// injected waveform
           wavearray<double>* wfREC = &vREC[n][vREC[n].size()-1];	// reconstructed waveform
           ComputeErrorWF(wfINJ, wfREC, n);
         }
#ifdef SAVE_WHT_PLOT
         PlotWaveform(NET->ifoName[n], wfINJ, wfREC, cfg, false, false);	// time	
         PlotWaveform(NET->ifoName[n], wfINJ, wfREC, cfg, true, false);		// fft
#endif
      }

      if(gLRETRY==0 && detected) {	// last time
         // final computation of the waveform errors
         for(int n=0;n<nIFO;n++) {
           for(int j=0;j<wAVR[n].size();j++) {
              if(wTRY[n][j]==0) continue;
              wAVR[n][j]/=wTRY[n][j];
	      wRMS[n][j]/=wTRY[n][j];
    	      wRMS[n][j]-=wAVR[n][j]*wAVR[n][j];
	      wRMS[n][j] =sqrt(wRMS[n][j]);
              //cout << n << " " << j << " wINJ : " << wINJ[n][j] << " wAVR : " << wAVR[n][j] 
              //                      << " wRMS : " << wRMS[n][j] << endl;
           }
#ifdef SAVE_PLOT_WERR
           PlotFinal(n); 
#endif
         } 
#ifdef SAVE_PLOT_WERR2
         PlotFinal2(NET, ID);
#endif
      }
#endif

      if(gLRETRY==0 && detected) {			// last time
#ifdef SAVE_MRA_PLOT					// monster event display
         PlotMRA(NET, ID, zlag, "last_wREC");
#endif
         if(WRC_RETRY>0) ComputeErrorStatistic(NET, cfg, ID); 	// compute final waveform error statistic
         DumpOutput(NET, EVT, cfg, ID, zlag, ofactor);	// dump event to output root file
      }

      // warning : must be done after DumpOutput (GetRInjWaveform destroy original pwc) 
      if(gLRETRY==WRC_RETRY) {				// first time
#ifdef SAVE_MRA_PLOT					// monster event display
         PlotMRA(NET, ID, zlag, "wREC");
#endif

         GetInjWaveform(NET, EVT, ID, zlag, factor);	// get & save injected waveform

         Wave2Sparse(NET,cfg,'v');			// save original sparse map to vSS
         Wave2Sparse(NET,cfg,'r');			// rec
         Wave2Sparse(NET,cfg,'j');			// inj
         Wave2Sparse(NET,cfg,'s');			// wht	
         Wave2Sparse(NET,cfg,'n');			// noise
         SaveSkyProb(NET,cfg,ID);

         GetRInjWaveform(NET, EVT, cfg, ID, zlag);	// get injected waveform in wREC TF sub-domain
      }

      // in the last loop the original event is reprocessed 
      // because must be saved on the final root file
      if(gLRETRY==1 && detected) {	
         WRC_SKYMASK=true;		
         WRC_NOISE=false;
         WRC_AMP_CAL_ERR=false;
         WRC_PHS_CAL_ERR=false;
         WRC_SKYMASK_REC_SKY=true; 

         // restore original sparse maps
         for(int n=0;n<nIFO;n++) {
           detector* pD = NET->getifo(n);
           pD->vSS=vSS[n];
         }
      }

      if(gLRETRY==0 && detected) {			// last time

         // restore original sparse maps
         for(int n=0;n<nIFO;n++) {
           detector* pD = NET->getifo(n);
           pD->vSS=vSS[n];
         }

         // clean all temporary WRC structures
         Clear();				
         detected=false;
      }
    }

    jfile->cd();

    if(EVT)   delete EVT;
  }

  return;
}

void ComputeResidualEnergy(wavearray<double>* wfINJ, wavearray<double>* wfREC, 
                           double& enINJ, double& enREC, double& xcorINJ_REC) {

  double bINJ = wfINJ->start();
  double eINJ = wfINJ->stop();
  double bREC = wfREC->start();
  double eREC = wfREC->stop();
  //cout.precision(14);
  //cout << "bINJ : " << bINJ << " eINJ : " << eINJ << endl;
  //cout << "bREC : " << bREC << " eREC : " << eREC << endl;

  double R=wfINJ->rate();

  int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
  int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);
  //cout << "oINJ : " << oINJ << " oREC : " << oREC << endl;

  double startXCOR = bINJ>bREC ? bINJ : bREC;
  double endXCOR   = eINJ<eREC ? eINJ : eREC;
  int sizeXCOR     = int((endXCOR-startXCOR)*R+0.5);
  //cout << "startXCOR : " << startXCOR << " endXCOR : " 
  //     << endXCOR << " sizeXCOR :" << sizeXCOR << endl;

  enINJ=0;
  enREC=0;
  xcorINJ_REC=0;
  if (sizeXCOR<=0) return;

  // the enINJ, enREC, xcorINJ_REC are computed in the INJ range
  for (int i=0;i<wfINJ->size();i++) enINJ+=wfINJ->data[i]*wfINJ->data[i];
  for (int i=0;i<sizeXCOR;i++) enREC+=wfREC->data[i+oREC]*wfREC->data[i+oREC];
  for (int i=0;i<sizeXCOR;i++) xcorINJ_REC+=wfINJ->data[i+oINJ]*wfREC->data[i+oREC];

  double erINJ_REC = enINJ+enREC-2*xcorINJ_REC;
  //cout << "enINJ : " << enINJ << " enREC : " << enREC 
  //     << " xcorINJ_REC : " << xcorINJ_REC << " enINJREC : " << enINJ+enREC-2*xcorINJ_REC << endl;
  //cout << "erINJ_REC/enINJ : " << erINJ_REC/enINJ << endl;

  return;
}

void ComputeErrorWF(wavearray<double>* wfINJ, wavearray<double>* wfREC, int ifoID) {

  double bINJ = wfINJ->start();
  double eINJ = wfINJ->stop();
  double bREC = wfREC->start();
  double eREC = wfREC->stop();
  //cout.precision(14);
  //cout << "bINJ : " << bINJ << " eINJ : " << eINJ << endl;
  //cout << "bREC : " << bREC << " eREC : " << eREC << endl;

  double R=wfINJ->rate();

  int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
  int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);
  //cout << "oINJ : " << oINJ << " oREC : " << oREC << endl;

  double startXCOR = bINJ>bREC ? bINJ : bREC;
  double endXCOR   = eINJ<eREC ? eINJ : eREC;
  int sizeXCOR     = int((endXCOR-startXCOR)*R+0.5);
  //cout << "startXCOR : " << startXCOR << " endXCOR : " 
  //     << endXCOR << " sizeXCOR :" << sizeXCOR << endl;

  for (int i=0;i<sizeXCOR;i++) {
     wTRY[ifoID][i+oINJ]+=1;
     wAVR[ifoID][i+oINJ]+=wfREC->data[i+oREC];
     wRMS[ifoID][i+oINJ]+=wfREC->data[i+oREC]*wfREC->data[i+oREC];
  }

  return;
}

void PlotMRA(network* NET, int ID, int k, TString tag) {

   if(NET->getwc(k)->sCuts[ID-1]!=-1) return; 	// no events 

   // set sCuts=1 for the events which must be not copied with cps to pwc
   std::vector<int> sCuts = NET->getwc(k)->sCuts; 	// save sCuts
   for(int i=0; i<(int)sCuts.size(); i++) if(i!=ID-1) NET->getwc(k)->sCuts[i]=1;  
   // after cpf, pwc contains only one event
   // ID can not be used to get the event, to get event use ID=1 (ex: for watplot)
   netcluster* pwc = new netcluster();  
   pwc->cpf(*(NET->getwc(k)));				// note: likelihood2G delete tdAmp 

   int nIFO = NET->ifoListSize();                       // number of detectors
   watplot WTS(const_cast<char*>("wts"));
   WTS.canvas->cd();
   char fname[1024];
   sprintf(fname, "%s/%s_l_tfmap_scalogram_%d.png",PLOT_DIR,tag.Data(),ID);
   cout << "write " << fname << endl;
   WTS.plot(pwc, 1, nIFO, 'L', 0, const_cast<char*>("COLZ"));
   WTS.canvas->Print(fname);
   WTS.clear();
   sprintf(fname, "%s/%s_n_tfmap_scalogram_%d.png",PLOT_DIR,tag.Data(),ID);
   cout << "write " << fname << endl;
   WTS.plot(pwc, 1, nIFO, 'N', 0, const_cast<char*>("COLZ"));
   WTS.canvas->Print(fname);
   WTS.clear();

   if(pwc) delete pwc;
   NET->getwc(k)->sCuts = sCuts; 			// restore sCuts
}

void PlotWaveform(TString ifo, wavearray<double>* wfINJ, wavearray<double>* wfREC, 
                  CWB::config* cfg, bool fft, bool strain) {

  watplot PTS(const_cast<char*>("ptswrc"),200,20,800,500);

  //cout << "Print " << fname << endl;
  double tmin = wfINJ->start()<wfREC->start() ? wfINJ->start() : wfREC->start();
  wfINJ->start(wfINJ->start()-tmin);
  wfREC->start(wfREC->start()-tmin);
  if(fft) {
    PTS.plot(wfINJ, const_cast<char*>("ALP"), 1, 0, 0, true, cfg->fLow, cfg->fHigh);
  } else {
    PTS.plot(wfINJ, const_cast<char*>("ALP"), 1, 0, 0);
  }
  PTS.graph[0]->SetLineWidth(1);
  if(fft) {
    PTS.plot(wfREC, const_cast<char*>("SAME"), 2, 0, 0, true, cfg->fLow, cfg->fHigh);
  } else {
    PTS.plot(wfREC, const_cast<char*>("SAME"), 2, 0, 0);
  }
  PTS.graph[1]->SetLineWidth(2);
  wfINJ->start(wfINJ->start()+tmin);
  wfREC->start(wfREC->start()+tmin);

  char label[64]="";
  if(fft) sprintf(label,"%s","fft");
  else    sprintf(label,"%s","time");
  if(strain) sprintf(label,"%s_%s",label,"strain");
  else       sprintf(label,"%s_%s",label,"white"); 

  char fname[1024];
  sprintf(fname, "%s/%s_wf_%s_inj_rec_gps_%d.root",PLOT_DIR,ifo.Data(),label,int(tmin));
  PTS.canvas->Print(fname); 
  cout << "write : " << fname << endl;
  //PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
}

void SetSkyMask(network* NET, CWB::config* cfg, int seed) {

  // NEW SKYMAP

  gRandom->SetSeed(seed);
  int isky = (int)hist_skyprob.GetRandom();
  WRC_SM_THETA = 90-rec_skyprob.getTheta(isky);
  WRC_SM_PHI  = rec_skyprob.getPhi(isky);

  if(WRC_SKYMASK_REC_SKY) {
    WRC_SM_THETA = 90-rec_theta;
    WRC_SM_PHI   = rec_phi;
  }

/*
  // test code : use injected sky location for skymask 
  WRC_SM_THETA = 90-inj_theta;
  WRC_SM_PHI   = inj_phi;
*/

  // --------------------------------------------------------
  // define SetSkyMask
  // --------------------------------------------------------
  sprintf(cfg->skyMaskFile,"--theta %f --phi %f --radius %f",WRC_SM_THETA,WRC_SM_PHI,WRC_SM_RADIUS);
  cout << cfg->skyMaskFile << endl;
  xCWB.SetSkyMask(NET,cfg,cfg->skyMaskFile,'e');

}

void Wave2Sparse(network* NET, CWB::config* cfg, char wave_type) {

   if(wave_type!='j' && wave_type!='r' && wave_type!='s' && wave_type!='n' && wave_type!='v') {
     cout << "Wave2Sparse - Error : wrong wave type" << endl; exit(1);
   }

   gRandom->SetSeed(0);

   //cout << "-----------------> Wave2Sparse" << endl; 

   // init waveform sparse maps
   int nIFO = NET->ifoListSize();			// number of detectors
   for(int n=0;n<nIFO;n++) {
      detector* pD = NET->getifo(n);
      if(wave_type=='v') vSS[n]=pD->vSS;
      if(wave_type=='r') rSS[n]=pD->vSS;
      if(wave_type=='j') jSS[n]=pD->vSS;
      if(wave_type=='s') sSS[n]=pD->vSS;
   }
   if(wave_type=='v') return;

   int nRES = cfg->l_high-cfg->l_low+1;     // number of frequency resolution levels

   // build waveform array
   wavearray<double> WF[NIFO_MAX];
   for(int n=0; n<nIFO; n++) {

     if(wave_type=='s') {
       WF[n].resize(sSS[n][0].wdm_nSTS);
       WF[n].rate(sSS[n][0].wdm_rate);
       WF[n].start(sSS[n][0].wdm_start);
       WF[n]=0;
       int wos = double(sINJ[n].start()-WF[n].start())*WF[n].rate();   
       for (int i=0;i<sINJ[n].size();i++) WF[n][wos+i] = sINJ[n][i];
     } 
     if(wave_type=='j') {
       WF[n].resize(jSS[n][0].wdm_nSTS);
       WF[n].rate(jSS[n][0].wdm_rate);
       WF[n].start(jSS[n][0].wdm_start);
       WF[n]=0;
       int wos = double(wINJ[n].start()-WF[n].start())*WF[n].rate();   
       for (int i=0;i<wINJ[n].size();i++) WF[n][wos+i] = wINJ[n][i];
     } 
     if(wave_type=='r') { 
       WF[n].resize(rSS[n][0].wdm_nSTS);
       WF[n].rate(rSS[n][0].wdm_rate);
       WF[n].start(rSS[n][0].wdm_start);
       WF[n]=0;
       int wos = double(wREC[n].start()-WF[n].start())*WF[n].rate();   
       for (int i=0;i<wREC[n].size();i++) WF[n][wos+i] = wREC[n][i];
     }
     if(wave_type=='n') { 	// noise
       detector* pD = NET->getifo(n);
       WF[n].resize(pD->vSS[0].wdm_nSTS);
       WF[n].rate(pD->vSS[0].wdm_rate);
       WF[n].start(pD->vSS[0].wdm_start);
       WF[n]=0;
       for (int i=0;i<WF[n].size();i++) WF[n][i] = gRandom->Gaus(0,1);
     }

#ifdef SAVE_WF_PLOT
     gwavearray<double> gw(&WF[n]);
     gw.Draw(GWAT_TIME); 
     watplot* plot = gw.GetWATPLOT();
     TString gfile;
     if(wave_type=='r') gfile=TString::Format("%s/REC_%s.png",PLOT_DIR,NET->ifoName[n]);
     if(wave_type=='j') gfile=TString::Format("%s/INJ_%s.png",PLOT_DIR,NET->ifoName[n]);
     if(wave_type=='s') gfile=TString::Format("%s/sINJ_%s.png",PLOT_DIR,NET->ifoName[n]);
     if(wave_type!='n') (*plot) >> gfile;
#endif
   }

   // fill waveform sparse maps
   WSeries<double> w;
   WDM<double>* pwdm = NULL;
   double r00,r90;
   for(int n=0; n<nIFO; n++) {
     for(int i=0; i<nRES; i++) {
        int level = i+cfg->l_low; 
        int layers = level>0 ? 1<<level : 0;
        pwdm = NET->getwdm(layers+1);
        w.Forward(WF[n],*pwdm);

        if((wave_type=='n')&&(i==nRES-1)) {		// time shuffled index for noise array vN
          sI.resize(w.sizeZero()); 
          for(int j=0;j<sI.size();j++) sI[j]=j; 
        }

        int size;
        if(wave_type=='r') size = rSS[n][i].sparseMap00.size();
        if(wave_type=='j') size = jSS[n][i].sparseMap00.size();
        if(wave_type=='s') size = sSS[n][i].sparseMap00.size();
        if(wave_type=='n') vN[n].push_back(w);

        for(int j=0; j<size; j++) {
          int index;
          if(wave_type=='r') index = rSS[n][i].sparseIndex[j];
          if(wave_type=='j') index = jSS[n][i].sparseIndex[j];
          if(wave_type=='s') index = sSS[n][i].sparseIndex[j];
          double v00 = w.pWavelet->pWWS[index];
          double v90 = w.pWavelet->pWWS[index+w.maxIndex()+1];
          if(wave_type=='r') {
            rSS[n][i].sparseMap00[j]=v00; 
            rSS[n][i].sparseMap90[j]=v90; 
          }
          if(wave_type=='j') {
            jSS[n][i].sparseMap00[j]=v00; 
            jSS[n][i].sparseMap90[j]=v90; 
          }
          if(wave_type=='s') {
            sSS[n][i].sparseMap00[j]=v00; 
            sSS[n][i].sparseMap90[j]=v90; 
          }
        }
     }
   }
}

void AddNoise2Sparse(network* NET, CWB::config* cfg, int seed) {

   //cout << "-----------------> AddNoise2Sparse" << endl; 

   gRandom->SetSeed(seed);

   TF2 *gaus2d = new TF2("gaus2d","exp(-x*x/2)*exp(-y*y/2)", -5,5,-5,5);

   double d2r = TMath::Pi()/180.;
   double C,S;

   // copy waveform sparse to detector sparse map
   int nIFO = NET->ifoListSize();			// number of detectors
   for(int n=0;n<nIFO;n++) {
      detector* pD = NET->getifo(n);
#ifdef USE_wWHT				// use wWHT (default wREC)
      if(vSS[n].size()>0) pD->vSS=vSS[n];
#else
      if(rSS[n].size()>0) pD->vSS=rSS[n];
#endif
   }

   int nRES = cfg->l_high-cfg->l_low+1;     // number of frequency resolution levels

   // random shuffle of vN[i] time indexes
   std::srand(time(0));
   random_shuffle(&(sI.data[0]),&(sI.data[sI.size()])); 

   // add noise/mis-calibration to sparse map
   double r00,r90;
   for(int n=0; n<nIFO; n++) {
     detector* pD = NET->getifo(n);
     if(WRC_AMP_CAL_ERR) {
       WRC_CAL_AMP = gRandom->Uniform(1-AMP_CAL_ERR,1+AMP_CAL_ERR);
       cout << "WRC_CAL_AMP : " << WRC_CAL_AMP << endl;
     }
     if(WRC_PHS_CAL_ERR) {
       WRC_CAL_PHS = gRandom->Uniform(1-PHS_CAL_ERR,1+PHS_CAL_ERR);
       cout << "WRC_CAL_PHS : " << WRC_CAL_PHS << endl;
       C = cos(-WRC_CAL_PHS*d2r);
       S = sin(-WRC_CAL_PHS*d2r);
     }
     for(int i=0; i<nRES; i++) {
        int size = pD->vSS[i].sparseMap00.size();
        for(int j=0; j<size; j++) {
          if(WRC_AMP_CAL_ERR) {
            pD->vSS[i].sparseMap00[j]*=WRC_CAL_AMP; 
            pD->vSS[i].sparseMap90[j]*=WRC_CAL_AMP; 
          }
          if(WRC_PHS_CAL_ERR) {
            double aa = pD->vSS[i].sparseMap00[j]; 
            double AA = pD->vSS[i].sparseMap90[j]; 

            pD->vSS[i].sparseMap00[j] =  C*aa + S*AA; 
            pD->vSS[i].sparseMap90[j] = -S*aa + C*AA ; 
	  }
          if(WRC_NOISE) {
            int index = pD->vSS[i].sparseIndex[j];

            int layers = vN[n][i].maxLayer()+1;   // numbers of frequency bins (first & last bins have df/2)
            int slices = vN[n][i].sizeZero();     // number of time bins
            //cout << "layers " << layers << " slices " << slices << endl;

            int levels = (1<<(nRES-1-i));

            // index = itime*layers + ifreq;
            int ifreq = index%layers;
            int itime = (index-ifreq)/layers;
            itime = levels*sI[itime/levels];
            index = itime*layers + ifreq;		// shuffled index

            r00 = vN[n][i].pWavelet->pWWS[index];
            r90 = vN[n][i].pWavelet->pWWS[index+vN[n][i].maxIndex()+1];
            //gaus2d->GetRandom2(r00,r90);

            pD->vSS[i].sparseMap00[j]+=r00; 
            pD->vSS[i].sparseMap90[j]+=r90; 
          }
        }
     }
   }
}

void SaveSkyProb(network* NET, CWB::config* cfg, int id) {

   std::vector<float>* vP;
   std::vector<int>*   vI;

   // save the probability skymap
   if(NET->wc_List[0].p_Map.size()) {

      vP = &(NET->wc_List[0].p_Map[id-1]);
      vI = &(NET->wc_List[0].p_Ind[id-1]);

      rec_skyprob = NET->getifo(0)->tau;
      rec_skyprob = 0.;

      for(int j=0; j<int(vP->size()); j++) {
         int i = (*vI)[j];
         double th = rec_skyprob.getTheta(i);
         double ph = rec_skyprob.getPhi(i);
         int k=rec_skyprob.getSkyIndex(th, ph);
         rec_skyprob.set(k,(*vP)[j]);
      }
   }

   hist_skyprob.SetBins(rec_skyprob.size(),0,rec_skyprob.size()-1);

   double prob=0;
   for(int l=0;l<rec_skyprob.size();l++) prob+=rec_skyprob.get(l);
   cout << "PROB : " << prob << endl;
   for(int l=0;l<rec_skyprob.size();l++) hist_skyprob.SetBinContent(l,rec_skyprob.get(l));

#ifdef SAVE_SKYPROB
   // plot rec_skyprob 
   gskymap gSM=rec_skyprob;
   //gSM.SetOptions("hammer","Celestial",2);
   gSM.SetZaxisTitle("probability");
   gSM.Draw(1);
   //TH2D* hsm = gSM.GetHistogram();
   //hsm->GetZaxis()->SetTitleOffset(0.85);
   //hsm->GetZaxis()->SetTitleSize(0.03);
   gSM.Print(TString::Format("%s/rec_skyprob.png",PLOT_DIR)); 
   exit(0);
#endif
}

void PlotFinal(int n) {

  // draw envelope
  gwavearray<double> geINJ=wINJ[n];
  CWB::mdc::PhaseShift(geINJ,90);
  for(int j=0;j<geINJ.size();j++) geINJ[j]=sqrt(wINJ[n][j]*wINJ[n][j]+geINJ[j]*geINJ[j]);
  geINJ.Draw(GWAT_TIME);

  gwavearray<double> geAVR=wINJ[n];
  for(int j=0;j<wAVR[n].size();j++) geAVR[j]=wAVR[n][j];
  CWB::mdc::PhaseShift(geAVR,90);
  for(int j=0;j<geAVR.size();j++) geAVR[j]=sqrt(wAVR[n][j]*wAVR[n][j]+geAVR[j]*geAVR[j]);
  geINJ.Draw(&geAVR,GWAT_TIME,"SAME",kRed);

  gwavearray<double> geRMS=wINJ[n];
  for(int j=0;j<wRMS[n].size();j++) geRMS[j]=wRMS[n][j];
  geINJ.Draw(&geRMS,GWAT_TIME,"SAME",kBlue);

  watplot* plot = geINJ.GetWATPLOT();
  TString gfile=TString::Format("%s/eAVR_%d.root",PLOT_DIR,n);
  (*plot) >> gfile;


  // draw frequency
  double dt=1./wINJ[n].rate();
  double Pi = TMath::Pi();
  wavearray<double> xx=wINJ[n];
  wavearray<double> XX=wINJ[n];
  CWB::mdc::PhaseShift(XX,90);
  gwavearray<double> gF=wINJ[n];
  for(int j=1;j<gF.size();j++) {
//double angle = atan2(XX[j],xx[j])*180/Pi;
//angle += ceil( -angle / 360 ) * 360;
//gF[j]=angle;
     gF[j]=atan2(XX[j],xx[j])-atan2(XX[j-1],xx[j-1]);
     if(gF[j]<-Pi) gF[j]+=2*Pi;
     if(gF[j]> Pi) gF[j]-=2*Pi;
     gF[j]/=2*Pi*dt;
  }
  gF.Draw(GWAT_TIME);

  wavearray<double> yy=wINJ[n];
  for(int j=0;j<yy.size();j++) yy[j]=wAVR[n][j];
  wavearray<double> YY=yy;
  CWB::mdc::PhaseShift(YY,90);
  gwavearray<double> gFF=wINJ[n];
  for(int j=1;j<gFF.size();j++) {
     gFF[j]=atan2(YY[j],yy[j])-atan2(YY[j-1],yy[j-1]);
     if(gFF[j]<-Pi) gFF[j]+=2*Pi;
     if(gFF[j]> Pi) gFF[j]-=2*Pi;
     gFF[j]/=2*Pi*dt;
  }
  gF.Draw(&gFF,GWAT_TIME,"SAME",kRed);

  plot = gF.GetWATPLOT();
  gfile=TString::Format("%s/FREQ_%d.root",PLOT_DIR,n);
  (*plot) >> gfile;


  // draw envelope difference
  gwavearray<double> geDIF=wINJ[n];
  for(int j=0;j<geDIF.size();j++) geDIF[j]=geAVR[j]-geINJ[j];
  geDIF.Draw(GWAT_TIME);
  geDIF.Draw(&geRMS,GWAT_TIME,"SAME",kBlue);

  plot = geDIF.GetWATPLOT();
  gfile=TString::Format("%s/eDIF_%d.root",PLOT_DIR,n);
  (*plot) >> gfile;


  // draw waveforms INJ vs AVR
  gwavearray<double> gwAVR=wINJ[n];
  gwAVR.Draw(GWAT_TIME);
  for(int j=0;j<wAVR[n].size();j++) gwAVR[j]=wAVR[n][j];
  gwAVR.Draw(GWAT_TIME,"SAME",kRed);
  for(int j=0;j<wAVR[n].size();j++) gwAVR[j]=wRMS[n][j];
  gwAVR.Draw(GWAT_TIME,"SAME",kBlue);

  plot = gwAVR.GetWATPLOT();
  gfile=TString::Format("%s/wAVR_%d.root",PLOT_DIR,n);
  (*plot) >> gfile;

  // draw waveforms REC vs AVR/INJ

  double bINJ = wINJ[n].start();
  double eINJ = wINJ[n].stop();
  double bREC = wREC[n].start();
  double eREC = wREC[n].stop();
  cout.precision(14);
  cout << "bINJ : " << bINJ << " eINJ : " << eINJ << endl;
  cout << "bREC : " << bREC << " eREC : " << eREC << endl;

  double R=wINJ[n].rate();

  int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
  int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);
  cout << "oINJ : " << oINJ << " oREC : " << oREC << endl;

  double startXCOR = bINJ>bREC ? bINJ : bREC;
  double endXCOR   = eINJ<eREC ? eINJ : eREC;
  int sizeXCOR     = int((endXCOR-startXCOR)*R+0.5);
  cout << "sizeXCOR : " << sizeXCOR << endl;
  cout << "wINJ[n].size() : " << wINJ[n].size() << endl;
  cout << "wREC[n].size() : " << wREC[n].size() << endl;

  gwavearray<double> gwREC = wINJ[n];
  for(int j=0;j<sizeXCOR;j++) gwREC[j+oINJ]=wREC[n][j+oREC];
  gwREC.Draw(GWAT_TIME);
  //for(int j=0;j<wAVR[n].size();j++) gwREC[j]=wAVR[n][j];
  for(int j=0;j<wINJ[n].size();j++) gwREC[j]=wINJ[n][j];
  gwREC.Draw(GWAT_TIME,"SAME",kRed);
  for(int j=0;j<wAVR[n].size();j++) gwREC[j]=wRMS[n][j];
  gwREC.Draw(GWAT_TIME,"SAME",kBlue);

  plot = gwREC.GetWATPLOT();
  gfile=TString::Format("%s/wREC_%d.root",PLOT_DIR,n);
  (*plot) >> gfile;
}

double GetSparseMap(SSeries<double>* SS, bool phase, int index) {

  int layer = SS->GetLayer(index);
  int start = SS->sparseLookup[layer];                // sparse table layer offset
  int end   = SS->sparseLookup[layer+1]-1;            // sparse table layer+1 offset

  int i = SS->binarySearch(SS->sparseIndex.data,start,end,index);

  if(i<0) return 0;

  return phase ? SS->sparseMap00[i] : SS->sparseMap90[i];
}

void PlotFinal2(network* NET, int ID) {

  double tMin,tMax;
  double start;	// use ifoID=0 as offset

  int nIFO = NET->ifoListSize();			// number of detectors

  for(int n=0;n<nIFO;n++) {

    // draw waveforms INJ vs AVR
    gwavearray<double> gw=wINJ[n];
    //gwavearray<double> gw=rINJ[n];
    
    if(n==0) {
      start = gw.start();
      gw.GetTimeRange(tMin, tMax, 0.998);
      cout << "tMin " << tMin << " tMax " << tMax << endl;
    }
    gw.SetTimeRange(start+tMin, start+tMax);
    gw.Draw(GWAT_TIME, "CUSTOM");
    for(int j=0;j<wAVR[n].size();j++) gw[j]=wAVR[n][j];
    gw.Draw(GWAT_TIME,"SAME",kRed);
    for(int j=0;j<wAVR[n].size();j++) gw[j]=wRMS[n][j];
    gw.Draw(GWAT_TIME,"SAME",kBlue);

    watplot* plot = gw.GetWATPLOT();

    char gtitle[256];
    sprintf(gtitle,"%s : %10.3f (gps) - injected(Black) - average(Red) - rms(Blue)",gw.start(),NET->ifoName[n]);
    plot->gtitle(gtitle,"time(sec)","amplitude");

    TString gfile=TString::Format("%s/inj_vs_avr_%s_%d.root",PLOT_DIR,NET->ifoName[n],ID);
    //TString gfile=TString::Format("%s/rinj_vs_avr_%s_%d.root",PLOT_DIR,NET->ifoName[n],ID);
    (*plot) >> gfile;

  }
}

void PlotFinal3(network* NET, int ID, int ifoID, wavearray<double>* w1, wavearray<double>* w2, TString tag, TString gtype) {

  double tMin,tMax;
  double start;	// use ifoID=0 as offset

  // draw waveforms INJ vs AVR
  gwavearray<double> gw=*w1;
    
  start = gw.start();
  gw.GetTimeRange(tMin, tMax, 0.998);
  cout << "tMin " << tMin << " tMax " << tMax << endl;
   
  gw.SetTimeRange(start+tMin, start+tMax);
  gw.Draw(GWAT_TIME, "CUSTOM");
  gw.Draw(w2,GWAT_TIME,"SAME",kRed);

  watplot* plot = gw.GetWATPLOT();

  char gtitle[256];
  sprintf(gtitle,"%s : %10.3f (gps) - %s",gw.start(),NET->ifoName[ifoID],tag.Data());
  plot->gtitle(gtitle,"time(sec)","amplitude");

  TString gfile=TString::Format("%s/%s_%s_%d.%s",PLOT_DIR,tag.Data(),NET->ifoName[ifoID],ID,gtype.Data());
  (*plot) >> gfile;

}

void ComputeErrorStatistic(network* NET, CWB::config* cfg, int ID) {

  int nIFO = NET->ifoListSize();			// number of detectors

  double enINJ[NIFO_MAX], enREC[NIFO_MAX], xcorINJ_REC[NIFO_MAX];

  // shift all data to the first ifo arrival time
  for(int n=0;n<nIFO;n++) {
    CWB::mdc::TimeShift(wREC[n],-wDREC[n]);	
    if(cfg->simulation) CWB::mdc::TimeShift(wINJ[n],-wDINJ[n]);	
    if(cfg->simulation) CWB::mdc::TimeShift(rINJ[n],-wDREC[n]);	
    for(int i=0;i<vREC[n].size();i++) {
      CWB::mdc::TimeShift(vREC[n][i],-vDREC[n][i]);	
    }
  }

  double rnum=0; 
  double rden=0; 
/*
  for(int n=0;n<nIFO;n++) {
    ComputeResidualEnergy(&wINJ[n], &wREC[n], enINJ[n], enREC[n], xcorINJ_REC[n]);
    //cout << n << " RESIDUAL : " << (enINJ[n]+enREC[n]-2*xcorINJ_REC[n])/enINJ[n] << endl;
    rnum+=(enINJ[n]+enREC[n]-2*xcorINJ_REC[n]);
    rden+=enINJ[n];
  }
  cout << "NET RESIDUAL : " << rnum/rden << endl; 
*/
  // copy wINJ (only the wREC time range) into winj
  wavearray<double> winj[NIFO_MAX];
  if(cfg->simulation) {
    for(int n=0;n<nIFO;n++) {
      winj[n]=wREC[n];
      winj[n]=0;

      double bINJ = wINJ[n].start();
      double eINJ = wINJ[n].stop();
      double bREC = wREC[n].start();
      double eREC = wREC[n].stop();

      double R=wINJ[0].rate();

      int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
      int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);

      double startXCOR = bINJ>bREC ? bINJ : bREC;
      double endXCOR   = eINJ<eREC ? eINJ : eREC;
      int sizeXCOR     = int((endXCOR-startXCOR)*R+0.5);

      if (sizeXCOR<=0) continue;

      for (int j=0;j<sizeXCOR;j++) winj[n][j+oREC]+=wINJ[n].data[j+oINJ];
    } 
  } 

  // copy vREC (only the wREC time range) into wavr
  wavearray<double> wavr[NIFO_MAX];
  for(int n=0;n<nIFO;n++) {
    wavr[n]=wREC[n];
    wavr[n]=0;
    double bINJ = wREC[n].start();
    double eINJ = wREC[n].stop();
    for(int i=0;i<vREC[n].size();i++) {

      double bREC = vREC[n][i].start();
      double eREC = vREC[n][i].stop();
      //cout.precision(14);
      //cout << "bINJ : " << bINJ << " eINJ : " << eINJ << endl;
      //cout << "bREC : " << bREC << " eREC : " << eREC << endl;

      double R=wREC[n].rate();

      int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
      int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);
      //cout << "oINJ : " << oINJ << " oREC : " << oREC << endl;

      double startXCOR = bINJ>bREC ? bINJ : bREC;
      double endXCOR   = eINJ<eREC ? eINJ : eREC;
      int sizeXCOR     = int((endXCOR-startXCOR)*R+0.5);
      //cout << "startXCOR : " << startXCOR << " endXCOR : " 
      //     << endXCOR << " sizeXCOR :" << sizeXCOR << endl;

      if (sizeXCOR<=0) continue;

      for (int j=0;j<sizeXCOR;j++) wavr[n][j+oINJ]+=vREC[n][i].data[j+oREC];
    }
    wavr[n]*=1./vREC[n].size();
  } 

  wavearray<float> vRES(vREC[0].size());

  for(int i=0;i<vREC[0].size();i++) {
    rnum=0; 
    rden=0; 
    for(int n=0;n<nIFO;n++) {
      ComputeResidualEnergy(&wavr[n], &vREC[n][i], enINJ[n], enREC[n], xcorINJ_REC[n]);
      //cout << i << " " << n << " RESIDUAL : " << (enINJ[n]+enREC[n]-2*xcorINJ_REC[n])/enINJ[n] << endl;
      rnum+=(enINJ[n]+enREC[n]-2*xcorINJ_REC[n]);
      rden+=enINJ[n];
    }
    //cout << "vRECvsAVR NET RESIDUAL : " << rnum/rden << endl; 
    vRES[i]=rnum/rden;
  }

  double jres=0;
  if(cfg->simulation) {
    rnum=0; 
    rden=0; 
    for(int n=0;n<nIFO;n++) {
      ComputeResidualEnergy(&wavr[n], &winj[n], enINJ[n], enREC[n], xcorINJ_REC[n]);
      //cout << n << " RESIDUAL : " << (enINJ[n]+enREC[n]-2*xcorINJ_REC[n])/enINJ[n] << endl;
      rnum+=(enINJ[n]+enREC[n]-2*xcorINJ_REC[n]);
      rden+=enINJ[n];
    }
    //cout << "INJvsAVR NET RESIDUAL : " << rnum/rden << endl; 

    rnum=0; 
    rden=0; 
    for(int n=0;n<nIFO;n++) {
      ComputeResidualEnergy(&wavr[n], &rINJ[n], enINJ[n], enREC[n], xcorINJ_REC[n]);
      //cout << n << " RESIDUAL : " << (enINJ[n]+enREC[n]-2*xcorINJ_REC[n])/enINJ[n] << endl;
      rnum+=(enINJ[n]+enREC[n]-2*xcorINJ_REC[n]);
      rden+=enINJ[n];
    }
    jres=rnum/rden;
    //cout << "rINJvsAVR NET RESIDUAL : " << rnum/rden << endl; 
  }

  rnum=0; 
  rden=0; 
  for(int n=0;n<nIFO;n++) {
    ComputeResidualEnergy(&wavr[n], &wREC[n], enINJ[n], enREC[n], xcorINJ_REC[n]);
    //cout << n << " RESIDUAL : " << (enINJ[n]+enREC[n]-2*xcorINJ_REC[n])/enINJ[n] << endl;
    rnum+=(enINJ[n]+enREC[n]-2*xcorINJ_REC[n]);
    rden+=enINJ[n];
  }
  //cout << "RECvsAVR NET RESIDUAL : " << rnum/rden << endl; 

#ifdef SAVE_WF_COMPARISON
  for(int n=0;n<nIFO;n++) {
//    PlotFinal3(NET, ID, n, &wavr[n], &wavr[n], "AVRvsAVR","root");
//    PlotFinal3(NET, ID, n, &winj[n], &winj[n], "INJvsINJ","root");
//    PlotFinal3(NET, ID, n, &wREC[n], &wREC[n], "RECvsREC","root");

    PlotFinal3(NET, ID, n, &rINJ[n], &winj[n], "rINJvsINJ","root");
    PlotFinal3(NET, ID, n, &rINJ[n], &wREC[n], "rINJvsREC","root");
    PlotFinal3(NET, ID, n, &rINJ[n], &wavr[n], "rINJvsAVR","root");

    PlotFinal3(NET, ID, n, &winj[n], &wavr[n], "INJvsAVR","root");
    PlotFinal3(NET, ID, n, &wREC[n], &wavr[n], "RECvsAVR","root");
  }
#endif

  wavearray<int> index(vRES.size());
  TMath::Sort(int(vRES.size()),const_cast<float*>(vRES.data),index.data,false);

  erR[0] = jres;

  for(int i=1;i<10;i++) erR[i]=vRES[index[i*int(vRES.size()/10.)]];

  erR[10]=0;
  for(int i=0;i<vRES.size();i++)  {
    if(vRES[i]<=jres) erR[10]+=1.;
  }
  erR[10]/=vRES.size();

  cout.precision(3);
  cout << "erR : ";
  for(int i=0;i<10;i++) cout << erR[i] << "/";
  cout << erR[10] << endl;

#ifdef SAVE_PLOT_RESIDUALS

  // plot residuals
  TCanvas* canvas2= new TCanvas("canvas2", "canvas2", 200, 20, 600, 600);
  TH1F* hres = new TH1F("residuals","residuals",10,0,1);
  for(int i=0;i<vRES.size();i++)  hres->Fill(vRES[i]);

  hres->Draw();
  canvas2->SetLogy();

  float ymax = hres->GetMaximum();
  TLine *line = new TLine(jres,0,jres,ymax);
  line->SetLineColor(kRed);
  line->Draw();

  gStyle->SetLineColor(kWhite);
  hres->SetTitle("Distribution of residuals"); 
  hres->SetStats(kTRUE); 
  canvas2->Print(TString::Format("%s/residuals_%d.png",PLOT_DIR,ID));

  delete hres;
  delete canvas2;

#endif
}

void Clear() {

  sI.resize(0);
  for(int n=0;n<NIFO_MAX;n++) {

    wINJ[n].resize(0);
    sINJ[n].resize(0);
    wREC[n].resize(0);

    while(!vSS[n].empty()) vSS[n].pop_back();
    vSS[n].clear(); std::vector<SSeries<double> >().swap(vSS[n]);
    while(!sSS[n].empty()) sSS[n].pop_back();
    sSS[n].clear(); std::vector<SSeries<double> >().swap(sSS[n]);
    while(!rSS[n].empty()) rSS[n].pop_back();
    rSS[n].clear(); std::vector<SSeries<double> >().swap(rSS[n]);
    while(!jSS[n].empty()) jSS[n].pop_back();
    jSS[n].clear(); std::vector<SSeries<double> >().swap(jSS[n]);

    while(!vN[n].empty()) vN[n].pop_back();
    vN[n].clear(); std::vector<WSeries<double> >().swap(vN[n]);

    while(!wTRY[n].empty()) wTRY[n].pop_back();
    wTRY[n].clear(); std::vector<int >().swap(wTRY[n]);

    while(!wAVR[n].empty()) wAVR[n].pop_back();
    wAVR[n].clear(); std::vector<double >().swap(wAVR[n]);

    while(!wRMS[n].empty()) wRMS[n].pop_back();
    wRMS[n].clear(); std::vector<double >().swap(wRMS[n]);

    while(!vDREC[n].empty()) vDREC[n].pop_back();
    vDREC[n].clear(); std::vector<double >().swap(vDREC[n]);

    while(!vREC[n].empty()) vREC[n].pop_back();
    vREC[n].clear(); std::vector<wavearray<double> >().swap(vREC[n]);

  }
}

void ReadUserOptions() {

  TString cwb_inet_options=TString(gSystem->Getenv("CWB_INET_OPTIONS"));
  if(cwb_inet_options.CompareTo("")!=0) {
    TString inet_options = cwb_inet_options;
    cout << inet_options << endl;
    if(!inet_options.Contains("--")) {  // parameters are used only by cwb_inet
  
      TObjArray* token = TString(inet_options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("wrc_id=")) {
          TString wrc_id=stok;
          wrc_id.Remove(0,wrc_id.Last('=')+1);
          if(wrc_id.IsDigit()) WRC_ID=wrc_id.Atoi();
          WRC_ID_CONFIG=WRC_ID;
          cout << "WRC_ID : " << WRC_ID << endl;
        }

        if(stok.Contains("wrc_retry=")) {
          TString wrc_retry=stok;
          wrc_retry.Remove(0,wrc_retry.Last('=')+1);
          if(wrc_retry.IsDigit()) WRC_RETRY=wrc_retry.Atoi();
          WRC_RETRY_CONFIG=WRC_RETRY;
          cout << "WRC_RETRY : " << WRC_RETRY << endl;
        }

        if(stok.Contains("wrc_ced_id=")) {
          TString wrc_ced_id=stok;
          wrc_ced_id.Remove(0,wrc_ced_id.Last('=')+1);
          if(wrc_ced_id.IsDigit()) WRC_CED_ID=wrc_ced_id.Atoi();
          WRC_CED_ID_CONFIG=WRC_CED_ID;
          cout << "WRC_CED_ID : " << WRC_CED_ID << endl;
        }

        if(stok.Contains("wrc_ced_retry=")) {
          TString wrc_ced_retry=stok;
          wrc_ced_retry.Remove(0,wrc_ced_retry.Last('=')+1);
          if(wrc_ced_retry.IsDigit()) WRC_CED_RETRY=wrc_ced_retry.Atoi();
          WRC_CED_RETRY_CONFIG=WRC_CED_RETRY;
          cout << "WRC_CED_RETRY : " << WRC_CED_RETRY << endl;
        }

        if(stok.Contains("wrc_skymask=")) {
          TString wrc_skymask=stok;
          wrc_skymask.Remove(0,wrc_skymask.Last('=')+1);
          if(wrc_skymask=="true")  WRC_SKYMASK=true;
          if(wrc_skymask=="false") WRC_SKYMASK=false;
          WRC_SKYMASK_CONFIG=WRC_SKYMASK;
        }

        if(stok.Contains("wrc_noise=")) {
          TString wrc_noise=stok;
          wrc_noise.Remove(0,wrc_noise.Last('=')+1);
          if(wrc_noise=="true")  WRC_NOISE=true;
          if(wrc_noise=="false") WRC_NOISE=false;
          WRC_NOISE_CONFIG=WRC_NOISE;
        }

      }
    }
  }
}

int _sse_mra_ps(network* NET, float* amp, float* AMP, float Eo, int K) {
// fast multi-resolution analysis inside sky loop                         
// select max E pixel and either scale or skip it based on the value of residual 
// pointer to 00 phase amplitude of monster pixels                               
// pointer to 90 phase amplitude of monster pixels                               
// Eo - energy threshold                                                         
//  K - max number of principle components to extract                                
// returns number of MRA pixels                                                      

   int j,n,mm,J;                                                                   
   int k = 0;                                                                      
   int m = 0;                                                                      
   int f = NIFO/4;                                                                 
   int V = (int)NET->rNRG.size();                                                 
   float*  ee = NET->rNRG.data;                            // residual energy     
   float*  c  = NULL;                                                              
   float   E2 = Eo/2;                                       // threshold           
   float   E;                                                                      
   float mam[NIFO];                                                                
   float mAM[NIFO];                                                                

   __m128* _m00 = (__m128*) mam;
   __m128* _m90 = (__m128*) mAM;
   __m128* _amp = (__m128*) amp;
   __m128* _AMP = (__m128*) AMP;
   __m128* _a00 = (__m128*) NET->a_00.data;
   __m128* _a90 = (__m128*) NET->a_90.data;

   while(k<K){

      for(j=0; j<V; ++j) if(ee[j]>ee[m]) m=j;               // find max pixel
      if(ee[m]<=Eo) break;  mm = m*f;                                        

      //cout<<k<<" "<<" V= "<<V<<" m="<<m<<" ee[m]="<<ee[m];

      float cc = 0.; 

      _sse_zero_ps(_m00);
      _sse_zero_ps(_m90);

      J = NET->wdmMRA.size()/7;                              
      c = NET->wdmMRA.getXTalk(m);                          // c1*c2+c3*c4=c1*c3+c2*c4=0
      for(j=0; j<J; j++) {                                                       
         n = int(c[0]+0.1);                                                      
         if(ee[n]>Eo) {                                                          
            _sse_rotadd_ps(_a00+n*f,c[1],_a90+n*f,c[3],_m00);    // construct 00 vector
            _sse_rotadd_ps(_a00+n*f,c[2],_a90+n*f,c[4],_m90);    // construct 90 vector
         }                                                                             
         if(ee[n]>0) cc += c[5];                                                       
         c += 7;                                                                           
      }                                                                                     
      _sse_mul_ps(_m00,cc>1?1./cc:0.7);                                                     
      _sse_mul_ps(_m90,cc>1?1./cc:0.7);                                                     
      E = _sse_abs_ps(_m00,_m90);                           // get PC energy                

      if(E > ee[m]) {                                       // correct overtuning PC
         _sse_cpf_ps(mam,_a00+mm);                          // store a00 for max pixel
         _sse_cpf_ps(mAM,_a90+mm);                          // store a90 for max pixel
      }                                                                               
      _sse_add_ps(_amp+mm,_m00);                            // update 00 PC           
      _sse_add_ps(_AMP+mm,_m90);                            // update 90 PC           

      J = NET->wdmMRA.size()/7;                              
      c = NET->wdmMRA.getXTalk(m);                          // c1*c2+c3*c4=c1*c3+c2*c4=0
      for(j=0; j<J; j++) {                                                              
         n = int(c[0]+0.1);                                                             
         if(E<E2 && n!=m) {c+=7; continue;}                                             
         if(ee[n]>Eo) {                                                                 
            ee[n] = _sse_rotsub_ps(_m00,c[1],_m90,c[2],_a00+n*f);    // subtract PC from a00
            ee[n]+= _sse_rotsub_ps(_m00,c[3],_m90,c[4],_a90+n*f);    // subtract PC from a90
            ee[n]+= 1.e-6;                                                                  
         }                                                                                  
         c += 7;                                                                            
      }                                                                                     
      //cout<<" "<<ee[m]<<" "<<k<<" "<<E<<" "<<EE<<" "<<endl;                                             
      //cout<<" "<<m<<" "<<cc<<" "<<ee[m]<<" "<<E<<" "<<pp[m]<<endl;                                      
      k++;                                                                                                
   }                                                                                                      
   return k;                                                                                              
}                                                                                        

void SetOutput(network* NET, netevent* &EVT, CWB::config* cfg) {

   // import slagShift
   float* gSLAGSHIFT=NULL; IMPORT(float*,gSLAGSHIFT)

   int nIFO = NET->ifoListSize();			// number of detectors

   // search output root file in the system list
   TFile* froot = NULL;                         
   TList *files = (TList*)gROOT->GetListOfFiles();
   outDump="";
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
          outDump=fname;
          outDump.ReplaceAll(".root.tmp",".txt");
          //cout << "output file name : " << fname << endl;
        }
     }                                                               
     if(!froot) {                                                    
       cout << "CWB_Plugin_xWRC.C : Error - output root file not found" << endl;
       gSystem->Exit(1);                                                                             
     }                                                                                      
   } else {                                                                                 
     cout << "CWB_Plugin_xWRC.C : Error - output root file not found" << endl;  
     gSystem->Exit(1);                                                                               
   }                                                                                        

   net_tree = (TTree *) froot->Get("waveburst");
   if(net_tree!=NULL) {
     EVT = new netevent(net_tree,nIFO);
     for(int j=0;j<nPAR;j++) net_tree->SetBranchAddress(parName[j],&parValue[j]);
     net_tree->SetBranchAddress("erR",erR);
   } else {
     EVT = new netevent(nIFO);
     net_tree = EVT->setTree();
     for(int j=0;j<nPAR;j++) 
        net_tree->Branch(parName[j],&parValue[j],TString::Format("%s/F",parName[j].Data()));
     net_tree->Branch("erR",erR,TString::Format("erR[%i]/F",11));
   }
   EVT->setSLags(gSLAGSHIFT);                          // set slags into netevent
   EVT->Psave=cfg->Psave;
}

void DumpOutput(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor) {

   for(int j=0;j<11;j++) NET->getwc(k)->sArea[ID-1][j]=rec_erA[j]; // restore erA
   if(cfg->dump) EVT->dopen(outDump.Data(),const_cast<char*>("a"),false);
   EVT->output2G(net_tree,NET,ID,k,factor);		// get reconstructed parameters
   if(cfg->dump) EVT->dclose();
   if(!cfg->cedDump) NET->getwc(k)->sCuts[ID-1]=1; 	// mark as processed                              
}

void GetRecWaveform(network* NET, netevent* EVT, int ID, int k, int factor) {

   // import retry
   int gLRETRY=-1; IMPORT(int,gLRETRY)

   int nIFO = NET->ifoListSize();			// number of detectors

   EVT->output2G(NULL,NET,ID,k,factor);		// get reconstructed parameters

   cout << "rec_theta : " << EVT->theta[0] << " rec_phi : " << EVT->phi[0] << endl;

   wavearray<double>** pwfREC = new wavearray<double>*[nIFO];
   for(int n=0; n<nIFO; n++) {

      detector* pd = NET->getifo(n);
      int idSize = pd->RWFID.size();
      int wfIndex=-1;
      for(int mm=0; mm<idSize; mm++) if (pd->RWFID[mm]==ID) wfIndex=mm;

      if(wfIndex<0) {
         cout << "CWB_Plugin_xWRC.C : Error : Reconstructed waveform not saved !!! : ID -> "
              << ID << " : detector " << NET->ifoName[n] << endl;
         continue;
      }
      if(wfIndex>=0) pwfREC[n] = pd->RWFP[wfIndex];
      wavearray<double>* wfREC = pwfREC[n];	// array of reconstructed waveforms

      double delay = EVT->time[n]-EVT->time[0];	 
      //cout << "DEB DELAY : " << n << " " << delay << endl;

      // save tst rec waveforms
      if(gLRETRY<WRC_RETRY) {	
         vREC[n].push_back(*wfREC);
	 vDREC[n].push_back(delay);
      }
      // save inj/rec waveforms
      if(gLRETRY==WRC_RETRY) {		// first time
	 wREC[n]=*wfREC;
	 wDREC[n]=delay;

         // save injected/reconstruced sky location
         rec_theta = EVT->theta[0];	// rec theta
         rec_phi   = EVT->phi[0];		// rec phi
              
         for(int j=0;j<11;j++) rec_erA[j] = EVT->erA[j];
      }

      // save wrc parameters
      parValue[0]=WRC_RETRY-gLRETRY;
      if(WRC_NOISE) parValue[1]=1; else parValue[1]=0;
      parValue[2]=WRC_SM_PHI;
      parValue[3]=WRC_SM_THETA;
      parValue[4]=WRC_SM_RADIUS;
      parValue[5]=WRC_CAL_AMP;
      parValue[6]=WRC_CAL_PHS;
   }
   delete [] pwfREC;
}

void GetInjWaveform(network* NET, netevent* EVT, int ID, int k, int factor) {

   // import retry
   int gLRETRY=-1; IMPORT(int,gLRETRY)

   int nIFO = NET->ifoListSize();			// number of detectors

   double recTime   = EVT->time[0];		// rec event time det=0

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
            cout << "CWB_Plugin_xWRC.C : Error : Injected waveform not saved !!! : detector "
                 << NET->ifoName[n] << endl;
            continue;
         }

         double rFactor = 1.;
         rFactor *= factor>0 ? factor : 1;
         wavearray<double>* wfINJ = pwfINJ[n];	// array of injected waveforms
         *wfINJ*=rFactor;			// injected wf is multiplied by the custom factor

         // save inj waveforms
         wINJ[n]=*wfINJ;

         //double delay = INJ.time[n]-INJ.time[0];
         double delay = EVT->time[n]-EVT->time[0];
   	 sINJ[n]=*wfINJ; CWB::mdc::TimeShift(sINJ[n],-delay);	
	 wDINJ[n]=delay;

         wTRY[n].resize(wINJ[n].size());
         wAVR[n].resize(wINJ[n].size());
         wRMS[n].resize(wINJ[n].size());
         for(int j=0;j<wINJ[n].size();j++) {wTRY[n][j]=0;wAVR[n][j]=0;wRMS[n][j]=0;};

         // save injected/reconstruced sky location
         inj_theta = EVT->theta[1];		// inj theta
         inj_phi   = EVT->phi[1];		// inj phi

         *wfINJ*=1./rFactor;			// restore injected amplitude
      }
      delete [] pwfINJ;
   }
}

void GetRInjWaveform(network* NET, netevent* EVT, CWB::config* cfg, int ID, int k) {

  int nIFO = NET->ifoListSize();			// number of detectors

  netcluster* pwc = NET->getwc(k);

  std::vector<int>* vint = &(pwc->cList[ID-1]);         // pixel list
  int V = vint->size();                                 // cluster size
  if(!V) return;                                                       

  int irate=0;
  if(cfg->optim) {                                      // extract optimal rate
    vector_int* pv = pwc->cRate.size() ? &(pwc->cRate[ID-1]) : NULL;          
    irate = pv!=NULL ? (*pv)[0] : 0;                                           
  }
  bool isPCs = !(cfg->optim&&std::isupper(cfg->search)); // are Principal Components ?

  detector* pD = NET->getifo(0);
  double R = pD->TFmap.rate();

  // ------------------------------------------------------------------
  // Find max pixel parameters                                          
  // ------------------------------------------------------------------ 

  int K=0;
  for(int n=0; n<V; n++) {                                       
    netpixel* pix = pwc->getPixel(ID,n);                        
    if(pix->layers%2==0) {                                       
      cout << "CWB_Plugin_eWRC.C - Error : is enabled only for WDM (2G)" << endl;
      exit(1);                                                                                     
    }                                                                                              
    for(int m=0; m<nIFO; m++) {        
#ifdef APPLY_INJ_MRA
      NET->a_00[n*NIFO+m]=0;                                                            
      NET->a_90[n*NIFO+m]=0;                                                            
#endif
      pix->setdata(0,'S',m);                                      
      pix->setdata(0,'P',m);                                      
    }                                                                                       
#ifdef APPLY_INJ_MRA
    NET->rNRG[n]=0;
#endif
    if(!pix->core) continue;                  // select only the principal components pixels       
    if((irate)&&(irate != int(pix->rate+0.5))) continue;  // if irate!=0 skip rate!=irate          
    for(int m=0; m<nIFO; m++) {                                                                    
      int index = pix->data[m].index;                                                              
      int ires = int(TMath::Log2(R/pix->rate))-cfg->l_low;                                          
      double aa = GetSparseMap(&sSS[m][ires], true, index);                                    
      double AA = GetSparseMap(&sSS[m][ires], false , index);                                    
      pix->setdata(aa,'S',m);                                      
      pix->setdata(AA,'P',m);                                      
#ifdef APPLY_INJ_MRA
      NET->a_00[n*NIFO+m]=aa;                                                            
      NET->a_90[n*NIFO+m]=AA;                                                            
      NET->rNRG[n]+=aa*aa+AA*AA;
#endif
    }                                                                                       
K++;
  }

#ifdef SAVE_MRA_PLOT					// monster event display
  PlotMRA(NET, ID, k, "rINJ");
#endif

#ifdef APPLY_INJ_MRA

  int V4  = V + (V%4 ? 4 - V%4 : 0);
  int V44 = V4 + 4;

  wavearray<float>  amp(NIFO*V44); amp=0;
  wavearray<float>  AMP(NIFO*V44); AMP=0;
  float Eo = 1e-3;	// must be > 1e-6
  //float Eo = 2*cfg->Acore*cfg->Acore*nIFO;

  _sse_mra_ps(NET, amp.data, AMP.data, Eo, K); 
  for(int n=0; n<V; n++) {                                       
    netpixel* pix = pwc->getPixel(ID,n);                        
    if(pix->layers%2==0) {                                       
      cout << "CWB_Plugin_eWRC.C - Error : is enabled only for WDM (2G)" << endl;
      exit(1);                                                                                     
    }                                                                                              
    for(int m=0; m<nIFO; m++) {        
      pix->setdata(0,'S',m);                                      
      pix->setdata(0,'P',m);                                      
    }                                                                                       
    if(!pix->core) continue;                  // select only the principal components pixels       
    if((irate)&&(irate != int(pix->rate+0.5))) continue;  // if irate!=0 skip rate!=irate          
    for(int m=0; m<nIFO; m++) {                                                                    
      double aa = amp[n*NIFO+m];                                    
      double AA = AMP[n*NIFO+m];                                    
      pix->setdata(aa,'S',m);                                      
      pix->setdata(AA,'P',m);                                      
    }                                                                                       
  }

#ifdef SAVE_MRA_PLOT					// monster event display
  PlotMRA(NET, ID, k, "mra_rINJ");
#endif

#endif

  //NET->getMRAwave(ID,0,'S',0,true);
  NET->getMRAwave(ID,0,'S',0,false);

  for(int n=0; n<nIFO; n++) {
    detector* pD = NET->getifo(n);
    rINJ[n]=pD->waveForm;
    rINJ[n].start(pwc->start+pD->waveForm.start()); 
    CWB::mdc::TimeShift(rINJ[n],wDINJ[n]);
  }
  return;
}
