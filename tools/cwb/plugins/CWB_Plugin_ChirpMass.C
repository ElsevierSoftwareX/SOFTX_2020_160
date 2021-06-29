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

#define _USE_EBBH

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
// DEFINES
// ---------------------------------------------------------------------------------

#define NCHIRP_MAX      4

#define CHI2_THR       -2.5    	// if negative the chirp pixels are tagged with pix->likelihood=0 after the selection
#define TMERGER_CUT     0.1
#define ZMAX_THR        0.2

#define DUMP_WAVE	false	// dump rec/inj/wht waveforms to output root file
#define DUMP_CLUSTER	false	// dump pixel cluster to output root file

// ---------------------------------------------------------------------------------
// USER CONFIG OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  float chi2_thr;
  float tmerger_cut;
  float zmax_thr;

  bool  dump_wave;
  bool  dump_cluster;
};

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void ClearWaveforms(detector* ifo);

void SetOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, bool dump_ebbh);
void DumpOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor);

std::vector<netpixel> GetCluster(network* NET, CWB::config* cfg, int lag, int id);

void SetEventWindow(CWB::config* cfg, double gps);

double GetCentralTime(wavearray<double> x);
double GetTimeBoundaries(wavearray<double> x, double P, double& bT, double& eT);
double GetFrequencyBoundaries(wavearray<double> x, double P, double& bF, double& eF);
double GetCentralFrequency(wavearray<double> x);

std::vector<netpixel> DoPCA(network* NET, CWB::config* cfg, int lag, int id);
void ResetPCA(network* NET, CWB::config* cfg, netcluster* pwc, std::vector<netpixel>* vPIX, int ID);

std::vector<wavearray<double> > GetInjWaveform(network* NET, netevent* EVT, int id, double factor);
std::vector<wavearray<double> > GetRecWaveform(network* NET, netevent* EVT, int id);

wavearray<double> GetWaveform(int ifoId, network* NET, int lag, int id, char type, bool shift=true);

wavearray<double> GetDifWaveform(wavearray<double>* wf1, wavearray<double>* wf2);
wavearray<double> GetAlignedWaveform(wavearray<double>* wf1, wavearray<double>* wf2);
wavearray<double> AddWaveforms(wavearray<double>* wf1, wavearray<double>* wf2);

void PlotWaveform(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wf1, 
                  wavearray<double>* wf2, wavearray<double>* wf3, wavearray<double>* wref, bool fft=false,TString pdir="", double P=0.99);

void GetEccentricityIndex(network* NET, int lag, int id);

void ReadUserOptions(CWB::config* cfg);
void PrintUserOptions();

// ----------------------------------------------------
// ROOT Output EBBH Parameters
// ----------------------------------------------------

struct _EBBH {				// structure for output estimated parameters

  int   nch; 				// number of analyzed chirp 
  float ei;				// eccentricity index
  float ch_mass[NCHIRP_MAX];		// chirp mass
  float ch_tmerger[NCHIRP_MAX];		// merger time
  float ch_merr[NCHIRP_MAX];		// mchirp error
  float ch_mchi2[NCHIRP_MAX];		// mchirp chi2
  float ch_energy[NCHIRP_MAX];		// chirp energy

  wavearray<double>* wDAT[NIFO_MAX];	// whitened data (in the vREC time range)
  wavearray<double>* wINJ[NIFO_MAX];	// whiten injected waveforms
  wavearray<double>* sINJ[NIFO_MAX];	// strain injected waveforms
  wavearray<double>* wREC[NIFO_MAX];	// whiten injected waveforms
  wavearray<double>* sREC[NIFO_MAX];	// strain injected waveforms

  double segGPS;			// segment start time
  float  cRATE;				// netcluste::rate param

  std::vector<netpixel> *pCLUSTER;	// netcluster
};

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions        gOPT;                   // global User Options

_EBBH		gEBBH;			// global output for estimated parameters
TTree* 		gTREE;			// output tree file name
TString 	gOUTPUT;		// output root file name
cwb 		gCWB;

wavearray<double> gHOT[NIFO_MAX];       // HoT time series used to save whitened data

double		gINRATE;		// input data sampling rate
int 		gRATEANA;          	// analysis data rate

wavearray<double> gwDAT[NIFO_MAX];	// whitened data (in the vREC time range)
wavearray<double> gwINJ[NIFO_MAX];	// whiten injected waveforms
wavearray<double> gsINJ[NIFO_MAX];	// strain injected waveforms
wavearray<double> gwREC[NIFO_MAX];	// whiten injected waveforms
wavearray<double> gsREC[NIFO_MAX];	// strain injected waveforms

double		gSEGGPS;		// segment start time
float		gCRATE;			// netcluster::rate
std::vector<netpixel> gCLUSTER;		// netcluster

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Plugin used for eBBH Parameter Estimation

  if(type==CWB_PLUGIN_CONFIG) {  
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)

    if(gCWB2G->istage!=CWB_STAGE_FULL) {cout<< "CWB_Plugin_ChirpMass - Error : EBBH can be executed only in CWB_STAGE_FULL mode" << endl;exit(1);}

    ReadUserOptions(cfg);       			// user config options : read from parPlugin

    if(cfg->pattern<=0) {
      cout << "CWB_Plugin_ChirpMass Error : EBBH enable only with pattern>0 !!!" << endl;
      exit(1);
    }
    cfg->outPlugin=true;  				// disable built-in output root file
    gINRATE=cfg->inRate;				// input data sampling rate
    gRATEANA=gCWB2G->rateANA;          			// analysis data rate

    gEBBH.ei = 1.e-20;
    gEBBH.nch = NCHIRP_MAX;
    for(int i=0;i<NCHIRP_MAX;i++) {
      gEBBH.ch_mass[i]    = 0;
      gEBBH.ch_tmerger[i] = 0;
      gEBBH.ch_merr[i]    = 0;
      gEBBH.ch_mchi2[i]   = 0;
      gEBBH.ch_energy[i]  = 0;
    }
  }

  if(type==CWB_PLUGIN_NETWORK) {
    PrintUserOptions();      				// print config options
  }

  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {
    // save whitened HoT
    int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==NET->getifo(n)->Name) {ifoID=n;break;}
    gHOT[ifoID] = *x;
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {			// AFTER EVENT RECONSTRUCTION 

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

      for(int j=0; j<(int)id.size(); j++) {  		// loop over cluster index

        int ID = size_t(id.data[j]+0.5);

        if(NET->getwc(k)->sCuts[ID-1]!=-1)  continue; 	// skip rejected/processed clusters

        EVT->output2G(NULL,NET,ID,k,ofactor);		// get reconstructed parameters

        // print event parameters
        cout << endl;
        cout << "event parameters : ID -> " << ID << endl;
        for(int n=0;n<nIFO;n++) printf("rec_time %s : %.4f\n",NET->ifoName[n], EVT->time[n]);
        cout << "rec_theta : " << EVT->theta[0] << " rec_phi : " << EVT->phi[0] << endl;
        cout << "SNRnet : " << sqrt(EVT->likelihood) << " netcc[0] : " << EVT->netcc[0] 
             << " rho[0] : " << EVT->rho[0] << " size : " << EVT->size[0] << endl;

        gSEGGPS = EVT->gps[0];				// segment start time

        // get waveforms
        if(cfg->simulation) {
          std::vector<wavearray<double> > vINJ;		// injected
          vINJ = GetInjWaveform(NET, EVT, ID, factor);  // get injected waveform
          if(vINJ.size()!=nIFO) {
            cout << "CWB_Plugin_ChirpMass Error : Injection Waveform Not Found !!!" << endl;
            exit(1);
          }
          for(int n=0;n<nIFO;n++) gwINJ[n] = vINJ[n];
        }
        for(int n=0;n<nIFO;n++) {
          gwREC[n] = GetWaveform(n, NET, k, ID, 'S');  		// get whiten reconstructed waveform
          gsREC[n] = GetWaveform(n, NET, k, ID, 's');  		// get strain reconstructed waveform
          gwDAT[n] = GetAlignedWaveform(&gHOT[n], &gwREC[n]); 	// get whitened data (in the vREC time range)
        }
        gCRATE = NET->getwc(k)->rate;
        gCLUSTER = GetCluster(NET, cfg, k, ID);		// get netcluster

        GetEccentricityIndex(NET, k, ID);		// get eccentricity index

        SetOutputFile(NET, EVT, cfg, true);             // set output root file
        DumpOutputFile(NET, EVT, cfg, ID, k, ofactor); 	// dump event to output root file

        for(int n=0;n<nIFO;n++) {
          detector* pD = NET->getifo(n);
          if(!cfg->simulation) ClearWaveforms(pD);     	// release waveform memory
        }
      }
    }

    jfile->cd();
    if(EVT) delete EVT;

    while(!gCLUSTER.empty()) gCLUSTER.pop_back();
    gCLUSTER.clear(); std::vector<netpixel>().swap(gCLUSTER);
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
           cout << "CWB_Plugin_ChirpMass.C : Error : Injected waveform not saved !!! : detector "
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
        cout << "CWB_Plugin_ChirpMass.C : Error : Reconstructed waveform not saved !!! : ID -> "
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

wavearray<double> GetWaveform(int ifoId, network* NET, int lag, int id, char type, bool shift) {

  if(type!='S' && type!='s' && type!='W' && type!='w' && type!='N' && type!='n') {
    cout << "GetWaveform : not valid type : Abort" << endl; exit(1);
  } 

  netcluster* pwc = NET->getwc(lag);
  detector* pd = NET->getifo(ifoId);
  wavearray<double>* wf;

  if(type=='S' || type=='s') {
    NET->getMRAwave(id,lag,type,0,shift);       // reconstruct whitened/strain shifted pd->waveForm
    wf = &pd->waveForm;
    wf->start(pwc->start+pd->waveForm.start());
  }

  if(type=='W' || type=='w') {
    NET->getMRAwave(id,lag,type,0,shift);       // reconstruct+noise whitened shifted pd->waveBand
    wf = &pd->waveBand;
    wf->start(pwc->start+pd->waveBand.start());
  }

  if(type=='N' || type=='n') {
    NET->getMRAwave(id,lag,'W',0,shift);       // reconstruct+noise whitened shifted pd->waveBand
    NET->getMRAwave(id,lag,'S',0,shift);       // reconstruct whitened shifted pd->waveForm
    pd->waveNull = pd->waveBand;
    pd->waveNull-= pd->waveForm;
    wf = &pd->waveNull;
    wf->start(pwc->start+pd->waveNull.start());
  }

  return *wf;
}

void PlotWaveform(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wf1, 
                  wavearray<double>* wf2, wavearray<double>* wf3, wavearray<double>* wref, bool fft, TString pdir, double P) {

  watplot PTS(const_cast<char*>("ptspe"),200,20,800,500);

  //cout << "Print " << ofname << endl;
  double tmin=1.e20;
  if(wref==NULL) return;
  if(wf1==NULL) return;
  else          tmin=wf1->start();
  if(wf2!=NULL) if(wf2->start()<tmin) tmin=wf2->start();
  if(wf3!=NULL) if(wf3->start()<tmin) tmin=wf3->start();

                                           wf1->start(wf1->start()-tmin);
  if(wf2!=NULL) if(wf2!=wf1)               wf2->start(wf2->start()-tmin);
  if(wf3!=NULL) if(wf3!=wf1 && wf3!=wf2)   wf3->start(wf3->start()-tmin);
  if(wref!=wf1 && wref!=wf2  && wref!=wf3) wref->start(wref->start()-tmin);

  if(fft) {
    PTS.plot(wf1, const_cast<char*>("AL"), kRed, 0, 0, true, cfg->fLow, cfg->fHigh);
  } else {
    double bT, eT;
    GetTimeBoundaries(*wref, P, bT, eT);
    PTS.plot(wf1, const_cast<char*>("AL"), kRed, bT, eT);
    PTS.graph[0]->GetXaxis()->SetRangeUser(bT, eT);
  }
  PTS.graph[0]->SetLineWidth(1);
  PTS.graph[0]->SetTitle(title);

  TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",tmin);
  PTS.graph[0]->GetXaxis()->SetTitle(xtitle);

  if(wf2!=NULL) {
    if(fft) {
      if(wf2!=NULL) PTS.plot(wf2, const_cast<char*>("SAME"), kBlue, 0, 0, true, cfg->fLow, cfg->fHigh);
    } else {
      if(wf2!=NULL) PTS.plot(wf2, const_cast<char*>("SAME"), kBlue, 0, 0);
    }
    if(wf2!=NULL) PTS.graph[1]->SetLineWidth(1);
  }

  if(wf3!=NULL) {
    if(fft) {
      if(wf3!=NULL) PTS.plot(wf3, const_cast<char*>("SAME"), kGreen+2, 0, 0, true, cfg->fLow, cfg->fHigh);
    } else {
      if(wf3!=NULL) PTS.plot(wf3, const_cast<char*>("SAME"), kGreen+2, 0, 0);
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

void SetOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, bool dump_ebbh) {

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
       cout << "CWB_Plugin_ChirpMass.C : Error - output root file not found" << endl;
       gSystem->Exit(1);                                                                             
     }                                                                                               
   } else {                                                                                          
     cout << "CWB_Plugin_ChirpMass.C : Error - output root file not found" << endl;                   
     gSystem->Exit(1);                                                                               
   }                                                                                                 

   gTREE = (TTree *) froot->Get("waveburst");
   if(gTREE!=NULL) {                         
     EVT = new netevent(gTREE,nIFO);         
   } else {                                                                      
     EVT = new netevent(nIFO);                                                   
     gTREE = EVT->setTree();                                                  
   }                                                                                        
   EVT->setSLags(gSLAGSHIFT);                          // set slags into netevent           
   EVT->Psave=cfg->Psave;                                                                   

   static bool ebbh_tree_init = false;
   if(dump_ebbh) { 

     gEBBH.segGPS = gSEGGPS;
     if(gOPT.dump_wave) {
       if(cfg->simulation) for(int n=0;n<nIFO;n++) gEBBH.wINJ[n] = &gwINJ[n];
       for(int n=0;n<nIFO;n++) gEBBH.wDAT[n] = &gwDAT[n];
       for(int n=0;n<nIFO;n++) gEBBH.wREC[n] = &gwREC[n];
       for(int n=0;n<nIFO;n++) gEBBH.sREC[n] = &gsREC[n];
     }
     if(gOPT.dump_cluster) {
       gEBBH.cRATE = gCRATE;
       gEBBH.pCLUSTER = &gCLUSTER;
     }
   
     if(!ebbh_tree_init) { 

       ebbh_tree_init = true;

       gTREE->Branch("ebbh_ei", &gEBBH.ei,"ebbh_ei/F");
       gTREE->Branch("ebbh_nch",&gEBBH.nch,"ebbh_nch/I");

       gTREE->Branch("ebbh_ch_mass",   gEBBH.ch_mass,TString::Format("ebbh_ch_mass[%i]/F",NCHIRP_MAX));
       gTREE->Branch("ebbh_ch_tmerger",gEBBH.ch_tmerger,TString::Format("ebbh_ch_tmerger[%i]/F",NCHIRP_MAX));
       gTREE->Branch("ebbh_ch_merr",   gEBBH.ch_merr,TString::Format("ebbh_ch_merr[%i]/F",NCHIRP_MAX));
       gTREE->Branch("ebbh_ch_mchi2",  gEBBH.ch_mchi2,TString::Format("ebbh_ch_mchi2[%i]/F",NCHIRP_MAX));
       gTREE->Branch("ebbh_ch_energy", gEBBH.ch_energy,TString::Format("ebbh_ch_energy[%i]/F",NCHIRP_MAX));

       gTREE->Branch("ebbh_seggps",&gEBBH.segGPS,"ebbh_seggps/D");
       if(gOPT.dump_wave) {
         for(int n=0;n<nIFO;n++) {
           if(cfg->simulation) gTREE->Branch(TString::Format("ebbh_wINJ_%d",n).Data(),"wavearray<double>",&gEBBH.wINJ[n],32000,0);
           gTREE->Branch(TString::Format("ebbh_wDAT_%d",n).Data(),"wavearray<double>",&gEBBH.wDAT[n],32000,0);
           gTREE->Branch(TString::Format("ebbh_wREC_%d",n).Data(),"wavearray<double>",&gEBBH.wREC[n],32000,0);
           gTREE->Branch(TString::Format("ebbh_sREC_%d",n).Data(),"wavearray<double>",&gEBBH.sREC[n],32000,0);
         }
       }
       if(gOPT.dump_cluster) {
         gTREE->Branch("ebbh_crate",&gEBBH.cRATE,"ebbh_crate/F");
         gTREE->Branch("ebbh_cluster","vector<netpixel>",&gEBBH.pCLUSTER,32000,0);
       }

     } else {

       gTREE->SetBranchAddress("ebbh_ei", &gEBBH.ei);
       gTREE->SetBranchAddress("ebbh_nch",&gEBBH.nch);

       gTREE->SetBranchAddress("ebbh_ch_mass",   gEBBH.ch_mass);
       gTREE->SetBranchAddress("ebbh_ch_tmerger",gEBBH.ch_tmerger);
       gTREE->SetBranchAddress("ebbh_ch_merr",   gEBBH.ch_merr);
       gTREE->SetBranchAddress("ebbh_ch_mchi2",  gEBBH.ch_mchi2);
       gTREE->SetBranchAddress("ebbh_ch_energy", gEBBH.ch_energy);

       gTREE->SetBranchAddress("ebbh_seggps",&gEBBH.segGPS);
       if(gOPT.dump_wave) {
         for(int n=0;n<nIFO;n++) {
           if(cfg->simulation) gTREE->SetBranchAddress(TString::Format("ebbh_wINJ_%d",n).Data(),&gEBBH.wINJ[n]);
           gTREE->SetBranchAddress(TString::Format("ebbh_wDAT_%d",n).Data(),&gEBBH.wDAT[n]);
           gTREE->SetBranchAddress(TString::Format("ebbh_wREC_%d",n).Data(),&gEBBH.wREC[n]);
           gTREE->SetBranchAddress(TString::Format("ebbh_sREC_%d",n).Data(),&gEBBH.sREC[n]);
         }
       }
       if(gOPT.dump_cluster) {
         gTREE->SetBranchAddress("ebbh_crate",&gEBBH.cRATE);
         gTREE->SetBranchAddress("ebbh_cluster",&gEBBH.pCLUSTER);
       }
     }
   }                                                                                        
}                                                                                           

void DumpOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor) {

   if(cfg->dump) EVT->dopen(gOUTPUT.Data(),const_cast<char*>("a"),false);        
   EVT->output2G(gTREE,NET,ID,k,factor);             // get reconstructed parameters
   if(cfg->dump) EVT->dclose();
   if(!cfg->cedDump) NET->getwc(k)->sCuts[ID-1]=1;   // mark as processed
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
  for(int j=0;j<size/2;j+=2) {
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

std::vector<netpixel> GetCluster(network* NET, CWB::config* cfg, int lag, int id) {

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

void GetEccentricityIndex(network* NET, int lag, int id) {

  netcluster* pwc = NET->getwc(lag);

  double chi2_thr    = gOPT.chi2_thr;
  double tmerger_cut = gOPT.tmerger_cut;
  double zmax_thr    = gOPT.zmax_thr;

  clusterdata CD = pwc->cData[id-1]; 	// save loudest chirp

  for(int i=0;i<NCHIRP_MAX;i++) {

    double echirp = pwc->mchirp(id,chi2_thr,tmerger_cut,zmax_thr);

    clusterdata* pcd = &(pwc->cData[id-1]);
    TF1* fit = &(pcd->fit);

    gEBBH.ch_mass[i]    = pcd->mchirp;
    gEBBH.ch_tmerger[i] = -fit->GetParameter(1)/fit->GetParameter(0);
    gEBBH.ch_merr[i]    = pcd->mchirperr;
    gEBBH.ch_mchi2[i]   = pcd->chi2chirp;
    gEBBH.ch_energy[i]  = echirp;
  }

  // compute the eccentrycity index
  int imax=0;
  float emax=0;
  // find chirp with max energy
  for(int i=0;i<NCHIRP_MAX;i++) if(gEBBH.ch_energy[i]>emax) {emax=gEBBH.ch_energy[i];imax=i;}
  int imax2=imax;
  // find chirp with highest chirp mass and energy > (max chirp energy)/4
  for(int i=0;i<NCHIRP_MAX;i++) {
    if(i==imax) continue;
    if(gEBBH.ch_mass[i]>0) {
      if(gEBBH.ch_mass[i]>gEBBH.ch_mass[imax]) if(gEBBH.ch_energy[i]>gEBBH.ch_energy[imax]/4.) imax2=i;
    }
  }
  imax=imax2;
  gEBBH.ei=0;
  for(int i=0;i<NCHIRP_MAX;i++) {
    if(gEBBH.ch_mass[i]>0) {
      if(gEBBH.ch_mass[i]<gEBBH.ch_mass[imax]) gEBBH.ei+=gEBBH.ch_energy[i];
      if(gEBBH.ch_mass[i]>gEBBH.ch_mass[imax]) gEBBH.ei-=gEBBH.ch_energy[i];
    }
  }
  gEBBH.ei/=gEBBH.ch_energy[imax];
  gEBBH.ei*=100;

  cout << endl;
  cout << "---------------------------------------------------" << endl;
  for(int i=0;i<NCHIRP_MAX;i++) cout << "mchirp " << i << " -> " << gEBBH.ch_mass[i] << "\t\t echirp -> " << gEBBH.ch_energy[i] << endl;
  cout << "Eccentricity Index : " << gEBBH.ei << endl;
  cout << "---------------------------------------------------" << endl;
  cout << endl;

  pwc->cData[id-1] = CD;		// restore loudest chirp
}

void ReadUserOptions(CWB::config* cfg) {

  TString options = cfg->parPlugin;

  gOPT.chi2_thr     = CHI2_THR;
  gOPT.tmerger_cut  = TMERGER_CUT;
  gOPT.zmax_thr     = ZMAX_THR;

  gOPT.dump_wave    = DUMP_WAVE;
  gOPT.dump_cluster = DUMP_CLUSTER;

  if(options.CompareTo("")!=0) {
    cout << options << endl;
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet

      TObjArray* token = TString(options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("ebbh_chi2_thr=")) {
          TString ebbh_chi2_thr=stok;
          ebbh_chi2_thr.Remove(0,ebbh_chi2_thr.Last('=')+1);
          if(ebbh_chi2_thr.IsFloat()) gOPT.chi2_thr=ebbh_chi2_thr.Atof();
        }

        if(stok.Contains("ebbh_tmerger_cut=")) {
          TString ebbh_tmerger_cut=stok;
          ebbh_tmerger_cut.Remove(0,ebbh_tmerger_cut.Last('=')+1);
          if(ebbh_tmerger_cut.IsFloat()) gOPT.tmerger_cut=ebbh_tmerger_cut.Atof();
        }

        if(stok.Contains("ebbh_zmax_thr=")) {
          TString ebbh_zmax_thr=stok;
          ebbh_zmax_thr.Remove(0,ebbh_zmax_thr.Last('=')+1);
          if(ebbh_zmax_thr.IsFloat()) gOPT.zmax_thr=ebbh_zmax_thr.Atof();
        }

        if(stok.Contains("ebbh_dump_wave=")) {
          TString opt=stok;
          opt.Remove(0,opt.Last('=')+1);
          if(opt=="false")  gOPT.dump_wave=false;
          if(opt=="true")   gOPT.dump_wave=true;
        }

        if(stok.Contains("ebbh_dump_cluster=")) {
          TString opt=stok;
          opt.Remove(0,opt.Last('=')+1);
          if(opt=="false")  gOPT.dump_cluster=false;
          if(opt=="true")   gOPT.dump_cluster=true;
        }
      }
    }
  }
}

void PrintUserOptions() {

  cout << "-----------------------------------------"     << endl;
  cout << "eBBH Plugin config options" << endl; 
  cout << "-----------------------------------------"     << endl << endl;
  cout << endl;
  cout << "chi2_thr      : " << gOPT.chi2_thr     << endl;
  cout << "tmerger_cut   : " << gOPT.tmerger_cut  << endl;
  cout << "zmax_thr      : " << gOPT.zmax_thr     << endl;
  cout << "dump_wave     : " << gOPT.dump_wave    << endl;
  cout << "dump_cluster  : " << gOPT.dump_cluster << endl;
  cout << endl;

}

