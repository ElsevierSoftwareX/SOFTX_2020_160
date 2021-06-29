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
#include "TMath.h"
#include "mdc.hh"
#include "regression.hh"
#include <vector>


#define REGRESSION_FILTER_LENGTH        8
#define REGRESSION_MATRIX_FRACTION      0.95
#define REGRESSION_SOLVE_EIGEN_THR      0.
#define REGRESSION_SOLVE_EIGEN_NUM      10
#define REGRESSION_SOLVE_REGULATOR      'h'
#define REGRESSION_APPLY_THR            0.8

#define LPR_LEVELD	8
#define LPR_LEVELF	6
#define LPR_TLPR	120.

#define USE_REGRESSION
//#define USE_LPR_FILTER

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!DATA_CONDITIONING
// This plugin shows how to implement a custom Data Conditioning (only 2G)

  cout << endl;
  cout << "-----> CWB_Plugin_DataConditioning.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->dcPlugin=true;   			// disable built in dc 
  }

  if(type==CWB_PLUGIN_DATA_CONDITIONING) {  

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_DAtaConditioning.C -> implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    wavearray<double> y;
    WSeries<double> wM;
    // import global variables
    size_t gRATEANA; IMPORT(size_t,gRATEANA)
    void* gMDC=NULL; IMPORT(void*,gMDC)
    WSeries<double>* xM = (WSeries<double>*)gMDC;
#ifdef USE_LPR_FILTER				      // used in 1G
    Meyer<double> S(1024,2);                          // set wavelet for production
#endif

    // get ifo id
    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    detector* pD = net->getifo(id);
    WSeries<double>* pTF = pD->getTFmap();
    wavearray<double>* hot = (wavearray<double>*)x;

    int layers = 0;
    for(int level=cfg->l_high; level>=cfg->l_low; level--) {
      for(int j=0;j<net->wdmMRA.nRes;j++)
         if(layers<net->wdmMRA.layers[j]) layers=net->wdmMRA.layers[j];
    }
    WDM<double> WDMwhite(layers,layers,6,10);         // set whitening WDM
    layers = gRATEANA/8;
    WDM<double> WDMlpr(layers,layers,6,10);           // set LPE filter WDM

#ifdef USE_LPR_FILTER				      // used in 1G
    pTF->Forward(*hot,S,LPR_LEVELD);
    pTF->lprFilter(2,0,LPR_TLPR,4.);
    pTF->Inverse();
    *hot = *pTF;      		                      // conditioned TS
#endif

#ifdef USE_REGRESSION				      // used in 2G
    // regression
    pTF->Forward(*hot,WDMlpr);
    regression rr(*pTF,const_cast<char*>("target"),1.,cfg->fHigh);
    rr.add(*hot,const_cast<char*>("target"));                    
    rr.setFilter(REGRESSION_FILTER_LENGTH);                         
    rr.setMatrix(net->Edge,REGRESSION_MATRIX_FRACTION);              
    rr.solve(REGRESSION_SOLVE_EIGEN_THR,REGRESSION_SOLVE_EIGEN_NUM,REGRESSION_SOLVE_REGULATOR);
    rr.apply(REGRESSION_APPLY_THR);                                                            
    *hot = rr.getClean();                                                                   
#endif

    // whitening
    pTF->Forward(*hot,WDMwhite);            
    pTF->setlow(cfg->fLow);                     
    pTF->sethigh(cfg->fHigh);                   
    pD->white(cfg->whiteWindow,0,cfg->segEdge,cfg->whiteStride);      	// calculate noise rms 
    pD->nRMS.bandpass(16.,0.,1);                                        // high pass filtering at 16Hz
    pTF->white(pD->nRMS,1);                                     	// whiten  0 phase WSeries
    pTF->white(pD->nRMS,-1);                                    	// whiten 90 phase WSeries
    if(cfg->simulation) {                                              	// estimated whitened MDC parms
      wM.Forward(*xM,WDMwhite);                                                                       
      wM.setlow(cfg->fLow);                                                                           
      wM.sethigh(cfg->fHigh);                                                                         
      // zero f<fLow to avoid whitening issues when psd noise is not well defined for f<fLow         
      int layers  = wM.maxLayer();                                                                   
      for(int j=0;j<layers;j++) if(wM.frequency(j)<cfg->fLow) {
        double layer = j+0.01;                                // -epsilon select 0 layer for 90 phase
        wM.getLayer(y, layer);y=0;wM.putLayer(y, layer);      //  0 phase
        wM.getLayer(y,-layer);y=0;wM.putLayer(y,-layer);      // 90 phase
      }
      // compute mdc snr and save whiten waveforms
      pD->setsim(wM,net->getmdcTime(),cfg->iwindow/2.,cfg->segEdge,true);
    }
    WSeries<double> wtmp = *pTF;
    pTF->Inverse();
    wtmp.Inverse(-2);
    *hot = *pTF;
    *hot += wtmp;
    *hot *= 0.5;

#ifdef USE_LPR_FILTER				              // used in 1G
    pTF->Forward(*hot,S,LPR_LEVELF);
    pTF->lprFilter(2,0,LPR_TLPR,4.);
    pTF->Inverse();
    *hot = *pTF;                                	      // conditioned TS
#endif

  }
  return;
}
