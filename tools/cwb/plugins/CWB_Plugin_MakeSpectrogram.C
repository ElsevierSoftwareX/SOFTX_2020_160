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
#include "Toolbox.hh"
#include "STFT.hh"

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!DISPLAY_SPECTRA
// Plugin to produce scalograms of input data

  cout << endl;
  cout << "-----> CWB_Plugin_MakeSpectrogram.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

//  CWB_PLUGIN stage = CWB_PLUGIN_STRAIN;
//  CWB_PLUGIN stage = CWB_PLUGIN_MDC;
  CWB_PLUGIN stage = CWB_PLUGIN_STRAIN_AND_MDC;
//  CWB_PLUGIN stage = CWB_PLUGIN_WHITE;

  if(type==stage) {
     TString data_type="unknown";
     if(stage==CWB_PLUGIN_STRAIN) data_type="strain";
     if(stage==CWB_PLUGIN_MDC) data_type="mdc";
     if(stage==CWB_PLUGIN_STRAIN_AND_MDC) data_type="strain_and_mdc";
     if(stage==CWB_PLUGIN_WHITE) data_type="white";
 
     bool isWavearray=false; 
     if(TString(cfg->analysis)=="1G") {
       if(stage==CWB_PLUGIN_STRAIN)         isWavearray=true;
       if(stage==CWB_PLUGIN_MDC)            isWavearray=true;
       if(stage==CWB_PLUGIN_STRAIN_AND_MDC) isWavearray=false;
       if(stage==CWB_PLUGIN_WHITE)          isWavearray=false;
     }
     if(TString(cfg->analysis)=="2G") {
       if(stage==CWB_PLUGIN_STRAIN)         isWavearray=true;
       if(stage==CWB_PLUGIN_MDC)            isWavearray=true;
       if(stage==CWB_PLUGIN_STRAIN_AND_MDC) isWavearray=true;
       if(stage==CWB_PLUGIN_WHITE)          isWavearray=true;
     }
 
     int level = 0;
     if(isWavearray) { // x is a wavearray 
       cout << "Input data is a wavearray" << endl;
     } else { // x is wseries
       cout << "Input data is a wseries" << endl;
       level = x->getLevel();
       x->Inverse();
     }

     cout << "Decomposition Level : " << level << endl;

     int nfact=4;
     int nfft=nfact*512;
     int noverlap=nfft-10;
     double fparm=nfact*6;

     CWB::STFT stft(*x,nfft,noverlap,"energy","gauss",fparm);

     double start = x->start()+cfg->segEdge;
     double stop  = x->start()+x->size()/x->rate()-cfg->segEdge;
     double flow  = 64;
     double fhigh = 2048;
     stft.Draw(start,stop,flow,fhigh,0,0,1);
     // dump spectrum
     char fname[1024];
     sprintf(fname,"%s/spectrogram_%s_%s_%d_%s_job%lu.root",cfg->dump_dir,data_type.Data(),ifo.Data(),
             int(x->start()),cfg->data_label,net->nRun);
     cout << endl << "Dump Spectrogram : " << fname << endl << endl;
     stft.Print(fname);

     int nIFO=net->ifoListSize();
     if(level) x->Forward(level);
     if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) gSystem->Exit(0);  // last ifo
  }

  return;
}
