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

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!DATA_CONDITIONING
// Implements the 2G frequency cuts after the whitening 

  cout << endl;
  cout << "-----> CWB_Plugin_WDM_freqCuts.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {
    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_freqCuts.C is implmenented only for 2G analysis" << endl;
      gSystem->Exit(1);
    }

    // cuts around the H1 S6A upconversion frequency ranges
    int nfreqCuts = 6;
    double freqCuts[100][2];
    freqCuts[0][0]=58; 		freqCuts[0][1]=62;
    freqCuts[1][0]=78; 		freqCuts[1][1]=81;
    freqCuts[2][0]=118; 	freqCuts[2][1]=122;
    freqCuts[3][0]=176; 	freqCuts[3][1]=184;
    freqCuts[4][0]=188; 	freqCuts[4][1]=191;
    freqCuts[5][0]=198; 	freqCuts[5][1]=201;

    // select max resolution level (the one used for whitening)
    int layers = 0;
    for(int level=cfg->l_high; level>=cfg->l_low; level--) {
      for(int j=0;j<net->wdmMRA.nRes;j++)
         if(layers<net->wdmMRA.layers[j]) layers=net->wdmMRA.layers[j];
    }
    WDM<double> WDMwhite(layers,layers,6,10);         

    // WDM transform
    WSeries<double> wM;             
    wM.Forward(*x,WDMwhite);

    wavearray<double> z; 
    layers  = wM.maxLayer();
    for(int j=0;j<layers;j++) {
      double freq = wM.frequency(j);
      double layer = j+0.01;		// -epsilon select 0 layer for 90 phase
      for(int k=0;k<nfreqCuts;k++) {
        // apply frequency cuts
        if((freq>freqCuts[k][0])&&(freq<freqCuts[k][1])) { 
          wM.getLayer(z, layer);z=0;wM.putLayer(z, layer);      //  0 phase
          wM.getLayer(z,-layer);z=0;wM.putLayer(z,-layer);      // 90 phase
        }
      }
    }

    // inverse WDM transform
    wavearray<double>* y = (wavearray<double>*)x;
    WSeries<double> wM2 = wM;
    wM.Inverse();
    wM2.Inverse(-2);
    *y  = wM;
    *y += wM2;
    *y *= 0.5;
  }

  return;
}
