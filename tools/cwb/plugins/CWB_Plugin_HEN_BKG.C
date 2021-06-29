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
#include <vector>


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!NOISE_MDC_SIMULATION
// Plugin to set SkyMaskCC 'on the fly' 

  if(type==CWB_PLUGIN_STRAIN) {  

    // ---------------------------------
    // read plugin config 
    // CWB_Plugin_HEN_BKG_Config.C
    // ---------------------------------

    if(TString(ifo).CompareTo(net->ifoName[0])==0) {

      char cmd[128]; 
      sprintf(cmd,"network* net = (network*)%p;",net);
      gROOT->ProcessLine(cmd);
      sprintf(cmd,"CWB::config* cfg = (CWB::config*)%p;",cfg);
      gROOT->ProcessLine(cmd);
      sprintf(cmd,"int xstart = %d;",(int)x->start());
      gROOT->ProcessLine(cmd);
      sprintf(cmd,"int xstop = %d;",(int)x->stop());
      gROOT->ProcessLine(cmd);
      int error=0;
      cfg->configPlugin.Exec(NULL,&error);
      if(error) {
        cout << "Error executing macro : " << cfg->configPlugin.GetTitle() << endl;
        cfg->configPlugin.Print(); 
        gSystem->Exit(1); 
      }
    }
  }

  return;
}
