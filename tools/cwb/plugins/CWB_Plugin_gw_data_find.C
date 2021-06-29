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


//!NOISE_MDC_SIMULATION
// Plugin used to retrive 'on the fly' the frame file paths using gw_data_find

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
#include "cwb2G.hh"
#include "mdc.hh"
#include <vector>

double Tb = 0;
double Te = 0;

#define FRAME_SCRATCH_TIME	10000	// seconds before and after the segment range

void
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  if(type==CWB_PLUGIN_CONFIG) {
    cfg->mdcPlugin=true;         // disable read mdc from frames
    cfg->dataPlugin=true;        // disable read data from frames
  }

  if(type==CWB_PLUGIN_INIT_JOB) {

    if(Tb!=0 || Te!=0) return;

    cout << "Execute CWB_Plugin_gw_data_find.C : Inject On The Fly MDC ..." << endl;

    // get data range
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
    Tb = gCWB2G->Tb-FRAME_SCRATCH_TIME;
    Te = gCWB2G->Te+FRAME_SCRATCH_TIME;
    //cout << "CWB_Plugin_gw_data_find.C - " << "Tb : " << int(Tb) << " Te : " << int(Te) << endl; 


    for(int n=0;n<cfg->nIFO;n++) {

       // get observatory
       TString observatory="";
       if(TString(net->getifo(n)->Name)=="L1") observatory="L";
       if(TString(net->getifo(n)->Name)=="H1") observatory="H";
       if(TString(net->getifo(n)->Name)=="V1") observatory="V";
       cout << "observatory : " << observatory << " type : " << type << endl; 
       if(observatory=="") {cout << "CWB_Plugin_gw_data_find.C : observatory not defined !!!" << endl;gSystem->Exit(1);}

       // get type
       TString type = cfg->frFiles[n];
       if(type=="") {cout << "CWB_Plugin_gw_data_find.C : type not defined !!!" << endl;gSystem->Exit(1);}

       // build gw_data_find command
       char cmd[1024];
       sprintf(cmd,"gw_data_find --observatory %s --type %s --gps-start-time %g --gps-end-time %g --url file > %s/%s_gw_data_find.frames",
               observatory.Data(), type.Data(), Tb, Te, cfg->input_dir, net->getifo(n)->Name); 

       // create frame files list
       cout << cmd << endl;
       gSystem->Exec(cmd);

       // set frame file list in the configuration object 
       sprintf(cfg->frFiles[n],"%s/%s_gw_data_find.frames", cfg->input_dir, net->getifo(n)->Name);

       // cout << "CWB_Plugin_gw_data_find.C : frFiles : " << cfg->frFiles[n] << endl;
    } 

    cfg->mdcPlugin=false;         // enable read mdc from frames
    cfg->dataPlugin=false;        // enable read data from frames

    gCWB2G->InitJob(); 
  }
}
