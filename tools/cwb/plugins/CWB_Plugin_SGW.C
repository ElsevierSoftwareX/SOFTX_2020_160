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
#include <vector>


// ---------------------------------------------------------------------------------
// DEFINES
// ---------------------------------------------------------------------------------

#define SGW_MDC_POLARIZATION	"SCALAR"
#define SGW_SEARCH_POLARIZATION	"SCALAR"

// ---------------------------------------------------------------------------------
// USER SGW PLUGIN OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  TString mdc_polarization;	// SCALAR/TENSOR polarizations used for MDC
  TString search_polarization;	// SCALAR/TENSOR polarizations used for SEARCH
};

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions gOPT;                           // global User Options

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  cout << endl;
  cout << "-----> macro/CWB_Plugin_SGW.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;


  if(type==CWB_PLUGIN_CONFIG) {  
    ResetUserOptions();                                 // set default config options
    ReadUserOptions(cfg->parPlugin);                    // user config options : read from parPlugin

    if(gOPT.search_polarization=="SCALAR") cfg->gamma = -1.0; // force hard constrain for SGW 
  }

  if(type==CWB_PLUGIN_NETWORK) {  
    // set polarization for MDC
    PrintUserOptions(cfg);      // print config options
    for(int n=0;n<(int)net->ifoListSize();n++) {
      detector *d = net->getifo(n);
      d->setPolarization(gOPT.mdc_polarization=="SCALAR" ? SCALAR : TENSOR);
    }
    // restore skymaps
    if(cfg->healpix) net->setSkyMaps(int(cfg->healpix));
    else             net->setSkyMaps(cfg->angle,cfg->Theta1,cfg->Theta2,cfg->Phi1,cfg->Phi2);
    net->setAntenna();
    net->setDelay(cfg->refIFO);
  }


  if(type==CWB_PLUGIN_ICOHERENCE) {  
    // set polarization for SEARCH
    PrintUserOptions(cfg);      // print config options
    for(int n=0;n<(int)net->ifoListSize();n++) {
      detector *d = net->getifo(n);
      d->setPolarization(gOPT.search_polarization=="SCALAR" ? SCALAR : TENSOR);
    }
    // restore skymaps
    if(cfg->healpix) net->setSkyMaps(int(cfg->healpix));
    else             net->setSkyMaps(cfg->angle,cfg->Theta1,cfg->Theta2,cfg->Phi1,cfg->Phi2);
    net->setAntenna();
    net->setDelay(cfg->refIFO);
  }

  return;
}

void ReadUserOptions(TString options) {

  // get plugin options 

  if(TString(options)!="") {

    //cout << "SGW options : " << options << endl;
    TObjArray* token = TString(options).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++) {

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();

      if(stok.Contains("--sgw_mdc_polarization=")) {
        gOPT.mdc_polarization=stok;
        gOPT.mdc_polarization.Remove(0,gOPT.mdc_polarization.Last('=')+1);
        gOPT.mdc_polarization.ToUpper();;
      }
      if(gOPT.mdc_polarization!="TENSOR" && gOPT.mdc_polarization!="SCALAR") {
        cout << "CWB_Plugin_SGW Error : parameter sgw_mdc_polarization must be : SCALAR/TENSOR" << endl;
        gSystem->Exit(1);
      }

      if(stok.Contains("--sgw_search_polarization=")) {
        gOPT.search_polarization=stok;
        gOPT.search_polarization.Remove(0,gOPT.search_polarization.Last('=')+1);
        gOPT.search_polarization.ToUpper();;
      }
      if(gOPT.search_polarization!="TENSOR" && gOPT.search_polarization!="SCALAR") {
        cout << "CWB_Plugin_SGW Error : parameter sgw_search_polarization must be : SCALAR/TENSOR" << endl;
        gSystem->Exit(1);
      }

    }
  }
}

void ResetUserOptions() {

    gOPT.mdc_polarization     = SGW_MDC_POLARIZATION;
    gOPT.search_polarization  = SGW_SEARCH_POLARIZATION;
} 


void PrintUserOptions(CWB::config* cfg) {

    cout << "-----------------------------------------"     << endl;
    cout << "SGW config options                       "     << endl;
    cout << "-----------------------------------------"     << endl << endl;
    cout << "SGW_MDC_POLARIZATION         " << gOPT.mdc_polarization    << endl;
    cout << "SGW_SEARCH_POLARIZATION      " << gOPT.search_polarization << endl;
    cout << endl;
}

