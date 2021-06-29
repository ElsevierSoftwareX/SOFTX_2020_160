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

// ---------------------------------------------------------------------------------
// HOW TO SETUP PLUGIN SimNoise IN CWB USER CONFIGURATION (EXAMPLE)
// ---------------------------------------------------------------------------------

/*
  TString optSN = "";                           // NOTE : add space at the end of each line

  optSN += "sn_strain_file_name=LIGO.txt ";     // add LIGO strain file name to first detector
  optSN += "sn_strain_file_name=Virgo.txt ";    // add Virgo strain file name to second detector

  optSN += "sn_seed=1 ";                        // add seed for random noise generation to first detector
  optSN += "sn_seed=2 ";                        // add seed for random noise generation to second detector

  strcpy(parPlugin,optSN.Data());               // set WF plugin parameters
*/

// ---------------------------------------------------------------------------------
// DEFINES
// ---------------------------------------------------------------------------------

#define SN_STRAIN_FILE_NAME    "" 
#define SN_SEED                150914

// ---------------------------------------------------------------------------------
// USER PLUGIN OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  TString  strain_file_name[NIFO_MAX];
  int  	   seed[NIFO_MAX];
};

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions        gOPT;                           // global User Options

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);

 
void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// Plugin to generate simulated gaussian noise

  cout << endl;
  cout << "-----> CWB_Plugin_SimNoise.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;


  if(type==CWB_PLUGIN_CONFIG) {
    cfg->dataPlugin=true; 				// disable read data from frames

    ResetUserOptions();                                 // set default config options
    ReadUserOptions(cfg->parPlugin);                    // user config options : read from parPlugin
  }

  if(type==CWB_PLUGIN_NETWORK) {
    PrintUserOptions(cfg);      // print config options
    // check user options
    for(int n=0;n<cfg->nIFO;n++) {
      for(int m=n+1;m<cfg->nIFO;m++) {
        if(gOPT.seed[n]==gOPT.seed[m]) {
          cout << "CWB_Plugin_SimNoise Error : seed must me unique !!!!" << endl;
          cout << endl;
          exit(1);
        }
      }
    }
    for(int n=0;n<cfg->nIFO;n++) CWB::Toolbox::checkFile(gOPT.strain_file_name[n]);
  }

  if(type==CWB_PLUGIN_DATA) {

    CWB::Toolbox TB;

    int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {ifoID=n;break;}

    int seed = gOPT.seed[ifoID];
    TString fName = gOPT.strain_file_name[ifoID];;
    cout << "CWB_Plugin_SimNoise " << net->getifo(ifoID)->Name << " -> seed " << seed << " - " << fName << endl;  

    int size=x->size();
    double start=x->start();
    TB.getSimNoise(*x, fName, seed, net->nRun);
    x->resize(size);
    x->start(start);
  }
  return;
}

void PrintUserOptions(CWB::config* cfg) {

    cout << "-----------------------------------------"     << endl;
    cout << "SN config options                        "     << endl;
    cout << "-----------------------------------------"     << endl << endl;

    for(int n=0;n<cfg->nIFO;n++) {
      cout << "SN_STRAIN_FILE_NAME  " << cfg->ifo[n] << "   " << gOPT.strain_file_name[n] << endl;
      cout << "SN_SEED              " << cfg->ifo[n] << "   " << gOPT.seed[n] << endl;
    }

    cout << endl;
}


void ReadUserOptions(TString options) {

  int n_strain_file_name=0;
  int n_seed=0;
  if(options.CompareTo("")!=0) {
    cout << options << endl;
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet

      TObjArray* token = TString(options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("sn_strain_file_name=")) {
          TString sn_strain_file_name=stok;
          sn_strain_file_name.Remove(0,sn_strain_file_name.Last('=')+1);
          gOPT.strain_file_name[n_strain_file_name]=sn_strain_file_name;
          if(n_strain_file_name<(NIFO_MAX-1)) n_strain_file_name++;
        }

        if(stok.Contains("sn_seed=")) {
          TString sn_seed=stok;
          sn_seed.Remove(0,sn_seed.Last('=')+1);
          if(sn_seed.IsDigit()) gOPT.seed[n_seed]=sn_seed.Atoi();
          if(n_seed<(NIFO_MAX-1)) n_seed++;
        }

      }
    }
  }
}


void ResetUserOptions() {

  for(int n=0;n<NIFO_MAX;n++) gOPT.strain_file_name[n] = SN_STRAIN_FILE_NAME;
  for(int n=0;n<NIFO_MAX;n++) gOPT.seed[n] = SN_SEED+n;
  
}
