/*
# Copyright (C) 2019 Sergey Klimenko
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
#include "cwb2G.hh"
#include "config.hh"
#include "network.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "Toolbox.hh"

// ---------------------------------------------------------------------------------
// HOW TO SETUP PLUGIN IN CWB USER CONFIGURATION (EXAMPLE)
// WORNING: This plugin works only with 3 detectors !!!
// ---------------------------------------------------------------------------------

/*
  TString opt_ulags = "";                  // NOTE : add space at the end of each line
  opt_ulags += "ulags_enabled=true ";      // true/false enable/disable unique lags (default=false)
  opt_ulags += "ulags_lagSize=599 ";       // if>0 force lagSize to be fixed
  strcpy(parPlugin,opt_ulags.Data());      // set plugin parameters
  lagSize = 0;				   // mandatory, used to manage command 'cwb_inet #job 0 ced 0 #lag' 
  lagStep = 1.;
  lagOff  = 0;
*/

// ---------------------------------------------------------------------------------
// USER UniqueLags PLUGIN OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  bool   enabled;	// true/false enable/disable unique lags (default=false)		
  int    lagSize; 	// if>0 force lagSize to be fixed
};

uoptions gOPT;          // global User Options
bool     gULAG;         // used to manage command 'cwb_inet #job 0 ced 0 #lag'

void ResetUserOptions(network* net, CWB::config* cfg);
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!CONFIG
// Plugin to setup dinamically the unique lags for 2/3/4 ifos network

  if(type==CWB_PLUGIN_CONFIG) {
    gULAG = (cfg->lagSize) ? true : false;      // single lag is select by the command 'cwb_inet #job 0 ced 0 #lag'
  }

  if(type==CWB_PLUGIN_NETWORK) {
    ResetUserOptions(net, cfg);                 // set default config options
    ReadUserOptions(cfg->parPlugin);            // user config options : read from parPlugin
    PrintUserOptions(cfg);                      // print config options
  }

  if(type==CWB_PLUGIN_INIT_JOB) {

    if(!gOPT.enabled) return;

    cout << endl;
    cout << "-----> CWB_Plugin_UniqueLags.C ";
    cout << "type " << type << endl;
    cout << endl;

    // get ifo id
    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) {
       if(ifo==net->ifoName[n]) {id=n; break;}
    }
    if(id>0) return;				// plugin is execute only for the first ifo

    if(cfg->simulation) {cout << "CWB_Plugin_UniqueLags.C : Error - works only for background" << endl; gSystem->Exit(1);} 
    if(cfg->nIFO<2 || cfg->nIFO>4) {cout << "CWB_Plugin_UniqueLags.C : Error - implemented only for 2/3/4 ifos network" << endl; gSystem->Exit(1);}

    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
    // get lagBuffer (used to store the contents of the user lagFile)
    TArrayC* lagBuffer = &(gCWB2G->lagBuffer);
    // get data range
    double segLen = gCWB2G->GetSegEnd()-gCWB2G->GetSegBegin();
    // get laxMax & lagSize
    int lagMax = int(segLen/cfg->lagStep);
    if(cfg->nIFO==2) {
      lagMax = lagMax-lagMax%2;	// force to be an even number
    }
    if(cfg->nIFO==3) {
      lagMax = lagMax-lagMax%3;
      lagMax = lagMax-lagMax%2;	// force to be an even number
    }
    int lagSize = lagMax-1;

    if(gOPT.lagSize>0) {	// if>0 force lagSize to be fixed
      lagSize = gOPT.lagSize;
      lagMax  = lagSize+1;
    } 

    cout << endl << "CWB_Plugin_UniqueLags.C : segLen -> " << segLen << " lagMax " << lagMax << " lagSize " << lagSize << endl << endl;

    // redirect cout to buffer and return contents
    std::stringstream buffer;
    std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());

    // generate unique lags
    if(cfg->nIFO==2) {
      cout << 0 << "\t" << 0 << "\t" << 0 << endl;
      std::vector<int> v[1];
      for(int i=1;i<lagSize;i++) {
        int lag0 = 0;
        int lag1 = i;
        cout << i << "\t" << lag0 << "\t" << lag1 << endl;
        v[0].push_back(lag1-lag0);
      }
      for(int k=0;k<1;k++) {
        sort(v[k].begin(), v[k].end());
        auto it = std::unique(v[k].begin(), v[k].end());
        if(!(it == v[k].end())) {cout << "CWB_Plugin_UniqueLags.C : Error - Duplicate Lags" << endl; gSystem->Exit(1);} 
      }
    }
    if(cfg->nIFO==3) {
      cout << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
      std::vector<int> v[3];
      for(int i=1;i<lagSize;i++) {
        int lag0 = 0;
        int lag1 = i;
        int lag2 = (i*2)%lagSize;
        cout << i << "\t" << lag0 << "\t" << lag1 << "\t" << lag2 << endl;
        v[0].push_back(lag1-lag0);
        v[1].push_back(lag2-lag0); v[2].push_back(lag2-lag1);
      }
      for(int k=0;k<3;k++) {
        sort(v[k].begin(), v[k].end());
        auto it = std::unique(v[k].begin(), v[k].end());
        if(!(it == v[k].end())) {cout << "CWB_Plugin_UniqueLags.C : Error - Duplicate Lags" << endl; gSystem->Exit(1);} 
      }
    }
    if(cfg->nIFO==4) {
      cout << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
      std::vector<int> v[6];
      for(int i=1;i<lagSize;i++) {
        int lag0 = 0;
        int lag1 = i;
        int lag2 = (i*2)%lagSize;
        int lag3 = (i*3)%lagSize;
        cout << i << "\t" << lag0 << "\t" << lag1 << "\t" << lag2 << "\t" << lag3 << endl;
        v[0].push_back(lag1-lag0);
        v[1].push_back(lag2-lag0); v[2].push_back(lag2-lag1);
        v[3].push_back(lag3-lag0); v[4].push_back(lag3-lag1); v[5].push_back(lag3-lag2);
      }
      for(int k=0;k<6;k++) {
        sort(v[k].begin(), v[k].end());
        auto it = std::unique(v[k].begin(), v[k].end());
        if(!(it == v[k].end())) {cout << "CWB_Plugin_UniqueLags.C : Error - Duplicate Lags" << endl; gSystem->Exit(1);} 
      }
    }

    // restore cout 
    std::cout.rdbuf(old);

    // set lagBuffer
    lagBuffer->Set(buffer.str().size(),buffer.str().c_str());
    // set end of string
    lagBuffer->Set(lagBuffer->GetSize()+1);
    (*lagBuffer)[lagBuffer->GetSize()-1]=0;
    gCWB2G->lagMode[0]='s';gCWB2G->lagMode[1]=0;      // change lagMode : setTimeShifts read lags from string
    if(!gULAG) cfg->lagSize=lagSize;
    cfg->lagMax=lagMax;
    //cout << gCWB2G->lagBuffer.GetArray() << endl; 
  }

  return;
}

void ReadUserOptions(TString options) {

  if(options.CompareTo("")!=0) {
    cout << options << endl;
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet

      TObjArray* token = TString(options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("ulags_enabled=")) {
          TString ulags_enabled=stok;
          ulags_enabled.Remove(0,ulags_enabled.Last('=')+1);
          if(ulags_enabled=="true")  gOPT.enabled=true;
          if(ulags_enabled=="false") gOPT.enabled=false;
        }

        if(stok.Contains("ulags_lagSize=")) {
          TString ulags_lagSize=stok;
          ulags_lagSize.Remove(0,ulags_lagSize.Last('=')+1);
          if(ulags_lagSize.IsDigit()) gOPT.lagSize=ulags_lagSize.Atoi();
          if((gOPT.lagSize%2)==0) {
            cout << "ReadUserOptions Error : ulags_lagSize = " << gOPT.lagSize << " must be odd" << endl;exit(1);
          }
        }
      }
    }
  }
}

void PrintUserOptions(CWB::config* cfg) {

    cout << "-----------------------------------------"     << endl;
    cout << "UniqueLags config options                "     << endl;
    cout << "-----------------------------------------"     << endl << endl;
    if(gOPT.enabled==false) cout << "ULAGS_ENABLED = false" << endl;
    if(gOPT.enabled==true)  cout << "ULAGS_ENABLED = true " << endl;
    cout << "ULAGS_LAG_SIZE = " << gOPT.lagSize             << endl;
    cout << endl;
}

void ResetUserOptions(network* net, CWB::config* cfg) {
    gOPT.enabled=false;
    gOPT.lagSize=0;
}

