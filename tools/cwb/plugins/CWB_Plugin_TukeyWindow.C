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
// HOW TO SETUP PLUGIN WF IN CWB USER CONFIGURATION (EXAMPLE)
// ---------------------------------------------------------------------------------

/*
  // Example for the L1 Glitch in GW170817

  TString opttw = "";                       // NOTE : add space at the end of each line

  opttw += "tw_gps_time=1187008881.38 ";    // TW central GPS time
  opttw += "tw_width=0.34 ";                // TW width
  opttw += "tw_range=0.14 ";                // TW range
  opttw += "tw_ced_spectrogram_zmax=1e-3 "; // TW set zmax in CED spectrograms
  opttw += "tw_dump=false ";                // TW dump ascii file

  strcpy(parPlugin,opttw.Data());           // set TW plugin parameters
  strcpy(comment,"TW");

*/

// ---------------------------------------------------------------------------------
// TUKEY WINDOW DEFINES
//
//             GPS_TIME
//                |
//                ^
// ----->       WIDTH       <------
//       \                 /
//        \               /
// -------->    RANGE    <--------
//
// ---------------------------------------------------------------------------------

#define TW_GPS_TIME    		0.0		
#define TW_WIDTH    		0.0			 
#define TW_RANGE		0.0
#define TW_CED_SPECTROGRAM_ZMAX	0.0
#define TW_DUMP			false	// true -> window is dumped to ascii file

// ---------------------------------------------------------------------------------
// USER TW PLUGIN OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  double gps_time; 
  double width;   
  double range;   
  double ced_spectrogram_zmax;   
  bool   dump;     
};

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions();

double TukeyWin(double x, double W, double R);
void DumpData(wavearray<double>* x , TString ofname);

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions gOPT;                           // global User Options



void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Plugin used to remove the L1 GW170817 glitch

  char ofname[256];

  if(type==CWB_PLUGIN_CONFIG) {

    ResetUserOptions();                                 // set default config options
    ReadUserOptions(cfg->parPlugin);                    // user config options : read from parPlugin
    PrintUserOptions();      				// print config options

    if(gOPT.gps_time==0) {
      cout << "CWB_Plugin_TukeyWindow Error : tw_gps_time not defined !!!" << endl;
      exit(1);
    }

    if(gOPT.width==0) {
      cout << "CWB_Plugin_TukeyWindow Error : tw_width not defined !!!" << endl;
      exit(1);
    }

    if(gOPT.range==0) {
      cout << "CWB_Plugin_TukeyWindow Error : tw_range not defined !!!" << endl;
      exit(1);
    }
  }

  if(type==CWB_PLUGIN_STRAIN) {
    double start = x->start();
    double stop  = x->start()+x->size()/x->rate();
    if(!(gOPT.gps_time>start && gOPT.gps_time<stop)) return;
    if(ifo=="L1") {
      cout << endl << "CWB_Plugin_TukeyWindow.C : Apply Tukey Window to GW170817 L1 Glitch !!!" << endl << endl;

      if(gOPT.dump) {
        sprintf(ofname,"%s_TukeyWindow.dat",ifo.Data());
        DumpData(x, ofname);
      }

      wavearray<double> TW = *x;
      double dt=1./x->rate();
      for(int i=0;i<x->size();i++) {
        double t = start+i*dt-(gOPT.gps_time-gOPT.width/2);
        double w = 1-TukeyWin(t,gOPT.width,gOPT.range);
        x->data[i]*=w;
        TW[i]=w;
      }

      if(gOPT.dump) {
        sprintf(ofname,"%s_TukeyWindow.dat",ifo.Data());
        DumpData(&TW, ofname);

        sprintf(ofname,"%s_TukeyWindow.dat",ifo.Data());
        DumpData(x, ofname);
      }
    }      
  }

  if(type==CWB_PLUGIN_CED) {
    // get object CED pointer 
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
    CWB::ced* ced = gCWB2G->GetCED();
    char ced_options[64];
    sprintf(ced_options,"ced_spectrogram_zmax=%f",gOPT.ced_spectrogram_zmax);
    if(ced!=NULL) ced->SetOptions(ced_options);
  }
}

double TukeyWin(double x, double W, double R) {

  double y=0;

  if(x<0 || x>W) return y;

  double Pi = TMath::Pi();

  if(x>=0       && x<R/2  ) y = 0.5*(1+cos(2*Pi/R*(x-R/2)));
  if(x>= R/2    && x<W-R/2) y = 1;
  if(x>=(W-R/2) && x<=W   ) y = 0.5*(1+cos(2*Pi/R*(x-W+R/2)));

  return y;
}

void DumpData(wavearray<double>* x , TString ofname) {

  ofstream out;
  out.open(ofname.Data(),ios::out);
  if (!out.good()) {cout << "CWB_Plugin_TukeyWindow.C : Error Opening Output File : " << ofname << endl;exit(1);}
  cout << "Create Output File : " << ofname << endl;
  out.precision(19);

  // write data
  int size = x->size();
  double dt=1./x->rate();
  for (int i=0; i<size; i++) {
    double start = gOPT.gps_time-x->start()-1;
    double stop  = start+2;
    double time = i*dt;
    if(time>start && time<stop) out << time << " " << x->data[i] << endl;
  }

  out.close();
}

void ReadUserOptions(TString options) {

  // get plugin options 

  if(TString(options)!="") {

    //cout << "SGW options : " << options << endl;
    TObjArray* token = TString(options).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++) {

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();

      if(stok.Contains("tw_dump=")) {
        TString tw_dump=stok;
        tw_dump.Remove(0,tw_dump.Last('=')+1);
        if(tw_dump=="root")  gOPT.dump=true;
      }

      if(stok.Contains("tw_gps_time=")) {
        TString tw_gps_time=stok;
        tw_gps_time.Remove(0,tw_gps_time.Last('=')+1);
        if(tw_gps_time.IsFloat()) gOPT.gps_time=tw_gps_time.Atof();
      }

      if(stok.Contains("tw_width=")) {
        TString tw_width=stok;
        tw_width.Remove(0,tw_width.Last('=')+1);
        if(tw_width.IsFloat()) gOPT.width=tw_width.Atof();
      }

      if(stok.Contains("tw_range=")) {
        TString tw_range=stok;
        tw_range.Remove(0,tw_range.Last('=')+1);
        if(tw_range.IsFloat()) gOPT.range=tw_range.Atof();
      }

      if(stok.Contains("tw_ced_spectrogram_zmax=")) {
        TString tw_ced_spectrogram_zmax=stok;
        tw_ced_spectrogram_zmax.Remove(0,tw_ced_spectrogram_zmax.Last('=')+1);
        if(tw_ced_spectrogram_zmax.IsFloat()) gOPT.ced_spectrogram_zmax=tw_ced_spectrogram_zmax.Atof();
      }
    }
  }
}

void ResetUserOptions() {

    gOPT.gps_time  = TW_GPS_TIME;
    gOPT.width     = TW_WIDTH;
    gOPT.range     = TW_RANGE;
    gOPT.range     = TW_CED_SPECTROGRAM_ZMAX;
    gOPT.dump      = TW_DUMP;
} 


void PrintUserOptions() {

    cout.precision(14);
    cout << "-----------------------------------------"     << endl;
    cout << "TW config options                        "     << endl;
    cout << "-----------------------------------------"     << endl << endl;
    cout << "TW_GPS_TIME         	" << gOPT.gps_time         	<< endl;
    cout << "TW_WIDTH            	" << gOPT.width            	<< endl;
    cout << "TW_RANGE            	" << gOPT.range            	<< endl;
    cout << "TW_CED_SPECTROGRAM_ZMAX    " << gOPT.ced_spectrogram_zmax  << endl;
    cout << "TW_DUMP             	" << gOPT.dump             	<< endl;
    cout << endl;
}

