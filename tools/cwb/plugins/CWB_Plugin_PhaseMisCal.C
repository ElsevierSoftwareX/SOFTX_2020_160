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
#include "mdc.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "gwavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "TComplex.h"

#define L1_PHS_CAL_ERR    10.     // 10 degrees in phase uncertainty
#define H1_PHS_CAL_ERR    10.     // 10 degrees in phase uncertainty
#define V1_PHS_CAL_ERR    10.     // 10 degrees in phase uncertainty

//#define SAVE_PLOT                 // plot data in ROOT format before/after mis-calibration

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!MISCALIBRATION
// Plugin to mis-calibrate in phase the MDC

// This macro implements the the test of the cWB offline analysis 
// robustness versus the official calibration systematics in phase.

  wavearray<double> y;

  cout << endl;
  cout << "-----> CWB_Plugin_PhaseMisCal.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_MDC) { // mdc
    cout << "APPLY PHASE MIS-CALIBRATION !!!" << endl;
#ifdef SAVE_PLOT
    gwavearray<double> gx(x);
    gx.Draw(GWAT_TIME);                 // draw signal before cfactor (BLACK)
#endif
    for (int nmdc=0; nmdc<(int)net->mdcListSize(); nmdc++) {
      TString mdcstring(net->getmdcList(nmdc));
      TObjArray* token = mdcstring.Tokenize(' ');
      TObjString* iname = (TObjString*)token->At(11);
      TString wavename = iname->GetString();
      TObjString* itime = (TObjString*)token->At(10);
      TString wavetime = itime->GetString();
      double mdctime = wavetime.Atof();
      if (mdctime<x->start()||mdctime>x->start()+x->size()/x->rate()) continue;
      //cout << mdcstring.Data() << endl;
      //printf(" Time : %s %f %f %f", wavetime.Data(), mdctime, x->start(), x->start()+x->size()/x->rate());
      //cout << "String: " << wavename.Data() << " time : " << wavetime.Data() << " " <<mdctime;
      int starti = (mdctime - x->start()-cfg->iwindow/2.)*x->rate();
      int stopi  = (mdctime - x->start()+cfg->iwindow/2.)*x->rate();
      if (starti<0) starti=0;
      if (stopi>(int)x->size()) stopi=x->size();
 
      double pShift = 0;
      if(ifo=="L1") pShift = L1_PHS_CAL_ERR;
      if(ifo=="H1") pShift = H1_PHS_CAL_ERR;
      if(ifo=="V1") pShift = V1_PHS_CAL_ERR;
      pShift *= (2*gRandom->Integer(2)-1.);	// +/- pShift

      //cout << " starti " << starti << " stopi " << stopi << " pShift : " << pShift << endl;
      int size = TMath::Nint(cfg->iwindow+0.5)*x->rate();
      wavearray<double> X(size);X=0;
      for (int jj=starti; jj<stopi; jj++) X[jj-starti] = x->data[jj];
      //cout << "Apply phase shift : " << pShift << " degrees" << endl;
      CWB::mdc::PhaseShift(X,pShift);
      for (int jj=starti; jj<stopi; jj++) x->data[jj] = X[jj-starti];
    }
#ifdef SAVE_PLOT
    gx.Draw(x,GWAT_TIME,"SAME",kRed);           // draw signal after cfactor (RED)
    TString gfile="CWB_Plugin_PhaseMisCal_Plot_MDC_"+ifo+".root";
    watplot* plot = gx.GetWATPLOT();            // get pointer to watplot object
    (*plot) >> gfile;
    cout << "CWB_Plugin_PhaseMisCal.C : created plot file name : " << gfile << endl;
#endif
  }

  return;
}
