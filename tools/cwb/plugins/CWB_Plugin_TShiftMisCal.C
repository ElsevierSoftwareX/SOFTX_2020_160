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

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!MISCALIBRATION
// Plugin to mis-calibrate in time the injeted MDC

  WSeries<double> y;

  cout << endl;
  cout << "-----> CWB_Plugin_TShiftMisCal.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_MDC) { // mdc

    // TEST PHASE CALIBRATION ERRORS  //GV

    // upconvertion 16kHz -> 64kHz

    x->FFTW(1);
    x->resize(4*x->size());
    x->rate(4*x->rate());
    x->FFTW(-1);

    y.resize(x->size());y.rate(x->rate());y.start(x->start());y=0;  // temporary array
    for (int nmdc=0; nmdc<(int)net->mdcListSize(); nmdc++) {
      TString mdcstring(net->getmdcList(nmdc));
      TObjArray* token = mdcstring.Tokenize(' ');
      TObjString* iname = (TObjString*)token->At(11);
      TString wavename = iname->GetString();
      TObjString* itime = (TObjString*)token->At(10);
      TString wavetime = itime->GetString();
      double mdctime = wavetime.Atof();
      if (mdctime<(int)x->start()||mdctime>x->start()+x->size()/x->rate()) continue;
      //cout << mdcstring.Data() << endl;
      //printf(" Time : %s %f %f %f", wavetime.Data(), mdctime, x->start(), x->start()+x->size()/x->rate());
      //cout << "String: " << wavename.Data() << " time : " << wavetime.Data() << " " <<mdctime<< endl;
      int cshift = 2*int(gRandom->Integer(3)-1);

      int starti = (mdctime - x->start()-1.)*x->rate();
      int stopi  = (mdctime - x->start()+1.)*x->rate();
      if(cshift>0) {
        if (starti<cshift) starti=cshift;
        if (stopi>(int)x->size()-cshift) stopi=x->size()-cshift;
      } else {
        if (starti<0) starti=0;
        if (stopi>(int)x->size()+cshift) stopi=x->size()+cshift;
      }

      cout << " Start: " << starti << " Stop : " << stopi << " cshift : " << cshift << endl;
      for (int jj=starti; jj<stopi; jj++) {
        y.data[jj] = x->data[jj-cshift];
      }
      //for (int jj=starti; jj<stopi; jj++) {
      //  if(x->data[jj]!=0) cout << jj << " " << cshift << " " << x->data[jj] << " " << y.data[jj] << endl;
      //}
    }
    (*x)=y;y.resize(0);

    // downconvertion 64kHz -> 16kHz
    x->FFTW(1);
    x->resize(x->size()/4);
    x->rate(x->rate()/4);
    x->FFTW(-1);
  }

  return;
}
