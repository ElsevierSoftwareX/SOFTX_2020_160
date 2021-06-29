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
#include "gwavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"

#define L1_AMP_CAL_ERR_SIGMA	0.2	// 20% amplitude uncertainty 
#define H1_AMP_CAL_ERR_SIGMA	0.2	// 20% amplitude uncertainty 
#define V1_AMP_CAL_ERR_SIGMA	0.1	// 10% amplitude uncertainty 

//#define SAVE_PLOT			// plot data in ROOT format before/after mis-calibration 

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!MISCALIBRATION
// Plugin to mis-calibrate in amplitude the MDC

// The jittering procedure depends on the model describing the calibration uncertainties distribution. 
// For each detector, jittering values are extracted, via MonteCarlo, from a gaussian distribution, 
// centered on 1 and with standard deviation equal to the interferometer nominal uncertainty 
// (i.e., 20% for LIGO observatories, 10% for Virgo). 
// Such values are then used to rescale the MDC waveforms amplitude before the injection in the detectors data. 

// The waveforms are injected in the detectors data after having been jittered, with different 
// values for the three detectors. Since cWB performs a coherent search, a fraction of events 
// is therefore rejected by the algorithm. How large such fraction is, depends on how extreme 
// the considered jittering models are. It is expected to get larger for increasing possible jittering values. 
// Given that part of the jittered injected signals is rejected by the algorithm, the detection efficiency 
// curves are expected to get lower. Such variations in the detection efficiency induce changes 
// of the 90% confidence upper limit on the total event rate. i

  cout << endl;
  cout << "-----> CWB_Plugin_AmplitudeMisCal.C" << endl;
  cout << "ifo  " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_MDC) { // mdc
    cout << "APPLY AMPLITUDE MIS-CALIBRATION !!!" << endl;
#ifdef SAVE_PLOT
    gwavearray<double> gx(x);
    gx.Draw(GWAT_TIME);			// draw signal before cfactor (BLACK)
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
      double cfactor = 0;
      if(ifo=="L1") cfactor = gRandom->Gaus(1,L1_AMP_CAL_ERR_SIGMA);
      if(ifo=="H1") cfactor = gRandom->Gaus(1,H1_AMP_CAL_ERR_SIGMA);
      if(ifo=="V1") cfactor = gRandom->Gaus(1,V1_AMP_CAL_ERR_SIGMA);
      //cout << " starti " << starti << " stopi " << stopi << " cfactor : " << cfactor << endl;
      //cout << "Apply amplitude factor : " << cfactor << endl;
      for (int jj=starti; jj<stopi; jj++) {
        //if(x->data[jj]!=0) cout << jj << " " << cfactor << " " << x->data[jj] << endl;
        x->data[jj] *= cfactor;
      }
    }
#ifdef SAVE_PLOT
    gx.Draw(x,GWAT_TIME,"SAME",kRed);		// draw signal after cfactor (RED)
    TString gfile="CWB_Plugin_AmplitudeMisCal_Plot_MDC_"+ifo+".root";
    watplot* plot = gx.GetWATPLOT();		// get pointer to watplot object
    (*plot) >> gfile;
    cout << "CWB_Plugin_AmplitudeMisCal.C : created plot file name : " << gfile << endl;
#endif
  }

  return;
}
