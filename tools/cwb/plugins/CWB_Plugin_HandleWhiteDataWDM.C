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
#include "watplot.hh"


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!DATA_HANDLING
// Plugin to handle whitened input data after its WDM transform  (only for 2G analysis !!!)

//  #define PLOT_WDM_TF           // Plot WDM Scalogram
//  #define DATA_COMPARISON       // Compare data before/after WDM data handling
  #define DATA_HANDLING         // example of data handling
  #define COPY_HANDLED_DATA     // copy handled data to into the original x array

  cout << endl;
  cout << "-----> CWB_Plugin_HandleWhiteDataWDM.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_WHITE) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_HandleWhiteDataWDM.C : Error - analysis type must be 2G !!!" << endl;
      gSystem->Exit(1);
    }

    wavearray<double> y = *x;  			// copy original data to y (contains the final handled data)
    WSeries<double> w;				// temporary array for data manipulation
    WDM<double>* pwdm = NULL;

    for(int level=cfg->l_high; level>=cfg->l_low; level--) {    // loop over levels

      int layers = level>0 ? 1<<level : 0;	// get layers
      cout << "level : " << level << " layers : " << layers << endl;

      pwdm = net->getwdm(layers+1);		// get pointer to wdm transform
      if(pwdm==NULL) {
        cout << "CWB_Plugin_HandleWhiteDataWDM.C : Error - WDM not defined !!!" << endl;
        gSystem->Exit(1);
      }
 
      w.Forward(y,*pwdm);			// apply wdm transformation 

#ifdef DATA_HANDLING
      // example of data handling
      double R = x->rate();

      double dF = R/(2*layers);     		// frequency resolution
      double dT = layers/R;          		// time resolution

      cout << "WDM Decomposition Level : " << w.getLevel()
           << " dF : " << dF << " Hz " << " dT : " << 1000*dT << " ms " <<endl;

      // Compute Energy
      int M = w.getLevel();
      int L = pwdm->maxLayer();
      cout << "WDM Decomposition Level : " << M << " Layers " << L << endl;

      int mF = int(w.size()/pwdm->nSTS);
      int nTC = w.size()/(M+1)/mF;		// # of Time Coefficients
      // Extract map00 & map90 : size = (M+1)*nTC
      double* map00 = pwdm->pWWS;
      double* map90 = map00 + (mF-1)*(M+1)*nTC;

      wavearray<double> xx;
      wavearray<double> XX;
      double E00=0;
      double E90=0;
      // this example show the corrispondence between map00,map90 and xx,XX
      for(int j=0;j<M+1;j++) {      		// loop over the layers
        w.getLayer(xx,j);           		// extract phase 00 @ layer j
        w.getLayer(XX,-j);          		// extract phase 90 @ layer j
        // for phase 00 -> (xx[i] = map00[i*(M+1)+j])
        // for phase 90 -> (XX[i] = map90[i*(M+1)+j])
        for(int i=0;i<xx.size();i++) {
          //cout << j << " P00 " << xx[i] << " " << map00[i*(M+1)+j]
          //<< " P90 " << XX[i] << " " << map90[i*(M+1)+j] << endl;
        }
        for(int i=0;i<xx.size();i++) E00+=xx[i]*xx[i];
        for(int i=0;i<XX.size();i++) E90+=XX[i]*XX[i];
      }
      cout << "E00 = " << E00 << endl;
      cout << "E90 = " << E90 << endl;
#endif

#ifdef PLOT_WDM_TF
      // Plot WDM Scalogram
      watplot WTS(const_cast<char*>("wdm"));
      //scalogram maps
      double start = w.start();
      double stop  = w.start()+w.size()/w.rate();
      double flow  = 64;
      double fhigh = 2048;
      WTS.plot(&w, 2, start, stop,const_cast<char*>("COLZ"));
      WTS.hist2D->GetYaxis()->SetRangeUser(flow, fhigh);
      // dump spectrum
      char fname[1024];
      sprintf(fname,"%s_wdm_scalogram_%d_lev%d.root", ifo.Data(), int(w.start()),w.getLevel());
      cout << endl << "Dump WDM Scalogram : " << fname << endl << endl;
      WTS.canvas->Print(fname);
#endif

      w.Inverse();				// inverse transform
      y = w;					// copy manipulated data to y
    }

#ifdef DATA_COMPARISON
    // Compare data before/after WDM data handling
    wavearray<double> z = y;

    watplot plot(const_cast<char*>("plot"),200,20,800,500);
    char gtitle[256];
    sprintf(gtitle,"WDM : Original(black) - After WDM Forward/Inverse Tranform(Red)");
    plot.gtitle(gtitle,"time(sec)","amplitude");
    double start = z.start()+cfg->segEdge;
    double stop  = z.start()+z.size()/z.rate()-cfg->segEdge;
    plot.goptions("alp", 1, start, stop);
    // draw signal
    z-=*x;       				// difference
    //*x >> plot;
    z >> plot;

    int nEdge = cfg->segEdge*y.rate();
    double ez=0;
    for(int i=nEdge;i<(int)y.size()-nEdge;i++) ez+=y[i]*y[i];
    double ex=0;
    for(int i=nEdge;i<(int)x->size()-nEdge;i++) ex+=x->data[i]*x->data[i];
    cout << "Energy of Difference / Energy of Input Data " << ez/ex << endl;

    //TString gfile=ifo+"_WDM_data_comparison.png"; 	// save plot the difference to png file
    TString gfile=ifo+"_WDM_data_comparison.root"; 	// save plot the difference to root file
    plot >> gfile;
#endif

#ifdef COPY_HANDLED_DATA
    // copy handled data to into the original x array
    for(int i=0;i<(int)x->size();i++) x->data[i]=y[i];
#endif
  }

  return;
}
