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

//!DATA_CONDITIONING

// Implements the 2G frequency cuts in pixel selection stage

inline bool  CheckRange(int i, double df, double flow, double fhigh)
             {
		double iflow=i*df-df/2.;
		double ifhigh=(i+1)*df-df/2.;
                if((iflow>=flow)&&(iflow<fhigh))   return true;
                if((ifhigh>flow)&&(ifhigh<=fhigh)) return true;
                return false;
             }

// fcut structure used to define the frequency cuts
#define _FCUT_MAX		100

int FCUT_MAX=_FCUT_MAX;

struct fcut {
  char  ifo[32];
  float flow;
  float fhigh;
  char  levels[1024];
};


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  cout << endl;
  cout << "-----> CWB_Plugin_fCuts.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_XCOHERENCE) {

    // WSeries contains energy (it is stored in the 0 phase amplitudes)

    // frequency cuts declaration
    fcut FCUT[_FCUT_MAX];

    // export to CINT 
    char cmd[128];
    sprintf(cmd,"fcut* FCUT = (fcut*)%p;",FCUT);
    gROOT->ProcessLine(cmd);

    // read plugin config
    cfg->configPlugin.Exec();
    // import resolution level
    int nFCUT=-1; IMPORT(int,nFCUT)
    // transform level list { 2, 6, 8} to {,2,6,8,}
    for(int i=0;i<nFCUT;i++) {
      TString Levels = FCUT[i].levels;
      Levels = ","+Levels+",";
      Levels.ReplaceAll(" ","");
      strcpy(FCUT[i].levels,Levels.Data());
    }

    // import resolution level
    size_t gILEVEL=-1; IMPORT(size_t,gILEVEL)
    char slevel[256];sprintf(slevel,",%d,",gILEVEL); 

    int nIFO = net->ifoListSize();
    detector* pD[NIFO_MAX];                  // pointers to detectors
    for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);

    for(int n=0; n<nIFO; n++) {              // produce TF maps with max over the sky energy

      WSeries<double>* WS = pD[n]->getTFmap();
      ifo = net->ifoName[n];

      if(WS->pWavelet->m_WaveType!=WDMT) {
        cout << "CWB_Plugin_fCuts.C - Error : works only for WDMT wavelet types" << endl;
        gSystem->Exit(1);                                                                   
      }                                                                                     


      int layers = WS->maxLayer()+1;  	// numbers of frequency bins (first & last bins have df/2)
      int slices = WS->sizeZero();    	// number of time bins

      float df = WS->resolution();      // frequency bin resolution (hz)
      float dt = 1./(2*df);             // time bin resolution (sec)

      int rate = int(1./dt);                                                                           

      if(n==0) cout << "layers : " << layers << "\t slices : " << slices << "\t rate : " << rate
                    << "\t dt : " << dt << "\t df : " << df << endl;                            

      for(int i=0;i<slices;i++) {
        for(int j=0;j<layers;j++) {
          for(int k=0;k<nFCUT;k++) { 
            if(FCUT[k].ifo==ifo && TString(FCUT[k].levels).Contains(slevel)) {
              if(CheckRange(j,df,FCUT[k].flow,FCUT[k].fhigh)) {
                WS->putSample(0,i,j);		// set to 0 the energy
              }
            }
          }
        }                                                                                                     
      }                                                                                                       
    }
  }
  return;
}
