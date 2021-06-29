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
#include "watplot.hh"
#include <vector>

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!NOISE_MDC_SIMULATION
// Plugin used in simulation mode to compute and save the mdc SNR (terminate after the read noise/mdc stage) 

  cout << endl;
  cout << "-----> CWB_Plugin_SNR.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  // set sim in simulation mode = 2
  if(type==CWB_PLUGIN_CONFIG) {	
    cfg->simulation = 2;
    cfg->nfactor = 1;
    cfg->factors[0] = 1;
  }

  if(type==CWB_PLUGIN_OREADDATA) {

    if(cfg->simulation!=2) { 
      cout << "CWB_Plugin_SNR.C - Error : this pluging is implemented only in simulation mode = 2" << endl;
      gSystem->Exit(1);
    }
 
    // create output root file to store snr infos
    char snrFile[1024]; 
    sprintf(snrFile,"%s/snr_%s_%d.root",cfg->output_dir,cfg->data_label,NET->nRun);
    TFile* froot = new TFile(snrFile, "RECREATE");
    if(froot==NULL) {
      cout << "CWB_Plugin_SNR.C - Error opening root file : " << snrFile << endl;
      gSystem->Exit(1);
    }
    froot->cd();

    // create mdc tree
    int nIFO = NET->ifoListSize();			// number of detectors
    injection mdc(nIFO);
    TTree* mdc_tree = mdc.setTree();
    // add iSNR branch to mdc_tree
    char ciSNR[16]; sprintf(ciSNR,"iSNR[%1d]/F",nIFO);
    float* iSNR = new float[nIFO];
    mdc_tree->Branch("iSNR", iSNR, ciSNR);

    NET->setVeto(cfg->iwindow);				// select mdc 
                                                        // in pipeline is done in the coherence stage

    int N = NET->mdc__ID.size();			// number of injectios

    std::vector<size_t> mdc__ID = NET->mdc__ID;		// save a copy mdc__ID

    NET->mdc__ID.resize(1);				// resize to 1 to save single event
    for (int k=0;k<N;k++) {				// loop over mdc
      int ID = mdc__ID[k];
      NET->mdc__ID[0] = ID;
      for(int i=0; i<nIFO; i++) {
        detector* pD = NET->getifo(i);
        iSNR[i] = pD->ISNR.data[ID];			// fill mdc with iSNR
        //cout << k << " SNR " << iSNR[i] << endl;
      }
      mdc.output(mdc_tree,NET,1);			// fill tree with mdc parameters
    }

    NET->mdc__ID = mdc__ID;				// restore mdc__ID

    mdc_tree->Write();					// save tree
    froot->Close();

    delete iSNR;

    gSystem->Exit(0); 					// terminate job

  }

  return;
}

