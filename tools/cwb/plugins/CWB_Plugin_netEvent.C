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
#include "gwavearray.hh"
#include <vector>

TFile* GetRootFile(network* NET);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Write user parameters to the output root file

  cout << endl;
  cout << "-----> CWB_Plugin_netEvent.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  float UserParms[NIFO_MAX];				// UserParms

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->outPlugin=true;  				// disable built-in output root file
  }

  if(type==CWB_PLUGIN_ILIKELIHOOD) {

    // search output root file in the system list
    TFile* froot = GetRootFile(NET);

    netevent* EVT;
    int nIFO = NET->ifoListSize();                      // number of detectors
    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree==NULL) {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("UserParms",UserParms,TString::Format("UserParms[%i]/F",cfg->nIFO));
    }
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_netEvent.C -> "
           << "CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_netEvent.C -> "
         << " gIFACTOR : " << gIFACTOR << endl;

    // import slagShift
    float* gSLAGSHIFT=NULL; IMPORT(float*,gSLAGSHIFT)

    int nIFO = NET->ifoListSize();                      // number of detectors
    int K = NET->nLag;                                  // number of time lag
    netevent* EVT;
    wavearray<double> id;
    //double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];
    double factor = cfg->factors[gIFACTOR];
    int rate = 0;                                       // select all resolutions

    // search output root file in the system list
    TFile* froot = GetRootFile(NET);

    TString outDump = froot->GetName();
    outDump.ReplaceAll(".root.tmp",".txt");

    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree!=NULL) {
      EVT = new netevent(net_tree,nIFO);
      net_tree->SetBranchAddress("UserParms",UserParms);
    } else {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("UserParms",UserParms,TString::Format("UserParms[%i]/F",cfg->nIFO));
    }
    EVT->setSLags(gSLAGSHIFT);        			// set slags into netevent

    for(int k=0; k<K; k++) {  				// loop over the lags

      id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);

      for(int j=0; j<(int)id.size(); j++) {  		// loop over cluster index

        int ID = size_t(id.data[j]+0.5);

        if(NET->getwc(k)->sCuts[ID-1]!=-1) continue;    // skip rejected/processed clusters

        double ofactor=0;
        if(cfg->simulation==4)      ofactor=-factor;
        else if(cfg->simulation==3) ofactor=-gIFACTOR;
        else                        ofactor=factor;

        EVT->output2G(NULL,NET,ID,k,ofactor);		// get reconstructed parameters

        netcluster* pwc = NET->getwc(k);

        UserParms[0]=111;
        UserParms[1]=222;

        std::vector<int> sCuts = NET->getwc(k)->sCuts;  // save cCuts
        // set sCuts=1 to the events which must be not copied with cps to pwc
        for(int i=0; i<(int)sCuts.size(); i++) if(i!=ID-1) NET->getwc(k)->sCuts[i]=1;

        // ID can not be used to get the event, to get event use ID=1 (ex: for watplot)
        NET->getwc(k)->sCuts = sCuts;                   // restore cCuts

        if(cfg->dump) EVT->dopen(outDump.Data(),const_cast<char*>("a"),false);
        EVT->output2G(net_tree,NET,ID,k,ofactor);       // get reconstructed parameters
        if(cfg->dump) {                                 
          // add UserParms to dump file
          fprintf(EVT->fP,"UserParms:      ");
          for(int i=0; i<nIFO; i++) fprintf(EVT->fP,"%f ",UserParms[i]);
          fprintf(EVT->fP,"\n");
        }
        if(cfg->dump) EVT->dclose();
        if(!cfg->cedDump) NET->getwc(k)->sCuts[ID-1]=1; // mark as processed
      }
    }

    jfile->cd();
    if(EVT) delete EVT;
  }
  return;
}

TFile* GetRootFile(network* NET) {

  // search output root file in the system list
  TFile* froot = NULL;                         
  TList *files = (TList*)gROOT->GetListOfFiles();
  TString outDump="";
  netevent* EVT;
  int nIFO = NET->ifoListSize();			// number of detectors
  if (files) {                                   
    TIter next(files);                           
    TSystemFile *file;                           
    TString fname;                               
    bool check=false;                            
    while ((file=(TSystemFile*)next())) {        
       fname = file->GetName();                  
       // set output root file as the current file
       if(fname.Contains("wave_")) {
         froot=(TFile*)file;froot->cd();
         outDump=fname;
         outDump.ReplaceAll(".root.tmp",".txt");
         //cout << "output file name : " << fname << endl;
       }
    }                                                               
    if(!froot) {                                                    
      cout << "CWB_Plugin_netEvent.C : Error - output root file not found" << endl;
      gSystem->Exit(1);                                                                             
    }                                                                                      
  } else {                                                                                 
    cout << "CWB_Plugin_netEvent.C : Error - output root file not found" << endl;  
    gSystem->Exit(1);                                                                               
  }                                                                                        

  return froot;
}

