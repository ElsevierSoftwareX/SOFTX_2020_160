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
#include <vector>


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!NOISE_MDC_SIMULATION
// Plugin to generate simulated gaussian noise and injected 'on the fly' splitted MDC
// The standard procedure inject signals only within the run
// This plugin allows to test injected signals splitted into two consecutive runs
// if simulation==0  the MDC are multiplied by factors[0] and added to the noise 
// WARNING!!! this plugin can not be used in snr mode (simulation=2)
//            because the events are splitted into two subevents 
//            and the subevents are differently normalized by the snr mode procedure 

  cout << endl;
  cout << "-----> CWB_Plugin_SplittedMDC.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->dataPlugin=true; // disable read data from frames
    cfg->mdcPlugin=true;  // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_DATA) {  

    CWB::Toolbox TB;

    int seed;
    if(ifo.CompareTo("L1")==0) seed=1000;
    if(ifo.CompareTo("H1")==0) seed=2000;
    if(ifo.CompareTo("V1")==0) seed=3000;
    if(ifo.CompareTo("J1")==0) seed=4000;
    if(ifo.CompareTo("A2")==0) seed=5000;
    if(ifo.CompareTo("Y2")==0) seed=6000;
    if(ifo.CompareTo("Y3")==0) seed=7000;

    TString fName;
    if(ifo.CompareTo("L1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("V1")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";
    if(ifo.CompareTo("J1")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";
    if(ifo.CompareTo("A2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("Y2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("Y3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";

    int size=x->size();
    double start=x->start();
    TB.getSimNoise(*x, fName, seed, net->nRun);
    x->resize(size);
    x->start(start);
  }

  // if simulation==0  the MDC are multiplied by factors[0] and added to the noise 
  if((type==CWB_PLUGIN_MDC)||((cfg->simulation==0)&&(type==CWB_PLUGIN_DATA))) {  

    if(cfg->simulation==2) {
      cout << "CWB_Plugin_SplittedMDC.C : Error - "
           << "Injections on the Edge are not allowed in snr mode (simulation=2)" << endl;
      gSystem->Exit(1);
    }

    char cmd[128];
    cout.precision(14);
    sprintf(cmd,"network* net = (network*)%p;",net);
    gROOT->ProcessLine(cmd);
    sprintf(cmd,"CWB::config* cfg = (CWB::config*)%p;",cfg);
    gROOT->ProcessLine(cmd);

    CWB::mdc MDC(net);

    // ---------------------------------
    // read plugin config 
    // ---------------------------------

    cfg->configPlugin.Exec();

    // ---------------------------------
    // set list of mdc waveforms
    // ---------------------------------
   
    IMPORT(CWB::mdc,MDC) 
    MDC.Print();

    // ---------------------------------
    // clear mdc list in the network class 
    // ---------------------------------

    if(cfg->simulation) {
      if(ifo.CompareTo(net->ifoName[0])==0) {
        net->mdcList.clear();
        net->mdcType.clear();
        net->mdcTime.clear();
      }
    }

    // --------------------------------------------------------------------------------------
    // get mdc data
    // --------------------------------------------------------------------------------------
    //
    // xxxx             segEdge            xxxx
    //                
    //     <--------segLen=2*tshift------->  
    //
    // xxxx---------------||---------------xxxx                                left shifted
    //                xxxx---------------||---------------xxxx                 current segment
    //                               xxxx---------------||---------------xxxx  right shifted
    //

    wavearray<double> y(x->size()); 
    y.rate(x->rate()); 
    y.start(x->start());
    double tshift = (x->size()/x->rate()-2*cfg->segEdge)/2.;
    int nedge = int(cfg->segEdge*x->rate());  

    // get buffer half left shifted 
    y.start(x->start()-tshift);						
    MDC.Get(y,ifo);

    if(cfg->simulation==0) {
      for(int i=0;i<x->size()/2;i++) x->data[i]+=cfg->factors[0]*y[i+x->size()/2-nedge];
    } else {
      for(int i=0;i<x->size()/2;i++) x->data[i]=y[i+x->size()/2-nedge];
      if(ifo.CompareTo(net->ifoName[0])==0) {
        for(int k=0;k<(int)MDC.mdcList.size();k++) { 
          net->mdcList.push_back(MDC.mdcList[k]); 
          net->mdcType.push_back(MDC.mdcType[k]); 
          net->mdcTime.push_back(MDC.mdcTime[k]); 
        }
      }
    }

    // get buffer half right shifted 
    y.start(x->start()+tshift);						
    MDC.Get(y,ifo);

    if(cfg->simulation==0) {
      for(int i=0;i<x->size()/2;i++) x->data[i+x->size()/2]+=cfg->factors[0]*y[i+nedge];
    } else {
      for(int i=0;i<x->size()/2;i++) x->data[i+x->size()/2]=y[i+nedge];
      if(ifo.CompareTo(net->ifoName[0])==0) {
        for(int k=0;k<(int)MDC.mdcList.size();k++) { 
          // if event is not duplicated it is added to the network mdc list
          bool duplicated=false;
          for(int m=0;m<(int)net->mdcList.size();m++) if(MDC.mdcTime[k]==net->mdcTime[m]) duplicated=true;
          if(!duplicated) { 
            net->mdcList.push_back(MDC.mdcList[k]); 
            net->mdcType.push_back(MDC.mdcType[k]); 
            net->mdcTime.push_back(MDC.mdcTime[k]); 
          }
        }
      }
    }

    if(cfg->simulation) {
      for(int k=0;k<(int)net->mdcList.size();k++) {
        cout << k << " mdcList " << net->mdcList[k] << endl;
        cout << k << " mdcTime " << net->mdcTime[k] << endl;
        cout << k << " mdcType " << net->mdcType[k] << endl;
      }
    }
  }

  return;
}
