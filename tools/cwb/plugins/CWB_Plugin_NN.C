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

//#define WATPLOT						// enable event monster plots 
#define MDC_OTF						// enable MDC On The Fly
#define NOISE_OTF					// enable NOISE On The Fly

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!NOISE_MDC_SIMULATION
// Plugin to include netcluster event stricture into the output wave root file (used in NN)

  cout << endl;
  cout << "-----> CWB_Plugin_NN.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
#ifdef NOISE_OTF
    cfg->dataPlugin=true; 				// disable read data from frames
#endif
#ifdef MDC_OTF
    cfg->mdcPlugin=true;  				// disable read mdc from frames
#endif
    cfg->outPlugin=true;  				// disable built-in output root file
  }


#ifdef NOISE_OTF
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
    TB.getSimNoise(*x, fName, seed, NET->nRun);
    x->resize(size);                           
    x->start(start);                           
  }                                            
#endif

#ifdef MDC_OTF
  if(type==CWB_PLUGIN_MDC) {  

    char cmd[128];
    sprintf(cmd,"network* net = (network*)%p;",NET);
    gROOT->ProcessLine(cmd);

    CWB::mdc MDC(NET);

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
    // get mdc data
    // ---------------------------------

    MDC.Get(*x,ifo);

    // ---------------------------------
    // set mdc list in the network class 
    // ---------------------------------

    if(ifo.CompareTo(NET->ifoName[0])==0) {
      NET->mdcList.clear();
      NET->mdcType.clear();
      NET->mdcTime.clear();
      NET->mdcList=MDC.mdcList;
      NET->mdcType=MDC.mdcType;
      NET->mdcTime=MDC.mdcTime;
    }

    cout.precision(14);
    for(int k=0;k<(int)NET->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)NET->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)NET->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;
  }
#endif

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_NN.C -> CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

#ifdef WATPLOT
    watplot WTS(const_cast<char*>("wts"));
#endif

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_NN.C -> " << " gIFACTOR : " << gIFACTOR << endl;

    // import slagShift
    float* gSLAGSHIFT=NULL; IMPORT(float*,gSLAGSHIFT)

    int nIFO = NET->ifoListSize();			// number of detectors
    int K = NET->nLag;  				// number of time lag          
    netevent* EVT;
    wavearray<double> id;
    double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];                 
    int rate = 0;					// select all resolutions

    // search output root file in the system list
    TFile* froot = NULL;                         
    TList *files = (TList*)gROOT->GetListOfFiles();
    TString outDump="";
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
           cout << "output file name : " << fname << endl;
         }
      }                                                               
      if(!froot) {                                                    
        cout << "CWB_Plugin_NN.C : Error - output root file not found" << endl;
        gSystem->Exit(1);                                                                             
      }                                                                                      
    } else {                                                                                 
      cout << "CWB_Plugin_NN.C : Error - output root file not found" << endl;  
      gSystem->Exit(1);                                                                               
    }                                                                                        

    netcluster* pwc = new netcluster();  
    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree!=NULL) {
      EVT = new netevent(net_tree,nIFO);
      net_tree->SetBranchAddress("netcluster",&pwc); 
    } else {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      // add new netcluster leaf
      net_tree->Branch("netcluster","netcluster",&pwc,32000,0); 
    }
    EVT->setSLags(gSLAGSHIFT);                          // set slags into netevent

    for(int k=0; k<K; k++) {  				// loop over the lags

      id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);

      for(int j=0; j<(int)id.size(); j++) {  		// loop over cluster index

        int ID = size_t(id.data[j]+0.5);

        if(NET->getwc(k)->sCuts[ID-1]!=-1) continue;    // skip rejected/processed clusters

        std::vector<int> sCuts = NET->getwc(k)->sCuts; 	// save cCuts
        // set sCuts=1 to the events which must be not copied with cps to pwc
        for(int i=0; i<(int)sCuts.size(); i++) if(i!=ID-1) NET->getwc(k)->sCuts[i]=1;  

        // after cpf, pwc contains only one event
        // ID can not be used to get the event, to get event use ID=1 (ex: for watplot)
        pwc->cpf(*(NET->getwc(k)));			// note: likelihood2G delete tdAmp 
        NET->getwc(k)->sCuts = sCuts; 			// restore cCuts

        double ofactor=0;
        if(cfg->simulation==4)      ofactor=-factor;
        else if(cfg->simulation==3) ofactor=-gIFACTOR;
        else                        ofactor=factor;
        if(cfg->dump) EVT->dopen(outDump.Data(),const_cast<char*>("a"),false);
        EVT->output2G(net_tree,NET,ID,k,ofactor);	// get reconstructed parameters
        if(cfg->dump) EVT->dclose();
        if(!cfg->cedDump) NET->getwc(k)->sCuts[ID-1]=1; // mark as processed

#ifdef WATPLOT						// monster event display
        WTS.canvas->cd();
        char fname2[1024];
        sprintf(fname2, "l_tfmap_scalogram_%d.png",ID);
        cout << "write " << fname2 << endl;
        WTS.plot(pwc, 1, nIFO, 'L', 0, const_cast<char*>("COLZ")); // use ID=1
        WTS.canvas->Print(fname2);
        WTS.clear();
#endif
      }
    }

    jfile->cd();

    if(pwc) delete pwc;
    if(EVT) delete EVT;
  }

  return;
}

