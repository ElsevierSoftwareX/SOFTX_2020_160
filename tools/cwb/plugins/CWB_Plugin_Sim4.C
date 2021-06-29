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

size_t SetFactors(wavearray<double> &wi, std::vector<double>* pT, std::vector<double>* pF, double dT);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  // Plugin to generate simulated MDC injected 'on the fly' using simulation=4

  cout << endl;
  cout << "-----> CWB_Plugin_Sim4.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  CWB::Toolbox TB;

  // import istage,jstage
  int gISTAGE=-1; IMPORT(int,gISTAGE) 
  int gJSTAGE=-1; IMPORT(int,gJSTAGE) 
  cout << "-----> CWB_Plugin_Sim4.C -> " << " gISTAGE : " << gISTAGE << " gJSTAGE : " << gJSTAGE << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->dataPlugin=false; // disable read data from frames
    cfg->mdcPlugin=true;  // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_MDC) {  

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR) 
    cout << "-----> CWB_Plugin_Sim4.C -> " << " gIFACTOR : " << gIFACTOR << endl;

    // export to CINT net,cfg
    char cmd[128];
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
   
    // import MDC form plugin config
    IMPORT(CWB::mdc,MDC) 
    MDC.Print();

    // import mdcHRSS form plugin config
    std::vector<double> mdcHRSS;
    IMPORT(std::vector<double>,mdcHRSS) 
    int M = mdcHRSS.size();
    for(int i=0;i<M;i++) cout << i << " mdcHRSS = " << mdcHRSS[i] << endl;

    // ---------------------------------
    // get mdc data
    // ---------------------------------

    MDC.Get(*x,ifo);

    // compute normalization factors
    int K = MDC.mdcList.size();
    int L = MDC.mdcType.size();
    std::vector<double> mdcFactor(K);
    double inj_hrss = MDC.GetInjHrss();
    for(int i=0;i<K;i++) {
      for(int j=0;j<L;j++) {
        if(TString(MDC.mdcList[i]).Contains(MDC.mdcType[j])) mdcFactor[i]=mdcHRSS[j]/inj_hrss;
      } 
    }
    for(int i=0;i<K;i++) cout << i << " mdcFactor = " << mdcFactor[i] << endl;

    // normalize x data
    SetFactors(*x, &MDC.mdcTime, &mdcFactor, cfg->iwindow/2.);	

    // rescale amplitudes stored in the mdcList
    for(int k=0; k<(int)K; k++) {
      int ilog[5] = {1,3,12,13,14};
      for(int l=0;l<5;l++) {
        double mfactor = l<2 ? mdcFactor[k] : mdcFactor[k]*mdcFactor[k];
        TString slog = TB.GetMDCLog(MDC.mdcList[k], ilog[l]);
        MDC.mdcList[k]=TB.SetMDCLog(MDC.mdcList[k], ilog[l], mfactor*slog.Atof());
      }
    }

    // ---------------------------------
    // set mdc list in the network class 
    // ---------------------------------

    if(ifo.CompareTo(net->ifoName[0])==0) {

      net->mdcList.clear();
      net->mdcType.clear();
      net->mdcTime.clear();
      net->mdcList=MDC.mdcList;
      net->mdcType=MDC.mdcType;
      net->mdcTime=MDC.mdcTime;
    }

    // for multi stage analysis store mdc infos to job root file
    if((gJSTAGE!=CWB_STAGE_FULL) && (ifo.CompareTo(net->ifoName[0])==0)) {

      int ifactor;
      std::vector<std::string>*  mdcList = new vector<std::string>;   // list of injections
      std::vector<std::string>*  mdcType = new vector<std::string>;   // list of injection types
      std::vector<double>*       mdcTime = new vector<double>;        // gps time of selected injections

      int inj_size = (int)MDC.mdcList.size(); 

      jfile->cd();

      TTree* mdc_tree = (TTree*)jfile->Get("injections");
      if(mdc_tree==NULL) {
        mdc_tree = new TTree("injections", "injections");
        mdc_tree->Branch("ifactor", &ifactor, "ifactor/I");
        mdc_tree->Branch("mdcList", &mdcList, inj_size, 0);
        mdc_tree->Branch("mdcType", &mdcType, inj_size, 0);
        mdc_tree->Branch("mdcTime", &mdcTime, inj_size, 0);
      } else {
        mdc_tree->SetBranchAddress("ifactor", &ifactor);
        mdc_tree->SetBranchAddress("mdcList", &mdcList);
        mdc_tree->SetBranchAddress("mdcType", &mdcType);
        mdc_tree->SetBranchAddress("mdcTime", &mdcTime);
      }

      ifactor=gIFACTOR;

      *mdcList=MDC.mdcList;
      *mdcType=MDC.mdcType;
      *mdcTime=MDC.mdcTime;

      mdc_tree->Fill();
     
      jfile->Write(); 
    }

    cout.precision(14);
    for(int k=0;k<(int)MDC.mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)MDC.mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)MDC.mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;

  }

  // for multi stage analysis read mdc infos to job root file
  if((gJSTAGE!=CWB_STAGE_FULL) && (type==CWB_PLUGIN_ILIKELIHOOD)) {  

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR) 
    cout << "-----> CWB_Plugin_Sim4.C -> " << " gIFACTOR : " << gIFACTOR << endl;

    int ifactor;
    std::vector<std::string>*  mdcList = new vector<std::string>;   // list of injections
    std::vector<std::string>*  mdcType = new vector<std::string>;   // list of injection types
    std::vector<double>*       mdcTime = new vector<double>;        // gps time of selected injections

    jfile->cd();
    TTree* mdc_tree = (TTree*)jfile->Get("injections");
    if(mdc_tree==NULL) {
      cout << "-----> CWB_Plugin_Sim4.C : Error-> " << " injection tree not present in job file" << endl;
      gSystem->Exit(1);
    } else {
      mdc_tree->SetBranchAddress("ifactor", &ifactor);
      mdc_tree->SetBranchAddress("mdcList", &mdcList);
      mdc_tree->SetBranchAddress("mdcType", &mdcType);
      mdc_tree->SetBranchAddress("mdcTime", &mdcTime);
    }

    int tree_size = (int)mdc_tree->GetEntries(); 

    bool ifactor_check=false;
    for(int i=0;i<tree_size;i++) {
      mdc_tree->GetEntry(i);
      if(ifactor==gIFACTOR) {
        ifactor_check=true;
        net->mdcList.clear();
        net->mdcType.clear();
        net->mdcTime.clear();
        net->mdcList=*mdcList;
        net->mdcType=*mdcType;
        net->mdcTime=*mdcTime;
      } 
    }
    if(!ifactor_check) {
      cout << "-----> CWB_Plugin_Sim4.C : Error-> " << " mdc infos not present in the injection tree " << endl;
      gSystem->Exit(1);
    }

    cout.precision(14);
    for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << net->mdcList[k] << endl;
    for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << net->mdcTime[k] << endl;
    for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << net->mdcType[k] << endl;
  }

  return;
}

// modify input signals (wi) at times pT according the factor pF
size_t 
SetFactors(wavearray<double> &wi, std::vector<double>* pT, std::vector<double>* pF, double dT) {

   int j,nstop,nstrt;
   size_t k;
   double F,T;
   size_t K    = pT->size();
   size_t N    = wi.size();
   size_t M    = 0;
   double rate = wi.rate();                        // simulation rate

   dT = fabs(dT);

   wavearray<double> w; w=wi;

   // isolate injections
   for(k=0; k<K; k++) {

     F = (*pF)[k];
     T = (*pT)[k] - w.start();

     nstrt = int((T - dT)*rate);
     nstop = int((T + dT)*rate);
     if(nstrt<=0) nstrt = 0;
     if(nstop>=int(N)) nstop = N;
     if(nstop<=0) continue;                     // outside of the segment
     if(nstrt>=int(N)) continue;                // outside of the segment

     for(j=nstrt; j<nstop; j++) {
       w.data[j]*=F;
     }

     M++;
   }
   wi = w;
   return M;
}

