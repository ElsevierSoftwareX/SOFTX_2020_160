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
//!EVENT_REPORT
// This plugin save skymap to output root file

  cout << endl;
  cout << "-----> CWB_Plugin_Skymap.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_OLIKELIHOOD) {  

    int i,j,k,n,m,l;
    int nIFO = net->ifoListSize();
    int K = net->nLag;            
    int M = net->mdc__ID.size();  
    int ID;                       
    char search = net->tYPe;                                   
                                                               
    wavearray<double> id;                                      
                                                               
    bool batch = gROOT->IsBatch();                             
    gROOT->SetBatch(true);                                     
                                                               
    watplot SMS(const_cast<char*>("sms"),200,20,800,400);      
                                                               
    char ifostr[20]  = "";                                     
    char strtime[1024]  = "";                                  
    char fname[1024];                                          
                                                               
    int    minTimeDet=nIFO;                                    
    double minTime=1e40;                                       
    double eventTime[NIFO];                                    
    double lagTime[NIFO];                                      
    int ifoid[NIFO],sortid[NIFO];                              
    double factor = 1;                                         


    //Fill in all skymaps
    double old_cc = net->netCC;
    double old_rho = net->netRHO;
    net->netCC = -1;
    net->netRHO = 0;

    // search output root file in the system list
    TFile* froot = NULL;                         
    TList *files = (TList*)gROOT->GetListOfFiles();
    if (files) {                                   
      TIter next(files);                           
      TSystemFile *file;                           
      TString fname;                               
      bool check=false;                            
      while ((file=(TSystemFile*)next())) {        
         fname = file->GetName();                  
         // set output root file as the current file
         if(fname.Contains("wave_")) {froot=(TFile*)file;froot->cd();}
      }                                                               
      if(!froot) {                                                    
        cout << "CWB_Plugin_skymap.C : Error - output root file not found" << endl;
        exit(1);                                                                   
      }                                                                            
    } else {                                                                       
      cout << "CWB_Plugin_skymap.C : Error - output root file not found" << endl;  
      exit(1);                                                                     
    }                                                                              

    TDirectory* cdskymap = (TDirectory*)froot->Get("skymap");
    if(cdskymap==NULL) cdskymap = froot->mkdir("skymap");
 
    netevent* EVT = new netevent(nIFO);
                                       
    int rate = int(2*net->getifo(0)->TFmap.resolution(0)+0.5); 
                                                               
    for(k=0; k<K; k++) {                                // loop over the lags
                                                                             
      id = net->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);        
                                                                             
      for(j=0; j<(int)id.size(); j++) {                 // loop over cluster index
                                                                                  
        ID = size_t(id.data[j]+0.5);                    // cluster id             
                                                                                  
        EVT->output(NULL,net,factor,ID,k);              // get reconstructed parameters
        if(EVT->rho[1] < cfg->cedRHO) continue;         // skip events under rho threshold
                                                                                          
        int masterDet=0;                                                                  
        int lagMin=2147483647;                                                            
        for(n=0; n<nIFO; n++) if(EVT->lag[n]<lagMin) {lagMin=int(EVT->lag[n]);masterDet=n;}
                                                                                           
        net->likelihood(search, net->acor, ID, k);      // exec likelihood search          
                                                                                           
        // get event network time                                                          
        double gps_start = EVT->time[masterDet]-EVT->duration[masterDet];                  
        double gps_stop  = EVT->time[masterDet]+EVT->duration[masterDet];                  

        // create a unique label
        for(n=0; n<nIFO; n++) sprintf(strtime, "%s_%.3f",strtime,EVT->start[n]); 

        //likelihood skymaps                                        
        TMarker inj(EVT->phi[1],90.-EVT->theta[1], 29);  // injected pos     (white star)
        TMarker rec(EVT->phi[0],90.-EVT->theta[0], 29);  // recostructed pos (black star)
        TMarker det(EVT->phi[3],90.-EVT->theta[3], 20);  // detected pos     (black dot )
        inj.SetMarkerSize(2.0); inj.SetMarkerColor(kWhite);                              
        rec.SetMarkerSize(2.0); rec.SetMarkerColor(kBlack);                              
        det.SetMarkerSize(1.0); det.SetMarkerColor(kBlack);                              

        char evtdir[64];
        sprintf(evtdir, "run_%d_id_%d",net->nRun,ID);
        TDirectory* cdevtdir = (TDirectory*)cdskymap->Get(evtdir);
        if(cdevtdir==NULL) cdevtdir = cdskymap->mkdir(evtdir);
        cdevtdir->cd();

        SMS.canvas->cd();
        sprintf(fname, "sensitivity_plus");
        SMS.plot(net->nSensitivity, 0);                         
        rec.Draw(); det.Draw(); if(M) inj.Draw();               
        SMS.canvas->Write(fname);                               
        SMS.clear();                                                                                        
        sprintf(fname, "sensitivity_cross");
        SMS.plot(net->nAlignment, 0);                                                                       
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "skystat");
        SMS.plot(net->nSkyStat, 0);                                                                         
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "likelihood");                                                                       
        SMS.plot(net->nLikelihood, 0);                                                                      
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "null_energy");                                                                      
        SMS.plot(net->nNullEnergy, 0);                                                                      
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "corr_energy");                                                                      
        SMS.plot(net->nCorrEnergy, 0);                                                                      
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "penalty");                                                                          
        SMS.plot(net->nPenalty, 0);                                                                         
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "disbalance");                                                                       
        SMS.plot(net->nDisbalance, 0);                                                                      
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "correlation");                                                                      
        SMS.plot(net->nCorrelation, 0);                                                                     
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "netindex");                                                                         
        SMS.plot(net->nNetIndex, 0);                                                                        
        rec.Draw(); det.Draw(); if(M) inj.Draw();                                                           
        SMS.canvas->Write(fname);                                                                           
        SMS.clear();                                                                                        
        sprintf(fname, "ellipticity");                                                                      
        SMS.plot(net->nEllipticity, 0);                                                                     
        rec.Draw(); det.Draw(); if(M) inj.Draw();
        SMS.canvas->Write(fname);
        SMS.clear();
        sprintf(fname, "polarisation");
        SMS.plot(net->nPolarisation, 0);
        rec.Draw(); det.Draw(); if(M) inj.Draw();
        SMS.canvas->Write(fname);
        SMS.clear();
        sprintf(fname, "probability");
        SMS.plot(net->nProbability, 0);
        rec.Draw(); det.Draw(); if(M) inj.Draw();
        SMS.canvas->Write(fname);
        SMS.clear();

      }                                                 // End loop on found events
    }                                                   // End loop on lags

    // restore NET thresholds
    net->netCC = old_cc;
    net->netRHO = old_rho;
  }

  return;
}

