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

//!COHERENCE

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// This plugin implements the built-in coherence stage (only 2G), user can provide its own coherence stage implementation

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_Coherence.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->cohPlugin=true;  	// disable built-in coherence stage
  }

  if(type==CWB_PLUGIN_ICOHERENCE) {

    cout << "type==CWB_PLUGIN_ICOHERENCE" << endl;

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_Sim4.C -> " << " gIFACTOR : " << gIFACTOR << endl;

    int rateANA=cfg->inRate>>cfg->levelR;
    int nIFO = net->ifoListSize();
    detector* pD[NIFO_MAX];                              // pointers to detectors
    for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);
    double mTau=net->getDelay(const_cast<char*>("MAX")); // maximum time delay

    int upN = rateANA/1024; if(upN<1) upN=1;              // calculate upsample factor
    int nRES = net->wdmListSize();			 // number of resolution levels

    double Eo;
    netcluster* pwc;

    for(int i=0; i<nRES; i++) {                                 // loop over TF resolutions
      // print level infos                                                                 
      int level=cfg->l_high-i;                                                              
      int layers = level>0 ? 1<<level : 0;                                                 
      int rate  = rateANA>>level;                                                          
      cout << "level : " << level << "\t rate(hz) : " << rate                              
           << "\t layers : " << layers << "\t df(hz) : " << rateANA/2./double(1<<level)    
           << "\t dt(ms) : " << 1000./rate << endl;                                        
  
      // produce TF maps with max over the sky energy
      for(int n=0; n<nIFO; n++) {                    
        WDM<double>* pwdm = net->wdmList[i];
        wavearray<double>* hot = net->getifo(n)->getHoT();
        net->getifo(n)->getTFmap()->maxEnergy(*hot,*pwdm,mTau,upN);
        //if(singleDetector) {                                            
        //  *(net->getifo(1)->getTFmap()) = *(net->getifo(0)->getTFmap());  
        //  break;                                                        
        //}
        // restore the frequency boundaries changed by the maxEnergy call
        net->getifo(n)->getTFmap()->setlow(cfg->fLow);                     
        net->getifo(n)->getTFmap()->sethigh(cfg->fHigh);                   
      }
  
      Eo = net->THRESHOLD(cfg->bpp);                        // threshold on pixel energy
      cout<<"thresholds in units of noise variance: Eo="<<Eo<<" Emax="<<Eo*2<<endl;   
  
      double TL = net->setVeto(cfg->iwindow);
      cout<<"live time in zero lag: "<<TL<<endl;          		// set veto array
      if(TL <= 0.) {cout<<"livetime is zero : exit"<<endl;exit(1);}	// exit if live time is zero
  
      // select pixels
      if(cfg->simulation) {cout<<"ifactor|clusters|pixels ";cout.flush();}
      else                {cout<<"lag|clusters|pixels ";    cout.flush();}
      int csize_tot=0;int psize_tot=0;                                   
      for(int j=0; j<(int)net->nLag; j++) {                               

         net->getNetworkPixels(j,Eo,Eo*2);
         net->cluster(1,1);               
         pwc = net->getwc(j);             

         // store cluster into temporary job file
         int cycle = cfg->simulation ? gIFACTOR : Long_t(pwc->shift);
         pwc->write(jfile,"coherence","clusters",0,cycle);         
         pwc->write(jfile,"coherence","clusters",-1,cycle);        
         cout<<cycle<<"|"<<pwc->csize()<<"|"<<pwc->size()<<" ";cout.flush();
         csize_tot+=pwc->csize(); psize_tot+=pwc->size();                   

         pwc->clear();
      }
    }
  }

  return;
}
