/*
# Copyright (C) 2019 Sergey Klimenko
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
#include "gwavearray.hh"
#include "regression.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "Toolbox.hh"

// ---------------------------------------------------------------------------------
// HOW TO SETUP PLUGIN IN CWB USER CONFIGURATION (EXAMPLE)
// ---------------------------------------------------------------------------------

/*

  TString opt_o3ac = "";                  // NOTE : add space at the end of each line

  opt_o3ac += "o3ac_ifo=enable ";         // enable 32hz of first ifo
  opt_o3ac += "o3ac_ifo=enable ";         // enable 32hz of seconf ifo
  opt_o3ac += "o3ac_ifo=disable ";        // disable 32hz of third ifo

  strcpy(parPlugin,opt_o3ac.Data());      // set plugin parameters

*/

// ---------------------------------------------------------------------------------
// USER O3aConditioning PLUGIN OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  bool   ifo[NIFO_MAX];
};

uoptions gOPT;                                  // global User Options

void ResetUserOptions(network* net, CWB::config* cfg);
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!32Hz_CONDITIONING
// Plugin to correct PSD variability in the 16-48Hz band  

  if(type==CWB_PLUGIN_NETWORK) {
    ResetUserOptions(net, cfg);                 // set default config options
    ReadUserOptions(cfg->parPlugin);            // user config options : read from parPlugin
    PrintUserOptions(cfg);                      // print config options
  }

  if(type==CWB_PLUGIN_IDATA_CONDITIONING) {
     // Marek put power 60 up-conversion regression here
  }
  
  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {
     cout << endl;
     cout << "-----> CWB_Plugin_O3aDataConditioning.C ";
     cout << "type " << type << endl;
     cout << endl;

    if(TString(cfg->analysis)=="1G") {
      cout << "CWB_Plugin_O3aConditioning Error - this plugin works only with 2G" << endl;
      gSystem->Exit(1);
    }
    
    // get ifo id
    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) {
       if(ifo==net->ifoName[n]) {id=n; break;}
    }
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    if(!gOPT.ifo[id]) return;                     // skip 32hz correction if ifo is not enabled

    wavearray<double> u,v,w,p,q,s,c;
    WSeries<double> ws, WS;
    detector* pD = net->getifo(id);
    wavearray<double>* hot = pD->getHoT();
    int M = int(hot->rate()/64+0.1);              // number of WDM layers with 32 Hz bandwidth    
    WDM<double> wdm(M,M,6,10);
    double edge = net->Edge;
     
    ws.Forward(*hot,wdm);
    ws.getLayer(c, 1); 
    ws.getLayer(s,-1);
    u = c; w = c;
    for(int i=0; i<u.size(); i++) {
       u.data[i]=sqrt(c.data[i]*c.data[i]+s.data[i]*s.data[i]);
       c.data[i]/=u.data[i]>0. ? u.data[i] : 1;
       s.data[i]/=u.data[i]>0. ? u.data[i] : 1;
    }
    v = u;
    double R = v.rate();                           // WDM layer rate (should be 64 Hz)
    double T = v.stop()-v.start();                 // data duration      
    v.lprFilter(4.,0,0.,edge,1);                   // apply lpe filter to clean layer 
    double um = u.median(edge*R,(T-edge)*R);       // median of u - original
    double vm = v.median(edge*R,(T-edge)*R);       // median of v - clean
    double vr,aa;
    
    v += um-vm;
    for(int i=0; i<u.size(); i++) {
       if(v.data[i]<0.) v.data[i]=0.;              // corrected magnitude can not be negative
       vr = u.data[i]-v.data[i];                   // rms variability
       aa = u.data[i]<2*um?(u.data[i]-um)/um:1.;   // amplitude correction factor
       if(vr<0 || aa<0) aa = 0;                    // no correction
       vr *= vr<1 ? vr : 1;                        // quadratic correction for vr<1
       v.data[i] = u.data[i]-vr*aa;                // apply all corrections
       w.data[i] = v.data[i]/u.data[i];
    }

    int n = int(edge*R);
    int m = int(0.6*R);

    p=v; q=v; q=0.;
    for(int i=0; i<u.size(); i++) {
     p.data[i] = 1-w.data[i];
     if(p.data[i]>0.9) p.data[i]=0.;
    }
    for(int i=n+m; i<u.size()-n-m; i++) {
       for(int j=0; j<m; j++) {
	  q.data[j] += p.data[i]*(p.data[i-j]+p.data[i+j]);
       }
    }
    q *= 1/q.data[0]; 
    
    pD->nVAR.resize(u.size());
    pD->nVAR = 1;
    
    for(int i=m+1; i<u.size()-m-1; i++) {
       double sp = 0;
       double sm = 0;
       for(int j=8; j<m; j++) {
	  aa = p.data[i-j]; sm += aa*aa*q.data[j];
	  aa = p.data[i+j]; sp += aa*aa*q.data[j];
       }
       aa = sp<sm ? sm : sp;                      
       aa = 1./(1.+aa);                             // multiplicity correction 
       v.data[i] *= aa;                             // data
       w.data[i] *= aa;                             // variability
       if(w.data[i]==0) w.data[i]=1.e-16;           // 23/08/19 fix nRMS = 0
       pD->nVAR.data[i] = w.data[i];                // noise rms correction
    }

    c*=v; ws.putLayer(c, 1); 
    s*=v; ws.putLayer(s,-1);
    WS = ws;
    ws.Inverse();   c = ws;
    WS.Inverse(-2); s = WS;
    *hot = c; *hot += s; *hot *= 0.5;              // whitened and corected time series
    pD->nVAR.rate(R);
    pD->nVAR.start(hot->start());
    pD->nVAR.setlow(16.);
    pD->nVAR.sethigh(48.);

    /*
    char ofname[128];
    sprintf(ofname,"%s_%d_variability.root",ifo.Data(),net->nRun);
    TFile *froot = new TFile(ofname, "RECREATE");
    gwavearray<double> gw=w; 
    gw.Write("32hz");
    froot->Close();
    */
  }

  return;
}

void ReadUserOptions(TString options) {

  int n_ifo=0;
  if(options.CompareTo("")!=0) {
    cout << options << endl;
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet

      TObjArray* token = TString(options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("o3ac_ifo=")) {
          TString o3ac_ifo=stok;
          o3ac_ifo.Remove(0,o3ac_ifo.Last('=')+1);
          if(o3ac_ifo=="enable")  gOPT.ifo[n_ifo]=true;
          if(o3ac_ifo=="disable") gOPT.ifo[n_ifo]=false;
          if(n_ifo<(NIFO_MAX-1)) n_ifo++;
        }
      }
    }
  }
}

void PrintUserOptions(CWB::config* cfg) {

    cout << "-----------------------------------------"     << endl;
    cout << "O3aConditioning config options           "     << endl;
    cout << "-----------------------------------------"     << endl << endl;
    for(int n=0;n<cfg->nIFO;n++) {
      if(gOPT.ifo[n]==true)  cout << "O3aC_IFO  " << cfg->ifo[n] << "   enabled"  << endl;
      if(gOPT.ifo[n]==false) cout << "O3aC_IFO  " << cfg->ifo[n] << "   disabled" << endl;
    }

    cout << endl;
}

void ResetUserOptions(network* net, CWB::config* cfg) {

    for(int n=0;n<cfg->nIFO;n++) {
       gOPT.ifo[n]=false;
       if(TString(net->ifoName[n])=="L1") gOPT.ifo[n]=true;
       if(TString(net->ifoName[n])=="H1") gOPT.ifo[n]=true;
    }
}

