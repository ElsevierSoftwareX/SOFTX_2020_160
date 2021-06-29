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
#include "cwb2G.hh"
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
#include "regression.hh"
#include <vector>


// regression parameters
#define REGRESSION_FILTER_LENGTH        8
#define REGRESSION_MATRIX_FRACTION      0.95
#define REGRESSION_SOLVE_EIGEN_THR      0.
#define REGRESSION_SOLVE_EIGEN_NUM      10
#define REGRESSION_SOLVE_REGULATOR      'h'
#define REGRESSION_APPLY_THR            0.8

//#define SAVE_RMS

// ---------------------------------------------------------------------------------
// WHAT IS?
// this plugin can be used for low latency analysis 
// HOW TO CONFIGURE THE Low Latency PLUGIN
// the following is an example : must be included in the config/user_parameters.C file
// see the DumpUserOptions function for the full description of parameters
// ---------------------------------------------------------------------------------
/*
  plugin = TMacro(gSystem->ExpandPathName("$HOME_CWB/plugins/CWB_Plugin_LowLatency.C"));        // Macro source

  TString optll = "";                           // NOTE : add space at the end of each line
  optll += "ll_rmsfile=rmsfile.root ";          // is the root file which contains the noise rms used for whitening data
  optll += "ll_grms=true ";                     // if grms==true  && rmsfile!="" -> RMS whitening is provided by the user by the rmsfile parameter 
  optll += "ll_fitsfile=fitsfile.fits ";        // is the name of file where the fits skymap is saved -> fitsfile.fits 

  strcpy(parPlugin,optll.Data());               // set LL plugin parameters
  strcpy(comment,"LL configuration example");

*/

// ---------------------------------------------------------------------------------
// DEFINES
// ---------------------------------------------------------------------------------

#define LL_RMSFILE              ""              // is the root file which contains the noise rms used for whitening data
                                                //              
#define LL_FITSFILE             ""              // if fitsfile!="" the cWB fits skymap is saved into fitsfile.fits
                                                //              
#define LL_GRMS                 false           // if grms==true  && rmsfile!="" -> RMS whitening is provided by the user by the rmsfile parameter 
                                                // if grms==false && rmsfile!="" -> RMS whitening is saved into rmsfile

// ---------------------------------------------------------------------------------
// USER CONFIG OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  TString rmsfile;
  bool    grms;
  TString fitsfile;
};

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);
void DumpUserOptions(TString odir, CWB::config* cfg);

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions        gOPT;                           // global User Options

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void DumpSkymap(network* NET, int lag, netevent* EVT, int ID);


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!DATA_CONDITIONING
// This plugin shows how to implement the low latency pipeline

  cout << endl;
  cout << "-----> CWB_Plugin_LowLatency.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  

    ResetUserOptions();                                 // set default config options
    ReadUserOptions(cfg->parPlugin);                    // user config options : read from parPlugin

    cfg->dcPlugin=true;   				// disable built in dc 
  }

  if(type==CWB_PLUGIN_NETWORK) {
    PrintUserOptions(cfg);                              // print config options
  }

  if(type==CWB_PLUGIN_DATA_CONDITIONING) {  

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_LowLatency.C -> implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // get data range
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)

    wavearray<double> y;
    WSeries<double> wM;
    // import global variables
    size_t gRATEANA; IMPORT(size_t,gRATEANA)

    // get ifo id
    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "CWB_Plugin_LowLatency.C : Error - bad ifo id" << endl; gSystem->Exit(1);}

    detector* pD = net->getifo(id);
    WSeries<double>* pTF = pD->getTFmap();
    wavearray<double>* hot = (wavearray<double>*)x;

    int layers = 0;
    for(int level=cfg->l_high; level>=cfg->l_low; level--) {
      for(int j=0;j<net->wdmMRA.nRes;j++)
         if(layers<net->wdmMRA.layers[j]) layers=net->wdmMRA.layers[j];
    }
    WDM<double> WDMwhite(layers,layers,6,10);         // set whitening WDM
    layers = gRATEANA/8;
    WDM<double> WDMlpr(layers,layers,6,10);           // set LPE filter WDM

    // regression
    pTF->Forward(*hot,WDMlpr);
    regression rr(*pTF,const_cast<char*>("target"),1.,cfg->fHigh);
    rr.add(*hot,const_cast<char*>("target"));
    rr.setFilter(REGRESSION_FILTER_LENGTH);
    rr.setMatrix(net->Edge,REGRESSION_MATRIX_FRACTION);
    rr.solve(REGRESSION_SOLVE_EIGEN_THR,REGRESSION_SOLVE_EIGEN_NUM,REGRESSION_SOLVE_REGULATOR);
    rr.apply(REGRESSION_APPLY_THR);
    *hot = rr.getClean();

    // whitening
    pTF->Forward(*hot,WDMwhite);
    pTF->setlow(cfg->fLow);
    pTF->sethigh(cfg->fHigh);

    if(gOPT.grms) {

      if(gOPT.rmsfile!="") {						// RMS whitening is loaded from root file

        cout << "CWB_Plugin_LowLatency.C -> LOAD WHITENING FILE: " << gOPT.rmsfile << endl; 

        TFile *frms = new TFile(gOPT.rmsfile);
        if(frms==NULL) {
          cout << "CWB_Plugin_LowLatency.C -> Failed to open RMS whitening file " <<  gOPT.rmsfile << endl;
          gSystem->Exit(1);
        }

        frms->ls();

        WSeries<double>* wrms = (WSeries<double>*)frms->Get(ifo);
        if(wrms==NULL) {
          cout << "CWB_Plugin_LowLatency.C -> Object RMS doesn't exist !!! " <<  endl;
          gSystem->Exit(1);
        }
        pD->nRMS = *wrms;

        int levels = pD->nRMS.getLevel()+1;             // number of levels
        int slices = pD->nRMS.size()/levels;            // number of nRMS samples 
        double length = slices*cfg->whiteStride;      	// nRMS len in sec

	double tcenter = (x->start()+x->stop())/2.; 
	double rms_start = tcenter - length/2.;

        pD->nRMS.start(tcenter);
//        pD->nRMS.start(rms_start);
        delete wrms;
      }

    } else {

      pD->white(cfg->whiteWindow,0,cfg->segEdge,cfg->whiteStride);   	// calculate noise rms 
      pD->nRMS.bandpass(16.,0.,1);                                   	// high pass filtering at 16Hz

      if(gOPT.rmsfile!="") {						// RMS whitening is saved to root file

        cout << "CWB_Plugin_LowLatency.C -> SAVE WHITENING FILE: " << gOPT.rmsfile << endl; 

        // save nRMS
        TFile *frms = id==0 ? new TFile(gOPT.rmsfile, "RECREATE") : new TFile(gOPT.rmsfile, "UPDATE");
        if(frms==NULL) {
          cout << "CWB_Plugin_LowLatency.C -> Failed to create file !!! " <<  endl;
          gSystem->Exit(1);
        }
        pD->nRMS.Write(ifo);
        frms->Close();
      }
    }

    pTF->white(pD->nRMS,1);                                        // whiten  0 phase WSeries
    pTF->white(pD->nRMS,-1);                                       // whiten 90 phase WSeries

    WSeries<double> wtmp = *pTF;
    pTF->Inverse();
    wtmp.Inverse(-2);
    *hot = *pTF;
    *hot += wtmp;
    *hot *= 0.5;
    // add infos to history
    char info[256];
    sprintf(info,"-IFO:%d-RMS:%g",id,hot->rms());
    gCWB2G->PrintAnalysisInfo(CWB_STAGE_CSTRAIN,"cwb2G::DataConditioning",info,false);
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {                    // AFTER EVENT RECONSTRUCTION 

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];
    double ofactor=0;
    if(cfg->simulation==4)      ofactor=-gIFACTOR;
    else if(cfg->simulation==3) ofactor=-gIFACTOR;
    else                        ofactor=factor;

    int nIFO = net->ifoListSize();                      // number of detectors
    int K = net->nLag;                                  // number of time lag
    int rate = 0;                                       // select all resolutions
    netevent* EVT = new netevent(nIFO);;
    wavearray<double> id;

    for(int k=0; k<K; k++) {                            // loop over the lags

      id = net->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);

      for(int j=0; j<(int)id.size(); j++) {             // loop over cluster index

        int ID = size_t(id.data[j]+0.5);

        if(net->getwc(k)->sCuts[ID-1]!=-1)  continue;   // skip rejected/processed clusters

        EVT->output2G(NULL,net,ID,k,ofactor);           // get reconstructed parameters

        // print event parameters
        cout << endl;
        cout << "CWB_Plugin_WF - event parameters : ID -> " << ID << endl;
        for(int n=0;n<nIFO;n++) printf("rec_time %s : %.4f\n",net->ifoName[n], EVT->time[n]);
        cout << "rec_theta : " << EVT->theta[0] << " rec_phi : " << EVT->phi[0] << endl;
        cout << "SNRnet : " << sqrt(EVT->likelihood) << " netcc[0] : " << EVT->netcc[0]
             << " rho[0] : " << EVT->rho[0] << " size : " << EVT->size[0] << endl;

        if(gOPT.fitsfile!="") DumpSkymap(net, k, EVT, ID);     // save skymap to fits file
      }
    }
  }

  return;
}

void PrintUserOptions(CWB::config* cfg) {

  cout << "-----------------------------------------"     << endl;
  cout << "LowLatency config options                "     << endl;
  cout << "-----------------------------------------"     << endl << endl;
  cout << "LL_RMSFILE         " << gOPT.rmsfile           << endl;
  cout << "LL_GRMS            " << gOPT.grms              << endl;
  cout << "LL_FITSFILE        " << gOPT.fitsfile          << endl;

  cout << endl;
}

void DumpUserOptions(TString odir, CWB::config* cfg) {

  TString ofName = odir+"/ll_config.txt";

  ofstream out;
  out.open(ofName,ios::out);
  if(!out.good()) {cout << "DumpUserOptions : Error Opening File : " << ofName << endl;exit(1);}

  TString info="";
  TString tabs="\t\t\t\t";

  char version[128];
  sprintf(version,"WAT Version : %s - SVN Revision : %s - Tag/Branch : %s",watversion('f'),watversion('r'),watversion('b'));

  out << endl;
  out << "--------------------------------"     << endl;
  out << "LL config options               "     << endl;
  out << "--------------------------------"     << endl;
  out << endl;

  out << "ll_grms            " << gOPT.grms << endl;
  out << tabs << info << endl;

  out << "ll_rmsfile         " << gOPT.rmsfile << endl;
  out << tabs << info << endl;

  out << "ll_fitsfile        " << gOPT.fitsfile << endl;
  out << tabs << info << endl;

  out.close();
}

void ResetUserOptions() {

  gOPT.rmsfile   = LL_RMSFILE;
  gOPT.grms      = LL_GRMS;
  gOPT.fitsfile  = LL_FITSFILE;
}

void ReadUserOptions(TString options) {

  if(options.CompareTo("")!=0) {
    cout << options << endl;
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet

      TObjArray* token = TString(options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("ll_rmsfile=")) {
          TString rmsfile=stok;
          rmsfile.Remove(0,rmsfile.Last('=')+1);
          gOPT.rmsfile=rmsfile;
        }

        if(stok.Contains("ll_grms=")) {
          TString grms=stok;
          grms.Remove(0,grms.Last('=')+1);
          gOPT.grms=grms;
        }

        if(stok.Contains("ll_fitsfile=")) {
          TString fitsfile=stok;
          fitsfile.Remove(0,fitsfile.Last('=')+1);
          gOPT.fitsfile=fitsfile;
        }
      }
    }
  }
}

void
DumpSkymap(network* NET, int lag, netevent* EVT, int ID) {

  // Dump2fits probability skymap  (healpix)
  skymap skyprobcc = NET->getifo(0)->tau;
  skyprobcc=0.;
  skymap skyprob = skyprobcc;
  skyprob=1.e-12;

  std::vector<float>* vP;
  std::vector<int>*   vI;

  vP = &(NET->wc_List[lag].p_Map[ID-1]);
  vI = &(NET->wc_List[lag].p_Ind[ID-1]);
  double th,ph,ra;
  int k;
  for(int j=0; j<int(vP->size()); j++) {
    int i = (*vI)[j];
    th = skyprob.getTheta(i);
    ph = skyprob.getPhi(i);

    k=skyprob.getSkyIndex(th, ph);
    skyprob.set(k,(*vP)[j]);

    ra = skyprob.getRA(i);
    k=skyprob.getSkyIndex(th, ra);
    skyprobcc.set(k,(*vP)[j]);
  }

  skyprobcc.Dump2fits(const_cast<char*>(gOPT.fitsfile.Data()),EVT->time[0],const_cast<char*>(""),const_cast<char*>("PROB"),const_cast<char*>("pix-1"),'C');

  return;
}

