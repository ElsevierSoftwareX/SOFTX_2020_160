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

//#define SAVE_WHT_PLOT			// enable event WHITE plots 

#define nPEAK	10			// number of peaks to be extracted from the reconstructed waveform

void GetGlitchParams(wavearray<double>* wf, int ifoID, float PF[][nPEAK], 
                     float PB[][nPEAK], float PA[][nPEAK], float PE[][nPEAK]);
void PlotWaveform(TString ifo, wavearray<double>* wfREC,
                  CWB::config* cfg, bool fft=false, bool strain=false);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Extract whitened reconstructed waveforms, and compute the first nPEAK parameters
// Save peak parameters to the output wave root file 

  cout << endl;
  cout << "-----> CWB_Plugin_WavePeaks.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  int   nP = nPEAK;			// number of peaks
  float PF[NIFO_MAX][nPEAK];		// peak frequency
  float PB[NIFO_MAX][nPEAK];		// peak bandwidth
  float PA[NIFO_MAX][nPEAK];		// peak amplitude/(max signal amplitude)
  float PE[NIFO_MAX][nPEAK];		// peak energy/(signal energy)

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->outPlugin=true;  				// disable built-in output root file
  }

  if(type==CWB_PLUGIN_ILIKELIHOOD) {
    NET->wfsave=true;                                   // enable waveform save
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_WavePeaks.C -> "
           << "CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_WavePeaks.C -> " 
         << " gIFACTOR : " << gIFACTOR << endl;

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
           //cout << "output file name : " << fname << endl;
         }
      }                                                               
      if(!froot) {                                                    
        cout << "CWB_Plugin_WavePeaks.C : Error - output root file not found" << endl;
        gSystem->Exit(1);                                                                             
      }                                                                                      
    } else {                                                                                 
      cout << "CWB_Plugin_WavePeaks.C : Error - output root file not found" << endl;  
      gSystem->Exit(1);                                                                               
    }                                                                                        

    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree!=NULL) {
      EVT = new netevent(net_tree,nIFO);
      net_tree->SetBranchAddress("nP",&nP);
      net_tree->SetBranchAddress("PF",PF);
      net_tree->SetBranchAddress("PB",PB);
      net_tree->SetBranchAddress("PA",PA);
      net_tree->SetBranchAddress("PE",PE);
    } else {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("nP",&nP,"nP/I");
      net_tree->Branch("PF",PF,TString::Format("PF[%i][%i]/F",cfg->nIFO,nPEAK));
      net_tree->Branch("PB",PB,TString::Format("PB[%i][%i]/F",cfg->nIFO,nPEAK));
      net_tree->Branch("PA",PA,TString::Format("PA[%i][%i]/F",cfg->nIFO,nPEAK));
      net_tree->Branch("PE",PE,TString::Format("PE[%i][%i]/F",cfg->nIFO,nPEAK));
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

        wavearray<double>** pwfREC = new wavearray<double>*[nIFO];
        detector* pd = NET->getifo(0);
        int idSize = pd->RWFID.size();

        int wfIndex=-1;
        for (int mm=0; mm<idSize; mm++) if (pd->RWFID[mm]==ID) wfIndex=mm;
        if(wfIndex==-1) continue;

        // extract whitened reconstructed waveforms
        for(int n=0; n<nIFO; n++) {

           pd = NET->getifo(n);

           pwfREC[n] = pd->RWFP[wfIndex];
           wavearray<double>* wfREC = pwfREC[n];	// array of reconstructed waveforms

#ifdef SAVE_WHT_PLOT
           //PlotWaveform(NET->ifoName[n], wfREC, cfg, false, false);
           PlotWaveform(NET->ifoName[n], wfREC, cfg, true, false);
#endif

           for(int i=0;i<nPEAK;i++) {
             GetGlitchParams(wfREC, n, PF, PB, PA, PE);
//             printf("%d -> \t %5.1f \t %5.1f \t %3.2f \t %3.2f\n",
//                    i,PF[n][i],PB[n][i],PA[n][i],PE[n][i]);
           }
        }
        delete [] pwfREC;

        std::vector<int> sCuts = NET->getwc(k)->sCuts;  // save cCuts
        // set sCuts=1 to the events which must be not copied with cps to pwc
        for(int i=0; i<(int)sCuts.size(); i++) if(i!=ID-1) NET->getwc(k)->sCuts[i]=1;

        // ID can not be used to get the event, to get event use ID=1 (ex: for watplot)
        NET->getwc(k)->sCuts = sCuts;                   // restore cCuts

        if(cfg->dump) EVT->dopen(outDump.Data(),const_cast<char*>("a"),false);
        EVT->output2G(net_tree,NET,ID,k,ofactor);       // get reconstructed parameters
        if(cfg->dump) EVT->dclose();
        if(!cfg->cedDump) NET->getwc(k)->sCuts[ID-1]=1; // mark as processed
      }
    }

    jfile->cd();
    if(EVT) delete EVT;
  }
  return;
}

void
GetGlitchParams(wavearray<double>* wf, int ifoID, float PF[][nPEAK], 
                float PB[][nPEAK], float PA[][nPEAK], float PE[][nPEAK]) {

  int size = wf->size()/2;
  double Fs=((double)wf->rate()/(double)(2*size));

  wavearray<float> ener(size); 

  wf->FFTW(1);
  for(int j=0;j<size;j++) ener[j]=pow(wf->data[2*j],2)+pow(wf->data[2*j+1],2);
  wf->FFTW(-1);
  wf->resetFFTW();

  // compute total signal energy
  double ET=0; 
  for(int j=0;j<size;j++) ET+=ener[j];
  // compute max signal amplitude
  double MA=0; 
  for(int j=0;j<size;j++) if(ener[j]>MA) MA=ener[j]; 
  MA=sqrt(MA);

  for(int k=0;k<nPEAK;k++) {

    // find frequency @ max energy
    float ME=0; 	// peak max energy
    int  jME=0; 	// peak max energy index
    for(int j=0;j<size;j++) if(ener[j]>ME) {ME=ener[j];jME=j;}

    // compute peak energy 
    float  fmean=ME*Fs*jME;
    float  frms=ME*pow(Fs*jME,2);
    double EP=ME;
    double ge=ME;
    ener[jME]=0;
    float finf=0;
    for(int j=jME-1;j>=0;j--) {
      float E = ener[j];
      float F = Fs*j;
      if(E<ge) {EP+=E;finf=F;fmean+=E*F;frms+=E*F*F;ener[j]=0;} else break;
      ge=E;
    }
    ge = ME;
    float fsup=size*Fs;
    for(int j=jME+1;j<size;j++) {
      float E = ener[j];
      float F = Fs*j;
      if(E<ge) {EP+=E;fsup=F;fmean+=E*F;frms+=E*F*F;ener[j]=0;} else break;
      ge=E;
    }

    fmean = fmean/EP;
    frms  = frms/EP-fmean*fmean;
    frms  = frms>0 ? sqrt(frms) : 0;

    PF[ifoID][k] = fmean;
    PA[ifoID][k] = sqrt(ME)/MA;
    PB[ifoID][k] = frms;  
    PE[ifoID][k] = ET>0 ? EP/ET : 0;
  }
}

void 
PlotWaveform(TString ifo, wavearray<double>* wfREC, 
             CWB::config* cfg, bool fft, bool strain) {

  watplot PTS(const_cast<char*>("ptswrc"),200,20,800,500);

  //cout << "Print " << fname << endl;
  double tmin = wfREC->start();
  wfREC->start(wfREC->start()-tmin);
  if(fft) {
    PTS.plot(wfREC, const_cast<char*>("ALP"), 1, 0, 0, true, cfg->fLow, cfg->fHigh);
  } else {
    PTS.plot(wfREC, const_cast<char*>("ALP"), 1, 0, 0);
  }
  PTS.graph[0]->SetLineWidth(1);
  wfREC->start(wfREC->start()+tmin);

  char label[64]="";
  if(fft) sprintf(label,"%s","fft");
  else    sprintf(label,"%s","time");
  if(strain) sprintf(label,"%s_%s",label,"strain");
  else       sprintf(label,"%s_%s",label,"white"); 

  char fname[1024];
  sprintf(fname, "%s_wf_%s_rec_gps_%d.root",ifo.Data(),label,int(tmin));
  PTS.canvas->Print(fname); 
  cout << "write : " << fname << endl;
  //PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
}
