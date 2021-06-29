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
#include "watplot.hh"
#include <stdio.h>


double GetBoundaries(wavearray<double> x, double P, double& bT, double& eT);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!EVENT_REPORT
// This plugin add to the output root file the whitened recunstructed waveforms 

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_waveform.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;


  if(type==CWB_PLUGIN_OLIKELIHOOD) {  

    int i,j,k,n,m,l;
    int nIFO = NET->ifoListSize();
    int K = NET->nLag;            
    int M = NET->mdc__ID.size();  
    int ID;                       
    char search = NET->tYPe;                                   
  
    wavearray<double> id;
  
    bool batch = gROOT->IsBatch();
    gROOT->SetBatch(true);        
  
    char ifostr[20]  = ""; 
    char strtime[1024]  = ""; 
    char fname[1024];   
  
    int    minTimeDet=nIFO;
    double minTime=1e40;   
    double eventTime[NIFO];
    double lagTime[NIFO];  
    int ifoid[NIFO],sortid[NIFO];   
    double factor = 1;
    detector *pD[NIFO];

 
    netevent* EVT = new netevent(nIFO);

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
        cout << "plugins/CWB_Plugin_waveform.C : Error - output root file not found" << endl;
        exit(1);
      }
    } else {
      cout << "plugins/CWB_Plugin_waveform.C : Error - output root file not found" << endl;
      exit(1);
    }

    int ndim=nIFO;
    int run=NET->nRun;
    int eventID;
    float rho,netcc,frequency,likelihood;
    float snr[NIFO];  
    double time[NIFO];  
    int lag[NIFO+1];  
    int slag[NIFO+1];  
    wavearray<double>* wave[NIFO];
    for(n=0; n<nIFO; n++) wave[n] = new wavearray<double>;
    TTree* wftree = (TTree *) froot->Get("waveform");
    if(wftree==NULL) {
      wftree = new TTree("waveform","waveform");
      wftree->Branch("ndim",&ndim,"ndim/I");
      wftree->Branch("run",&run,"run/I");
      wftree->Branch("eventID",&eventID,"eventID/I");
      wftree->Branch("rho",&rho,"rho/F");
      wftree->Branch("netcc",&netcc,"netcc/F");
      wftree->Branch("likelihood",&likelihood,"likelihood/F");
      char csnr[16];sprintf(csnr,"snr[%1d]/F",nIFO);
      wftree->Branch("snr",snr,csnr);
      wftree->Branch("frequency",&frequency,"frequency/F");
      char ctime[16];sprintf(ctime,"time[%1d]/D",nIFO);
      wftree->Branch("time",time,ctime);
      char clag[16];sprintf(clag,"lag[%1d]/F",nIFO+1);
      char cslag[16];sprintf(cslag,"slag[%1d]/F",nIFO+1);
      wftree->Branch("lag",lag,clag);
      wftree->Branch("slag",slag,cslag);
      char label[256];
      for(n=0; n<nIFO; n++) {
        sprintf(label, "%s_waveform", NET->ifoName[n]);
        wftree->Branch(label,"wavearray<double>",&wave[n],32000,0);
      }
    } else {
      wftree->SetBranchAddress("ndim",&ndim);
      wftree->SetBranchAddress("run",&run);
      wftree->SetBranchAddress("eventID",&eventID);
      wftree->SetBranchAddress("rho",&rho);
      wftree->SetBranchAddress("netcc",&netcc);
      wftree->SetBranchAddress("frequency",&frequency);
      wftree->SetBranchAddress("likelihood",&likelihood);
      wftree->SetBranchAddress("snr",snr);
      wftree->SetBranchAddress("time",time);
      wftree->SetBranchAddress("lag",lag);
      wftree->SetBranchAddress("slag",slag);
      char label[256];
      for(n=0; n<nIFO; n++) {
        wave[n] = new wavearray<double>;
        sprintf(label, "%s_waveform", NET->ifoName[n]);
        wftree->SetBranchAddress(label,&wave[n]);
      }
    }

    for(n=0; n<nIFO; n++) pD[n] = NET->getifo(n);

    int rate = int(2*NET->getifo(0)->TFmap.resolution(0)+0.5); 

    for(k=0; k<K; k++) {  				// loop over the lags
  
      id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);

      for(j=0; j<(int)id.size(); j++) {  		// loop over cluster index
  
        ID = size_t(id.data[j]+0.5);			// cluster id

        eventID = ID;
  
        EVT->output(NULL,NET,factor,ID,k); 		// get reconstructed parameters
        if(EVT->rho[1] < cfg->cedRHO) continue;    	// skip events under rho threshold
 
        rho = EVT->rho[1]; 
        netcc = EVT->netcc[0]; 
        frequency = EVT->frequency[0]; 
        likelihood = EVT->likelihood; 
        for(n=0; n<nIFO; n++)   snr[n] = EVT->sSNR[n];
        for(n=0; n<nIFO; n++)  time[n] = EVT->time[n];
        for(n=0; n<=nIFO; n++)  lag[n] = EVT->lag[n];
        for(n=0; n<=nIFO; n++) slag[n] = EVT->slag[n];
 
        int masterDet=0;
        int lagMin=2147483647;
        for(n=0; n<nIFO; n++) if(EVT->lag[n]<lagMin) {lagMin=int(EVT->lag[n]);masterDet=n;}
  
        NET->likelihood(search, NET->acor, ID, k);	// exec likelihood search
 
        // get event network time
        double gps_start = EVT->time[masterDet]-EVT->duration[masterDet];
        double gps_stop  = EVT->time[masterDet]+EVT->duration[masterDet];

        // create a unique label
        //for(n=0; n<nIFO; n++) sprintf(strtime, "%s_%.3f",strtime,EVT->start[n]); 
        for(n=0; n<nIFO; n++) sprintf(strtime, "%s_%d",strtime,(int)EVT->start[n]); 
 
        // compute event time & lags time
        for(n=0; n<nIFO; n++) eventTime[n]=(EVT->start[n]+EVT->stop[n])/2.;
        for(n=0; n<nIFO; n++) minTime = (eventTime[n]<minTime) ? eventTime[n] : minTime;
        for(n=0; n<nIFO; n++) lagTime[n]=eventTime[n]-minTime;                          
        for(n=0; n<nIFO; n++) minTimeDet = (eventTime[n]<=minTime) ? n : minTimeDet;    

        NET->getwave(ID, k, 'W');
        wavearray<double> w;
        for(n=0; n<nIFO; n++) {
          w.start(0);
          w.rate(pD[n]->waveForm.rate());
          w.resize(pD[n]->waveForm.size());
          for(int i=0;i<w.size();i++) w[i]=pD[n]->waveForm.data[i];

          double bT,eT;
          GetBoundaries(w, 0.99, bT, eT);
          //cout << w.start() << " " << bT << " " << eT << " " << eT-bT << endl;

          // select data which contains the 0.99 of energy
          double dt = 1./w.rate();
          int size = (eT-bT)/dt;
          int os = bT/dt;
          wave[n]->resize(size);
          wave[n]->start(w.start());
          wave[n]->rate(w.rate());
          for(int i=0;i<size;i++) wave[n]->data[i]=w[i+os];
        } 

        wftree->Fill();
      } 						// End loop on found events
    }

    wftree->Write("", TObject::kOverwrite);

    gROOT->SetBatch(batch);  // restore batch status

    jfile->cd(); 

    for(n=0; n<nIFO; n++) delete wave[n];
    delete EVT;  
  }

  if(type==CWB_PLUGIN_CLOSE_JOB) {  
  }

  return;
}


double
GetBoundaries(wavearray<double> x, double P, double& bT, double& eT) {

  if(P<0) P=0;
  if(P>1) P=1;

  int N = x.size();

  double E = 0;                                                 // signal energy
  double avr = 0;                                               // average      
  for(int i=0;i<N;i++) {avr+=i*x[i]*x[i]; E+=x[i]*x[i];}
  int M=int(avr/E);                                             // central index

  // search range which contains percentage P of the total energy E
  int jB=0;
  int jE=N-1;
  double a,b;
  double sum = ((M>=0)&&(M<N)) ? x[M]*x[M] : 0.;
  for(int j=1; j<N; j++) {
    a = ((M-j>=0)&&(M-j<N)) ? x[M-j] : 0.;
    b = ((M+j>=0)&&(M+j<N)) ? x[M+j] : 0.;
    if(a) jB=M-j;
    if(b) jE=M+j;
    sum += a*a+b*b;
    if(sum/E > P) break;
  }

  bT = x.start()+jB/x.rate();
  eT = x.start()+jE/x.rate();

  return eT-bT;
}


