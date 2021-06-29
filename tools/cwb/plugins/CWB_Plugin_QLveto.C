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

//#define PLOT_LIKELIHOOD
//#define PLOT_WHITENED_WAVEFORMS

#define NTHR 1   		
#define ATHR 7.58859   		

double GetQfactor(wavearray<double>* wf, double frequency, bool fixAmax);
double GetPeakFrequency(wavearray<double>* wf);
void   GetGaussFitPars(wavearray<double>* wf, double& mean, double& sigma, bool doenv=true);
void   GetGaussFitPars2(wavearray<double>* wf, double& mean, double& sigma, bool fixAmax);

void  GetQveto(wavearray<double>* wf, float &Qveto, float &Qfactor);
void  GetLveto(netcluster* pwc, int cid, int nifo, float* Lveto);
void  PlotWaveform(TString ifo, wavearray<double>* wfREC,
                   CWB::config* cfg, bool fft=false, bool strain=false);
void ClearWaveforms(detector* ifo);

std::vector<netpixel> DoPCA(network* NET, CWB::config* cfg, int lag, int id);
void ResetPCA(network* NET, CWB::config* cfg, netcluster* pwc, std::vector<netpixel>* vPIX, int ID);


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Extract whitened reconstructed waveforms, and compute the Qveto, Lveto parameters

//  cout << endl;
//  cout << "-----> CWB_Plugin_QLveto.C" << endl;
//  cout << "ifo " << ifo.Data() << endl;
//  cout << "type " << type << endl;
//  cout << endl;

  float Qveto[4];					// Qveto
  float Lveto[3];                               	// Lveto

  if(type==CWB_PLUGIN_CONFIG) {  
    cout << endl;
    cout << "-----> CWB_Plugin_QLveto.C" << endl;
    cout << endl;
    cfg->outPlugin=true;  				// disable built-in output root file
  }

  if(type==CWB_PLUGIN_ILIKELIHOOD) {
    NET->wfsave=true;                                   // enable waveform save

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
        cout << "CWB_Plugin_QLveto.C : Error - output root file not found" << endl;
        gSystem->Exit(1);                                                                             
      }                                                                                      
    } else {                                                                                 
      cout << "CWB_Plugin_QLveto.C : Error - output root file not found" << endl;  
      gSystem->Exit(1);                                                                               
    }                                                                                        

    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree==NULL) {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("Qveto",Qveto,TString::Format("Qveto[%i]/F",4));
      net_tree->Branch("Lveto",Lveto,TString::Format("Lveto[%i]/F",3));
    } else {
      TBranch* branch;
      bool qveto_exists=false;
      bool lveto_exists=false;
      TIter next(net_tree->GetListOfBranches());
      while ((branch=(TBranch*)next())) {
        if(TString("Qveto").CompareTo(branch->GetName())==0) qveto_exists=true;
        if(TString("Lveto").CompareTo(branch->GetName())==0) lveto_exists=true;
      }
      next.Reset();
      if(!qveto_exists) net_tree->Branch("Qveto",Qveto,TString::Format("Qveto[%i]/F",4));
      if(!lveto_exists) net_tree->Branch("Lveto",Lveto,TString::Format("Lveto[%i]/F",3));
    }
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_QLveto.C -> "
           << "CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    //cout << "-----> CWB_Plugin_QLveto.C -> " 
    //     << " gIFACTOR : " << gIFACTOR << endl;

    // import slagShift
    float* gSLAGSHIFT=NULL; IMPORT(float*,gSLAGSHIFT)

    int nIFO = NET->ifoListSize();			// number of detectors
    int K = NET->nLag;  				// number of time lag          
    netevent* EVT;
    wavearray<double> id;
    //double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];                 
    double factor = cfg->factors[gIFACTOR];                 
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
        cout << "CWB_Plugin_QLveto.C : Error - output root file not found" << endl;
        gSystem->Exit(1);                                                                             
      }                                                                                      
    } else {                                                                                 
      cout << "CWB_Plugin_QLveto.C : Error - output root file not found" << endl;  
      gSystem->Exit(1);                                                                               
    }                                                                                        

    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree!=NULL) {
      EVT = new netevent(net_tree,nIFO);
      net_tree->SetBranchAddress("Qveto",Qveto);
      net_tree->SetBranchAddress("Lveto",Lveto);
    } else {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("Qveto",Qveto,TString::Format("Qveto[%i]/F",4));
      net_tree->Branch("Lveto",Lveto,TString::Format("Lveto[%i]/F",3));
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

        netcluster* pwc = NET->getwc(k);
        cout << endl << "----------------------------------------------------------------" << endl;

        // extract whitened reconstructed waveforms
        Qveto[0]=Qveto[1]=1.e20;
        double Qfactor1[NIFO_MAX];
        double Qfactor2[NIFO_MAX];
        for(int n=0; n<nIFO; n++) {

           pd = NET->getifo(n);

           pwfREC[n] = pd->RWFP[wfIndex];
           wavearray<double>* wfREC = pwfREC[n];	// array of reconstructed waveforms

#ifdef PLOT_WHITENED_WAVEFORMS
           //PlotWaveform(NET->ifoName[n], wfREC, cfg, false, false);
           PlotWaveform(NET->ifoName[n], wfREC, cfg, true, false);
#endif
           // store in Qveto[0/1] the minimum velue of Qveto/Qfactor between the reconstructed,whitened waveforms in all ifos
           float qveto,qfactor;

	   // reconstructed whitened waveform
           NET->getMRAwave(ID,k,'S',0,true);
           GetQveto(&(pd->waveForm), qveto, qfactor);
           if(qveto<Qveto[0]) Qveto[0]=qveto; 
           if(qfactor<Qveto[1]) Qveto[1]=qfactor; 

	   Qfactor1[n] = qfactor;
           double fpeak = GetPeakFrequency(&(pd->waveForm));
           Qfactor2[n] = GetQfactor(&(pd->waveForm), fpeak, true);      // get qfactor with amax fixed;
           cout << endl;
           cout << "Qfactor1/2 : " << pd->Name << " frequency[0] = " << EVT->frequency[0] << " fpeak = " << fpeak 
		<< " -> " << Qfactor1[n] << " / " << Qfactor2[n] << endl;
           cout << endl;

	   // whitened waveform
           NET->getMRAwave(ID,k,'W',0,true);
           GetQveto(&(pd->waveBand), qveto, qfactor);
           if(qveto<Qveto[0]) Qveto[0]=qveto; 
           if(qfactor<Qveto[1]) Qveto[1]=qfactor; 

           cout << "Qveto : " << pd->Name << " Qveto[0] = " << Qveto[0] 
                                          << " Qveto[1] = " << Qveto[1] << endl;

           if(!cfg->simulation) ClearWaveforms(pd);	// release waveform memory
        }
        delete [] pwfREC;

        // compute Qveto[2] using the weighted of Qfactor1 values
        Qveto[2]=0.;
        for(int n=0; n<nIFO; n++) Qveto[2]+=EVT->sSNR[n]*Qfactor1[n];
        Qveto[2]/=EVT->likelihood;
        // compute Qveto[3] using the weighted of Qfactor2 values
        Qveto[3]=0.;
        for(int n=0; n<nIFO; n++) Qveto[3]+=EVT->sSNR[n]*Qfactor2[n];
        Qveto[3]/=EVT->likelihood;
        cout << "Qveto   :" << " Qveto[2] = " << Qveto[2]
                            << " Qveto[3] = " << Qveto[3] << endl;

        std::vector<netpixel> vPIX;
        if(cfg->pattern>0) vPIX = DoPCA(NET, cfg, k, ID);            // do PCA analysis
        GetLveto(pwc, ID, nIFO, Lveto);
        if(cfg->pattern>0) ResetPCA(NET, cfg, pwc, &vPIX, ID);	     // restore WP pwc

	cout << endl;
        cout << "Lveto : " << "fmean : " << Lveto[0] << " frms : " << Lveto[1] 
             << " Energy Ratio : " << Lveto[2] << endl << endl;
        cout << "----------------------------------------------------------------" << endl;

        std::vector<int> sCuts = NET->getwc(k)->sCuts;  // save cCuts
        // set sCuts=1 to the events which must be not copied with cps to pwc
        for(int i=0; i<(int)sCuts.size(); i++) if(i!=ID-1) NET->getwc(k)->sCuts[i]=1;

        // ID can not be used to get the event, to get event use ID=1 (ex: for watplot)
        NET->getwc(k)->sCuts = sCuts;                   // restore cCuts

        if(cfg->dump) EVT->dopen(outDump.Data(),const_cast<char*>("a"),false);
        EVT->output2G(net_tree,NET,ID,k,ofactor);       // get reconstructed parameters
        if(cfg->dump) {                                 
          // add Qveto to dump file
          fprintf(EVT->fP,"Qveto:      ");
          for(int i=0; i<2*nIFO; i++) fprintf(EVT->fP,"%f ",Qveto[i]);
          fprintf(EVT->fP,"\n");
          // add Lveto to dump file
          fprintf(EVT->fP,"Lveto:      ");
          for(int i=0; i<3; i++) fprintf(EVT->fP,"%f ",Lveto[i]);
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

void
GetQveto(wavearray<double>* wf, float &Qveto, float &Qfactor) { 

  wavearray<double> x = *wf;; 

  // resample data by a factor 4
  int xsize=x.size();
  x.FFTW(1);
  x.resize(4*x.size());
  x.rate(4*x.rate());
  for(int j=xsize;j<x.size();j++) x[j]=0;
  x.FFTW(-1);

  // extract max/min values and save the absolute values in the array 'a'
  wavearray<double> a(x.size());
  int size=0;
  double dt = 1./x.rate();
  double prev=x[0];
  double xmax=0;
  for (int i=1;i<x.size();i++) {
    if(fabs(x[i])>xmax) xmax=fabs(x[i]);
    if(prev*x[i]<0) {
      a[size]=xmax;
      size++;
      xmax=0;
    }
    prev=x[i];
  }

  // find max value/index ans save on  amax/imax  
  int imax=-1;
  double amax=0;
  for (int i=1;i<size;i++) {
    if(a[i]>amax) {amax=a[i];imax=i;}
  }

/*
  cout << endl;
  cout << "a[imax-2] " << a[imax-2] << endl;
  cout << "a[imax-1] " << a[imax-1] << endl;
  cout << "a[imax]   " << a[imax] << endl;
  cout << "a[imax+1] " << a[imax+1] << endl;
  cout << "a[imax+2] " << a[imax+2] << endl;
  cout << endl;
*/

  // compute Qveto 
  double ein=0;	// energy of max values inside NTHR
  double eout=0;	// energy of max values outside NTHR
  for (int i=0;i<size;i++) {
    if(abs(imax-i)<=NTHR) {
      ein+=a[i]*a[i];
      //cout << i << " ein " << a[i] << " " << amax << endl;
    } else {
      if(a[i]>amax/ATHR) eout+=a[i]*a[i];
      //if(a[i]>amax/ATHR) cout << i << " eout " << a[i] << " " << amax << endl;
    }
  }
  Qveto = ein>0 ? eout/ein : 0.;
  //cout << "Qveto : " << Qveto << " ein : " << ein << " eout : " << eout << endl;

  // compute Qfactor 
  float R = (a[imax-1]+a[imax+1])/a[imax]/2.;
  Qfactor = sqrt(-pow(TMath::Pi(),2)/log(R)/2.);
  //cout << "Qfactor : " << Qfactor << endl;

  return;
}

void 
GetLveto(netcluster* pwc, int cid, int nifo, float* Lveto) {
//          
// input                                                                                      
//        pwc    : pointer to netcluster object                                                          
//        cid    : cluster id                                                                            
//        nifo   : number of detectors      
// output                                                             
//     Lveto[0]  : line frequency
//     Lveto[1]  : line bandwitdh
//     Lveto[2]  : line enery ratio (line_energy / total_energy)
//   
                                                                                             
  Lveto[0] = Lveto[1] = Lveto[2] = 0;

  std::vector<int>* vint = &(pwc->cList[cid-1]);        // pixel list
  int V = vint->size();                                 // cluster size
  if(!V) return;                                                       

  // ------------------------------------------------------------------
  // Find max pixel parameters                                          
  // ------------------------------------------------------------------

  double likeMax=0;	// maximum pixel's energy 
  double likeTot=0;	// total cluster energy
  double freqMax;	// frequency of the pixel with max energy
  double dfreqMax;	// df of the pixel with max energy
  for(int n=0; n<V; n++) {
    netpixel* pix = pwc->getPixel(cid,n);
    if(pix->layers%2==0) {        
      cout << "CWB_Plugin_QLveto.C - Error : is enabled only for WDM (2G)" << endl;
      exit(1);                                                                                     
    }                                                                                              
    if(!pix->core) continue;                  // select only the principal components pixels

    double likePix=0;                                                                       
    for(int m=0; m<nifo; m++) {                                                          
      likePix += pow(pix->getdata('S',m),2);  // like whitened reconstructed signal 00
      likePix += pow(pix->getdata('P',m),2);  // like whitened reconstructed signal 90
    }                                                                                     

    double freq = pix->frequency*pix->rate/2.; 
    double df = pix->rate/2.;

    likeTot+=likePix;
    if(likePix>likeMax) {likeMax=likePix;freqMax=freq;dfreqMax=df;}
  }
  //cout << "likeMax : " << likeMax << " likeTot : " << likeTot 
  //     << " freqMax : " << freqMax << " dfreqMax : " << dfreqMax << endl;

  // ------------------------------------------------------------------
  // Compute Lveto parameters                                          
  // ------------------------------------------------------------------

  double fmean=0;	// line mean frequency
  double frms=0;	// line bandwidth	
  double likeLin=0;	// line energy
  for(int n=0; n<V; n++) {
    netpixel* pix = pwc->getPixel(cid,n);
    if(!pix->core) continue;                  // select only the principal components pixels

    double likePix=0;                                                                       
    for(int m=0; m<nifo; m++) {                                                          
      likePix += pow(pix->getdata('S',m),2);  // like whitened reconstructed signal 00
      likePix += pow(pix->getdata('P',m),2);  // like whitened reconstructed signal 90
    }                                                                                     

    // the estimation is done for all pixels 
    // with freq in the range [freqMax-dfreqMax, freqMax+dfreqMax]
    double freq = pix->frequency*pix->rate/2.; 
    if(fabs(freq-freqMax)<=dfreqMax) {
      likeLin += likePix;
      fmean   += likePix*freq;
      frms    += likePix*freq*freq;
    }
  }                                                                   
   
  fmean = fmean/likeLin;
  frms  = frms/likeLin-fmean*fmean;
  frms  = frms>0 ? sqrt(frms) : 0.;

  if(frms<dfreqMax/2.) frms=dfreqMax/2.;

  // ------------------------------------------------------------------
  // Save Lveto parameters                                          
  // ------------------------------------------------------------------

  Lveto[0] = fmean;     			// line mean frequency     
  Lveto[1] = frms;        			// line bandwidth   
  Lveto[2] = likeTot>0. ? likeLin/likeTot : 0.;	// energy ratio energy inside_line/total

  // ------------------------------------------------------------------
  // plot time-frequency energy                                             
  // ------------------------------------------------------------------

#if defined PLOT_LIKELIHOOD 
  watplot WTS(const_cast<char*>("wts"));
  WTS.plot(pwc, cid, nifo, 'L', 0, const_cast<char*>("COLZ"));
  WTS.canvas->Print("l_tfmap_scalogram.png");
#endif

  return;
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
  //sprintf(fname, "%s_wf_%s_rec_gps_%d.root",ifo.Data(),label,int(tmin));
  sprintf(fname, "%s_wf_%s_rec_gps_%d.png",ifo.Data(),label,int(tmin));
  PTS.canvas->Print(fname); 
  cout << "write : " << fname << endl;
  //PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
}

void 
ClearWaveforms(detector* ifo) {

  int n;

  n = ifo->IWFP.size();
  for (int i=0;i<n;i++) {
    wavearray<double>* wf = (wavearray<double>*)ifo->IWFP[i];
    delete wf;
  }
  ifo->IWFP.clear();
  ifo->IWFID.clear();

  n = ifo->RWFP.size();
  for (int i=0;i<n;i++) {
    wavearray<double>* wf = (wavearray<double>*)ifo->RWFP[i];
    delete wf;
  }
  ifo->RWFP.clear();
  ifo->RWFID.clear();
}

std::vector<netpixel> DoPCA(network* NET, CWB::config* cfg, int lag, int id) {

  double ee;

  size_t nIFO = NET->ifoList.size();

  float  En = 2*NET->acor*NET->acor*nIFO;     // network energy threshold in the sky loop

  int size = NET->a_00.size();
  int f_ = NIFO/4;

  netpixel* pix;
  netcluster* pwc = NET->getwc(lag);
  std::vector<netpixel> vPIX;

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, id, false);   // buffer for pixel IDs
  int V = pI.size();
  if(V>cfg->BATCH) return vPIX;                                 // attach TD amp to pixels < V

  wavearray<float>  xi(size); xi=0;           // PC 00 array
  wavearray<float>  XI(size); XI=0;           // PC 90 array

  __m128* _xi  = (__m128*) xi.data;           // set pointer to PC 00 array
  __m128* _XI  = (__m128*) XI.data;           // set pointer to PC 90 array

  __m128* _aa  = (__m128*) NET->a_00.data;    // set pointer to 00 array
  __m128* _AA  = (__m128*) NET->a_90.data;    // set pointer to 90 array

  int nPC = 0;
  for(int j=0; j<V; j++) {
    int jf = j*f_;                            // source sse pointer increment
    _sse_zero_ps(_xi+jf);                     // zero MRA amplitudes
    _sse_zero_ps(_XI+jf);                     // zero MRA amplitudes
    ee = _sse_abs_ps(_aa+jf,_AA+jf);          // total pixel energy / quadrature
    if(ee>En) nPC++; else ee=0.;              // count core pixels
    NET->rNRG.data[j] = ee;                   // init residual energy array
    NET->pNRG.data[j] = NET->rNRG.data[j];    // init residual energy array
  }

  nPC = NET->_sse_mra_ps(xi.data,XI.data,En,nPC);  // get principal components

  for(int j=0; j<V; j++) {                    // loop over principal components
     pix = pwc->getPixel(id,pI[j]);
     vPIX.push_back(*pix);		      // save original pixels
     pix->core = false;
     ee = NET->pNRG.data[j];                  // total pixel energy
     if(ee<En) continue;
     pix->core = true;
     for(int i=0; i<nIFO; i++) {
        pix->setdata(double(xi.data[j*NIFO+i]),'S',i);    // store 00 whitened response PC
        pix->setdata(double(XI.data[j*NIFO+i]),'P',i);    // store 90 whitened response PC
     }
  }

  return vPIX;
}

void ResetPCA(network* NET, CWB::config* cfg, netcluster* pwc, std::vector<netpixel>* vPIX, int ID) {

  std::vector<int> pI = NET->wdmMRA.getXTalk(pwc, ID, false);   // buffer for pixel IDs
  int V = pI.size();
  if(V>cfg->BATCH) return;                                 	// attach TD amp to pixels < V
  for(int j=0; j<V; j++) {
    netpixel* pix = pwc->getPixel(ID,pI[j]);
    *pix = (*vPIX)[j];
  }

  while(!vPIX->empty()) vPIX->pop_back();
  vPIX->clear(); std::vector<netpixel>().swap(*vPIX);
}

double
GetQfactor(wavearray<double>* wf, double frequency, bool fixAmax) {
//
// wf        -> input waveform
// frequency -> input peak frequency
// fixAmax   -> true/false : fix/not-fix amplitude of the gaussian fit
//
// -> the input waveform wf is resampled 4x
// -> get mean and sigma from a gaussian fit of the waveform envelope
// -> compute the range [xmin,xmax] = [mean-3*sigma,mean-3*sigma]]
// -> compute area of the envelope in the rage [xmin,xmax]
// -> compute qfactor using area, max_amplitude and frequency
// -> return qfactor
//

  wavearray<double> x = *wf;;

  // resample data by a factor 4
  int xsize=x.size();

  x.FFTW(1);
  x.resize(4*x.size());
  x.rate(4*x.rate());
  for(int j=xsize;j<x.size();j++) x[j]=0;
  x.FFTW(-1);

  // compute range [xmin,xmax] used for final fit using the square of x
  double gmean, gsigma;
  GetGaussFitPars2(&x, gmean, gsigma, fixAmax);
//  cout << "GetQfactor -> gmean " << gmean << " gsigma " << gsigma << endl;

  double xmin = gmean-3.0*gsigma;
  double xmax = gmean+3.0*gsigma;
//  cout << "GetQfactor -> xmin " << xmin << " xmax " << xmax << endl;

  x = CWB::mdc::GetEnvelope(&x);      // get envelope

  // extract max/min values and save the absolute values in the array 'a'
  double dt = 1./x.rate();
  // find max value/index and save on  amax/imax  
  int imax=-1;
  double amax=0;
  for (int i=1;i<x.size();i++) {
    if(x[i]>amax) {amax=x[i];imax=i;}
  }

  double area=0.;
  for (int i=1;i<x.size();i++) {
    double t=i*dt;
    if(t>xmin && t<xmax) area+=x[i]*dt;
  }

  double sigma = area/amax/sqrt(2*TMath::Pi());
//  cout << "sigma integral is " << sigma << endl;
  double qfactor = sigma*(TMath::TwoPi()*frequency);    // compute qfactor

  return qfactor;
}

double
GetPeakFrequency(wavearray<double>* wf) {
//
// wf        -> input waveform
//
// -> multiply the input waveform wf by its envelope
// -> perform FFT transform
// -> compute max frequency -> peak
// -> return frequency peak
//

  wavearray<double> x  = *wf;
  wavearray<double> ex = CWB::mdc::GetEnvelope(&x);      // get envelope
  int N = x.size();
  for(int i=0;i<N;i++) x[i]*=ex[i];

  x.FFTW(1);

  // fill time/amp graph
  double df=(x.rate()/(double)N)/2.;
  double amax=0.;
  double fmax=0.;
  for(int i=0;i<N;i+=2) {
    double freq=i*df;
    double amp=x[i]*x[i]+x[i+1]*x[i+1];
    if(amp>amax) {fmax=freq;amax=amp;}
  }

  return fmax;
}

void
GetGaussFitPars(wavearray<double>* wf, double& mean, double& sigma, bool doenv) {
//
// wf          -> input waveform
// doenv       -> true/false do/not-do waveform envelope
//
// mean, sigma -> output gaussian fit parameters
//
// -> get evelope from the input waveform wf and perform fit with gaussian function
// -> return mean and gigma of gaussian fit
//

  wavearray<double> x = doenv ? CWB::mdc::GetEnvelope(wf) : *wf;      // get envelope

  // fill time/amp graph
  int N = x.size();
  double* time = new double[N];
  double* amp = new double[N];
  double dt = 1./x.rate();
  for(int i=0;i<N;i++) {
    time[i]=i*dt;
    amp[i]=x[i];
  }
  TGraph gr(N,time,amp);

  // fit signal envelope with gaussian function
  gr.Fit("gaus","q0");
  TF1* gaus = (TF1*)gr.GetFunction("gaus");

  mean = gaus->GetParameter("Mean");           // extract sigma from gauss fit
  sigma = gaus->GetParameter("Sigma");         // extract sigma from gauss fit

  delete [] time;
  delete [] amp;

  return;
}

void
GetGaussFitPars2(wavearray<double>* wf, double& mean, double& sigma, bool fixAmax) {
//
// wf          -> input waveform
// fixAmax     -> true/false : fix/not-fix amplitude of the gaussian fit
//
// mean, sigma -> output gaussian fit parameters
//
//
// -> get the square of the evelope from the input waveform wf 
// -> perform fit with gaussian function and get mean and gigma of gaussian fit
// -> compute the range [xmin,xmax] = [mean-3*sigma,mean-3*sigma]]
//
// -> get the evelope from the input waveform wf 
// -> perform fit with gaussian function in the range [xmin,xmax] fixing the amplitude to max envelope (according to fixAmax value)
// -> return mean and gigma of gaussian fit
//

  wavearray<double> x = CWB::mdc::GetEnvelope(wf);      // get envelope

  // compute range [xmin,xmax] used for final fit using the square of x
  // x*x highlight the main peak
  wavearray<double> x2 = x;                             // getx*x
  for(int i=0;i<x.size();i++) x2[i]=x[i]*x[i];

  double gmean, gsigma;
  GetGaussFitPars(&x2, gmean, gsigma, false);
//  cout << "GetGaussFitPars2 -> gmean x2 " << gmean << " gsigma x2 " << gsigma << endl;

  double xmin = gmean-3.0*gsigma;
  double xmax = gmean+3.0*gsigma;
//  cout << "GetGaussFitPars2 -> xmin " << xmin << " xmax " << xmax << endl;

  // fill time/amp graph
  int N = x.size();
  double* time = new double[N];
  double* amp = new double[N];
  double dt = 1./x.rate();
  double amax=0;
  for(int i=0;i<N;i++) {
    time[i]=i*dt;
    amp[i]=x[i];
    if(fabs(x[i])>amax) amax=fabs(x[i]);
  }
  TGraph gr(N,time,amp);
//  cout << "GetGaussFitPars2 -> amax " << amax << endl;

  // fit signal envelope with gaussian function
  TF1 *xgaus = new TF1("xgaus","gaus",xmin,xmax);
  if(fixAmax) xgaus->SetParLimits(0,amax*0.999,amax*1.001);     // fix gauss constant par to amax
  gr.Fit("xgaus","RQ");

  mean = xgaus->GetParameter("Mean");           // extract sigma from gauss fit
  sigma = xgaus->GetParameter("Sigma");         // extract sigma from gauss fit

  delete xgaus;
  delete [] time;
  delete [] amp;

  return;
}

