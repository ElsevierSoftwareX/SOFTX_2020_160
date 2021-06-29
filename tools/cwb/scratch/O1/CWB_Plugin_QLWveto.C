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

float GetQveto(wavearray<double>* wf);
void  GetLveto(netcluster* pwc, int cid, int nifo, float* Lveto);
void  GetWveto(netcluster* pwc, int cid, int nifo, float* Wveto);
void  PlotWaveform(TString ifo, wavearray<double>* wfREC,
                   CWB::config* cfg, bool fft=false, bool strain=false);


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Extract whitened reconstructed waveforms, and compute the Qveto, Lveto & Wveto parameters

  cout << endl;
  cout << "-----> CWB_Plugin_QLWveto.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  float Qveto[2*NIFO_MAX];				// Qveto
  float Lveto[3];                               	// Lveto
  float Wveto[2];                               	// Wveto

  if(type==CWB_PLUGIN_CONFIG) {  
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
        cout << "CWB_Plugin_QLWveto.C : Error - output root file not found" << endl;
        gSystem->Exit(1);                                                                             
      }                                                                                      
    } else {                                                                                 
      cout << "CWB_Plugin_QLWveto.C : Error - output root file not found" << endl;  
      gSystem->Exit(1);                                                                               
    }                                                                                        

    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree==NULL) {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("Qveto",Qveto,TString::Format("Qveto[%i]/F",2*cfg->nIFO));
      net_tree->Branch("Lveto",Lveto,TString::Format("Lveto[%i]/F",3));
      net_tree->Branch("Wveto",Wveto,TString::Format("Wveto[%i]/F",2));
    }
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_QLWveto.C -> "
           << "CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_QLWveto.C -> " 
         << " gIFACTOR : " << gIFACTOR << endl;

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
        cout << "CWB_Plugin_QLWveto.C : Error - output root file not found" << endl;
        gSystem->Exit(1);                                                                             
      }                                                                                      
    } else {                                                                                 
      cout << "CWB_Plugin_QLWveto.C : Error - output root file not found" << endl;  
      gSystem->Exit(1);                                                                               
    }                                                                                        

    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree!=NULL) {
      EVT = new netevent(net_tree,nIFO);
      net_tree->SetBranchAddress("Qveto",Qveto);
      net_tree->SetBranchAddress("Lveto",Lveto);
      net_tree->SetBranchAddress("Wveto",Wveto);
    } else {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("Qveto",Qveto,TString::Format("Qveto[%i]/F",2*cfg->nIFO));
      net_tree->Branch("Lveto",Lveto,TString::Format("Lveto[%i]/F",3));
      net_tree->Branch("Wveto",Wveto,TString::Format("Wveto[%i]/F",2));
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
        GetLveto(pwc, ID, nIFO, Lveto);
        cout << "Lveto : " << "fmean : " << Lveto[0] << " frms : " << Lveto[1] 
             << " Energy Ratio : " << Lveto[2] << endl << endl;
        GetWveto(pwc, ID, nIFO, Wveto);
        cout << "Wveto : " << " Slope : " << Wveto[0] << " Correlation : " << Wveto[1] << endl << endl;

        // extract whitened reconstructed waveforms
        for(int n=0; n<nIFO; n++) {

           pd = NET->getifo(n);

           pwfREC[n] = pd->RWFP[wfIndex];
           wavearray<double>* wfREC = pwfREC[n];	// array of reconstructed waveforms

#ifdef PLOT_WHITENED_WAVEFORMS
           //PlotWaveform(NET->ifoName[n], wfREC, cfg, false, false);
           PlotWaveform(NET->ifoName[n], wfREC, cfg, true, false);
#endif
	   // reconstructed whitened waveform
           NET->getMRAwave(ID,k,'S',0,true);
           Qveto[n] = GetQveto(&(pd->waveForm));
	   // whitened waveform
           NET->getMRAwave(ID,k,'W',0,true);
           Qveto[n+nIFO] = GetQveto(&(pd->waveBand));

           //Qveto[n] = GetQveto(wfREC);
           cout << "Qveto : " << pd->Name << " Qveto[R] = " << Qveto[n] 
                                          << " Qveto[W] = " << Qveto[n+nIFO] << endl;
        }
        cout << "----------------------------------------------------------------" << endl;
        delete [] pwfREC;

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
          // add Wveto to dump file
          fprintf(EVT->fP,"Wveto:      ");
          for(int i=0; i<3; i++) fprintf(EVT->fP,"%f ",Wveto[i]);
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

float
GetQveto(wavearray<double>* wf) { 

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
  float Qveto = ein>0 ? eout/ein : 0.;
  //cout << "Qveto : " << Qveto << " ein : " << ein << " eout : " << eout << endl;

  return Qveto;
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
      cout << "CWB_Plugin_QLWveto.C - Error : is enabled only for WDM (2G)" << endl;
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
GetWveto(netcluster* pwc, int cid, int nifo, float* Wveto) {
//          
// input                                                                                      
//        pwc    : pointer to netcluster object                                                          
//        cid    : cluster id                                                                            
//        nifo   : number of detectors      
// output                                                             
//     Wveto[0]  : whistle slope
//     Wveto[1]  : whistle correlation
//   
                                                                                             
  Wveto[0] = Wveto[1] = Wveto[2] = Wveto[3] = Wveto[4] = 0;

  std::vector<int>* vint = &(pwc->cList[cid-1]);        // pixel list
  int V = vint->size();                                 // cluster size
  if(!V) return;                                                       

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> w;
  
  // extract pixels                                             
  for(int j=0; j<V; j++) {                  
    netpixel* pix = pwc->getPixel(cid,j);
    if(pix->layers%2==0) {        
      cout << "CWB_Plugin_QLWveto.C - Error : is enabled only for WDM (2G)" << endl;
      exit(1);                                                                                     
    }                                                                                              
    if(pix->likelihood<1. || pix->frequency==0) continue;
    //if(!pix->core) continue;                  // select only the principal components pixels

    double time = int(double(pix->time)/pix->layers)/pix->rate; // time in seconds from the start
    double freq = pix->frequency*pix->rate/2.;                                    

    x.push_back(time);
    y.push_back(freq);
    w.push_back(pix->likelihood);
  }                                                                                         
  int size = x.size();
  if(size<5) return;                                                                       

  double xcm, ycm, qxx, qyy, qxy, ew;                                                           
  xcm = ycm  = qxx = qyy = qxy = ew = 0;                                                         
                                                                                             
  for(int i=0; i<size; ++i) {                                                                  
    xcm += x[i]*w[i];                                                                           
    ycm += y[i]*w[i];                                                                           
    ew  += w[i]; 
  }                                                                                         
  xcm /= ew;                                                                                
  ycm /= ew;                                                                                
                                                                                             
  for(int i=0; i<size; ++i) {  
    qxx += (x[i] - xcm)*(x[i] - xcm)*w[i];                                                      
    qyy += (y[i] - ycm)*(y[i] - ycm)*w[i];                                                      
    qxy += (x[i] - xcm)*(y[i] - ycm)*w[i];                                                      
  }                                                                                         

  double beta  = qxy/qxx;		// slope
  double alpha = ycm-beta*xcm;		// intercept
  double corr  = qxy/sqrt(qxx*qyy);	// correlation
  double duration = sqrt(qxx/ew);	// duration
  double bandwidth = sqrt(qyy/ew);	// bandwidth

  //printf("alpha : %lf  , beta : %lf , corr : %lf,  dur : %lf , bw : %lf \n", 
  //       alpha, beta, corr, duration, bandwidth);                                                 

  Wveto[0] = beta;			// slope
  Wveto[1] = fabs(corr);		// correlation 

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
