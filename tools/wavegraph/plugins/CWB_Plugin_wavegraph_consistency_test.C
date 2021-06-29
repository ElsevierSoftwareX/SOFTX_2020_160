/*
# Copyright (C) 2019 Eric Chassande-Mottin, Philippe Bacon, Gayathri V, Archana Pai
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
#include "TPolyMarker.h"
#include "mdc.hh"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iterator>
#include "wavegraph.hh"
#include "wavepath.hh"
#include "watplot.hh"

//!COHERENCE

//#define WAVEGRAPH_FILE  "config/wavegraph_cbc_1d4x1d4_to_25x25.txt"
//#define WAVEGRAPH_FILE  "config/wavegraph_2.5x2.5_25x25_newton_log_r290-fixed.txt"
//#define WAVEGRAPH_FILE  "config/wavegraph_2.5x2.5_25x25_newton_log_3481_0d2e.txt"
//#define WAVEGRAPH_FILE  "config/wavegraph_2.5x2.5_10x10_newton_log_strip_AdvVirgo.txt"
#define WAVEGRAPH_FILE  "config/wavegraph_2.5x2.5_25x25_newton_log_20150310.txt"
#define WAVEGRAPH_PLOT_WAVEPATH   // plot wavepath on top of the TF maps
#define WAVEGRAPH_THR 1.0 
#define WAVEGRAPH_WIDTH 0 
#define WAVEGRAPH_PENAL_FACTOR 1.0
#define WAVEGRAPH_WDMDUMP "mydumpWDM" // dump folder from user's home (to do not add path to home!)
// #undef WAVEGRAPH_WDMDUMP // uncomment to disable dump of WDM data to file

// wavecraft coherence
void wcoherence(TFile* jfile, CWB::config* cfg, network* net);	
void getNetworkEnergyMap(int LAG, network* net, WSeries<double>& MAP, int& nM, int* IN, int& jB, int& jE);   
void PlotWDM(WSeries<double>* WS, wavearray<double>* t, wavearray<double>* f);

WSeries<double> edMAP[NIFO_MAX][NRES_MAX]; 	// detector energy TF maps
vector<double> energy_thresholds; // energy thresholds at each resolution level

#define NWG_PARAMS 1
float WavegraphParams[NWG_PARAMS];  // Wavegraph parameters sent to cWB output

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// This plugin implements the wavegraph in coherence stage (only 2G)

  cout << endl;
  cout << "-----> CWB_Plugin_wavegraph.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->cohPlugin=true;  	// disable built-in coherence stage
    cfg->outPlugin=true;  	// disable built-in output root file

  }

  if(type==CWB_PLUGIN_ICOHERENCE) {

    cout << "type==CWB_PLUGIN_ICOHERENCE" << endl;

    int rateANA=cfg->inRate>>cfg->levelR;
    int nIFO = net->ifoListSize();
    detector* pD[NIFO_MAX];                              // pointers to detectors
    for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);
    double mTau=net->getDelay(const_cast<char*>("MAX")); // maximum time delay

    int upN = rateANA/1024; if(upN<1) upN=1;              // calculate upsample factor
    int nRES = net->wdmListSize();			 // number of resolution levels

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

        // copy TF energy map to edMAP
        edMAP[n][i] = *(net->getifo(n)->getTFmap()); 
        edMAP[n][i].stop(net->getifo(n)->getTFmap()->stop());        
      } 

      double Eo = net->THRESHOLD(cfg->bpp);                      // threshold on pixel energy
      energy_thresholds.push_back(Eo);
      cout.precision(5);
      cout<<"thresholds in units of noise variance: Eo="<<Eo;
    }

    wcoherence(jfile, cfg, net);	// wavegraph coherence	

    // New sparse map handler
    jfile->Close();			// flush clusters to jfile
    TString jname = jfile->GetName();	// get job file name
    cout << "jname : " << jname << endl;
    jfile = new TFile(jname, "UPDATE");
    // import gCWB2G
    cwb2G* gCWB2G=NULL; IMPORT(cwb2G*,gCWB2G)
    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_wavegraph.C -> gIFACTOR : " << gIFACTOR << endl;
    // write sparse table to job file
    gCWB2G->WriteSparseTFmap(jfile, gIFACTOR, "csparse", "coherence");
    jfile->Close();			// flush sparse maps to jfile

  }

  return;
}

  if(type==CWB_PLUGIN_ILIKELIHOOD) {

    // search output root file in the system list
    TFile* froot = GetRootFile(NET);

    netevent* EVT;
    int nIFO = NET->ifoListSize();                      // number of detectors
    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree==NULL) {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("WavegraphParams",WavegraphParams,TString::Format("WavegraphParams[%i]/F",NWG_PARAMS));
    }
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_netEvent.C -> "
           << "CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_netEvent.C -> "
         << " gIFACTOR : " << gIFACTOR << endl;

    // import slagShift
    float* gSLAGSHIFT=NULL; IMPORT(float*,gSLAGSHIFT)

    int nIFO = NET->ifoListSize();                      // number of detectors
    int K = NET->nLag;                                  // number of time lag
    netevent* EVT;
    wavearray<double> id;
    //double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];
    double factor = cfg->factors[gIFACTOR];
    int rate = 0;                                       // select all resolutions

    // search output root file in the system list
    TFile* froot = GetRootFile(NET);

    TString outDump = froot->GetName();
    outDump.ReplaceAll(".root.tmp",".txt");

    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree!=NULL) {
      EVT = new netevent(net_tree,nIFO);
      net_tree->SetBranchAddress("WavegraphParams",WavegraphParams);
    } else {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
      net_tree->Branch("WavegraphParams",WavegraphParams,TString::Format("WavegraphParams[%i]/F",NWG_PARAMS));
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

        netcluster* pwc = NET->getwc(k);

        // ** WavegraphParams[0]=111; **//
        // ** WavegraphParams[1]=222; **//

        std::vector<int> sCuts = NET->getwc(k)->sCuts;  // save cCuts
        // set sCuts=1 to the events which must be not copied with cps to pwc
        for(int i=0; i<(int)sCuts.size(); i++) if(i!=ID-1) NET->getwc(k)->sCuts[i]=1;

        // ID can not be used to get the event, to get event use ID=1 (ex: for watplot)
        NET->getwc(k)->sCuts = sCuts;                   // restore cCuts

        if(cfg->dump) EVT->dopen(outDump.Data(),const_cast<char*>("a"),false);
        EVT->output2G(net_tree,NET,ID,k,ofactor);       // get reconstructed parameters
        if(cfg->dump) {                                 
          // add WavegraphParams to dump file
          fprintf(EVT->fP,"WavegraphParams:      ");
          for(int i=0; i<NWG_PARAMS; i++) fprintf(EVT->fP,"%f ",WavegraphParams[i]);
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



//**************************************************************************
// wavegraph coherence stage
//**************************************************************************
void wcoherence(TFile* jfile, CWB::config* cfg, network* net) {

  // import ifactor
  int gIFACTOR=-1; IMPORT(int,gIFACTOR)
  cout << "-----> CWB_Plugin_wavegraph.C -> " << " gIFACTOR : " << gIFACTOR << endl;

  int rateANA=cfg->inRate>>cfg->levelR; //  New sparse map handler

  int nIFO = net->ifoListSize();		// number of detectors
  int nRES = net->wdmListSize();		// number of resolution levels

  netcluster* pwc;
  netpixel pix(nIFO); 

  // network energy maps (sum over the detectors)
  WSeries<double> enMAP[NRES_MAX]; 	
  // store energy TF maps pointers to data vector
  vector< WSeries<double>* > data;
  for(int i=0; i<nRES; i++) data.push_back(&enMAP[i]);

  // check that graph file exists
  if (access( WAVEGRAPH_FILE, F_OK ) == -1) {
    std::cerr << "CWB_Plugin_wavegraph.C : Error - graph file " << WAVEGRAPH_FILE << " not found\n";
    gSystem->Exit(1);
  }
    
  // load graph
  wavegraph graph;
  graph.create(WAVEGRAPH_FILE);
  if(!graph.is_topologically_sorted()) {
    std::cerr << "CWB_Plugin_wavegraph.C : Error - the graph is not topologically sorted\n";
    gSystem->Exit(1);
  }

  // set injection time window (Tinj +/- iwindow/2)
  double TL = net->setVeto(cfg->iwindow);
  cout<<"live time in zero lag: "<<TL<<endl;          		
  // exit if live time is zero
  if(TL <= 0.) {cout<<"livetime is zero : exit"<<endl;gSystem->Exit(1);}	

  if(cfg->simulation) {cout<<"ifactor|clusters|pixels|paths ";cout.flush();}
  else                {cout<<"lag|clusters|pixels|paths ";    cout.flush();}
  
  // select pixels
  int csize_tot=0;int psize_tot=0;                                   
  for(int j=0; j<(int)net->nLag; j++) {                     // loop over time lags

     int nM[NRES_MAX];
     int IN[NRES_MAX][NIFO_MAX]; 			
     int jB[NRES_MAX]; 			
     int jE[NRES_MAX]; 			
 
     for(int i=0; i<nRES; i++) {                            // loop over TF resolutions
       for(int n=0; n<nIFO; n++) {                    
         *(net->getifo(n)->getTFmap()) = edMAP[n][i]; 
         net->getifo(n)->getTFmap()->stop(edMAP[n][i].stop());
       }
       int ii = nRES-i-1;
       getNetworkEnergyMap(j,net,enMAP[ii],nM[ii],IN[ii],jB[ii],jE[ii]);
     }
     
#ifdef WAVEGRAPH_WDMDUMP
     // quick and dirty dump of data cube for further inspection
     char wname[1024]; static int count=0;
     
     // set path and name of the output file 
     sprintf(wname,"%s/%s/WDM_%10.0f_%d.root",getenv("HOME"),WAVEGRAPH_WDMDUMP,enMAP->start(),count); count++;
     cout << "XXX write raw WDM to file:" << wname;

     TFile *wroot = new TFile(wname, "CREATE");         // create output root file for wdm TF maps
     
     for(int i=0; i<nRES; i++) {                            // loop over TF resolutions                                                     
       int level=cfg->l_low+i;
       char wwlabel[32];sprintf(wwlabel,"map:%d",level);
       enMAP[i].Write(wwlabel);
     }

     cout << " -- done XXX" << endl; cout.flush();
     
     wroot->Close();
#endif

     // Apply wavegraph
     std::vector<cluster> clusters = graph.clustering(WAVEGRAPH_THR,data,cfg->segEdge,WAVEGRAPH_WIDTH,WAVEGRAPH_PENAL_FACTOR,energy_thresholds);

     cout << "Number of paths found by wavegraph : " << clusters.size() << endl;

#ifdef WAVEGRAPH_PLOT_WAVEPATH
     // get best cluster
     cluster best_cluster=clusters.front();
     int M = best_cluster.size();

     // plot best wavepath
     wavearray<double> t(M);
     wavearray<double> f(M);
     wavearray<double> scale(M);
     for (int n = M; n-- > 0;) {
       t[n]=best_cluster[n].time;
       f[n]=best_cluster[n].freq;
       scale[n]=best_cluster[n].log2scale;
       //cout << "t=" << t[n] << " sec, f=" << f[n] << " Hz, scale=" << scale[n] << endl;
       //cout << "ixt=" << best_cluster[n].timeix << " , ixf=" << best_cluster[n].freqix 
       //     << " ixscale=" << best_cluster[n].scaleix << endl;
     }
     //PlotWDM(data.at(4), &t, &f);
#endif

     // fill netwok pixels list
     int nclusters=0;
     int npixels=0;
     for(int i=0; i<nRES; i++) {                            // loop over TF resolutions

       int nPix=0;
       for (int k=0; k<clusters.size(); k++) {		      // loop over clusters 

         int N = clusters[k].size();

         for (int n = N; n-- > 0;) {

           int timeix  = clusters[k][n].timeix;
           int freqix  = clusters[k][n].freqix;
           int scaleix = clusters[k][n].scaleix; 
  
           if(i!=scaleix) continue;

           double R  = enMAP[scaleix].wrate();               // pixel layer rate        
           int I  = enMAP[scaleix].maxLayer()+1;             // number of layers        

           //cout << "ixt=" << timeix << " , ixf=" << freqix << " ixscale=" 
           //     << scaleix << " R " << R << " I " << I << endl;

           pix.core = true;    
           pix.rate = R;       
           pix.layers = I;     

           double E=0;
           for(int m=0; m<nIFO; m++) {
              //cout << i << " " << m << " mN " << nM[i] << " IN[i] " << IN[i][m] 
              //     << " jE " << jE[i] << " jB " << jB[i] << " " << R << endl;
              int kk = IN[i][m]+timeix-jB[i];
              if(kk >= jE[i]) kk -= (jE[i]-jB[i]); // circular buffer
              int jj = I*kk+freqix; 
              pix.data[m].index = jj;                                 
              pix.data[m].asnr = edMAP[m][nRES-scaleix-1].data[jj];                        
              E += edMAP[m][nRES-scaleix-1].data[jj];                                      
           }                                                         
           int kk = IN[i][nM[i]]+timeix-jB[i];
           if(kk >= jE[i]) kk -= (jE[i]-jB[i]);	 // circular buffer
           int jj = I*kk+freqix; 
           pix.time = jj;                                                       
           pix.frequency = freqix;                                                  
           pix.likelihood = E;                                                 
           net->wc_List[j].append(pix);          // save pixels in wc_List
           nPix++;
         }
       } 
 
       if(nPix) net->setRMS();		 // associare noise rms to the selected pixels 

       net->cluster(1,1);               	 // cluster pixels
       pwc = net->getwc(j);             

       // store cluster into temporary job file
       int cycle = cfg->simulation ? gIFACTOR : Long_t(pwc->shift);
       pwc->write(jfile,"coherence","clusters",0,cycle);         
       pwc->write(jfile,"coherence","clusters",-1,cycle,-(rateANA>>(cfg->l_low+i)));   // New sparse map handler
       //       pwc->write(jfile,"coherence","clusters",-1,cycle);        
       nclusters+=pwc->csize();
       npixels+=pwc->size();
       csize_tot+=pwc->csize(); psize_tot+=pwc->size();                   
       
       pwc->clear();				 // clean pixels
     }
     
     // New sparse map handler
     // write clusters to the job file & free trees memory
     jfile->Write();
     TList* fList = gDirectory->GetList();
     TObject *obj;
     TIter nextobj(fList);
     while ((obj = (TObject *) nextobj())) {
       if(TString(obj->GetName()).Contains("clusters")) delete obj;
     }

    int cycle = cfg->simulation ? gIFACTOR : Long_t(pwc->shift);
    cout<<cycle<<"|"<<nclusters<<"|"<<npixels<<"|"<<clusters.size()<<" ";cout.flush();
  }
  cout << endl;
}

//**************************************************************************
//:select TF samples using the wavegraph algorithm   
//**************************************************************************
void getNetworkEnergyMap(int LAG, network* net, WSeries<double>& MAP, int& nM, int* IN, int& jB, int& jE)   
{                                                                           
// return the incoherent energy sum over the time shifted detectors network        
// works with WDM/wavelet energy TF maps                                    
//
// LAG - time shift lag defining how detectors are shifted wrt each other.  
// net - network pointer                
// MAP - incoherent energy sum over the time shifted detectors network              
// nM  - index of the first detector
// IN  - array (size = nIFO) of the time indices  
// jB  - number of samples in the edges 
// jE  - last good sample in the layer
//

   size_t nIFO = net->ifoList.size();       		// number of detectors

   if(nIFO>NIFO) {
      cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(): " 
          <<"invalid number of detectors or\n";
      gSystem->Exit(1);                                                                          
   }                                           
   if(net->getifo(0)->getTFmap()->w_mode != 1) {    
      cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(): invalid whitening mode.\n"; 
      gSystem->Exit(1);                                                                          
   }                                                                  

   WSeries<double>* pTF = net->getifo(0)->getTFmap();  // pointer to first TF map
   MAP = *pTF; MAP=0.;                                 // initialize TF map      
   wavearray<double>* hTS = net->getifo(0)->getHoT();  // pointer to first TS data
                                                                                  
   int i,j,k,m,n,NN,jj,jb,je,J,K;

   double R  = pTF->wrate();                      // pixel layer rate        
   double r  = hTS->rate();                       // TS rate                 
   int N  = pTF->size();                          // size of TF array        
   int M  = hTS->size();                          // size of TS array        
   int I  = pTF->maxLayer()+1;                    // number of layers        
   int II = pTF->maxLayer()-1;                    // number of layers - 2    
       jB = int(net->Edge*R+0.001);               // number of samples in the edges
   if(jB&1) {
      cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(1): WDM parity violation\n"; 
      gSystem->Exit(1);
   }             
                                                                                        
   if(jB < 3) {                                                                         
      cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(): insufficient data edge length.\n";            
      gSystem->Exit(1);                                                                          
   }                                                                                    

   netpixel pix(nIFO); 
   pix.core = true;    
   pix.rate = R;       
   pix.layers = I;     
                       
   int     in[NIFO];                              // pixel time index
   double* PDATA;                                                          
   double* pmap;                                                           
   double* pdata[NIFO];                           // pointers to data
   double* pp[5];                                 // pointers to sorted F-arrays
   for(n=0; n<nIFO; n++) {                        // pointers to data           
      pdata[n] = net->getifo(n)->getTFmap()->data;                                         
   }                                                                                  

   size_t count = 0;                              // live pixel counter  
   double a,b,E,Ct,Cb,Ht,Hb;                                             

   if(net->veto.size() != M) {                    // set veto array if it is not set
      net->veto.resize(M); net->veto = 1;                                                     
   }                                                                                
   short* pveto = net->veto.data;                 // pointer to veto                
                                                                                    
   net->wc_List[LAG].clear();                     // clear netcluster structure     
   net->livTime[LAG] = 0.;                        // clear live time counters       
   net->wc_List[LAG].setlow(pTF->getlow());                                        
   net->wc_List[LAG].sethigh(pTF->gethigh());                                      

   a  = 1.e10; nM = 0;                            // master detector    
   for(n=0; n<nIFO; n++) {                                              
      b = net->getifo(n)->lagShift.data[LAG];     // shift in seconds   
      if(a>b) { a = b; nM = n; }                                        
   }                                                                    
                                                                        
   for(n=0; n<nIFO; n++) {                                              
      b = net->getifo(n)->lagShift.data[LAG];     // shift in seconds   
      K = int((b-a)*R+0.001);                     // time shift wrt reference
      if(K&1) {
        cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(2): WDM parity violation\n"; 
        gSystem->Exit(1);
      }
      in[n] = IN[n] = K+jB;                       // time index of first pixel in the layer 
   }                                                                                        
                                                                                            
   int ib=1;                                                                                
   int ie=I;                                                                                
   for(i=0; i<I; i++) {                           // select bandwidth                       
      if(pTF->frequency(i) <= pTF->gethigh()) ie=i;                                         
      if(pTF->frequency(i) <= pTF->getlow())  ib=i+1;                                       
   }                                                                                        
   if(ie>I-1) ie = I-1;                           // required by catalog                    
   if(ib<1)   ib = 1;                             // required by catalog                    

   slice S = pTF->getSlice(0);
   jE = S.size()-jB;                              // last good sample in the layer
   NN = jE-jB;                                    // #of good samples in the layer
   if(jE&1) {
      cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(3): WDM parity violation\n"; 
      gSystem->Exit(1);
   }       

   //cout<<r<<" "<<R<<" "<<I<<" "<<jB<<" "<<net->veto.size()<<endl;
   //cout<<ib<<" "<<ie<<" "<<NN<<" "<<jB<<" "<<jE<<endl;            

   for(jj=0; jj<NN; jj++) {                       // loop over time stamps

      double VETO = 1.;
      pmap = MAP.data+(jj+jB)*I;                  // pointer to 0 F sample in MAP
      for(n=0; n<nIFO; n++) {                                                    
         if(in[n] >= jE) in[n] -= NN;             // go to jB sample             
         jb = int(in[n]*r/R+0.01);                // first veto index            
         je = int((in[n]+1)*r/R+0.01);            // last veto index             
         while(jb<je) if(!pveto[jb++]) VETO=0.;   // set veto value              
         PDATA = &(pdata[n][in[n]*I]);            // pointer to 0 F sample       
         for(i=0; i<I; i++) pmap[i]+=*PDATA++;    // sum energy                  
         in[n]++;                                 // increment index pointer     
      }                                                                          
                                                                                 
      for(i=0; i<I; i++) {                                                       
         pmap[i] *= VETO;                                                        
      }                                                                          
      count += VETO;                              // count live time             
   }                                                                             

// set metadata in wc_List
   net->wc_List[LAG].start = pTF->start();  
   net->wc_List[LAG].stop  = pTF->stop();
   net->wc_List[LAG].rate  = pTF->rate();
   net->livTime[LAG] = count/R;                   // live time depends on resolution
}

TFile* GetRootFile(network* NET) {

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
      cout << "CWB_Plugin_netEvent.C : Error - output root file not found" << endl;
      gSystem->Exit(1);                                                                             
    }                                                                                      
  } else {                                                                                 
    cout << "CWB_Plugin_netEvent.C : Error - output root file not found" << endl;  
    gSystem->Exit(1);                                                                               
  }                                                                                        

  return froot;
}


void PlotWDM(WSeries<double>* WS, wavearray<double>* t, wavearray<double>* f) {

  // Plot WDM Scalogram
  watplot WTS(const_cast<char*>("WTS"));

  int layers = WS->maxLayer()+1;  // numbers of frequency bins (first & last bins have df/2)
  int slices = WS->sizeZero();    // number of time bins                                    

  float df = WS->resolution();    // frequency bin resolution (hz)
  float dt = 1./(2*df);           // time bin resolution (sec)    

  //scalogram maps

  double start = WS->start();
  double stop  = WS->start()+slices*dt;
  double flow  = 0;                    
  double fhigh = (layers-1)*df;        
  cout.precision(14);                  
  cout << "start " << start << " stop " << stop << " flow " << flow << " fhigh " << fhigh << endl;
  WTS.plot(*WS, 2, start, stop,const_cast<char*>("COLZ"));                                        
  // set frequency range                                                                          
  WTS.hist2D->GetYaxis()->SetRangeUser(flow, fhigh);                                              

  // draw white dashed polyline (example)
  TPolyLine cluster_lines(t->size(),t->data,f->data);

  cluster_lines.SetLineColor(kWhite);
  cluster_lines.SetLineStyle(2);
  cluster_lines.Draw();

  TPolyMarker cluster_points(t->size(),t->data,f->data);
  cluster_points.SetMarkerSize(2);
  cluster_points.SetMarkerStyle(kDot);
  cluster_points.SetMarkerColor(kWhite);
  cluster_points.Draw();

  // write plots and data to file
  char fname[1024];
  sprintf(fname,"wavegraph_demo_plots.root");
  WTS.canvas->cd();
  WTS.canvas->Print(fname);
}

