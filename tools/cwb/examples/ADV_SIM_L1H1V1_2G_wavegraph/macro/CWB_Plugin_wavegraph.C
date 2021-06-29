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

#define WAVEGRAPH 	// uncomment to enable the wavegraph analysis
#define WAVEGRAPH_FILE  "config/wavegraph_cbc_1d4x1d4_to_25x25.txt"
//#define PLOT_WAVEPATH   // plot wavepath on top of the TF maps
#define WAVEGRAPH_THR 400

//#define TEST_SHIFTED_DATA	// noise+mdc are produced in simulation=0 

// wavecraft coherence
void wcoherence(TFile* jfile, CWB::config* cfg, network* net);	
void getNetworkEnergyMap(int LAG, network* net, WSeries<double>& MAP, int& nM, int* IN, int& jB, int& jE);   
void PlotWDM(WSeries<double>* WS, wavearray<double>* t, wavearray<double>* f);

// standard coherence
void scoherence(TFile* jfile, CWB::config* cfg, network* net);	
long getNetworkPixels(int LAG, double Eo, double Em, network* net);

WSeries<double> edMAP[NIFO_MAX][NRES_MAX]; 	// detector energy TF maps

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// This plugin implements the wavegraph in coherence stage (only 2G)

  cout << endl;
  cout << "-----> CWB_Plugin_wavegraph.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
#ifdef TEST_SHIFTED_DATA
    cfg->dataPlugin=true; 	// disable read data from frames
#else
    cfg->mdcPlugin=true;  	// disable read mdc from frames   
#endif
    cfg->cohPlugin=true;  	// disable built-in coherence stage
    //cfg->scPlugin=true;  	// disable built-in supercluster function
  }

#ifdef TEST_SHIFTED_DATA

  if(type==CWB_PLUGIN_DATA) {  

    CWB::Toolbox TB;

    int seed;
    if(ifo.CompareTo("L1")==0) seed=1000;
    if(ifo.CompareTo("H1")==0) seed=2000;
    if(ifo.CompareTo("V1")==0) seed=3000;
    if(ifo.CompareTo("J1")==0) seed=4000;

    TString fName;
    if(ifo.CompareTo("L1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("V1")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";
    if(ifo.CompareTo("J1")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";            

    int size=x->size();
    double start=x->start();
    TB.getSimNoise(*x, fName, seed, net->nRun);
    x->resize(size);                           
    x->start(start);                           

    wavearray<double> y = *x;
    wavearray<double> z = y; 

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
                                        
    IMPORT(CWB::mdc,MDC)                
    MDC.Print();                        

    // ---------------------------------
    // get mdc data                     
    // ---------------------------------

    MDC.Get(y,ifo);

    int binLen   = (cfg->segLen+2*cfg->segEdge)  * y.rate();
    int binEdge  = cfg->segEdge * y.rate();                 

    int segShift=0;
    if(ifo=="L1") segShift=0;
    if(ifo=="H1") segShift=50;
    if(ifo=="V1") segShift=20; 
    int binShift = segShift * y.rate();

    cout << ifo << " " << binEdge+binShift << " " << binLen-binEdge << endl;
    cout << ifo << " " << binEdge << " " << binEdge+binShift << endl;

    z=0;
    int j=binEdge;
    for(int i=binLen-binEdge-binShift;i<binLen-binEdge;i++) z[j++]=y[i];
    for(int i=binEdge;i<binLen-binEdge-binShift;i++) z[j++]=y[i];
    y=z;

    // ---------------------------------
    // set mdc list in the network class
    // ---------------------------------

    cout.precision(14);
    for(int k=0;k<(int)MDC.mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)MDC.mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)MDC.mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;

    y *= cfg->factors[0];	// normalization

    *x+=y;
  }

#else

  if(type==CWB_PLUGIN_MDC) {  

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
                                        
    IMPORT(CWB::mdc,MDC)                
    MDC.Print();                        

    // ---------------------------------
    // get mdc data                     
    // ---------------------------------

    MDC.Get(*x,ifo);

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

    cout.precision(14);
    for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;
  }                                                                                                

#endif

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
      }
    } 

#ifdef WAVEGRAPH
    wcoherence(jfile, cfg, net);	// wavegraph coherence	
#else
    scoherence(jfile, cfg, net);	// standard coherence
#endif
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

  int nIFO = net->ifoListSize();		// number of detectors
  int nRES = net->wdmListSize();		// number of resolution levels

  netcluster* pwc;
  netpixel pix(nIFO); 

  // network energy maps (sum over the detectors)
  WSeries<double> enMAP[NRES_MAX]; 	
  // store energy TF maps pointers to data vector
  vector< WSeries<double>* > data;
  for(int i=0; i<nRES; i++) data.push_back(&enMAP[i]);

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
       }
       int ii = nRES-i-1;
       getNetworkEnergyMap(j,net,enMAP[ii],nM[ii],IN[ii],jB[ii],jE[ii]);
     }

     // apply wavegraph
     std::vector<cluster> clusters = graph.clustering(WAVEGRAPH_THR,data,cfg->segEdge);
     //cout << "Number of paths found by wavegraph : " << clusters.size() << endl;

#ifdef PLOT_WAVEPATH
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
     PlotWDM(data.at(4), &t, &f);
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
       pwc->write(jfile,"coherence","clusters",-1,cycle);        
       nclusters+=pwc->csize();
       npixels+=pwc->size();
       csize_tot+=pwc->csize(); psize_tot+=pwc->size();                   

       pwc->clear();				 // clean pixels
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

//**************************************************************************
// standard coherence stage
//**************************************************************************
void scoherence(TFile* jfile, CWB::config* cfg, network* net) {

  // import ifactor
  int gIFACTOR=-1; IMPORT(int,gIFACTOR)
  cout << "-----> CWB_Plugin_wavegraph.C -> " << " gIFACTOR : " << gIFACTOR << endl;

  int nIFO = net->ifoListSize();
  int nRES = net->wdmListSize();			 // number of resolution levels

  double Eo;
  netcluster* pwc;
 
  for(int i=0; i<nRES; i++) {                            // loop over TF resolutions

    for(int n=0; n<nIFO; n++) {                    
      // restore edMAP energy map 
      *(net->getifo(n)->getTFmap()) = edMAP[n][i]; 
    }

    Eo = net->THRESHOLD(cfg->bpp);                        // threshold on pixel energy
    cout<<"thresholds in units of noise variance: Eo="<<Eo<<" Emax="<<Eo*2<<endl;   

    double TL = net->setVeto(cfg->iwindow);
    cout<<"live time in zero lag: "<<TL<<endl;          		// set veto array
    if(TL <= 0.) {cout<<"livetime is zero : exit"<<endl;gSystem->Exit(1);}	// exit if live time is zero
  
    // select pixels
    if(cfg->simulation) {cout<<"ifactor|clusters|pixels ";cout.flush();}
    else                {cout<<"lag|clusters|pixels ";    cout.flush();}
    int csize_tot=0;int psize_tot=0;                                   
    for(int j=0; j<(int)net->nLag; j++) {                               

       //net->getNetworkPixels(j,Eo,Eo*2);
       getNetworkPixels(j,Eo,Eo*2,net);
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

//**************************************************************************
//:select TF samples by value of the network excess energy: 2-8 detectors   
//**************************************************************************
long getNetworkPixels(int LAG, double Eo, double Em, network* net)   
{                                                                           
// 2G analysis algorithm for selection of significant network pixles        
// works with WDM/wavelet energy TF maps                                    
// LAG - time shift lag defining how detectors are shifted wrt each other.  
// net - network pointer                                 

   size_t nIFO = net->ifoList.size();       		// number of detectors

   if(nIFO>NIFO) {
      cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(): " 
          <<"invalid number of detectors or\n";
      return 0;                                
   }                                           
   if(net->getifo(0)->getTFmap()->w_mode != 1) {    
      cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(): invalid whitening mode.\n"; 
      return 0;                                                       
   }                                                                  

   WSeries<double>* pTF = net->getifo(0)->getTFmap();  // pointer to first TF map
   WSeries<double> MAP; MAP = *pTF; MAP=0.;            // initialize TF map      
   wavearray<double>* hTS = net->getifo(0)->getHoT();  // pointer to first TS data
                                                                                  
   int i,j,k,m,n,NN,jj,nM,jE,jb,je,J,K;                                           

   double Eh = Em*Em;                                  // halo energy^2           
   double R  = pTF->wrate();                           // pixel layer rate        
   double r  = hTS->rate();                            // TS rate                 
   int N  = pTF->size();                               // size of TF array        
   int M  = hTS->size();                               // size of TS array        
   int I  = pTF->maxLayer()+1;                         // number of layers        
   int II = pTF->maxLayer()-1;                         // number of layers - 2    
   int jB = int(net->Edge*R+0.001);                    // number of samples in the edges
   if(jB&1) {cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(1): WDM parity violation\n"; exit(1);}             
                                                                                        
   if(jB < 3) {                                                                         
      cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(): insufficient data edge length.\n";            
      exit(1);                                                                          
   }                                                                                    

   netpixel pix(nIFO); 
   pix.core = true;    
   pix.rate = R;       
   pix.layers = I;     
                       
   int     in[NIFO];                                    // pixel time index
   int     IN[NIFO];                                    // pixel time index
   double* PDATA;                                                          
   double* pmap;                                                           
   double* pdata[NIFO];                                 // pointers to data
   double* pp[5];                                       // pointers to sorted F-arrays
   for(n=0; n<nIFO; n++) {                              // pointers to data           
      pdata[n] = net->getifo(n)->getTFmap()->data;                                         
   }                                                                                  

   long nPix = 0;
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
      if(K&1) {cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(2): WDM parity violation\n"; exit(1);}
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
   if(jE&1) {cout<<"CWB_Plugin_wavegraph.C : getNetworkPixels(3): WDM parity violation\n"; exit(1);}       

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
         if(pmap[i]<Eo || i<ib) pmap[i]=0.;        // zero sub-threshold pixels   
         if(pmap[i]>Em) pmap[i]=Em+0.1;            // degrade loud pixels         
      }                                                                          
      count += VETO;                               // count live time             
   }                                                                             

   for(jj=0; jj<NN; jj++) {                        // loop over time stamps

      pmap = MAP.data+(jj+jB)*I;                   // pointer to 0 F sample in MAP
      for(n=0; n<nIFO; n++) {                                                     
         if(IN[n] >= jE) IN[n] -= NN;              // go to jB sample             
      }                                                                           
      for(n=0; n<5; n++) pp[n]=pmap+(n-2)*I;       // initialize set of pointers  
      for(i=ib; i<ie; i++) {                                                      
         if((E=pp[2][i])<Eo) continue;             // skip subthreshold pixels    
         Ct = pp[2][i+1]+pp[3][ i ]+pp[3][i+1];    // top core                    
         Cb = pp[2][i-1]+pp[1][ i ]+pp[1][i-1];    // bottom core                 
         Ht = pp[4][i+1];                          // top halo                    
         Ht+= i<II? pp[4][i+2]+pp[3][i+2] : 0.;    // top halo                    
         Hb = pp[0][i-1];                          // bottom halo                 
         Hb+= i>1 ? pp[0][i-2]+pp[1][i-2] : 0.;    // bottom halo                 

         if((Ct+Cb)*E<Eh && 
            (Ct+Ht)*E<Eh && 
            (Cb+Hb)*E<Eh && 
            E<Em) continue; 
                            
         E = 0;             
         for(n=0; n<nIFO; n++) {
            j = IN[n]*I+i;                         // sample index
cout << n << " DEB1 " << " I " << I << " IN[n] " << IN[n] << " i " << i << " j " << j << endl;
            pix.data[n].index = j;                                 
            pix.data[n].asnr = pdata[n][j];                        
            E += pdata[n][j];                                      
         }                                                         
         j = IN[nM]*I+i;                           // reference sample index
         pix.time = j;                                                       
         pix.frequency = i;                                                  
cout << "DEB2 " << " I " << I << " IN[nM] " << IN[nM] << " i " << i << " j " << j << endl;
//cout << "DEB " << "pix.time " << pix.time << " pix.frequency " << pix.frequency << " I " << I << " IN[nM] " << IN[nM] << " i " << i << endl;
         pix.likelihood = E;                                                 
         net->wc_List[LAG].append(pix);            // save pixels in wc_List
         nPix++;                                                             
      }                                                                      
      for(n=0; n<nIFO; n++) IN[n]++;               // increment IN          
   }                                                                         

// set metadata in wc_List
   net->wc_List[LAG].start = pTF->start();  
   net->wc_List[LAG].stop  = pTF->stop();
   net->wc_List[LAG].rate  = pTF->rate();
   net->livTime[LAG] = count/R;                    // live time depends on resolution

   if(nPix) net->setRMS();

   return nPix;
}

