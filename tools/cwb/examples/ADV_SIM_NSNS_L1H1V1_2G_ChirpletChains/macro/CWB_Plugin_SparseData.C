#define XIFO 4
#define _USE_HEALPIX

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
#include "TProfile.h"
#include "mdc.hh"
#include "watplot.hh"
#include <vector>

// Plugin to Test Chirplet Chains

#define GET_PIXELS

//#define DUMP_WDM_TF_IFO
//#define DUMP_WDM_TF_NET

#define MIN_SKYRES_HEALPIX      4
#define MIN_SKYRES_ANGLE        3

void PlotData(WSeries<double>* WS, TString tag, TString odir);
vector<bool> SelectUnique(network* net);
void GetPixels(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString scycle, vector<SSeries<double> > &vSS);
void PlotMRData(vector<SSeries<double> > &vSS, TString odir);
Double_t fchirp(Double_t* x, Double_t* par);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString scycle, int type)  {

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_SparseData.C" << endl;
  cout << "cycle " << scycle.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_ILIKELIHOOD) {  
     cout << "type==CWB_PLUGIN_ILIKELIHOOD" << endl;
     int nIFO = net->ifoListSize();
     detector* pD[NIFO_MAX];                       //! pointers to detectors
     for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);
     for(int n=0;n<nIFO;n++) cout << n << " detector Name : " << pD[n]->Name << endl;

#ifdef GET_PIXELS

     int nSS = pD[0]->vSS.size();
     SSeries<double> ss;
     vector<SSeries<double> > vSS;
     for(int i=0; i<nSS; i++) vSS.push_back(ss);	          // fill vector
     for(int i=0;i<nSS;i++) vSS[i] = pD[0]->vSS[i];               // copy sparse map
     for(int i=0;i<(int)vSS.size();i++) vSS[i].Expand(true);      // expand sparse map (only core pixels)

     GetPixels(jfile, cfg, net, x, scycle, vSS);

     PlotMRData(vSS, cfg->dump_dir);				// plot event energy integrated over the resolutions	

/*
     // plot data at different resolutions
     for(int i=0;i<nSS;i++) {			// loop over the resolutions
       WSeries<double>* WS = (WSeries<double>*)(&vSS[i]); 
       PlotData(WS, "EE",cfg->dump_dir); 
       vSS[i].Shrink();				// resize 0 wseries : leave only sparse table 
     }
*/
#endif 

#ifdef DUMP_WDM_TF_IFO
     for(int n=0;n<nIFO;n++) {
       int nSS = pD[n]->vSS.size();
       cout << n << " -> " << nSS << endl;
       for(int i=0;i<nSS;i++) {

         pD[n]->vSS[i].Expand(true);	// rebuild wseries from sparse table 
         cout << "REBUILD SIZE " << pD[n]->vSS[i].size() << endl; 

         WSeries<double>* WS = (WSeries<double>*)(&(pD[n]->vSS[i])); 
         PlotData(WS, pD[n]->Name,cfg->dump_dir); 

//         pD[n]->vSS[i].Shrink();	// resize 0 wseries : leave only sparse table 
       } 
     }
#endif

#ifdef DUMP_WDM_TF_NET
     int nSS = pD[0]->vSS.size();
     for(int i=0;i<nSS;i++) {
       SSeries<double> SS = pD[0]->vSS[i];
       for(int j=0;j<(int)SS.sparseIndex.size();j++) {
         SS.sparseMap00[j]=0;
         SS.sparseMap90[j]=0;
         if(!SS.sparseType.TestBitNumber(j)) continue;
         for(int n=0;n<nIFO;n++) {
           SS.sparseMap00[j]+=pow(pD[n]->vSS[i].sparseMap00[j],2);
           SS.sparseMap90[j]+=pow(pD[n]->vSS[i].sparseMap90[j],2);
         }
         SS.sparseMap00[j]=sqrt(SS.sparseMap00[j]);
         SS.sparseMap90[j]=sqrt(SS.sparseMap90[j]);
       }

       SS.Expand(true);	// rebuild wseries from sparse table 
       cout << "REBUILD SIZE " << SS.size() << endl; 

       WSeries<double>* WS = (WSeries<double>*)(&SS); 
       char strnet[1024]="";
       for(int n=0;n<nIFO;n++) sprintf(strnet,"%s%s",strnet, pD[n]->Name);
       PlotData(WS, strnet,cfg->dump_dir); 

//       pD[n]->vSS[i].Shrink();	// resize 0 wseries : leave only sparse table 
     } 
#endif
  }

  return;
}

void PlotData(WSeries<double>* WS, TString tag, TString odir) {

  // Plot WDM Scalogram
  watplot WTS(const_cast<char*>("WTS"));

  //scalogram maps
  double start = WS->start();
  double stop  = WS->start()+WS->size()/WS->rate();
  double flow  = 0;
  double fhigh = 2048;
  WTS.plot(*WS, 2, start, stop,const_cast<char*>("COLZ"));
  WTS.hist2D->GetYaxis()->SetRangeUser(flow, fhigh);

  // fill WTS.hist2D with nearest pixels
  int xsize=WTS.hist2D->GetNbinsX();
  int ysize=WTS.hist2D->GetNbinsY();
  cout << "xsize : " << xsize << endl;
  cout << "ysize : " << ysize << endl;

  // dump spectrum
  char fname[1024];
  sprintf(fname,"%s/%s_wdm_scalogram_%d_lev%d.root", odir.Data(), tag.Data(), int(WS->start()),WS->getLevel());
  cout << endl << "Dump WDM Scalogram : " << fname << endl << endl;
  WTS.canvas->cd();
  WTS.canvas->Print(fname);

  return;
}

void PlotMRData(vector<SSeries<double> > &vSS, TString odir) {

  TCanvas* canvas;

  canvas= new TCanvas("SS", "SS", 200, 20, 800, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);
  canvas->SetRightMargin(0.10);
  canvas->SetLeftMargin(0.10);
  canvas->SetBottomMargin(0.13);
  canvas->SetBorderMode(0);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);


  TH2F* hist2D;

  int nSS = vSS.size();
  for(int n=0;n<nSS;n++) {                      // loop over the resolutions
    int layers = vSS[n].GetLayers();            // numbers of frequency bins (first & last bins have df/2)
    int slices = vSS[n].GetSlices();            // number of time bins
    float dt = vSS[n].GetTimeResolution();      // time bin resolution (sec)
    float df = vSS[n].GetFreqResolution();      // frequency bin resolution (Hz)
    int rate = int(1./dt);
    cout << n << " ->\t layers : " << layers << "\t slices : " << slices << "\t rate : " << rate
              << "\t dt : " << dt << "\t df : " << df << endl;
  }

  float dt = vSS[0].GetTimeResolution();
  float df = vSS[nSS-1].GetFreqResolution()/2;

  int nT = vSS[0].GetSlices();
  int nF = 2*(vSS[nSS-1].GetLayers()-1);

  cout << endl;
  cout << "  nT=" << nT << "  nF=" << nF 
       << "  dt=" << dt << "  df=" << df  
       << "  nT*dt=" << nT*dt << "  nF*df=" << nF*df << endl;

  hist2D=new TH2F("vSS","", nT, 0, nT*dt, nF, 0., nF*df);

  for(int n=0;n<nSS;n++) {                      // loop over the resolutions
    int slices = vSS[n].GetSlices();            // number of time bins
    int layers = vSS[n].GetLayers();            // numbers of frequency bins (first & last bins have df/2)

    int mT = nT/slices;
    int mF = nF/(2*(layers-1));

    for(int i=0;i<slices;i++) {
      for(int j=0;j<layers;j++) {
        int S = (j==0&&j==nF-1) ? 1 : 2;
//if(layers!=9) continue;
        float ee = vSS[n].GetMap00(i,j);
        ee/=(mT*S*mF);  // normalize energy
        for(int p=0;p<mT;p++) {
          for(int q=0;q<S*mF;q++) {
            int slice = mT*i+p+1;
            int layer = (j==0&&j==nF-1) ? mF*j+q+1 : S*mF*(j-1)+q+1;
            double EE = hist2D->GetBinContent(slice,layer);
            hist2D->SetBinContent(slice,layer,ee+EE);
          }
        }
      }
    }
  }

  hist2D->SetXTitle("time, sec");
  hist2D->SetYTitle("frequency, Hz");

  hist2D->SetStats(kFALSE);
  hist2D->SetTitleFont(12);
  hist2D->SetFillColor(kWhite);

  char title[256];
  sprintf(title,"Scalogram (%s)","energy");
  hist2D->SetTitle(title);

  hist2D->GetXaxis()->SetNdivisions(506);
  hist2D->GetXaxis()->SetLabelFont(42);
  hist2D->GetXaxis()->SetLabelOffset(0.014);
  hist2D->GetXaxis()->SetTitleOffset(1.4);
  hist2D->GetYaxis()->SetTitleOffset(1.2);
  hist2D->GetYaxis()->SetNdivisions(506);
  hist2D->GetYaxis()->SetLabelFont(42);
  hist2D->GetYaxis()->SetLabelOffset(0.01);
  hist2D->GetZaxis()->SetLabelFont(42);
  hist2D->GetZaxis()->SetNoExponent(false);
  hist2D->GetZaxis()->SetNdivisions(506);

  hist2D->GetXaxis()->SetTitleFont(42);
  hist2D->GetXaxis()->SetTitle("Time (sec)");
  hist2D->GetXaxis()->CenterTitle(true);
  hist2D->GetYaxis()->SetTitleFont(42);
  hist2D->GetYaxis()->SetTitle("Frequency (Hz)");
  hist2D->GetYaxis()->CenterTitle(true);

  hist2D->GetZaxis()->SetTitleOffset(0.6);
  hist2D->GetZaxis()->SetTitleFont(42);
  //hist2D->GetZaxis()->SetTitle(ztitle);
  hist2D->GetZaxis()->CenterTitle(true);

  hist2D->GetXaxis()->SetLabelSize(0.03);
  hist2D->GetYaxis()->SetLabelSize(0.03);
  hist2D->GetZaxis()->SetLabelSize(0.03);

  canvas->cd();

  hist2D->Draw("COLZ");

/*
  // draw profile
  TProfile *prof = hist2D->ProfileX(" ");
  prof->SetErrorOption("g");  // s,i,g
  prof->SetLineColor(kBlack);
  prof->Draw("LP");
  //prof->Draw("same");

  // fit chirp
  TF1* chirp = new TF1("fchirp",fchirp,0, nT*dt,3);
  chirp->SetLineColor(kRed);
  chirp->SetLineWidth(1);
  //chirp->SetParameter(1,xmax);
  chirp->SetParameter(2,3./8.);
  chirp->SetParLimits(2,3./8.,3./8.);

  double chi2;
  double par[3];
  double epar[3];
  prof->Fit("fchirp","+rob=0.75");
  chi2 = chirp->GetChisquare();
  chirp->GetParameters(par);
  for (int k=0;k<3;k++) epar[k] = chirp->GetParError(k);
  cout << "chi2 " << chi2 << endl;
  for (int k=0;k<3;k++) cout << "Par  " << k << " : " << par[k] << endl;
  for (int k=0;k<3;k++) cout << "ePar " << k << " : " << epar[k] << endl;
  int ndf = chirp->GetNDF();
  double plevel = TMath::Prob(chi2,ndf);
  cout << "ndf " << ndf << endl;
*/ 

  // dump spectrum
  char fname[1024];
  sprintf(fname,"%s/SS_wdm_scalogram.root", odir.Data());
  cout << endl << "Dump SS Scalogram : " << fname << endl << endl;
  canvas->Print(fname);

  delete hist2D;
  delete canvas;
}

vector<bool> SelectUnique(network* net) {

  int nIFO = net->ifoListSize();

  vector<int> v[NIFO];
  short*  ml[NIFO];

  // init delay index array 
  for(int i=0; i<NIFO; i++) {
    if(i<nIFO) ml[i] = net->getifo(i)->index.data;
    else       ml[i] = net->getifo(0)->index.data;
  }

  int N = net->getifo(0)->index.size();
  vector<bool> skip(N);

  for(int n=0;n<N;n++) {
    skip[n]=false;
    int equals=0;
    for(int m=0;m<v[0].size();m++) {
      equals=0;
      for(int i=0;i<nIFO;i++) if(ml[i][n]==v[i][m]) equals++;
      if(equals==nIFO) break;
    }
    if(equals!=nIFO) {
      for(int i=0;i<nIFO;i++) v[i].push_back(ml[i][n]);
    } else skip[n]=true;
  }
  cout << "Unique sky positions : " << v[0].size() << "/" << N << endl;
  return skip;
}

void
GetPixels(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString scycle, vector<SSeries<double> > &vSS)  {

  std::vector<int>* vint;
  std::vector<int>* vtof;
  wavearray<double> cid;                   // get cluster ID
  netpixel* pix;                              
  double TDRate;                              

  int nIFO = net->ifoListSize();
  netcluster* pwc = net->getwc(scycle.Atoi());

  // decrease skymap resolution to improve performances
  double skyres=0;                                     
  if(cfg->healpix) skyres = cfg->healpix>MIN_SKYRES_HEALPIX ? MIN_SKYRES_HEALPIX : 0;
  else             skyres = cfg->angle<MIN_SKYRES_ANGLE ? MIN_SKYRES_ANGLE : 0;      
  if(skyres) {                                                                       
    if(cfg->healpix) net->setSkyMaps(int(skyres));                                   
    else             net->setSkyMaps(skyres,cfg->Theta1,cfg->Theta2,cfg->Phi1,cfg->Phi2);
    net->setAntenna();                                                                   
    net->setDelay(cfg->refIFO);                                                          
  }                                                                                      
  TDRate = (cfg->inRate>>cfg->levelR)*1;  // time-delay filter rate                      
  cout << "TDRate : " << TDRate << endl;                                                 
  net->setDelayIndex(TDRate);                                                            

  // init delay index array 
  short*  ml[NIFO];         
  for(int i=0; i<NIFO; i++) {
    if(i<nIFO) ml[i] = net->getifo(i)->index.data;
    else       ml[i] = net->getifo(0)->index.data;
  }                                               

  int L = net->skyMask.size();             // sky grid size                   

  //vector<bool> skip(net->getifo(0)->index.size());
  //for(int i=0;i<skip.size();i++) skip[i]=false;   
  vector<bool> skip = SelectUnique(net);        

  // build lookup table  ;  rate -> vSS resolution index                                                     
  size_t rateANA=cfg->inRate>>cfg->levelR;                                                            
  int* luTable = new int[rateANA];                                                                    
  for(int i=0;i<(int)vSS.size();i++) {          // loop over the resolutions
    int layers = vSS[i].GetLayers();		// numbers of frequency bins (first & last bins have df/2)
    int slices = vSS[i].GetSlices();            // number of time bins
    float dt = vSS[i].GetTimeResolution();      // time bin resolution (sec) 
    float df = vSS[i].GetFreqResolution();      // frequency bin resolution (Hz) 
    int rate = int(1./dt);                                                                      
    cout << i << " ->\t layers : " << layers << "\t slices : " << slices << "\t rate : " << rate
              << "\t dt : " << dt << "\t df : " << df << endl;
    luTable[rate]=i;                                                                                
  }                                                                                                   

  // read cluster list & metadata netcluster object
  TString cldir = (jfile->Get("likelihood;1")!=NULL) ? "likelihood" : "supercluster";
  int cycle = cfg->simulation ? scycle.Atoi() : Long_t(net->wc_List[scycle.Atoi()].shift);
  vector<int> clist = pwc->read(jfile,cldir.Data(),"clusters",0,cycle);                   
  pwc->print();                                                                           

  // read pixels & tdAmp into netcluster pwc
  for(int m=0;m<(int)clist.size();m++) {                      // loop over the cluster list
    int nmax = -1;                                            // load all tdAmp                                
    //int nmax = 2000;        // read no more than maxPix pixels from a cluster
    pwc->read(jfile,cldir.Data(),"clusters",nmax,cycle,0,clist[m]);            
    cid = pwc->get((char*)"ID",  0,'S',0);                    // get cluster ID                 
    pwc->setcore(false,size_t(cid.data[cid.size()-1]+0.1));                    
    pwc->loadTDampSSE(*net, 'a', cfg->BATCH, cfg->BATCH);     // attach TD amp to pixels

    int K = cid.size();
    cout << "Number of Clusters " << K << endl;
    for(int k=0; k<K; k++) {                                  // loop over clusters

      int id = size_t(cid.data[k]+0.1);
      cout << k << " Cluster ID " << id << endl;

      if(pwc->sCuts[id-1] != -2) continue;                    // skip rejected/processed clusters

      vint = &(pwc->cList[id-1]);                             // pixel list

      int V = vint->size();
      cout << k << " Cluster Size : " << V << endl;
      if(!V) continue;                             

      for(int j=0; j<V; j++) {                                // loop over pixels
        pix = pwc->getPixel(id,j);                                           
        if(!pix) {                                                           
          cout<<"Plugin error: NULL pointer"<<endl;                          
          gSystem->Exit(1);                                                  
        }                                                                    

        int tsize = pix->tdAmp[0].size();
        tsize/=2;                        

        int rate = int(pix->rate+0.5);
//        skymap sm(int(MIN_SKYRES_HEALPIX));                                     
        double maxEE=0;                              // max energy over the sky loop
        for(int l=0;l<L;l++) {                       // sky loop 
          if(skip[l]) continue;                      // skip duplicate ifo delays                             
          double EE=0;                                                            
          for(int n=0; n<nIFO; n++) {                                             
            int q=tsize/2+ml[n][l];                                               
            double aa = pix->tdAmp[n].data[q];       // copy TD 00 data           
            double AA = pix->tdAmp[n].data[q+tsize]; // copy TD 90 data           
            EE+=(aa*aa+AA*AA)/2;                                                  
          }                                                                                                     
//          if(j==30) sm.set(l,EE);                                                                             
          if(EE>maxEE) maxEE=EE;                                                                                
        }                                                                                                       
        //cout << luTable[rate] << " " << id << " " << j << " " << maxEE << endl;                                 
        char fitsName[256];sprintf(fitsName,"EE.fits");                                                         
//        if(j==30) {sm.Dump2fits(fitsName);gSystem->Exit(0);}                                                  
        int index = pix->data[0].index;                                                                         
        vSS[luTable[rate]].SetMap00(index,maxEE);                                                        
        vSS[luTable[rate]].SetMap90(index,maxEE);                                                            
      }                                                                                                         
      pwc->sCuts[id-1]=1;                                                                                       
    }                                                                                                           
  }                                                                                                             
  pwc->clear();                                      // unload clusters from network class 

  // restore skymap resolution
  if(skyres) {                
    if(cfg->healpix) net->setSkyMaps(int(cfg->healpix));
    else            net->setSkyMaps(cfg->angle,cfg->Theta1,cfg->Theta2,cfg->Phi1,cfg->Phi2);
    net->setAntenna();                                                                      
    net->setDelay(cfg->refIFO);                                                             
    if(cfg->mask>0.) net->setSkyMask(cfg->mask,cfg->skyMaskFile);                           
    if(strlen(cfg->skyMaskCCFile)>0) net->setSkyMaskCC(cfg->skyMaskCCFile);                 
  }                                                                                         
  TDRate = (cfg->inRate>>cfg->levelR)*cfg->upTDF;   // time-delay filter rate                
  net->setDelayIndex(TDRate);                                                               

  delete [] luTable;                                                                    
}                                                             

Double_t fchirp(Double_t* x, Double_t* par) {

  if(par[1]-x[0]<0) {TF1::RejectPoint();return 0;}

  double fitval = par[0]*pow(par[1]-x[0],-par[2]);
  return fitval;
}

