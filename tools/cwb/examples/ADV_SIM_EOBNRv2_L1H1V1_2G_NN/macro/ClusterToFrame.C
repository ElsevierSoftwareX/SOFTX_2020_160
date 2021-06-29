//
// Convert netcluster principal components to NN (Neural Network) frames (NDIMxNDIM pixels)
// Author : Gabriele Vedovato
// Note : run with ACLiC -> root -l 'ClusterToFrame.C+("ifile.root","ofile.root")'
//

#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPaletteAxis.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TMath.h"
#include "mdc.hh"
#include "watplot.hh"
#include <vector>

TH2F* hist2D=NULL;
TCanvas* canvas;

void PlotFrame(std::vector<double>* nn_frame, int cid);
void ConvertToFrame(netcluster* pwc, int cid, int nifo, char type, int irate, std::vector<double>* nn_frame);

#define NDIM 6

//#define WATPLOT         // if defined the event monster plots are produced
//#define PLOT_FRAME      // frame display

void ClusterToFrame(TString ifName, TString ofName) {

#ifdef WATPLOT
    watplot WTS(const_cast<char*>("wts"));
#endif

   // open input root file
   TFile* ifile = new TFile(ifName);
   if(ifile==NULL) {
     cout << "ClusterToFrame - Error : error opening file : " << ifName.Data() << endl;
     gSystem->Exit(1);
   } 

   // get waveburst tree
   TTree* nn_itree = (TTree *) ifile->Get("waveburst");
   if(nn_itree==NULL) {
     cout << "ClusterToFrame - Error : tree waveburst not found !!!" << endl;
     gSystem->Exit(1);
   } 

   netcluster* pwc = new netcluster();
   nn_itree->SetBranchAddress("netcluster",&pwc);

   // get number of detectors
   int nIFO = nn_itree->GetUserInfo()->GetSize();
   if(nIFO==0) {
     cout<<"ClusterToFrame - Error : wrong user detector list in header tree"<<endl; 
     gSystem->Exit(1);
   }

   // open output root file
   TFile ofile(ofName,"RECREATE");
   // create output tree 
   TTree* nn_otree = (TTree*)nn_itree->CloneTree(0);
   //nn_otree->SetName("nn");
   //nn_otree->SetTitle("nn");
   // add the NDIM*NDIM pixels vector for NN analysis
   std::vector<double>* nn_frame = new vector<double>;    
   nn_otree->Branch("nnframe", &nn_frame, NDIM*NDIM, 0);

   // read nn_itree
   int nn_isize = nn_itree->GetEntries();
   cout << "nn_itree size : " << nn_isize << endl;

   for(int m=0;m<nn_isize;m++) {

     if(m%100==0) cout << "ClusterToFrame -> " << m << "/" << nn_isize << endl; 

     nn_itree->GetEntry(m);

#ifdef WATPLOT                          // monster event display
     WTS.canvas->cd();
     char fname[1024];
     sprintf(fname, "monster_event_display_%d.png",m);
     cout << "write " << fname << endl;
     WTS.plot(pwc, 1, nIFO, 'L', 0, const_cast<char*>("COLZ")); // use ID=1
     WTS.canvas->Print(fname);
     WTS.clear();
#endif

     int ID = 1;			// ID=1 because each netcluster contains only 1 event
     ConvertToFrame(pwc, ID, nIFO, 'L', 0, nn_frame);
     //for(int i=0;i<nn_frame->size();i++) cout << i << " " << (*nn_frame)[i] << endl;
#ifdef PLOT_FRAME                       // frame display
     PlotFrame(nn_frame, m);
#endif

     nn_otree->Fill();
     nn_frame->clear();
   }

   nn_otree->Write();
   ofile.Close();

   ifile->Close();

   gSystem->Exit(0);
}

void                                                                        
ConvertToFrame(netcluster* pwc, int cid, int nifo, char type, int irate, std::vector<double>* nn_frame) {

  double FRAME[NDIM][NDIM];

  if(type != 'L' && type != 'N') return;

  double RATE = pwc->rate;                              // original rate

  std::vector<int>* vint = &(pwc->cList[cid-1]);        // pixel list

  int V = vint->size();                                 // cluster size
  if(!V) return;                                                       

  bool optrate=false;
  if(irate==1) {                                        // extract optimal rate
    optrate=true;                                                              
    std::vector<int>* pv = pwc->cRate.size() ? &(pwc->cRate[cid-1]) : NULL;          
    irate = pv!=NULL ? (*pv)[0] : 0;                                           
  }                                                                            

  int minLayers=1000;
  int maxLayers=0;   
  double minTime=1e20;
  double maxTime=0.;  
  double minFreq=1e20;
  double maxFreq=0.;  
  for(int j=0; j<V; j++) {                      // loop over the pixels
    netpixel* pix = pwc->getPixel(cid,j);                               
    if(!pix->core) continue;                                            

    if((irate)&&(irate != int(pix->rate+0.5))) continue;  // if irate!=0 skip rate!=irate

    if(pix->layers<minLayers) minLayers=pix->layers;
    if(pix->layers>maxLayers) maxLayers=pix->layers;

    double rate = RATE/(pix->layers-1); 
    double time=((double)pix->time/rate)/pix->layers;
    if(time<minTime) minTime=time;                   
    if(time>maxTime) maxTime=time;                   

    double freq = pix->frequency*rate/2.; 
    if(freq<minFreq) minFreq=freq;        
    if(freq>maxFreq) maxFreq=freq;        
  }                                       

  int minRate=RATE/(maxLayers-1);
  int maxRate=RATE/(minLayers-1);

  double mindt = 1./maxRate;
  double maxdt = 1./minRate;
  double mindf = minRate/2.;
  double maxdf = maxRate/2.;

  //cout << "minRate : " << minRate << "\t\t\t maxRate : " << maxRate << endl;
  //cout << "minTime : " << minTime << "\t\t\t maxTime : " << maxTime << endl;  
  //cout << "minFreq : " << minFreq << "\t\t\t maxFreq : " << maxFreq << endl;  
  //cout << "mindt   : " << mindt   << "\t\t\t maxdt   : " << maxdt << endl;  
  //cout << "mindf   : " << mindf   << "\t\t\t maxdf   : " << maxdf << endl;  

  int iminTime = int(minTime);
  int imaxTime = int(maxTime+1);
  int nTime = (imaxTime-iminTime)*maxRate;

  if(hist2D) { delete hist2D; hist2D=NULL; }
  hist2D=new TH2F("WTF", "WTF", nTime, iminTime, imaxTime, 2*maxLayers, 0, RATE/2);
  hist2D->SetXTitle("time, sec");                                                  
  hist2D->SetYTitle("frequency, Hz");                                              

  double dFreq = (maxFreq-minFreq)/10.>2*maxdf ? (maxFreq-minFreq)/10. : 2*maxdf ;
  double mFreq = minFreq-dFreq<0 ? 0 : minFreq-dFreq;                             
  double MFreq = maxFreq+dFreq>RATE/2 ? RATE/2 : maxFreq+dFreq;                   
  hist2D->GetYaxis()->SetRangeUser(mFreq, MFreq);                                 

  double dTime = (maxTime-minTime)/10.>2*maxdt ? (maxTime-minTime)/10. : 2*maxdt ;
  double mTime = minTime-dTime<iminTime ? iminTime : minTime-dTime;               
  double MTime = maxTime+dTime>imaxTime ? imaxTime : maxTime+dTime;               
  hist2D->GetXaxis()->SetRangeUser(mTime,MTime);                                  

  int npix=0;
  double Null=0;
  double Likelihood=0;
  for(int n=0; n<V; n++) {
    netpixel* pix = pwc->getPixel(cid,n);
    if(!pix->core) continue;             

    if((irate)&&(irate != int(pix->rate+0.5))) continue;  // if irate!=0 skip rate!=irate

    npix++;
           
    double sSNR=0;
    double wSNR=0;
    for(int m=0; m<nifo; m++) {                 
      sSNR += pow(pix->getdata('S',m),2);          // snr whitened reconstructed signal 00
      sSNR += pow(pix->getdata('P',m),2);          // snr whitened reconstructed signal 90
      wSNR += pow(pix->getdata('W',m),2);          // snr whitened at the detector 00     
      wSNR += pow(pix->getdata('U',m),2);          // snr whitened at the detector 90     
    }                                                                                     
    double null = wSNR-sSNR;                                                              

    int iRATE = RATE/(pix->layers-1); 
    int M=maxRate/iRATE;              
    int K=2*(maxLayers-1)/(pix->layers-1);
    double itime=((double)pix->time/(double)iRATE)/pix->layers;
    int i=(itime-iminTime)*maxRate;                            
    int j=pix->frequency*K;                                    
    if(iRATE!=irate && irate!=0) continue;                     
    Null+=null;                                                
    Likelihood+=sSNR;                                          
    int L=0;int R=1;while (R < iRATE) {R*=2;L++;}              
    for(int m=0;m<M;m++) {                                     
      for(int k=0;k<K;k++) {                                   
        if(type=='L') hist2D->SetBinContent(i+1+m,j+1+k-K/2,sSNR);
        else          hist2D->SetBinContent(i+1+m,j+1+k-K/2,null);           
      }                                                                      
    }                                                                        
  }                                                                          

  int imin=hist2D->GetNbinsX();
  int imax=0;
  int jmin=hist2D->GetNbinsY();
  int jmax=0;

  for (int i=0;i<=hist2D->GetNbinsX();i++) {
    for (int j=0;j<=hist2D->GetNbinsY();j++) {
      double X = hist2D->GetXaxis()->GetBinCenter(i)+hist2D->GetXaxis()->GetBinWidth(i)/2.;
      double Y = hist2D->GetYaxis()->GetBinCenter(j)+hist2D->GetYaxis()->GetBinWidth(j)/2.;
      double Z = hist2D->GetBinContent(i,j);
      if(Z) {
        if(i<imin) imin=i;
        if(i>imax) imax=i;
        if(j<jmin) jmin=j;
        if(j>jmax) jmax=j;
      }
      //if(Z) cout << i << " " << j << " " << X << " " << Y << " " << Z << endl;
    }
  }

  //cout << "imin " << imin << " imax " << imax << endl;
  //cout << "jmin " << jmin << " jmax " << jmax << endl;


  int di = imax-imin+1;
  int dj = jmax-jmin+1;

  //cout << "di " << di << " dj " << dj << endl;

  int ri = NDIM-di%NDIM;
  int rj = NDIM-dj%NDIM;
 
  di = di+ri; 
  dj = dj+rj; 

  //cout << "di " << di << " dj " << dj << endl;

  int lri = ri/2;
  int rri = ri-lri;
  //cout << "lri " << lri << " rri " << rri << endl;

  imin-=lri;
  imax+=rri;

  int lrj = rj/2;
  int rrj = rj-lrj;
  //cout << "lrj " << lrj << " rrj " << rrj << endl;

  jmin-=lrj;
  jmax+=rrj;

  //cout << "imin " << imin << " imax " << imax << endl;
  //cout << "jmin " << jmin << " jmax " << jmax << endl;

  int xdim = (imax-imin+1)/NDIM;
  int ydim = (jmax-jmin+1)/NDIM;
  //cout << "xdim " << xdim << " ydim " << ydim << endl;

  for(int i=0;i<NDIM;i++) for(int j=0;j<NDIM;j++) FRAME[i][j]=0;

  for (int i=0;i<=hist2D->GetNbinsX();i++) {
    for (int j=0;j<=hist2D->GetNbinsY();j++) {
      double X = hist2D->GetXaxis()->GetBinCenter(i)+hist2D->GetXaxis()->GetBinWidth(i)/2.;
      double Y = hist2D->GetYaxis()->GetBinCenter(j)+hist2D->GetYaxis()->GetBinWidth(j)/2.;
      double Z = hist2D->GetBinContent(i,j);
      if(i<imin || i>imax) continue;
      if(j<jmin || j>jmax) continue;
      int I = (i-imin)/xdim;
      int J = (j-jmin)/ydim;
      FRAME[I][J]+=Z;
      //if(Z) cout << i << " " << j << " " << X << " " << Y << " " << Z << endl;
      //if(Z) cout << I << " " << J << " " << FRAME[I][J] << endl;
    }
  }

  // fill nn frame
  for(int i=0;i<NDIM;i++) for(int j=0;j<NDIM;j++) nn_frame->push_back(FRAME[i][j]);
  // normalization
  double norm=0;
  int nn_size = nn_frame->size();
  for(int i=0;i<nn_size;i++) norm+=(*nn_frame)[i];
  for(int i=0;i<nn_size;i++) (*nn_frame)[i]/=norm;

  return;
}

void 
PlotFrame(std::vector<double>* nn_frame, int cid) {

  // get dimension of the frame 
  int nframe = nn_frame->size();
  if(nframe != pow(int(sqrt(nframe)),2)) {
    cout << "PlotFrame - Error : size is not a square number " << endl;
    gSystem->Exit(1);
  }
  nframe = sqrt(nframe);

  if(hist2D) { delete hist2D; hist2D=NULL; }
  hist2D=new TH2F("frame", "frame", nframe, 0, nframe, nframe, 0, nframe);

  hist2D->SetStats(kFALSE);
  hist2D->SetTitleFont(12);
  hist2D->SetFillColor(kWhite);

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
  hist2D->GetXaxis()->SetTitle("Time");
  hist2D->GetXaxis()->CenterTitle(true);     
  hist2D->GetYaxis()->SetTitleFont(42);      
  hist2D->GetYaxis()->SetTitle("Frequency");
  hist2D->GetYaxis()->CenterTitle(true);         

  hist2D->GetZaxis()->SetTitleOffset(0.6);
  hist2D->GetZaxis()->SetTitleFont(42);   
  hist2D->GetZaxis()->CenterTitle(true);  

  hist2D->GetXaxis()->SetLabelSize(0.03);
  hist2D->GetYaxis()->SetLabelSize(0.03);
  hist2D->GetZaxis()->SetLabelSize(0.03);

  for(int i=0;i<nframe;i++) {
    for(int j=0;j<nframe;j++) {
      //cout << i << " " << j << " " << (*nn_frame)[i*nframe+j] << endl;
      hist2D->SetBinContent(i+1,j+1,(*nn_frame)[i*nframe+j]);
    }
  }


  canvas= new TCanvas("mra", "mra", 200, 20, 600, 600);
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

  hist2D->Draw("COLZ");

  // change palette's width
  canvas->Update();        
  TPaletteAxis *palette = (TPaletteAxis*)hist2D->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.91);                                                                   
  palette->SetX2NDC(0.933);                                                                  
  palette->SetTitleOffset(0.92);                                                             
  palette->GetAxis()->SetTickSize(0.01);                                                     
  canvas->Modified();                                                                        

  char fname[1024];
  sprintf(fname, "nn_frame_display_%d.png",cid);
  cout << "write " << fname << endl;
  canvas->Print(fname);

}

