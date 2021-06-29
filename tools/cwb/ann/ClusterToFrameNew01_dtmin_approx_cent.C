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


//
// Convert netcluster principal components to NN (Neural Network) frames (NDIMxNDIM pixels)
// Author : Gabriele Vedovato
// Note : run with ACLiC -> root -l 'ClusterToFrameNew01_dtmin_approx_cent.C+("ifile.root","ofile.root")'
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
#include <algorithm>
#include "TTreeFormula.h"
TH2F* hist2D=NULL;
TCanvas* canvas;
#define REF 0.1
#define NPIX 20
#define NDIM 8
#define NBIN 20
#define MAXL 10
#define MINL 4
#define maxLayers 1025

void PlotFrame(std::vector<double>* nn_frame, int cid);

double ConvertToFrame(netcluster* pwc, int cid, int nifo, char type, int irate, double &dT,double &dF,double &Fc,int &ninf,std::vector<double>* nn_frame, std::vector<double>* dt);

//#define WATPLOT         // if defined the event monster plots are produced
//#define PLOT_FRAME      // frame display
//#define distribution

void ClusterToFrameNew01_dtmin_approx_cent(TString ifName, TString ofName) {
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
   int Ntot=nn_itree->GetEntries();

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

// add the NDIM*NDIM pixels vector (nn_frame) for NN analysis and the time vector (dt) which contains the time-duration of each selected pixel
   std::vector<double>* nn_frame = new vector<double>;    
   std::vector<double>* dt = new vector<double>;
   nn_frame->resize(NDIM*NDIM);
   dt->resize(NPIX);

   int ind=0;    
   double dT=0.;    
   double dF=0.;    
   double Fc=0.;    
   double deltat[NPIX];
 
   for (int iu=0; iu<NPIX; iu++) deltat[iu]=0.;
   int ninf=0;

//defining brannces for adding other information
   nn_otree->Branch("nnframe", &nn_frame, NDIM*NDIM, 0);
   nn_otree->Branch("dt_corepixels",&dt,NPIX,0);
   nn_otree->Branch("index",&ind,"ind/I");
   nn_otree->Branch("t_duration",&dT,"duration/D");
   nn_otree->Branch("f_duration",&dF,"f_duration/D");
   nn_otree->Branch("central_f",&Fc,"central_f/D");

// read nn_itree
   int nn_isize = nn_itree->GetEntries();
   cout << "nn_itree size : " << nn_isize << endl;
     double SNRf=0.;

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

// using the real conversion-function (ConvertToFrame)

     ind=m;
     double SNRfi=0.;
     int ID = 1;			// ID=1 because each netcluster contains only 1 event

     for (int kk=0; kk<NDIM*NDIM; kk++) (*nn_frame)[kk]=0.;
     for (int kk=0; kk<NPIX; kk++) (*dt)[kk]=0.;
     SNRfi=ConvertToFrame(pwc, m, nIFO, 'L', 0, dT, dF, Fc,ninf, nn_frame,dt);
     SNRf=SNRf+SNRfi; //evaluating the avarage weight of the latest pixel

#ifdef PLOT_FRAME                       // frame display
     PlotFrame(nn_frame, m);
#endif

//filling  the output tree
     nn_otree->Fill();
   }

//write and close the output file
   cout<<"nn_size: "<<nn_isize<<endl;
   nn_otree->Write();
   ofile.Close();
   cout<<ofName<<endl;
         
   ifile->Close();
   cout<<"number events with SNRfi<"<<REF<<": "<<ninf<<endl;
   cout<<"media importanza ultimo bin: "<<SNRf/nn_isize<<endl;
   gSystem->Exit(0);
}

double                                                                 
ConvertToFrame(netcluster* pwc, int n_ev, int nifo, char type, int irate,double &dT, double &dF, double &Fc,int &ninf, std::vector<double>* nn_frame, std::vector<double>* dt) {
  int cid=1;
  double FRAME[NDIM][NDIM];
  if(type != 'L' && type != 'N') return 0.;

  double RATE = pwc->rate;                              // original rate>4096

  std::vector<int>* vint = &(pwc->cList[cid-1]);        // pixel list

  int V = vint->size();                                 // cluster size
  if(!V) return 0.;                                                       

  bool optrate=false;
  if(irate==1) {                                        // extract optimal rate
    optrate=true;                                                              
    std::vector<int>* pv = pwc->cRate.size() ? &(pwc->cRate[cid-1]) : NULL;          
    irate = pv!=NULL ? (*pv)[0] : 0;                                           
  }                                                                            

 //Layers: number of pixels pixels contained in the frequency range (64-2048) in a particular level of decomposition

  // int minLayers=1000;
  // int maxLayers=0;   
  double minTime=1e20;
  double maxTime=0.;  
  double minFreq=1e20;
  double maxFreq=0.; 

  //in freq sta prendendo i centri dei pixels mentre nei tempi prendi i limite inf
  int num_core=0; 
  
//loop over the pixels
  for(int j=0; j<V; j++) {                      // loop over the pixels
    netpixel* pix = pwc->getPixel(cid,j);                               
    if(!pix->core) continue;                    // selection of the only core-pixels                        
    
    if((irate)&&(irate != int(pix->rate+0.5))) continue;  // if irate!=0 skip rate!=irate

//    if(pix->layers<minLayers) minLayers=pix->layers;
//    if(pix->layers>maxLayers) maxLayers=pix->layers;
    
//   int minLayers=pow(2,MINL)+1;
//   int maxLayers=pow(2,MAXL)+1;

    
    //deltat[j]=((pix->layers-1)/RATE);  
    int max_ind=0;
    double rate = RATE/(pix->layers-1); 
    double time=((double)pix->time/rate)/pix->layers;
    if(time<minTime) minTime=time;                   
    if(time>maxTime) {maxTime=time;                   
    //max_ind=j;
    //max_amp=pow(pix->getdata(),2);
    }
    num_core=num_core+1;
   //pix->frequency: nindice in freq
   //rate/2=df
   //siccome in freq inizio con pixels mezzi così(pix->frequency*rate/2.) arrivo alla freq centrale-> layers sono quindi sempre dispari
   //pix->time=itime*(pix->layer)+ifreq {ifreq e itime idici interi di array}
   //time=itime*dt[s]
   //freq=ifreq*df[Hz]

    double freq = pix->frequency*rate/2.; 
    if(freq<minFreq) minFreq=freq;        
    if(freq>maxFreq) maxFreq=freq;        
  }
  int const ncorepix=num_core;
  double deltat[ncorepix];
  for (int u=0;u<ncorepix;u++) deltat[u]=0.;
  //cout<<"nPIX: "<<NPIX<<" size: "<<dt->size()<<endl;
  //std::sort(dt->begin(),dt->begin()+NPIX);
 // for (int u=0;u<NPIX;u++) cout<<"dt: "<<(*dt)[u+n_ev*NPIX]<<endl;
 // for (int u=0;u<NPIX;u++) cout<<"evento numero "<<n_ev<<" pixel numero: "<<u<<" dt="<<deltat[u]<<endl;
  //int minRate=RATE/(maxLayers-1);
 // int maxRate=RATE/(minLayers-1);
  int minRate=RATE/(pow(2,MAXL));
  int maxRate=RATE/(pow(2,MINL));
  dT=maxTime-minTime;
 // cout<<" dT: "<<dT<<endl;
  dF=maxFreq-minFreq;
 // cout<<" dF: "<<dF<<endl;
  Fc=(maxFreq+minFreq)/2;
 // cout<<" Fc: "<<Fc<<endl;
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
  hist2D=new TH2F("WTF", "WTF", nTime, iminTime, imaxTime, 2*maxLayers-2, 0, RATE/2);
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
  int ncore=0;
  int npix=0;
  double Null2=0;
  double Likelihood2=0;
  double Null=0;
  double Likelihood=0;
 double amp=0.;
 double ttime=0.;
 TTree* tS=new TTree("Time_SNR","Time_SNR");
 tS->Branch("amplitude",&amp,"amplitude/D");
 tS->Branch("central_time",&ttime,"central_time/D");
 for(int n=0; n<V; n++) {
    netpixel* pix = pwc->getPixel(cid,n);
    if(!pix->core) continue;             
    double rate = RATE/(pix->layers-1); 
    deltat[ncore]=1./rate;
    ncore=ncore+1;
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
    amp=sSNR;
    //cout<<"pixel n: "<<n<<" sSNR: "<<sSNR;
    int iRATE = RATE/(pix->layers-1); 
    int M=maxRate/iRATE;              
    int K=2*(maxLayers-1)/(pix->layers-1);
    double itime=((double)pix->time/(double)iRATE)/pix->layers;
    ttime=itime;
    //cout<<" ttime: "<<ttime<<endl;
    tS->Fill();
    int i=(itime-iminTime)*maxRate;                            
    int j=pix->frequency*K;                                    
    if(iRATE!=irate && irate!=0) continue;                     
    Null+=null;                                                
    Likelihood+=sSNR;                                          
    int L=0;int R=1;while (R < iRATE) {R*=2;L++;}              
  } 
  //cout<<"ncore: "<<ncore<<endl;                                                                         
  //cout<<"minTime: "<<minTime<<" maxTime "<<maxTime<<endl;
  double deltaBIN=(maxTime-minTime)/(NBIN-1);
  //cout<<"deltaBIN: "<<deltaBIN<<endl;
  TH1D* h=new TH1D("SNR_distribution","SNR_distribution",NBIN,1,NBIN+1);
  //double timec[NPIX];
  double timec2[ncorepix];
 // for (int u=0; u<NPIX;u++) timec[u]=0.;
  char cuts[1024];
  double time_bin[NBIN];
  for (int nbin=0;nbin<NBIN;nbin++){
	double mintime_bin=maxTime-(nbin)*deltaBIN;
        if (nbin==0||nbin==NBIN-1) mintime_bin=mintime_bin-0.001;
        time_bin[nbin]=mintime_bin;
 //       cout<<"numero bin "<<nbin<<" mintime_bin "<<mintime_bin<<endl;  
 	 sprintf(cuts,"central_time>%5.5f",mintime_bin);
  	 TTreeFormula* treeformula=new TTreeFormula("cuts",cuts,tS);
 	 int err = treeformula->Compile(cuts);
 	 if(err) {
 		 cout << "DisplayROC.C - wrong input cuts " << cuts << endl;
 		 //return -1;
	    }
         double ybin=0.;
	 tS->SetBranchAddress("central_time",&ttime);
	 tS->SetBranchAddress("amplitude",&amp);
	 for(int n=0;n<tS->GetEntries();n++) {                     // loop over the detected events
	    tS->GetEntry(n);                          // load event #n
//	    timec[n]=ttime;
            timec2[n]=ttime;
    //        cout<<"ttime: "<<ttime<<" amp: "<<amp<<" di evento nume: "<<n+1<<" su "<<tS->GetEntries()<<endl;
             if(treeformula!=NULL) {
	      if(treeformula->EvalInstance()==0) {      // skip entry if it does not satisfy selection cuts
	    //    cout << "Skip entry number : " << n << endl;
	        continue;
	      }
	      ybin=ybin+amp;
  //            cout<<"n "<<n<<" min time bin: "<<mintime_bin<<" central time: "<<ttime<<" amp: "<<amp<<" xBIN: "<<mintime_bin<<" yBIN: "<<ybin<<endl;
	    }
	 }
         h->Fill(nbin+1,ybin);
	 ybin=0.;
  }
std::sort(timec2,timec2+ncorepix);
//std::sort(timec,timec+NPIX);
int index[ncorepix];                                       
TMath::Sort(ncorepix,timec2,index); 
  for (int u=0;u<NPIX;u++){
	if(u<ncorepix) {
		//dt->push_back(deltat[index[u]]);
		(*dt)[u] = deltat[index[u]];
		//cout<<"dt: "<<deltat[index[u]]<<" central_time "<<timec2[u]<<endl;
		}
	//else dt->push_back(0.);
	else (*dt)[u] = 0.;
//	cout<<""<<endl;
//		 cout<<"deltat "<<deltat[u]<<" dt: "<<(*dt)[u]<<endl;
		}
#ifdef distribution
TCanvas* h_canv=new TCanvas("SNR","SNR",0,0,500,500);
h->Draw();
char h_char[512];
sprintf(h_char,"SNRisto/event%i.png",n_ev);
h_canv->SaveAs(h_char);
#endif
double SNRfi0=0.;
double SNRfi=0.;
int js=1;
while (h->GetBinContent(js)==0&&js<NBIN) js=js+1;
SNRfi0=(h->GetBinContent(js)/h->GetBinContent(NBIN));
int indtemp=0;
double ampTOT=0.;
double ampFin=0.;
double amp2[ncorepix];
for (int u=0;u<ncorepix;u++) amp2[u]=0.;
double t_limit=0.;
    for(int iu=0;iu<ncorepix;iu++){
	 tS->SetBranchAddress("central_time",&ttime);
	 tS->SetBranchAddress("amplitude",&amp);
	 char cuts2[1024];
 	 sprintf(cuts2,"central_time>%5.5f",timec2[ncorepix-1-iu]-0.001);
  	 TTreeFormula* treeformula2=new TTreeFormula("cuts2",cuts2,tS);
 	 int err2 = treeformula2->Compile(cuts2);
 	 if(err2) {
 		 cout << "DisplayROC.C - wrong input cuts " << cuts2 << endl;
 		 //return -1;
		break;
	    }
	 ampTOT=0.;
	 for(int n=0;n<tS->GetEntries();n++) {                     // loop over the detected events
	    tS->GetEntry(n);                          // load event #n
            //cout<<"tempo limite analizzato "<<timec2[ncorepix-1-iu]<<" ttime "<<ttime<<" evento nume: "<<n+1<<" su "<<tS->GetEntries()<<endl;
             if(treeformula2!=NULL) {
	      if(treeformula2->EvalInstance()==0) {      // skip entry if it does not satisfy selection cuts
	       // cout << "Skip entry number : " << n << endl;
	        continue;
	      }
	      ampTOT=ampTOT+amp;
              //cout<<"timec "<<timec2[ncorepix-1-iu]<< " ampTOT "<<ampTOT<<" amp "<<amp<<endl;
	  }
	  }
	amp2[iu]=ampTOT;
    //if(ampTOT>REF&&iu!=0) t_limit= timec2[ncorepix-iu];//così t_limit è l'ultimo da togliere
	      if((ampTOT/h->GetBinContent(NBIN))>REF) {
			if(iu>0) ampFin=amp2[iu-1];
			//cout<<"fial amplitude"<<ampFin<<endl;
			t_limit= timec2[ncorepix-1-iu];//così t_limit è il primo da lasciare
			//if(iu>0) cout<<"ultimo tempo da buttare:"<<timec2[ncorepix-iu];
			break;
			}
	     // if((ampTOT/h->GetBinContent(NBIN))>REF&&iu==0) t_limit=0.;
 
	}
//cout<<" primo timo da tenere:"<<t_limit<<endl;
//SNRfi=ampFin/(h->GetBinContent(NBIN-1));
SNRfi=ampFin/(h->GetBinContent(NBIN));
//cout<<"SNRfi0: "<<SNRfi0<<" SNRfi: "<<SNRfi<<endl;

 for(int n=0; n<V; n++) {
    netpixel* pix = pwc->getPixel(cid,n);
    if(!pix->core) continue;             
    if((irate)&&(irate != int(pix->rate+0.5))) continue;  // if irate!=0 skip rate!=irate
    double rate = RATE/(pix->layers-1); 
    double sSNR2=0;
    double wSNR2=0;
    //cout<<"limit time: "<<t_limit<<" time from direct pixel: "<<((double)pix->time/rate)/pix->layers<<endl;
    //cout<<"likelihhod: "<<Likelihood2<<endl;
    if (SNRfi>0.&&SNRfi<=REF&&(((double)pix->time/rate)/pix->layers)>t_limit) continue;
    for(int m=0; m<nifo; m++) {                 
      sSNR2 += pow(pix->getdata('S',m),2);          // snr whitened reconstructed signal 00
      sSNR2 += pow(pix->getdata('P',m),2);          // snr whitened reconstructed signal 90
      wSNR2 += pow(pix->getdata('W',m),2);          // snr whitened at the detector 00     
      wSNR2 += pow(pix->getdata('U',m),2);          // snr whitened at the detector 90     
    }                                                                                     
    double null2 = wSNR2-sSNR2;                                                              
    int iRATE = RATE/(pix->layers-1); 
    int M=maxRate/iRATE;              
    int K=2*(maxLayers-1)/(pix->layers-1);
    double itime=((double)pix->time/(double)iRATE)/pix->layers;
    ttime=itime;
    int i=(itime-iminTime)*maxRate;                            
    int j=pix->frequency*K;                                    
    if(iRATE!=irate && irate!=0) continue;                     
    Null2+=null2;                                                
    Likelihood2+=sSNR2;
    //cout<<"likelihood dopo aggiornamento: "<<Likelihood2<<endl;                                       
    int L=0;int R=1;while (R < iRATE) {R*=2;L++;}              
    for(int m=0;m<M;m++) {                                     
      for(int k=0;k<K;k++) {                                   
        double SNRsum=hist2D->GetBinContent(i+1+m,j+1+k-K/2);
	if(type=='L') hist2D->SetBinContent(i+1+m,j+1+k-K/2,sSNR2+SNRsum);
        else          hist2D->SetBinContent(i+1+m,j+1+k-K/2,null2+SNRsum);           
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
  if (di%NDIM==0) ri=0;
  if (dj%NDIM==0) rj=0;
/*    int lri = ri/2;
  int rri = ri-lri;
  //cout << "lri " << lri << " rri " << rri << endl;

  imin-=lri;
  imax+=rri;

  int lrj = rj/2;
  int rrj = rj-lrj;
  //cout << "lrj " << lrj << " rrj " << rrj << endl;

  jmin-=lrj;
  jmax+=rrj;
*/
//cout << "di: " << di << " ri: " << ri << endl;
//  cout << "dj: " << dj << " rj: " << rj << endl;
//  cout << "imin " << imin << " imax " << imax << endl;
//  cout << "jmin " << jmin << " jmax " << jmax << endl;
  bool ripi=false;
  bool ripj=false;
//cout<<"resto diverso da zero"<<endl; 

 int lri = ri/2;
 int rri = ri-lri;
 int lrj = rj/2;
 int rrj = rj-lrj;


  if (ri<4 || di<=4){
          di = di+ri; 
	  imin-=lri;
	  imax+=rri;
	  //imin-=ri;
	}

  else {
	  imin+=di%NDIM;
          di = di-di%NDIM;
          ripi=true;
	}

   if (rj<4 || dj<=4){
	  dj = dj+rj; 
	  jmin-=lrj;
	  jmax+=rrj;

	  //jmax+=rj;
	}
   else {
	  jmax-=dj%NDIM;
          dj = dj-dj%NDIM;
          ripj=true;
	}

//  cout << "imin " << imin << " imax " << imax << endl;
//  cout << "jmin " << jmin << " jmax " << jmax << endl;
//cout<<" def limiti"<<endl;
  int xdim = (imax-imin+1)/NDIM;
  int ydim = (jmax-jmin+1)/NDIM;
  //cout << "xdim " << xdim << " ydim " << ydim << endl;
//cout<<" def xdim,ydim"<<endl;

  for(int i=0;i<NDIM;i++) for(int j=0;j<NDIM;j++) FRAME[i][j]=0;

  for (int i=0;i<=hist2D->GetNbinsX();i++) {
//cout<<"i, j "<<i<<" "<<" i min "<<imin<<xmax<<" xdim "<<xdim<<" jmin "<<jmin<<" ydim "<<ydim<<endl;
    for (int j=0;j<=hist2D->GetNbinsY();j++) {

      bool ripi2=false;
      bool ripj2=false;

      double X = hist2D->GetXaxis()->GetBinCenter(i)+hist2D->GetXaxis()->GetBinWidth(i)/2.;
      double Y = hist2D->GetYaxis()->GetBinCenter(j)+hist2D->GetYaxis()->GetBinWidth(j)/2.;
      double Z = hist2D->GetBinContent(i,j);
      if(i>imax) continue;
      if(j<jmin) continue;

     // if((i<imin)&&(ri<4|| di<=4)) continue;
     // if((j>jmax)&&(rj<4|| dj<=4)) continue;

     // if((i<imin)&&!(ri<4|| di<=4)) i=i+xdim;
     // if((j>jmax)&&!(rj<4|| dj<=4)) j=j-ydim; 

      if (i<imin && ripi==false) continue;
      if (j<jmin && ripj==false) continue;

      if (i<imin && ripi==true) {i=i+xdim;ripi2=true;}
      if (j<jmin && ripj==true) {j=j-ydim;ripj2=true;}

//cout<<"i, j "<<i<<" "<<j<<" i min "<<imin<<" xdim "<<xdim<<" jmin "<<jmin<<" ydim "<<ydim<<endl;
      int I = (i-imin)/xdim;
      int J = (j-jmin)/ydim;
//cout << I << " " << J << endl;
      FRAME[I][J]+=Z;
//cout << "sono 	ui" << endl;

	//if((i-xdim<imin)&&!(ri<4|| di<=4)) i=i-xdim;
	//if((j+ydim>jmax)&&!(rj<4|| dj<=4)) j=j+ydim;

        if (ripi2==true) i=i-xdim;
        if (ripj2==true) j=j+ydim;
      //if(Z) cout << i << " " << j << " " << X << " " << Y << " " << Z << endl;
      //if(Z) cout << I << " " << J << " " << FRAME[I][J] << endl;
    }
  }

  // fill nn frame
  //for(int i=0;i<NDIM;i++) for(int j=0;j<NDIM;j++) nn_frame->push_back(FRAME[i][j]);
  for(int i=0;i<NDIM;i++) for(int j=0;j<NDIM;j++) (*nn_frame)[i*NDIM+j] = FRAME[i][j];
  // normalization
  double norm=0;
  int nn_size = nn_frame->size();
  for(int i=0;i<nn_size;i++) norm+=(*nn_frame)[i];
  for(int i=0;i<nn_size;i++) (*nn_frame)[i]/=norm;
  if (SNRfi<=REF) ninf=ninf+1;
  delete h;
  return SNRfi0;

//cout<<" fine funzione" <<endl;
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

