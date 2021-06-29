// this macro shows how to read the eBBH results from output ROOT file produced with CWB_Plugin_eBBH.C
// Author : G.Vedovato

#define NCHIRP_MAX	4

// ------------------------------------------------------------
// ebbh structure
// ------------------------------------------------------------

struct EBBH {                           // structure for output for estimated parameters

  TString name; 

  std::vector<TString> ifo; 

  float isnr;
  float osnr;

  double seggps;                        // segment start time

  float parms[NIFO_MAX];                // null pixel statistic, for each detector 

  float mchirp;
  float mass[2];
  float spin[6];                        // spin
  float e0;
  float rp0;
  float dist;
  float redshift;

  float ch_mass[NCHIRP_MAX];            // chirp mass
  float ch_tmerger[NCHIRP_MAX];         // merger time
  float ch_merr[NCHIRP_MAX];            // mchirp error
  float ch_mchi2[NCHIRP_MAX];           // mchirp chi2
  float ch_energy[NCHIRP_MAX];          // chirp energy

  float ei;				// eccentricity index

  netcluster nc;      			// netcluster

  wavearray<double> wdat[NIFO_MAX];     // whitened data (in the vREC time range)
  wavearray<double> winj[NIFO_MAX];     // whiten injected waveforms
  wavearray<double> sinj[NIFO_MAX];     // strain injected waveforms
  wavearray<double> wrec[NIFO_MAX];     // whiten reconstructed waveforms
  wavearray<double> srec[NIFO_MAX];     // strain injected waveforms

  wavearray<double> wsim[2];     	// strain simulated waveforms hp,hc
};

// ------------------------------------------------------------
// functions
// ------------------------------------------------------------

int ReadDataFromROOT(TString ifile, EBBH* ebbh);, 

void PlotMultiChirp(EBBH* ebbh);
void PlotMultiChirp2(EBBH* ebbh, int ifo=0, TString ptype="", TString wtype="");
void PlotMChirp(EBBH* ebbh);
void PlotChirp(EBBH* ebbh, int mch_order, double& tmerger, double& xmin, double& xmax, EColor color);

double GetTimeBoundaries(wavearray<double> x, double P, double& bT, double& eT);
wavearray<double> GetAlignedWaveform(wavearray<double>* wf, wavearray<double>* wref);
wavearray<double> GetDifWaveform(wavearray<double>* wf1, wavearray<double>* wf2);

double GetFitParameter(double mchirp); 
void PlotCluster(watplot* WTS, netcluster* pwc, double fLow, double fHigh);

void FindChirp(EBBH* ebbh, int mch_order, double& tmerger, double& xmin, double& xmax, int ifo, TString ptype, TString wtype);

void FillHistogram(EBBH* ebbh, TH2F&* hN, TH2F&* hL, TH2F&* heT, TH2F&* heF, double zmax_thr=0.);
void PlotHistogram(EBBH* ebbh, TH2F&* hH);
double ChirpIntegralHistogram(EBBH* ebbh, TH2F&* h2d, double dMc);

void DrawMChirpInstantaneousFrequency(EBBH* ebbh, TString wtype, int ifo);

void GetCBC(EBBH* ebbh);
void GetEBBH(EBBH* ebbh);

void Init();
void Clear();

// -----------------------------------------------------------------------------------------------------------
// defines
// -----------------------------------------------------------------------------------------------------------

#define FLOW	16
#define FHIGH	1024

// default threshold (used in network::likelihoodWP)
/*
#define CHI2_THR	2.0
#define TMERGER_CUT	1.e20
#define ZMAX_THR	0.0
*/


//#define CHI2_THR       -2.5	// if negative the chirp pixels are tagged with pix->likelihood=0 after the selection
//#define TMERGER_CUT	0.1
#define CHI2_THR       -4.5	// if negative the chirp pixels are tagged with pix->likelihood=0 after the selection
#define TMERGER_CUT	0.01
#define ZMAX_THR	0.2



// ------------------------------------------------------------
// global variables
// ------------------------------------------------------------

int  gTYPE;
int  gENTRY;
bool gPLOT;
TString gOFILE;

CWB::mdc* MDC;
watplot* WTS;
watplot* PCH;
CWB::STFT* stft;
TF1* fchirp;

TH2F* hN;
TH2F* hL;
TH2F* heT;
TH2F* heF;

TF1* xfit[NCHIRP_MAX];
TGraphErrors* xgr[NCHIRP_MAX];

TCanvas* canvas;

int EccentricityIndex(TString ifroot, int gtype=0, int entry=0, int ifo=0, EBBH* xebbh=NULL) {

  // if gtype<0 the plot is saved with name with the following format : gTYPE_gtype_file_name.png

  // gtype =  0  -> no plots, only computation of ei
  // gtype =  1  -> plot chirp fit function with rec spectrogram
  // gtype =  2  -> plot chirp fit function with likelihood
  // gtype =  3  -> plot chirp fit function with likelihood without the main chirp
  //                show residuals after chirp1 removal ( only with PlotMultiChirp2 & "like" )
  // gtype =  4  -> plot with inj spectrogram
  // gtype =  5  -> plot likelihood TF after with pixels selected with zmax_thr 
  // gtype =  6  -> plot chirp mass fit in f^-3/8 vs time plane
  // gtype =  7  -> plot chirp mass fit in f^-3/8 vs time plane for all chirp
  // gtype =  8  -> plot ifo waveform in time
  // gtype =  9  -> plot ifo waveform in frequency
  // gtype =  10 -> plot whitened data spectrogram
  // gtype =  11 -> plot inj waveform  spectrogram
  // gtype =  12 -> plot rec waveform  spectrogram
  // gtype =  13 -> plot inj waveform  instantaneous frequency mchirp 
  // gtype =  14 -> plot rec waveform  instantaneous frequency mchirp 
  // gtype =  15 -> plot sim waveform (GetCBC)  instantaneous frequency mchirp  
  // gtype =  16 -> plot sim waveform (GetEBBH) instantaneous frequency mchirp  

  Init();

  if(gtype<-16 || gtype>16) {cout << "Error, gtype parameter value not allowed : [-1:11]" << endl;exit(1);}
  if(ifroot=="")            {cout << "Error, input root file name not defined" << endl;exit(1);}

  gPLOT = gtype<0 ? true : false;
  gTYPE = abs(gtype);

  if(gPLOT) gROOT->SetBatch(true);

  gENTRY = entry;

  TString fDir  = gSystem->DirName(ifroot);
  TString fName = gSystem->BaseName(ifroot);

  if(fName.Contains("*")) {
    TString fTag = fName;
    fTag.ReplaceAll("*",""); 
    vector<TString> fList = CWB::Toolbox::getFileListFromDir(fDir, ".root", "wave_", fTag, true);
    if(fList.size()==0) {cout << "Error, no file selected" << endl;exit(1);}
    if(fList.size()>1)  {cout << "Error, multiple files selected" << endl;exit(1);}

    ifroot = fList[0];
    fDir  = gSystem->DirName(ifroot);
    fName = gSystem->BaseName(ifroot);
  }

  CWB::Toolbox::checkFile(ifroot);

  gOFILE = fName;
  gOFILE.ReplaceAll(".root",".png");	// output file name for plot
  gOFILE = TString::Format("gTYPE_%d_%s",gTYPE,gOFILE.Data());	
  if(gPLOT) cout << endl << "Output Plot File Name : " << gOFILE << endl << endl;

  CWB::Toolbox::checkFile(ifroot);

  EBBH ebbh;
  MDC = new CWB::mdc();

  // read ebbh from input root file
  int err = ReadDataFromROOT(ifroot, &ebbh);
  if(err) return -100000;

  if(gTYPE==15) GetCBC(&ebbh);     	// get simulated CBC waveforms
  if(gTYPE==16) GetEBBH(&ebbh);     	// get simulated EBBH waveforms

  // set waveforms start time at start segment time
  ebbh.wdat[ifo].start(ebbh.seggps ? ebbh.wdat[ifo].start()-ebbh.seggps : 0);
  ebbh.winj[ifo].start(ebbh.seggps ? ebbh.winj[ifo].start()-ebbh.seggps : 0);
  ebbh.wrec[ifo].start(ebbh.seggps ? ebbh.wrec[ifo].start()-ebbh.seggps : 0);
  ebbh.srec[ifo].start(ebbh.seggps ? ebbh.srec[ifo].start()-ebbh.seggps : 0);

  if(gTYPE==0)             {PlotMultiChirp2(&ebbh); if(xebbh) {*xebbh=ebbh; Clear(); return 0;} else exit(0);}
  if(gTYPE==1)             {PlotMultiChirp2(&ebbh, ifo, "stft", "wrec");return  0;}
  if(gTYPE==2 || gTYPE==3) {PlotMultiChirp2(&ebbh, ifo, "like", "wrec");return  0;}
  if(gTYPE==4)             {PlotMultiChirp2(&ebbh, ifo, "stft", "winj");return  0;}
  if(gTYPE==5) {
                           FillHistogram(&ebbh, hN, hL, heT, heF, ZMAX_THR);
                           PlotHistogram(&ebbh, hL); return 0;
/*
                           FillHistogram(&ebbh, hN, hL, heT, heF, 0);
                           double cenergy = ChirpIntegralHistogram(&ebbh, hL, 0);
                           PlotHistogram(&ebbh, hL); return 0;
                           TH1F* hC = new TH1F("hC","hC",41,-20,20);
                           for(int i=0;i<=40;i++) {
                             double dMc = i-20;
                             double cenergy = ChirpIntegralHistogram(&ebbh, hL, dMc);
                             cout << "CHIRP ENERGY : " << dMc << " -> " << cenergy << endl;
                             hC->SetBinContent(i,cenergy);
                           }
                           hC->Draw();
                           return 0; 
*/
                           //PlotHistogram(&ebbh, hL); return 0;
                           //PlotHistogram(&ebbh, heF); return 0;
                           //PlotHistogram(&ebbh, heT); return 0;
                           //PlotHistogram(&ebbh, hN); return 0;
  }
  if(gTYPE==6)             {PlotMChirp(&ebbh); return 0;}
  if(gTYPE==7)             {PlotMultiChirp(&ebbh);return 0;}
  if(gTYPE==8) {
  	                   MDC->SetZoom(0.999);
                           watplot* plot=MDC->Draw(ebbh.winj[ifo],MDC_TIME,"ALP ZOOM",kRed);
                           MDC->Draw(ebbh.wrec[ifo],MDC_TIME,"SAME",kBlack);
                           if(plot) plot->graph[0]->SetTitle(TString::Format("%s : Rec(black), Inj(red))",ebbh.name.Data()));
                           if(gPLOT) {stft->Print(gOFILE);exit(0);} else return 0;
  }
  if(gTYPE==9) {
  	                   MDC->SetZoom(0.999);
                           watplot* plot=MDC->Draw(ebbh.winj[ifo],MDC_FFT,"ALP ZOOM",kRed);
                           MDC->Draw(ebbh.wrec[ifo],MDC_FFT,"SAME",kBlack);
                           if(plot) plot->graph[0]->SetTitle(TString::Format("%s : Rec(black), Inj(red))",ebbh.name.Data()));
                           if(gPLOT) {MDC->GetWatPlot()->canvas->Print(gOFILE);exit(0);} else return 0;
  }
  if(gTYPE==10)            {MDC->Draw(ebbh.wdat[ifo],MDC_TF);if(gPLOT) {MDC->GetSTFT()->Print(gOFILE);exit(0);} else return 0;}
  if(gTYPE==11)            {MDC->Draw(ebbh.winj[ifo],MDC_TF);if(gPLOT) {MDC->GetSTFT()->Print(gOFILE);exit(0);} else return 0;}
  if(gTYPE==12)            {MDC->Draw(ebbh.wrec[ifo],MDC_TF);if(gPLOT) {MDC->GetSTFT()->Print(gOFILE);exit(0);} else return 0;}

  if(gTYPE==13)            {DrawMChirpInstantaneousFrequency(&ebbh, "winj", ifo); return;}
  if(gTYPE==14)            {DrawMChirpInstantaneousFrequency(&ebbh, "wrec", ifo); return;}
  if(gTYPE==15)            {DrawMChirpInstantaneousFrequency(&ebbh, "wsim", ifo); return;}

  return 0;
}

void Init() {

  MDC = NULL;
  WTS = NULL;
  PCH = NULL;
  stft = NULL;

  hN = NULL;
  hL = NULL;
  heT = NULL;
  heF = NULL;

  for(int i=0;i<NCHIRP_MAX;i++) {
    xfit[i] = NULL;
    xgr[i] = NULL;
  }

  fchirp = NULL;

  canvas = NULL;
}

void Clear() {

  if(MDC != NULL)  delete MDC;
  if(WTS != NULL)  delete WTS;
  if(PCH != NULL)  delete PCH;
  if(stft != NULL) delete stft;

  if(hN != NULL)   delete hN;
  if(hL != NULL)   delete hL;
  if(heT != NULL)  delete heT;
  if(heF != NULL)  delete heF;

  for(int i=0;i<NCHIRP_MAX;i++) {
    if(xfit[i] != NULL) delete xfit[i];
    if(xgr[i] != NULL)  delete xgr[i];
  }

  if(fchirp != NULL) delete fchirp;

  if(canvas != NULL) delete canvas;
}

int ReadDataFromROOT(TString ifile, EBBH* ebbh) { 
//
// ----------------------------------------------------
// ROOT Output PE Parameters
// ----------------------------------------------------

  TFile* froot = new TFile(ifile,"READ");
  TTree* itree = froot->Get("waveburst");
  if(itree==NULL) {cout << "ReadDataFromROOT - Error : no waveburst present in the file" << endl;return 1;}

  // get detector list
  TList* list = itree->GetUserInfo();
  int nIFO=list->GetSize();
  if (nIFO==0) {cout << "ReadDataFromROOT - Error : no ifo present in the tree" << endl;return 1;}
  for (int n=0;n<list->GetSize();n++) {
    detector* pDetector;
    pDetector = (detector*)list->At(n);
    ebbh->ifo.push_back(pDetector->Name);
    detectorParams dParams = pDetector->getDetectorParams();
    //pDetector->print();                                                       
  }

  itree->SetBranchAddress("ndim",&nIFO);
  itree->GetEntry(gENTRY);

  double seggps=0;			   // segment start time
  float crate;				   // netcluster::rate
  float parms[NIFO_MAX];                   // ebbh parameters
  float eBBH[4];                           // ebbh parameters
  float range[2];                          // distance (Mpc)
  float iSNR[NIFO_MAX];
  float oSNR[NIFO_MAX];	
  std::vector<netpixel>* cluster;
  wavearray<double>* wDAT[NIFO_MAX];
  wavearray<double>* wINJ[NIFO_MAX];
  wavearray<double>* wREC[NIFO_MAX];
  wavearray<double>* sREC[NIFO_MAX];

  for(int n=0;n<nIFO;n++) {
    wDAT[n] = new wavearray<double>;
    wINJ[n] = new wavearray<double>;
    wREC[n] = new wavearray<double>;
    sREC[n] = new wavearray<double>;
  }

  cluster = new std::vector<netpixel>;

  itree->SetBranchAddress("mass",ebbh->mass);
  itree->SetBranchAddress("eBBH",eBBH);
  itree->SetBranchAddress("range",range);
  itree->SetBranchAddress("spin",ebbh->spin);
  itree->SetBranchAddress("iSNR",iSNR);
  itree->SetBranchAddress("oSNR",oSNR);

  itree->SetBranchAddress("ebbh_parms",parms);
  for(int n=0;n<nIFO;n++) {
    itree->SetBranchAddress(TString::Format("ebbh_wDAT_%d",n).Data(),&wDAT[n]);
    itree->SetBranchAddress(TString::Format("ebbh_wINJ_%d",n).Data(),&wINJ[n]);
    itree->SetBranchAddress(TString::Format("ebbh_wREC_%d",n).Data(),&wREC[n]);
    itree->SetBranchAddress(TString::Format("ebbh_sREC_%d",n).Data(),&sREC[n]);
  }

  itree->SetBranchAddress("ebbh_ch_mass",   ebbh->ch_mass);
  itree->SetBranchAddress("ebbh_ch_tmerger",ebbh->ch_tmerger);
  itree->SetBranchAddress("ebbh_ch_merr",   ebbh->ch_merr);
  itree->SetBranchAddress("ebbh_ch_mchi2",  ebbh->ch_mchi2);
  itree->SetBranchAddress("ebbh_ch_energy", ebbh->ch_energy);

  itree->SetBranchAddress("ebbh_seggps",&seggps);
  itree->SetBranchAddress("ebbh_crate",&crate);
  itree->SetBranchAddress("ebbh_cluster",&cluster);

  itree->GetEntry(gENTRY);

  ebbh->e0 = eBBH[1];
  ebbh->rp0 = eBBH[2];
  ebbh->redshift = eBBH[3];
  ebbh->dist = range[1];
  ebbh->isnr = 0;
  for(int n=0;n<nIFO;n++) ebbh->isnr += iSNR[n];
  ebbh->isnr = sqrt(ebbh->isnr);
  ebbh->osnr = 0;
  for(int n=0;n<nIFO;n++) ebbh->osnr += oSNR[n];
  ebbh->osnr = sqrt(ebbh->osnr);

  for(int n=0;n<nIFO;n++) {
    ebbh->wdat[n] = GetAlignedWaveform(wDAT[n],wREC[n]);
    ebbh->winj[n] = GetAlignedWaveform(wINJ[n],wREC[n]);
    ebbh->wrec[n] = GetAlignedWaveform(wREC[n],wREC[n]);
    ebbh->srec[n] = GetAlignedWaveform(sREC[n],wREC[n]);
  }

  cout << "segment start time " << seggps << endl;
  cout << "cluster rate       " << crate << endl;
  cout << "cluster size       " << cluster->size() << endl;
 
  if(cluster->size()==0) {cout << "ReadDataFromROOT - Error : no clusters in the root file !!!" << endl;return 1;}

  ebbh->nc.rate=crate;
  std::vector<int> pList;
  for(int i=0;i<cluster->size();i++) {
    ebbh->nc.pList.push_back((*cluster)[i]);
    pList.push_back(i);
  }
  ebbh->nc.cList.push_back(pList);
  clusterdata cd;  // dummy cluster data
  ebbh->nc.cData.push_back(cd);

  ebbh->seggps = seggps;

  ebbh->mchirp = pow(ebbh->mass[0]*ebbh->mass[1],3./5.)/pow(ebbh->mass[0]+ebbh->mass[1],1./5.);

  char name[1024];
  sprintf(name,"source : m1=%0.1f m2=%0.1f mchirp=%0.1f e0=%0.1f rp0=%0.1f dist=%0.1f (Mpc) redshift=%0.1f isnr=%0.1f osnr=%0.1f",
                ebbh->mass[0], ebbh->mass[1], ebbh->mchirp, ebbh->e0, ebbh->rp0, ebbh->dist, ebbh->redshift, ebbh->isnr, ebbh->osnr);
  ebbh->name = name;

  for(int n=0;n<nIFO;n++) {
    if(wDAT[n]) delete wDAT[n];
    if(wINJ[n]) delete wINJ[n];
    if(wREC[n]) delete wREC[n];
    if(sREC[n]) delete sREC[n];
  }
  if(cluster) delete cluster;
  if(itree)   delete itree;
  if(froot)   delete froot;

  return 0;
}


double GetTimeBoundaries(wavearray<double> x, double P, double& bT, double& eT) {

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

wavearray<double> GetAlignedWaveform(wavearray<double>* wf1, wavearray<double>* wf2) {

   wavearray<double> wf = *wf2;
   wf=0;

   if(wf1==NULL)      return wf;
   if(wf1->size()==0) return wf;

   double R=wf1->rate();

   double b_wf1 = wf1->start();
   double e_wf1 = wf1->start()+wf1->size()/R;
   double b_wf2 = wf2->start();
   double e_wf2 = wf2->start()+wf2->size()/R;

   int o_wf1 = b_wf1>b_wf2 ? 0 : int((b_wf2-b_wf1)*R+0.5);
   int o_wf2 = b_wf1<b_wf2 ? 0 : int((b_wf1-b_wf2)*R+0.5);

   double startXCOR = b_wf1>b_wf2 ? b_wf1 : b_wf2;
   double endXCOR   = e_wf1<e_wf2 ? e_wf1 : e_wf2;
   int sizeXCOR  = int((endXCOR-startXCOR)*R+0.5);

   for(int i=0;i<sizeXCOR;i++) wf[i+o_wf2] = wf1->data[i+o_wf1];

   return wf;
}

wavearray<double> GetDifWaveform(wavearray<double>* wf1, wavearray<double>* wf2) {

   double R=wf1->rate();

   double b_wf1 = wf1->start();
   double e_wf1 = wf1->start()+wf1->size()/R;
   double b_wf2 = wf2->start();
   double e_wf2 = wf2->start()+wf2->size()/R;

   int o_wf1 = b_wf1>b_wf2 ? 0 : int((b_wf2-b_wf1)*R+0.5);
   int o_wf2 = b_wf1<b_wf2 ? 0 : int((b_wf1-b_wf2)*R+0.5);

   double startXCOR = b_wf1>b_wf2 ? b_wf1 : b_wf2;
   double endXCOR   = e_wf1<e_wf2 ? e_wf1 : e_wf2;
   int sizeXCOR  = int((endXCOR-startXCOR)*R+0.5);

   wavearray<double> wfdif(sizeXCOR);
   wfdif=0.;
   wfdif.rate(R);
   wfdif.start(b_wf1+double(o_wf1)/R);

   for(int i=0;i<sizeXCOR;i++) wfdif[i] = wf1->data[i+o_wf1] - wf2->data[i+o_wf2];

   return wfdif;
}

double GetFitParameter(double mchirp) {

  const double G  = watconstants::GravitationalConstant();
  const double SM = watconstants::SolarMass();
  const double C  = watconstants::SpeedOfLightInVacuo();
  const double Pi = TMath::Pi();


  double A = (96./5.) * pow(Pi,8./3.) * pow(G*SM*mchirp/C/C/C, 5./3);
  double B = 3./8.;

  double p0 = -A/B;

  return p0;
}

void PlotCluster(watplot* WTS, netcluster* pwc, double fLow, double fHigh) {

  double RATE = pwc->rate;                        // original rate
  std::vector<int>* vint = &(pwc->cList[0]);        // pixel list
  int V = vint->size();                           // cluster size
  for(int j=0; j<V; j++) {                      // loop over the pixels
    netpixel* pix = pwc->getPixel(1,j);
    if(!pix->core) continue;
    if(pix->likelihood==0) {
      for(int m=0; m<2; m++) {
        pix->setdata(0,'S',m);          // snr whitened reconstructed signal 00
        pix->setdata(0,'P',m);          // snr whitened reconstructed signal 90
        pix->setdata(0,'W',m);          // snr whitened at the detector 00
        pix->setdata(0,'U',m);          // snr whitened at the detector 90
      }
    }
  }
  WTS->plot(pwc, 1, 2, 'L', 0, const_cast<char*>("COLZ"));
  WTS->hist2D->GetYaxis()->SetRangeUser(fLow, fHigh);
  WTS->canvas->SetLogy();
  return;
}

void PlotMultiChirp(EBBH* ebbh) {

  EColor color[NCHIRP_MAX] = {kBlack, kBlue, kGreen, kMagenta};

  double mchirp[NCHIRP_MAX];
  double energy[NCHIRP_MAX];

  canvas = new TCanvas;
  canvas->SetGridx();
  canvas->SetGridy();

  // plot multi chirp likelihood tf  

  double tmerger=0;
  double xmin=0;
  double xmax=0;

  for(int i=0; i<NCHIRP_MAX; i++) {
    PlotChirp(ebbh, i, tmerger, xmin, xmax, color[i]); 
    mchirp[i] = ebbh->ch_mass[i];
    energy[i] = ebbh->ch_energy[i];
  }

  char title[256];
  sprintf(title,"chirp mass : rec =  %3.3f  /   %3.3f  /   %3.3f  /   %3.3f - chirp_energy =  %3.2f  /  %3.2f  /  %3.2f  /  %3.2f",
              mchirp[0],mchirp[1],mchirp[2],mchirp[3],energy[0],energy[1],energy[2],energy[3]);
  xgr[0]->SetTitle(title);

  if(gPLOT) {canvas->Print(gOFILE);exit(0);}

  return;
}

void PlotChirp(EBBH* ebbh, int mch_order, double& tmerger, double& xmin, double& xmax, EColor color) {

  double chi2_thr    = CHI2_THR;
  double tmerger_cut = TMERGER_CUT;
  double zmax_thr    = ZMAX_THR;

  double echirp = ebbh->nc.mchirp5(1,zmax_thr,chi2_thr,tmerger_cut);

  clusterdata* pcd = &ebbh->nc.cData[0];

  ebbh->ch_mass[mch_order]   = pcd->mchirp; 
  ebbh->ch_merr[mch_order]   = pcd->mchirperr; 
  ebbh->ch_mchi2[mch_order]  = pcd->chi2chirp; 
  ebbh->ch_energy[mch_order] = echirp;

  TGraphErrors* gr = &(pcd->chirp);

  xgr[mch_order] = new TGraphErrors(gr->GetN(),gr->GetX(),gr->GetY(),gr->GetEX(),gr->GetEY());

  char title[256];
  sprintf(title,"chirp mass : rec = %3.3f [%3.2f] , chi2 = %3.2f",
              pcd->mchirp,pcd->mchirperr,pcd->chi2chirp);

  xgr[mch_order]->SetTitle(title);
  xgr[mch_order]->GetHistogram()->SetStats(kFALSE);
  xgr[mch_order]->GetHistogram()->SetTitleFont(12);
  xgr[mch_order]->SetFillColor(kWhite);
  xgr[mch_order]->SetLineColor(kBlack);
  xgr[mch_order]->GetXaxis()->SetNdivisions(506);
  xgr[mch_order]->GetXaxis()->SetLabelFont(42);
  xgr[mch_order]->GetXaxis()->SetLabelOffset(0.014);
  xgr[mch_order]->GetXaxis()->SetTitleOffset(1.4);
  xgr[mch_order]->GetYaxis()->SetTitleOffset(1.2);
  xgr[mch_order]->GetYaxis()->SetNdivisions(506);
  xgr[mch_order]->GetYaxis()->SetLabelFont(42);
  xgr[mch_order]->GetYaxis()->SetLabelOffset(0.01);
  xgr[mch_order]->GetXaxis()->SetTitleFont(42);
  xgr[mch_order]->GetXaxis()->SetTitle("Time (sec)");
  xgr[mch_order]->GetXaxis()->CenterTitle(true);
  xgr[mch_order]->GetYaxis()->SetTitleFont(42);
  xgr[mch_order]->GetYaxis()->SetTitle("Frequency^{-8/3}");
  xgr[mch_order]->GetYaxis()->CenterTitle(true);
  xgr[mch_order]->GetXaxis()->SetLabelSize(0.03);
  xgr[mch_order]->GetYaxis()->SetLabelSize(0.03);

  xgr[mch_order]->SetLineColor(kGreen);
  xgr[mch_order]->SetMarkerStyle(20);
  xgr[mch_order]->SetMarkerColor(color);
  xgr[mch_order]->SetMarkerSize(1);
  mch_order==0 ? xgr[mch_order]->Draw("XAP") : xgr[mch_order]->Draw("XPsame");

  TF1* fit = &(pcd->fit);
  double mchirp = pcd->mchirp;
  if(tmerger==0) tmerger = -fit->GetParameter(1)/fit->GetParameter(0);
  //if(xmin==0)    xmin = fit->GetXmin();
  //if(xmax==0)    xmax = tmerger-0.001;
  xmin = fit->GetXmin();
  xmax = tmerger-0.001;
  double p0 = GetFitParameter(mchirp);
  double p1 = tmerger;

  xfit[mch_order] = new TF1(fit->GetName(), fit->GetTitle(), xmin, xmax);
  xfit[mch_order]->SetLineColor(kRed);
  xfit[mch_order]->SetParameter(0, p0);
  xfit[mch_order]->SetParameter(1, p1);
  xgr[mch_order]->Fit(xfit[mch_order],"Q");
  xfit[mch_order]->Draw("same");

  ebbh->ch_tmerger[mch_order] = -xfit[mch_order]->GetParameter(1)/xfit[mch_order]->GetParameter(0);; 

  return;
}

void PlotMChirp(EBBH* ebbh) {

  double chi2_thr    = CHI2_THR;
  double tmerger_cut = TMERGER_CUT;
  double zmax_thr    = ZMAX_THR;

  ebbh->nc.cData[0].mchirp = 0;
  ebbh->nc.cData[0].mchirperr = 0;
  ebbh->nc.cData[0].tmrgr = 0;
  ebbh->nc.cData[0].tmrgrerr = 0;
  ebbh->nc.cData[0].chi2chirp = 0;

  PCH = new watplot(const_cast<char*>("chirp"));

  for(int i=0; i<NCHIRP_MAX; i++) {

    ebbh->nc.mchirp5(1,zmax_thr,chi2_thr,tmerger_cut);

    clusterdata* pcd = &ebbh->nc.cData[0];
    printf("mchirp : %.2e %.3f %.3f %.3f %.3f \n\n",
           pcd->mchirp, pcd->mchirperr, pcd->tmrgr,
           pcd->tmrgrerr,  pcd->chi2chirp);

    // draw chirp : f^(-8/3) vs time
    if(gPLOT) {
      TString ofile = gOFILE;
      ofile.ReplaceAll("gTYPE_6",TString::Format("gTYPE_6_%d",i+1).Data());

      PCH->canvas->cd();
      PCH->plot(pcd, 0);
      PCH->canvas->Print(ofile);
    } else {
      cout << "WARNING : Not Implemented : use gtype = -6" << endl;exit(0);
      char title[256];
      sprintf(title,"chirp mass : rec = %3.3f [%3.2f] , chi2 = %3.2f",
                  pcd->mchirp,pcd->mchirperr,pcd->chi2chirp);
      TGraphErrors* gr = &(pcd->chirp);
      xgr[i] = new TGraphErrors(gr->GetN(),gr->GetX(),gr->GetY(),gr->GetEX(),gr->GetEY());
      xgr[i]->SetTitle(title);
      xgr[i]->Draw("ALP");
      return;
      TF1* fit = &(pcd->fit);
      double tmerger = -fit->GetParameter(1)/fit->GetParameter(0);
      double xmin = fit->GetXmin()-0.1;
      double xmax = tmerger-0.001;
      xfit[i] = new TF1(fit->GetName(), "pow([0]*(x-[1]),-3./8.)", xmin, xmax);
      double p0 = GetFitParameter(pcd->mchirp);
      double p1 = tmerger;
      xfit[i]->SetLineColor(kRed);
      xfit[i]->SetParameter(0, p0);
      xfit[i]->SetParameter(1, p1);
      xfit[i]->Draw("same");
      return; 
    }
  }

  exit(0);
}

void PlotMultiChirp2(EBBH* ebbh, int ifo, TString ptype, TString wtype) { 

  double tmerger=0;
  double xmin=0;
  double xmax=0;

  // plot chirp

  cout << endl;
  cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << ebbh->name << endl;
  cout << endl << "         ";
  cout <<   "spin1x=" << ebbh->spin[0] << "\tspin1y=" << ebbh->spin[1] << "\tspin1z=" << ebbh->spin[2]  
       << "\tspin2x=" << ebbh->spin[3] << "\tspin2y=" << ebbh->spin[4] << "\tspin2z=" << ebbh->spin[2] << endl; 
  cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;

  for(int i=0;i<NCHIRP_MAX;i++) FindChirp(ebbh, i, tmerger, xmin, xmax, ifo, ptype, wtype);

  cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;

  // print final 

  int imax=0;
  float emax=0;
  // find chirp with max energy
  for(int i=0;i<NCHIRP_MAX;i++) if(ebbh->ch_energy[i]>emax) {emax=ebbh->ch_energy[i];imax=i;}
  int imax2=imax;
  // find chirp with highest chirp mass and energy > (max chirp energy)/4
  for(int i=0;i<NCHIRP_MAX;i++) {
    if(i==imax) continue;
    if(ebbh->ch_mass[i]>0) {
      if(ebbh->ch_mass[i]>ebbh->ch_mass[imax]) if(ebbh->ch_energy[i]>ebbh->ch_energy[imax]/4.) imax2=i;
    }
  }
  imax=imax2;
  float ebbh->ei=0;
  for(int i=0;i<NCHIRP_MAX;i++) {
    if(ebbh->ch_mass[i]>0) {
      if(ebbh->ch_mass[i]<ebbh->ch_mass[imax]) ebbh->ei+=ebbh->ch_energy[i];
      if(ebbh->ch_mass[i]>ebbh->ch_mass[imax]) ebbh->ei-=ebbh->ch_energy[i];
    }
  }
  ebbh->ei/=ebbh->ch_energy[imax];
  ebbh->ei*=100;

  cout << "CHIRP_M ->   " << " MCH_MASS : " << " rec :" << ebbh->ch_mass[imax] << " [" << ebbh->ch_merr[imax] << "]" 
       << " / "<< " inj : " << ebbh->mchirp << "\t\t\t\tMCH_EI : " << ebbh->ei << "\t\tMCH_TIME : " << ebbh->ch_tmerger[imax] << " (sec)" << endl;
  cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;

  return;
}

void FindChirp(EBBH* ebbh, int mch_order, double& tmerger, double& xmin, double& xmax, int ifo, TString ptype, TString wtype) { 

  if((ptype!="like")&&(ptype!="stft")&&(ptype!="")) {
    cout << "PlotMultiChirp2 Error : wrong ptype, available value are : like/stft" << endl;
    exit(1); 
  }

  double chi2_thr    = CHI2_THR;
  double tmerger_cut = TMERGER_CUT;
  double zmax_thr    = ZMAX_THR;

  EColor lcolor = (ptype=="like") ? kBlack : kWhite;
  int lwidth=3;
  int lstyle=2;

  if(mch_order==0 && ptype=="like") {
    WTS = new watplot(const_cast<char*>("wts"));
    WTS->plot(&ebbh->nc, 1, 2, 'L', 0, const_cast<char*>("COLZ"));
    WTS->hist2D->GetYaxis()->SetRangeUser(FLOW, FHIGH);
    WTS->canvas->SetLogy();
  }

  double echirp = ebbh->nc.mchirp5(1,zmax_thr,chi2_thr,tmerger_cut);

  if(mch_order==0 && ptype=="like") {
    if(gTYPE==3) PlotCluster(WTS, &ebbh->nc, FLOW, FHIGH);  
  }

  clusterdata* pcd = &ebbh->nc.cData[0];
  TF1* fit = &(pcd->fit);

  double mchirp = pcd->mchirp;

  if(tmerger==0) tmerger = -fit->GetParameter(1)/fit->GetParameter(0);
  double Xmin = fit->GetXmin()-0.1;
  double Xmax = tmerger-0.001;
  if(mch_order==0) {xmin=Xmin;xmax=Xmax;}
  if(Xmin<xmin) xmin = Xmin;
  if(Xmax>xmax) xmax = Xmax;

  double p0 = GetFitParameter(mchirp);
  double p1 = tmerger;

  ebbh->ch_mass[mch_order]    = pcd->mchirp; 
  ebbh->ch_tmerger[mch_order] = -fit->GetParameter(1)/fit->GetParameter(0); 
  ebbh->ch_merr[mch_order]    = pcd->mchirperr; 
  ebbh->ch_mchi2[mch_order]   = pcd->chi2chirp; 
  ebbh->ch_energy[mch_order]  = echirp;

  // fit function -> freq = p0 * (time-p1);
  xfit[mch_order] = new TF1(fit->GetName(), "pow([0]*(x-[1]),-3./8.)", xmin, xmax);
  xfit[mch_order]->SetLineColor(lcolor);
  xfit[mch_order]->SetLineWidth(lwidth);
  xfit[mch_order]->SetLineStyle(lstyle);
  xfit[mch_order]->SetParameter(0, p0);
  xfit[mch_order]->SetParameter(1, p1);
  //mch_order==0 ? xfit[mch_order]->SetParameter(1, p1) : xfit[mch_order]->FixParameter(1, p1);

  if(mch_order==0 && ptype=="stft") {
    if(wtype=="wdat") MDC->Draw(ebbh->wdat[ifo],MDC_TF);
    if(wtype=="winj") MDC->Draw(ebbh->winj[ifo],MDC_TF);
    if(wtype=="wrec") MDC->Draw(ebbh->wrec[ifo],MDC_TF);
    stft = MDC->GetSTFT();
    stft->GetHistogram()->SetTitle(ebbh->name);
    stft->GetHistogram()->GetYaxis()->SetRangeUser(FLOW,FHIGH);
  } 
  if(stft) stft->GetHistogram()->GetXaxis()->SetRangeUser(xmin,xmax+0.3);

  if(ptype!="") xfit[mch_order]->Draw("same");

  if(gPLOT && mch_order==NCHIRP_MAX-1 && ptype=="stft") {stft->Print(gOFILE);exit(0);}
  if(gPLOT && mch_order==NCHIRP_MAX-1 && ptype=="like") {WTS->canvas->Print(gOFILE);exit(0);}

  // format standard output
  cout.precision(1);
  cout.setf(ios::fixed, ios::floatfield);

  cout << "CHIRP_"   << mch_order+1 << " ->   "   
       << " MCH_MASS : "    << ebbh->ch_mass[mch_order] << "\t\tMCH_ERR : " << ebbh->ch_merr[mch_order] 
       << "\t\tMCH_CHI2 : " << ebbh->ch_mchi2[mch_order] << "\t\tMCH_EN : "  << ebbh->ch_energy[mch_order] 
       << "\t\tMCH_TIME : " << setprecision(2) << ebbh->ch_tmerger[mch_order] << " (sec)" << endl;

  return;
}

void FillHistogram(EBBH* ebbh, TH2F&* hN, TH2F&* hL, TH2F&* heT, TH2F&* heF, double zmax_thr) {

  netcluster* pwc = &(ebbh->nc);

  int cid = 1;
 int nifo = ebbh->ifo.size();

  double RATE = pwc->rate;                              // original rate

  std::vector<int>* vint = &(pwc->cList[cid-1]);        // pixel list

  int V = vint->size();                                 // cluster size
  if(!V) return;                                               

  double zmax = -1e100;
  double etot = 0;
  for(int j=0; j<V; j++) {
     netpixel* pix = pwc->getPixel(cid, j);
     if(pix->likelihood<1. || pix->frequency==0 ) continue;
     if(pix->likelihood>zmax) zmax = pix->likelihood;
     etot+=pix->likelihood;
  }
  double zthr = zmax_thr*zmax;

  int minLayers=1000;
  int maxLayers=0;   
  double minTime=1e20;
  double maxTime=0.;  
  double minFreq=1e20;
  double maxFreq=0.;  
  for(int j=0; j<V; j++) {                      // loop over the pixels
    netpixel* pix = pwc->getPixel(cid,j);                               
    if(!pix->core) continue;                                           

    if(pix->layers<minLayers) minLayers=pix->layers;
    if(pix->layers>maxLayers) maxLayers=pix->layers;

    double dt = 1./pix->rate;
    double time = int(pix->time/pix->layers)/double(pix->rate);         // central bin time
    time -= dt/2.;                                                      // begin bin time
    if(time<minTime) minTime=time;                   
    if(time+dt>maxTime) maxTime=time+dt;                   

    double freq = pix->frequency*pix->rate/2.; 
    if(freq<minFreq) minFreq=freq;     
    if(freq>maxFreq) maxFreq=freq;     
  }                                    

  int minRate=RATE/(maxLayers-1);
  int maxRate=RATE/(minLayers-1);

  double mindt = 1./maxRate;
  double maxdt = 1./minRate;
  double mindf = minRate/2.;
  double maxdf = maxRate/2.;

  cout.precision(1);
  cout.setf(ios::fixed, ios::floatfield);

/*
  cout << "minRate : " << minRate << "\t\t\t maxRate : " << maxRate << endl;
  cout << "minTime : " << minTime << "\t\t\t maxTime : " << maxTime << endl;
  cout << "minFreq : " << minFreq << "\t\t\t maxFreq : " << maxFreq << endl;
  cout << "mindt   : " << mindt   << "\t\t\t maxdt   : " << maxdt << endl;
  cout << "mindf   : " << mindf   << "\t\t\t maxdf   : " << maxdf << endl;
*/

  double iminTime = minTime-maxdt;
  double imaxTime = maxTime+maxdt;
  int nTime = (imaxTime-iminTime)*maxRate;

  if(hN) { delete hN; hN=NULL; }
  hN = new TH2F("hN", "hN", nTime, iminTime, imaxTime, 2*(maxLayers-1), 0, RATE/2);
  if(hL) { delete hL; hL=NULL; }
  hL = new TH2F("hL", "hL", nTime, iminTime, imaxTime, 2*(maxLayers-1), 0, RATE/2);
  if(heT) { delete heT; heT=NULL; }
  heT = new TH2F("heT", "heT", nTime, iminTime, imaxTime, 2*(maxLayers-1), 0, RATE/2);
  if(heF) { delete heF; heF=NULL; }
  heF = new TH2F("heF", "heF", nTime, iminTime, imaxTime, 2*(maxLayers-1), 0, RATE/2);

  double dFreq = (maxFreq-minFreq)/10.>2*maxdf ? (maxFreq-minFreq)/10. : 2*maxdf ;
  double mFreq = minFreq-dFreq<0 ? 0 : minFreq-dFreq;
  double MFreq = maxFreq+dFreq>RATE/2 ? RATE/2 : maxFreq+dFreq;

  double dTime = (maxTime-minTime)/10.>2*maxdt ? (maxTime-minTime)/10. : 2*maxdt ;
  double mTime = minTime-dTime<iminTime ? iminTime : minTime-dTime;
  double MTime = maxTime+dTime>imaxTime ? imaxTime : maxTime+dTime;

  int npix=0;
  double Likelihood=0;
  for(int n=0; n<V; n++) {
    netpixel* pix = pwc->getPixel(cid,n);
    if(!pix->core) continue;            

    if(pix->likelihood<zthr) continue;

    double sSNR=0;
    for(int m=0; m<nifo; m++) {                 
      sSNR += pow(pix->getdata('S',m),2);          // snr whitened reconstructed signal 00
      sSNR += pow(pix->getdata('P',m),2);          // snr whitened reconstructed signal 90
    }                                                                      

    int iRATE = int(pix->rate+0.5); 
    int M=maxRate/iRATE;              
    int K=2*(maxLayers-1)/(pix->layers-1);
    double dt = 1./pix->rate;
    double itime = int(pix->time/pix->layers)/double(pix->rate);        // central bin time
    itime -= dt/2.;                                                     // begin bin time
    int i=(itime-iminTime)*maxRate;                            
    int j=pix->frequency*K;                                    
    Likelihood+=sSNR; 
//    sSNR /= M*K;                                // Normalization (watplot not use normalization !!!)
//sSNR = pix->likelihood;                                // Use original oriinal likelihood before WP
    for(int m=0;m<M;m++) {                                     
      for(int k=0;k<K;k++) {
        double L = hL->GetBinContent(i+1+m,j+1+k-K/2);
        hL->SetBinContent(i+1+m,j+1+k-K/2,sSNR+L);

        double N = hN->GetBinContent(i+1+m,j+1+k-K/2);
        hN->SetBinContent(i+1+m,j+1+k-K/2,N+1);

        double eT = heT->GetBinContent(i+1+m,j+1+k-K/2);
        //heT->SetBinContent(i+1+m,j+1+k-K/2,eT+M*M*L);
        heT->SetBinContent(i+1+m,j+1+k-K/2,eT+M);

        double eF = heF->GetBinContent(i+1+m,j+1+k-K/2);
        //heF->SetBinContent(i+1+m,j+1+k-K/2,eF+K*K*L);
        heF->SetBinContent(i+1+m,j+1+k-K/2,eF+K);
      }                                                               
    }                                                                 
    npix++;
  }                                                                   
  cout << "Likelihood : " << Likelihood << " npix : " << npix << endl;

  for(int i=0;i<=hL->GetNbinsX();i++) {
    for(int j=0;j<=hL->GetNbinsY();j++) {

      double L = hL->GetBinContent(i,j);
      if(L<=0) continue;

      double N = hN->GetBinContent(i,j);

      double eT = heT->GetBinContent(i,j);
      //heT->SetBinContent(i,j,mindt*sqrt(eT/L));
      //heT->SetBinContent(i,j,mindt);
      heT->SetBinContent(i,j,mindt*eT/N);

      double eF = heF->GetBinContent(i,j);
      //heF->SetBinContent(i,j,(mindf/2.)*sqrt(eF/L));
      //heF->SetBinContent(i,j,(mindf/2.));
      heF->SetBinContent(i,j,(mindf/2.)*eF/N);
    }
  }

  return;
}

void PlotHistogram(EBBH* ebbh, TH2F&* h2d) {

  if(!h2d) return;

  if(!canvas) canvas = new TCanvas;
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLogy();

  h2d->SetTitle(ebbh->name);

  h2d->GetYaxis()->SetRangeUser(FLOW,FHIGH);

  h2d->SetStats(kFALSE);
  h2d->SetTitleFont(12);
  h2d->SetFillColor(kWhite);

  h2d->GetXaxis()->SetNdivisions(506);
  h2d->GetXaxis()->SetLabelFont(42);
  h2d->GetXaxis()->SetLabelOffset(0.014);
  h2d->GetXaxis()->SetTitleOffset(1.4);
  h2d->GetYaxis()->SetTitleOffset(1.2);
  h2d->GetYaxis()->SetNdivisions(506);
  h2d->GetYaxis()->SetLabelFont(42);
  h2d->GetYaxis()->SetLabelOffset(0.01);
  h2d->GetZaxis()->SetLabelFont(42);
  h2d->GetZaxis()->SetNoExponent(false);
  h2d->GetZaxis()->SetNdivisions(506);

  h2d->GetXaxis()->SetTitleFont(42);
  h2d->GetXaxis()->SetTitle("Time (sec)");
  h2d->GetXaxis()->CenterTitle(true);
  h2d->GetYaxis()->SetTitleFont(42);
  h2d->GetYaxis()->SetTitle("Frequency (Hz)");
  h2d->GetYaxis()->CenterTitle(true);

  h2d->GetZaxis()->SetTitleOffset(0.6);
  h2d->GetZaxis()->SetTitleFont(42);
  h2d->GetZaxis()->CenterTitle(true);

  h2d->GetXaxis()->SetLabelSize(0.03);
  h2d->GetYaxis()->SetLabelSize(0.03);
  h2d->GetZaxis()->SetLabelSize(0.03);

  h2d->Draw("colz");
  if(fchirp) fchirp->Draw("same");

  if(gPLOT) {canvas->Print(gOFILE);exit(0);}

  return;
}

double ChirpIntegralHistogram(EBBH* ebbh, TH2F&* h2d, double dMc) {

  if(!h2d) return;

  const double G  = watconstants::GravitationalConstant();
  const double SM = watconstants::SolarMass();
  const double C  = watconstants::SpeedOfLightInVacuo();
  const double Pi = TMath::Pi();

  double Mc = ebbh->ch_mass[0];
  //cout << "Estimated Mc : " << Mc << endl;

  Mc+=dMc;

  double p0 = 256.*Pi/5*pow(G*Mc*SM*Pi/C/C/C, 5./3);
  double p1 = ebbh->ch_tmerger[0];

  double xmin = h2d->GetXaxis()->GetXmin();
  double xmax = h2d->GetXaxis()->GetXmax();
  if(p1<xmax) xmax=p1;
  cout.precision(3);
  //cout << "xmin : " << xmin << " xmax : " << xmax << endl;

  // fit function -> freq = p0 * (time-p1);
  fchirp = new TF1("fchirp", "pow([0]*([1]-x),-3./8.)", xmin, xmax);
  fchirp->SetLineColor(kBlack);
  fchirp->SetLineWidth(3);
  fchirp->SetLineStyle(2);
  fchirp->SetParameter(0, p0);
  fchirp->SetParameter(1, p1);

  double dT = hL->GetXaxis()->GetBinWidth(0);
  double dF = hL->GetYaxis()->GetBinWidth(0);
  //cout << "dT = " << dT << " dF = " << dF << endl;

  double cenergy=0.; 
  for(int i=1;i<=hL->GetNbinsX();i++) {
    double T = hL->GetXaxis()->GetBinCenter(i)+hL->GetXaxis()->GetBinWidth(i)/2.;
    double DF = fchirp->Eval(T)-fchirp->Eval(T-dT); 
    for(int j=0;j<=hL->GetNbinsY();j++) {

      double L = hL->GetBinContent(i,j);
      if(L<=0) continue;

      double F = hL->GetYaxis()->GetBinCenter(j)+hL->GetYaxis()->GetBinWidth(j)/2.;

      //if(F-fchirp->Eval(T-dT)<DF && (F-fchirp->Eval(T-dT)>0))    hL->SetBinContent(i,j,0);
      if(F-fchirp->Eval(T-dT)<DF && (F-fchirp->Eval(T-dT)>0))    cenergy+=L;
    }
  }

  //cout << "CHIRP ENERGY : " << cenergy << endl;

  delete fchirp; fchirp=NULL;

  return cenergy;
}

void DrawMChirpInstantaneousFrequency(EBBH* ebbh, TString wtype, int ifo) {

  wavearray<double> wf; 

  if(wtype=="winj") wf = ebbh->winj[ifo];     // whiten injected waveforms
  if(wtype=="sinj") wf = ebbh->sinj[ifo];     // strain injected waveforms
  if(wtype=="wrec") wf = ebbh->wrec[ifo];     // whiten reconstructed waveforms
  if(wtype=="srec") wf = ebbh->srec[ifo];     // strain injected waveforms
  if(wtype=="wsim") wf = ebbh->wsim[0];       // simulated waveforms -> hp

  double start = wf.start();
  double imchirp = ebbh->ch_mass[0];
  //double tmerger = ebbh->ch_tmerger[0]-start;
  double tmerger = ebbh->ch_tmerger[0];
  cout << "mchirp5 - mchirp  : " << imchirp << endl;
  cout << "mchirp5 - tmerger : " << tmerger << endl;

  //wf.start(0);

  // Get hp envelope
  wavearray<double> ewf = CWB::Toolbox::getHilbertEnvelope(wf);
  MDC->SetZoom(0.999);
  MDC->Draw(wf,MDC_TIME,"ALP ZOOM",kGray);
  MDC->Draw(ewf,MDC_TIME,"same",kRed);

  watplot* pts = MDC->GetWatPlot();
  TGraph* gr = pts->getGraph(0);
  gr->SetTitle(ebbh->name);

  // Get wf instantaneous frequency
  wavearray<double> fwf = CWB::Toolbox::getHilbertIFrequency(wf);
  //MDC->Draw(fwf,MDC_TIME,"ALP ZOOM",kBlack);return;

  // F -> F^-3/8
  double sF = 128;
  for(int i=0;i<fwf.size();i++) {
    fwf[i] = fwf[i]>FLOW ? pow(fwf[i]/sF,-8./3.) : 0;
  }
  //MDC->Draw(fwf,MDC_TIME,"ALP ZOOM",kBlack);return;

  const double G  = watconstants::GravitationalConstant();
  const double SM = watconstants::SolarMass();
  const double C  = watconstants::SpeedOfLightInVacuo();
  const double Pi = TMath::Pi();
 
  double M1=ebbh->mass[0]; 
  double M2=ebbh->mass[1]; 
  double Mc = pow(M1*M2,3./5.)/pow(M1+M2,1./5.);
  double MC = Mc;
  double mc = Mc;
  cout << "Mc " << Mc << endl;
  Mc *= SM;
  double kk = 256.*Pi/5*pow(G*SM*Pi/C/C/C, 5./3);
  double K  = 256.*Pi/5*pow(G*Mc*Pi/C/C/C, 5./3);
  // scaling of frequency: in units of 128Hz
  kk *= pow(sF, 8./3);
  K  *= pow(sF, 8./3);

  wavearray<double> tt = fwf;
  wavearray<double> ff = fwf;
  tt=0;
  ff=0;

  double dt = 1./ff.rate();
  double Tc = tmerger;
  for(int i=0;i<ff.size();i++) {
    double T = ff.start()+i*dt;
    //double F = pow(Tc-K*T,-3./8.);
    double F = T<Tc ? pow(K*(Tc-T),-3./8.) : 0;
    //double F = pow(K*(T),-3./8.);
    if(F<0 || F>1000) F=0; 
    F = F>0 ? pow(F,-8./3.) : 0;
    tt.data[i] = T;
    ff.data[i] = F;
    if(F>1) F=0; 
  }
  //MDC->Draw(ff,MDC_TIME,"ALP",kBlack);return;

  // remove data where the waveform amplitude (ewf) is negligeable
  double emax=ewf.max();
  int SS = 0;
  for(int i=0;i<fwf.size();i++) {
    if(ewf[i]/emax<0.1) SS++; else break;
  }
  int EE=0;
  for(int i=fwf.size()-1;i>=0;i--) {
    if(ewf[i]/emax<0.2) EE++; else break;
  }

  // ff errors
  wavearray<double> eff = ff;
  //eff=1.e20;
  //for(int i=0;i<fwf.size();i++) eff[i] = ewf[i]/emax>0.01 ? 1./ewf[i] : 1.e20;
  eff=0;

  TGraphErrors *grr = new TGraphErrors(fwf.size()-SS-EE, tt.data+SS, fwf.data+SS, 0, eff.data+SS);
  grr->GetHistogram()->SetTitle(ebbh->name);  
  grr->GetHistogram()->GetXaxis()->SetTitle("time (sec)");  
  grr->GetHistogram()->GetYaxis()->SetTitle("freq^-3/8 ");  
  grr->GetHistogram()->GetXaxis()->SetTitleOffset(1.1);  
  grr->GetHistogram()->GetYaxis()->SetTitleOffset(1.1);  
  grr->SetMarkerStyle(7);
  grr->SetMarkerColor(kGray);
  grr->SetLineColor(kGray);

  TF1 *func = new TF1("func", "[0]*(x-[1])", 0.05,Tc);
  //func->FixParameter(0,-kk*pow(MC,1./0.6));
  //func->FixParameter(0,-kk*pow(imchirp,1./0.6));
  //func->FixParameter(1,Tc);
  func->SetParameter(0,-kk*pow(imchirp,1./0.6));
  func->SetParameter(1,Tc);
  func->SetLineColor(kRed);
  TCanvas *canvas = new TCanvas("canvas", "Linear and robust linear fitting");
  canvas->SetGrid();
  grr->Draw("ALP");
  printf("Ordinary least squares:\n");
  grr->Fit(func);
  func->Draw("same");

  double mch = -func->GetParameter(0)/kk;
  double mchirp = mch>0 ? pow(mch, 0.6) : -pow(-mch, 0.6);

  cout << "mchirp : " << mchirp << endl;
  cout << "mchirp-MC : " << mchirp-MC << endl;
  cout << "m1 m2 ichirp ochirp ochirp-ichirp : " << M1 << " " << M2 << " " << MC << " " << " " << mchirp << " " << mchirp-MC << endl;


  double mm = (pow(80.,-8./3.)-pow(30.,-8./3.))/(0.13);
  double mch1 = -mm/kk;
  double mchirp1 = mch1>0 ? pow(mch1, 0.6) : -pow(-mch1, 0.6);
  cout << "mchirp1 : " << mchirp1 << endl;

  return;
}

void GetCBC(EBBH* ebbh) {

  double To = 931072130;

  double Tc=ebbh->ch_tmerger[0]; 

  double M1 = ebbh->mass[0]<ebbh->mass[1] ? ebbh->mass[0] : ebbh->mass[1]; 
  double M2 = ebbh->mass[0]>ebbh->mass[1] ? ebbh->mass[0] : ebbh->mass[1]; 

  double S1 = fabs(ebbh->spin[0])<fabs(ebbh->spin[1]) ? fabs(ebbh->spin[0]) : fabs(ebbh->spin[1]); 
  double S2 = fabs(ebbh->spin[0])>fabs(ebbh->spin[1]) ? fabs(ebbh->spin[0]) : fabs(ebbh->spin[1]); 

  double DIST=ebbh->dist; 

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------
  TString inspOptions="";
  inspOptions+= TString::Format("--gps-start-time %f --gps-end-time %f ",Tc+To,Tc+To);
  inspOptions+= "--taper-injection start --seed 962488 ";
  inspOptions+= "--dir ./ ";
  inspOptions+= "--f-lower 24.000000 ";
  inspOptions+= "--time-step 60 ";

  inspOptions+= "--waveform IMRPhenomBthreePointFivePN ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--i-distr uniform ";

  inspOptions+= TString::Format("--min-mass1 %f --max-mass1 %f ",M1,M1);
  inspOptions+= TString::Format("--min-mass2 %f --max-mass2 %f ",M2,M2);
  inspOptions+= TString::Format("--min-mtotal %f --max-mtotal %f ",M1+M1,M2+M2);

  inspOptions+= "--m-distr componentMass ";
  inspOptions+= "--enable-spin --aligned ";

  inspOptions+= TString::Format("--min-spin1 %f --max-spin1 %f ",S1,S1);
  inspOptions+= TString::Format("--min-spin2 %f --max-spin2 %f ",S2,S2);

  inspOptions+= "--amp-order 0 ";
  inspOptions+= "--d-distr volume ";
  inspOptions+= TString::Format("--min-distance %f --max-distance %f ",DIST,DIST);
  //inspOptions+= "--output inspirals.xml ";              // set output xml file

  MDC->SetInspiral("IMRPhenomBthreePointFivePN",inspOptions.Data());

  // Get the first waveform hp,hx components starting from gps = 931072130
  ebbh->wsim[0] = MDC->GetInspiral("hp",Tc+To-200,Tc+To+200);
  ebbh->wsim[1] = MDC->GetInspiral("hx",Tc-200,Tc+200);
  ebbh->wsim[0].start(ebbh->wsim[0].start()-To);
  ebbh->wsim[0].resample(2*FHIGH);
  ebbh->wsim[1].start(ebbh->wsim[1].start()-To);
  ebbh->wsim[1].resample(2*FHIGH);
  cout << "size : " << ebbh->wsim[0].size() << " rate : " << ebbh->wsim[0].rate() << " start : " << (int)ebbh->wsim[0].start() << endl;

  return;
}

void GetEBBH(EBBH* ebbh) {

  double To = ebbh->wrec[0].start();

  double Tc=ebbh->ch_tmerger[0]; 

  double M1 = ebbh->mass[0]<ebbh->mass[1] ? ebbh->mass[0] : ebbh->mass[1]; 
  double M2 = ebbh->mass[0]>ebbh->mass[1] ? ebbh->mass[0] : ebbh->mass[1]; 

  double S1 = fabs(ebbh->spin[0])<fabs(ebbh->spin[1]) ? fabs(ebbh->spin[0]) : fabs(ebbh->spin[1]); 
  double S2 = fabs(ebbh->spin[0])>fabs(ebbh->spin[1]) ? fabs(ebbh->spin[0]) : fabs(ebbh->spin[1]); 

  double DIST=ebbh->dist; 

  double E0  = ebbh->e0;
  double RP0 = ebbh->rp0>00 ? ebbh->rp0 : 14;
  double REDSHIFT = ebbh->redshift;

  // get ebbh parameters from input command line
  TString ebbh_parms = TString::Format("1 %f %f %f %f %f %f",M1,M2,RP0,E0,DIST,REDSHIFT);
  cout << "ebbh_parms : " << ebbh_parms << endl;

  int gps = To+Tc;

  vector<mdcpar> par;
  par.resize(1);
  par[0].name=ebbh_parms;
  MDC->AddWaveform(MDC_EBBH, par);
  MDC->Print(0);

  // Get the first waveform hp,hx components starting from gps = 931072130

  waveform wf = MDC->GetWaveform("eBBH");
  ebbh->wsim[0] = wf.hp;
  ebbh->wsim[1] = wf.hx;
  ebbh->wsim[0].start(ebbh->wrec[0].start()-ebbh->seggps);
  ebbh->wsim[0].resample(2*FHIGH);
  ebbh->wsim[1].start(ebbh->wrec[1].start()-ebbh->seggps);
  ebbh->wsim[1].resample(2*FHIGH);
  cout << "size : " << ebbh->wsim[0].size() << " rate : " << ebbh->wsim[0].rate() << " start : " << (int)ebbh->wsim[0].start() << endl;

  return;
}
