// this macro shows how to read the PE results from ROOT and ASCII file
// plot the waveforms with tghe same format presented in the PE CED
// Author : G.Vedovato

#include "constants.hh"

void ReadDataFromASCII(TString ipath, int ifo, wavearray<double>* wrec, wavearray<double>* wmed, 
                                wavearray<double>* wl50, wavearray<double>* wu50,
                                wavearray<double>* wl90, wavearray<double>* wu90,
                                wavearray<double>* frec, wavearray<double>* fmed, 
                                wavearray<double>* fl50, wavearray<double>* fu50,
                                wavearray<double>* fl90, wavearray<double>* fu90);

std::vector<TString> ReadDataFromROOT(TString ifile, int ifo, wavearray<double> *winj, wavearray<double> *wrec, wavearray<double> *wwht,
                               wavearray<double> *wmed, wavearray<double> *wl50, wavearray<double> *wu50, 
                               wavearray<double> *wl90, wavearray<double> *wu90);

void PlotWaveformAsymmErrors(  TString ofname, TString title, wavearray<double>* wrec,
                               wavearray<double>* wmed, wavearray<double>* wl50, wavearray<double>* wu50,
                               wavearray<double>* wl90, wavearray<double>* wu90, wavearray<double>* wref, 
                               TString pdir, double P, bool freq=false, bool showerr=true);

void DumpWaveform(int ifo,     wavearray<double>* wrec, wavearray<double>* wmed, 
                               wavearray<double>* wl50, wavearray<double>* wu50,
                               wavearray<double>* wl90, wavearray<double>* wu90);


double GetTimeBoundaries(wavearray<double> x, double P, double& bT, double& eT);
wavearray<double> GetAlignedWaveform(wavearray<double>* wf, wavearray<double>* wref);
wavearray<double> GetDifWaveform(wavearray<double>* wf1, wavearray<double>* wf2);

void FrequencyCut(wavearray<double>* x, double bF, double eF);
void PlotSpectrogram(wavearray<double>* x, double tstart=0, double tstop=0, TString title="", TString ofname="", double tchirp=0., double mchirp=0.);


#define IN_CWB_ASCII_FILE_L1  "L1_pe_wave.dat"
#define IN_CWB_ASCII_FILE_H1  "H1_pe_wave.dat"

#define OUT_CWB_ASCII_FILE_L1 "OUT_L1_pe_wave.dat"
#define OUT_CWB_ASCII_FILE_H1 "OUT_H1_pe_wave.dat"


//#define FREQUENCY_CUT
#define FLOW	16
#define FHIGH	256

CWB::STFT* stft;
TF1* fchirp;

void DrawWaveformPE(TString ipath, int gtype=0, int ifo=0, double tshift=0, TString label="PLOT", double P=0.99) {
//
// ipath- > if(gtype<0) input ced idir 
//          if(gtype>0) input root file name
//
// gtype -> plot type
//        < 0 -> plot instantaneous frequency : must be used only when ifo[0]=L1 and ifo[1]=H1
//          0 -> plot reconstructed vs whitened
//          1 -> plot reconstructed vs median
//          2 -> plot rec vs inj (only if injection is present in the input file)
//          3 -> plot the difference between med/rec and inj waveforms
//          4 -> plot spectrogram ifo
//          5 -> plot spectrogram ifo_0 + ifo_1(inv+tshift) : must be used only when ifo[0]=L1 and ifo[1]=H1
//          6 -> plot spectrogram ifo_0 - ifo_1(inv+tshift) : must be used only when ifo[0]=L1 and ifo[1]=H1
//        >=7 -> dump data to ascii file
//
// ifo   -> detector id : 0,1,...
//
// tshift-> the second detector is shifted by tshift and sign inverted : must be used only when ifo[0]=L1 and ifo[1]=H1
//          GW150914  -> tshift =-0.0069
//          GW151226  -> tshift =-0.0011
//          LVT151012 -> tshift = 0.001
//
// label -> used for title and output plot files
//
// P     -> [0,1] is used to restrict the plot time interval around the wave energy percentage P
//

  bool bexit = label.Contains("EXIT") ? true : false;
  label.ReplaceAll("EXIT","");

  fchirp=NULL;
  stft=NULL;

  double tchirp=0.;
  double mchirp=0.;

  if(label=="GW150914")  {tshift=-0.0069;tchirp=8.42;mchirp=30.0;}
  if(label=="GW151226")  {tshift=-0.0011;tchirp=7.65;mchirp=8.9;}
  if(label=="LVT151012") {tshift= 0.0010;tchirp=7.43;mchirp=24.0;}	// ?

  wavearray<double> winj[NIFO_MAX];

  wavearray<double> wrec[NIFO_MAX];
  wavearray<double> wwht[NIFO_MAX];
  wavearray<double> wmed[NIFO_MAX];
  wavearray<double> wl50[NIFO_MAX];
  wavearray<double> wu50[NIFO_MAX];
  wavearray<double> wl90[NIFO_MAX];
  wavearray<double> wu90[NIFO_MAX];

  wavearray<double> frec[NIFO_MAX];
  wavearray<double> fmed[NIFO_MAX];
  wavearray<double> fl50[NIFO_MAX];
  wavearray<double> fu50[NIFO_MAX];
  wavearray<double> fl90[NIFO_MAX];
  wavearray<double> fu90[NIFO_MAX];

  int nifo=0;
  std::vector<TString> ifoname;

  if(gtype<0) {
    // read data from PE ascii output files
    ReadDataFromASCII(ipath, ifo, &wrec[ifo], &wmed[ifo], &wl50[ifo], &wu50[ifo], &wl90[ifo], &wu90[ifo], 
                      &frec[ifo], &fmed[ifo], &fl50[ifo], &fu50[ifo], &fl90[ifo], &fu90[ifo]);
    PlotWaveformAsymmErrors("", "", &frec[ifo], &fmed[ifo], &fl50[ifo], &fu50[ifo], &fl90[ifo], 
                      &fu90[ifo], &wrec[ifo], ".", P, true); 
    return;
  } else {
    // read data from PE root output file
    for(int n=0;n<NIFO_MAX;n++) {
      ifoname = ReadDataFromROOT(ipath, n, &winj[n], &wrec[n], &wwht[n], &wmed[n], &wl50[n], &wu50[n], &wl90[n], &wu90[n]);
      nifo = ifoname.size();
      if(n>=nifo-1) break;
    }
    if(ifo>nifo-1) {cout << "DrawWaveformPE Error : input ifo=" << ifo << " parameter must be < nifo=" << nifo << endl;exit(1);}
#ifdef FREQUENCY_CUT
    // frequency cut for whitened data
    for(int n=0;n<nifo;n++) FrequencyCut(&wwht[n], FLOW, FHIGH);
#endif
  }

  // compute med,start,stop times for input wavearray
  double tmed   = wrec[0].start()+(wrec[0].size()/wrec[0].rate())/2.;
  double tstart = tmed-1;
  double tstop  = tmed+1;
  GetTimeBoundaries(wrec[0], P, tstart, tstop);
  //cout << "Time Boundaries : " << tstart-wrec[0].start() << " " << tstop-wrec[0].start() << endl;exit(0);

  if(gtype==0) {
    // plot reconstructed vs whitened
    PlotWaveformAsymmErrors("", "", &wrec[ifo], &wwht[ifo], &wl50[ifo], &wu50[ifo], &wl90[ifo], &wu90[ifo], &wrec[ifo], "", P, false, false);
    return;
  }
  if(gtype==1) {
    // plot reconstructed vs median
    PlotWaveformAsymmErrors("", "", &wrec[ifo], &wmed[ifo], &wl50[ifo], &wu50[ifo], &wl90[ifo], &wu90[ifo], &wrec[ifo], "", P);
    return;
  }
  if(gtype==2) {
    // plot reconstructed vs injected
    PlotWaveformAsymmErrors("", "", &winj[ifo], &wmed[ifo], &wl50[ifo], &wu50[ifo], &wl90[ifo], &wu90[ifo], &wrec[ifo], "", P);
    return;
  }
  if(gtype==3) {  	
    // plot difference between med/rec and inj
    wavearray<double> wdif[NIFO_MAX];
    wdif[ifo] = GetDifWaveform(&wmed[ifo], &winj[ifo]);
    PlotWaveformAsymmErrors("", "", &wdif[ifo], &wdif[ifo], &wl50[ifo], &wu50[ifo], &wl90[ifo], &wu90[ifo], &wrec[ifo], "", P);
    return;
  } 
  if(gtype==4) {  	
    // plot spectrogram ifo
    TString title = TString::Format("%s whitened data %s",label.Data(),ifoname[ifo].Data());
    TString ofile = TString::Format("%s_whitened_data_%s_Spectrogram.png",label.Data(),ifoname[ifo].Data());
    PlotSpectrogram(&wwht[ifo], tstart, tstop,title,ofile,tchirp,mchirp);
    if(bexit) exit(0); else return;
  } 
  if(gtype==5) {  	
    // plot spectrogram ifo_0 + ifo_1(inv+tshift)
    TString title = TString::Format("%s whitened data %s + %s(inv/tshift)",label.Data(),ifoname[0].Data(),ifoname[1].Data());
    TString ofile = TString::Format("%s_whitened_data_%s_plus_%s_Spectrogram.png",label.Data(),ifoname[0].Data(),ifoname[1].Data());
    // invert sign of the second ifo
    wwht[1]*=-1;
    // shift time of the second ifo
    CWB::mdc::TimeShift(wwht[1],tshift);
    // add first and second ifo
    wwht[0]+=wwht[1];
    PlotSpectrogram(&wwht[ifo], tstart, tstop,title,ofile,tchirp,mchirp);
    if(bexit) exit(0); else return;
  } 
  if(gtype==6) {  	
    // plot spectrogram ifo_0 - ifo_1(inv+tshift)
    TString title = TString::Format("%s whitened data %s - %s(inv/tshift)",label.Data(),ifoname[0].Data(),ifoname[1].Data());
    TString ofile = TString::Format("%s_whitened_data_%s_minus_%s_Spectrogram.png",label.Data(),ifoname[0].Data(),ifoname[1].Data());
    // invert sign of the second ifo
    wwht[1]*=-1;
    // shift time of the second ifo
    CWB::mdc::TimeShift(wwht[1],tshift);
    // subctract second from the first ifo
    wwht[0]-=wwht[1];
    PlotSpectrogram(&wwht[ifo], tstart, tstop,title,ofile,tchirp,mchirp);
    if(bexit) exit(0); else return;
  } 
  if(gtype>=7) {  	
    // dump data to ascii file
    DumpWaveform(ifo, &wrec[ifo], &wmed[ifo], &wl50[ifo], &wu50[ifo], &wl90[ifo], &wu90[ifo]);
    exit(0);
  }

  return;
}

std::vector<TString> ReadDataFromROOT(TString ifile, int ifo, wavearray<double> *winj, wavearray<double> *wrec, wavearray<double> *wwht,
                                       wavearray<double> *wmed, wavearray<double> *wl50, wavearray<double> *wu50, 
                                       wavearray<double> *wl90, wavearray<double> *wu90) {
//
// ----------------------------------------------------
// ROOT Output PE Parameters
// ----------------------------------------------------

  float pe_erR[11];                       // probability distribution of residuals
  float pe_erF[11];                       // probability distribution of frequency residuals
  int   pe_trials;                           // number of effective trials
  float pe_pnul[2*NIFO_MAX];              // null pixel statistic, for each detector pnul[0]=avr, pnul=rms
  float pe_snet[2];                       // SNRnet statistic, 0 -> avr, 1 -> rms  
  float pe_ff[2];                         // Fitting Factor statistic, 0 -> avr, 1 -> rms  
  float pe_of[2];                         // Overlap Factor statistic, 0 -> avr, 1 -> rms  
  float pe_mch[2];                        // chirp mass statistic, 0 -> avr, 1 -> rms  

  wavearray<double>* pe_wINJ[NIFO_MAX];
  wavearray<double>* pe_wREC[NIFO_MAX];
  wavearray<double>* pe_wWHT[NIFO_MAX];
  wavearray<double>* pe_wMED[NIFO_MAX];
  wavearray<double>* pe_wL50[NIFO_MAX];
  wavearray<double>* pe_wU50[NIFO_MAX];
  wavearray<double>* pe_wL90[NIFO_MAX];
  wavearray<double>* pe_wU90[NIFO_MAX];
  wavearray<double>* pe_wAVR[NIFO_MAX];
  wavearray<double>* pe_wRMS[NIFO_MAX];

  TFile* froot = new TFile(ifile,"READ");
  if(froot==NULL) {cout << "ReadWaveformPE Error : opening input root file" << endl;exit(1);}
  TTree* gTREE = (TTree*)froot->Get("waveburst");
  if(gTREE==NULL) {cout << "ReadWaveformPE Error : no waveburst present in the file" << endl;exit(1);}

  // get detector list
  TList* list = gTREE->GetUserInfo();
  int nIFO=list->GetSize();
  std::vector<TString> ifoname;
  if (nIFO==0) {cout << "ReadWaveformPE Error : no ifo present in the tree" << endl;exit(1);}
  for (int n=0;n<list->GetSize();n++) {
    detector* pDetector;
    pDetector = (detector*)list->At(n);
    ifoname.push_back(pDetector->Name);
    detectorParams dParams = pDetector->getDetectorParams();
    //pDetector->print();                                                       
  }

  for(int n=0;n<nIFO;n++) {
    pe_wINJ[n] = new wavearray<double>;
    pe_wREC[n] = new wavearray<double>;
    pe_wWHT[n] = new wavearray<double>;
    pe_wMED[n] = new wavearray<double>;
    pe_wL50[n] = new wavearray<double>;
    pe_wU50[n] = new wavearray<double>;
    pe_wL90[n] = new wavearray<double>;
    pe_wU90[n] = new wavearray<double>;
    pe_wAVR[n] = new wavearray<double>;
    pe_wRMS[n] = new wavearray<double>;
  }

  gTREE->SetBranchAddress("pe_trials",&pe_trials);
  gTREE->SetBranchAddress("pe_erR",pe_erR);
  gTREE->SetBranchAddress("pe_erF",pe_erF);
  gTREE->SetBranchAddress("pe_pnul",pe_pnul);
  gTREE->SetBranchAddress("pe_snet",pe_snet);
  gTREE->SetBranchAddress("pe_ff",pe_ff);
  gTREE->SetBranchAddress("pe_of",pe_of);
  gTREE->SetBranchAddress("pe_mch",pe_mch);

  for(int n=0;n<nIFO;n++) {
    gTREE->SetBranchAddress(TString::Format("pe_wINJ_%d",n).Data(),&pe_wINJ[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wREC_%d",n).Data(),&pe_wREC[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wWHT_%d",n).Data(),&pe_wWHT[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wMED_%d",n).Data(),&pe_wMED[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wL50_%d",n).Data(),&pe_wL50[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wU50_%d",n).Data(),&pe_wU50[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wL90_%d",n).Data(),&pe_wL90[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wU90_%d",n).Data(),&pe_wU90[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wAVR_%d",n).Data(),&pe_wAVR[n]);
    gTREE->SetBranchAddress(TString::Format("pe_wRMS_%d",n).Data(),&pe_wRMS[n]);
  }

  gTREE->GetEntry(0);

  cout.precision(4); 
  cout << "---------------------> size " << pe_wMED[0]->size() << endl;
  cout << "---------------------> pe_trials " << pe_trials << endl;

  int n = ifo;

  *winj = GetAlignedWaveform(pe_wINJ[n],pe_wREC[n]);
  *wrec = GetAlignedWaveform(pe_wREC[n],pe_wREC[n]);
  *wwht = GetAlignedWaveform(pe_wWHT[n],pe_wREC[n]);
  *wmed = GetAlignedWaveform(pe_wMED[n],pe_wREC[n]);
  *wl50 = GetAlignedWaveform(pe_wL50[n],pe_wREC[n]);
  *wu50 = GetAlignedWaveform(pe_wU50[n],pe_wREC[n]);
  *wl90 = GetAlignedWaveform(pe_wL90[n],pe_wREC[n]);
  *wu90 = GetAlignedWaveform(pe_wU90[n],pe_wREC[n]);

  for(int n=0;n<nIFO;n++) {
    delete pe_wINJ[n];
    delete pe_wREC[n];
    delete pe_wWHT[n];
    delete pe_wMED[n];
    delete pe_wL50[n];
    delete pe_wU50[n];
    delete pe_wL90[n];
    delete pe_wU90[n];
    delete pe_wAVR[n];
    delete pe_wRMS[n];
  }

  return ifoname;
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

void PlotWaveformAsymmErrors(TString ofname, TString title, wavearray<double>* wrec,
                             wavearray<double>* wmed, wavearray<double>* wl50, wavearray<double>* wu50,
                             wavearray<double>* wl90, wavearray<double>* wu90, wavearray<double>* wref, 
                             TString pdir, double P, bool freq, bool showerr) {

  int size = wrec->size();

  wavearray<double> time(size);
  wavearray<double> etime(size); etime=0;
  for (int i=0; i<size; i++) time[i] = i/wrec->rate();

  double bT, eT;
  GetTimeBoundaries(*wref, P, bT, eT);
  bT-=wref->start();
  eT-=wref->start();

  // set to 0 the frequency values outside the time range -> fix the y scale autoscale
  // info : this procedure modify the frequency input data but it is not relevant
  if(freq) {
    for(int i=0;i<wrec->size();i++) {
      if(time[i]>bT && time[i]<eT) continue;
      wrec->data[i]=0; wmed->data[i]=0; wl50->data[i]=0; wu50->data[i]=0; wl90->data[i]=0; wu90->data[i]=0;
    }
  }

  TGraphAsymmErrors* egr90 = new TGraphAsymmErrors(size,time.data,wmed->data,etime.data,etime.data,wl90->data,wu90->data);
  egr90->SetLineColor(17);
  egr90->SetFillStyle(1001);
  egr90->SetFillColor(17);
  egr90->GetXaxis()->SetTitle("time (s)");
  if(freq) egr90->GetYaxis()->SetTitle("frequency (hz)");
  else     egr90->GetYaxis()->SetTitle("magnitude");
  egr90->SetTitle(title);

  TGraphAsymmErrors* egr50 = new TGraphAsymmErrors(size,time.data,wmed->data,etime.data,etime.data,wl50->data,wu50->data);
  egr50->SetLineColor(15);
  egr50->SetFillStyle(1001);
  egr50->SetFillColor(15);

  TGraph* agr = new TGraph(size,time.data,wmed->data);
  agr->SetLineWidth(1);
  agr->SetLineColor(kWhite);
  agr->SetLineStyle(1);

  TGraph* gr = new TGraph(size,time.data,wrec->data);
  gr->SetLineWidth(1);
  gr->SetLineColor(2);

  TCanvas* canvas = new TCanvas("wrec", "wrec",200,20,800,500);
  canvas->cd();
  canvas->SetGridx();
  canvas->SetGridy();

  egr90->GetXaxis()->SetRangeUser(bT, eT);
  if(showerr) {
    egr90->Draw("A4");
    egr50->Draw("4same");
    agr->Draw("Lsame");
    gr->Draw("Lsame");
  } else {
    agr->GetXaxis()->SetTitle("time (s)");
    if(freq) agr->GetYaxis()->SetTitle("frequency (hz)");
    else     agr->GetYaxis()->SetTitle("magnitude");
    agr->SetTitle(title);
    agr->GetXaxis()->SetRangeUser(bT, eT);
    agr->SetLineColor(16);
    gr->SetLineWidth(2);
    agr->Draw("AL");
    gr->Draw("Lsame");
  }

  if(ofname!="") {
    ofname = TString(pdir)+TString("/")+ofname;
    canvas->Print(ofname);
    cout << "write : " << ofname << endl;
  }

/*
  delete canvas;
  delete egr50;
  delete egr90;
  delete agr;
  delete gr;
*/
  return;
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

void PlotWaveformAsymmErrors(TString ofname, TString title, CWB::config* cfg, wavearray<double>* wrec,
                             wavearray<double>* wmed, wavearray<double>* wl50, wavearray<double>* wu50,
                             wavearray<double>* wl90, wavearray<double>* wu90, wavearray<double>* wref, TString pdir, double P, bool freq) {

  int size = wrec->size();

  wavearray<double> time(size);
  wavearray<double> etime(size); etime=0;
  for (int i=0; i<size; i++) time[i] = i/wrec->rate();

  double bT, eT;
  GetTimeBoundaries(*wref, P, bT, eT);
  bT-=wref->start();
  eT-=wref->start();

  // set to 0 the frequency values outside the time range -> fix the y scale autoscale
  // info : this procedure modify the frequency input data but it is not relevant
  if(freq) {
    for(int i=0;i<wrec->size();i++) {
      if(time[i]>bT && time[i]<eT) continue;
      wrec->data[i]=0; wmed->data[i]=0; wl50->data[i]=0; wu50->data[i]=0; wl90->data[i]=0; wu90->data[i]=0;
    }
  }

  TGraphAsymmErrors* egr90 = new TGraphAsymmErrors(size,time.data,wmed->data,etime.data,etime.data,wl90->data,wu90->data);
  egr90->SetLineColor(17);
  egr90->SetFillStyle(1001);
  egr90->SetFillColor(17);
  egr90->GetXaxis()->SetTitle("time (s)");
  if(freq) egr90->GetYaxis()->SetTitle("frequency (hz)");
  else     egr90->GetYaxis()->SetTitle("magnitude");
  egr90->SetTitle(title);

  TGraphAsymmErrors* egr50 = new TGraphAsymmErrors(size,time.data,wmed->data,etime.data,etime.data,wl50->data,wu50->data);
  egr50->SetLineColor(15);
  egr50->SetFillStyle(1001);
  egr50->SetFillColor(15);

  TGraph* agr = new TGraph(size,time.data,wmed->data);
  agr->SetLineWidth(1);
  agr->SetLineColor(kWhite);
  agr->SetLineStyle(1);

  TGraph* gr = new TGraph(size,time.data,wrec->data);
  gr->SetLineWidth(1);
  gr->SetLineColor(2);

  TCanvas* canvas = new TCanvas("wrec", "wrec",200,20,800,500);
  canvas->cd();
  canvas->SetGridx();
  canvas->SetGridy();

  egr90->GetXaxis()->SetRangeUser(bT, eT);
  egr90->Draw("A4");
  egr50->Draw("4same");
  agr->Draw("Lsame");
  gr->Draw("Lsame");

  ofname = TString(pdir)+TString("/")+ofname;
  if(ofname!="") {
    canvas->Print(ofname);
    cout << "write : " << ofname << endl;
  }

  delete canvas;
  delete egr50;
  delete egr90;
  delete agr;
  delete gr;
}

// Dumps reconstructed waveform/time/errors array in ASCII format.
void DumpWaveform(int ifo, wavearray<double>* wrec, wavearray<double>* wmed, 
                           wavearray<double>* wl50, wavearray<double>* wu50,
                           wavearray<double>* wl90, wavearray<double>* wu90) {

  char ofname[256];

  if(ifo==0) sprintf(ofname,OUT_CWB_ASCII_FILE_L1);
  if(ifo==1) sprintf(ofname,OUT_CWB_ASCII_FILE_H1);

  ofstream out;
  out.open(ofname,ios::out);
  if (!out.good()) {cout << "Error Opening Output File : " << ofname << endl;exit(1);}
  cout << "Create Output File : " << ofname << endl;
  out.precision(19);

  // write header
  out << "#whitened data : time, amp_point, amp_median, amp_lower_50_perc, amp_lower_90_perc, amp_upper_50_perc, amp_upper_90_perc" << endl;

  // write data
  int size = wrec->size();
  double dt=1./wrec->rate();
  for (int i=0; i<size; i++) {
    double time = i*dt+wrec->start();

    double vrec = wrec->data[i];
    double vmed = wmed->data[i];
    double vl50 = fabs(wl50->data[i]);
    double vu50 = fabs(wu50->data[i]);
    double vl90 = fabs(wl90->data[i]);
    double vu90 = fabs(wu90->data[i]);

    double l50 = vmed-vl50;
    double u50 = vmed+vu50;
    double l90 = vmed-vl90;
    double u90 = vmed+vu90;

    // dump full percentiles
    out << time << " " << vrec << " " << vmed << " " << l50 << " " << l90 << " " << u50 << " " << u90 << endl;
    // dump dif percentiles from median
    //out << time << " " << vrec << " " << vmed << " " << vl50 << " " << vl90 << " " << vu50 << " " << vu90 << endl;

  }

  out.close();
}

// Read reconstructed waveform/time/errors array from ascii file
void ReadDataFromASCII(TString ipath, int ifo, wavearray<double>* wrec, wavearray<double>* wmed, 
                                wavearray<double>* wl50, wavearray<double>* wu50,
                                wavearray<double>* wl90, wavearray<double>* wu90,
                                wavearray<double>* frec, wavearray<double>* fmed, 
                                wavearray<double>* fl50, wavearray<double>* fu50,
                                wavearray<double>* fl90, wavearray<double>* fu90) {

// #whitened data : time, amp_point, amp_mean, amp_rms, amp_median, amp_lower_50_perc, amp_lower_90_perc, amp_upper_50_perc, amp_upper_90_perc, 
//                        frq_point, frq_mean, frq_rms, frq_median, frq_lower_50_perc, frq_lower_90_perc, frq_upper_50_perc, frq_upper_90_perc


  char ifname[256];

  if(ifo==0) sprintf(ifname,"%s/%s",ipath.Data(),IN_CWB_ASCII_FILE_L1);
  if(ifo==1) sprintf(ifname,"%s/%s",ipath.Data(),IN_CWB_ASCII_FILE_H1);

  ifstream in(ifname);
  if(!in.good()) {cout << "Error Opening File : " << ifname << endl;exit(1);}

  int size=0;
  char str[1024];
  int fpos=0;
  in.getline(str,1024); // skip first line : header
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    size++;
  }
  cout << "size " << size << endl;
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);
  if (size==0) {cout << "Error : File " << ifname<< " is empty" << endl;gSystem->Exit(1);}

  double time;
  double vrec, vavr, vrms, vmed, vl50, vl90, vu50, vu90;  
  double grec, gavr, grms, gmed, gl50, gl90, gu50, gu90;  

  wrec->rate(2048.);
  wmed->rate(2048.);
  wl50->rate(2048.);
  wu50->rate(2048.);
  wl90->rate(2048.);
  wu90->rate(2048.);

  wrec->resize(size);
  wmed->resize(size);
  wl50->resize(size);
  wu50->resize(size);
  wl90->resize(size);
  wu90->resize(size);

  frec->rate(2048.);
  fmed->rate(2048.);
  fl50->rate(2048.);
  fu50->rate(2048.);
  fl90->rate(2048.);
  fu90->rate(2048.);

  frec->resize(size);
  fmed->resize(size);
  fl50->resize(size);
  fu50->resize(size);
  fl90->resize(size);
  fu90->resize(size);

  // skip header
  in.getline(str,1024); 

  int i=0;
  while(1) {

    in >> time 
       >> vrec >> vavr >> vrms >> vmed >> vl50 >> vl90 >> vu50 >> vu90
       >> grec >> gavr >> grms >> gmed >> gl50 >> gl90 >> gu50 >> gu90;  

    if (!in.good()) break;

    wrec->data[i] = vrec;
    wmed->data[i] = vmed;

    // read full percentile -> reduce to relative percentiles
    wl50->data[i] = vmed-vl50;
    wu50->data[i] = vu50-vmed;
    wl90->data[i] = vmed-vl90;
    wu90->data[i] = vu90-vmed;

    if(i==0) {
      wrec->start(time);
      wmed->start(time);
      wl50->start(time);
      wu50->start(time);
      wl90->start(time);
      wu90->start(time);
    }

    frec->data[i] = grec;
    fmed->data[i] = gmed;

    // read full percentile -> reduce to relative percentiles
    fl50->data[i] = gmed-gl50;
    fu50->data[i] = gu50-gmed;
    fl90->data[i] = gmed-gl90;
    fu90->data[i] = gu90-gmed;

    if(i==0) {
      frec->start(time);
      fmed->start(time);
      fl50->start(time);
      fu50->start(time);
      fl90->start(time);
      fu90->start(time);
    }

    i++;
  }

  wrec->resize(i);
  wmed->resize(i);
  wl50->resize(i);
  wu50->resize(i);
  wl90->resize(i);
  wu90->resize(i);

  frec->resize(i);
  fmed->resize(i);
  fl50->resize(i);
  fu50->resize(i);
  fl90->resize(i);
  fu90->resize(i);

}

void FrequencyCut(wavearray<double>* x, double bF, double eF) {

  // cut frequency range bF,eF
  double F=0.;
  double dF=(x->rate()/(double)x->size())/2.;
  x->FFTW(1);
  for(int j=0;j<x->size()/2;j+=2) {
    F = j*dF;
    if(F<bF || F>eF) {x->data[j]=0;x->data[j+1]=0;}
  }
  x->FFTW(-1);
}

void PlotSpectrogram(wavearray<double>* x, double tstart, double tstop, TString title, TString ofname, double tchirp, double mchirp) {

  if(tstart==0) tstart=x->start();
  if(tstop==0)  tstop=x->stop();

  // START CHIRP FUNCTION
  if(mchirp>0) {
    double Mc = mchirp;     
    double Tc = tchirp;

    const double G  = watconstants::GravitationalConstant();
    const double SM = watconstants::SolarMass();
    const double C  = watconstants::SpeedOfLightInVacuo();
    const double Pi = TMath::Pi();

    double p0 = 256.*Pi/5*pow(G*Mc*SM*Pi/C/C/C, 5./3);
    double p1 = Tc;

    double tmin=Tc-0.9;
    double tmax=Tc;

    // chirp function -> freq = p0 * (time-p1);
    fchirp = new TF1("fchirp", "pow([0]*([1]-x),-3./8.)", tmin, tmax);
    fchirp->SetLineColor(kWhite);
    fchirp->SetLineWidth(1);
    fchirp->SetLineStyle(2);
    fchirp->SetParameter(0, p0);
    fchirp->SetParameter(1, p1);
  } 
  // END CHIRP FUNCTION

  TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",x->start());

  int nfact=4;
  int nfft=nfact*512;
  int noverlap=nfft-10;
  double fparm=nfact*6;
  int ystart = int((tstart-x->start()-1)*x->rate());
  int ystop  = int((tstop-x->start()+1)*x->rate());
  ystart-=nfft;
  ystop+=nfft;
  int ysize=ystop-ystart;
  wavearray<double> y;y.resize(ysize);y.rate(x->rate());y.start(ystart/x->rate());

  // stft use dt=y.rate() to normalize data but whitened data are already normalized by dt 
  // so before stft data must be divided by 1./sqrt(dt)
  for(int i=0;i<(int)y.size();i++) y.data[i]=x->data[i+ystart]/sqrt(y.rate());

  stft = new CWB::STFT(y,nfft,noverlap,"energy","gauss",fparm);

  TCanvas* canvas;
  double xtstart = nfft/x->rate()+ystart/x->rate();
  double xtstop = (ysize-nfft)/x->rate()+ystart/x->rate();

  xtstart+=0.9;xtstop-=0.9;
  stft->Draw(xtstart,xtstop,16,1024,0,0,1);
  if(fchirp!=NULL) fchirp->Draw("SAME");
  stft->GetHistogram()->SetTitle(title);
  stft->GetHistogram()->GetXaxis()->SetTitle(xtitle);
  canvas = stft->GetCanvas();
  canvas->SetLogy(true);
  stft->GetHistogram()->GetXaxis()->SetTitle(xtitle);
  if(ofname!="") stft->Print(ofname);

  y.resize(0);
}

