//
// Read data from frame and display content 
// Author : Gabriele Vedovato

//#define FRLIST   "input/L1.lst"
//#define CHNAME   "L1:FAKE-GLITCHES"
//#define ADV_PSD  "$HOME_WAT/tools/cwb/plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt"

//#define FRLIST   "input/H1.lst"
//#define CHNAME   "H1:FAKE-GLITCHES"
//#define ADV_PSD  "$HOME_WAT/tools/cwb/plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt"

#define FRLIST  "input/V1.lst"
#define CHNAME  "V1:FAKE-GLITCHES"
#define ADV_PSD "$HOME_WAT/tools/cwb/plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt"

#define GPS 1040739010
#define LENGTH 10

#define RESAMPLE

{
  char gtitle[256];
  TString gfile;

  TString psdName  = gSystem->ExpandPathName(ADV_PSD);

  CWB::Toolbox TB;
  CWB::mdc MDC;   

  // read target data
  wavearray<double>  x;
  x.start(GPS-LENGTH); x.stop(GPS+LENGTH);

  CWB::frame frt(FRLIST,CHNAME);
  frt >> x;

  cout << "start " << x.start() << endl;
  x.start(0);

#ifdef RESAMPLE
  // resample glitch channel
  int levelR = 2;
  WSeries<double> wM;              // WSeries
  Meyer<double> B(1024);           // set wavelet for resampling
  wM.Forward(x,B,levelR);
  wM.getLayer(x,0);
#endif

//  for(int i=x.size()/2;i<x.size();i++) x[i]=0;
  //MDC.Draw(x,MDC_FFT);return;        // Draw Waveform in frequency domain
  MDC.Draw(x);return;                // Draw Waveform in time domain     

  // compute x energy in time domain
  double Et=0;                      
  double dt=1./x.rate();            
  for(int i=0;i<x.size();i++) Et+=x[i]*x[i]*dt;
  cout << "Et = " << Et << endl;               

  x.FFTW(1);
  double df=(double)x.rate()/(double)x.size();
  x*=1./df;     // double side spectrum       

  // compute x energy in frequency domain
  double Ef=0;                           
  for(int i=0;i<x.size();i++) Ef+=x[i]*x[i]*df;
  cout << "Ef = " << 2*Ef << endl;               

  double fWidth = x.rate()/2;           // bandwidth
  double dFreq  = df;                   // frequency resolution
  wavearray<double> psd = TB.GetDetectorPSD(psdName,fWidth,dFreq);      // psd is one side PSD
  //double df=psd.rate()/psd.size();                                                          
  cout << "PSD : " << " rate : " << psd.rate() << " size : " << psd.size() << " df " << df << endl;
  psd*=1./sqrt(2); // single -> double side PSD

  // compute SNR
  double fLow=64;
  double fHigh=2048;
  double snr=0; 
  for(int i=0;i<psd.size();i++) {
    double freq=i*df;
    if(freq<fLow || freq>fHigh) continue;
    snr+=(x[2*i]*x[2*i]+x[2*i+1]*x[2*i+1])/pow(psd[i],2)*df;
  }
  snr*=2;       // add negative side integration
  snr=sqrt(snr);                                                                        
  cout << "SNR " << snr << endl;                                                        

  // define frequency array
  wavearray<double> f=psd; 
  for(int i=0;i<f.size();i++) f[i]=i*df;

  // define mdc amplitude array
  wavearray<double> a=psd;
  for(int i=0;i<a.size();i++) a[i]=sqrt(x[2*i]*x[2*i]+x[2*i+1]*x[2*i+1]);

  TObjArray* token = TString(psdName).Tokenize(TString('/'));
  TObjString* sname = (TObjString*)token->At(token->GetEntries()-1);
  TString Title = sname->GetString();                               

  // display PSD
  TCanvas canvas;
  canvas.SetLogx();
  canvas.SetLogy();
  canvas.SetGridx();
  canvas.SetGridy();
  TGraph grMDC(psd.size(),f.data,a.data);
  grMDC.SetMarkerColor(kBlue);
  grMDC.SetLineColor(kBlue);
  TGraph grPSD(psd.size(),f.data,psd.data);
  grPSD.SetTitle(Title);
  grPSD.SetMarkerColor(kRed);
  grPSD.SetLineColor(kRed);


  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle("Advanced Detectors : Design Sensitivity Curves");
  //mg->SetTitle("");

  mg->Add(&grPSD);
  mg->Add(&grMDC);

  mg->Paint("APL");

  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);

  mg->GetHistogram()->GetXaxis()->SetRangeUser(8,8192);
  mg->GetHistogram()->GetYaxis()->SetRangeUser(1e-25,1e-21);

  //mg->GetXaxis()->SetTitle(gr[0]->GetXaxis()->GetTitle());
  mg->GetXaxis()->SetLabelFont(42);
  mg->GetYaxis()->SetLabelFont(42);
  mg->GetXaxis()->SetTitleFont(42);
  mg->GetYaxis()->SetTitleFont(42);
  mg->GetXaxis()->SetTitleOffset(1.20);
  mg->GetYaxis()->SetTitleOffset(1.05);
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitle("Frequency (Hz)      ");
  mg->GetYaxis()->SetTitle("#frac{1}{#sqrt{Hz}}          ");

  mg->Draw("APL");

/*
  grPSD.GetHistogram()->GetXaxis()->SetRangeUser(8,8192);
  grPSD.GetHistogram()->GetYaxis()->SetRangeUser(1e-24,1e-21);
  grPSD.GetHistogram()->GetXaxis()->SetTitle("Hz");
  grPSD.GetHistogram()->GetYaxis()->SetTitle("PSD");
  grPSD.GetHistogram()->GetYaxis()->SetTitleOffset(1.2);

  grPSD.Draw("ALP");
  grMDC.Draw("SAME");
*/
}

