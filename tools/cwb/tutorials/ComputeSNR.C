// 
// generate MDC & Read PSD from file & compute SNR
// Author : Gabriele Vedovato

{
  #define ADV_LIGO_PSD    "$HOME_WAT/tools/cwb/plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt"
  #define DELTA_FREQ	0.125	// Hz

  //#define WNB_MDC
  #define RD_MDC
  //#define SGQ_MDC

  #define AMP_FACTOR 1e-21

  #define FLOW	32
  #define FHIGH	2048

  #include <vector>

  TString psdName  = gSystem->ExpandPathName(ADV_LIGO_PSD);

  CWB::Toolbox TB;
  CWB::mdc MDC;

#ifdef WNB_MDC
  wavearray<double> x = MDC.GetWNB(1000., 1000., 0.01);
#endif
#ifdef RD_MDC
  wavearray<double> x = MDC.GetRD(1000., 0.2, 10.);
#endif
#ifdef SGQ_MDC
  wavearray<double> x = MDC.GetSGQ(1053., 9);
#endif

  x*=AMP_FACTOR;			 
  //MDC.Draw(x,MDC_FFT);return;        // Draw Waveform in frequency domain
  //MDC.Draw(x);return;                // Draw Waveform in time domain

  // compute x energy in time domain
  double Et=0;
  double dt=1./x.rate();
  for(int i=0;i<x.size();i++) Et+=x[i]*x[i]*dt;
  cout << "Et = " << Et << endl;

  x.FFT(1);
  double df=(double)x.rate()/(double)x.size();
  x*=1./df;	// double side spectrum

  // compute x energy in frequency domain
  double Ef=0;
  for(int i=0;i<x.size();i++) Ef+=x[i]*x[i]*df;	
  cout << "Ef = " << 2*Ef << endl;

  double fWidth = x.rate()/2;		// bandwidth
  double dFreq  = df;			// frequency resolution
  wavearray<double> psd = TB.GetDetectorPSD(psdName,fWidth,dFreq);	// psd is one side PSD
  //double df=psd.rate()/psd.size();
  cout << "PSD : " << " rate : " << psd.rate() << " size : " << psd.size() << " df " << df << endl;
  psd*=1./sqrt(2); // single -> double side PSD

  // compute SNR
  double fLow=FLOW;
  double fHigh=FHIGH;
  double xsnr=0;
  for(int i=0;i<psd.size();i++) {
    double freq=i*df;
    if(freq<fLow || freq>fHigh) continue;
    xsnr+=(x[2*i]*x[2*i]+x[2*i+1]*x[2*i+1])/pow(psd[i],2)*df;
  }
  xsnr*=2;       // add negative side integration
  xsnr=sqrt(xsnr);
  cout << "SNR " << xsnr << endl;

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
  grPSD.GetHistogram()->GetXaxis()->SetRangeUser(8,8192);
  grPSD.GetHistogram()->GetYaxis()->SetRangeUser(1e-24,1e-21);
  grPSD.GetHistogram()->GetXaxis()->SetTitle("Hz");
  grPSD.GetHistogram()->GetYaxis()->SetTitle("PSD");
  grPSD.GetHistogram()->GetYaxis()->SetTitleOffset(1.2);

  grPSD.Draw("ALP");
  grMDC.Draw("SAME");
}

