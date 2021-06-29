//
// Display Spetrogram of white noise
// Author : Gabriele Vedovato

{

  #define SAMPLERATE 4096

  using namespace CWB;

  wavearray<double> x(70*SAMPLERATE);
  x.rate(SAMPLERATE);
  x.start(0);

  double sigma=sqrt(SAMPLERATE);
  for(int i=0;i<x.size();i++) x.data[i]=gRandom->Gaus(0,sigma);

  double dt=1./x.rate();
  double rms=0;
  for(int i=0;i<x.size();i++) rms+=dt*pow(x.data[i],2);
  rms/=x.size();
  rms=sqrt(rms);
  cout << "rms : " << rms << endl;

  int nfft=4*512;
  int noverlap=4*512-10;
  double fparm=24.;
  STFT stft(x,nfft,noverlap,"amplitude","gauss",fparm); 
  //STFT stft(x,nfft,noverlap,"energy","gauss",fparm); 
  //STFT stft(x,nfft,"energy","hamming"); 
  //STFT stft(x,nfft,"energy","hann"); 
  double Tmin=53.2;
  double Tmax=54.1;
  double Fmin=0.0;
  double Fmax=600.0;
  //double Fmin=1.0;
  //double Fmax=1024.0+512.0;
  double Zmin=0.0;
  double Zmax=0.0;
  char title[256];
  sprintf(title,"Spectrogram"); 
  stft.SetTitle(title);
  TH2D* h2 = stft.GetHistogram();
  stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax,1);
  TCanvas* canvas = stft.GetCanvas();
  //canvas->SetLogy(true);
  //canvas->SetWindowSize(800,600);

  //stft.Print("plots/H1_BD_amplitude_spectrogram.png");
  //stft.Print("plots/H1_BD_energy_spectrogram.png");
}
