//
// Display Spetrogram from formula
// Author : Gabriele Vedovato

{

  #define SAMPLERATE 4096

  using namespace CWB;

  wavearray<double> x;
  x.resize(6*SAMPLERATE);
  x.rate(SAMPLERATE);
  x.start(0);
  x=0;
  x[3*SAMPLERATE]=1;

  double pi = TMath::Pi();
  double dt=1./x.rate();
  //for(int i=0;i<x.size();i++) x[i]=sin(2*pi*200*i*dt);
  for(int i=0;i<x.size();i++) x[i]=sin(2*pi*100*(i*dt)*(i*dt));

  int nfact=4;
  int nfft=nfact*512;
  int noverlap=nfft-10;
/*
  int nfact=1;
  int nfft=nfact*128;
  int noverlap=nfft-nfft/4;
*/
  //int noverlap=1;
  double fparm=nfact*6;
  STFT stft(x,nfft,noverlap,"energy","gauss",fparm); 
  //STFT stft(x,nfft,noverlap,"energy","rectangular",fparm); 
  double Tmin=4;
  double Tmax=5;
  double Fmin=800.0;
  double Fmax=1000.0;
  double Zmin=0.0;
  double Zmax=0.0;
  char title[256]="rectangular window";
  stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax,1);
}
