//
// Display Spetrogram of a delta signal
// Author : Gabriele Vedovato

{

  #define SAMPLERATE 4096

  using namespace CWB;

  wavearray<double> x;
  x.resize(76*SAMPLERATE);
  x.rate(SAMPLERATE);
  x.start(0);
  x=0;
  x[54*SAMPLERATE]=1;

  int nfact=1;
  int nfft=nfact*128;
  int noverlap=nfft-nfft/4;
  //int noverlap=1;
  double fparm=nfact*6;
  STFT stft(x,nfft,noverlap,"energy","gauss",fparm); 
  //STFT stft(x,nfft,noverlap,"energy","rectangular",fparm); 
  double Tmin=53.95;
  double Tmax=54.05;
  double Fmin=0.0;
  double Fmax=350.0;
  double Zmin=0.0;
  double Zmax=0.0;
  char title[256]="rectangular window";
  stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax);
}
