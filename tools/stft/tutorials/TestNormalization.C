//
// Display Spetrogram (Test)
// Author : Gabriele Vedovato

{

  #define SAMPLERATE 4096

  #define IFILE_NAME "data/noiseH1_968654512_S6_online_1_848615060_D.dat"
  //#define IFILE_NAME "data/noiseL1_968654512_S6_online_1_848615060_D.dat"
  //#define IFILE_NAME "data/noiseV1_968654512_S6_online_1_848615060_D.dat"

  using namespace CWB;

  wavearray<double> x;
  x.ReadBinary(IFILE_NAME);
  x.rate(SAMPLERATE);
  x.start(100);

  int nscratch=8*x.rate();
  double dt=1./x.rate();

/*
  wavearray<double> t;
  t.resize(x.size());
  for(int i=0;i<t.size();i++) t.data[i]=i*dt;

  TGraph gr(t.size()-2*nscratch,&t.data[nscratch],&x.data[nscratch]);
  gr.Draw("ALP");
  return;
*/
/*
  //for(int i=0;i<x.size();i++) x.data[i]+=100;
  x.FFT(1);
  double df=x.rate()/x.size();
  for(int i=0;i<x.size();i++) {
    double f=i*df/2;
    if(f<1) x.data[i]=0;
  }
  x.FFT(-1);
*/

  //for(int i=0;i<x.size();i++) x.data[i]*=sqrt(x.rate());

  int N=0; 
  double avr=0;
  for(int i=nscratch;i<x.size()-nscratch;i++) {avr+=x.data[i];N++;}
  avr/=N;
  cout << "avr : " << avr << endl;
  N=0;
  double rms=0;
  for(int i=nscratch;i<x.size()-nscratch;i++) {rms+=dt*pow(x.data[i]-avr,2);N++;}
  rms/=N;
  rms=sqrt(rms);
  cout << "rms : " << rms << endl;
  //exit(0);

  int nfact=4;
  int nfft=nfact*512;
  int noverlap=nfft-10;
  double fparm=nfact*6;
  //STFT stft(x,nfft,noverlap,"amplitude","gauss",fparm); 
  STFT stft(x,nfft,noverlap,"energy","gauss",fparm); 
  //STFT stft(x,nfft,"energy","hamming"); 
  //STFT stft(x,nfft,"energy","hann"); 
  double Tmin=53.5+100;
  double Tmax=54.5+100;
  double Fmin=64.0;
  double Fmax=2048.0;
  //double Fmin=1.0;
  //double Fmax=1024.0+512.0;
  double Zmin=0.0;
  double Zmax=0.0;
  char title[256];
  //sprintf(title,"Spectrogram : %s",IFILE_NAME); 
  //stft.SetTitle(title);
  TH2D* h2 = stft.GetHistogram();
  stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax,1);
  TCanvas* canvas = stft.GetCanvas();
  //canvas->SetLogy(true);
  //canvas->SetWindowSize(800,600);

  //stft.Print("plots/H1_BD_amplitude_spectrogram.png");
  //stft.Print("plots/H1_BD_energy_spectrogram.png");
}
