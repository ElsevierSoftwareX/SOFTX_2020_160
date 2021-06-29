//
// Display Spetrogram of data from binary files
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
  x.start(0);

  int offset=53.2*x.rate();
  for(int i=0;i<2*x.rate();i++) x.data[i]=x.data[i+offset];
  x.resize(2*x.rate());

  int nfft=4*512;
  int noverlap=4*512-10;
  double fparm=24.;
  STFT stft(x,nfft,noverlap,"gauss",fparm); 
  double Tmin=0.0;
  double Tmax=1.0;
  double Fmin=0.0;
  double Fmax=600.0;
  double Zmin=0.0;
  double Zmax=0.0;
  char title[256];
  sprintf(title,"Spectrogram : %s",IFILE_NAME); 
  stft.SetTitle(title);
  stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax,1);
}
