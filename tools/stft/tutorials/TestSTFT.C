//
// Display Spetrogram of data from binary file
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
  cout << x.size() << endl;

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
  sprintf(title,"Spectrogram : %s",IFILE_NAME); 
  stft.SetTitle(title);
  TH2D* h2 = stft.GetHistogram();
  stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax,1);
  TCanvas* canvas = stft.GetCanvas();
  //canvas->SetLogy(true);
  canvas->SetWindowSize(800,600);

  //stft.Print("plots/H1_BD_amplitude_spectrogram.png");
  //stft.Print("plots/H1_BD_energy_spectrogram.png");
}
