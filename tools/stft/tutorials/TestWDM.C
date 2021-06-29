//
// Display Spetrogram of data from binary file
// Author : Gabriele Vedovato

{

  #define SAMPLERATE (4*4096)

  #define IFILE_NAME "/home/vedovato/Y2/coherent/cwbstft-1.0.0/tst/data/noiseH1_968654512_S6_online_1_848615060_D.dat"
  //#define IFILE_NAME "/home/vedovato/Y2/coherent/cwbstft-1.0.0/tst/data/noiseL1_968654512_S6_online_1_848615060_D.dat"
  //#define IFILE_NAME "/home/vedovato/Y2/coherent/cwbstft-1.0.0/tst/data/noiseV1_968654512_S6_online_1_848615060_D.dat"

  using namespace CWB;

  wavearray<double> x;
  x.ReadBinary(IFILE_NAME);
  x.rate(SAMPLERATE);
  x.start(0);

x.resize(512*SAMPLERATE);
cout << x.size() << endl;
//double df = x.rate()/x.size();
//cout << "df : " << df << endl; 

//  int nfact=4;
  int nfact=1;
//  int nfft=nfact*512;
  int nfft=x.rate()/16;
  int noverlap=nfft/4;
  double fparm=nfact*6;
x=0;
x[x.rate()/2]=1;
  //STFT stft(x,nfft,noverlap,"amplitude","gauss",fparm); 
  STFT stft(x,nfft,noverlap,"energy","gauss",fparm); 
  //STFT stft(x,nfft,"energy","hamming"); 
  //STFT stft(x,nfft,"energy","hann"); 
//  double Tmin=53.5;
//  double Tmax=54.5;
  double Tmin=0;
  double Tmax=1;
  double Fmin=0.0;
  double Fmax=600.0;
  double Zmin=0.0;
  double Zmax=0.0;
  char title[256];
  //sprintf(title,"Spectrogram : %s",IFILE_NAME); 
  //stft.SetTitle(title);
  //TH2D* h2 = stft.GetHistogram();
  //stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax,1);
  stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax);
  //TCanvas* canvas = stft.GetCanvas();
  //canvas->SetLogy(true);

  //stft.Print("H1_BD_amplitude_spectrogram.png");
  //stft.Print("H1_BD_energy_spectrogram.png");
}
