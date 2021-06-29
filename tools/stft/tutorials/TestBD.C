//
// Display Big Dog Spetrogram from binary file
// Author : Gabriele Vedovato

{

  #define SAMPLERATE 4096

  #define IFILE_NAME "/home/vedovato/Y2/coherent/cwbstft-1.0.0/tst/data/noiseH1_968654512_S6_online_1_848615060_D.dat"
  //#define IFILE_NAME "/home/vedovato/Y2/coherent/cwbstft-1.0.0/tst/data/noiseL1_968654512_S6_online_1_848615060_D.dat"
  //#define IFILE_NAME "/home/vedovato/Y2/coherent/cwbstft-1.0.0/tst/data/noiseV1_968654512_S6_online_1_848615060_D.dat"

  using namespace CWB;

  wavearray<double> x;
  x.ReadBinary(IFILE_NAME);
  x.rate(SAMPLERATE);
  x.start(0);
  cout << x.size() << " " << x.size()/x.rate() << endl;
//exit(0);
/*
  int M=512;
  int N=16;
  int P=1;

  M/=P; N/=P;
  int nfact=1;
  int nfft=2*M;
  int noverlap=nfft-nfft/N;
*/

  int nfact=4;
  int nfft=nfact*512;
  int noverlap=nfft-10;

  //int noverlap=nfft/sqrt(2);
  double fparm=nfact*6;
  //STFT stft(x,nfft,noverlap,"amplitude","gauss",fparm); 
  STFT stft(x,nfft,noverlap,"energy","gauss",fparm); 
  //STFT stft(x,nfft,noverlap,"energy","hamming",); 
  //STFT stft(x,nfft,noverlap,"energy","hann"); 
  //STFT stft(x,nfft,noverlap,"energy","rectangular"); 
  double Tmin=53.5;
  double Tmax=54.5;
  double Fmin=0.0;
  double Fmax=350.0;
  double Zmin=0.0;
  double Zmax=0.0;
  char title[256];
  //TH2D* h2 = stft.GetHistogram();
  //sprintf(title,"BigDog STFTScan : M=%d - N=%d (Hann)",M,N);
//  sprintf(title,"BigDog STFTScan : M=%d - N=%d (Gauss)",M,N);
  stft.SetTitle("BigDog STFTScan : M=1024 - N=102 (Gauss)");
  stft.SetTitle(title);
  stft.Draw(Tmin,Tmax,Fmin,Fmax,Zmin,Zmax,1);
  //TCanvas* canvas = stft.GetCanvas();
  //canvas->SetLogy(true);

  //stft.Print("H1_BD_amplitude_spectrogram.png");
  //stft.Print("H1_BD_energy_spectrogram.png");
  //char ofile[256];sprintf(ofile,"BD_STFT_HANN_M%d_N%d.png",M,N);

//  char ofile[256];sprintf(ofile,"BD_STFT_GAUSS_M%d_N%d.png",M,N);
  char ofile[256];sprintf(ofile,"plots/BD_STFT_GAUSS_M%d_N%d.png",1024,102);
  stft.Print(ofile);
}
