//
// Sim Data & Produce Bicoherence
// Author : Gabriele Vedovato

#define SRATE 1024.
#define BSIZE (1024*10)
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 2
#define BCM_Y_NUM_SLICES 16
#define BCM_X_SLICE_INDEX 0
#define BCM_Y_SLICE_INDEX 0

#define NLOOP 100

{

  int bsize = BSIZE;
  double srate = SRATE;

  CWB::Bicoherence* blc = new CWB::Bicoherence("CH1", "CH2", srate, bsize,
                                               BCM_X_NUM_SLICES, BCM_Y_NUM_SLICES,
                                               BCM_X_SLICE_INDEX, BCM_Y_SLICE_INDEX,
                                               BCM_ORDER);


  double pi = TMath::Pi();
  TRandom normal;

  wavearray<double> x(bsize);
  x.rate(srate);
  wavearray<double> y(bsize/8); 
  y.rate(srate/8);
  for (int cnt=0;cnt<NLOOP;cnt++) {
    //Test with fake injected signal
    for (int j=0;j<x.size();j++) {
      int n=x.size()*cnt+j;
      double s0=sin(150*pi*2*n/x.rate());
      double s1=sin(15*pi*2*n/x.rate());
      double K=10;
      double value = s0+s1+K*s0*s1+normal.Gaus(); 
      x[j] = value;
    }
    for (int j=0;j<y.size();j++) {
      int n=y.size()*cnt+j;
      double s1=sin(15*pi*2*n/y.rate());
      double value = s1+normal.Gaus(); 
      y[j] = value;
    }
    cout << "Loop : " << blc->GetAverages() << endl;
    blc->MakeBicoherence(x,y);
  }
  blc->DrawBicoherence();
//  blc->Reset();
}
