//
// Sim Data & Produce Bicoherence
// Author : Gabriele Vedovato


#define SRATE 1024.
#define BSIZE (1024*10)
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 2
#define BCM_Y_NUM_SLICES 2
#define BCM_X_SLICE_INDEX 0
#define BCM_Y_SLICE_INDEX 0


{

  int bsize = BSIZE;
  double srate = SRATE;

  CWB::Bicoherence* blc = new CWB::Bicoherence("CH1", "CH2", srate, bsize,
                                               BCM_X_NUM_SLICES, BCM_Y_NUM_SLICES,
                                               BCM_X_SLICE_INDEX, BCM_Y_SLICE_INDEX,
                                               BCM_ORDER);


  double pi = TMath::Pi();
  double dt = 1./srate;
  TRandom normal;

  wavearray<double> x(bsize);
  wavearray<double> y; y.resize(0);
  x.rate(srate);
  for (int cnt=0;cnt<100;cnt++) {
    //Test with fake injected signal
    for (int j=0;j<x.size();j++) {
      int n=x.size()*cnt+j;
      double s0=sin(150*pi*2*n*dt);//-cos(10*pi*2*n*dt);
      double s1=sin(15*pi*2*n*dt);//-cos(250*pi*2*n*dt);
      double K=10;
//     double value = 1+s0+s1+K*s0*s1+normal.Gaus(); //GV
      double value = s0+s1+K*s0*s1+normal.Gaus(); //vir
      x[j] = value;
    }
    cout << "Loop : " << blc->GetAverages() << endl;
    blc->MakeBicoherence(x,y);
  }
  blc->DrawBicoherence();
//  blc->Reset();

}
