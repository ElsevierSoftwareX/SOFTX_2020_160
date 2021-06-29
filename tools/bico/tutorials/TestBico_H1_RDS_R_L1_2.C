//
// Read Frame & Produce Bicoherence
// Author : Gabriele Vedovato

#define SRATE 16384.
#define BSIZE (16384*20)
#define BCM_ORDER 2
/*
#define BCM_X_NUM_SLICES 32
#define BCM_Y_NUM_SLICES 32
#define BCM_X_SLICE_INDEX 0
#define BCM_Y_SLICE_INDEX 0
*/
//#define BCM_X_NUM_SLICES 256
#define BCM_X_NUM_SLICES 64
//#define BCM_Y_NUM_SLICES 256
#define BCM_Y_NUM_SLICES 64
#define BCM_X_SLICE_INDEX 1
#define BCM_Y_SLICE_INDEX 0


//#define FRLIST_NAME "input/S6A_R3_H1_LDAS_C02_L2.frames"
//#define CHNAME "H1:LDAS-STRAIN"
//#define CHNAME "L1:LDAS-STRAIN"
//#define CHNAME "V1:h_16384Hz"

#define FRLIST_NAME "input/H1_RDS_R_L1.frl"
//#define CHNAME_1 "H1:LSC-DARM_CTRL"
#define CHNAME_1 "H1:LSC-DARM_ERR"
#define CHNAME_2 "H0:PEM-LVEA2_V1"
//#define CHNAME_2 "H0:PEM-LSC1_MAGX"

#define START 942449664

#define NLOOP 100

{
  // cwb toolbox
  CWB::Toolbox frTB;

  // set & get frame file list
  int nfrFiles=frTB.frl2FrTree(FRLIST_NAME);
  cout << "nfrFiles : " << nfrFiles << endl;

  int bsize = BSIZE;
  double srate = SRATE;

  CWB::Bicoherence* blc = new CWB::Bicoherence(CHNAME_1, CHNAME_2, srate, bsize,
                                               BCM_X_NUM_SLICES, BCM_Y_NUM_SLICES,
                                               BCM_X_SLICE_INDEX, BCM_Y_SLICE_INDEX,
                                               BCM_ORDER);

  blc->SetOutputGraph(0);

  double pi = TMath::Pi();
  double dt = 1./srate;
  TRandom normal;

  double start = START;

  wavearray<double> x(bsize);
  wavearray<double> y; y.resize(0);
  x.rate(srate);
  for (int cnt=0;cnt<NLOOP;cnt++) {
    double stop = start+bsize/srate;
    frfile FRF = frTB.getFrList(start, stop, 0);
    frTB.readFrames(FRF,CHNAME_1,x);  
    frTB.readFrames(FRF,CHNAME_2,y);  
/*   
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
*/
    cout << "Loop : " << blc->GetAverages() << endl;
    blc->MakeBicoherence(x,y);
    start=stop; 
  }
  blc->DrawBicoherence();
//  blc->Reset();

}
