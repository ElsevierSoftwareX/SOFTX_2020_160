//
// Read Frame & Produce Bicoherence
// Author : Gabriele Vedovato


#define SRATE 16384.
#define BSIZE (16384*20)
#define BCM_ORDER 1
/*
#define BCM_X_NUM_SLICES 32
#define BCM_Y_NUM_SLICES 32
#define BCM_X_SLICE_INDEX 0
#define BCM_Y_SLICE_INDEX 0
*/
#define BCM_X_NUM_SLICES 256
#define BCM_Y_NUM_SLICES 256
//#define BCM_Y_NUM_SLICES 1024
#define BCM_X_SLICE_INDEX 1
#define BCM_Y_SLICE_INDEX 1


//#define FRLIST_NAME "input/S6A_R3_H1_LDAS_C02_L2.frames"
//#define CHNAME "H1:LDAS-STRAIN"
//#define CHNAME "L1:LDAS-STRAIN"
//#define CHNAME "V1:h_16384Hz"

#define FRLIST_NAME "input/H1_RDS_R_L1.frl"
//#define CHNAME_1 "H1:LSC-DARM_CTRL"
#define CHNAME_1 "H1:LSC-DARM_ERR"
#define CHNAME_2 "H0:PEM-LVEA_MIC"
//#define CHNAME_2 "H1:OMC-QPD3_SUM_IN1_DAQ"

#define START 942449664

#define NLOOP 100

{
  // cwb toolbox
  CWB::Toolbox frTB;

  CWB::frame frl(FRLIST_NAME,"","README",true);

  frl.setChName(CHNAME_1);
//  frl.setVerbose(true);
  frl.setSRIndex(14);

  int bsize = BSIZE;
  double srate = SRATE;

  CWB::Bicoherence* blc = new CWB::Bicoherence(CHNAME_1, CHNAME_2, srate, bsize,
                                               BCM_X_NUM_SLICES, BCM_Y_NUM_SLICES,
                                               BCM_X_SLICE_INDEX, BCM_Y_SLICE_INDEX,
                                               BCM_ORDER);


  double pi = TMath::Pi();
  double dt = 1./srate;
  TRandom normal;

  double start = START;

  wavearray<double> x(bsize);
  wavearray<double> y(bsize);
  x.rate(srate);
//  for (int cnt=0;cnt<NLOOP;cnt++) {
  int cnt_rejected=0;
  int cnt=0;
  while (cnt<NLOOP) {
    double stop = start+bsize/srate;
    x.start(start); x.stop(stop);
    y.start(start); y.stop(stop);
    if(blc->segListCheck(start,stop)) {
      cout << "READ Data Channel " << CHNAME_1 << endl;
      frl >> x;
      cout << "READ Data Channel " << CHNAME_2 << endl;
      frl >> y;
      cout << "X Channel rate : " << x.rate() << " " << x.rms() << endl;
      cout << "Y Channel rate : " << y.rate() << " " << y.rms() << endl;
/*   
      //Test with fake injected signal
      for (int j=0;j<x.size();j++) {
        int n=x.size()*cnt+j;
        double s0=sin(150*pi*2*n*dt);//-cos(10*pi*2*n*dt);
        double s1=sin(15*pi*2*n*dt);//-cos(250*pi*2*n*dt);
        double K=10;
//       double value = 1+s0+s1+K*s0*s1+normal.Gaus(); //GV
        double value = s0+s1+K*s0*s1+normal.Gaus(); //vir
        x[j] = value;
      }
*/
      cout << "Loop : " << blc->GetAverages() << endl;
      if(blc->MakeBicoherence(x,y)) cnt++; else cout << "Warning - buffer rejected : " << cnt_rejected++ << endl;
    } else cout << "Warning - buffer rejected : " << cnt_rejected++ << endl;
    start=stop;
  }

  cout << "Warning - buffer rejected : " << cnt_rejected << endl;

//  blc->Reset();

}
