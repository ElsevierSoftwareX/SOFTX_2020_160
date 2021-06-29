//
// Read Frame & Produce Bicoherence
// Author : Gabriele Vedovato

#include <vector>

#define SRATE 16384.
#define BSIZE (16384*100)
#define BCM_ORDER 2
/*
#define BCM_X_NUM_SLICES 32
#define BCM_Y_NUM_SLICES 32
#define BCM_X_SLICE_INDEX 0
#define BCM_Y_SLICE_INDEX 0
*/
#define BCM_X_NUM_SLICES 128
//#define BCM_X_NUM_SLICES 64
#define BCM_Y_NUM_SLICES 512
//#define BCM_Y_NUM_SLICES 256
#define BCM_X_SLICE_INDEX 2
#define BCM_Y_SLICE_INDEX 0


//#define FRLIST_NAME "input/S6A_R3_H1_LDAS_C02_L2.frames"
//#define CHNAME "H1:LDAS-STRAIN"
//#define CHNAME "L1:LDAS-STRAIN"
//#define CHNAME "V1:h_16384Hz"

/*
#define FRLIST_NAME_1 "input/H1_LDAS_C02_L2.frl"
#define FRLIST_NAME_2 "input/H1_LDAS_C02_L2.frl"
#define CHNAME_1 "H1:LDAS-STRAIN"
#define CHNAME_2 "H1:LDAS-STRAIN"
*/
#define FRLIST_NAME_1 "input/H1_RDS_R_L1.frl"
#define FRLIST_NAME_2 "input/H1_RDS_R_L1.frl"
#define CHNAME_1 "H1:LSC-DARM_ERR"
#define CHNAME_2 "H0:PEM-LSC1_MAGY"
//#define CHNAME_1 "H1:LSC-DARM_CTRL"
//#define CHNAME_2 "H0:PEM-LVEA2_V1"
//#define CHNAME_2 "H0:PEM-LSC1_MAGZ"
//#define CHNAME_2 "H0:PEM-HAM6_ACCX"
//#define CHNAME_2 "H0:PEM-COIL_MAGX"
//#define CHNAME_2 "H0:PEM-ISCT10_ACCX"
//#define CHNAME_2 "H0:PEM-LVEA_SEISX"
//#define CHNAME_2 "H0:PEM-BSC1_MAG1X"
//#define CHNAME_2 "H0:PEM-BSC1_MAG1Y"
//#define CHNAME_2 "H0:PEM-BSC1_MAG1Z"

#define START 942449664

#define REBIN 80
//#define REBIN 10
#define NLOOP 100

#define BLC_THRESHOLD 0.5

#define ODIR_NAME "bico_files"

//#define WRITE_BICO_INFOS

void
TestBico_H1_RDS_R_L1_3(TString chname_2="") {

  if(chname_2.Sizeof()==1) chname_2=CHNAME_2;

  CWB::frame frl1(FRLIST_NAME_1,"","README",true);
  CWB::frame frl2(FRLIST_NAME_2,"","README",true);

  frl1.setChName(CHNAME_1);
  frl2.setChName(CHNAME_2);
//  frl1.setVerbose(true);
//  frl2.setVerbose(true);
  frl1.setSRIndex(14);
  frl2.setSRIndex(14);

  int bsize = BSIZE;
  double srate = SRATE;

  CWB::Bicoherence* blc = new CWB::Bicoherence(CHNAME_1, chname_2.Data(), srate, bsize,
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
      frl1 >> x;
      cout << "READ Data Channel " << CHNAME_2 << endl;
      frl2 >> y;
      cout << "X Channel rate : " << x.rate() << " " << x.rms() << endl;
      cout << "Y Channel rate : " << y.rate() << " " << y.rms() << endl;
      cout << "Loop : " << blc->GetAverages() << endl;
      if(blc->MakeBicoherence(x,y)) cnt++; else cout << "Warning - buffer rejected : " << cnt_rejected++ << endl;
    } else cout << "Warning - buffer rejected : " << cnt_rejected++ << endl;
    start=stop;
  }

  cout << "Warning - buffer rejected : " << cnt_rejected << endl;


#ifdef WRITE_BICO_INFOS
  vector<bico> listBico = blc->GetBicoherence(BLC_THRESHOLD,REBIN);
  cout << "list size : " << listBico.size() << endl;
  char ofile[256];
  sprintf(ofile,"%s/bico_%s_%s.txt",ODIR_NAME,CHNAME_1,chname_2.Data());
  ofstream out;
  out.open(ofile,ios::out);

  for(int i=0;i<listBico.size();i++) {
    printf("%d -> %3.2f - %3.2f : %1.2f\n",i,listBico[i].x,listBico[i].y,listBico[i].c);
    char ostr[256];
    sprintf(ostr,"%d -> %3.2f - %3.2f : %1.2f\n",i,listBico[i].x,listBico[i].y,listBico[i].c);
    out << ostr <<endl;
  }
  out.close();
#else
  blc->DrawBicoherence(REBIN);
#endif
//  blc->Reset();

}
