//
// Read Strain from Frame & Produce Bicoherence
// Author : Gabriele Vedovato


#include <vector>

#define SRATE 16384.
#define BSIZE (16384*40)
#define BCM_ORDER 2
/*
#define BCM_X_NUM_SLICES 32
#define BCM_Y_NUM_SLICES 32
#define BCM_X_SLICE_INDEX 0
#define BCM_Y_SLICE_INDEX 0
*/
//#define BCM_X_NUM_SLICES 512
#define BCM_X_NUM_SLICES 128
//#define BCM_Y_NUM_SLICES 512
#define BCM_Y_NUM_SLICES 2048
#define BCM_X_SLICE_INDEX 2
#define BCM_Y_SLICE_INDEX 0

/*
#define FRLIST_NAME_1 "input/S6A_R3_H1_LDAS_C02_L2.frames"
#define CHNAME_1 "H1:LDAS-STRAIN"
#define FRLIST_NAME_2 "input/S6A_R3_H1_LDAS_C02_L2.frames"
#define CHNAME_2 "H1:LDAS-STRAIN"
#define START 931035776
*/

//#define CHNAME "L1:LDAS-STRAIN"
//#define CHNAME "V1:h_16384Hz"

/*
#define FRLIST_NAME_1 "input/VSR2_V1_HrecV3.frames"
#define FRLIST_NAME_2 "input/VSR2_V1_HrecV3.frames"
#define CHNAME_1 "V1:h_16384Hz"
#define CHNAME_2 "V1:h_16384Hz"
#define START 942449664
*/

#define FRLIST_NAME_1 "input/H1_RDS_R_L1.frl"
#define CHNAME_1 "H1:LSC-DARM_ERR"
#define FRLIST_NAME_2 "input/H1_RDS_R_L1.frl"
#define CHNAME_2 "H1:LSC-DARM_ERR"
#define START 942449664

//#define CHNAME_1 "H1:LSC-DARM_CTRL"


#define REBIN 1
//#define REBIN 10
#define NLOOP 100

#define BLC_THRESHOLD 0.5

#define ODIR_NAME "bico_files"

//#define WRITE_BICO_INFOS

void
TestBico_STRAIN(TString chname_2="") {

  if(chname_2.Sizeof()==1) chname_2=CHNAME_1;

  // cwb toolbox
  CWB::Toolbox frTB_1;
  CWB::Toolbox frTB_2;
  // set & get frame file list
  int nfrFiles_1=frTB_1.frl2FrTree(FRLIST_NAME_1);
  cout << "nfrFiles_1 : " << nfrFiles_1 << endl;
  int nfrFiles_2=frTB_2.frl2FrTree(FRLIST_NAME_2);
  cout << "nfrFiles_2 : " << nfrFiles_2 << endl;

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
  wavearray<double> y; y.resize(0);
  x.rate(srate);
  for (int cnt=0;cnt<NLOOP;cnt++) {
    double stop = start+bsize/srate;
    frfile FRF_1 = frTB_1.getFrList(start, stop, 0);
    frTB_1.readFrames(FRF_1,CHNAME_1,x);  
    frfile FRF_2 = frTB_2.getFrList(start, stop, 0);
    frTB_2.readFrames(FRF_2,chname_2.Data(),y);  
    cout << "X Channel rate : " << x.rate() << endl;
    cout << "Y Channel rate : " << y.rate() << endl;
    cout << "Loop : " << blc->GetAverages() << endl;
    blc->MakeBicoherence(x,y);
    start=stop; 
  }

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
