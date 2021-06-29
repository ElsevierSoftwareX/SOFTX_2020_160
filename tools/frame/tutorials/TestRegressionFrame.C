#define FRLIST_TARGET  "input/H1_LDAS_C02_L2.frl"
#define FRLIST_WITNESS "input/H1_RDS_R_L1.frl"

#define nWITNESS 9
#define SCRATCH 8.
#define dFREQ 0.5
//#define dFREQ 1.0

//#define LINE_60HZ
//#define LINE_120HZ
#define LINE_180HZ

#define GPS 942449664
#define LENGTH 616

#define nFILTER 2

using namespace CWB;
watplot* plot;

void TestRegressionFrame() {

  TString target = "H1:LDAS-STRAIN";

  TString witness[nWITNESS] = { 
                               "H0:PEM-BSC10_MAGX",
                               "H1:SUS-ETMX_COIL_LL",
                               "H1:SUS-ETMX_COIL_LR",
                               "H1:SUS-ETMX_COIL_UL",
                               "H1:SUS-ETMX_COIL_UR",
                               "H1:SUS-ITMX_COIL_LL",
                               "H1:SUS-ITMX_COIL_LR",
                               "H1:SUS-ITMX_COIL_UL",
                               "H1:SUS-ITMX_COIL_UR"
                              };

#ifdef LINE_60HZ
  double fLow=50; double fHigh=70;
#endif
#ifdef LINE_120HZ
  double fLow=110; double fHigh=130;
#endif
#ifdef LINE_180HZ
  double fLow=170; double fHigh=190;
#endif

  // declare watplot 
  plot = new watplot(const_cast<char*>("plot"),200,20,800,500);
  plot->goptions(const_cast<char*>("alp logy"), 1, 0., 0., true, fLow, fHigh, true, 64);

  // read target data
  wavearray<double>  xx;
  xx.start(GPS); xx.stop(GPS+LENGTH);

  frame frt(FRLIST_TARGET,target);
  frt >> xx;

  // read witness data
  wavearray<double>* yy[nWITNESS];
  frame frw(FRLIST_WITNESS);
  for(int n=0;n<nWITNESS;n++) {
    yy[n] = new wavearray<double>;
    yy[n]->start(GPS); yy[n]->stop(GPS+LENGTH);
    frw.setChName(witness[n]);
    frw >> *yy[n];
    cout << "yy : " << n << " " << yy[n]->size() << " " << yy[n]->rate() << endl;
  }

  // resample target 16384Hz -> 2048Hz
  Biorthogonal<double> Bio(512);
  WSeries<double> wT(Bio);
  wT.Forward(xx,3);
  wT.getLayer(xx,0);

  int level=xx.rate()/(2*dFREQ);
  WDM<double> WD(level, level, 4, 12);

  wT.Forward(xx,WD);
  wT.setlow(fLow);
  wT.sethigh(fHigh);
  regression rr(wT,const_cast<char*>("target"));
  for(int n=0;n<nWITNESS;n++)  rr.add(*yy[n],const_cast<char*>("witness"));
  for(int n=2;n<=nWITNESS;n++) rr.add(1,n,const_cast<char*>("bi-witness"));
  for(int n=2;n<=nWITNESS;n++) rr.mask(n);

  rr.setFilter(nFILTER);
  rr.setMatrix(SCRATCH,1.);
  rr.solve(0.,0,'h');
  //rr.solve(0.,nWITNESS,'h');
  rr.apply();

  wavearray<double> nn = rr.getNoise();
  wavearray<double> cc = rr.getClean();

  xx >> *plot;
  cc >> *plot;

#ifdef LINE_60HZ
  plot->gtitle("H1 regression @ 60Hz");
  *plot >> const_cast<char*>("H1_regression_60Hz.png");
#endif
#ifdef LINE_120HZ
  plot->gtitle("H1 regression @ 120Hz");
  *plot >> const_cast<char*>("H1_regression_1200Hz.png");
#endif
#ifdef LINE_180HZ
  plot->gtitle("H1 regression @ 180Hz");
  *plot >> const_cast<char*>("H1_regression_180Hz.png");
#endif

} 
