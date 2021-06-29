// 
// How to apply regression to subtract power line at 180Hz and related sidebands

//Target channel
#define HCHANNEL "H1:LSC-STRAIN"
#define FRLIST_NAME "S5_H1_RDS_C03_L2.frl"
//Frame list file for auxliary channels
#define FRLIST_NAME2 "H-R-8159.frl"
//Time
#define START 815974196 
#define LENGHT 1200
//Frequency near 180Hz
#define F1 175
#define F2 185

{
  //Auxiliary channels for carrier
  #define NA 1
  TString auch[]={"H0:PEM-COIL_MAGZ","H0:PEM-BSC10_MAGX","H0:PEM-BSC10_MAGY","H0:PEM-BSC10_MAGZ"};

  //Auxiliary channels for side-bands
  #define NB 12
  TString channels[]={"H1:SUS-ITMX_COIL_UL","H1:SUS-ITMX_COIL_LL",
                      "H1:SUS-ETMX_COIL_UR","H1:SUS-ETMY_COIL_UR",
                      "H1:SUS-ETMY_SENSOR_UL","H1:SUS-ETMY_SENSOR_LL",
                      "H1:SUS-ETMY_SENSOR_UR","H1:SUS-ETMY_SENSOR_LR",
                      "H1:SUS-ETMX_SENSOR_UL","H1:SUS-ETMX_SENSOR_LL",
                      "H1:SUS-ETMX_SENSOR_UR","H1:SUS-ETMX_SENSOR_LR"};

  using namespace CWB;
  int totalscratch=32;

  //Fill wavearray h with target channel
  wavearray<double> h;
  h.start(START-totalscratch);
  h.stop(LENGHT+START+totalscratch);
  frame frt(FRLIST_NAME,HCHANNEL);
  frt >> h;

  //Resample to 2048 Hz
  Meyer<double> B(1024);           // set wavelet for resampling
  WSeries<double> ww;
  ww.Forward(h,B,2);
  ww.getLayer(h,0);

  //Make WDM transform, resolution = 1Hz
  int lev=h.rate()/2;
  WDM<double> wdtf(lev, 2*lev, 6, 10);
  WSeries<double> tfmap;
  tfmap.Forward(h, wdtf);

  //Adding target channel on regression object
  regression r;
  r.add(tfmap,const_cast<char*>("hchannel"));
  r.mask(0);
  r.unmask(0,F1,F2);

  //Adding auxiliary channels
  wavearray<double> x;
  x.start(START-totalscratch);
  x.stop(LENGHT+START+totalscratch);
  frame frw(FRLIST_NAME2);
  for (int i=0; i<NA; i++)
  {
   x=0;
   frw.setChName(auch[i]);
   frw >> x;
   r.add(x,const_cast<char*>(auch[i].Data()));
   cout << auch[i].Data() << endl;
  }
  for (int j=0; j<NB; j++)
  {
   x=0;
   frw.setChName(channels[j]);
   frw >> x;
   r.add(x,const_cast<char*>(channels[j].Data()),0,10);
   cout << channels[j].Data() << endl;
  }
  //Multiply carrier channels for side-bands channels
  for (int i=0; i<NA; i++) for (int j=0; j<NB; j++) r.add(i+1,NA+1+j,const_cast<char*>(channels[j].Data()));
  //Mask not useful channels (side-bands) 
  //for (int i=0; i<NA; i++) r.mask(i+1);
  for (int j=0; j<NB; j++) r.mask(NA+1+j);

  //Calculate prediction
  r.setFilter(10);
  r.setMatrix(totalscratch,.95);
  r.solve(0.2, 0, 'h');
  r.apply(0.2);

  //Draw plot of target and target-prediction
  watplot plot;
  plot.goptions(const_cast<char*>("alp logy"), 1, START, START+LENGHT, true, F1, F2,true, 50);
  h >> plot;
  wavearray<double> hh = r.getClean();
  hh >> plot;
  TString ofile = "Regression_H1_bic.png";
  plot >> ofile;

  //Write ranking for each frequency layer
  wavearray<double> freq=r.vfreq;
  wavearray<double> rank=r.getRank(0);
  for (int i=0; i<freq.size(); i++) cout << freq.data[i] << " " << rank.data[i] << endl;
}
