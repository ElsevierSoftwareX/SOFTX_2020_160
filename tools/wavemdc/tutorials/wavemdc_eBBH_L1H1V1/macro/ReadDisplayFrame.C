//
// Read data from frame and display content 
// Author : Gabriele Vedovato

#define FRLIST "frames/L1H1V1-wavemdc_eBBH-Log.frl"
#define CHNAME "H1:EBBH-STRAIN_INJ"

// first frame
#define GPS 931158000
#define LENGTH 1000

// first event
//#define GPS 931158090
//#define LENGTH 10

//#define SAVE_PLOTS
//#define DISPLAY_FFT

{
  char gtitle[256];
  TString gfile;

  // read target data
  wavearray<double>  x;
  x.start(GPS); x.stop(GPS+LENGTH);

  CWB::frame frt(FRLIST,CHNAME);
  frt >> x;

  cout << "start " << x.start() << endl;
 
  x.start(0);
  watplot tplot(const_cast<char*>("tplot"),200,20,800,500);
  sprintf(gtitle,"Reconstructed Signal");
  tplot.gtitle(gtitle,"time(sec)","amplitude");
  tplot.goptions("alp", 1, 0., 0.);
  // draw signal
  x >> tplot;

#ifdef SAVE_PLOTS
  // save plot to file
  gfile="hp_time.png";
  tplot >> gfile;
  gfile="hp_time.root";
  tplot >> gfile;
#endif

#ifdef DISPLAY_FFT
  watplot fplot(const_cast<char*>("fplot"),200,20,800,500);
  fplot.gtitle(gtitle,"frequency (Hz)","strain/#sqrt{Hz}");
  fplot.goptions("alp logy", 1, 0., 0., true, 0,0);
  x >> fplot;

#ifdef SAVE_PLOTS
  // save plot to file
  gfile="hp_freq.png";
  fplot >> gfile;
  gfile="hp_freq.root";
  fplot >> gfile;
#endif
#endif
 
}
