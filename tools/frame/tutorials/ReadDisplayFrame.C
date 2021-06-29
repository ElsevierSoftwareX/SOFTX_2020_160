//
// Read data from frame and display content 
// Author : Gabriele Vedovato

{
  #define FRLIST "oframes.frl"
  #define CHNAME "H1:FAKE-STRAIN_INJ"
  #define GPS 931072138
  #define LENGTH 64

  //#define SAVE_PLOTS
  #define DISPLAY_FFT

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
  tplot.goptions(const_cast<char*>("alp"), 1, 0., 0.);
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
  fplot.goptions(const_cast<char*>("alp logy"), 1, 0., 0., true, 0,0);
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
