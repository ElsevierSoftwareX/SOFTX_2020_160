
#define FRLIST_1 "frames/L1H1V1-wavemdc_EOBNRv2HMpseudoFourPN-Log.frl"
#define CHNAME_1 "H1:FAKE-STRAIN_INJ"
#define GPS_1 931072130

#define FRLIST_2 "input/chris.frl"
#define CHNAME_2 "H1:FAKE-STRAIN_INJ"
#define GPS_2 931072130

#define LENGTH 64

//#define SAVE_PLOTS
//#define DISPLAY_FFT

{
  char gtitle[256];
  TString gfile;

  // read frame_1 
  wavearray<double>  x1;
  x1.start(GPS_1); x1.stop(GPS_1+LENGTH);
  CWB::frame frt_1(FRLIST_1,CHNAME_1);
  frt_1 >> x1;

  // read frame_2 
  wavearray<double>  x2;
  x2.start(GPS_2); x2.stop(GPS_2+LENGTH);
  CWB::frame frt_2(FRLIST_2,CHNAME_2);
  frt_2 >> x2;

  //x2-=x1;   // diff

  x1.start(0);
  x2.start(0);
  watplot tplot(const_cast<char*>("tplot"),200,20,800,500);
  sprintf(gtitle,"Reconstructed Signal");
  tplot.gtitle(gtitle,"time(sec)","amplitude");
  tplot.goptions("alp", 1, 0., 0.);
  // draw signal
  x1 >> tplot;
  x2 >> tplot;

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
