// This macro use getEBBH & display eBBH hp,hx 

#define DRAW_TIME
//#define DRAW_FFT
//#define SAVE_PLOT

{  
  wavearray<double> hp, hx;
  getEBBH(16.395803, 19.401641, 25.317980, 0.24891893, hp, hx);

  // the distance of source is G*M/c^2  meters 
  double G = watconstants::GravitationalConstant();
  double M = watconstants::SolarMass();
  double c = watconstants::SpeedOfLightInVacuo();
  double pc = watconstants::Parsec();
  double distance_source_Kpc = G*M/(c*c)/pc/1.e3;

  // rescale hp,hx to 10 Kpc
  hp*=distance_source_Kpc/10;
  hx*=distance_source_Kpc/10;

  gwavearray<double> gw(&hp);

#ifdef DRAW_TIME
  gw.Draw();
  gw.Draw(&hx,GWAT_TIME,"SAME",kRed);
#endif

#ifdef DRAW_FFT
  gw.Draw(GWAT_FFT);
  gw.Draw(&hx,GWAT_FFT,"SAME",kRed);
#endif

#ifdef SAVE_PLOT
  watplot* plot = gw.GetWATPLOT();

  char gtitle[256];
  sprintf(gtitle,"eBBH : hp(black) - hx(red)");
#ifdef DRAW_TIME
  plot->gtitle(gtitle,"time(sec)","amplitude");
  TString gfile="eBHH_time_plot.png";
#endif
#ifdef DRAW_FFT
  plot->gtitle(gtitle,"freq(hz)","amplitude");
  TString gfile="eBHH_freq_plot.png";
#endif

  // save plot to file
  (*plot) >> gfile;

  exit(0);
#endif

}
