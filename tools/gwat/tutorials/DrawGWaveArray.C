// This example shows how to produce time,freq,time-freq plots of
// a signal stored in a wavearray structure

{
  //#define SAVE_PLOT	// uncomment to save plot to file

  #define PLOT_TYPE	GWAT_TIME 	// draw time shifted signal (black)	
  //#define PLOT_TYPE       GWAT_FFT	// draw on the same plot the original signal (red)
  //#define PLOT_TYPE       GWAT_TF		// get pointer to watplot object

  #define SIG_FREQ	10		// signal frequency
  #define TSHIFT	0.05		// time shift (sec)


  gwavearray<double> gw(16384);                 // instantiation with size 16384 samples
  gw.rate(16384);                               // set sample rate = 16384 Hz
  gw.start(0);                                  // set start time = 0

  // fill with a sinusoid @ 100 Hz
  double dt=1/gw.rate();
  for(int i=0;i<gw.size();i++) gw[i]=sin(2*PI*SIG_FREQ*dt*i);

  // draw methods

  //gw.Draw(GWAT_TIME);                           // draw signal in time domain
  //gw.Draw(GWAT_FFT);                            // draw signal in frequency domain
  //gw.Draw(GWAT_TF);                             // draw signal in time-frequency domain


  // apply a time shift of 0.005 sec
  //gw.TimeShift(TSHIFT);
  //gw.Draw(GWAT_TIME,"SAME",kRed);

  // create a new time shifted signal @ SIG_FREQ Hz
  wavearray<double> u;
  u.resize(16384);
  u.rate(16384);
  u.start(0);
  for(int i=0;i<u.size();i++) u[i]=sin(2*PI*SIG_FREQ*(dt*i+TSHIFT));

  // create a new gwavearray using the wavearray u
  gwavearray<double> gu(&u);
  //gwavearray<double> gu = u;

  gu.Draw(PLOT_TYPE);				// draw time shifted signal (black)	
  gu.Draw(&gw,PLOT_TYPE,"SAME",kRed);		// draw on the same plot the original signal (red)

#ifdef SAVE_PLOT
  if(PLOT_TYPE==GWAT_TF) {
    CWB::STFT* stft = gu.GetSTFT();
    stft->Print("gwavearray_plot.png");
  }
#endif

  watplot* plot = gu.GetWATPLOT();		// get pointer to watplot object

  if((PLOT_TYPE==GWAT_TIME) || (PLOT_TYPE==GWAT_FFT)) {

    // set x,y axis titles and plot title
    char gtitle[256];
    sprintf(gtitle,"Original(Red) - Time Shifted(Black)");
    if(PLOT_TYPE==GWAT_TIME) plot->gtitle(gtitle,"time(sec)","amplitude");
    if(PLOT_TYPE==GWAT_FFT)  plot->gtitle(gtitle,"frequency(hz)","amplitude");

    // save plot to file
#ifdef SAVE_PLOT
    TString gfile="gwavearray_plot.png";
    (*plot) >> gfile;

    cout << "created plot file name : " << gfile << endl;
    exit(0);
#endif

  }

  return plot->canvas;  // used by THtml
}
