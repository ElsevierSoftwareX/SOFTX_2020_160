// This example shows how to produce a scalogram plot of
// a signal stored in a wavearray structure

{
  //#define SAVE_PLOT			// uncomment to save plot to file

  #define SIG_FREQ	200		// centra signal frequency

  // instantiation
  Meyer<double> S(1024,2);         	// set wavelet for production

  wavearray<double> x(16384);
  x.rate(16384);
  x.start(0);

  // add sin gaussian signal
  double dt=1/x.rate();
  double f=SIG_FREQ;			// frequency (Hz)
  double s=0.01;			// gaussian RMS
  for(int i=0;i<x.size();i++) {
    int j=i-x.size()/2;
    x[i]=exp(-pow(dt*j,2)/2/s/s)*sin(2*PI*f*dt*j);
  }

  WSeries<double> w(x,S);

  // do transformation level=8 and assign to gWSeries class
  gWSeries<double> gw(w);
  gw.Forward(8);
  // or 
//  w.Forward(8);
//  gWSeries<double> gw(&w);

  cout << "level : " << gw.getLevel() << endl;


  // plot scalogram
  gw.DrawSG("FULL");			// full time range
//  gw.DrawSG();			// zoom time range


  watplot* plot = gw.GetWATPLOT();

  // set title
  char gtitle[256];
  sprintf(gtitle,"Scalogram");
  plot->gtitle(gtitle,"time(sec)","amplitude");

  // set the frequency range to be displayed
  plot->hist2D->GetYaxis()->SetRangeUser(SIG_FREQ-100, SIG_FREQ+100);

#ifdef SAVE_PLOT
  // save plot to file
  TString gfile="gwavearray_plot.png";
  plot->canvas->Print(gfile);

  cout << "created plot file name : " << gfile << endl;

  exit(0);
#endif

  return plot->canvas;  // used by THtml 
}
