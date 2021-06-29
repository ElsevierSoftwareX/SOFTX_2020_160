//
// Test  Resampling Sata to Power of 2
// Author : Gabriele Vedovato

{
  #define N_IFO	1

  CWB::frame fr;

  TString ifo[N_IFO]={"L1"};
  CWB::mdc MDC(N_IFO,ifo);

//  wavearray<double> w = MDC.GetSGQ(554., 8.9);
//  wavearray<double> w = MDC.GetWNB(3500., 1000., 0.01);
  wavearray<double> w = MDC.GetWNB(500., 4000., 0.01);
  cout << endl;
  cout << "w srate " << w.rate() << endl;
  cout << "w start " << w.start() << endl;
  cout << "w stop  " << w.stop() << endl;
  cout << "w size  " << w.size() << endl;
  cout << endl;

  // resampling
  double srate = 2*4300;
  wavearray<double> ww = w; 

  ww.FFTW(1);
  ww.resize(srate);
  ww.FFTW(-1);
  ww.rate(srate);

  cout << endl;
  cout << "ww srate " << ww.rate() << endl;
  cout << "ww start " << ww.start() << endl;
  cout << "ww stop  " << ww.stop() << endl;
  cout << "ww size  " << ww.size() << endl;
  cout << endl;

  fr.resample(ww,11);
  cout << endl;
  cout << "resampled ww srate " << ww.rate() << endl;
  cout << "resampled ww start " << ww.start() << endl;
  cout << "resampled ww stop  " << ww.stop() << endl;
  cout << "resampled ww size  " << ww.size() << endl;
  cout << endl;

  watplot plot(const_cast<char*>("plot"),200,20,800,500);

  char gtitle[256];  
  sprintf(gtitle,"Original PSD (Black) Reasampled PSD (Red)");
  plot.gtitle(gtitle,"frequency (Hz)","strain/#sqrt{Hz}");
  plot.goptions(const_cast<char*>("alp logy"), 1, 0, 0, true, 0,w.rate()/2., true, 32);

  w >> plot; ww >> plot; 

/*
  char ofName[256];  
  sprintf(ofName,"ResamplePSD.png");
  cout << "write results to " << ofName << endl;
  TString gfile=ofName;
  plot >> gfile;

  exit(0);
*/
}
