//
// This example show how to apply the band filters
// Author : Gabriele Vedovato

{
  //#define WAVELET		// comment to select WDM
  #define FLOW 	90.	// low frequency cut
  #define FHIGH	150.    // high frequency cut

  TStopwatch watchJob;  // job benchmark

  double flow = FLOW;
  double fhigh = FHIGH;	

  detector D;

  double rate = 4096.;	// sample rate
  double length = 600; 	// sec

  wavearray<double> x(length*rate); x.rate(rate);
  for(int i=0;i<x.size();i++) x[i]=gRandom->Gaus(0,1);	// fill with gaussian noise
  gwavearray<double> gx(&x);

  watchJob.Start();

  // decomposition	
#ifdef WAVELET
  int level = 11;			// df = 1Hz
  Meyer<double> S(1024,2);         	// set wavelet for production
  D.getTFmap().Forward(x,S,level);
#else
  int layers = x.rate()/2;
  WDM<double> wdm(layers,layers,6,10);  // df = (x.rate()/2)/layers = 1Hz
  D.getTFmap()->Forward(x,wdm); 	
#endif

  if(D.getTFmap()->pWavelet->m_WaveType==WDMT)  cout << "WDM" << endl;
  else						cout << "WAVELET" << endl;
 
  D.getTFmap()->setlow(flow);
  D.getTFmap()->sethigh(fhigh);

  D.bandPass(flow,fhigh);		// band pass filter
//  D.bandCut(flow,fhigh);		// band cut  filter
//  D.lowPass(flow);			// low  pass filter
//  D.highPass(fhigh);			// high pass filter

  D.getTFmap()->Inverse();		// return to time domain

  watchJob.Stop();
  cout << "Job   Elapsed Time - " << watchJob.RealTime() << " sec" << endl;

  wavearray<double> y = *D.getTFmap();  // get band filtered data

  gx.Draw(GWAT_FFT,"ALP NOLOGX NOLOGY");		// original data
  gx.Draw(&y,GWAT_FFT,"SAME NOLOGX NOLOGY",kRed);	// band filtered data
}
