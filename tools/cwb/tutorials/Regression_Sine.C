//
//How to apply regression to subtract sinusoidal line at 100Hz
//Generate time series w(t) and x(t) each containing Gaussian white noise (unit variance)
//Construct the witness channel: A(t) = 0.1*w(t) + 0.25*sin(2pi*100*t)
//Construct the target channel:  H(t) = x(t) + sin(2pi*100*t)
//Plot Fourier Transform of target and cleaned
//

{
  //Time
  #define LENGHT 1200
  #define RATE 2048
  //Frequency near 180Hz
  #define F1 95
  #define F2 105

  double omega = TMath::TwoPi()*100;

  using namespace CWB;
  int totalscratch=32;

  //Fill wavearray h with target channel
  int N = RATE*(LENGHT+2*totalscratch);
  wavearray<double> h;
  h.rate(RATE);
  h.resize(N);
  h.start(0);
  h.stop(LENGHT+2*totalscratch);

  // time series is filled with white noise data: 
  TRandom3 rnd(0);   
  for(int i=0; i<N; i++) h[i] = rnd.Gaus();

  // add a sine wave @ 100Hz:
  double T = 1./RATE;
  for(int i=0; i<N; i++) h[i] += sin(omega*i*T);

  //Make WDM transform, resolution = 1Hz
  int lev=h.rate()/2;
  WDM<double> wdtf(lev, 2*lev, 6, 10);
  WSeries<double> tfmap;
  tfmap.Forward(h, wdtf);

  //Defining auxiliary channels
  wavearray<double> x;
  x.rate(RATE);
  x.resize(N);
  x.start(0);
  x.stop(LENGHT+2*totalscratch);

  // time series is filled with white noise data: 
  for(int i=0; i<N; i++) x[i] = 0.1*rnd.Gaus();

  // add a sine wave @ 100Hz:
  T = 1./RATE;
  for(int i=0; i<N; i++) x[i] += 0.25*sin(omega*i*T);

  //Adding target channel
  regression r;
  r.add(tfmap,const_cast<char*>("hchannel"));

  //Considering only the interested frequency band
  r.mask(0);
  r.unmask(0,F1,F2);

  //Adding witness channel
  r.add(x,const_cast<char*>("witness"));

  //Calculate prediction
  r.setFilter(5);
  r.setMatrix(totalscratch);
  r.solve(0., 4, 'h');
  r.apply(0.2);

  //Draw plot of target and target-prediction
  watplot plot;
  gPad->SetLogy();
  plot.plot(h,const_cast<char*>("alp"), 1, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  wavearray<double> clean = r.getClean();
  plot.plot(clean,const_cast<char*>("same"), 2, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.graph[0]->SetTitle("Black : noisy data | Red : cleaned data");
  plot >> const_cast<char*>("Regression_Sine.png");

  //Write ranking for each frequency layer
  wavearray<double> freq=r.vfreq;
  wavearray<double> rank=r.getRank(0);
  cout << "Ranking:" << endl;
  for (int i=0; i<freq.size(); i++) cout << freq.data[i] << " " << rank.data[i] << endl;
  cout << endl;

  //Write eigenvalues for each frequency layer
  cout << "Eigen-values" << endl;
  wavearray<double>* eigen = new wavearray<double>[freq.size()];
  for (int i=0; i<freq.size(); i++) eigen[i]=r.getVEIGEN(i);
  for (int j=0; j<eigen[0].size(); j++)
   {
    for (int i=0; i<freq.size(); i++)
     printf("%.3f\t",eigen[i].data[j]);
    cout << endl;
   }
}
