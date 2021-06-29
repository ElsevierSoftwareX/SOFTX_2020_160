//
//How to apply regression to subtract sinusoidal line at 100Hz
//Generate time series w(t) and x(t) each containing Gaussian white noise (unit variance)
//Construct the witness channel: A(t) = 0.1*w(t) + 0.25*sin(2pi*100*t)
//Construct the target channel:  H(t) = x(t) + sin(2pi*100*t)
//Compare different regression parameters
//Plot Fourier Transform of target and cleaned
//


{
  //Time
  #define LENGHT 1200
  #define RATE 2048
  //Frequency near 180Hz
  #define F1 95
  #define F2 105

  //Auxiliary channels

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
  double omega = TMath::TwoPi()*100;
  for(int i=0; i<N; i++) h[i] += sin(omega*i*T);

  //Make WDM transform, resolution = 1Hz
  int lev=h.rate()/2;
  WDM<double> wdtf(lev, 2*lev, 6, 10);
  WSeries<double> tfmap;
  tfmap.Forward(h, wdtf);

  //Adding auxiliary channels
  wavearray<double> x;
  x.rate(RATE);
  x.resize(N);
  x.start(0);
  x.stop(LENGHT+2*totalscratch);

  for(int i=0; i<N; i++) x[i] = 0.1*rnd.Gaus();
  // add a sine wave @ 100Hz:
  T = 1./RATE;
  omega = TMath::TwoPi()*100;
  for(int i=0; i<N; i++) x[i] += 0.25*sin(omega*i*T);

  //Adding channels
  regression r;
  r.add(tfmap,const_cast<char*>("hchannel"));
  r.mask(0);
  r.unmask(0,F1,F2);
  r.add(x,const_cast<char*>("witness"));

  regression r2=r;
  regression r3=r;
  regression r4=r;

  //Calculate prediction
  //r.setFilter(10);
  //r.setMatrix(totalscratch,.95);
  //r.solve(0.2, 0, 'h');
  //r.apply(0.2);
  r.setFilter(5);
  r.setMatrix(totalscratch);
  r.solve(0., 4, 'h');
  r.apply(0.2);

  r2.setFilter(5);
  r2.setMatrix(totalscratch,.95);
  r2.solve(0, 4, 'h');
  r2.apply(0.8);

  r3.setFilter(5);
  r3.setMatrix(totalscratch,.95);
  r3.solve(0, 2, 'h');
  r3.apply(0.2);

  r4.setFilter(2);
  r4.setMatrix(totalscratch,.95);
  r4.solve(0., 4, 'h');
  r4.apply(0.2);

  //Draw plot of target and target-prediction
  TCanvas *c1 = new TCanvas("comp","comp",0,0,800,600);
  c1->Divide(2,2);

  watplot plot(const_cast<char*>("plot"));

  c1->cd(1);
  gPad->SetLogy();
  plot.plot(h,const_cast<char*>("alp"), 1, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  wavearray<double> clean = r.getClean();
  plot.plot(clean,const_cast<char*>("same"), 2, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.graph[0]->SetTitle("Black : noisy data | Red : cleaned data");

  c1->cd(2);
  gPad->SetLogy();
  plot.plot(h,const_cast<char*>("alp"), 1, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  clean = r2.getClean();
  plot.plot(clean,const_cast<char*>("same"), 2, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.graph[2]->SetTitle("Black : noisy data | Red : cleaned data");

  c1->cd(3);
  gPad->SetLogy();
  plot.plot(h,const_cast<char*>("alp"), 1, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  clean = r3.getClean();
  plot.plot(clean,const_cast<char*>("same"), 2, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.graph[4]->SetTitle("Black : noisy data | Red : cleaned data");

  c1->cd(4);
  gPad->SetLogy();
  plot.plot(h,const_cast<char*>("alp"), 1, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  clean = r4.getClean();
  plot.plot(clean,const_cast<char*>("same"), 2, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.graph[6]->SetTitle("Black : noisy data | Red : cleaned data");

  c1->SaveAs("Regression_Sine_parameters.png");

}
