//
// How to apply regression to subtract upconverted 2Hz noise at 100 +/- 2 Hz
// before the upconversion the 2Hz noise is modify by a transfer function
// the transfer function take into account the coupling between witness and target channels
// Note: in this example the data are cleaned without removing the lines at 100 Hz and 102Hz
// 

{
  //Time
  #define LENGHT 1200
  #define RATE 2048
  //Frequency near 180Hz
  #define F1 95
  #define F2 105

  //Draw plot of target and target-prediction
  watplot plot(const_cast<char*>("plot"));
  TCanvas *c1 = plot.canvas;
  c1->Divide(2,2);

  double omega = TMath::TwoPi()*100;
  double omega_lf = TMath::TwoPi()*2;

  int totalscratch=32;
  int N = RATE*(LENGHT+2*totalscratch);
  double dt = 1./double(RATE);
  double df = double(RATE)/N/2.;

  TRandom3 rnd(1);   

  using namespace CWB;

  //Fill wavearray h with target channel
  wavearray<double> target;
  target.rate(RATE);
  target.resize(N);
  target.start(0);
  target.stop(LENGHT+2*totalscratch);
  for(int i=0; i<N; i++) target[i] = rnd.Gaus();	// white gaussian noise

  // generate witness channel : gauss noise @ 2Hz with gauss shape 
  wavearray<double> witness = target;
  for(int i=0; i<N; i++) witness[i] = rnd.Gaus();	// white gaussian noise
  witness.FFT(1);
  for(int i=0;i<N;i++) {				// select freq range [1,3] with gauss shape
    double f = double(i)*df;
    witness[i]*=1000.*exp(-pow(2*(f-2.)/1.,2));
    witness[i]/=RATE;
  }
  witness.FFT(-1);
  for(int i=0; i<N; i++) witness[i] += 0.01*rnd.Gaus();// add white noise

  // plot witness channel
  c1->cd(1);
  gPad->SetLogy();
  plot.plot(witness, const_cast<char*>("alp logy"), 1, totalscratch, totalscratch+LENGHT, true, 0, 5, true, 50);

  // define 100Hz line 
  wavearray<double> line_100hz = target;
  for(int i=0; i<N; i++) line_100hz[i] = sin((omega)*i*dt);

  // define 102Hz line 
  wavearray<double> line_102hz = target;
  for(int i=0; i<N; i++) line_102hz[i] = 0.2*sin((omega+omega_lf)*i*dt);

  // apply to witness channel a transfer function (coupling between witness and target channels)
  wavearray<double> witness_tf = witness;
  witness_tf.FFT(1);
  for(int i=0;i<N;i++) {
    double f = double(i)*df;
    if(f<=1) witness_tf[i]*=0.35*(pow(1-1,4)+1);
    if(f>1&&f<3) witness_tf[i]*=0.35*(pow(f-1,4)+1); 
    if(f>=3) witness_tf[i]*=0.35*(pow(3-1,4)+1);
  }
  witness_tf.FFT(-1);
  plot.plot(witness_tf, const_cast<char*>("same"), 2, totalscratch, totalscratch+LENGHT, true, 0, 5, true, 50);
  plot.graph[0]->SetTitle("Black : Witness Channel - Red : After Transfer Function");;

  // add to target channel a sine wave @ 100Hz + sine wave @ 102Hz + bicoherence 100Hz & 2Hz:
  for(int i=0; i<N; i++) target[i] += line_100hz[i] + line_102hz[i] + 50.*line_100hz[i]*witness_tf[i];
  // plot target channel
  c1->cd(2);
  gPad->SetLogy();
  plot.plot(target, const_cast<char*>("alp logy"), 1, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.graph[2]->SetTitle("Gauss + 100Hz Line + bicoherence 100Hz line & 2Hz line + 102Hz line");;

  //Make the target channel WDM transform, resolution = 0.5Hz
  int lev=target.rate()/2;
  lev*=2;
  WDM<double> wdtf(lev, 2*lev, 6, 10);
  WSeries<double> tfmap;
  tfmap.Forward(target, wdtf);

  //Adding target channel
  regression r;
  r.add(tfmap,const_cast<char*>("hchannel"));

  //Considering only the interested frequency band
  r.mask(0);
  r.unmask(0,F1,F2);

  //Adding witness channels
  r.add(target,const_cast<char*>("witness"),99.7,100.2);	// select target channel at range [99.7,100.2]
  r.add(witness,const_cast<char*>("witness_lf"),0,6);	// select witness channel at range [0,6]
  //Constructing channel describing multi-linear coupling
  r.add(1,2,const_cast<char*>("bic_witness"));
  //masking carrier channels and 
  r.mask(1);
  r.mask(2);

  //Calculate prediction
  r.setFilter(4);
  r.setMatrix(totalscratch);
  r.solve(0., 20, 's');
  r.apply(0.2);

  // Draw cleaned data 
  wavearray<double> clean = r.getClean();
  c1->cd(3);
  gPad->SetLogy();
  plot.plot(target, const_cast<char*>("alp logy"), 1, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.plot(clean, const_cast<char*>("same"), 2, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.graph[3]->SetTitle("Black : noisy data | Red : cleaned data");;

  // Draw cleaned data after subtraction of injected 100Hz & 102Hz lines
  c1->cd(4);
  gPad->SetLogy();
  plot.plot(clean, const_cast<char*>("alp logy"), 1, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  for(int i=0; i<N; i++) clean[i] -= line_100hz[i];
  for(int i=0; i<N; i++) clean[i] -= line_102hz[i];
  plot.plot(clean, const_cast<char*>("same"), 2, totalscratch, totalscratch+LENGHT, true, F1, F2, true, 50);
  plot.graph[5]->GetHistogram()->GetYaxis()->SetRangeUser(0.015,2.8);;
  plot.graph[5]->SetTitle("cleaned data, black: original, red: after subtraction of injected 100Hz & 102Hz lines");;

  // save plots
  c1->SaveAs("Regression_Sine_Gaus_Bic.png");

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
  for (int j=0; j<eigen[0].size(); j++) {
    for (int i=0; i<freq.size(); i++) printf("%.3f\t",eigen[i].data[j]);
    cout << endl;
  }

}
