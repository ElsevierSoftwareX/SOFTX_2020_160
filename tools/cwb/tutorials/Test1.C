//Test 1: Removal of broadband noise using a single, highly correlated witness channel.
//You can prepare the channels with the desired correlation this way, for instance:
//Generate time series w(t) and x(t) each containing Gaussian white noise (unit variance)
//Construct the witness channel: A(t) = w(t)
//Construct the target channel:  H(t) = x(t) + 0.8*w(t)
//Compare the ASD (or simply the RMS) of H(t) before and after regression.
//The regression should be able to reduce the noise amplitude by a factor of
//     1/sqrt(1^2+0.8^2) = 1/sqrt(1.64) ~ 0.78

{
  //Time
  #define LENGHT 1200
  #define SCRATCH 32
  #define RATE 2048

  using namespace CWB;

  //define x channel properties -> gauss 1
  int N = RATE*(LENGHT+2*SCRATCH);
  wavearray<double> x;
  x.rate(RATE);
  x.resize(N);
  x.start(0);
  x.stop(LENGHT+2*SCRATCH);

  //define w channel properties -> gauss 2
  wavearray<double> w = x;

  // time series is filled with white noise data: 
  TRandom3 rnd(0);   
  for(int i=0; i<N; i++) x[i] = rnd.Gaus();
  for(int i=0; i<N; i++) w[i] = rnd.Gaus();

  //Fill target and witness
  wavearray<double> A = w;
  wavearray<double> H = x; w*= 0.8; H += w;

  //Make WDM transform, resolution = 1Hz
  int lev=H.rate()/2;
  WDM<double> wdtf(lev, 2*lev, 6, 10);
  WSeries<double> tfmap;
  tfmap.Forward(H, wdtf);

  //Adding target channel
  regression r;
  r.add(tfmap,const_cast<char*>("hchannel"));

  //Adding witness channel
  r.add(A,const_cast<char*>("witness"));

  //Calculate prediction
  r.setFilter(5);
  r.setMatrix(SCRATCH);
  r.solve(0.5, 0, 'h');
  r.apply(0.2);

  wavearray<double> clean=r.getClean();
  cout << "Ratio rms: (" << clean.rms() << "/" << H.rms() << ")= " << clean.rms()/H.rms() << endl;

  cout << "x         : " << x.mean() << " " << x.rms() << endl;
  cout << "clean     : " << clean.mean() << " " << clean.rms() << endl;
  clean -= x;
  cout << "clean-x   : " << clean.mean() << " " << clean.rms() << endl;

  wavearray<double> eigen=r.getVEIGEN(-1);
  eigen.rate(11);
  watplot plot;
  TCanvas *c1 = plot.canvas;
  c1->Divide(1,2);

  c1->cd(1);
  plot.plot(eigen,const_cast<char*>("alp"),1);
  plot.graph[0]->SetTitle("Eigen-values of all layers");
   
  wavearray<double> freq=r.vfreq;
  wavearray<double> rank=r.getRank(0);
  c1->cd(2);
  TGraph* g=new TGraph(freq.size(),freq.data,rank.data);
  g->SetLineColor(1);
  g->Draw("alp");
  g->SetTitle("Ranking for all layers");
  
  c1->SaveAs("Test1.png");
  
  exit(0);
}
