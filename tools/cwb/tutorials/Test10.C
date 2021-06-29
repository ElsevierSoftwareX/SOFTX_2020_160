//Test 10: Two good witness channels plus one slightly correlated one
//Witness channels: A(t) = w(t)
//                  B(t) = v(t)
//                  C(t) = y(t) + 0.5*z(t)
//Target channel:   H(t) = x(t) + 0.8*w(t) + 0.6*v(t) + z(t)
//Compare the ASD (or simply the RMS) of H(t) before and after regression.
//Although C(t) is correlated, using it in the regression would add more noise than
//it removes.  So hopefully the code will use A and B to regress but ignore C.
//The result, then, should be the same as Test 9.
//* Try this with hard, soft, and mild regulators * 

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
  wavearray<double> v = x;
  wavearray<double> y = x;
  wavearray<double> z = x;

  // time series is filled with white noise data: 
  TRandom3 rnd(0);   
  for(int i=0; i<N; i++) x[i] = rnd.Gaus();
  for(int i=0; i<N; i++) w[i] = rnd.Gaus();
  for(int i=0; i<N; i++) v[i] = rnd.Gaus();
  for(int i=0; i<N; i++) y[i] = rnd.Gaus();
  for(int i=0; i<N; i++) z[i] = rnd.Gaus();

  //Fill target and witness
  wavearray<double> A = w;
  wavearray<double> B = v;
  wavearray<double> H = x; w*= 0.8; H += w; v*= 0.6; H += v; H += z;
  wavearray<double> C = y; z*= 0.5; C += z;

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
  r.add(B,const_cast<char*>("witness"));
  r.add(C,const_cast<char*>("witness"));

  regression r_s = r;  
  regression r_m = r;  

  //Calculate prediction
  r.setFilter(5);
  r.setMatrix(SCRATCH);
  r.solve(0, 20, 'h');
  r.apply(0.2);

  wavearray<double> clean=r.getClean();

  cout << "HARD:" << endl;
  cout << "Ratio rms: (" << clean.rms() << "/" << H.rms() << ")= " << clean.rms()/H.rms() << endl;

  cout << "x         : " << x.mean() << " " << x.rms() << endl;
  cout << "clean     : " << clean.mean() << " " << clean.rms() << endl;
  clean -= x;
  cout << "clean-x   : " << clean.mean() << " " << clean.rms() << endl;
  cout << "-------------------" << endl;

  //Calculate prediction
  r_s.setFilter(5);
  r_s.setMatrix(SCRATCH);
  r_s.solve(0, 20, 's');
  r_s.apply(0.2);

  clean=r_s.getClean();

  cout << "SOFT:" << endl;
  cout << "Ratio rms: (" << clean.rms() << "/" << H.rms() << ")= " << clean.rms()/H.rms() << endl;

  cout << "x         : " << x.mean() << " " << x.rms() << endl;
  cout << "clean     : " << clean.mean() << " " << clean.rms() << endl;
  clean -= x;
  cout << "clean-x   : " << clean.mean() << " " << clean.rms() << endl;
  cout << "-------------------" << endl;

  //Calculate prediction
  r_m.setFilter(5);
  r_m.setMatrix(SCRATCH);
  r_m.solve(0, 20, 'm');
  r_m.apply(0.2);

  clean=r_m.getClean();

  cout << "MILD:" << endl;
  cout << "Ratio rms: (" << clean.rms() << "/" << H.rms() << ")= " << clean.rms()/H.rms() << endl;

  cout << "x         : " << x.mean() << " " << x.rms() << endl;
  cout << "clean     : " << clean.mean() << " " << clean.rms() << endl;
  clean -= x;
  cout << "clean-x   : " << clean.mean() << " " << clean.rms() << endl;
  cout << "-------------------" << endl;

  wavearray<double> eigen=r.getVEIGEN(-1);
  eigen.rate(22);
  watplot plot;
  TCanvas *c1 = plot.canvas;
  c1->SetCanvasSize(1600,800);
  c1->Divide(1,2);

  c1->cd(1);
  plot.plot(eigen,const_cast<char*>("alp"),1);
  plot.graph[0]->SetTitle("Eigen-values of all layers");

  c1->cd(2);
  TPad* P=(TPad*)c1->GetPad(2);
  P->Divide(3,1);

  wavearray<double> freq=r.vfreq;
  wavearray<double> rank=r.getRank(0);
  wavearray<double> rank1=r.getRank(1);
  wavearray<double> rank2=r.getRank(2);
  wavearray<double> rank3=r.getRank(3);
  P->cd(1);
  TMultiGraph* mg=new TMultiGraph();
  TGraph* g=new TGraph(freq.size(),freq.data,rank.data);
  g->SetLineColor(1);
  mg->Add(g);
  TGraph* g1=new TGraph(freq.size(),freq.data,rank1.data);
  g1->SetLineColor(2);
  mg->Add(g1);
  TGraph* g2=new TGraph(freq.size(),freq.data,rank2.data);
  g2->SetLineColor(3);
  mg->Add(g2);
  TGraph* g3=new TGraph(freq.size(),freq.data,rank3.data);
  g3->SetLineColor(4);
  mg->Add(g3);
  mg->Draw(const_cast<char*>("alp"));
  mg->SetTitle("HARD: Ranking for all layers");

  wavearray<double> freq_s=r_s.vfreq;
  wavearray<double> rank_s=r_s.getRank(0);
  wavearray<double> rank1_s=r_s.getRank(1);
  wavearray<double> rank2_s=r_s.getRank(2);
  wavearray<double> rank3_s=r_s.getRank(3);
  P->cd(2);
  TMultiGraph* mg_s=new TMultiGraph();
  TGraph* g_s=new TGraph(freq_s.size(),freq_s.data,rank_s.data);
  g_s->SetLineColor(1);
  mg_s->Add(g_s);
  TGraph* g1_s=new TGraph(freq_s.size(),freq_s.data,rank1_s.data);
  g1_s->SetLineColor(2);
  mg_s->Add(g1_s);
  TGraph* g2_s=new TGraph(freq_s.size(),freq_s.data,rank2_s.data);
  g2_s->SetLineColor(3);
  mg_s->Add(g2_s);
  TGraph* g3_s=new TGraph(freq_s.size(),freq_s.data,rank3_s.data);
  g3_s->SetLineColor(4);
  mg_s->Add(g3_s);
  mg_s->Draw(const_cast<char*>("alp"));
  mg_s->SetTitle("SOFT: Ranking for all layers");

  wavearray<double> freq_m=r_m.vfreq;
  wavearray<double> rank_m=r_m.getRank(0);
  wavearray<double> rank1_m=r_m.getRank(1);
  wavearray<double> rank2_m=r_m.getRank(2);
  wavearray<double> rank3_m=r_m.getRank(3);
  P->cd(3);
  TMultiGraph* mg_m=new TMultiGraph();
  TGraph* g_m=new TGraph(freq_m.size(),freq_m.data,rank_m.data);
  g_m->SetLineColor(1);
  mg_m->Add(g_m);
  TGraph* g1_m=new TGraph(freq_m.size(),freq_m.data,rank1_m.data);
  g1_m->SetLineColor(2);
  mg_m->Add(g1_m);
  TGraph* g2_m=new TGraph(freq_m.size(),freq_m.data,rank2_m.data);
  g2_m->SetLineColor(3);
  mg_m->Add(g2_m);
  TGraph* g3_m=new TGraph(freq_m.size(),freq_m.data,rank3_m.data);
  g3_m->SetLineColor(4);
  mg_m->Add(g3_m);
  mg_m->Draw(const_cast<char*>("alp"));
  mg_m->SetTitle("MILD: Ranking for all layers");

  c1->SaveAs("Test10.png");


  exit(0);

}
