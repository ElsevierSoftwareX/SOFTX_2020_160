{
  //
  // Write & Read gwavearray object to/from root file
  // Author : Gabriele Vedovato

  gwavearray<double> x1(1024*4);
  x1.rate(1024*4);
  double dt=1/x1.rate();
  for(int i=0;i<x1.size();i++) x1[i]=sin(2*PI*100*dt*i);

  wavearray<double> x2(1024*4);
  x2.rate(1024*4);
  for(int i=0;i<x2.size();i++) x2[i]=sin(2*PI*100*dt*i);

  gwavearray<double> x3(&x2);

  TFile *froot = new TFile("gwavearray_test.root", "RECREATE");
x2.resize(1024);
  x2.Write("clusters-lev:1-lag:678");
x2.resize(2*1024);
  x2.Write("clusters-lev:2-lag:678");
x2.resize(4*1024);
  x2.Write("clusters-lev:3-lag:678");
  froot->Close();

  TFile *f = new TFile("gwavearray_test.root");

  f->ls();

  gwavearray<double>* w1 = (gwavearray<double>*)f->Get("clusters-lev:1-lag:678");
  cout << w1->size() << endl;

  gwavearray<double>* w2 = (gwavearray<double>*)f->Get("clusters-lev:2-lag:678");
  cout << w2->size() << endl;

  gwavearray<double>* w3 = (gwavearray<double>*)f->Get("clusters-lev:3-lag:678");
  cout << w3->size() << endl;

  f->Close();

  //exit(0);
}
