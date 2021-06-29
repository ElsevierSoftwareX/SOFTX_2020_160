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
  x1.Write("X1");
  x2.Write("X2");
  x3.Write("X3");
  froot->Close();

  TFile *f = new TFile("gwavearray_test.root");

  f->ls();

  gwavearray<double>* w3 = (gwavearray<double>*)f->Get("X3");
  cout << w3->size() << endl;
  for(int i=0;i<10;i++) cout << i << " " << w3->data[i] << endl;
  w3->DrawFFT();
  //w3->Draw(GWAT_FFT);
  gwavearray<double>* w1 = (gwavearray<double>*)f->Get("X1");
  w1->DrawFFT();

  f->Close();

  //exit(0);
}
