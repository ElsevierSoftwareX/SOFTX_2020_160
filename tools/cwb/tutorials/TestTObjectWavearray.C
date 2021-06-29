{
  //
  // Write & Read wavearray object to/from root file
  // Author : Gabriele Vedovato


  wavearray<double> x1(1024*128);
  x1.rate(1024*4);
  for(int i=0;i<x1.size();i++) x1[i]=i;

  wavearray<double> x2(1024*128);
  x2.rate(1024*4);
  for(int i=0;i<x2.size();i++) x2[i]=i+1000;

  TFile *froot = new TFile("test.root", "RECREATE");
  if(froot==NULL) {
    cout << "Failed to create file !!! " <<  endl;
    gSystem->Exit(1);
  }

  x1.Write("wavearray<double>");
  x2.Write("wavearray<double>");
  froot->Close();

  TFile *f = new TFile("test.root");
  if(f==NULL) {
    cout << "Failed to open file test.root !!! " <<  endl;
    gSystem->Exit(1);
  }

  f->ls();

  wavearray<double>* w = (wavearray<double>*)f->Get("wavearray<double>;2");
  if(w==NULL) {
    cout << "Object wavearray not exist !!! " <<  endl;
    gSystem->Exit(1);
  }

  cout << w->size() << endl;
  for(int i=0;i<10;i++) cout << i << " " << w->data[i] << endl;

  f->Close();

  exit(0);
}
