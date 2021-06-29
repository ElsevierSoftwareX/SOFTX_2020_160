{
  //
  // Write & Read wavearray, Wseries objects to/from root file
  // Author : Gabriele Vedovato


  wavearray<double> x1(1024*128);
  x1.rate(1024*4);
  for(int i=0;i<x1.size();i++) x1[i]=i;

  Meyer<double> S(512,2);
  WSeries<double> w1(x1,S);
  w1.Forward(1);
  cout << "w1 level " << w1.getLevel() << endl;

  TFile *froot = new TFile("test.root", "RECREATE");
  if(froot==NULL) {
    cout << "Failed to create file !!! " <<  endl;
    gSystem->Exit(1);
  }

  w1.Write("WSeries<double>");
  froot->Close();

  TFile *f = new TFile("test.root");
  if(f==NULL) {
    cout << "Failed to open file test.root !!! " <<  endl;
    gSystem->Exit(1);
  }

  f->ls();

  WSeries<double>* w = (WSeries<double>*)f->Get("WSeries<double>;1");
  if(w==NULL) {
    cout << "Object WSeries not exist !!! " <<  endl;
    gSystem->Exit(1);
  }

  cout << w->size() << endl;
  cout << "w level " << w->getLevel() << endl;
  for(int i=0;i<10;i++) cout << i << " " << w->data[i] << endl;

  f->Close();

  exit(0);
}
