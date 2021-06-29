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

  Meyer<double> S(512,2);
  WSeries<double> w1(x1,S);
  WSeries<double> w2(x2,S);
  w2.Forward(2);
  cout << "S level " << S.getLevel() << endl;
  cout << "w1 level " << w1.getLevel() << endl;
  cout << "w2 level " << w2.getLevel() << endl;

  TFile *froot = new TFile("test.root", "RECREATE");
  S.Write("S");
  wavearray<double>(w1).Write("w1");
  w2.Write(TString("mdc")+"w2");
  froot->Close();

  TFile *f = new TFile("test.root");

  f->ls();

//  Meyer<double>* S2 = (Meyer<double>*)f->Get("Meyer<double>;1");
//  cout << "S2 level " << S2->getLevel() << endl;

  //WSeries<double>* w = (WSeries<double>*)f->Get("w2");
  //cout << "w level " << w->getLevel() << endl;
  //w->Inverse(2);

  //wavearray<double>* w = (wavearray<double>*)f->Get("w1");
  //cout << w->size() << endl;
  //for(int i=0;i<10;i++) cout << i << " " << w->data[i] << endl;

  for(int i=0;i<w2.size();i++) w2.data[i]=0;
  w2 = *(wavearray<double>*)f->Get(TString("mdc")+"w2");
  cout << "w2 level " << w2.getLevel() << endl;
  w2.Inverse(2);
  for(int i=0;i<10;i++) cout << i << " " << w2.data[i] << endl;

/*
  TIter nextkey(f.GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
     WSeries<double> *w = (WSeries<double>*)key->ReadObj();
     if(w==NULL) {cout << "Null object !!!" << endl;exit(1);}
     cout << w.size() << endl;
     for(int i=0;i<10;i++) cout << i << " " << w.data[i] << endl;
  }
*/

  f->Close();

  exit(0);
}
