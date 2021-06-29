{
  //
  // Write & Read WaveDMT object to/from root file
  // Author : Gabriele Vedovato


  WaveDWT<double> S(512,2);
  S.setLevel(2);
  cout << "S " << S.getLevel() << endl;

  TFile *froot = new TFile("test.root", "RECREATE");
  if(froot==NULL) {
    cout << "Failed to create file !!! " <<  endl;
    gSystem->Exit(1);
  }

  S.Write("WaveDWT<double>");
  froot->Close();

  TFile *f = new TFile("test.root");
  if(f==NULL) {
    cout << "Failed to open file test.root !!! " <<  endl;
    gSystem->Exit(1);
  }

  f->ls();

  WaveDWT<double>* S2 = (WaveDWT<double>*)f->Get("WaveDWT<double>;1");
  if(S2==NULL) {
    cout << "Object WaveDWT not exist !!! " <<  endl;
    gSystem->Exit(1);
  }

  cout << "S2 " << S2->getLevel() << endl;

  f->Close();

  exit(0);
}
