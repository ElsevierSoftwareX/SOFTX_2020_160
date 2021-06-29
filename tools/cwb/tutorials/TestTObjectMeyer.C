{
  //
  // Write & Read Meyer object to/from root file
  // Author : Gabriele Vedovato


  Meyer<double> S(512,2);
  cout << "S level " << S.getLevel() << endl;

  TFile *froot = new TFile("test.root", "RECREATE");
  if(froot==NULL) {
    cout << "Failed to create file !!! " <<  endl;
    gSystem->Exit(1);
  }
  S.Write("S");
  froot->Close();

  TFile *f = new TFile("test.root");
  if(f==NULL) {
    cout << "Failed to open file test.root !!! " <<  endl;
    gSystem->Exit(1);
  }

  f->ls();

  Meyer<double>* S2 = (Meyer<double>*)f->Get("S");
  if(S2==NULL) {
    cout << "Object S not exist !!! " <<  endl;
    gSystem->Exit(1);
  }
  cout << "S2 level " << S2->getLevel() << endl;
  f->Close();

  exit(0);
}
