{
  //
  // Write & Read skymap object to/from root file
  // Author : Gabriele Vedovato

  skymap sm(3);

  for(int l=0;l<sm.size();l++) sm.set(l,l);

  TFile *froot = new TFile("test.root", "RECREATE");
  if(froot==NULL) {
    cout << "Failed to create file !!! " <<  endl;
    gSystem->Exit(1);
  }

  sm.Write("skymap");
  froot->Close();

  TFile *f = new TFile("test.root");
  if(f==NULL) {
    cout << "Failed to open file test.root !!! " <<  endl;
    gSystem->Exit(1);
  }

  f->ls();

  skymap* sm2 = (skymap*)f->Get("skymap");
  if(sm2==NULL) {
    cout << "Object skymap not exist !!! " <<  endl;
    gSystem->Exit(1);
  }

  for(int l=0;l<sm2->size();l++) cout << l << " " << sm2->get(l) << endl;

  f->Close();

  exit(0);
}
