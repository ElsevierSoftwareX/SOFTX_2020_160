{
  //
  // Write & Read config object to/from root file
  // Author : Gabriele Vedovato

  skymap sm(int(7));

  skymap sm1=sm;

  TFile *froot = new TFile("test.root", "RECREATE");
  sm1.Write("sm1");
  froot->Close();

  TFile *f = new TFile("test.root");

  f->ls();

  skymap *ism = (skymap*)f->Get("sm1");

  f->Close();

  exit(0);
}
