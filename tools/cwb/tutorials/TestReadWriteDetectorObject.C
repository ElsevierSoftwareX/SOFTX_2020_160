{
  //
  // Write & Read config object to/from root file
  // Author : Gabriele Vedovato

  detector D(const_cast<char*>("H1"));
  D.setPolarization(SCALAR);

  detector D1=D;

  TFile *froot = new TFile("test.root", "RECREATE");
  D1.Write(const_cast<char*>("D"));
  froot->Close();

  TFile *f = new TFile("test.root");

  f->ls();

  detector *iD = (detector*)f->Get(const_cast<char*>("D"));
  iD->print();

  f->Close();

  exit(0);
}
