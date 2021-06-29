{
  //
  // Write & Read netcluster object to/from root file
  // Author : Gabriele Vedovato

  netcluster nc;

  netpixel pix;
  pix.frequency=100;

  nc.append(pix);

  TFile *froot = new TFile("test.root", "RECREATE");
  if(froot==NULL) {
    cout << "Failed to create file !!! " <<  endl;
    gSystem->Exit(1);
  }

  nc.Write("cluster");
  froot->Close();

  TFile *f = new TFile("test.root");
  if(f==NULL) {
    cout << "Failed to open file test.root !!! " <<  endl;
    gSystem->Exit(1);
  }

  f->ls();

  netcluster* nc2 = (netcluster*)f->Get("cluster");
  if(nc2==NULL) {
    cout << "Object cluster not exist !!! " <<  endl;
    gSystem->Exit(1);
  }

  netpixel* pix2 = nc2->getPixel(0,0);
  cout << "pix2 frequency " << pix2->frequency << endl;

  f->Close();

  exit(0);
}
