{
  //
  // Write & Read netpixel object to/from root file
  // Author : Gabriele Vedovato

  netpixel pix;
  pix.frequency=100;

  TFile *froot = new TFile("test.root", "RECREATE");
  if(froot==NULL) {
    cout << "Failed to create file !!! " <<  endl;
    gSystem->Exit(1);
  }

  pix.Write("netpixel");
  froot->Close();

  TFile *f = new TFile("test.root");
  if(f==NULL) {
    cout << "Failed to open file test.root !!! " <<  endl;
    gSystem->Exit(1);
  }

  f->ls();

  netpixel* pix2 = (netpixel*)f->Get("netpixel");
  if(pix2==NULL) {
    cout << "Object pixel not exist !!! " <<  endl;
    gSystem->Exit(1);
  }

  cout << "pix2 frequency " << pix2->frequency << endl;

  f->Close();

  exit(0);
}
