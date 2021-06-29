
void DrawCanvas(TString iFile) {

  TFile *ifile = TFile::Open(iFile.Data());
  if(ifile==NULL) {
    cout << "Open Error : " << iFile.Data() << endl;
    gSystem->Exit(1);
  }
  TCanvas* canvas = NULL;
  if(canvas==NULL) canvas = (TCanvas*)ifile->Get("WSeries");
  if(canvas==NULL) canvas = (TCanvas*)ifile->Get("WTS");
  if(canvas==NULL) {
    cout << "Canvas WSeries/WTS not present " << endl;
    gSystem->Exit(1);
  }

  canvas->Draw();
}
