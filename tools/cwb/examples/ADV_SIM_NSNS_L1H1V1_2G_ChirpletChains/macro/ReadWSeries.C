
void ReadWSeries(TString fname, TString ifo, int level) {

  gROOT->LoadMacro("macro/PlotWSeries.C");

  TFile froot(fname);	// open fname
  froot.ls();		// list objects in saved the file

  char label[32];sprintf(label,"%s:%d",ifo.Data(),level);
  WSeries<double>* WS = (WSeries<double>*)froot.Get(label);
  cout << WS->size() << endl;

  char ofname[32];sprintf(ofname,"wdm_%s_%d.root",ifo.Data(),level);
  PlotWSeries(WS,ofname); 

  exit(0);
}
