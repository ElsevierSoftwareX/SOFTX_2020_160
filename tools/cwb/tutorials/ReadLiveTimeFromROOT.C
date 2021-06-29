//
// Read livetime tree from input root file 
// Author : Gabriele Vedovato
//
// Example: root -l -b 'ReadLiveTimeFromROOT.C("merge/live_O2_K02_C02c_LH_BBH_BKG_run1.M1.root",2)'
//

void ReadLiveTimeFromROOT(TString liveFileName, int nIFO) {
//
// liveFileName: input livetime root file
// nIFO        : number of detectors
//

  TFile *_file0 = TFile::Open(liveFileName);

  TTree* tree = (TTree *) gROOT->FindObject("liveTime");
  if(tree==NULL) {cout << "ScanLIVE : liveTime tree not found !!!" << endl;exit(1);}

  char selection[1024];
  sprintf(selection,"lag[%d]==0 && slag[%d]==0",nIFO,nIFO);

  tree->Draw("start[0]:stop[0]:live",selection,"goff");	// select zero lag segments

  int size = (Int_t)tree->GetSelectedRows();		// get the selected entries
  double* start = tree->GetV1();
  double* stop = tree->GetV2();
  double* live = tree->GetV3();

  double live_dq2=0;
  double live_dq1=0;
  cout.precision(14);
  for(int i=0;i<size;i++)  {
    cout << i << "\tstart: " << start[i] << "\tstop: " << stop[i] << endl;
    live_dq1+=stop[i]-start[i];
    live_dq2+=live[i];
  }

  cout << "livetime zero lag after DQ CAT1         : " << live_dq1 << " sec\t" << live_dq1/(24*3600.) << " days" << endl;
  cout << "livetime zero lag after DQ CAT2 & GATING: " << live_dq2 << " sec\t" << live_dq2/(24*3600.) << " days" << endl;

  exit(0);
}
