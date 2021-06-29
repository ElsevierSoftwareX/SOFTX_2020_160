//
// Test set/get detector info to/from tree
// Author : Gabriele Vedovato

{

  #define N_IFO 3
  #define IFILE_NAME "test_streamer.root"

  char ifo[N_IFO][4]={"L1","H1","V1"};

  CWB::Toolbox::checkFile(IFILE_NAME);
  TFile ofile(IFILE_NAME,"RECREATE"); 
  TTree *otree = new TTree("tree", "tree");

  detector* pD[N_IFO];  

  for(int n=0;n<N_IFO;n++) {
    pD[n] = new detector(ifo[n]);
    otree->GetUserInfo()->Add(pD[n]);
  }
  otree->Write();
  ofile.Close();

  TFile ifile("test_streamer.root"); 
  TTree* itree = (TTree *) gROOT->FindObject("tree");
  TList* list = itree->GetUserInfo();
  for (int n=0;n<list->GetSize();n++) {
    detector* pDetector = (detector*)list->At(n);
    detectorParams dParams = pDetector->getDetectorParams();
    if(n==0) pDetector->print();
    cout << dParams.name << endl;
    //cout << pDetector->Name << endl;
cout << "----> DEB1" << endl;
    detector* D = new detector(*pDetector);
cout << "----> DEB2" << endl;
    D->print();
cout << "----> DEB3" << endl;
exit(0);
  }
  ifile.Close();

  detector* D = new detector(*pDetector);

  //detector V1("V1");
  //V1.print();

  exit(0);
}

