//
// Test set/get detector info to/from tree
// Author : Gabriele Vedovato

{

  #define N_IFO 9

                        // V1
  detectorParams dP[N_IFO] = {{"X1", 43.6314,  10.5045, 0.0, 0,   90, 0,   30},
                          {"X2", 43.6314,  10.5045, 0.0, 0,  150, 0, -150},
                          {"X3", 43.6314,  10.5045, 0.0, 0,  -90, 0,  -30},

                         // L1
                          {"Y1", 30.5629, -90.7742, 0.0, 0,   90, 0,   30},
                          {"Y2", 30.5629, -90.7742, 0.0, 0,  150, 0, -150},
                          {"Y3", 30.5629, -90.7742, 0.0, 0,  -90, 0,  -30},

                         // H1
                          {"Z1", 46.4551, -119.408, 0.0, 0,   90, 0,   30},
                          {"Z2", 46.4551, -119.408, 0.0, 0,  150, 0, -150},
                          {"Z3", 46.4551, -119.408, 0.0, 0,  -90, 0,  -30},
                         };

  TFile ofile("test_streamer.root","RECREATE"); 
  TTree *otree = new TTree("tree", "tree");

  detector* pD[N_IFO];  

  for(int n=0;n<N_IFO;n++) {
    pD[n] = new detector(dP[n]);
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
  }
  ifile.Close();

  //detector V1("V1");
  //V1.print();

  exit(0);
}

