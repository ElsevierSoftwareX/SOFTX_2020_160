//
// Test get detector info to/from tree
// Author : Gabriele Vedovato

{
  #define IFILE_NAME "../../cwbpipe/tutorials/S6A_R4_SIM_BRST_L1H1V1/data/wave_931081200_600_S6A_R4_SIM_BRST_L1H1V1_4.8_job8.root"

  TFile ifile(IFILE_NAME); 
  CWB::Toolbox::checkFile(IFILE_NAME);
  TTree* itree = (TTree *) gROOT->FindObject("waveburst");
  TList* list = itree->GetUserInfo();
  for (int n=0;n<list->GetSize();n++) {
    detector* pDetector = (detector*)list->At(n);
    detectorParams dParams = pDetector->getDetectorParams();
    if(n==0) pDetector->print();
    //cout << dParams.name << endl;
    //cout << pDetector->Name << endl;
  }
  ifile.Close();

  exit(0);
}

