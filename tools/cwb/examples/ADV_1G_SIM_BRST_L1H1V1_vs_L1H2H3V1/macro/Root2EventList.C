void Root2EventList(TString wFile, TString oFileName) {

  netevent event(wFile.Data());    // open wave root file
  int nevt = event.GetEntries();   // get number of detected events
  if(nevt==0) {
    cout << "no events are presents in the file : " << WAVE_FILE << endl;
    gSystem->Exit(1);
  } else cout << "nevt : " << nevt << endl;

  // print detectors infos
  TList* ifoList = event.fChain->GetUserInfo();
  int nIFO = ifoList->GetSize();
  for (int n=0;n<nIFO;n++) {
    detector* pDetector = (detector*)ifoList->At(n);
    detectorParams dParams = pDetector->getDetectorParams();
    cout << dParams.name << endl;
    pDetector->print();
  }


  ofstream out;
  out.open(oFileName.Data(), ios::out);
  if (!out.good()) {cout << "Error Opening File : " << oFileName.Data() << endl;exit(1);}
  cout << "Write file : " << oFileName.Data() << endl;
  out << "#\t" << "gps\t" << "name\t" <<  "theta\t" <<  "phi\t" <<  "psi\t" << "rho\t" << "iota\t" << "hrss" << endl;
  out.precision(14);

  for(int n=0;n<nevt;n++) {       // loop over the detected events

    event.GetEntry(n);            // load event #n

    double gps   = event.time[nIFO];   
    TString name = "SG849Q8d9";   
    float theta  = event.theta[1]; 
    float phi    = event.phi[1];    
    float psi    = event.psi[1];    
    float rho    = 0;     
    float iota   = 0;   
    double hrss  = event.strain[1];

    out << gps << "\t" << name << "\t" << theta << "\t" << phi << "\t" << psi << "\t" << rho << "\t" << iota << "\t" << hrss;

  }

  out.close();

  gSystem->Exit(0);
}
