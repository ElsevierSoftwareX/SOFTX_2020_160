//
// This example show how to read event parameters from root mdc file
// The full parameter list is defined in trunk/wat/injection.hh
// Author : Gabriele Vedovato


void ReadFileMDC(TString fName, int nIFO=0) {

  injection event(fName.Data(),nIFO);   // open wave root file
  int nevt = event.GetEntries();        // get number of detected events
  if(nevt==0) {
    cout << "events are not present in the file : " << fName.Data() << endl;
    gSystem->Exit(1);
  } else cout << "nevt : " << nevt << endl;

  // if detector objects are present then print detector infos
  TList* ifoList = event.fChain->GetUserInfo();
  detectorParams dParams[NIFO_MAX];
  for (int n=0;n<ifoList->GetSize();n++) {
    detector* pDetector = (detector*)ifoList->At(n);
    dParams[n] = pDetector->getDetectorParams();
    cout << dParams[n].name << endl;
    pDetector->print();
  }

  // if nIFO is not provided (nIFO=0) then it is retrieved from the input root file
  nIFO = nIFO==0 ? ifoList->GetSize() : nIFO;
  if(nIFO==0) {
    cout << "nIFO is not contained in the root file, must be declared as the input macro parameter" << endl;
    gSystem->Exit(1);
  }

  for(int n=0;n<nevt;n++) {       // loop over the detected events

    event.GetEntry(n);		  // load event #n

    cout.precision(14);  
    cout << "-----------------------------------" << endl;
    cout << " Event Injection Parameters   " << n << endl;
    cout << "-----------------------------------" << endl;
    cout << endl;
    cout << "mdc index   : " << event.type << endl;  
    cout << "phi         : " << event.phi[0] << endl;
    cout << "theta       : " << event.theta[0] << endl;
    cout << endl;
    for(int i=0;i<nIFO;i++) {
      cout << "-----------------------------------" << endl;
      if(ifoList->GetSize()==nIFO)
        cout << "ifo : " << dParams[i].name << endl;
      else
        cout << "ifo : " << i << endl;
      cout << "-----------------------------------" << endl;
      cout << endl;
      cout << "time        : " << event.time[i] << endl;
      cout << "hrss        : " << event.hrss[i] << endl; 
      cout << endl;
    }
  }

  exit(0);
}
