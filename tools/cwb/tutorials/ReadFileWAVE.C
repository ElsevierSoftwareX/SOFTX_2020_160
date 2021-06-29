//
// This example show how to read event parameters from root wave file
// The full parameter list is defined in trunk/wat/netevent.hh
// Author : Gabriele Vedovato


void ReadFileWAVE(TString fName, TString cuts="", int nIFO=0) {

  netevent event(fName.Data(),nIFO);   // open wave root file
  int nevt = event.GetEntries();       // get number of detected events
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

  // check cuts
  TTreeFormula* treeformula=NULL;
  if(cuts!="") {
    treeformula = new TTreeFormula("cuts", cuts.Data(), event.fChain);
    int err = treeformula->Compile(cuts.Data());
    if(err) {
      cout << "ReadFileWAVE.C - wrong input cuts " << cuts << endl;
      return -1;
    }
  }

  for(int n=0;n<nevt;n++) {       		// loop over the detected events

    event.GetEntry(n);		  		// load event #n
    
    if(treeformula!=NULL) {
      if(treeformula->EvalInstance()==0) {	// skip entry if it does not satisfy selection cuts
        cout << "Skip entry number : " << n << endl;
        continue;
      }
    }

    cout.precision(14);  
    cout << "-----------------------------------" << endl;
    cout << " Event Parameters             " << n << endl;
    cout << "-----------------------------------" << endl;
    cout << endl;
    cout << "rho         : " << "rec = " << event.rho[0] << endl;
    cout << "cc          : " << "rec = " << event.netcc[1] << endl;
    cout << "subnet      : " << "rec = " << event.netcc[3] << endl;
    cout << "network snr : " << "rec = " << sqrt(event.likelihood) << endl;
    cout << "frequency   : " << "rec = " << event.frequency[0] 
                             << " fLow = " << event.low[0] << " fHigh = " << event.high[0] << endl;
    cout << "phi         : " << "rec = " << event.phi[0] << " inj = " << event.phi[1] << endl;
    cout << "theta       : " << "rec = " << event.theta[0] << " inj = " << event.theta[1] << endl;
    cout << "error region: " << event.erA[0] << endl;
    cout << endl;
    for(int i=0;i<nIFO;i++) {
      cout << "-----------------------------------" << endl;
      if(ifoList->GetSize()==nIFO)
        cout << "ifo : " << dParams[i].name << endl;
      else
        cout << "ifo : " << i << endl;
      cout << "-----------------------------------" << endl;
      cout << endl;
      cout << "time        : " << "rec = " << event.time[i] << " inj = " << event.time[nIFO+i] << endl;
      cout << "time range  : " << "tStart = " << event.start[i] << " tStop = " << event.stop[i] << endl;
      cout << "snr         : " << event.snr[i] << endl; 
      cout << "hrss        : " << event.hrss[i] << endl; 
      cout << endl;
    }
  }

  exit(0);
}
