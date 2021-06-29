//
// This example show how to read/analyze skymap probability produced by CWB_Plugin_SkyProb.C plugin
// Author : Gabriele Vedovato

{

  #define COORDINATES "cWB"
  //#define COORDINATES "Geographic"
  //#define COORDINATES "Celestial"

  #define PROJECTION ""
  //#define PROJECTION "hammer"
  //#define PROJECTION "sinusoidal"

  #define EC_THR	0.5		// Euler Characteristic Threshold
  //#define DISPLAY_SKYMAP	// display skyprob skymap
  //#define ANTENNA_PATTERN	// tests with antenna pattern
  //#define EC_HIST		// display Euler Characteristic histogram
  //#define FILL_HOLE		// fill hole with the average of 8 neighbors pixels
  //#define SKYPROB_HIS		// display skyprob values @ injection dir histogram
  
  #include <vector>

  // read cwb parameters file
  CWB::Toolbox TB;
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TString cwb_parameter_name = gSystem->Getenv("CWB_PARAMETERS_FILE");
  gROOT->Macro(cwb_parameter_name);

  // define an alternative output_dir 
  strcpy(output_dir,"data");

  // read root file list from output directory
  vector<TString> fileList = TB.getFileListFromDir(output_dir, ".root", "wave_");
  if(fileList.size()==0) {
    cout << "No root files in the directory : " << output_dir << endl;
    gSystem->Exit(1);
  }

#ifdef EC_HIST
  TH1F* hec = new TH1F("hec","hec",50,0,10);
#endif
#ifdef SKYPROB_HIST
  TH1F* hskyprob = new TH1F("hskyprob","hskyprob",50,0,1);
#endif

  for(int n=0;n<fileList.size();n++) {		// loop over the root output files

    cout << n << " " << fileList[n].Data() << endl;
    netevent event(fileList[n].Data());      	// open wave root file
    int nevt = event.GetEntries();           	// get number of detected events
    if(nevt==0) continue;
    else cout << "nevt : " << nevt << endl;

    // get number of ifos (nIFO) and ifo name list
    TList* ifoList = event.fChain->GetUserInfo();
    detectorParams dParams[NIFO_MAX];
    int nIFO=ifoList->GetSize();
    TString IFO[NIFO_MAX];
    for (int k=0;k<nIFO;k++) {
      detector* pDetector = (detector*)ifoList->At(k);
      dParams[k] = pDetector->getDetectorParams();
      IFO[k]=dParams[k].name; 
      cout << dParams[k].name << endl;
      //pDetector->print(); 			// print detector infos
    }
    gnetwork* gNET = new gnetwork(nIFO,IFO);

    wavearray<int> index;

    for(int m=0;m<nevt;m++) {       	// loop over the detected events
      event.GetEntry(m);		// load event #n
      skymap* Psm = event.Psm;      	// probability skymap
      if(Psm->size()!=0) {
        if(Psm->getType()==1) {		// check if it is healpix
          if(fabs(event.time[0]-event.time[nIFO])<0.1) {        // filter output events

            gskymap* gSM = new gskymap(*(event.Psm));		

#ifdef FILL_HOLE	
            for(int l=0;l<Psm->size();l++) {
              index = Psm->neighbors(l);
              int M=0;double a=0;
              for(int k=0;k<index.size();k++) if(Psm->get(index[k])) {M++;a+=Psm->get(index[k]);}
              if(M) a/=M;		// average
              if(M==8) gSM->set(l,a);	// fill hole with the average of 8 neighbors pixels
            }
#endif
		 
            int ec = gSM->getEulerCharacteristic(EC_THR);	// compute the euler characteristic
            cout << n << " " << m << " EC -> " << ec << " erA[0] -> " << event.erA[0] << endl;
#ifdef EC_HIST
            hec->Fill(ec);
#endif

#ifdef DISPLAY_SKYMAP
            gStyle->SetPalette(1);
            gSM->SetOptions(PROJECTION,COORDINATES);
            char title[256];
            sprintf(title,"probability skymap gps = %d",int(event->time[0]));
            gSM->SetTitle(title);
            //gSM->GetHistogram()->GetZaxis()->SetRangeUser(TH, 1.);
            //for(int l=0;l<gSM->size();l++) if(gSM->get(l)<=EC_THR) gSM->set(l,0); else gSM->set(l,1);
            for(int l=0;l<gSM->size();l++) if(gSM->get(l)<=EC_THR) gSM->set(l,0); // set to 0 the pixels < EX_THR 

#ifdef ANTENNA_PATTERN
`	    // check the sky antenna pattern
            for(int l=0;l<gSM->size();l++) {
              double theta = gSM->getTheta(l);
              double phi   = gSM->getPhi(l);
              int polarization=2;	// |Fx|/|F+| DPF
              double ac = gNET->GetAntennaPattern(phi,theta,0,polarization);
              //gSM->set(l,ac); 
              //if(ac<0.2) gSM->set(l,ac); else gSM->set(l,0);
              if(ac<0.1) gSM->set(l,0); 
            }
#endif

            gSM->Draw();				// draw skymap
            double RA  = event->phi[1];
            double DEC = event->theta[1];
            gSM->DrawMarker(RA,DEC, 29, 1.5, kBlack);	// draw injected direction
            return;
#endif

#ifdef SKYPROB_HIST
            int l = gSM->getSkyIndex(event->theta[1],event->phi[1]);
            double skyprob = gSM->get(l) ;
            cout << "SkyProb @ injection pixel : " << l << " " << skyprob << endl;
            hskyprob->Fill(skyprob);
#endif

          }
        }
      }
    }
    delete gNET;
  }

#ifdef EC_HIST
  hec->Draw();
#endif
#ifdef SKYPROB_HIST
  hskyprob->Draw();
#endif
}
