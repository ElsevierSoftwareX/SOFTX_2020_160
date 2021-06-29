//
// This example show how to convert skymap probability from root to fits files in celestial coordinates
// Author : Gabriele Vedovato

{
  #define WAVE_FILE "merge/wave_bbh_low_lh_wp10_run2.M1.C_salvo.root"
  #define NETRHO 0
  #define NETCC  0.0
  #define ODIR "fits"

  netevent event(WAVE_FILE);      // open wave root file
  int nevt = event.GetEntries();  // get number of detected events
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

  for(int n=0;n<nevt;n++) {       // loop over the detected events

    event.GetEntry(n);		  // load event #n
 
    skymap* Psm = event.Psm;      // probability skymap

    cout << "-----------------------------------" << endl;
    cout << "run         : " << event.run << endl;
    cout << "gps         : " << int(event.time[0]) << endl;
    cout << "rho         : " << event.rho[1] << endl;
    cout << "cc          : " << event.netcc[0] << endl;
    cout << "network snr : " << sqrt(event.likelihood) << endl;
    cout << "-----------------------------------" << endl<<endl;

    if(Psm->size()!=0) {

      if(Psm->getType()==1) {		// check if it is healpix
        // dump skymap to fits (event.time[0] is the gps time of the first detector)
        if((event.rho[1]>NETRHO) && (event.netcc[0]>NETCC)) {    // filter output events

          // Dump2fits probability skymap in celestial coordinates (healpix)
          skymap skyprobcc = *Psm;
          skyprobcc=0.;

          double th,ph,ra;
          int k;
          for(j=0; j<int(Psm->size()); j++) {
            th = Psm->getTheta(j);
            ph = Psm->getPhi(j);

            k=Psm->getSkyIndex(th, ph);

            ra = Psm->getRA(j);
            k=Psm->getSkyIndex(th, ra);
            skyprobcc.set(k,Psm->get(j));
          }

          char fits_name[256];sprintf(fits_name,"%s/probability_skymap_%d.fits.gz",ODIR,int(event.time[0]));
          cout << "Dump fits file : " << fits_name << endl;
          skyprobcc.Dump2fits(fits_name,event.time[0],"2G:MRA un-modeled",const_cast<char*>("PROB"),const_cast<char*>("pix-1"),'C');   

          double prob=0;
          for(int l=0;l<skyprobcc.size();l++) prob+=skyprobcc.get(l);
          cout << "prob : " << prob << endl;
        }
      }
    }
  }

  exit(0);
}
