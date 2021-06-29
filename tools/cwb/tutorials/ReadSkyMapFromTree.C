//
// This example show how to read/save to fits/display skymap probability from output root file 
// Author : Gabriele Vedovato

//#define COORDINATES "cWB"
//#define COORDINATES "Geographic"
#define COORDINATES "Celestial"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define WAVE_FILE "data.1/wave_966383954_600_uniform_in_snr_dbg1_1_job1.root"

#define XANALYSIS "2G"
#define XSEARCH   'r'

void ReadSkyMapFromTree() {

  CWB::Toolbox TB;

  // load configuration parameters
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TString cwb_parameter_name = gSystem->Getenv("CWB_PARAMETERS_FILE");
  gROOT->Macro(cwb_parameter_name);
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  TString cwb_uparameter_name = gSystem->Getenv("CWB_UPARAMETERS_FILE");
  gROOT->Macro(cwb_uparameter_name);

  netevent* event =  new netevent(WAVE_FILE);      // open wave root file
  int nevt = event->GetEntries();  // get number of detected events
  if(nevt==0) {
    cout << "no events are presents in the file : " << WAVE_FILE << endl;
    gSystem->Exit(1);
  } else cout << "nevt : " << nevt << endl;

  // print detectors infos
  TList* ifoList = event->fChain->GetUserInfo();
  int nIFO = ifoList->GetSize();
  for (int n=0;n<nIFO;n++) {
    detector* pDetector = (detector*)ifoList->At(n);
    detectorParams dParams = pDetector->getDetectorParams();
    cout << dParams.name << endl;
    pDetector->print();
  }

  for(int n=0;n<nevt;n++) {       // loop over the detected events

    event->GetEntry(n);		  // load event #n
 
    skymap* Psm = event->Psm;      // probability skymap

    cout << "-----------------------------------" << endl;
    cout << "gps         : " << int(event->time[0]) << endl;
    cout << "rho         : " << event->rho[1] << endl;
    cout << "cc          : " << event->netcc[0] << endl;
    cout << "network snr : " << sqrt(event->likelihood) << endl;
    cout << "rec RA      : " << event->phi[2] << endl;
    cout << "rec DEC     : " << event->theta[2] << endl;
    cout << "-----------------------------------" << endl<<endl;

    if(Psm->size()!=0) {

      if(Psm->getType()==1) {		// check if it is healpix
        // dump skymap to fits (event->time[0] is the gps time of the first detector)
        char fits_name[256];sprintf(fits_name,"probability_skymap_%d.fits.gz",int(event->time[0]));
        // build configur info 
        char configur[64];
        if (XSEARCH=='r')                 sprintf(configur,"%s un-modeled",XANALYSIS);
        if((XSEARCH=='i')||(XSEARCH=='I')) sprintf(configur,"%s elliptical",XANALYSIS);
        if((XSEARCH=='s')||(XSEARCH=='S')) sprintf(configur,"%s linear",XANALYSIS);
        if((XSEARCH=='g')||(XSEARCH=='G')) sprintf(configur,"%s circular",XANALYSIS);
        cout << "Dump fits file : " << fits_name << endl;
        Psm->Dump2fits(fits_name,event->time[0],configur);   
      }

      // plot probability skymap 
      gskymap* gSM = new gskymap(*(event->Psm));
      gSM->SetOptions(PROJECTION,COORDINATES);
      char title[256];sprintf(title,"probability skymap gps = %d",int(event->time[0]));
      gSM->SetTitle(title);
      gSM->Draw();
      double RA  = event->phi[2];
      double DEC = event->theta[2];
      RA=360-RA;
      gSM->DrawMarker(RA,DEC, 29, 1.5, kWhite);  

      double prob=0;
      for(int l=0;l<Psm->size();l++) prob+=Psm->get(l);
      cout << "prob : " << prob << endl;

    }

    break;  // terminate after the first event
  }
}
