//
// This example show how to convert skymap probability from a list of output/root files to fits files
// Author : Gabriele Vedovato

{
  #include <vector>

  #define WAVE_FILE "data/wave_931158008_600_ADV_SIM_BRST_LF_GWGC_L1H1_1G_run3_14.1421_job1.root"
  #define NETRHO 4
  #define NETCC  0.6
  #define ODIR "fits"
  //#define PRINT_DETECTORS


  CWB::Toolbox TB;
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TString cwb_parameter_name = gSystem->Getenv("CWB_PARAMETERS_FILE");
  gROOT->Macro(cwb_parameter_name);

  // read root file list from output directory
  vector<TString> fileList = TB.getFileListFromDir(output_dir, ".root", "wave_");
  if(fileList.size()==0) {
    cout << "No root files in the directory : " << output_dir << endl;
    gSystem->Exit(1);
  }

  // create output directory
  char cmd[256];
  sprintf(cmd,"mkdir -p %s",ODIR);
  cout << cmd << endl;
  gSystem->Exec(cmd);

  for(int n=0;n<fileList.size();n++) {

    cout << fileList[n].Data() << endl;
    netevent event(fileList[n].Data());      // open wave root file
    int nevt = event.GetEntries();           // get number of detected events
    if(nevt==0) continue;
    else cout << "nevt : " << nevt << endl;

#ifdef PRINT_DETECTORS
    // print detectors infos
    TList* ifoList = event.fChain->GetUserInfo();
    int nIFO = ifoList->GetSize();
    for (int m=0;m<nIFO;m++) {
      detector* pDetector = (detector*)ifoList->At(m);
      detectorParams dParams = pDetector->getDetectorParams();
      cout << dParams.name << endl;
      pDetector->print();
    }
#endif

    for(int i=0;i<nevt;i++) {       // loop over the detected events

      event.GetEntry(i);		  // load event #n
 
      skymap* Psm = event.Psm;      // probability skymap

      cout << "--------------------------------------------------" << endl;
      cout << "gps         : " << int(event.time[0]) << endl;
      cout << "rho         : " << event.rho[1] << endl;
      cout << "cc          : " << event.netcc[0] << endl;
      cout << "network snr : " << sqrt(event.likelihood) << endl;
      cout << "--------------------------------------------------" << endl<<endl;

      if(Psm->size()!=0) {

        if(Psm->getType()==1) {		// check if it is healpix
          // dump skymap to fits (event.time[0] is the gps time of the first detector)
          if((event.rho[1]>NETRHO) && (event.netcc[0]>NETCC)) {    // filter output events
            char fits_name[256];sprintf(fits_name,"%s/probability_skymap_%d.fits.gz",ODIR,int(event.time[0]));
            cout << "--------------------------------------------------" << endl<<endl;
            cout << "Dump fits file : " << fits_name << endl;
            Psm->Dump2fits(fits_name,event.time[0]);   

            // compute total probability (Psm->get(l) return the probability for index l)
            int ml=0;
            double prob=0;   // probability
            double mprob=0;  // max probability
            double cprob=0;  // cumulative probability
            for(int l=0;l<Psm->size();l++) {prob=Psm->get(l);cprob+=prob;if(prob>mprob) {mprob=prob;ml=l;}}
            cout << "cprob : " << cprob << endl;
            cout << "mprob : " << mprob << " ml " << ml << endl;

            // print RA,Dec form maximum probability
            // mRA & event.phi[2] are not equal because of approximations
            double mRA  = Psm->getPhi(ml);
            double mDec = Psm->getTheta(ml);
            cout << "max prob RA      : " << mRA << " RA (root) : " << event.phi[2] << endl;
            cout << "max prob Dec     : " << 90-mDec << " Dec (root) : " << event.theta[2] << endl;
            cout << "--------------------------------------------------" << endl<<endl;
          }
        }
      }
    }
  }

  exit(0);
}
