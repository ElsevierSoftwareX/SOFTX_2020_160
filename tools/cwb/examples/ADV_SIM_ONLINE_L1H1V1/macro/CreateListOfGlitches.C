#define NUM_GLITCHES		2033		// number of glitches to be generated
#define GPS_START		1040739000.0	// start gps time

#define GLITCH_JITTER		10 		// sec
#define GLITCH_TIME_STEP	30 		// sec

#define SNR_GLITCH_LOW 		40		// low value of the glitch snr 
#define SNR_GLITCH_HIGH		40		// high value of the glitch snr

#define SEED			1234		// seed for random generation

{
  double rad2deg = 180./TMath::Pi();

  int nIFO=3;
  TString ifo[nIFO]={"L1","H1","V1"};
  network* net = new network();
  for(int n=0;n<nIFO;n++) net->add(new detector(const_cast<char*>(ifo[n].Data())));
  net->setRunID(SEED);				// run number is used as seed in configPlugin
 
  TString mdc_type="glitch"; 

  // load MDC setup  (define MDC set)
  TString IFO="";				// disable SetSkyDistribution in configPlugin
  gROOT->Macro("macro/CWB_Plugin_GOnlineFrame_Config.C");

  // create Glitch lists for each detector
  TString fname[nIFO];
  FILE *fp[nIFO];
  for(int k=0;k<nIFO;k++) { 
    fname[k]=TString("input/")+ifo[k]+"_glitches.lst";
    if((fp[k]=fopen(fname[k].Data(),"w")) == NULL ) {
      cout << "CreateListOfGlitches.C error : cannot open file " << fname[k].Data() << endl;
      gSystem->Exit(1);
    }
    cout << "open file " << fname[k].Data() << endl;
    fprintf(fp[k],"#          gps                       name       theta");
    fprintf(fp[k],"             phi              psi           rho            iota");
    fprintf(fp[k],"             snr		ID              id\n");
  }

  cout << endl;
  cout.precision(14);
  for(int n=0;n<NUM_GLITCHES;n++) {
    if(n%100==0) cout << "CreateListOfGlitches.C : " << n << "/" << NUM_GLITCHES << endl;
    double gps = GPS_START+n*GLITCH_TIME_STEP+gRandom->Uniform(-GLITCH_JITTER,GLITCH_JITTER);

    std::vector<int> vID;
    for (int i=0;i<MDC.wfList.size();i++) vID.push_back(i);

    for(int k=0;k<nIFO;k++) {
      double theta = rad2deg*acos(gRandom->Uniform(-1,1));
      double phi   = gRandom->Uniform(0,360);
      double psi   = gRandom->Uniform(0,180);
      double rho   = 1;
      double iota  = 0;	// inclination angle : the line of sight is perpendicular to the orbit plane
      double snr   = gRandom->Uniform(SNR_GLITCH_LOW,SNR_GLITCH_HIGH);

      // get a random waveform name from the MDC set
      int xID = (int)gRandom->Uniform(0,vID.size());
      int ID = vID[xID];
      int id = (int)gRandom->Uniform(0,1+MDC.wfList[ID].list.size());
      // erase the IDth element
      vID.erase(vID.begin()+xID);
      waveform wf = MDC.GetWaveform(ID,id);
      fprintf(fp[k],"%10.3f\t% 25s\t% 0.2f\t\t%3.2f\t\t%3.2f\t\t%3.2f\t\t%3.2f\t\t%3.2f\t\t%d\t\t%d\n",
              gps,wf.name.Data(),theta,phi,psi,rho,iota,snr,ID,id);
    }
  }
  for(int k=0;k<nIFO;k++) fclose(fp[k]);
  gSystem->Exit(0);
}
