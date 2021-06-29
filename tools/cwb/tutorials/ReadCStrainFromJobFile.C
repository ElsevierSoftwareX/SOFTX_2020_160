//
// read & display conditioned strain from job cstrain file
// Author : Gabriele Vedovato

{
  //#define DISPLAY_TIME
  #define DISPLAY_PSD
  #define JOB_CSTRAIN_FILE "data/cstrain_XXXX_jobYY.root"	// define cstrain root file
  #define FACTOR_INDEX 0

  #define SCRATCH_TIME	10					// is the segEdge parameter (sec)

  int nIFO = 2;
  TString ifo[2] = {"L1","H1"}; 
  wavearray<double> hot[2];  //! temporary time series
  wavearray<double> x;
  wavearray<double>* px;
  double factor=10.*sqrt(nIFO);
  TString jname = JOB_CSTRAIN_FILE;

  // open input job file
  TFile* jfile = new TFile(jname);
  if(jfile==NULL||!jfile->IsOpen())
    {cout << "Error : file " << jname << " not found" <<  endl;exit(1);}

  char cstrain_dir[32];sprintf(cstrain_dir,"cstrain/cstrain-f%d",FACTOR_INDEX);
  for(int i=0; i<nIFO; i++) {
    // read ifo strain from temporary job file
    px = (wavearray<double>*)jfile->Get(TString(cstrain_dir)+"/"+ifo[i]);
    if(px==NULL) {cout<<"Error : "<<cstrain_dir<<" not present"<<endl;exit(1);}
    hot[i] = *px; delete px;

    cout.precision(14);
    cout << hot[i].size() << " " << hot[i].rate() << " " << hot[i].start() << endl; 

    int jS = SCRATCH_TIME*hot[i].rate();
    int jE = hot[i].size()-jS;
    cout << "DEB SIZE " << hot[i].size() << " " << hot[i].rate() << " " << int(hot[i].start()) << endl;

    double rms=0;
    for(int j=jS; j<jE; j++)  {
      //if(fabs(hot[i].data[j])>1) cout << "DEB0 " << i << " " << n << " " << hot[i].data[j] << endl;
      rms+=pow(hot[i].data[j],2);
    }
    rms/=(jE-jS);
    cout << "DEB RMS = " << sqrt(rms) << endl;

  }

#if defined (DISPLAY_TIME) || defined (DISPLAY_PSD)
  for(int n=0;n<nIFO;n++) hot[n].start(0);
  watplot plot(const_cast<char*>("plot"),200,20,800,500);
  char gtitle[256];
  sprintf(gtitle,"Network : ");
  for(int n=0;n<nIFO;n++) {
    if(n==0) sprintf(gtitle,"%s %s (black) ", gtitle, ifo[n].Data());
    if(n==1) sprintf(gtitle,"%s %s (red)   ", gtitle, ifo[n].Data());
    if(n==2) sprintf(gtitle,"%s %s (green) ", gtitle, ifo[n].Data());
  }

#ifdef DISPLAY_TIME
  plot.gtitle(gtitle,"time (sec)","strain");
  plot.goptions(const_cast<char*>("alp"), 1, 0., 0.);
#endif

#ifdef DISPLAY_PSD
  double flow=64;
  double fhigh=2048;
  double tstart = hot[0].start();
  double tstop  = tstart+hot[0].size()/hot[0].rate();
  plot.gtitle(gtitle,"frequency (Hz)","strain/#sqrt{Hz}");
  plot.goptions(const_cast<char*>("alp logx logy"), 1, tstart, tstop, true, flow,fhigh, true, 4);
#endif

  for(int n=0;n<nIFO;n++) hot[n] >> plot;	
#endif
  
}
