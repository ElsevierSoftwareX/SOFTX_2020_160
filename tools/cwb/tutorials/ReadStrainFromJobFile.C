//
// Draw read & display strain & mdc from job strain file 
// Author : Gabriele Vedovato

{
  //#define DISPLAY_TIME
  #define DISPLAY_PSD
  //#define DISPLAY_MDC
  #define JOB_STRAIN_FILE "data/strain_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_MSA_job1.root"

  int nIFO = 3;
  TString ifo[3] = {"L1","H1","V1"}; 
  wavearray<double> hot[3];  //! temporary time series
  wavearray<double> x;
  wavearray<double>* px;
  double factor=10.*sqrt(nIFO);
  TString jname = JOB_STRAIN_FILE;

  // open input job file
  TFile* jfile = new TFile(jname);
  if(jfile==NULL||!jfile->IsOpen())
    {cout << "Error : file " << jname << " not found" <<  endl;exit(1);}

  for(int i=0; i<nIFO; i++) {
    // read ifo strain from temporary job file
    px = (wavearray<double>*)jfile->Get(TString("strain/")+ifo[i]);
    hot[i] = *px; delete px;
#ifdef DISPLAY_MDC
    hot[i]=0;	
#endif

    px = (wavearray<double>*)jfile->Get(TString("mdc/")+ifo[i]);
    (*px)*=factor;
    hot[i].add(*px);
    delete px;

    cout.precision(14);
    cout << hot[i].size() << " " << hot[i].rate() << " " << hot[i].start() << endl; 
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
