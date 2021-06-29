CWB::mdc* MDC;

void DrawSGE() {
  //
  // Show how to use mdc class to get & draw SG elliptical waveforms
  // Author : Gabriele Vedovato

  #define N_IFO 3

  TString ifo[N_IFO]={"L1","H1","V1"};
  MDC = new CWB::mdc(N_IFO,ifo);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par;

  par.resize(2);

  par[0].name="frequency"; par[0].value=100.;
  par[1].name="Q"; par[1].value=8.9;
  MDC->AddWaveform(MDC_SGE, par);

  MDC->Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjHrss(2.5e-21);
  MDC->SetInjRate(0.02);
  MDC->SetInjJitter(10.0);

  // -----------------------------------------------------------------------------------
  // define sky distribution : fixed earth/celestial direction & ellipticity : iota = 45
  // -----------------------------------------------------------------------------------
  float DEC  = -30.0;		// [-90:90]
  float RA   = 80.0;		// [0:360]
  float PSI  = 0.0;		// [0:180]
  float IOTA = 45;		// 0/90	-> circular/linear
  float GPS  = 10000000;
  float seed = 0;
  par.resize(6);
  par[0].name="theta";   par[0].value=DEC;
  par[1].name="phi";     par[1].value=RA;
  par[2].name="psi";     par[2].value=PSI;
  par[3].name="entries"; par[3].value=1000;
  par[4].name="iota";    par[4].value=IOTA;     // ellipticity [0:180] deg
                                                // if<0 || >180 -> random
  par[5].name="gps";     par[5].value=GPS;      // fixed injection time

  //MDC->SetSkyDistribution(MDC_CELESTIAL_FIX,par,seed); // theta,phi are in celestial coordinates	
  MDC->SetSkyDistribution(MDC_EARTH_FIX,par,seed);	 // theta,phi are in earth coordinates

  // Get the first waveform hp,hx components 
  wf = MDC->GetWaveform("SGE100Q8d9",0);
  if((wf.hp.size()==0)||(wf.hx.size()==0)) {
    cout << "Error : Waveform not present in the MDC set !!!" << endl;
    gSystem->Exit(1); 
  }

  cout << "hp : size : " << wf.hp.size() << " rate : " << wf.hp.rate() << " start : " << (int)wf.hp.start() << endl;
  cout << "hx : size : " << wf.hx.size() << " rate : " << wf.hx.rate() << " start : " << (int)wf.hx.start() << endl;
  wf.hp.start(0);		// set start to 0 (needed by draw Method)
  wf.hx.start(0);

  // uncomment the following lines to draw hp,hx components
  
  //MDC->Draw(wf.hp,MDC_TIME,"ALP ZOOM");
  //MDC->Draw(wf.hx,MDC_TIME,"same",kRed);

  //MDC->Draw(wf.hp,MDC_FFT,"ALP ZOOM"); // draw hp in frequency domain
  //MDC->Draw(wf.hp,MDC_TF,"ZOOM");     // draw hp in time-frequency domain;

  // build elliptical waveform
/*
  double deg2rad = TMath::Pi()/180.;
  double ePlus  = (1+cos(IOTA*deg2rad)*cos(IOTA*deg2rad))/2;
  double eCross = cos(IOTA*deg2rad);
  cout << "ePlus : " << ePlus << " eCross : " << eCross << endl;
  wavearray<double> ef=wf.hp;
  for(int i=0;i<ef.size();i++) ef[i]=ePlus*wf.hp[i]+eCross*wf.hx[i];
  MDC->DrawTime(ef,"ALP ZOOM");
*/

  // --------------------------------------------------------
  // get and plot the waveforms on L1,H1,V1
  // print the mdc parameters
  // --------------------------------------------------------

  // get mdc buffer draw waveform
  wavearray<double> x(4*16384);
  x.rate(16384);

  // draw V1
  x.start(GPS-2);
  TString log = MDC->Get(x,ifo[0]);	// get log
  x.start(0);
  MDC->DrawTime(x,"ALP ZOOM");

  // get central time
  double To = MDC->GetCentralTime(x);
  cout << "Central Time :" << To << " (sec)" << endl;

  // get central frequency
  double Fo = MDC->GetCentralFrequency(x);
  cout << "Central Frequency : " << Fo << " (Hz)" << endl;

  // draw H1
  x.start(GPS-2);
  MDC->Get(x,ifo[1]);
  x.start(0);
  MDC->DrawTime(x,"SAME",kRed);

  // draw L1
  x.start(GPS-2);
  MDC->Get(x,ifo[2]);
  x.start(0);
  MDC->DrawTime(x,"SAME",kBlue);

  // print mdc parameters
  cout << endl;
  cout << log << endl;
}
