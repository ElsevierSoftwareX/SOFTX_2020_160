CWB::mdc* MDC;

void DrawSGC() {
  //
  // Show how to use mdc class to get & draw SGC waveforms
  // Author : Gabriele Vedovato

  #include <vector>
  #define N_IFO 3

  TString ifo[N_IFO]={"L1","H1","V1"};
  MDC = new CWB::mdc(N_IFO,ifo);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par;

  par.resize(2);

  par[0].name="frequency"; par[0].value=100.;
  par[1].name="Q"; par[1].value=8.9;
  MDC->AddWaveform(MDC_SGC, par);

  par[0].name="frequency"; par[0].value=153.;
  par[1].name="Q"; par[1].value=8.9;
  MDC->AddWaveform(MDC_SGC, par);

  par[0].name="frequency"; par[0].value=1053.;
  par[1].name="Q"; par[1].value=9;
  MDC->AddWaveform(MDC_SGC, par);

  par[0].name="frequency"; par[0].value=1304.;
  par[1].name="Q"; par[1].value=9;
  MDC->AddWaveform(MDC_SGC, par);

  MDC->Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjHrss(2.5e-21);
  MDC->SetInjRate(0.02);
  MDC->SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  // 934361014 ra : 3.4909999370575 dec : -27.687000274658 radius : 1.3999999761581
  float DEC = -27.687000274658;
  float RA = 3.4909999370575;
  float seed = 0;
  par.resize(3);
  par[0].name="theta";   par[0].value=DEC;
  par[1].name="phi";     par[1].value=RA;
  par[2].name="entries"; par[2].value=1000;
  MDC->SetSkyDistribution(MDC_CELESTIAL_FIX,par,seed);

  // Get the first waveform hp,hx components 
  wf = MDC->GetWaveform("SGC1304Q9",0);

  cout << "size : " << wf.hp.size() << " rate : " << wf.hp.rate() << " start : " << (int)wf.hp.start() << endl;
  wf.hp.start(0);		// set start to 0 (needed by draw Method)
  wf.hx.start(0);

  //MDC->Draw(wf.hp,MDC_TIME,"ALP ZOOM");
  //MDC->Draw(wf.hx,MDC_TIME,"same",kRed);

  //MDC->Draw(wf.hp,MDC_FFT,"ALP ZOOM");	// draw hp in frequency domain
  //MDC->Draw(wf.hp,MDC_TF,"ZOOM");     // draw hp in time-frequency domain;

}
