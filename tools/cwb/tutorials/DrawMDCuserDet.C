CWB::mdc* MDC;
 
void DrawMDCuserDet() {
  //
  // Draw Builtin Waveform with the mdc class (for user defined detectors)
  // Author : Gabriele Vedovato

  #include <vector>

  // -------------------------------------------------------- 
  // define network
  // -------------------------------------------------------- 

  int nIFO = 4;
  char ifo[NIFO_MAX][4];
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"");    // if empty then select user detectors
  strcpy(ifo[2],"");    // if empty then select user detectors
  strcpy(ifo[3],"V1");

  detectorParams detParms[4] = {
                          {"L1", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, (    -197.716)},

                          {"Y2", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, ( +45-125.999)},  // H2 LCI
                          {"Y3", 46.4551, -119.408, 0.0, 0, ( +45-125.999), 0, (    -125.999)},  // H3 LCI

                          {"V1", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},
                         };

  network NET;
  detector* pD[NIFO_MAX];
  for(int i=0; i<nIFO; i++) {
    if(strlen(ifo[i])>0) pD[i] = new detector(ifo[i]);        // built in detector
    else                 pD[i] = new detector(detParms[i]);   // user define detector
  }
  for(int i=0; i<nIFO; i++) NET.add(pD[i]);

  MDC = new CWB::mdc(&NET); 

/*
  network* net = MDC->GetNetwork();
  int nIFO=net->ifoListSize();
  for(int n=0;n<nIFO;n++) cout << n << " " << net->ifoName[n] << endl;
  detector *ppD[NIFO_MAX];
  for(int n=0; n<nIFO; n++) ppD[n] = net->getifo(n);
  double theta=100.;
  double phi=50.;
  double psi=10.;
  wavecomplex F;

  for(int n=0; n<nIFO; n++) {
    F = pD[n]->antenna(theta,phi,psi);
    cout << n << " -> F+ : " << F.real() << " Fx : " << F.imag() << " " << pD[n]->getTau(theta,phi) << endl;
    F = ppD[n]->antenna(theta,phi,psi);
    cout << n << " -> F+ : " << F.real() << " Fx : " << F.imag() << " " << ppD[n]->getTau(theta,phi) << endl;
    cout << endl;
  }
  exit(0);
*/

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par;
/*
  par.resize(3);
  par[0].name="frequency"; par[0].value=250.;
  par[1].name="bandwidth"; par[1].value=100.;
  par[2].name="duration";  par[2].value=0.1;
  MDC->AddWaveform(MDC_WNB, par);
*/

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=3.;
  MDC->AddWaveform(MDC_SG, par);

/*
  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=8.9;
  MDC->AddWaveform(MDC_SG, par);
*/
/*
  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=9.;
  MDC->AddWaveform(MDC_SGC, par);
*/
  MDC->Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjHrss(2.5e-21);
  MDC->SetInjRate(0.0333333);
  MDC->SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  par.resize(4);
  par[0].name="theta"; par[0].value=30;
  par[1].name="phi";   par[1].value=60;
  par[2].name="psi";   par[2].value=0;
  par[3].name="gps";   par[3].value=968650050;
  MDC->SetSkyDistribution(MDC_EARTH_FIX,par);


  int id=1;
  MDC->Draw("Y2",968650000,968650100,id,MDC_TIME);
  MDC->Draw("Y3",968650000,968650100,id,MDC_TIME,"SAME",2);
//  MDC->Draw("L1",968650000,968650100,id,MDC_TIME,"SAME",3);
//  MDC->Draw("V1",968650000,968650100,id,MDC_TIME,"SAME",5);

  //MDC->Draw("Y3",968650000,968650100,id,MDC_FFT);
  //MDC->Draw("Y3",968650000,968650100,id,MDC_TF);
}
