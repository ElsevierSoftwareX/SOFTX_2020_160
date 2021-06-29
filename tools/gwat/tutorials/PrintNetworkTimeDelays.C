//
// Print Network Time Delays
// Author : Gabriele Vedovato


void PrintNetworkTimeDelays(TString network, double phi, double theta, bool burstmdc=true) {

  int nIFO=0;
  TString ifo[10];
  if(network.Contains("H1")) ifo[nIFO++]="H1";   // LHO1
  if(network.Contains("L1")) ifo[nIFO++]="L1";   // LLO
  if(network.Contains("G1")) ifo[nIFO++]="G1";   // GEO
  if(network.Contains("V1")) ifo[nIFO++]="V1";   // VIRGO
  if(network.Contains("T1")) ifo[nIFO++]="T1";   // TAMA
  if(network.Contains("H2")) ifo[nIFO++]="H2";   // LHO2
  if(network.Contains("A1")) ifo[nIFO++]="A1";   // AIGO
  if(network.Contains("O1")) ifo[nIFO++]="O1";   // AURIGA
  if(network.Contains("N1")) ifo[nIFO++]="N1";   // NAUTILUS
  if(network.Contains("E1")) ifo[nIFO++]="E1";   // EXPLORER
  if(network.Contains("A2")) ifo[nIFO++]="A2";   // AUSTRALIAN 90°
  if(network.Contains("J1")) ifo[nIFO++]="J1";   // JAPANESE
  if(network.Contains("I1")) ifo[nIFO++]="I1";   // INDINGO
  if(network.Contains("I2")) ifo[nIFO++]="I2";   // INDINGO 45 deg

  if(nIFO==0) {cout << "No detectors defined !!! " << endl;exit(1);}

  char ifostr[32]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,ifo[n].Data());
  }
  cout << "Network : " << ifostr << endl;

  TString title;

  gnetwork* gNET = new gnetwork(nIFO,ifo);

/* ------------------------------
  earthCenterTime = _EarthCtrGPS;
  phi = _External_phi;
  theta = _External_x;
  psi = _External_psi;
--------------------------------- */
  double Pi = TMath::Pi();
  if(burstmdc) { 
    theta = acos(theta);
    theta*= 180/Pi;
    phi = phi > 0 ? phi : 2*Pi+phi;
    phi*= 180/Pi;
//    phi = sm.phi2RA(phi,_EarthCtrGPS);
//    psi*= 180/Pi;
  }

  for(int i=0;i<nIFO;i++) {
    for(int j=i+1;j<nIFO;j++) {
      cout << ifo[i].Data() << " " << ifo[j].Data() << " -> " 
           << gNET->GetDelay(ifo[i].Data(),ifo[j].Data(),phi,theta) << " sec " << endl;
    }
  }
  exit(0);
}

