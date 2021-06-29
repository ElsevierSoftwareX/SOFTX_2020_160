//
// Convert Delay to Sky Coordinates
// Author : Gabriele Vedovato

#define ODIR_NAME "plots"

void Delays2Directions(TString network="L1H1V1") {

  int nIFO=0;
  TString ifo[10];
  if(network.Contains("V1")) ifo[nIFO++]="V1";   // VIRGO
  if(network.Contains("H1")) ifo[nIFO++]="H1";   // LHO1
  if(network.Contains("L1")) ifo[nIFO++]="L1";   // LLO
  if(network.Contains("G1")) ifo[nIFO++]="G1";   // GEO
  if(network.Contains("T1")) ifo[nIFO++]="T1";   // TAMA
  if(network.Contains("H2")) ifo[nIFO++]="H2";   // LHO2
  if(network.Contains("A1")) ifo[nIFO++]="A1";   // AIGO
  if(network.Contains("O1")) ifo[nIFO++]="O1";   // AURIGA
  if(network.Contains("N1")) ifo[nIFO++]="N1";   // NAUTILUS
  if(network.Contains("E1")) ifo[nIFO++]="E1";   // EXPLORER
  if(network.Contains("A2")) ifo[nIFO++]="A2";   // AUSTRALIAN 90°
  if(network.Contains("J1")) ifo[nIFO++]="J1";   // JAPANESE

  if(nIFO==0) {cout << "No detectors defined !!! " << endl;exit(1);}

  char ifostr[32]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,ifo[n].Data());
  }
  cout << "Network : " << ifostr << endl;

  gnetwork gNET(nIFO,ifo);

  TString ifoName[3];
  for(int n=0; n<nIFO; n++) ifoName[n]=TString(gNET.getifo(n)->Name);
  for(int n=0; n<nIFO; n++) cout << n << " " << ifoName[n].Data() << endl;

/*
H-V : 874649008.2204 - 874649008.2001 =  0.0203
L-V : 874649008.2172 - 874649008.2008 =  0.0165
L-H : 874649008.2257 - 874649008.2286 = -0.0029
*/

  double ph1,ph2,th1,th2;

  gNET.Delay2Coordinates(0,0.0203,0.0165,ph1,th1,ph2,th2);
  cout << "ph1 : " << ph1 << " th1 : " << th1 << " ph2 : " << ph2 << " th2 : " << th2 << endl;
  CwbToGeographic(ph1,th1,ph1,th1);
  CwbToGeographic(ph2,th2,ph2,th2);
  cout << "ph1 : " << ph1 << " th1 : " << th1 << " ph2 : " << ph2 << " th2 : " << th2 << endl;

  exit(0);
}

