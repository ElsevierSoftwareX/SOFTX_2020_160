//
// Draw Network Detector Index skymap
// Author : Gabriele Vedovato

#define ODIR_NAME "plots"
#define OFILE_EXT "png"
//#define OLABEL "_OLD"
//#define OLABEL "_NEW"
//#define OLABEL "_SGW"
#define OLABEL "_SGW_hammer"
//#define OFILE_EXT "root"
//#define WRITE_PLOT

#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

//#define PROJECTION ""
#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define DISPLAY_WORLD_MAP
#define WORLD_MAP_DIR "$CWB_GWAT/data/"

#define PALETTE 0

void DrawNetworkDetectorIndexScalar(TString network="L1H1V1", double gamma = 0.0, double snr=-1.0, bool btitle = true) {

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
  if(network.Contains("I1")) ifo[nIFO++]="I1";   // INDIGO
  if(network.Contains("I2")) ifo[nIFO++]="I2";   // INDIGO 45 deg

  if(nIFO==0) {cout << "No detectors defined !!! " << endl;exit(1);}

  char ifostr[32]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,ifo[n].Data());
  }
  cout << "Network : " << ifostr << endl;

  TString title;

  gnetwork* gNET = new gnetwork(nIFO,ifo);

  gNET->setSkyMaps(0.4,0,180,0,360);
  gNET->setAntenna();
  gNET->setDelay(const_cast<char*>(ifo[0].Data()));

  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

#ifdef DISPLAY_WORLD_MAP
  gSM->SetWorldMap();
#endif

  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleFont(42);

  for(int n=0;n<(int)gNET->ifoListSize();n++) {
    detector *d = gNET->getifo(n);
    d->setPolarization(SCALAR);
  }

  cout << "gamma : " << gamma << endl;
  gNET->DrawNetworkDetectorIndex(gamma,20,snr,PALETTE,btitle);
  gNET->DrawSitesShortLabel(kBlack);
  gNET->DrawSites(kBlack,2.0);
  gNET->DrawSitesArms(1000000,kWhite,3.0);


#ifdef WRITE_PLOT
  char ofileName[128]=ODIR_NAME;
  sprintf(ofileName,"%s/",ofileName);
  for(int n=0; n<nIFO; n++) {
    sprintf(ofileName,"%s%s",ofileName,ifo[n].Data());
  }
  sprintf(ofileName,"%s%s%s.%s",ofileName,"_NDI",OLABEL,OFILE_EXT);
  cout << "Write : " << ofileName << endl;
  gNET->Print(ofileName);
  exit(0);
#endif
}

