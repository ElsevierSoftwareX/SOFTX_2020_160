//
// Test Site J1
// Author : Gabriele Vedovato


#define L1_ENABLED
#define H1_ENABLED
#define V1_ENABLED
//#define H2_ENABLED
//#define G1_ENABLED
//#define T1_ENABLED
//#define A1_ENABLED
//#define A2_ENABLED
//#define O1_ENABLED
//#define N1_ENABLED
//#define E1_ENABLED
#define J1_ENABLED

#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define DISPLAY_WORLD_MAP

//#define WRITE_PLOT

#define INJ_PHI    9.064681e-01
#define INJ_THETA -2.885531e-01
#define INJ_PSI    1.427236

//#define TITLE "LHO1 LLO VIRGO LCGT Network Sensitivity |F+|"
//#define TITLE "LHO1 LLO VIRGO LCGT Network Alignment |Fx|/|F+|"

//#define TITLE "LHO1 LLO VIRGO LCGT AIGO2 Network Sensitivity |F+|"
//#define TITLE "LHO1 LLO VIRGO LCGT AIGO2 Network Alignment |Fx|/|F+|"

//#define TITLE "LHO1 LLO VIRGO LCGT AIGO2 Network Sensitivity |F+|"
#define TITLE "LHO1 LLO VIRGO LCGT AIGO1 Network Alignment |Fx|/|F+|"

void TestSiteJ1(bool draw_sensitivity = true) {

  int nIFO=0;
  TString ifo[10];
#ifdef H1_ENABLED
   ifo[nIFO++]="H1";   // LHO1
#endif
#ifdef L1_ENABLED
   ifo[nIFO++]="L1";   // LLO
#endif
#ifdef G1_ENABLED
   ifo[nIFO++]="G1";  // GEO
#endif
#ifdef V1_ENABLED
   ifo[nIFO++]="V1";   // VIRGO
#endif
#ifdef T1_ENABLED
   ifo[nIFO++]="T1";  // TAMA
#endif
#ifdef H2_ENABLED
   ifo[nIFO++]="H2";   // LHO2
#endif
#ifdef A1_ENABLED
   ifo[nIFO++]="A1";   // AIGO
#endif
#ifdef A2_ENABLED
   ifo[nIFO++]="A2";   // AIGO 45°
#endif
#ifdef O1_ENABLED
   ifo[nIFO++]="O1";   // AURIGA
#endif
#ifdef N1_ENABLED
   ifo[nIFO++]="N1";   // NAUTILUS
#endif
#ifdef E1_ENABLED
   ifo[nIFO++]="E1";   // EXPLORER
#endif
#ifdef J1_ENABLED
   ifo[nIFO++]="J1";   // LCGT
#endif

  char ifostr[32]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,ifo[n].Data());
  }

  TString title;
  TString ofileName;

  gnetwork* gNET = new gnetwork(nIFO,ifo);
  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

  gNET->setSkyMaps(0.4,0,180,0,360);
  gNET->setAntenna();
  gNET->setDelay(const_cast<char*>(ifo[0].Data()));

#ifdef DISPLAY_WORLD_MAP
  gSM->SetWorldMap();
#endif

  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  gNET->SetTitle(ifostr);
#ifdef TITLE
  gNET->SetTitle(TITLE);
#endif

  gNET->DrawAntennaPattern(draw_sensitivity,true);  // true alignment
  gNET->DrawSites(kWhite,0.6);
  gNET->DrawSitesLabel(kWhite,0.025);

  gNET->GetSite("H1");
  gNET->GetSite("L1");
  gNET->GetSite("V1");
  gNET->GetSite("J1");

  double Pi=TMath::Pi();

  double theta = 0;
  theta = acos(INJ_THETA);
  theta*= 180/Pi;

  double phi = 0;
  phi = INJ_PHI > 0 ? INJ_PHI : 2*Pi+INJ_PHI;
  phi*= 180/Pi;

  double psi = 0;
  psi = INJ_PSI;
  psi*= 180/Pi;

  cout << "phi : " << phi << " theta : " << theta << " psi : " << psi << endl;

  double Fp = gNET->GetAntennaPattern(phi,theta,psi,true);
  double Fc = gNET->GetAntennaPattern(phi,theta,psi,false);

  cout << "DPF -> " << " Fp " << Fp << " Fc : " << Fc << endl;

  detector J1((char*)"J1");
  detector L1((char*)"L1");
  cout << endl;
  wavecomplex fJ1,fL1;
  fJ1 = J1.antenna(theta,phi,psi);
  fL1 = L1.antenna(theta,phi,psi);
  cout << endl;
  cout << " fJ1p : " << fJ1.real() << " fJ1x : " << fJ1.imag() << endl;
  cout << " fL1p : " << fL1.real() << " fL1x : " << fL1.imag() << endl;

#ifdef WRITE_PLOT
  cout << "Write : " << ofileName << endl;
  gSM->Print(ofileName);
#endif
}

