//
// Draw Bar Detectors Antenna Patterns
// Author : Gabriele Vedovato

//#define L1_ENABLED
//#define H1_ENABLED
//#define V1_ENABLED
//#define H2_ENABLED
//#define G1_ENABLED
//#define T1_ENABLED
//#define A1_ENABLED
#define O1_ENABLED
#define N1_ENABLED
//#define E1_ENABLED

#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"
//#define PROJECTION "parabolic"

#define DISPLAY_WORLD_MAP

//#define WRITE_PLOT

#define DRAW_SGRA
#define SGRA_NAME         "SgrA*"
#define SGRA_DEC          -29.11667
#define SGRA_RA           266.41667

#define DRAW_SGR1806
#define SGR1806_NAME "SGR_1806_20"
#define SGR1806_DATE "2004-12-27 21:30:26.68 UTC Mon"
#define SGR1806_GPS_TIME 788218239.68
#define SGR1806_RA  270.095
#define SGR1806_DEC (90-(-20.41666))

void DrawO1N1E1(bool draw_sensitivity = true) {

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
#ifdef O1_ENABLED
   ifo[nIFO++]="O1";   // AURIGA
#endif
#ifdef N1_ENABLED
   ifo[nIFO++]="N1";   // NAUTILUS
#endif
#ifdef E1_ENABLED
   ifo[nIFO++]="E1";   // EXPLORER
#endif

  cout << ifo[0] << " " << ifo[1] << endl;

  char ifostr[32]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,ifo[n].Data());
  }

  char title[256];
  TString ofileName;

  //sprintf(title,"AURIGA NAUTILUS EXPLORER - Network Sensitivity");
  sprintf(title,"AURIGA NAUTILUS - Network Sensitivity (%s)",SGR1806_DATE);

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

  gSM->SetGalacticDiskColor(kBlack);
  gSM->SetGalacticDisk(788218239.68);

  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  //h2->GetZaxis()->SetRangeUser(0,1.4);
  gNET->SetTitle(title);

  gNET->DrawAntennaPattern(draw_sensitivity,57);
  gNET->DrawSites(kBlack,1.0);
  //gNET->DrawSitesLabel(kBlack,0.025);

  for(int n=0;n<nIFO;n++) gNET->GetSite(ifo[n]);
  double theta[3];
  double phi[3];
  for(int n=0;n<nIFO;n++) theta[n]=gNET->GetSite(ifo[n],"theta");
  for(int n=0;n<nIFO;n++) phi[n]=gNET->GetSite(ifo[n],"theta");

  double Tphi[3]={22,-4,-25};
  double Ttheta[3]={44.25,30.25,53.58};

  gSM->DrawText(Tphi[0], Ttheta[0], "AURIGA", 0.04, kBlack);
  gSM->DrawText(Tphi[1], Ttheta[1], "NAUTILUS", 0.04, kBlack);
  //gSM->DrawText(Tphi[2], Ttheta[2], "EXPLORER", 0.04, kBlack);

#ifdef DRAW_SGRA
  if(COORDINATES=="Geographic") {
    gSM->DrawMarker(SGRA_RA-420, SGRA_DEC, 29, 2.0, kBlack);
    gSM->DrawText(SGRA_RA-8-420, SGRA_DEC-8, SGRA_NAME, 0.04, kBlack);
  }
#endif

#ifdef DRAW_SGR1806
  double sphi,stheta;
  sphi=SGR1806_RA;
  stheta=SGR1806_DEC;
  skymap sm;
  sphi = sm.RA2phi(sphi,SGR1806_GPS_TIME);
  if(COORDINATES=="Geographic") {
    CwbToGeographic(sphi,stheta,sphi,stheta);
  } else {
    gSM->DrawMarker(sphi,stheta, 29, 2.0, kBlack);
    gSM->DrawText(sphi-14, stheta+8, SGR1806_NAME, 0.04, kBlack);
  }
#endif

  double mTau=gNET->GetDelay(ifo[0],ifo[1],"MAX");  // maximum time delay
  double dTau=gNET->GetDelay(ifo[0],ifo[1]);       // time delay difference
  cout<<"maximum time delay between detectors: "<<mTau<<endl;
  cout<<"       maximum time delay difference: "<<dTau<<endl;

#ifdef WRITE_PLOT
  cout << "Write : " << ofileName << endl;
  gSM->Print(ofileName);
#endif
}

