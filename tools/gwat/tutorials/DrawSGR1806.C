//
// Draw SGR1806 sky position & bar detectors antenna pattern
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
#define E1_ENABLED

#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define DISPLAY_WORLD_MAP

//#define WRITE_PLOT

// ARG 1806-20 SOURCE DEFINITION (RA 18h 05.7m  Dec -20d 25m) Distance : 14.5 ± 1.4 kpc
// http://www.journals.uchicago.edu/ApJ/journal/issues/ApJ/v478n2/33968/
 
#define DRAW_SGR1806
#define SGR1806_NAME "SGR_1806_20"
#define SGR1806_DATE "2004-12-27 21:30:26.68 UTC Mon" 
#define SGR1806_GPS_TIME 788218239.68
//#define SGR1806_RA  270.095
//#define SGR1806_DEC (90-(-20.41666))
#define SGR1806_RA  (270.095-360.0)
#define SGR1806_DEC -20.41666
 
#define DRAW_SGRA
#define SGRA_NAME "SgrA*"
#define SGRA_DATE "2004-12-27 21:30:26 UTC Mon" 
#define SGRA_GPS_TIME 788218239.68
// cWB
//#define SGRA_RA  266.41667
//#define SGRA_DEC (90-(-29.11667))
// Geographic
#define SGRA_RA  (266.41667-360.0)
#define SGRA_DEC -29.11667

//#define DRAW_DELAY

void DrawSGR1806(bool draw_sensitivity = true) {

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

  if(draw_sensitivity) sprintf(title,"%s     %s      DPF F+ - Network : %s",SGR1806_NAME,SGR1806_DATE,ifostr);
  else sprintf(title,"%s     %s      DPF Fx - Network : %s",SGR1806_NAME,SGR1806_DATE,ifostr);

  gnetwork* gNET = new gnetwork(nIFO,ifo);
  gskymap* gSM = gNET->GetGskymap();

  gNET->setSkyMaps(0.4,0,180,0,360);
  gNET->setAntenna();
  gNET->setDelay(const_cast<char*>(ifo[0].Data()));

  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

#ifdef DISPLAY_WORLD_MAP
  gSM->SetWorldMap();
#endif

  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  //h2->GetZaxis()->SetRangeUser(0,1.4);
  gSM->SetTitle(title);

  double phi,theta;
#ifdef DRAW_SGR1806
  phi=SGR1806_RA;
  theta=SGR1806_DEC;
#endif

#ifdef DRAW_DELAY
  //gNET->DrawDelay(ifo[0],ifo[1]);
  //gNET->DrawDelay(ifo[0],ifo[1],phi,theta);
  gNET->DrawDelay(ifo[0],ifo[1],phi,theta,930.);
#else
  gNET->DrawAntennaPattern(draw_sensitivity);
#endif
  gNET->DrawSites(kWhite,0.6);
  gNET->DrawSitesLabel(kWhite,0.025);

  gNET->GetSite(ifo[0]);

#ifdef DRAW_SGR1806
  gSM->DrawMarker(phi,theta, SGR1806_GPS_TIME, 29, 2.0, kYellow);  
  gSM->DrawText(phi-14, theta+8, SGR1806_GPS_TIME, SGR1806_NAME, 0.04, kWhite);
  cout << "Time Delay -> " << gNET->GetDelay(ifo[0],ifo[1],phi,theta) << endl;
  gNET->DrawCircles(phi,theta,SGR1806_GPS_TIME,kWhite);
#endif

#ifdef DRAW_SGRA
  phi=SGRA_RA;
  theta=SGRA_DEC;
  //gSM->DrawMarker(phi,theta, 29, 2.0, kBlack);  
  //gSM->DrawText(phi-14, theta-14, SGRA_NAME, 0.04, kBlack);
  gSM->DrawMarker(phi,theta, SGR1806_GPS_TIME, 29, 2.0, kBlack);  
  gSM->DrawText(phi-14, theta-14, SGR1806_GPS_TIME, SGRA_NAME, 0.04, kBlack);
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

