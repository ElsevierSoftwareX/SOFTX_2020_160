//
// CreateCelestialSkyMask 
// Create a cwb sky mask : use the source coordinates
// Author : Gabriele Vedovato

// ---------------------------------------------
// source definitions
// ---------------------------------------------
#define SOURCE_NAME "SOURCE"

#define SOURCE_GPS -931158395

// if SOURCE_GPS>0 : SOURCE_RA and SOURCE_DEC are in Celestial coordinated
// if SOURCE_GPS=0 : SOURCE_RA and SOURCE_DEC are in Geographic coordinated (used for earth skymask)
// if SOURCE_GPS<0 : SOURCE_RA and SOURCE_DEC are in Geographic coordinated 
//                   SOURCE_GPS is used to convert them into Celestial coordinates
// Geographic : theta [-90,90] , phi [0,360]
#define SOURCE_RA  60
#define SOURCE_DEC 30

// if uncommented then GPS time is computed from UTC time
//#define SOURCE_UTC "2004-12-27 21:30:26.68 UTC Mon"

// SGR_1806_20
/*
#define SOURCE_NAME "SGR_1806_20"
#define SOURCE_GPS 788218239.68  // SGR_1806_20
#define SOURCE_RA  270.095
#define SOURCE_DEC -20.41666
#define SOURCE_UTC "2004-12-27 21:30:26.68 UTC Mon"
*/

// ---------------------------------------------
// draw definitions
// ---------------------------------------------
#define RESOLUTION  2
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define DISPLAY_WORLD_MAP
#define WORLD_MAP_DIR "$CWB_GWAT/data/"

// ---------------------------------------------
// skymask definitions
// ---------------------------------------------
#define SKYMASK_RADIUS 20  // degrees
//#define SAVE_SKYMASK
#define DRAW_SKYMASK

using namespace ROOT::Math;

void CreateCelestialSkyMask() {

  gSystem->Load("libMathCore");

  gskymap* gSM = new gskymap(0.4,0,180,0,360);
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
  int L = gSM.size();

  double gps = SOURCE_GPS;
  double ph=SOURCE_RA;
  double th=SOURCE_DEC;
  double th=90-th;  // geographic -> CWB

#if (SOURCE_GPS!=0)
#if (SOURCE_GPS<0)
  gps=-gps;
  ph = gSM->phi2RA(ph,gps);
#else
  wat::Time time(gps);
#ifdef SOURCE_UTC	// if SOURCE_UTC is defined then it is used as source gps time
  gps = time.GetGPS();
#endif
  // RA -> phi conversion
  ph = gSM->RA2phi(ph,gps);
#endif
#endif

  cout.precision(14);
  cout << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << "Source Coordinates                                   " << endl; 
  cout << "-----------------------------------------------------" << endl;
#if (SOURCE_GPS!=0)
  cout << "Time                                                 " << endl;
  cout << "UTC                   -  " << wat::Time(gps).GetDateString() << endl;
  cout << "GPS                   -  " << gps << endl;
  cout << endl;
  cout << "Celestial coordinates                                " << endl;
  cout << "DEC (deg)             -  " << SOURCE_DEC << endl;
  cout << "RA  (deg)             -  " << SOURCE_RA << endl;
  cout << endl;
#endif
  cout << "CWB coodinates                                       " << endl;
  cout << "THETA (deg)           -  " << th << endl;
  cout << "PHI (deg)             -  " << ph << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << endl;

  Polar3DVector ov1(1, th*TMath::Pi()/180, ph*TMath::Pi()/180);

  int n=0;
  for (int l=0;l<L;l++) {
    double phi = gSM->getPhi(l);
    double theta = gSM->getTheta(l);

    Polar3DVector ov2(1, theta*TMath::Pi()/180, phi*TMath::Pi()/180 );
    double Dot = ov1.Dot(ov2);
    double dOmega = 180.*TMath::ACos(Dot)/TMath::Pi();
    //cout << "dOmega : " << dOmega << endl;

    if(dOmega<SKYMASK_RADIUS) gSM->set(l,1); else gSM->set(l,0);
  }

#ifdef SAVE_SKYMASK
  char oFile[256];
#if (SOURCE_GPS<0)
  sprintf(oFile,"CelestialSkyMask_DEC_%3.1f_RA_%3.1f_GPS_N%3.1f_RADIUS_%3.1f",
          SOURCE_DEC,ph,gps,SKYMASK_RADIUS);
#else
  sprintf(oFile,"CelestialSkyMask_DEC_%3.1f_RA_%3.1f_GPS_%3.1f_RADIUS_%3.1f",
          SOURCE_DEC,SOURCE_RA,gps,SKYMASK_RADIUS);
#endif
  TString oFileName = oFile;
  oFileName.ReplaceAll(".","d");
  oFileName=oFileName+".txt";
  cout << "Save File : " << oFileName.Data() << endl;
  ofstream out;
  out.open(oFileName.Data(), ios::out);
  if (!out.good()) {cout << "Error Opening File : " << oFileName.Data() << endl;exit(1);}
  for (int l=0;l<L;l++) out << l << " " << gSM->get(l) << endl;
  out.close();
  exit(0);
#endif

// Draw SkyMask in geographical coordinates
#ifdef DRAW_SKYMASK

#ifdef DISPLAY_WORLD_MAP
  TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);
  gSM->SetWorldMapPath(world_map.Data());
  gSM->SetWorldMap();
#endif

  double th1,ph1;
  CwbToGeographic(ph,th,ph1,th1);

  char title[256];
  sprintf(title,"Sky Mask : radius = %d^{0} - center = (DEC=%3.3f,RA=%3.3f)",SKYMASK_RADIUS,th1,ph1);
  gSM->SetTitle(title);

  gSM->Draw(0);

  gSM->DrawMarker(ph1,th1, 29, 2.0, kYellow);
  gSM->DrawText(ph1-14, th1+8, SOURCE_NAME, 0.04, kBlack);
  //gSM->DrawCircles(ph1,th1,(Color_t)kWhite);
#endif

}
