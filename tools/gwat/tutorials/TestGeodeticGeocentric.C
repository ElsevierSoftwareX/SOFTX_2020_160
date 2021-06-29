//
// Test Geodetic to Geocentric 
// Author : Gabriele Vedovato


void printangles(double lt, double lg);


void TestGeodeticGeocentric() {

  double X,Y,Z;

  //L1
/*
  X=-7.427604192e4; Y=-5.496283721e6; Z=3.224257016e6;
  double Az = 252.2835;  // xArm L1
  double Az = 162.2835;  // yArm L1
*/
  //H1

  X=-2.161414928e6; Y=-3.834695183e6; Z=4.600350224e6;
  double Az = 324;  // xArm H1
//  double Az = 234;  // yArm H1

  //V1
/*
  X=4.5463741e6; Y=8.429897e5; Z=4.378577e6;
  double Az = 19.4326;  // xArm V1
  double Az = 289.4326;  // yArm V1
*/
  //J1
//  X=-3.7685777e6; Y=3.49218552e6; Z=3.76722999e6;
  //A1
//  X=-2.3567784e6; Y=4.8970238e6; Z=-3.3173147e6;



  int nIFO=4;
  TString ifo[4]={"H1","L1","V1","J1"};
  gnetwork gNET(nIFO,ifo);
  gNET.GetSite("H1");
  gNET.GetSite("L1");
  gNET.GetSite("V1");
  gNET.GetSite("J1");


  double latitude,longitude,elevation;

  GeocentricToGeodetic(X,Y,Z,latitude,longitude,elevation); 

  double rad2grad = 180./TMath::Pi();
  double grad2rad = TMath::Pi()/180.;
  cout << latitude*rad2grad << " " << longitude*rad2grad << " " << elevation << endl;

  printangles(latitude*rad2grad, longitude*rad2grad);

  double uN[3];
  double Alt = 0;
  //double Az = 90;
  GetCartesianComponents( uN, Alt*grad2rad, Az*grad2rad, latitude, longitude);
  cout << "uN : " << uN[0] << " " << uN[1] << " " << uN[2] << endl;
  cout << uN[0]*uN[0]+uN[1]*uN[1]+uN[2]*uN[2] << endl;

  cout << X << " " << Y << " " << Z << endl;
  GeodeticToGeocentric(latitude,longitude,elevation,X,Y,Z); 
  cout << X << " " << Y << " " << Z << endl;

  exit(0);
}

void printangles(double lt, double lg) {
  char LAT;
  double lt_t=lt;
  if(lt_t>0) LAT='N'; else {LAT='S';lt_t=-lt_t;}
  int lt_d = int(lt_t);
  int lt_m = int((lt_t-lt_d)*60);
  float lt_s = (lt_t-lt_d-lt_m/60.)*3600.;

  char LON;
  double lg_t=lg;
//  if(lg_t>180) LON='E'; else {LON='W';lg_t=180-lg_t;}
  int lg_d = int(lg_t);
  int lg_m = int((lg_t-lg_d)*60);
  float lg_s = (lg_t-lg_d-lg_m/60.)*3600.;

  cout << "LAT: " << LAT << " " << lt_d << ", " << lt_m << ", " << lt_s << " LONG: " << LON << " " << lg_d   << ", " << lg_m   << ", " << lg_s   << endl;
}

