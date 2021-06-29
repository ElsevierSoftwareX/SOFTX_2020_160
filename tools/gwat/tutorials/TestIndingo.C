//
// Test CartesianComponents with INDIGO
// Author : Gabriele Vedovato


#define LATITUDE  14.4
#define LONGITUDE 76.4

void TestIndingo() {

  double rad2grad = 180./TMath::Pi();
  double grad2rad = TMath::Pi()/180.;

  double latitude=LATITUDE;
  double longitude=LONGITUDE;
  double elevation=0.0;
  double X,Y,Z;

  cout << latitude << " " << longitude << " " << elevation << endl;
  latitude*=grad2rad;
  longitude*=grad2rad;
  elevation*=grad2rad;
  GeodeticToGeocentric(latitude,longitude,elevation,X,Y,Z); 
  cout.precision(10);
  cout << X << " " << Y << " " << Z << endl;

  double uN[3];
  double AltN = 0;
  double AzN = 0;
  GetCartesianComponents( uN, AltN*grad2rad, AzN*grad2rad, latitude, longitude);
  cout << "uN : " << uN[0] << " " << uN[1] << " " << uN[2] << endl;

  double uE[3];
  double AltE = 0;
  double AzE = 90;
  GetCartesianComponents( uE, AltE*grad2rad, AzE*grad2rad, latitude, longitude);
  cout << "uE : " << uE[0] << " " << uE[1] << " " << uE[2] << endl;

  exit(0);
}
