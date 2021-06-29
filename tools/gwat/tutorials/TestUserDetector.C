//
// Test User Detector
// Author : Gabriele Vedovato

{

  detectorParams dP = {"X1", 0.0,  0.0, 0.0, 0,   90, 0,   0};
  //detectorParams dP = {"X1", 0.0,  0.0, 0.0, 0,   90, 0,   30};
  //detectorParams dP = {"X1", 0.0,  0.0, 0.0, 0,   80, 0,   20};
  //detectorParams dP = {"X1", 0.0,  0.0, 0.0, 0,   20, 0,   -40};

  detector D(dP);

  detectorParams udP = D.getDetectorParams();
  cout << "AzX : " << udP.AzX << endl;
  cout << "AzY : " << udP.AzY << endl;

  D.rotate(-60);

  udP = D.getDetectorParams();
  cout << "AzX : " << udP.AzX << endl;
  cout << "AzY : " << udP.AzY << endl;

  detector J1((char*)"J1");
  detectorParams jdP = J1.getDetectorParams();
  cout << "AzX : " << jdP.AzX << endl;
  cout << "AzY : " << jdP.AzY << endl;

  exit(0);
}
