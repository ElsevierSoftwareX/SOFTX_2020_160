// Macro is producing ring skymask for SN2007gr
// Author : Marek Szczepanczyk
//

void RingSkyMask() {

  //cWB  uses spherical coordinates
  
  // The position of HL detectors in 
  // 1) geographical coordinates:
  // H (46deg 72'18.53"N, 119deg24'27.56"W)
  // L (30deg 33'46.42"N,  90deg46'27.27"W)
  double tgH = +( 46 + 72/60 + 18.53/3600);
  double pgH = -(119 + 24/60 + 27.56/3600);
  double tgL = +( 30 + 33/60 + 46.42/3600);
  double pgL = -( 90 + 46/60 + 27.27/3600);
  // 2) cwb coordinates
  double thH =  90 - tgH;
  double phH = 360 + pgH;
  double thL =  90 - tgL;
  double phL = 360 + pgL;

  // Unit vectors
  TVector3 vX ( 1,0,0 );
  TVector3 vY ( 0,1,0 );
  TVector3 vZ ( 0,0,1 );
  
  // Unit vectors along detectors and unit vectors in cWB coordinates
  TVector3 vH ( sin(thH*PI/180) * cos(phH*PI/180), sin(thH*PI/180) * sin(phH*PI/180), cos(thH*PI/180) );
  TVector3 vL ( sin(thL*PI/180) * cos(phL*PI/180), sin(thL*PI/180) * sin(phL*PI/180), cos(thL*PI/180) );
  TVector3 vHL = vH - vL;
  //TVector3 vHL = vZ;
  // Position of HL detector on the sky:
  Double_t thHL = vHL.Theta() * 180 / PI;
  Double_t phHL = vHL.Phi()   * 180 / PI;
  cout << "The HL baseline postion on the sky is (thHL,phHL) = ( " << thHL << ", " << phHL << " )" << endl;

  // Position of SN2007gr on the sky:
  double RASN  =  40.8666;
  double DECSN = +37.3458;
  double gpss = 870772910; // GPS start time of the supernova
  double gpse = 871215278; // GPS end   time of the supernova
  double gps  = gpss;
  cout << "UTC of gps start: " << wat::Time(gps).GetDateString() << endl;
  // The location of SN cwb cooridinates
  double thSN = 90 - DECSN;
  cout << gps << endl;

  // Define a skymask
  gskymap* gSM = new gskymap(0.4,0,180,0,360);
  //#define COORDINATES "cwb"
  #define COORDINATES "geographic"
  #define RESOLUTION  1
  TString projection="cartesian";
  gSM->SetOptions(projection,COORDINATES,RESOLUTION);
  gSystem->Unlink("./plot/SN2007gr_ringskymask.gif"); // delete old file 

  double ap = 5;
  double dgps = 30 * 60;
  if (gpse - gpss > 48 * 60 * 60)
    double N = 48;
  else
    double N = floor( (gpse - gpss) / ( dgps ) );
  
  cout << "Number of frames in the gif file: " << N << endl;

  for(int j=0; j<N; ++j)
  {

    double gps = gps + dgps;
    //cout << gps << endl;
    double phSN = gSM->RA2phi(RASN,gps);
    TVector3 vSN(sin(thSN*PI/180)*cos(phSN*PI/180), sin(thSN*PI/180)*sin(phSN*PI/180), cos(thSN*PI/180));
    Double_t angSN = vHL.Angle(vSN)*180/PI;

    int L = gSM->size();	// skymap size

    int n=0;
    for (int l=0;l<L;l++) {
      // get theta,phi in CWB coordinates
      double phi = gSM->getPhi(l);
      double theta = gSM->getTheta(l);
      TVector3 vpx(sin(theta*PI/180)*cos(phi*PI/180), sin(theta*PI/180)*sin(phi*PI/180), cos(theta*PI/180));
      Double_t angpx = vHL.Angle(vpx)*180/PI;
      //double phid = fabs(phi-phHL);
      //double thetad = fabs(theta-thHL);

      //if( ( ( angpx < angSN + ap/2 ) && ( angpx > angSN - ap/2 ) ) || (thetad*thetad + phid*phid < ap*ap)  ) 
      //if( (thetad*thetad + phid*phid < ap*ap) || ((thetad-180)*(thetad-180) + phid*phid < ap*ap) ||
      //    (thetad*thetad + (phid-360)*(phid-360) < ap*ap) || ((thetad-180)*(thetad-180) + (phid-360)*(phid-360)< ap*ap) ||
      //    ( angpx < angSN + ap/2 ) && ( angpx > angSN - ap/2 ) 
      //    ) 
      if    ( ( angpx < angSN + ap/2 ) && ( angpx > angSN - ap/2 ) )
        gSM->set(l,1); 
      else
        gSM->set(l,0);
    }

    cout << "Iter j = " << j << " (thSN, phSN) = (" << thSN << ", " << phSN << ")" << " "<< "angSN = " << angSN <<endl;
    //cout << "The position of the sky is (thSN,phSN) = ( " << thSN << ", " << phSN << ")" << endl;
    //cout << "angSN = " << angSN << endl;

    //TString utc(wat::Time(gps).GetDateString());
    //TString utc(gps);
    //cout << utc << endl;
    TString title = Form("SN2007gr  gps: %9.0f ", gps );
    gSM->Draw();
    gSM->SetTitle(title);
    double pgHL = phHL;
    double tgHL = 90 - thHL;
    if (phSN<180) 
        double pgSN = phSN;
    else 
        double pgSN = phSN - 360;
    double tgSN = 90 - thSN;
    gSM->DrawMarker(pgHL,tgHL, 20, 1.0, kBlack);
    gSM->DrawMarker(pgSN,tgSN, 29, 2.0, kYellow);
    gSM->DrawText(pgHL-12, tgHL+6, "pHL", 0.05, kBlack);
    gSM->DrawText(pgSN-12, tgSN+6, "SN" , 0.05, kBlack);
    gSM->Print("./plot/SN2007gr_ringskymask.gif+40");
  }

  gSM->Print("./plot/SN2007gr_ringskymask.gif++");
}
