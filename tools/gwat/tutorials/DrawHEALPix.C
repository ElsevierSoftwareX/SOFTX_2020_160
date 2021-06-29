//
// Draw HEALPix skymap
// Author : Gabriele Vedovato

#define RESOLUTION 2

char COORDINATES[256]="cWB";
//char COORDINATES[256]="Geographic";

char PROJECTION[256]="";
//char PROJECTION[256]="hammer";
//char PROJECTION[256]="sinusoidal";

void DrawHEALPix(int order=3, bool save=false) {

  double grad2rad = TMath::Pi()/180.;
  double rad2grad = 180./TMath::Pi();

  gskymap* gSM = new gskymap(order);
  int L = gSM->size();
  cout << "Sky Patches : " << L << endl;

/*
  for(int l=0;l<L;l++) {
    double ph = gSM->getPhi(l);
    double th = gSM->getTheta(l);
    double sph = gSM->getPhiStep(l);
    double sth = gSM->getThetaStep(l);
    cout << "th : " << th << " ph : " << ph << " " <<  sth << " " << sph << endl;  
  }
*/
/*
  int ring=0; int startpix; int ringpix; double costheta; double sintheta; bool shifted=false;
  for(int i=1;i<=20;i++) {
    ring=i;
    gSM->healpix->get_ring_info(ring, startpix, ringpix, costheta, sintheta, shifted);
    cout << "ring : " << ring << " startpix : " << startpix << " ringpix : " << ringpix << endl;
  }
*/

  gSM->SetOptions(PROJECTION,COORDINATES);

  TH2D* h2 = (TH2D*)gSM->GetHistogram();

  int resolution = RESOLUTION;

  int size=2*180*resolution*2*360*resolution;

  double* phi = new double[size];
  double* theta = new double[size];
  double* binc = new double[size];

  int k=0;
  for(int i=0;i<2*180*resolution;i++) {
    for(int j=0;j<2*360*resolution;j++) {
      phi[k]=(double)j/(2*resolution);
      theta[k]=(double)i/(2*resolution);
      int l = gSM->getSkyIndex(theta[k],phi[k]);
      binc[k]=l%4+1;
      k++;
    }
  }
/*
  if (coordinate.CompareTo("GEOGRAPHIC")==0) {
    for (int i=0;i<k;i++) CwbToGeographic(phi[i],theta[i],phi[i],theta[i]);
  }
*/
  gSM->FillData(k, phi, theta, binc);
  char title[256];
  sprintf(title,"HEALPix skymap order = %d : Sky Patches = %d",order,L);
  gSM->SetTitle(title);

  delete [] phi;
  delete [] theta;
  delete [] binc;

  gSM->Draw();

  for(int l=0;l<L;l++) {
    double ph = gSM->getPhi(l);
    double th = gSM->getTheta(l);
    gSM->DrawMarker(ph,th, 20, 0.5, kBlack);  
  }

  char ofileName[256];
  sprintf(ofileName,"HEALPix_skymap_order_%d.png",order);
  if(save) gSM->Print(ofileName);
}

