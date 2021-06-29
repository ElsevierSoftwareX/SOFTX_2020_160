//
// Read and Draw Gravitational Wave Galaxy Catalog
// Author : Gabriele Vedovato

//#define RESOLUTION  2 
#define RESOLUTION  4

#define GWGCCatalog "$CWB_GWAT/data/GWGCCatalog_Rev1d8.txt"
#define DISTANCE_THR 50 // Mpc

// celestial coordinates
#define SGRA_NAME         "SgrA*"
#define SGRA_DEC          -29.11667
#define SGRA_RA           266.41667

#define NGC0224_NAME      "M31" 
#define NGC0224_DEC       41.2687 
#define NGC0224_RA        10.6846 

#define NGC0292_NAME      "SMC"
#define NGC0292_DEC       -72.8002
#define NGC0292_RA        13.1583

#define ESO056_115_NAME   "LMC"
#define ESO056_115_DEC    -69.7561
#define ESO056_115_RA     80.8941

//#define COORDINATES "Geographic"
#define COORDINATES ""
#define PROJECTION "hammer"

//#define WRITE_PLOT

#define NMAX 2136800  

//#define SKYMASK_FILE_NAME  "SkyMaskCC_GWGC50MPC_Rev1d7_R0d4.txt"
//#define SKYMASK_FILE_NAME  "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d4.txt"

//#define SKYMASK_FILE_NAME  "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d4_PROB_INV_DIST.txt"
#define SKYMASK_FILE_NAME  "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d4_NSOURCES.txt"
#define SKYMASK_RESOLUTION 0.4

//#define SKYMASK_FILE_NAME  "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d4_PHPAD.txt"
//#define SKYMASK_RESOLUTION 0.4

//#define SKYMASK_FILE_NAME  "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d2_PHPAD.txt"
//#define SKYMASK_RESOLUTION 0.2

//#define SKYMASK_FILE_NAME  "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d1_PHPAD.txt"
//#define SKYMASK_RESOLUTION 0.1

#define REJECTED_SKY_PIXEL_PERCENTAGE 1.0
//#define ADD_PHI_PADDING
#define OUTPUT_PROBABILITY_SKYMAP

int ReadGWGCCatalogToSkymap(skymap& sm);
int ReadGWGCCatalog(double*& x, double*& y, double*& z, int& nGWGC);
int ReadBRST50MPC_S6_ToSkymask(skymap& sm);

void CreateProbSkyMaskGWGC() {

  double *x,*y,*z;

  int nGWGC=0;
  int size = ReadGWGCCatalog(x,y,z,nGWGC);
  //int size = ReadBRST50MPC_S6(x,y,z,nGWGC);

#ifdef SKYMASK_FILE_NAME

  skymap sm(SKYMASK_RESOLUTION,0,180,0,360);
  int L = sm.size();
  //nGWGC=ReadGWGCCatalogToSkymap(sm);
  nGWGC=ReadBRST50MPC_S6_ToSkymask(sm);

#ifdef ADD_PHI_PADDING
  skymap sm2(sm);
  for (int l=0;l<L;l++) {
    if(sm.get(l)) {
      double th = sm.getTheta(l);
      double ph = sm.getPhi(l);
      double sp = sm.getPhiStep(l);

      int ll = sm.getSkyIndex(th,ph-sp);    // neighbour sky index
      int lr = sm.getSkyIndex(th,ph+sp);    // neighbour sky index

      sm2.set(ll,sm.get(ll)+1);
      sm2.set(lr,sm.get(lr)+1);
    } 
  } 
  for (int l=0;l<L;l++) sm.set(l,sm2.get(l));
#endif

  ofstream out;
  out.open(SKYMASK_FILE_NAME, ios::out);
  if (!out.good()) {cout << "Error Opening File : " << SKYMASK_FILE_NAME << endl;exit(1);}
  cout << "Write File : " << SKYMASK_FILE_NAME << endl;
#ifdef OUTPUT_PROBABILITY_SKYMAP
  double norm=0;
  for (int l=0;l<L;l++) norm+=sm.get(l);
  //for (int l=0;l<L;l++) out << l << " " << sm.get(l)/norm << endl;
  for (int l=0;l<L;l++) out << l << " " << sm.get(l) << endl;
  out.close();
#else
  double perc=0;
  for (int l=0;l<L;l++) {
    int mask = sm.get(l)>0 ? 1 : 0;
    out << l << " " << mask << endl;
    perc+=mask;
  }
  out.close();
  perc/=L;
  cout << "SkyMaskCC perc : " << perc*100 << endl; 
#endif

#endif

  char title[256];
  //sprintf(title,"Statistics of sources in the GWGC catalog out to a %d Mpc measured distance (# %d)",DISTANCE_THR,nGWGC);
  sprintf(title,"Statistics of sources in the BRST50MPC_S6 MDC SET out to a %d Mpc measured distance (# %d)",DISTANCE_THR,nGWGC);

  gskymap* gSM = new gskymap(SKYMASK_RESOLUTION,0,180,0,360);
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
  gSM->SetTitle(title);
  gSM->SetGridxColor(kWhite);
  gSM->SetGridyColor(kWhite);
  gSM->SetGalacticDisk(0);
  gSM->SetGalacticDiskColor(kYellow);
  gSM->FillData(size, x, y, z);
  //gSM->FillData(sm);
  gSM->SetZaxisTitle("Mpc");
  gSM->Draw(-2);

  gSM->DrawMarker(SGRA_RA-180, SGRA_DEC, 4, 2.0, kYellow);
  gSM->DrawText(SGRA_RA-5-180, SGRA_DEC+5, SGRA_NAME, 0.04, kYellow);
/*
  gSM->DrawMarker(NGC0292_RA-180, NGC0292_DEC, 4, 2.0, kYellow);
  gSM->DrawText(NGC0292_RA-10-180, NGC0292_DEC+10, NGC0292_NAME, 0.04, kYellow);

  gSM->DrawMarker(ESO056_115_RA-180, ESO056_115_DEC, 4, 2.0, kYellow);
  gSM->DrawText(ESO056_115_RA-0-180, ESO056_115_DEC+10, ESO056_115_NAME, 0.04, kYellow);

  gSM->DrawMarker(NGC0224_RA-180, NGC0224_DEC, 4, 2.0, kYellow);
  gSM->DrawText(NGC0224_RA+5-180, NGC0224_DEC-10, NGC0224_NAME, 0.04, kYellow);
*/

#ifdef WRITE_PLOT
  TString ofileLabel = gSystem->ExpandPathName(GWGCCatalog);
  ofileLabel.ReplaceAll(".txt","");
  ofileLabel.ReplaceAll(".","");
  ofileLabel.ReplaceAll("/","");
  char ofileName[256];
  sprintf(ofileName,"%s_Distance%dMpc%s.png",ofileLabel.Data(),DISTANCE_THR,PROJECTION);
  cout << "ofileName : " << ofileName << endl;

  cout << "Write : " << ofileName << endl;
  skyplot->Print(ofileName);
#endif
}

int
ReadGWGCCatalogToSkymap(skymap& sm) {

  ifstream in;
  in.open(gSystem->ExpandPathName(GWGCCatalog),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << gSystem->ExpandPathName(GWGCCatalog) << endl;exit(1);}

  char iline[1024];
  in.getline(iline,1024);  // ski first line (header)
  nGWGC=0;
  while (1) {

    in.getline(iline,1024);
    if (!in.good()) break;
    //cout << iline << endl;
    TObjArray* tok = TString(iline).Tokenize(TString('|'));

    TObjString* tname = (TObjString*)tok->At(1);
    TObjString* tra   = (TObjString*)tok->At(2);
    TObjString* tdec  = (TObjString*)tok->At(3);
    TObjString* tdist = (TObjString*)tok->At(14);

    TString name = tname->GetString();
    double ra    = tra->GetString().Atof();
    double dec   = tdec->GetString().Atof();
    double dist  = tdist->GetString().Atof();
    if (dist<DISTANCE_THR) {

      double th = -(dec-90);
      double ph = ra*360./24.;
      int ind = sm.getSkyIndex(th,ph);
      sm.set(ind,sm.get(ind)+1);

      nGWGC++;
    }
  }
  in.close();

  return nGWGC;
}


int
ReadGWGCCatalog(double*& x, double*& y, double*& z, int& nGWGC) {


  x = new double[NMAX];
  y = new double[NMAX];
  z = new double[NMAX];

  ifstream in;
  in.open(gSystem->ExpandPathName(GWGCCatalog),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << gSystem->ExpandPathName(GWGCCatalog) << endl;exit(1);}

  char iline[1024];
  in.getline(iline,1024);  // ski first line (header)
  nGWGC=0;
  while (1) {

    in.getline(iline,1024);
    if (!in.good()) break;
    //cout << iline << endl;
    TObjArray* tok = TString(iline).Tokenize(TString('|'));

    TObjString* tname = (TObjString*)tok->At(1);
    TObjString* tra   = (TObjString*)tok->At(2);
    TObjString* tdec  = (TObjString*)tok->At(3);
    TObjString* tdist = (TObjString*)tok->At(14);

    TString name = tname->GetString();
    double ra    = tra->GetString().Atof();
    double dec   = tdec->GetString().Atof();
    double dist  = tdist->GetString().Atof();
    if (dist<DISTANCE_THR) {

      //Geographic coordinates
      double ph = ra*360./24.+180.;
      double th = dec;
//if(ph<0) cout << "ph  " << ph << endl;
//if(th<0) cout << "th  " << th << endl;
      //double ph = ra*360./24.;
      //double th = -(dec-90);
      //double th = -(-dec-90);

      x[nGWGC]=ph;
      y[nGWGC]=th;
      z[nGWGC]=dist;
      nGWGC++;
if (name.CompareTo("NGC0224")==0) {  // M31
//if (name.CompareTo("NGC0292")==0) {  // SMC
//if (name.CompareTo("ESO056-115")==0) {  // LMC
//if (fabs(dist-0.05)<0.01) {
  cout << name.Data() << " ra : " << ra << " dec : " << dec << " " << dist << endl;
  cout << "PH : " << ph << " TH : " << th << " " << dist << endl;
}
    }
  }
  in.close();
  cout << "nGWGC : " << nGWGC << endl;

  int size=nGWGC;
  for (int i=0;i<360*RESOLUTION;i++) {
    for (int j=0;j<180*RESOLUTION;j++) {
      double ph = i/(double)RESOLUTION;
      double th = j/(double)RESOLUTION;
      x[size]=ph;
      y[size]=th;
      z[size]=DISTANCE_THR;
      size++;
    }
  }

  Int_t *index = new Int_t[size];
  TMath::Sort(size,z,index,true);
  double T[NMAX];
  for (int i=0;i<size;i++) T[i]=x[index[i]];
  for (int i=0;i<size;i++) x[i]=T[i];
  for (int i=0;i<size;i++) T[i]=y[index[i]];
  for (int i=0;i<size;i++) y[i]=T[i];
  for (int i=0;i<size;i++) T[i]=z[index[i]];
  for (int i=0;i<size;i++) z[i]=T[i];
  delete [] index;   

  cout << "size : " << size << endl; 
  return size;
}


int
ReadBRST50MPC_S6_ToSkymask(skymap& sm) {

  double Pi = TMath::Pi();

  TFile *ifile = TFile::Open(SKYMASK_FILE_NAME);
  TTree* itree = (TTree *) gROOT->FindObject("MDC");
  int isize = itree->GetEntries();
  //itree->Draw("15*2.5e-21/SimHrss:External_x:External_phi:External_psi","","goff");
  itree->Draw("15*2.5e-21/SimHrss:External_x:External_phi:EarthCtrGPS","","goff");
  double* idistance = itree->GetV1(); 
  double* itheta    = itree->GetV2(); 
  double* iphi      = itree->GetV3(); 
  //double* ipsi      = itree->GetV4(); 
  double* igps      = itree->GetV4(); 

  int nGWGC=0;
  for(int n=0;n<isize;n++) { 

    double dist  = idistance[n];

    if (dist<DISTANCE_THR) {

      double th = 0;
      th = acos(itheta[n]);
      th*= 180/Pi;

      double ph = 0;
      ph = iphi[n] > 0 ? iphi[n] : 2*Pi+iphi[n];
      ph*= 180/Pi;
      ph = sm.phi2RA(ph, igps[n]);

      int ind = sm.getSkyIndex(th,ph);
      //sm.set(ind,sm.get(ind)+1./dist);  // weihtd with the distance
      sm.set(ind,sm.get(ind)+1);  // weihtd with the distance

      nGWGC++;
    }
  }

  return nGWGC;
}


int
ReadBRST50MPC_S6(double*& x, double*& y, double*& z, int& nGWGC) {

  double Pi = TMath::Pi();

  x = new double[NMAX];
  y = new double[NMAX];
  z = new double[NMAX];

  TFile *ifile = TFile::Open(SKYMASK_FILE_NAME);
  TTree* itree = (TTree *) gROOT->FindObject("MDC");
  int isize = itree->GetEntries();
  //itree->Draw("15*2.5e-21/SimHrss:External_x:External_phi:External_psi","","goff");
  itree->Draw("15*2.5e-21/SimHrss:External_x:External_phi:EarthCtrGPS","","goff");
  double* idistance = itree->GetV1(); 
  double* itheta    = itree->GetV2(); 
  double* iphi      = itree->GetV3(); 
  //double* ipsi      = itree->GetV4(); 
  double* igps      = itree->GetV4(); 

  skymap sm;

  nGWGC=0;
  for(int n=0;n<isize;n++) { 

    double dist  = idistance[n];

    if (dist<DISTANCE_THR) {

      double th = 0;
      th = acos(itheta[n]);
      th*= 180/Pi;
      th = -(th-90);

      double ph = 0;
      ph = iphi[n] > 0 ? iphi[n] : 2*Pi+iphi[n];
      ph*= 180/Pi;
      //ph = sm.RA2phi(ph, igps[n]);
      ph = sm.phi2RA(ph, igps[n]);
      ph-=180;
/*
      double psi = 0;
      psi = ipsi[n];
      psi*= 180/Pi;
*/
      x[nGWGC]=ph;
      y[nGWGC]=th;
      z[nGWGC]=dist;
      nGWGC++;
    }
  }
  cout << "nGWGC : " << nGWGC << endl;

  int size=nGWGC;
  for (int i=0;i<360*RESOLUTION;i++) {
    for (int j=0;j<180*RESOLUTION;j++) {
      double ph = i/(double)RESOLUTION;
      double th = j/(double)RESOLUTION;
      x[size]=ph;
      y[size]=th;
      z[size]=DISTANCE_THR;
      size++;
    }
  }

  Int_t *index = new Int_t[size];
  TMath::Sort(size,z,index,true);
  double T[NMAX];
  for (int i=0;i<size;i++) T[i]=x[index[i]];
  for (int i=0;i<size;i++) x[i]=T[i];
  for (int i=0;i<size;i++) T[i]=y[index[i]];
  for (int i=0;i<size;i++) y[i]=T[i];
  for (int i=0;i<size;i++) T[i]=z[index[i]];
  for (int i=0;i<size;i++) z[i]=T[i];
  delete [] index;   

  cout << "size : " << size << endl; 
  return size;
}

