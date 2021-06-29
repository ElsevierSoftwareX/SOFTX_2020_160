//
// Generate GWGC SkyMask 
// Author : Gabriele Vedovato


// -----------------------------------------
// defines GWGC catalog
// -----------------------------------------

#define GWGCCatalog "$CWB_GWAT/data/GWGCCatalog_Rev1d8.txt"
#define DISTANCE_THR 50 // define the maximum distance of galaxies in Mpc

#define NMAX 2136800  // Maximum number of GWGC entries

// -----------------------------------------
// defines skymask options
// -----------------------------------------

// define where to store the output skymask
#define SKYMASK_FILE_NAME  "input/SkyMaskCC_GWGC_Rev1d8_XXX_PAD.txt"

#define SKYMASK_RESOLUTION 0.4	// define the sky resolution of the built-in sky segm
#define HEALPIX_ORDER 7		// if commented the built-in sky segmentaion is used

#define ADD_PADDING	// padding skymask 
//#define OUTPUT_PROBABILITY_SKYMAP

// -----------------------------------------
// defines draw options
// -----------------------------------------
//#define PLOT_GWGC		// plot GWGC skymap distribution
//#define SAVE_GWGC_PLOT	// save plot GWGC
//#define PLOT_SKYMASK		// plot skymask

#define COORDINATES "Geographic"
#define PROJECTION ""
//#define PROJECTION "hammer"
#define RESOLUTION  4


void CreateSkyMaskGWGC() {

  double *x,*y,*z;

  x = new double[NMAX];
  y = new double[NMAX];
  z = new double[NMAX];

  int nGWGC=0;
  int size = ReadGWGCCatalog(x,y,z,nGWGC);


#ifdef HEALPIX_ORDER
  skymap sm(int(HEALPIX_ORDER));
#else
  skymap sm(SKYMASK_RESOLUTION,0,180,0,360);
#endif
  int L = sm.size();
  nGWGC=ReadGWGCCatalogToSkymap(sm);

#ifdef ADD_PADDING
  skymap smpad(sm);
  wavearray<int> neighbors;
  for (int l=0;l<L;l++) {
    if(sm.get(l)) {
      double th = sm.getTheta(l);
      double ph = sm.getPhi(l);
      double sp = sm.getPhiStep(l);
      double st = sm.getThetaStep(l);

#ifdef HEALPIX_ORDER
      // pixels are rhombus
      // all neighbors pixels are used for padding (8 pixels)
      neighbors=sm.neighbors(l);
      for(int k=0;k<neighbors.size();k++) {
        int ll = neighbors[k];
        if(ll>=0) smpad.set(ll,sm.get(ll)+1);
      }
#else
      // pixels are rectangular with base parallel to latitude circles
      // only left/right pixels are used for padding
      int ll = sm.getSkyIndex(th,ph-sp);    // neighbour sky index
      int lr = sm.getSkyIndex(th,ph+sp);    // neighbour sky index
      smpad.set(ll,sm.get(ll)+1);
      smpad.set(lr,sm.get(lr)+1);
#endif
    } 
  } 
  for (int l=0;l<L;l++) sm.set(l,smpad.get(l));
#endif

#ifdef PLOT_SKYMASK
  gskymap* gSMPAD = new gskymap(smpad);
  gSMPAD->SetOptions(PROJECTION,COORDINATES);
  gSMPAD->SetTitle("GWGC skymask");
  gSMPAD->Draw(0);
#endif

#ifdef SKYMASK_FILE_NAME

  TString ofile = SKYMASK_FILE_NAME;
  char skyres[256];
#ifdef HEALPIX_ORDER
  sprintf(skyres,"HPX%d",HEALPIX_ORDER);
#else
  sprintf(skyres,"R%1.2f",SKYMASK_RESOLUTION);
#endif
  ofile.ReplaceAll("XXX",TString(skyres).ReplaceAll(".","d").Data());

  ofstream out;
  out.open(ofile.Data(), ios::out);
  if (!out.good()) {cout << "Error Opening File : " << ofile.Data() << endl;exit(1);}
  cout << "Write File : " << ofile.Data() << endl;
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

#ifdef PLOT_GWGC
  char title[256];
  sprintf(title,"Statistics of sources in the GWGC catalog out to a %d Mpc measured distance (# %d)",DISTANCE_THR,nGWGC);

#ifdef HEALPIX_ORDER
  gskymap* gSM = new gskymap(int(HEALPIX_ORDER));
#else
  gskymap* gSM = new gskymap(SKYMASK_RESOLUTION,0,180,0,360);
#endif
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
  return;

#ifdef SAVE_GWGC_PLOT
  char ofileName[256];
  sprintf(ofileName,"GWGCCatalog_Rev1d7_Distance%dMpc%s.PNG",DISTANCE_THR,PROJECTION);
  cout << "Write : " << ofileName << endl;
  gSM->Print(ofileName);
#endif
#endif

  delete [] x;
  delete [] y;
  delete [] z;
}

int
ReadGWGCCatalogToSkymap(skymap& sm) {

  ifstream in;
  in.open(gSystem->ExpandPathName(GWGCCatalog),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << gSystem->ExpandPathName(GWGCCatalog) << endl;exit(1);}

  char iline[1024];
  in.getline(iline,1024);  // ski first line (header)
  int nGWGC=0;
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

      double th = dec;
      double ph = ra*360./24.;
      CelestialToCwb(ph,th,ph,th);
      int ind = sm.getSkyIndex(th,ph);
      sm.set(ind,sm.get(ind)+1);

      nGWGC++;
    }
  }
  in.close();

  return nGWGC;
}

int
ReadGWGCCatalog(double* x, double* y, double* z, int& nGWGC) {

  ifstream in;
  in.open(gSystem->ExpandPathName(GWGCCatalog),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << gSystem->ExpandPathName(GWGCCatalog) << endl;exit(1);}

  char iline[1024];
  in.getline(iline,1024);  // ski first line (header)
  int nGWGC=0;
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

      x[nGWGC]=ph;
      y[nGWGC]=th;
      z[nGWGC]=dist;
      nGWGC++;
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

