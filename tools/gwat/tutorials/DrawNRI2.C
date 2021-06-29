//
// Draw Network Response Index skymap
// Author : Gabriele Vedovato

#define L1_ENABLED
#define H1_ENABLED
#define V1_ENABLED
#define H2_ENABLED
//#define G1_ENABLED
//#define T1_ENABLED
//#define A1_ENABLED
//#define A2_ENABLED
//#define J1_ENABLED

#define H1_SIGMA 1
#define L1_SIGMA 1
#define G1_SIGMA 1
#define V1_SIGMA 1
#define T1_SIGMA 1
//#define H2_SIGMA 2
#define H2_SIGMA 1
#define A1_SIGMA 1
#define A2_SIGMA 1
#define J1_SIGMA 1


//#define NTRIES  1000000
#define NTRIES  100000

//#define USE_NOISE

//#define WRITE_PLOT
#define PLOT_POSTFIX  "_WM"

#define SNR 10

#define RESOLUTION  1
//#define RESOLUTION  2
//#define RESOLUTION  4


//#define WRONG_DIRECTION


//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define DISPLAY_WORLD_MAP

//#define DISPLAY_PERC_UNDER2

//#define NET_GAMMA 0.02
#define NET_GAMMA 0.2
//#define NET_GAMMA 0.0


static inline    int net9(double*, double, double*, double, double);
int GetSkyMapSensitivity(double*& x, double*& y, double*& z, double*& t, 
                         TString& title, TString& ofileName, int& ndet);


void DrawNRI2() {

  int ndet=0;
  double *x,*y,*z,*t;
  TString title;
  TString ofileName;

  int size = GetSkyMapSensitivity(x,y,z,t,title,ofileName,ndet);
  if (ndet==0) {cout << "No detector !!!" << endl;exit(0);}

  TH2D* Fx = new TH2D("Fx","Fx", 360*RESOLUTION, 0, 360, 180*RESOLUTION, 0, 180);
  TH2D* h2 = new TH2D("angle","angle", 360*RESOLUTION, 0, 360, 180*RESOLUTION, 0, 180);
  h2->SetStats(kFALSE);

  h2->GetXaxis()->SetNdivisions(70318);
  h2->GetXaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetLabelOffset(0.012);
  h2->GetXaxis()->SetTitleOffset(1.1);
  h2->GetXaxis()->SetTitleFont(72);
  h2->GetYaxis()->SetNdivisions(409);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelOffset(0.01);
  h2->GetZaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetXaxis()->SetTitle("Phi");
  h2->GetXaxis()->CenterTitle(true);
  h2->GetYaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitle("Theta");
  h2->GetYaxis()->CenterTitle(true);
  h2->GetZaxis()->SetLabelOffset(0.00001);
  h2->GetZaxis()->SetNoExponent(false);

  h2->SetTitle(title);
#ifdef DISPLAY_PERC_UNDER2
  int ncnt_under2[360*RESOLUTION][180*RESOLUTION];
  for(int i=0;i<360*RESOLUTION;i++) for(int j=0;j<180*RESOLUTION;j++) ncnt_under2[i][j]=0; 
#endif
  int ncnt[360*RESOLUTION][180*RESOLUTION];
  for(int i=0;i<360*RESOLUTION;i++) for(int j=0;j<180*RESOLUTION;j++) ncnt[i][j]=0; 
  for (int i=0;i<size;i++) {
    int ii=int(x[i]); if(ii>359) ii=359;
    int jj=int(y[i]); if(jj>179) jj=179;
    ii*=RESOLUTION;
    jj*=RESOLUTION;
    ncnt[ii][jj]++;
#ifdef DISPLAY_PERC_UNDER2
    if(z[i]<2) ncnt_under2[ii][jj]++;
    Fx->SetBinContent(ii+1,jj+1,t[i]);
#else
    double binc = h2->GetBinContent(ii+1,jj+1);
    h2->SetBinContent(ii+1,jj+1,binc+z[i]);
#endif
  }
  for(int i=0;i<360*RESOLUTION;i++) for(int j=0;j<180*RESOLUTION;j++) { 
#ifdef DISPLAY_PERC_UNDER2
    if(ncnt[i][j]>0) h2->SetBinContent(i+1,j+1,100.*(double(ncnt_under2[i][j])/double(ncnt[i][j])));
#else
    double binc = h2->GetBinContent(i+1,j+1);
    if(ncnt[i][j]>0) binc/=ncnt[i][j];
    h2->SetBinContent(i+1,j+1,binc);
#endif
  }
  //h2->Draw("colfz");

  size=180*RESOLUTION*360*RESOLUTION;
  double* xx = new double[size];
  double* yy = new double[size];
  double* zz = new double[size];

  int cnt=0;
  for(int i=0;i<360*RESOLUTION;i++) for(int j=0;j<180*RESOLUTION;j++) { 
#ifdef DISPLAY_PERC_UNDER2
    double perc_under2 = h2->GetBinContent(i+1,j+1);
    double fx = Fx->GetBinContent(i+1,j+1);
#else
    double nri = h2->GetBinContent(i+1,j+1);
#endif
    xx[cnt]=i;
    yy[cnt]=j;
#ifdef DISPLAY_PERC_UNDER2
    zz[cnt]=perc_under2;
    //if(perc_under2>90) zz[cnt]=perc_under2;
    //if(perc_under2>0) zz[cnt]=fx;
    //if(perc_under2>0) zz[cnt]=fx*(100.-perc_under2);
    //zz[cnt]=fx;
#else
    zz[cnt]=nri;
#endif
    cnt++;
  }

  gskymap* gSM = new gskymap(int(7));
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION);

#ifdef DISPLAY_WORLD_MAP
  gSM->SetWorldMap();
#endif
  TH2D* hh2 = (TH2D*)gSM->GetHistogram();
//  hh2->GetZaxis()->SetRangeUser(0,1);
  gSM->SetTitle(title);
  if(COORDINATES=="Geographic") {
    for (int i=0;i<size;i++) {
      CwbToGeographic(xx[i],yy[i],xx[i],yy[i]);
    }
  }
  gSM->FillData(size, xx, yy, zz);
//  gSM->SetPalette(1);
  gSM->Draw(0);


#ifdef WRITE_PLOT
  cout << "Print File : " << ofileName.Data() << endl;
  canvas->Print(ofileName);
#endif

}

#define NDETECTORS  9

int
GetSkyMapSensitivity(double*& x, double*& y, double*& z, double*& t, 
                     TString& title, TString& ofileName, int& ndet) {

   bool detector_selected[NDETECTORS];
#ifdef H1_ENABLED
   detector_selected[0]=true;   // LHO1
#else
   detector_selected[0]=false;   // LHO1
#endif
#ifdef L1_ENABLED
   detector_selected[1]=true;   // LLO
#else
   detector_selected[1]=false;   // LLO
#endif
#ifdef G1_ENABLED
   detector_selected[2]=true;  // GEO
#else
   detector_selected[2]=false;  // GEO
#endif
#ifdef V1_ENABLED
   detector_selected[3]=true;   // VIRGO
#else
   detector_selected[3]=false;   // VIRGO
#endif
#ifdef T1_ENABLED
   detector_selected[4]=true;  // TAMA
#else
   detector_selected[4]=false;  // TAMA
#endif
#ifdef H2_ENABLED
   detector_selected[5]=true;   // LHO2
#else
   detector_selected[5]=false;   // LHO2
#endif
#ifdef A1_ENABLED
   detector_selected[6]=true;   // AIGO
#else
   detector_selected[6]=false;   // AIGO
#endif
#ifdef A2_ENABLED
   detector_selected[7]=true;   // AIGO2
#else
   detector_selected[7]=false;   // AIGO2
#endif
#ifdef J1_ENABLED
   detector_selected[8]=true;   // LCGT
#else
   detector_selected[8]=false;   // LCGT
#endif

   TString detector_name[NDETECTORS];
   detector_name[0]="H1";       // LHO1
   detector_name[1]="L1";       // LLO
   detector_name[2]="G1";       // GEO
   detector_name[3]="V1";       // VIRGO
   detector_name[4]="T1";       // TAMA
   detector_name[5]="H2";       // LHO2
   detector_name[6]="A1";       // AIGO
   detector_name[7]="A2";       // AIGO2
   detector_name[8]="J1";       // LCGT

   double detector_sensitivity[NDETECTORS];
   detector_sensitivity[0]=H1_SIGMA;   // LHO1
   detector_sensitivity[1]=L1_SIGMA;   // LLO
   detector_sensitivity[2]=G1_SIGMA;   // GEO
   detector_sensitivity[3]=V1_SIGMA;   // VIRGO
   detector_sensitivity[4]=T1_SIGMA;   // TAMA
   detector_sensitivity[5]=H2_SIGMA;   // LHO2
   detector_sensitivity[6]=A1_SIGMA;   // AIGO
   detector_sensitivity[7]=A2_SIGMA;   // AIGO2
   detector_sensitivity[8]=J1_SIGMA;   // LCGT

   TString sdetectors;
   if (detector_selected[0]) sdetectors+="LHO1 ";
   if (detector_selected[1]) sdetectors+="LLO ";
   if (detector_selected[2]) sdetectors+="GEO ";
   if (detector_selected[3]) sdetectors+="VIRGO ";
   if (detector_selected[4]) sdetectors+="TAMA ";
   if (detector_selected[5]) sdetectors+="LHO2 ";
   if (detector_selected[6]) sdetectors+="AIGO ";
   if (detector_selected[7]) sdetectors+="AIGO2 ";
   if (detector_selected[8]) sdetectors+="LCGT ";

   for (int k=0;k<NDETECTORS;k++) {
     if (detector_selected[k]) {
       ndet++;
     }
   }

   char hacc_title[256];
#ifdef USE_NOISE
   sprintf(hacc_title,"%s Network Response Index (gamma=%2.2f)",sdetectors.Data(),NET_GAMMA);
#else
   sprintf(hacc_title,"%s Network Response Index (gamma=%2.2f)",sdetectors.Data(),NET_GAMMA);
#endif
   title = hacc_title;


   char ofile_title[256];
#ifdef USE_NOISE
   sprintf(ofile_title,"%sNRI%s",sdetectors.Data(),PLOT_POSTFIX);
#else
   sprintf(ofile_title,"%sNRI%s",sdetectors.Data(),PLOT_POSTFIX);
#endif
   if(PROJECTION=="hammer") sprintf(ofile_title,"%s_hammer.png",ofile_title);
   else                     sprintf(ofile_title,"%s.png",ofile_title);

   ofileName = ofile_title;
   ofileName.ReplaceAll(" ","_");

   detector* D[NDETECTORS];

   for(int n=0;n<9;n++) D[n] = new detector((char*)detector_name[n].Data());

   skymap sm(0.4,0,180,0,360);
   int L = sm.size();

   double pi = TMath::Pi();
   wavecomplex sF[NDETECTORS];
   wavecomplex F[NDETECTORS];
   double      X[NDETECTORS];
   int size = NTRIES;
   x = new double[size];
   y = new double[size];
   z = new double[size];
   t = new double[size];
   int cnt=0;
   gRandom->SetSeed(0);
   for (int i=0;i<NTRIES;i++) {
      if(i%100000==0) cout << i << endl;
      int idx = int(gRandom->Uniform(0,L));

      //double phi=sm.getPhi(idx);
      //double theta=sm.getTheta(idx);
      //double psi=0;

      double phi=gRandom->Uniform(0,360);
      double theta=gRandom->Uniform(0,180);
      double psi=gRandom->Uniform(0,180);

      double Hp=gRandom->Uniform(-1,1);
      double Hx=gRandom->Uniform(-1,1);

      for (int k=0;k<NDETECTORS;k++) {
        if (detector_selected[k]) {
          F[k] = D[k]->antenna(theta,phi,psi);
          X[k] = Hp*F[k].real()+Hx*F[k].imag();
        }
      } 
         
      double gp=0;
      double gx=0;
      double gI=0;
      for (int k=0;k<NDETECTORS;k++) {
        if (detector_selected[k]) {
          gp+=F[k].real()*F[k].real();
          gx+=F[k].imag()*F[k].imag();
          gI+=F[k].real()*F[k].imag();
        }
      }
      double gR = (gp-gx)/2.;
      double gr = (gp+gx)/2.;
      double gc = sqrt(gR*gR+gI*gI);           // norm of complex antenna pattern

      double Xp=0;
      double Xx=0;
      for (int k=0;k<NDETECTORS;k++) {
        if (detector_selected[k]) {
          Xp+=X[k]*F[k].real();
          Xx+=X[k]*F[k].imag();
        }
      } 
      double uc = (Xp*gx - Xx*gI);            // u cos of rotation to PCF
      double us = (Xx*gp - Xp*gI);            // u sin of rotation to PCF
      double vc = (gp*uc + gI*us);            // (u*f+)/|u|^2 - 'cos' for v
      double vs = (gx*us + gI*uc);            // (u*fx)/|u|^2 - 'sin' for v
      double u[NDETECTORS];
      double v[NDETECTORS];
      double um=0; 
      double vm=0; 
      for (int k=0;k<NDETECTORS;k++) {
        u[k]=0; 
        v[k]=0; 
        if (detector_selected[k]) {
          u[k]=F[k].real()*uc+F[k].imag()*us;
          um+=u[k]*u[k];
          v[k]=-F[k].imag()*vc+F[k].real()*vs;
          vm+=v[k]*v[k];
        }
      }

// Compute Network Detector Index
      double gamma=NET_GAMMA;
      double GAMMA = 1.-gamma*gamma;                // network regulator
      double ndi=net9(u,um,v,vm,GAMMA);

      x[cnt]=phi;
      y[cnt]=theta;
      z[cnt]=ndi;
      cnt++;
   }
   return size;
}

inline int net9(double* u, double um, double* v, double vm, double g) {
  double q = (1.-g)*um;

  return int(u[0]*u[0]>q) - int((u[0]*u[0]/um+v[0]*v[0]/vm)>g) +
         int(u[1]*u[1]>q) - int((u[1]*u[1]/um+v[1]*v[1]/vm)>g) +
         int(u[2]*u[2]>q) - int((u[2]*u[2]/um+v[2]*v[2]/vm)>g) +
         int(u[3]*u[3]>q) - int((u[3]*u[3]/um+v[3]*v[3]/vm)>g) +
         int(u[4]*u[4]>q) - int((u[4]*u[4]/um+v[4]*v[4]/vm)>g) +
         int(u[5]*u[5]>q) - int((u[5]*u[5]/um+v[5]*v[5]/vm)>g) +
         int(u[6]*u[6]>q) - int((u[6]*u[6]/um+v[6]*v[6]/vm)>g) +
         int(u[7]*u[7]>q) - int((u[7]*u[7]/um+v[7]*v[7]/vm)>g) +
         int(u[8]*u[8]>q) - int((u[8]*u[8]/um+v[8]*v[8]/vm)>g);
}

