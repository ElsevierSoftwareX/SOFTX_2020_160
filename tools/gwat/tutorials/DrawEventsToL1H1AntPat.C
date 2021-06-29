#define XIFO 4
#pragma GCC system_header

#include "constants.hh"
#include "gnetwork.hh"

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TComplex.h>
#include <TStyle.h>
#include <TRandom.h>
#include <complex.h>
#include "TVector3.h"
#include "TRotation.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Rotation3D.h"
#include "Math/Polar3D.h"

//
// Draw Events & Antenna Pattern for L1H1 network
// Author : Gabriele Vedovato


#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

//#define WRITE_PLOT


#define IFILE_NAME "merge/wave_sinegaussian_noVirgo.M4.root"
#define DRAW_EVENTS

//#define DRAW_ANGULAR_DISTANCE
#define ANGULAR_DISTANCE 90

//using namespace CWB;

void DrawEventsToL1H1AntPat(bool binj=true, int polarization=3) {

  // binj=true  -> select injected directions
  // binj=false -> select reconstructed directions

  // polarization is the antenna pattern type 

  if (!gROOT->GetClass("Polar3DVector")) gSystem->Load("libMathCore");

  using namespace ROOT::Math;

  double r2d = 180./TMath::Pi();
  double d2r = TMath::Pi()/180.;

  double speed_light = watconstants::SpeedOfLightInVacuo();

  // define L1H1 network
  int nIFO=2;
  TString ifo[2]={"L1","H1"};
  
  gnetwork* gNET = new gnetwork;

  detector* pD[2];
  for(int i=0; i<nIFO; i++) pD[i] = new detector((char*)ifo[i].Data()); // built in detector
  for(int i=0; i<nIFO; i++) gNET->add(pD[i]);

  // setup skymap options
  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);


  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleFont(42);
  h2->GetZaxis()->SetRangeUser(0,1.0);

  // draw antenna pattern
  gNET->DrawAntennaPattern(polarization);

  // compute direction of D1 D2 axis -> vector v12
  XYZVector D1(pD[0]->Rv[0],pD[0]->Rv[1],pD[0]->Rv[2]);
  XYZVector D2(pD[1]->Rv[0],pD[1]->Rv[1],pD[1]->Rv[2]);
  D1=D1/speed_light;
  D2=D2/speed_light;
  XYZVector D12 = D1-D2;
  TVector3 vD12(D12.X(),D12.Y(),D12.Z());

  double th12 = r2d*vD12.Theta();
  double ph12 = r2d*vD12.Phi();
  cout << "coordinates D12 " << ph12 << " " << th12 << endl;
  Polar3DVector v12(1, d2r*th12, d2r*ph12);

  // open reconstructed event file
#ifdef DRAW_EVENTS
  TFile *ifile = TFile::Open(IFILE_NAME);
  if(ifile==NULL) {cout<<"Error opening file : "<<IFILE_NAME<<endl;exit(1);}
  TTree* itree = (TTree *) gROOT->FindObject("waveburst");
  if(itree==NULL) {cout<<"Error opening tree : "<<"waveburst"<<endl;exit(1);}
  itree->Draw("theta[0]:phi[0]:theta[1]:phi[1]","abs(time[0]-time[2])<0.1 && netcc[0]>0.7 && rho[0]>7 ","goff");
  int isize=itree->GetSelectedRows();
  cout << "isize : " << isize << endl;

  double* th0 = itree->GetV1();
  double* ph0 = itree->GetV2();
  double* th1 = itree->GetV3();
  double* ph1 = itree->GetV4();

#ifdef DRAW_ANGULAR_DISTANCE
  TH1F* h1 = new TH1F("hist","hist",100,0,180);
#endif

  // draw events
  double ph,th;
  for (int i=0;i<isize;i++) {

    // compute distance (dOmega) between injected and reconstructed directions
    Polar3DVector v0(1, d2r*th0[i], d2r*ph0[i]);
    Polar3DVector v1(1, d2r*th1[i], d2r*ph1[i]);
    double dot = v0.Dot(v1);
    double dOmega = r2d*TMath::ACos(dot);

    // compute distance (dOmega12) between injected and D1 D2 axis
    double dot12 = v12.Dot(v1);
    double dOmega12 = r2d*TMath::ACos(dot12);

    // discart events outside the ring (width=5degrees) at distance ANGULAR_DISTANCE 
    if(fabs(dOmega12-ANGULAR_DISTANCE)>5) continue;

#ifdef DRAW_ANGULAR_DISTANCE
    h1->Fill(dOmega);
#endif

    if(binj) {ph=ph1[i]; th=th1[i];}
    else     {ph=ph0[i]; th=th0[i];}

    if(COORDINATES=="cWB") {
      gSM->DrawMarker(ph, th, 20, 0.4, kBlack);  	  // cWB
    }
    if(COORDINATES=="Geographic") {
      double phi,theta;
      CwbToGeographic(ph,th,phi,theta);
      gSM->DrawMarker(phi,theta, 20, 0.4, kBlack);  // Geographic
    }
  }
#endif

  // draw circle at distance -ANGULAR_DISTANCE from D1 D2 axis
  double phi=ph12;
  double theta=th12-ANGULAR_DISTANCE;
  if(COORDINATES=="Geographic") CwbToGeographic(phi,theta,phi,theta);
  gNET->DrawCircles(phi,theta,(Color_t)kWhite,1,1,true);

  gSM->GetCanvas()->Update();

  // draw angular distance
#ifdef DRAW_ANGULAR_DISTANCE
  gStyle->SetLineColor(kBlack);
  h1->Draw();
#endif

  // write plot
#ifdef WRITE_PLOT
  TObjArray* token = TString(IFILE_NAME).Tokenize(TString('/'));
  TObjString* sfile = (TObjString*)token->At(token->GetEntries()-1);
  TString TITLE = sfile->GetString();
  TString ofile = sfile->GetString();
  ofile.ReplaceAll(".root","_EventsVsAntPat.png");
  cout << "Write : " << ofile << endl;
  gSM->Print(ofile);
  gSystem->Exit(0);
#endif
}

