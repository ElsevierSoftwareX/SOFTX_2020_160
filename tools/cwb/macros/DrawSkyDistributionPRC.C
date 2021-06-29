/*
# Copyright (C) 2019 Gabriele Vedovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


//
// Draw injected/reconstructed locations of the detected events
// Note : this macro is used to generate the PRC report (cwb_report merge_label prc)
// Author : Gabriele Vedovato

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
#include "TVector3.h"
#include "TRotation.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Rotation3D.h"
#include "constants.hh"

#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

//#define WRITE_PLOT

#define DRAW_EVENTS

//#define DRAW_ANGULAR_DISTANCE
#define ANGULAR_DISTANCE 90

//#define DRAW_PP

#define DISPLAY_WORLD_MAP
#define WORLD_MAP_DIR "$CWB_GWAT/data/"

using namespace CWB;

void 
DrawSkyDistributionPRC(TString data_label, TString odir, TString merge_label, vector<TString> ifo, 
                       detectorParams* dP, float T_win, int pp_inetcc, float T_cor, int pp_irho, float T_cut, 
                       float T_vED, float T_pen, float T_ifar, bool binj=true, TString polarization="TENSOR", bool save_antpat=false) {

  // binj=true  -> select injected directions
  // binj=false -> select reconstructed directions

  if (!gROOT->GetClass("Polar3DVector")) gSystem->Load("libMathCore");

  using namespace ROOT::Math;

  double r2d = 180./TMath::Pi();
  double d2r = TMath::Pi()/180.;

  double speedLight = watconstants::SpeedOfLightInVacuo();

  int nIFO = ifo.size();

  detector* pD[NIFO_MAX];
  for(int n=0;n<nIFO;n++) {
    if(ifo[n]!="") pD[n] = new detector((char*)ifo[n].Data()); 	// built in detector
    else           pD[n] = new detector(dP[n]);			// user detector
  }
  for(int n=0;n<nIFO;n++) cout << n << " " << pD[n]->Name << endl;

  gnetwork* gNET = new gnetwork;
  for(int i=0; i<nIFO; i++) gNET->add(pD[i]);

  if(polarization=="SCALAR") {
    // setting for SGW
    for(int n=0;n<(int)gNET->ifoListSize();n++) {
      detector *d = gNET->getifo(n);
      d->setPolarization(SCALAR);
    }
  }

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

  if(save_antpat) { 				// draw antenna patterns

    char ofileName[1024];
    TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);

    gSM->SetWorldMapPath(world_map.Data());
    gSM->SetWorldMap();
    gNET->DrawAntennaPattern(2);		// |fx|/|f+|
    gNET->DrawSitesShortLabel(kBlack);
    gNET->DrawSites(kBlack,2.0);
    gNET->DrawSitesArms(1000000,kWhite,3.0);
    sprintf(ofileName,"%s/antpat_alignment.png",odir.Data());
    cout << "Write : " << ofileName << endl;
    gSM->Print(ofileName);

    gSM->SetWorldMapPath(world_map.Data());
    gSM->SetWorldMap();
    gNET->DrawAntennaPattern(3);		// sqrt(|F+|^2+|Fx|^2)
    gNET->DrawSitesShortLabel(kBlack);
    gNET->DrawSites(kBlack,2.0);
    gNET->DrawSitesArms(1000000,kWhite,3.0);
    sprintf(ofileName,"%s/antpat_sensitivity.png",odir.Data());
    cout << "Write : " << ofileName << endl;
    gSM->Print(ofileName);
  }

  gSM->SetWorldMap(false);
#ifdef DRAW_PP
  gNET->DrawAntennaPattern(3,0,true,3);
  cout << "order : " << gSM->getOrder() << endl;
  cout << "size : " << gSM->size() << endl;
  int apsize = gSM->size(); 
  double* ap = new double[apsize]; 
  for(int i=0;i<apsize;i++) ap[i] = gSM->get(i);
  //for(int i=0;i<apsize;i++) cout << i << " " << gSM->get(i) << endl;
  Int_t *iap = new Int_t[apsize];
  TMath::Sort(apsize,ap,iap,true);
  //for(int i=0;i<apsize;i++) cout << i << " " << ap[iap[i]] << endl;
#else
  gNET->DrawAntennaPattern(3,0,true);
#endif

  // compute direction of D1 D2 axis -> vector v12
  XYZVector D1(pD[0]->Rv[0],pD[0]->Rv[1],pD[0]->Rv[2]);
  XYZVector D2(pD[1]->Rv[0],pD[1]->Rv[1],pD[1]->Rv[2]);
  D1=D1/speedLight;
  D2=D2/speedLight;
  XYZVector D12 = D1-D2;
  TVector3 vD12(D12.X(),D12.Y(),D12.Z());

  double th12 = r2d*vD12.Theta();
  double ph12 = r2d*vD12.Phi();
  cout << "coordinates D12 " << ph12 << " " << th12 << endl;
  Polar3DVector v12(1, d2r*th12, d2r*ph12);

  // open reconstructed event file
#ifdef DRAW_EVENTS

  char wave_file_name[1024];
  sprintf(wave_file_name,"merge/wave_%s.%s.root",data_label.Data(),merge_label.Data());

  TFile* ifile = TFile::Open(wave_file_name);
  if(ifile==NULL) {cout<<"Error opening file : "<<wave_file_name<<endl;exit(1);}
  TTree* itree = (TTree *) gROOT->FindObject("waveburst");
  if(itree==NULL) {cout<<"Error opening tree : "<<"waveburst"<<endl;exit(1);}

  char cut[1024];
  char tmp[1024];
  sprintf(cut,"abs(time[0]-time[%d])<%f && netcc[%d]>%f && rho[%d]>%f",
          nIFO,T_win,pp_inetcc,T_cor,pp_irho,T_cut);
  if(T_vED>0)  {strcpy(tmp,cut);sprintf(cut,"%s && neted[0]/ecor<%f",tmp,T_vED);}
  if(T_pen>0)  {strcpy(tmp,cut);sprintf(cut,"%s && penalty>%f",tmp,T_pen);}
  if(T_ifar>0) {strcpy(tmp,cut);sprintf(cut,"%s && ifar>(24.*3600.*365.)*%f",tmp,T_ifar);}

  itree->Draw("theta[0]:phi[0]:theta[1]:phi[1]",cut,"goff");
  int isize=itree->GetSelectedRows();
  cout << "isize : " << isize << endl;

  double* th0 = itree->GetV1();
  double* ph0 = itree->GetV2();
  double* th1 = itree->GetV3();
  double* ph1 = itree->GetV4();

#ifdef DRAW_ANGULAR_DISTANCE
  TH1F* h1 = new TH1F("hist","hist",100,0,180);
#endif

  // create gskymap to store events 
  gskymap* gEVT = new gskymap((int)3);
  gEVT->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
  *gEVT = 0;

  // draw events
  double markerSize = isize>10000 ? 0.3 : 0.4;
  double ph,th;
  for (int i=0;i<isize;i++) {

    // compute distance (dOmega) between injected and reconstructed directions
    Polar3DVector v0(1, d2r*th0[i], d2r*ph0[i]);
    Polar3DVector v1(1, d2r*th1[i], d2r*ph1[i]);
    //cout << th0[i] << " " << ph0[i] << " " << th1[i] << " " << ph1[i] << endl;
    double dot = v0.Dot(v1);
    double dOmega = r2d*TMath::ACos(dot);

    // compute distance (dOmega12) between injected and D1 D2 axis
    double dot12 = v12.Dot(v1);
    double dOmega12 = r2d*TMath::ACos(dot12);

    // discart events outside the ring (width=5degrees) at distance ANGULAR_DISTANCE 
//    if(fabs(dOmega12-ANGULAR_DISTANCE)>5) continue;

#ifdef DRAW_ANGULAR_DISTANCE
    h1->Fill(dOmega);
#endif

    if(binj) {ph=ph1[i]; th=th1[i];}
    else     {ph=ph0[i]; th=th0[i];}

    if(COORDINATES=="cWB") {
      gSM->DrawMarker(ph, th, 20, markerSize, kBlack);  	  // cWB
    }
    if(COORDINATES=="Geographic") {
      double phi,theta;
      CwbToGeographic(ph,th,phi,theta);
      gSM->DrawMarker(phi,theta, 20, markerSize, kBlack);  // Geographic

      int index = gEVT->getSkyIndex(th,ph);
      gEVT->set(index,gEVT->get(index)+1);
    }
  }
#endif

  // draw circle at distance -ANGULAR_DISTANCE from D1 D2 axis
  double phi=ph12;
  double theta=th12-ANGULAR_DISTANCE;
  if(COORDINATES=="Geographic") {
    CwbToGeographic(phi,theta,phi,theta);  
  }
  gNET->DrawCircles(phi,theta,(Color_t)kWhite,1,1,true);

  gSM->GetCanvas()->Update();

  // draw angular distance
#ifdef DRAW_ANGULAR_DISTANCE
  gStyle->SetLineColor(kBlack);
  h1->Draw("HIST");
#endif

  // write plot
  char ofname[1024];
  if(binj) {
    sprintf(ofname,"%s/inj_detected_antpat.png",odir.Data());
  } else {
    sprintf(ofname,"%s/rec_detected_antpat.png",odir.Data());
  }
  cout << "Write : " << ofname << endl;
  gSM->Print(ofname);
  //gSystem->Exit(0);

#ifdef DRAW_PP
  TCanvas* canvas = new TCanvas;
  canvas->SetGridx();
  canvas->SetGridy();

  double* x = new double[apsize];
  double* y = new double[apsize];
  for(int i=0;i<apsize;i++) {
    //x[i]=ap[iap[i]];
    x[i]=pow(ap[iap[i]],4);
    y[i]=gEVT->get(iap[i]);
//x[i]=pow(ap[i],4);
//y[i]=gEVT->get(i);
  }
  for(int i=1;i<apsize;i++) {
    x[i]+=x[i-1];
    y[i]+=y[i-1];
  }
  for(int i=0;i<apsize;i++) {
    x[i]/=x[apsize-1];
    y[i]/=y[apsize-1];
  }
  TGraph* gr = new TGraph(apsize,x,y);
  gr->Draw("ALP");
  gr->SetLineColor(kRed);
  gr->SetMarkerColor(kRed);

  TH1F* hist = gr->GetHistogram();
  hist->SetStats(kFALSE);
  char title[256];
  sprintf(title,"%s - pp plot",TITLE);
  hist->SetTitle(title);
  //hist->GetXaxis()->SetTitle("(F+^{2}+Fx^{2}) Probability");
  hist->GetXaxis()->SetTitle("(F+^{2}+Fx^{2})^{2} Probability");
  hist->GetYaxis()->SetTitle("Frequentist Probability");
  hist->GetXaxis()->SetRangeUser(0,1);
  hist->GetYaxis()->SetRangeUser(0,1);

//  for(int i=0;i<apsize;i++) cout << i << " " << ap[iap[i]] << " " << gEVT->get(iap[i]) << endl;
//  gEVT->Draw();
#endif

}

