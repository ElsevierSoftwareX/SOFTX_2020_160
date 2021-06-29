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
// Draw PRC Cos Omega distribution
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

TH1F* co_hist;
TCanvas* co_canvas;

void
DrawCosOmegaPRC(TString data_label, TString odir, TString merge_label, int nIFO, float T_win, 
                int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  if (!gROOT->GetClass("Polar3DVector")) gSystem->Load("libMathCore");

  using namespace ROOT::Math;

  char wave_file_name[1024];
  sprintf(wave_file_name,"merge/wave_%s.%s.root",data_label.Data(),merge_label.Data());

  TFile* ifile = TFile::Open(wave_file_name);
  TTree* itree = (TTree *) gROOT->FindObject("waveburst");   

  char cut[1024];
  char tmp[1024];
  sprintf(cut,"abs(time[0]-time[%d])<%f && netcc[%d]>%f && rho[%d]>%f",
          nIFO,T_win,pp_inetcc,T_cor,pp_irho,T_cut);
  if(T_vED>0)  {strcpy(tmp,cut);sprintf(cut,"%s && neted[0]/ecor<%f",tmp,T_vED);}
  if(T_pen>0)  {strcpy(tmp,cut);sprintf(cut,"%s && penalty>%f",tmp,T_pen);}
  if(T_ifar>0) {strcpy(tmp,cut);sprintf(cut,"%s && ifar>(24.*3600.*365.)*%f",tmp,T_ifar);}

  itree->Draw("theta[0]:phi[0]:theta[1]:phi[1]",cut,"goff");
  int size = itree->GetSelectedRows();
  cout << "detected events : " << size << endl;

  double* th0 = itree->GetV1();
  double* ph0 = itree->GetV2();
  double* th1 = itree->GetV3();
  double* ph1 = itree->GetV4();

  co_hist = new TH1F("co_hist","co_hist",24,-1,1);
  //co_hist = new TH1F("co_hist","co_hist",100,0,180);

  double deg2rad = TMath::Pi()/180.;

  for(int i=0;i<size;i++) {

    //cout << i << " " << th0[i] << " " << ph0[i] << endl;
    //cout << i << " " << th1[i] << " " << ph1[i] << endl;

    Polar3DVector v0(1, th0[i]*deg2rad, ph0[i]*deg2rad);
    Polar3DVector v1(1, th1[i]*deg2rad, ph1[i]*deg2rad);

    double Dot = v0.Dot(v1);
    double dOmega = 180.*TMath::ACos(Dot)/TMath::Pi();

/*
    double inj_phi = ph1[i]*deg2rad;
    double est_phi = ph0[i]*deg2rad;
    double inj_theta = th1[i]*deg2rad;
    double est_theta = th0[i]*deg2rad;
    cos_delta_theta = cos(inj_theta)*cos(est_theta) + sin(inj_theta)*sin(est_theta)*cos(inj_phi-est_phi);
    //co_hist->Fill(cos_delta_theta);
*/

    //cout << i << "Dot : " << Dot << " dOmega : " << dOmega << endl;
    co_hist->Fill(Dot);
    //co_hist->Fill(dOmega);
  }

  // Normalization
  double integral = 0;
  for (int i=0;i<=co_hist->GetNbinsX();i++) integral+=co_hist->GetBinContent(i);
  for (int i=0;i<=co_hist->GetNbinsX();i++) co_hist->SetBinContent(i,co_hist->GetBinContent(i)/integral);

  char title[256];
  sprintf(title,"angle offset : #events = %d",size);
  co_hist->SetTitle(title);
  co_hist->GetXaxis()->SetTitle("cos(angle offset)");
  co_hist->GetYaxis()->SetTitle("probability density");
  co_hist->GetXaxis()->SetRangeUser(-1,1);
  co_hist->GetYaxis()->SetRangeUser(0.001,1);
//  co_hist->GetYaxis()->SetRangeUser(1e-3,1);
  co_hist->GetXaxis()->SetTitleOffset(1.35);
  co_hist->GetYaxis()->SetTitleOffset(1.35);
  co_hist->GetXaxis()->CenterTitle(true);
  co_hist->GetYaxis()->CenterTitle(true);
  co_hist->GetXaxis()->SetLabelFont(42);
  co_hist->GetXaxis()->SetTitleFont(42);
  co_hist->GetYaxis()->SetLabelFont(42);
  co_hist->GetYaxis()->SetTitleFont(42);
  co_hist->SetLineColor(kRed);
  co_hist->SetFillColor(kRed);

  co_hist->SetStats(kFALSE);
  co_canvas = new TCanvas("fom", "PRC", 300,40, 600, 500);
  co_canvas->SetGridx();
  co_canvas->SetGridy();
  co_canvas->SetLogy();
  co_hist->Draw("HIST");  

  char ofname[1024];
  sprintf(ofname,"%s/cos_theta.gif",odir.Data());
  co_canvas->Print(ofname);
  char cmd[1024];
  TString pfname(ofname);
  pfname.ReplaceAll(".gif",".png");
  sprintf(cmd,"convert %s %s",ofname,pfname.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"rm %s",ofname);
  cout << cmd << endl;
  gSystem->Exec(cmd);
  //exit(0);

}

