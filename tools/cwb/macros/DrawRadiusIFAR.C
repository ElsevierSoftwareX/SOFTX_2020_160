/*
# Copyright (C) 2019 Francesco Salemi
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
/*
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
*/

TCanvas* co_canvas;

void DrawRadiusIFAR(std::vector<double> vR, std::vector<double> veR, std::vector<double> vsifar, std::vector<double> vseifar, char* netdir){

  
  co_canvas = new TCanvas("sd", "SD", 3, 47, 1000, 802);
  co_canvas->SetGridx();
  co_canvas->SetGridy();
  co_canvas->SetLogx();

  TGraphErrors *co_gr = new TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, &veR[0]);
  co_gr->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
  co_gr->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
  co_gr->GetXaxis()->SetLimits(0.01,1000.);
  co_gr->GetYaxis()->SetRangeUser(500,800);
  co_gr->GetXaxis()->SetTitleOffset(1.3);
  co_gr->GetYaxis()->SetTitleOffset(1.25);
  co_gr->GetXaxis()->SetTickLength(0.01);
  co_gr->GetYaxis()->SetTickLength(0.01);
  co_gr->GetXaxis()->CenterTitle(kTRUE);
  co_gr->GetYaxis()->CenterTitle(kTRUE);
  co_gr->GetXaxis()->SetTitleFont(42);
  co_gr->GetXaxis()->SetLabelFont(42);
  co_gr->GetYaxis()->SetTitleFont(42);
  co_gr->GetYaxis()->SetLabelFont(42);
  co_gr->SetMarkerStyle(20);
  co_gr->SetMarkerSize(1.0);
  co_gr->SetMarkerColor(1);
  co_gr->SetLineColor(kBlue);
  co_gr->SetTitle("");
  co_gr->SetFillColor(kBlue);
  co_gr->SetFillStyle(3001);
 // h->Draw("e2same"); 
  co_gr->Draw("aple3");
 

 

  char fname[1024];  
  sprintf(fname, "%s/ROC_IFAR.png", netdir);
  co_canvas->Update();
  co_canvas->SaveAs(fname);

    //exit(0);

}

