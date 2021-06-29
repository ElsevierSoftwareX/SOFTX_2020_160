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
// Draw PRC Search Area
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

TH1F* sa_hist;
TCanvas* sa_canvas;

void
DrawSearchAreaPRC(TString data_label, TString odir, TString merge_label, int nIFO, float T_win, 
                  int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

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

  itree->Draw("erA[0]",cut,"goff");
  int size = itree->GetSelectedRows();
  cout << "detected events : " << size << endl;

  double* erA0 = itree->GetV1();
  Int_t *index = new Int_t[size];

  TMath::Sort(size,itree->GetV1(),index,false);

  double area_min = pow(erA0[index[0]],2);
  double area_max = pow(erA0[index[size-1]],2);
  area_min = area_min ? area_min : 0.1;
  int nbin = TMath::Ceil(area_max/area_min);

  sa_hist = new TH1F("sa_hist","sa_hist",nbin,area_min,area_max);

  for(int i=0;i<size;i++) {
    double sarea = erA0[index[i]]*erA0[index[i]];
    sa_hist->Fill(sarea);
  }

  // Normalization
  double integral = 0;
  for (int i=0;i<=sa_hist->GetNbinsX();i++) integral+=sa_hist->GetBinContent(i);
  for (int i=0;i<=sa_hist->GetNbinsX();i++) sa_hist->SetBinContent(i,sa_hist->GetBinContent(i)/integral);
  for (int i=2;i<=sa_hist->GetNbinsX();i++) sa_hist->SetBinContent(i,sa_hist->GetBinContent(i)+sa_hist->GetBinContent(i-1));
  // find area at 50%
  double search_area50=0;
  for (int i=0;i<=sa_hist->GetNbinsX();i++) if(sa_hist->GetBinContent(i)>0.5) {search_area50=sa_hist->GetBinCenter(i);break;}
  cout << "search_area50 : " << search_area50 << endl;
  // find area at 90%
  double search_area90=0;
  for (int i=0;i<=sa_hist->GetNbinsX();i++) if(sa_hist->GetBinContent(i)>0.9) {search_area90=sa_hist->GetBinCenter(i);break;}
  cout << "search_area90 : " << search_area90 << endl;


  sa_canvas = new TCanvas("fom", "PRC", 300,40, 600, 600);
  sa_canvas->SetGridx();
  sa_canvas->SetGridy();
  sa_canvas->SetLogx();

  // dump search area

  char ofname[1024];
  ofstream out;
  sprintf(ofname,"%s/search_area.txt",odir.Data());
  out.open(ofname,ios::out);
  if (!out.good()) {cout << "DrawSearchAreaPRC.C: Error Opening Output File : " << ofname << endl;exit(1);}
  cout << "Create Dump Search Area Output File : " << ofname << endl;
  int hsize = sa_hist->GetNbinsX()+1;
  double* fraction = new double[hsize];
  double* perc = new double[hsize];
  double* eperc = new double[hsize];
  char ostring[256];
  for (int i=0;i<=sa_hist->GetNbinsX();i++) {
    // binomial error
    double N = integral;
    double k = integral*sa_hist->GetBinContent(i);
    if(k>N) k=N;
    fraction[i] = sa_hist->GetBinCenter(i);
    perc[i] = sa_hist->GetBinContent(i);
    eperc[i] = (1/N)*sqrt(k*(1-k/N));
    sprintf(ostring,"%.6f\t%.6f\t%.6f\t%d\n",sa_hist->GetBinCenter(i), perc[i], eperc[i], int(N));
    out << ostring;
  }
  out.close();

  TGraphErrors* sa_graph = new TGraphErrors(hsize, fraction, perc, 0, eperc);
  sa_graph->GetXaxis()->SetTitleOffset(1.3);
  sa_graph->GetYaxis()->SetTitleOffset(1.25);
  sa_graph->GetXaxis()->SetTickLength(0.01);
  sa_graph->GetYaxis()->SetTickLength(0.01);
  sa_graph->GetXaxis()->CenterTitle(kTRUE);
  sa_graph->GetYaxis()->CenterTitle(kTRUE);
  sa_graph->GetXaxis()->SetTitleFont(42);
  sa_graph->GetXaxis()->SetLabelFont(42);
  sa_graph->GetYaxis()->SetTitleFont(42);
  sa_graph->GetYaxis()->SetLabelFont(42);
  sa_graph->SetMarkerStyle(20);
  sa_graph->SetMarkerSize(0.2);
  sa_graph->SetMarkerColor(1);
  sa_graph->SetLineColor(kOrange);
  sa_graph->SetLineWidth(3);
  sa_graph->SetTitle("");
  sa_graph->SetFillColor(kOrange);
  sa_graph->SetFillStyle(3001);

  char title[256];
  sprintf(title,"Search Area - Ndet = %d",size); 
  sa_graph->SetTitle(title);
  sa_graph->Paint("alp");
  sa_graph->GetXaxis()->SetTitle("Searched area deg^{2}");
  sa_graph->GetYaxis()->SetTitle("Cumulative fraction of events");
  sa_graph->GetXaxis()->SetLimits(1, 10000);
  sa_graph->GetYaxis()->SetRangeUser(0, 1);
  sa_graph->GetXaxis()->SetTitleOffset(1.3);
  sa_graph->GetYaxis()->SetTitleOffset(1.25);
  sa_graph->GetXaxis()->SetTickLength(0.01);
  sa_graph->GetYaxis()->SetTickLength(0.01);
  sa_graph->GetXaxis()->CenterTitle(kTRUE);
  sa_graph->GetYaxis()->CenterTitle(kTRUE);
  sa_graph->GetXaxis()->SetTitleFont(132);
  sa_graph->GetXaxis()->SetLabelFont(132);
  sa_graph->GetYaxis()->SetTitleFont(132);
  sa_graph->GetYaxis()->SetLabelFont(132);

  sa_graph->Draw("aple3");

  // draw the legend
  TLegend *legend = new TLegend(0.1, 0.75, 0.5, 0.9,NULL,"brNDC");
  legend->SetBorderSize(1);
  legend->SetTextFont(132);
  legend->SetTextSize(0.03);
  legend->SetLineColor(kBlack);
  legend->SetFillColor(kWhite);
  char label[256];
  sprintf(label,"50%% : %0.2f (deg^{2})",search_area50);
  legend->AddEntry(sa_graph,label,"");
  sprintf(label,"90%% : %0.2f (deg^{2})",search_area90);
  legend->AddEntry(sa_graph,label,"");
  legend->Draw();

  sprintf(ofname,"%s/search_area.gif",odir.Data());
  sa_canvas->Update();
  sa_canvas->SaveAs(ofname);
  char cmd[1024];
  TString pfname(ofname);
  pfname.ReplaceAll(".gif",".png");
  sprintf(cmd,"convert %s %s",ofname,pfname.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"rm %s",ofname);
  cout << cmd << endl;
  gSystem->Exec(cmd);

  delete sa_graph;
  delete [] fraction;
  delete [] perc;
  delete [] eperc;
}

