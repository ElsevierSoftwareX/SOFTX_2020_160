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


// draw sensitivity curves

TCanvas* canvas;
TGraph* gr;
TMultiGraph* mg;
  
void cwb_draw_sensitivity(TString ifname="", bool save_plot=false, bool range_fix=true) {

  if(ifname.CompareTo("")==0) {
    if(gSystem->Getenv("CWB_SENSITIVITY_FILE_NAME")!=NULL) {
      ifname=TString(gSystem->Getenv("CWB_SENSITIVITY_FILE_NAME"));
    }
    if(gSystem->Getenv("CWB_SENSITIVITY_SAVE_PLOT")!=NULL) {
      TString cwb_save_plot=TString(gSystem->Getenv("CWB_SENSITIVITY_SAVE_PLOT"));
      if(cwb_save_plot.CompareTo("")!=0) {
        if(!cwb_save_plot.IsFloat()) {cout<< "Error : CWB_SENSITIVITY_SAVE_PLOT is not a number" << endl;exit(1);}
        if((cwb_save_plot.Atoi()==0)||(cwb_save_plot.Atoi()==1)) save_plot=cwb_save_plot.Atoi();
      }
    }
    if(gSystem->Getenv("CWB_SENSITIVITY_RANGE_FIX")!=NULL) {
      TString cwb_range_fix=TString(gSystem->Getenv("CWB_SENSITIVITY_RANGE_FIX"));
      if(cwb_range_fix.CompareTo("")!=0) {
        if(!cwb_range_fix.IsFloat()) {cout<< "Error : CWB_SENSITIVITY_RANGE_FIX is not a number" << endl;exit(1);}
        if((cwb_range_fix.Atoi()==0)||(cwb_range_fix.Atoi()==1)) range_fix=cwb_range_fix.Atoi();
      }
    }
  }
  if(ifname.CompareTo("")==0) {cout << "cwb_draw_sensitivity - Error : File not exist ! " << endl;exit(1);}

  TString stitle=ifname;
  stitle.ReplaceAll(".txt","");
  char title[512];
  TObjArray* token = TString(stitle).Tokenize(TString("/"));
  sprintf(title,"Sensitivity One Side - %s",((TObjString*)token->At(token->GetEntries()-1))->GetString().Data());

  ifstream in;
  in.open(ifname.Data(), ios::in); 
  if (!in.good()) {cout << "Error Opening File : " << ifname.Data() << endl;exit(1);}

  int size=0;
  char str[1024];
  int fpos=0;
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') size++;
  }
  //cout << "size " << size << endl;
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);

  double* freq = new double[size];
  double* sh = new double[size];
  int lines=0;
  while (1) {
    in >> freq[lines] >> sh[lines];
    if (!in.good()) break;
    size=lines;
    lines++;
  }
  in.close();

 
  gStyle->SetTitleH(0.032);
  gStyle->SetTitleW(0.98);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(72);
  gStyle->SetMarkerColor(50);
  gStyle->SetLineColor(kWhite);

  gr = new TGraph(size,freq,sh);
  gr->SetLineColor(kBlue);
  gr->SetLineWidth(1);
  gr->SetMarkerColor(kBlue);
 
  canvas = new TCanvas("Sensitivity", "ShOneSide", 300,40, 1000, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetLogx();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetMarkerStyle(7);
  gStyle->SetMarkerSize(2);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);

  mg = new TMultiGraph;
  mg->SetTitle(title);
  gStyle->SetTitleW(0.98);
  gStyle->SetTitleH(0.03);
  gStyle->SetFillColor(kWhite);
  gStyle->SetLineColor(kWhite);

  mg->Add(gr);
  mg->Paint("APL");

  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.03);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.03);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(80);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(82);
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);

  if(range_fix) {
    mg->GetHistogram()->GetXaxis()->SetRangeUser(16,8*1024);
    mg->GetHistogram()->GetYaxis()->SetRangeUser(5e-24,1e-21);
  }

  mg->GetXaxis()->SetLabelFont(22);
  mg->GetYaxis()->SetLabelFont(22);
  mg->GetXaxis()->SetTitleFont(22);
  mg->GetYaxis()->SetTitleFont(22);
  mg->GetXaxis()->SetTitleSize(0.03);
  mg->GetYaxis()->SetTitleSize(0.03);
  mg->GetXaxis()->SetTitle("Frequency (Hz)");
  mg->GetYaxis()->SetTitle("#frac{1}{#sqrt{Hz}}");
  mg->Draw("APL");

  if(save_plot) {
    TString ofname(ifname);
    ofname.ReplaceAll(".txt",".gif"); 
    cout << ofname.Data() << endl;
    canvas->Print(ofname.Data());
    char cmd[256];
    TString pfname(ofname);
    pfname.ReplaceAll(".gif",".png"); 
    sprintf(cmd,"convert %s %s",ofname.Data(),pfname.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",ofname.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    exit(0);
  }
}
