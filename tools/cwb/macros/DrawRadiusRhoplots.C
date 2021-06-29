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
// Draw Radius vs IFAR
// Note : this macro is used to generate the CBC report
// Author : Francesco Salemi


#define MINRADIUS 0
#define MAXRADIUS 1000
#define MINRHO 5
#define MAXRHO 100

//TCanvas* co_canvas;
//TGraphErrors* co_gr1,co_gr2,co_gr3,co_gr4;
//TMultiGraph* mg, mgr, mg2, mg2r, mg3, mg3r;
//TLegend* leg, leg2, leg3;
//void DrawRadiusIFAR2(char* sim_file_name, char* mdc_file_name, float shell_volume, char* netdir, double liveTot, float T_ifar=0.0, float T_win=0.2, float internal_volume=0.0, int TRIALS = 1){
void DrawRadiusRhoplots(char* sim_file_name, char* mdc_file_name, float shell_volume, TString opt){

  TCanvas* co_canvas = new TCanvas("sd2", "SD2", 3, 47, 1000, 802);
  co_canvas->SetGridx();
  co_canvas->SetGridy();
  co_canvas->SetLogx();

  //gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/CreateGraphRadiusIFAR.C"));


  //========== Sensitive Distance vs rho[1] for varius Mtot bins ================
  
  TMultiGraph* mgr = new TMultiGraph();
  TGraphErrors* co_gr1 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, "(mass[0]+mass[1] > 4.0) && (mass[0]+mass[1] < 27.25)", shell_volume,kBlue, opt+"rho");
  mgr->Add(co_gr1);
  TGraphErrors* co_gr2 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, "(mass[0]+mass[1] > 27.25) && (mass[0]+mass[1] < 51.50)", shell_volume,kCyan, opt+"rho");
  mgr->Add(co_gr2);
  TGraphErrors* co_gr3 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, "(mass[0]+mass[1] > 51.50) && (mass[0]+mass[1] < 75.75)", shell_volume,kGreen-9, opt+"rho");
  mgr->Add(co_gr3);
  TGraphErrors* co_gr4 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, "(mass[0]+mass[1] > 75.75) && (mass[0]+mass[1] < 100.0)", shell_volume,kOrange, opt+"rho");
  mgr->Add(co_gr4);
  
  mgr->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
  mgr->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
  mgr->GetXaxis()->SetLimits(MINRHO,MAXRHO);
 // mgr->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
  mgr->GetXaxis()->SetTitleOffset(1.3);
  mgr->GetYaxis()->SetTitleOffset(1.25);
  mgr->GetXaxis()->SetTickLength(0.01);
  mgr->GetYaxis()->SetTickLength(0.01);
  mgr->GetXaxis()->CenterTitle(kTRUE);
  mgr->GetYaxis()->CenterTitle(kTRUE);
  mgr->GetXaxis()->SetTitleFont(42);
  mgr->GetXaxis()->SetLabelFont(42);
  mgr->GetYaxis()->SetTitleFont(42);
  mgr->GetYaxis()->SetLabelFont(42);

  mgr->Draw("aple3");
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9, "", "brNDC");
 leg->AddEntry(co_gr1,"M_{total} #in [4.00, 27.25] M_{#odot}", "l");
 leg->AddEntry(co_gr2,"M_{total} #in [27.25, 51.50] M_{#odot}", "l");
 leg->AddEntry(co_gr3,"M_{total} #in [51.50, 75.75] M_{#odot}", "l");
 leg->AddEntry(co_gr4,"M_{total} #in [75.75, 100.00] M_{#odot}", "l");
 leg->SetFillColor(0);
  leg->SetFillColorAlpha(0, 1.0);
  leg->Draw();
  
  sprintf(fname, "%s/ROC_rho1_Mtot.png", netdir);
  co_canvas->SetLogy();
  co_canvas->Update();
  co_canvas->SaveAs(fname);
  co_canvas->SetLogy(0);
  
  //========== Sensitive Distance vs rho[1] for varius Xeff bins ================
  TMultiGraph* mg2r = new TMultiGraph();
  co_gr1 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, 
	"(spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])>-1 && (spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])<-0.4", 		shell_volume,kViolet-1, opt+"rho");
  mg2r->Add(co_gr1);
  co_gr2 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, 
	"(spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])>-0.4 && (spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])<0.4", 		shell_volume,kGreen-9, opt+"rho");
  mg2r->Add(co_gr2);
  co_gr3 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, 
	"(spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])>0.4 && (spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])<1", 		shell_volume,kRed, opt+"rho");
  mg2r->Add(co_gr3);

  mg2r->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
  mg2r->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
  mg2r->GetXaxis()->SetLimits(MINRHO,MAXRHO);
  //mg2r->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
  mg2r->GetXaxis()->SetTitleOffset(1.3);
  mg2r->GetYaxis()->SetTitleOffset(1.25);
  mg2r->GetXaxis()->SetTickLength(0.01);
  mg2r->GetYaxis()->SetTickLength(0.01);
  mg2r->GetXaxis()->CenterTitle(kTRUE);
  mg2r->GetYaxis()->CenterTitle(kTRUE);
  mg2r->GetXaxis()->SetTitleFont(42);
  mg2r->GetXaxis()->SetLabelFont(42);
  mg2r->GetYaxis()->SetTitleFont(42);
  mg2r->GetYaxis()->SetLabelFont(42);

  mg2r->Draw("aple3");

  TLegend* leg2 = new TLegend(0.6, 0.8, 0.9, 0.9, "", "tlNDC");
  leg2->AddEntry(co_gr1,"#chi_{eff} #in [-1, 0.4]", "l");
  leg2->AddEntry(co_gr2,"#chi_{eff} #in [-0.4, 0.4]", "l");
  leg2->AddEntry(co_gr3,"#chi_{eff} #in [0.4, 1]", "l");
  leg2->SetFillColor(0);

  leg2->SetFillColorAlpha(0, 1.0);
  leg2->Draw();
  
  sprintf(fname, "%s/ROC_rho1_chieff.png", netdir);
  co_canvas->SetLogy();
  co_canvas->Update();
  co_canvas->SaveAs(fname);
  co_canvas->SetLogy(0);


	//========== Sensitive Distance vs rho[1] for varius Mchirp bins ================
  TMultiGraph* mg3r = new TMultiGraph();
  co_gr1 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, 
	"chirp[0]>1.74 && chirp[0]<8.07",shell_volume,kBlue, opt+"rho");
  mg3r->Add(co_gr1);
  co_gr2 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, 
	"chirp[0]>8.07 && chirp[0]<14.92",shell_volume,kCyan, opt+"rho");
  mg3r->Add(co_gr2);
  co_gr3 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, 
	"chirp[0]>14.92 && chirp[0]<21.77",shell_volume,kGreen-9, opt+"rho");
  mg3r->Add(co_gr3);
  co_gr4 = CreateGraphRadiusIFAR(sim_file_name, mdc_file_name, 
	"chirp[0]>21.77 && chirp[0]<100.0",shell_volume,kOrange, opt+"rho");
  mg3r->Add(co_gr4);

  mg3r->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
  mg3r->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
  mg3r->GetXaxis()->SetLimits(MINRHO,MAXRHO);
 // mg3r->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
  mg3r->GetXaxis()->SetTitleOffset(1.3);
  mg3r->GetYaxis()->SetTitleOffset(1.25);
  mg3r->GetXaxis()->SetTickLength(0.01);
  mg3r->GetYaxis()->SetTickLength(0.01);
  mg3r->GetXaxis()->CenterTitle(kTRUE);
  mg3r->GetYaxis()->CenterTitle(kTRUE);
  mg3r->GetXaxis()->SetTitleFont(42);
  mg3r->GetXaxis()->SetLabelFont(42);
  mg3r->GetYaxis()->SetTitleFont(42);
  mg3r->GetYaxis()->SetLabelFont(42);

  mg3r->Draw("aple3");
 
  //char lab[256];
  TLegend* leg3 = new TLegend(0.6, 0.7, 0.9, 0.9, "", "tlNDC");
  leg3->AddEntry(co_gr1,"M_{chirp} #in [1.74, 8.07] M_{#odot}", "l");
  leg3->AddEntry(co_gr2,"M_{chirp} #in [8.07, 14.92] M_{#odot}", "l");
  leg3->AddEntry(co_gr3,"M_{chirp} #in [14.92, 21.77] M_{#odot}", "l");
  leg3->AddEntry(co_gr4,"M_{chirp} #in [21.77, 100.00] M_{#odot}", "l");
  leg3->SetFillColor(0);
  leg3->SetFillColorAlpha(0, 1.0);
  leg3->Draw();
  
  sprintf(fname, "%s/ROC_rho1_chirp.png", netdir);
  co_canvas->SetLogy();
  co_canvas->Update();
  co_canvas->SaveAs(fname);
  delete co_canvas;
  delete mgr,mg2r,mg3r;
  delete leg, leg2, leg3;
  delete co_gr1,co_gr2,co_gr3,co_gr4; 
   // exit(0);

}
