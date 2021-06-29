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
// Draw injected vs detected events distributions vs SNRnet
// Note : this macro is used to generate the PRC report (cwb_report merge_label prc)
// Author : Gabriele Vedovato

namespace DrawRECvsINJ_namespace {   // used to avoid conflict with global variables

#include <vector>

void RECvsINJ(int mtype, TString gtype, TString ptitle, TString xtitle, TString ytitle, int nIFO, float T_win,
              int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar);

// --------------------------------------------------------
// Global variables
// --------------------------------------------------------
TCanvas*  gCANVAS = NULL;		// canvas object
TTree*    gTRWAVE = NULL;		// wave tree
TTree*    gTRMDC  = NULL;		// mdc tree


void 
DrawRECvsINJ(TString gtype, TString data_label, TString odir, TString merge_label, int nIFO, float T_win, 
             int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  if(gtype!="distance" && gtype!="snr") {
    cout << "DrawRECvsINJ : Error - wrong input gtype (snr/distance)" << endl;
    gSystem->Exit(1);
  }

  // ------------------------------------------------
  // create canvas
  // ------------------------------------------------
  gCANVAS = new TCanvas("fom", "PRC", 300,40, 600, 500);
  gCANVAS->Range(-19.4801,-9.25,-17.4775,83.25);
  gCANVAS->SetBorderSize(2);
  gCANVAS->SetFrameFillColor(0);
  gCANVAS->SetGridx();
  gCANVAS->SetGridy();

  // ------------------------------------------------
  // open wave/mdc trees
  // ------------------------------------------------
  // open wave file
  char sim_file_name[1024];
  sprintf(sim_file_name,"merge/wave_%s.%s.root",data_label.Data(),merge_label.Data());
  TFile* fwave = TFile::Open(sim_file_name);
  gTRWAVE = (TTree*) gROOT->FindObject("waveburst");
  // open mdc file
  char mdc_file_name[1024];
  sprintf(mdc_file_name,"merge/mdc_%s.%s.root",data_label.Data(),merge_label.Data());
  TFile *fmdc = TFile::Open(mdc_file_name);
  gTRMDC = (TTree*) gROOT->FindObject("mdc");

  TString ptitle;
  TString xtitle;
  TString ytitle;

  if(gtype=="snr") gCANVAS->SetLogx(true);
  gCANVAS->SetLogy(true);
  ptitle="Reconstructed events vs Injected events";
  gStyle->SetOptStat(0);
  //if(nIFO==2) xtitle = "sqrt(hrss[0]^2+hrss[1]^2)";
  //if(nIFO==3) xtitle = "sqrt(hrss[0]^2+hrss[1]^2+hrss[2]^2)";
  xtitle = gtype=="snr" ? "Injected SNRnet" : "distance (Mpc)";
  ytitle = "counts";

  RECvsINJ(0, gtype, ptitle, xtitle, ytitle, nIFO, T_win, pp_inetcc, T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  char ofname[1024];
  if(gtype=="snr") {
    sprintf(ofname,"%s/inj_vs_rec_snr.gif",odir.Data());
  } else {
    sprintf(ofname,"%s/inj_vs_rec_distance.gif",odir.Data());
  }
  gCANVAS->Print(ofname);
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

void RECvsINJ(int mtype, TString gtype, TString ptitle, TString xtitle, TString ytitle, int nIFO, float T_win,
              int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  char sel[256]="";
  char tmp[256]="";
  char cut[256]="";
  char title[256]="";

  // INJECTED
  sprintf(cut,"");
  if(gtype=="snr") {
    strcpy(tmp,"snr[0]*snr[0]");
    for(int i=1;i<nIFO;i++) {sprintf(sel,"%s+snr[%d]*snr[%d]",tmp,i,i);strcpy(tmp,sel);}
    strcpy(tmp,sel);sprintf(sel,"sqrt(%s)>>hist(10000,1,10000)",tmp);
  } else {
    double dmin = gTRMDC->GetMinimum("distance");
    double dmax = gTRMDC->GetMaximum("distance");
    sprintf(sel,"distance>>hist(30,%f,%f)",dmin,dmax);
  }
  gTRMDC->Draw(sel,cut,"");
  int nmdc = gTRMDC->GetSelectedRows();
  cout << "nmdc : " << nmdc << endl; 

  TH2F *htemp = (TH2F*)gPad->GetPrimitive("hist");
  htemp->GetXaxis()->SetTitle(xtitle);
  htemp->GetYaxis()->SetTitle(ytitle);
  //htemp->GetXaxis()->SetRangeUser(4e-25,1e-21);
  if(gtype=="snr") {
    htemp->GetYaxis()->SetRangeUser(0.1, pow(10.,TMath::Ceil(TMath::Log10(htemp->GetMaximum()))));
  } else {
    htemp->GetXaxis()->SetRangeUser(gTRMDC->GetMinimum("distance"),gTRMDC->GetMaximum("distance"));
    htemp->SetMinimum(0.9);
  }
  htemp->GetXaxis()->SetTitleOffset(1.35);
  htemp->GetYaxis()->SetTitleOffset(1.35);
  htemp->GetXaxis()->CenterTitle(true);
  htemp->GetYaxis()->CenterTitle(true);
  htemp->GetXaxis()->SetLabelFont(42);
  htemp->GetXaxis()->SetTitleFont(42);
  htemp->GetYaxis()->SetLabelFont(42);
  htemp->GetYaxis()->SetTitleFont(42);


  // DETECTED
  if(gtype=="snr") {
    strcpy(tmp,"iSNR[0]");
    for(int i=1;i<nIFO;i++) {sprintf(sel,"%s+iSNR[%d]",tmp,i);strcpy(tmp,sel);}
    strcpy(tmp,sel);sprintf(sel,"sqrt(%s)",tmp);
  } else {
    sprintf(sel,"range[1]");
  }
  sprintf(cut,"abs(time[0]-time[%d])<%g && netcc[%d]>%g && rho[%d]>%g", 
          nIFO,T_win,pp_inetcc,T_cor,pp_irho,T_cut);
  if(T_vED>0)  {strcpy(tmp,cut);sprintf(cut,"%s && neted[0]/ecor<%f",tmp,T_vED);}
  if(T_pen>0)  {strcpy(tmp,cut);sprintf(cut,"%s && penalty>%f",tmp,T_pen);}
  if(T_ifar>0) {strcpy(tmp,cut);sprintf(cut,"%s && ifar>(24.*3600.*365.)*%f",tmp,T_ifar);}
  gTRWAVE->SetLineColor(kRed);
  gTRWAVE->SetFillColor(kRed);
  gTRWAVE->Draw(sel,cut,"same");
  int nwave = gTRWAVE->GetSelectedRows();
  cout << "nwave : " << nwave << endl; 
  sprintf(title,"%s",cut);

  sprintf(title,"%s : inj = %d : rec : %d",ptitle.Data(),nmdc,nwave);
  htemp->SetTitle(title);

/*
  // print plot
  char ofname[256];
  if(mtype==0) {
    sprintf(ofname,"%s/%s_%s.png",
            gODIR.Data(),gPTYPE.Data(),"ALL");
  } else {
    sprintf(ofname,"%s/%s_%s.png",
            gODIR.Data(),gPTYPE.Data(),gMDC->wfList[mtype-1].name.Data());
  }
  cout << ofname << endl;
  gCANVAS->Print(ofname);
*/

  return;
}

} // end namespace

void 
DrawRECvsINJ(TString gtype, TString data_label, TString odir, TString merge_label, int nIFO, float T_win, 
             int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  DrawRECvsINJ_namespace::DrawRECvsINJ(gtype, data_label, odir, merge_label, nIFO, T_win, 
                                       pp_inetcc, T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar); 
  return;
}
