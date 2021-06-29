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
// Draw the Normalized Residual Errors & Reaconstructed SNRnet vs Injected SNRnet
// Note : this macro is used to generate the PRC report (cwb_report merge_label prc)
// Author : Gabriele Vedovato

namespace DrawWRC_namespace {	// used to avoid conflict with global variables

#include <vector>

void Plot(TString gtype, TString ptitle, TString xtitle, TString ytitle, int nIFO, float T_win,
          int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar);

// --------------------------------------------------------
// Global variables
// --------------------------------------------------------
TCanvas*  gCANVAS = NULL;		// canvas object
TTree*    gTRWAVE = NULL;		// wave tree


void 
DrawWRC(TString gtype, TString data_label, TString odir, TString merge_label, int nIFO, float T_win, 
             int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  if(gtype!="nre" && gtype!="snr" && gtype!="ff" && gtype!="of" && gtype!="ofnre" && gtype!="offf" &&
     gtype!="mch1" && gtype!="mch2" && gtype!="mch3" && gtype!="mch4" && 
     gtype!="mch5" && gtype!="mch6" && gtype!="mch7" && gtype!="mch8" &&
     gtype!="mch9" && gtype!="mch10" && gtype!="mch11" && gtype!="mch12") {
    cout << "DrawWRC : Error - wrong input gtype (snr/nre/ff/of/ofnre/offf/mch1:12)" << endl;
    gSystem->Exit(1);
  }

  // ------------------------------------------------
  // create canvas
  // ------------------------------------------------
  gCANVAS = new TCanvas("fom", "PRC", 300,40, 600, 600);
  gCANVAS->Range(-19.4801,-9.25,-17.4775,83.25);
  if(gtype=="ofnre" || gtype=="offf") gCANVAS->SetRightMargin(0.13);
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

  TString ptitle;
  TString xtitle;
  TString ytitle;

  gCANVAS->SetLogx(true);
  if(gtype=="snr")  gCANVAS->SetLogy(true);
  if(gtype=="snr")  ptitle="Reconstructed SNRnet vs Injected SNRnet";
  if(gtype=="nre")  ptitle="Normalized Residual Energy (NRE)";
  if(gtype=="ff")   ptitle="Fitting Factor (FF)";
  if(gtype=="of")   ptitle="Overlap Factor (OF)";
  if(gtype=="ofnre") {
    ptitle="Overlap Factor vs Normalized Residual Energy (OFNRE)";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="offf") {
    ptitle="Overlap Factor vs Fitting Factor";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch1") {
    ptitle="Rec Chirp Mass vs Inj Chirp Mass (CM)";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch2") {
    ptitle="Rec-Inj Chirp Mass vs Injected SNRnet";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch3") {
    ptitle="Chirp Ellipticity vs Inj Chirp Mass";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch4") {
    ptitle="Chirp Ellipticity vs Injected SNRnet";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch5") {
    ptitle="Chirp Pixel Fraction vs Inj Chirp Mass";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch6") {
    ptitle="Chirp Pixel Fraction vs Injected SNRnet";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch7") {
    ptitle="Chirp Energy Fraction vs Inj Chirp Mass";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch8") {
    ptitle="Chirp Energy Fraction vs Injected SNRnet";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch9") {
    ptitle="Chirp Mass Error vs Inj Chirp Mass";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch10") {
    ptitle="Chirp Mass Error vs Injected SNRnet";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch11") {
    ptitle="PE Rec-Inj Chirp Mass vs Injected SNRnet";
    gCANVAS->SetLogx(false);
  }
  if(gtype=="mch12") {
    ptitle="PE Chirp Mass Error vs Injected SNRnet";
    gCANVAS->SetLogx(false);
  }
  gStyle->SetOptStat(0);
  //if(nIFO==2) xtitle = "sqrt(hrss[0]^2+hrss[1]^2)";
  //if(nIFO==3) xtitle = "sqrt(hrss[0]^2+hrss[1]^2+hrss[2]^2)";
  xtitle = "Injected SNRnet";
  if(gtype=="snr")  ytitle = "Reconstructed SNRnet";
  if(gtype=="nre")  ytitle = "NRE";
  if(gtype=="ff")   ytitle = "FF";
  if(gtype=="of")   ytitle = "OF";
  if(gtype=="ofnre") {
    xtitle = "1-NRE";
    ytitle = "OF";
  }
  if(gtype=="offf") {
    xtitle = "FF";
    ytitle = "OF";
  }
  if(gtype=="mch1") {
    xtitle = "Injected Chirp Mass";
    ytitle = "Reconstructed Chirp Mass";
  }
  if(gtype=="mch2") ytitle = "Rec-Inj Chirp Mass";
  if(gtype=="mch3") {
    xtitle = "Injected Chirp Mass";
    ytitle = "Chirp Ellipticity";
  }
  if(gtype=="mch4") ytitle = "Chirp Ellipticity";
  if(gtype=="mch5") {
    xtitle = "Injected Chirp Mass";
    ytitle = "Chirp Pixel Fraction";
  }
  if(gtype=="mch6") ytitle = "Chirp Pixel Fraction";
  if(gtype=="mch7") {
    xtitle = "Injected Chirp Mass";
    ytitle = "Chirp Energy Fraction";
  }
  if(gtype=="mch8") ytitle = "Chirp Energy Fraction";
  if(gtype=="mch9") {
    xtitle = "Injected Chirp Mass";
    ytitle = "Chirp Mass Error";
  }
  if(gtype=="mch10") ytitle = "Chirp Mass Error";
  if(gtype=="mch11") ytitle = "PE Rec-Inj Chirp Mass";
  if(gtype=="mch12") ytitle = "PE Chirp Mass Error";
  Plot(gtype, ptitle, xtitle, ytitle, nIFO, T_win, pp_inetcc, T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);

  char ofname[1024];
  if(gtype=="snr")  sprintf(ofname,"%s/osnr_vs_isnr.gif",odir.Data());
  if(gtype=="nre")  sprintf(ofname,"%s/nre_vs_isnr.gif",odir.Data());
  if(gtype=="ff")   sprintf(ofname,"%s/ff_vs_isnr.gif",odir.Data());
  if(gtype=="of")   sprintf(ofname,"%s/of_vs_isnr.gif",odir.Data());
  if(gtype=="ofnre") sprintf(ofname,"%s/of_vs_nre.gif",odir.Data());
  if(gtype=="offf") sprintf(ofname,"%s/of_vs_ff.gif",odir.Data());
  if(gtype=="mch1") sprintf(ofname,"%s/omch_vs_imch.gif",odir.Data());
  if(gtype=="mch2") sprintf(ofname,"%s/dmch_vs_isnr.gif",odir.Data());
  if(gtype=="mch3") sprintf(ofname,"%s/elch_vs_imch.gif",odir.Data());
  if(gtype=="mch4") sprintf(ofname,"%s/elch_vs_isnr.gif",odir.Data());
  if(gtype=="mch5") sprintf(ofname,"%s/pfch_vs_imch.gif",odir.Data());
  if(gtype=="mch6") sprintf(ofname,"%s/pfch_vs_isnr.gif",odir.Data());
  if(gtype=="mch7") sprintf(ofname,"%s/efch_vs_imch.gif",odir.Data());
  if(gtype=="mch8") sprintf(ofname,"%s/efch_vs_isnr.gif",odir.Data());
  if(gtype=="mch9") sprintf(ofname,"%s/emch_vs_imch.gif",odir.Data());
  if(gtype=="mch10") sprintf(ofname,"%s/emch_vs_isnr.gif",odir.Data());
  if(gtype=="mch11") sprintf(ofname,"%s/opemch_vs_isnr.gif",odir.Data());
  if(gtype=="mch12") sprintf(ofname,"%s/epemch_vs_isnr.gif",odir.Data());
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

void Plot(TString gtype, TString ptitle, TString xtitle, TString ytitle, int nIFO, float T_win,
          int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  if(gtype=="mch11" || gtype=="mch12") {  // check if pe_mch leaf exist
    TBranch* branch; 
    bool pe_mch_exists=false;
    TIter next(gTRWAVE->GetListOfBranches());
    while ((branch=(TBranch*)next())) {
      if(TString("pe_mch").CompareTo(branch->GetName())==0) pe_mch_exists=true;
    }
    next.Reset();
    if(!pe_mch_exists) return;
  }

  char sel[256]="";
  char tmp[256]="";
  char cut[256]="";
  char title[256]="";
  char num[256]="";
  char den[256]="";
  char ffnum[256]="";
  char iffden[256]="";
  char offden[256]="";

  // DETECTED
  sprintf(cut,"");
  // fitting factor
  strcpy(tmp,"ioSNR[0]");
  for(int i=1;i<nIFO;i++) {sprintf(ffnum,"%s+ioSNR[%d]",tmp,i);strcpy(tmp,ffnum);}
  strcpy(tmp,"iSNR[0]");
  for(int i=1;i<nIFO;i++) {sprintf(iffden,"%s+iSNR[%d]",tmp,i);strcpy(tmp,iffden);}
  strcpy(tmp,"oSNR[0]");
  for(int i=1;i<nIFO;i++) {sprintf(offden,"%s+oSNR[%d]",tmp,i);strcpy(tmp,offden);}
  // residual energy
  strcpy(tmp,"(iSNR[0]+oSNR[0]-2*ioSNR[0])");
  for(int i=1;i<nIFO;i++) {sprintf(num,"%s+(iSNR[%d]+oSNR[%d]-2*ioSNR[%d])",tmp,i,i,i);strcpy(tmp,num);}
  strcpy(tmp,"iSNR[0]");
  for(int i=1;i<nIFO;i++) {sprintf(den,"%s+iSNR[%d]",tmp,i);strcpy(tmp,den);}
  if(gtype=="snr")  sprintf(sel,"sqrt(likelihood):sqrt(%s)>>htemp",den);
  if(gtype=="nre")  sprintf(sel,"(%s)/(%s):sqrt(%s)>>htemp",num,den,den);
  if(gtype=="ff")   sprintf(sel,"(%s)/sqrt((%s)*(%s)):sqrt(%s)>>htemp",ffnum,iffden,iffden,den);
  if(gtype=="of")   sprintf(sel,"(%s)/sqrt((%s)*(%s)):sqrt(%s)>>htemp",ffnum,iffden,offden,den);
  if(gtype=="ofnre") sprintf(sel,"(%s)/sqrt((%s)*(%s)):1-(%s)/(%s)>>htemp",ffnum,iffden,offden,num,den);
  if(gtype=="offf") sprintf(sel,"(%s)/sqrt((%s)*(%s)):(%s)/sqrt((%s)*(%s))>>htemp",ffnum,iffden,offden,ffnum,iffden,iffden);
  if(gtype=="mch1") sprintf(sel,"chirp[1]:chirp[0]>>htemp");
  if(gtype=="mch2") sprintf(sel,"chirp[1]-chirp[0]:sqrt(%s)>>htemp",den);
  if(gtype=="mch3") sprintf(sel,"chirp[3]:chirp[0]>>htemp");
  if(gtype=="mch4") sprintf(sel,"chirp[3]:sqrt(%s)>>htemp",den);
  if(gtype=="mch5") sprintf(sel,"chirp[4]:chirp[0]>>htemp");
  if(gtype=="mch6") sprintf(sel,"chirp[4]:sqrt(%s)>>htemp",den);
  if(gtype=="mch7") sprintf(sel,"chirp[5]:chirp[0]>>htemp");
  if(gtype=="mch8") sprintf(sel,"chirp[5]:sqrt(%s)>>htemp",den);
  if(gtype=="mch9") sprintf(sel,"chirp[2]:chirp[0]>>htemp");
  if(gtype=="mch10") sprintf(sel,"chirp[2]:sqrt(%s)>>htemp",den);
  if(gtype=="mch11") sprintf(sel,"pe_mch[0]-chirp[0]:sqrt(%s)>>htemp",den);
  if(gtype=="mch12") sprintf(sel,"pe_mch[1]:sqrt(%s)>>htemp",den);
  sprintf(cut,"abs(time[0]-time[%d])<%g && netcc[%d]>%g && rho[%d]>%g", 
          nIFO,T_win,pp_inetcc,T_cor,pp_irho,T_cut);
  if(T_vED>0)  {strcpy(tmp,cut);sprintf(cut,"%s && neted[0]/ecor<%f",tmp,T_vED);}
  if(T_pen>0)  {strcpy(tmp,cut);sprintf(cut,"%s && penalty>%f",tmp,T_pen);}
  if(T_ifar>0) {strcpy(tmp,cut);sprintf(cut,"%s && ifar>(24.*3600.*365.)*%f",tmp,T_ifar);}
  gTRWAVE->SetMarkerColor(kRed);
  gTRWAVE->SetMarkerStyle(20);
  gTRWAVE->SetMarkerSize(0.5);
  gTRWAVE->SetLineColor(kRed);
  gTRWAVE->SetFillColor(kRed);
  if(gtype=="ofnre" || gtype=="offf") gTRWAVE->Draw(sel,cut,"colz");
  else                                gTRWAVE->Draw(sel,cut);
  int nwave = gTRWAVE->GetSelectedRows();
  cout << "nwave : " << nwave << endl; 
  sprintf(title,"%s",cut);

  if(gtype=="nre") {
    double* inre = gTRWAVE->GetV1(); 
    for(int i=0;i<nwave;i++) if(inre[i]>1) inre[i]=1; 
    for(int i=0;i<nwave;i++) if(inre[i]<0.001) inre[i]=0.001; 
  }

  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  sprintf(title,"%s : rec : %d",ptitle.Data(),nwave);
  htemp->SetTitle(title);
  htemp->GetXaxis()->SetTitle(xtitle);
  htemp->GetYaxis()->SetTitle(ytitle);
  htemp->GetXaxis()->SetMoreLogLabels();
  if(gtype=="snr") htemp->GetYaxis()->SetMoreLogLabels();
  if(gtype=="nre") htemp->GetYaxis()->SetRangeUser(0.0,1.0);
  if(gtype=="ff")  htemp->GetYaxis()->SetRangeUser(0.0,1.0);
  if(gtype=="of")  htemp->GetYaxis()->SetRangeUser(0.0,1.0);
  if(gtype=="ofnre") {htemp->GetXaxis()->SetRangeUser(-4.0,1.0);htemp->GetYaxis()->SetRangeUser(-1.0,1.0);}
  if(gtype=="offf") {htemp->GetXaxis()->SetRangeUser(-2.0,2.0);htemp->GetYaxis()->SetRangeUser(-1.0,1.0);}
  htemp->GetXaxis()->SetTitleOffset(1.35);
  htemp->GetYaxis()->SetTitleOffset(1.50);
  htemp->GetXaxis()->CenterTitle(true);
  htemp->GetYaxis()->CenterTitle(true);
  htemp->GetXaxis()->SetLabelFont(42);
  htemp->GetXaxis()->SetTitleFont(42);
  htemp->GetYaxis()->SetLabelFont(42);
  htemp->GetYaxis()->SetTitleFont(42);

  // draw reference line
  if(gtype=="snr" || gtype=="mch1") {
    double lmin=htemp->GetYaxis()->GetXmin();
    double lmax=htemp->GetYaxis()->GetXmax();
    if(htemp->GetXaxis()->GetXmin()>lmin) lmin=htemp->GetXaxis()->GetXmin();
    if(htemp->GetXaxis()->GetXmax()<lmax) lmax=htemp->GetXaxis()->GetXmax();
    TLine *line = new TLine(lmin,lmin,lmax,lmax);
    line->SetLineColor(kBlue);
    line->Draw();
  } else {
  }

  return;
}

} // end namespace

void 
DrawWRC(TString gtype, TString data_label, TString odir, TString merge_label, int nIFO, float T_win, 
             int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  DrawWRC_namespace::DrawWRC(gtype, data_label, odir, merge_label, nIFO, T_win, 
                             pp_inetcc, T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar); 
  return;
}
