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
// Draw Error Regions Percentage vs Coverage 
// Note : this macro is used to generate the PRC report (cwb_report merge_label prc)
// Author : Gabriele Vedovato


#define NMDC 1

//#define DUMP_ASCII

double GetPercentage(TString gtype, int iderA, TString fname, int nIFO, float T_win, int pp_inetcc,
                     float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar);

void 
DrawCoverageVsPercentagePRC(TString gtype, TString data_label, TString odir, TString merge_label, int nIFO, float T_win, 
                            int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  if(gtype!="nre" && gtype!="prc" && gtype!="fre") {
    cout << "DrawCoverageVsPercentagePRC : Error - wrong input gtype (prc/nre/fre)" << endl;
    gSystem->Exit(1);
  }


  char wave_file_name[1024];
  sprintf(wave_file_name,"merge/wave_%s.%s.root",data_label.Data(),merge_label.Data());

  TObjArray* token = TString(wave_file_name).Tokenize(TString('/'));
  TObjString* sfile = (TObjString*)token->At(token->GetEntries()-1);
  TString ofile = sfile->GetString();
  ofile.ReplaceAll(".root",".txt");
  //TString TITLE = sfile->GetString();
  //TITLE.ReplaceAll(".root","");
  TString TITLE = "";
  if(gtype=="prc") TITLE = "pp-plot-prc";
  if(gtype=="nre") TITLE = "pp-plot-nre";
  if(gtype=="fre") TITLE = "pp-plot-fre";

#ifdef DUMP_ASCII
  char ofname[1024];
  if(gtype=="prc") sprintf(ofname,"%s/pp_plot_prc.txt",odir.Data());
  if(gtype=="nre") sprintf(ofname,"%s/pp_plot_nre.txt",odir.Data());
  if(gtype=="fre") sprintf(ofname,"%s/pp_plot_fre.txt",odir.Data());
  ofstream out;
  out.open(ofname,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << ofname << endl;exit(1);}
  cout << "Create file : " << ofname << endl;
#endif

  // if coverage > percentage then erA is too large, the true erA is smaller (erA is conservative) 
  // if coverage < percentage then erA is too small, the true erA is greater (erA is optimistic) 

  double coverage[NMDC][11];
  for (int i=0;i<NMDC;i++) {
    coverage[i][0]=0.;
    coverage[i][10]=100.;
    for (int j=1;j<10;j++) {
      coverage[i][j]=GetPercentage(gtype, j,wave_file_name, nIFO, T_win, pp_inetcc,
                                   T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);
      if(coverage[i][j]<0) return;	// tree do not contains the branch (for example the nre erR array)
#ifdef DUMP_ASCII
      out << j*10 << " " << coverage[i][j] << endl;
#endif
      cout << j*10 << " " << coverage[i][j] << endl;
    }
  }
#ifdef DUMP_ASCII
  out.close();
#endif

  TCanvas* canvas = new TCanvas("fom", "PRC", 300,40, 600, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetFillColor(0);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLogy(false);
  canvas->SetLogx(false);

  gStyle->SetTitleH(0.062);
  gStyle->SetTitleW(0.98);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(72);
  gStyle->SetLineColor(kWhite);
  gStyle->SetPalette(1,0);
  gStyle->SetNumberContours(256);

  Style_t markers[32]= {20,21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,
                        21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,20 };

  Color_t colors[32] = { 6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7,
                         6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7 };

  double perc[11]={0,10,20,30,40,50,60,70,80,90,100};
  TGraph* gr[NMDC+11];
  for(int i=0;i<NMDC;i++) {
    gr[i] = new TGraph(11,perc,coverage[i]);
    gr[i]->SetLineColor(colors[i]);
    gr[i]->SetLineWidth(1);
    gr[i]->SetMarkerSize(1);
    gr[i]->SetMarkerColor(colors[i]);
    gr[i]->SetMarkerStyle(markers[i]);
  }
  double x[2]={0,100};
  double y[2]={0,100};
  gr[NMDC+1] = new TGraph(2,x,y);
  gr[NMDC+1]->SetLineColor(1);
  gr[NMDC+1]->SetLineWidth(2);
  gr[NMDC+1]->SetLineStyle(2);

  TMultiGraph* mg = new TMultiGraph();
  for(int i=0;i<NMDC;i++) mg->Add(gr[i],"lp");
  mg->Add(gr[NMDC+1],"lp");
  mg->Paint("ap");
  char title[256];
  sprintf(title,"%s",TITLE.Data());
  mg->GetHistogram()->SetTitle(title);
  mg->GetHistogram()->GetXaxis()->SetTitle("Percentage");
  mg->GetHistogram()->GetYaxis()->SetTitle("Coverage");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(0,100);
  mg->GetHistogram()->GetYaxis()->SetRangeUser(0,100);
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetXaxis()->CenterTitle(true);
  mg->GetHistogram()->GetYaxis()->CenterTitle(true);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetXaxis()->SetTitleFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetTitleFont(42);
  mg->GetHistogram()->GetXaxis()->SetNdivisions(511);
  mg->Draw("a");

  TLegend* leg;
  double hleg = 0.8-NMDC*0.05;
  leg = new TLegend(0.1291513,hleg,0.6244966,0.8738739,NULL,"brNDC");

  leg->SetBorderSize(1);
  leg->SetTextAlign(22);
  leg->SetTextFont(12);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.04);
  leg->SetLineColor(kBlack);
  leg->SetFillColor(kWhite);

  char label[256];
  for(int i=0;i<NMDC;i++) {
    sprintf(label,"%s ",data_label.Data());
    leg->AddEntry(gr[i],label,"lp");
  }
//  leg->Draw();

  char ofname[1024];
  if(gtype=="prc") sprintf(ofname,"%s/pp_plot_prc.gif",odir.Data());
  if(gtype=="nre") sprintf(ofname,"%s/pp_plot_nre.gif",odir.Data());
  if(gtype=="fre") sprintf(ofname,"%s/pp_plot_fre.gif",odir.Data());
  TString gfileName=ofname;
  canvas->Print(gfileName);
  TString pfileName=gfileName;
  pfileName.ReplaceAll(".gif",".png");
  char cmd[1024];
  sprintf(cmd,"convert %s %s",gfileName.Data(),pfileName.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"rm %s",gfileName.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);
  //exit(0);
}

double GetPercentage(TString gtype, int iderA, TString fname, int nIFO, float T_win, int pp_inetcc,
                     float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  TFile *ifile = TFile::Open(fname.Data());
  TTree* itree = (TTree *) gROOT->FindObject("waveburst");
  itree->SetEstimate(itree->GetEntries());

  // check if erR exist
  if(gtype=="nre") {
    TString nre_name="erR";
    TBranch* branch;
    bool nre_exist=false;
    TIter next(itree->GetListOfBranches());
    while ((branch=(TBranch*)next())) {
      if (nre_name.CompareTo(branch->GetName())==0) nre_exist=true;
    }
    if(!nre_exist) return -1;
  }

  // check if erF exist
  if(gtype=="fre") {
    TString fre_name="erF";
    TBranch* branch;
    bool fre_exist=false;
    TIter next(itree->GetListOfBranches());
    while ((branch=(TBranch*)next())) {
      if (fre_name.CompareTo(branch->GetName())==0) fre_exist=true;
    }
    if(!fre_exist) return -1;
  }

  char selection[1024];
  char cut[1024];
  char tmp[1024];
  sprintf(cut,"abs(time[0]-time[%d])<%f && netcc[%d]>%f && rho[%d]>%f",
          nIFO,T_win,pp_inetcc,T_cor,pp_irho,T_cut);
  if(T_vED>0)  {strcpy(tmp,cut);sprintf(cut,"%s && neted[0]/ecor<%f",tmp,T_vED);}
  if(T_pen>0)  {strcpy(tmp,cut);sprintf(cut,"%s && penalty>%f",tmp,T_pen);}
  if(T_ifar>0) {strcpy(tmp,cut);sprintf(cut,"%s && ifar>(24.*3600.*365.)*%f",tmp,T_ifar);}
  //cout << cut << endl;

  if(gtype=="prc") sprintf(selection,"erA[0]-erA[%d]>>hist_cumulative_erA1(2000,-100,100)",iderA);
  if(gtype=="nre") sprintf(selection,"erR[0]-erR[%d]>>hist_cumulative_erA1(2000,-1,1)",iderA);
  if(gtype=="fre") sprintf(selection,"erF[0]-erF[%d]>>hist_cumulative_erA1(2000,-1,1)",iderA);

  itree->Draw(selection,cut,"goff");
  TH1F *hist_cumulative_erA1 = (TH1F*)gDirectory->Get("hist_cumulative_erA1");
  int size = itree->GetSelectedRows();
  //cout << "size : " << size << endl;
  if(size==0) {
    delete hist_cumulative_erA1;
    delete itree;
    delete ifile;
    return 0;
  }

  double integral = hist_cumulative_erA1->ComputeIntegral();
  if (integral==0) {cout << "Empty histogram !!!" << endl;exit(0);}
  double* cumulative = hist_cumulative_erA1->GetIntegral();
  for (int i=0;i<hist_cumulative_erA1->GetNbinsX();i++) hist_cumulative_erA1->SetBinContent(i,cumulative[i]/integral);
  //cout << "integral " << integral << endl;


  double perc = 100.*hist_cumulative_erA1->GetBinContent(1001);
  delete hist_cumulative_erA1;
  delete itree;
  delete ifile;
  //ifile->Close();

  return perc;
}

