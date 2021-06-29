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


#define ALPHA 0.1

#define MAX_GW	100


void DrawMatchPE(TString ifname, TString odir, TString label, bool title=true) {

  // init blind colors 
  Color_t color[4];
  color[0] = CWB::Toolbox::getTableau10BlindColor("DarkOrange1");
  color[1] = CWB::Toolbox::getTableau10BlindColor("DeepSkyBlue4");
  color[2] = CWB::Toolbox::getTableau10BlindColor("DarkGray");
  color[3] = CWB::Toolbox::getTableau10BlindColor("SandyBrown");

  ifstream in;
  in.open(ifname.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << ifname.Data() << endl;exit(1);}

  int N=0;
  double match[MAX_GW],pvalue[MAX_GW],X[MAX_GW],zero[MAX_GW];
  double median[MAX_GW],l99[MAX_GW],l90[MAX_GW],l50[MAX_GW],u50[MAX_GW],u90[MAX_GW],u99[MAX_GW];
  char gw_name[MAX_GW][256];
  char dummy[256];
  while(1) {
    in >> gw_name[N] >> dummy >> dummy >> match[N] >> dummy >> pvalue[N] >> dummy >> median[N] 
       >> dummy >> l99[N] >> dummy >> l90[N] >> dummy >> l50[N] >> dummy >> u50[N] >> dummy >> u90[N] >> dummy >> u99[N];
    if(!in.good()) break;
    cout <<" "<< gw_name[N] <<" "<< match[N] <<" "<< pvalue[N] <<" "<< median[N] 
         <<" "<< l99[N] <<" "<< l90[N] <<" "<< l50[N] <<" "<< u50[N] <<" "<< u90[N] <<" "<< u99[N] << endl;;
    zero[N]=0.;

    l50[N] = fabs(median[N]-l50[N]);
    u50[N] = fabs(u50[N]-median[N]);

    l90[N] = (median[N]-l90[N]);
    u90[N] = (u90[N]-median[N]);

    l99[N] = fabs(median[N]-l99[N]);
    u99[N] = fabs(u99[N]-median[N]);

    X[N]=N+1;
    N++;
  }
  in.close();

  int *index = new int[N];
  TMath::Sort(N,match,index,false);

  double smatch[MAX_GW];
  double smedian[MAX_GW],sl99[MAX_GW],sl90[MAX_GW],sl50[MAX_GW],su50[MAX_GW],su90[MAX_GW],su99[MAX_GW];
  double xmax=0,ymax=0;
  TString sgw_name[MAX_GW];
  for(int n=0;n<N;n++) {
   X[n]/=N;
   smatch[n]=match[index[n]];
   smedian[n]=median[index[n]];
   sl50[n]=l50[index[n]];
   su50[n]=u50[index[n]];
   sl90[n]=l90[index[n]];
   su90[n]=u90[index[n]];
   sl99[n]=l99[index[n]];
   su99[n]=u99[index[n]];
   sgw_name[n]=gw_name[index[n]];
   if(smatch[n]>xmax) xmax=smatch[n];
   if(smedian[n]+su99[n]>ymax) ymax=smedian[n]+su99[n];
   cout << n << "\t" << sgw_name[n] << "\tmatch " << smatch[n] << "\tmedian " << smedian[n] << "\tlow90 " << sl90[n] << "\tup90 " << su90[n] << endl;
  }

  TCanvas* canvas = new TCanvas("fom", "fom", 300,40, 800, 800); 
  canvas->SetGrid();
  canvas->SetLeftMargin(0.15);
  canvas->SetBottomMargin(0.15);

  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleTextColor(kBlack);
  gStyle->SetTitleFont(12,"D");

  // draw frame
  TH2F *frame = new TH2F("frame","",1000,0.,1,1000,0.,1);
  if(title) frame->SetTitle(TString("Matching")+" ( "+label+" )"); else frame->SetTitle("");
  frame->SetStats(0);
  frame->GetXaxis()->SetTitle("OnSource Match");
  frame->GetYaxis()->SetTitle("OffSource Match");
  frame->GetXaxis()->SetTitleOffset(1.50);
  frame->GetYaxis()->SetTitleOffset(1.50);
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetYaxis()->CenterTitle(kTRUE);
  frame->GetXaxis()->SetTitleFont(132);
  frame->GetXaxis()->SetLabelFont(132);
  frame->GetYaxis()->SetTitleFont(132);
  frame->GetYaxis()->SetLabelFont(132);
  frame->GetXaxis()->SetTitleSize(0.04);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelOffset(0.01);
  frame->GetXaxis()->SetNdivisions(9);
  frame->GetYaxis()->SetNdivisions(9);
  frame->LabelsOption("x");
  frame->Draw();

  // draw match
  TGraphAsymmErrors* gmatch = new TGraphAsymmErrors(N,match,median,zero,zero,l90,u90);
  gmatch->SetMarkerStyle(20);
  gmatch->SetMarkerSize(1.0);
  gmatch->SetLineWidth(1);
  gmatch->SetLineColor(kBlack);
  gmatch->SetMarkerColor(color[0]);
  gmatch->Draw("PSAME");

  // Draw diagonal
  double xdiag[2]={0,1};
  double ydiag[2]={0,1};
  TGraph* gdiag = new TGraph(2,xdiag,ydiag);
  gmatch->SetLineColor(color[1]);
  gdiag->SetLineStyle(9);
  gdiag->Draw("SAME");

  // draw legend
  TLegend *leg = new TLegend(0.4047619,0.1614987,0.8897243,0.2635659,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextAlign(22);
  leg->SetTextFont(12);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.03);
  leg->SetLineColor(kBlack);
  leg->SetFillColor(kWhite);
  leg->AddEntry(gmatch,"median ( 90% confidence )","lp");
  leg->AddEntry(gdiag,"null hypothesis","lp");
  leg->Draw();

  if(odir!="") odir=odir+"/";
  TString ofname=odir+gSystem->BaseName(ifname);;
  ofname.ReplaceAll(".txt","_match.png"); 

  canvas->Print(ofname);
}

