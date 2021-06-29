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

// This macro is used to find the correspondence between local and global coverage

// N=20 -> 98.950% global = 90.0089% local
// N=23 -> 99.032% global = 90.0013% local

#define NTRIALS		100000		// trials

#define NPVALUES	20		// number of pvalues
#define COVERAGE 	98.950

//#define NPVALUES	23		// number of pvalues
//#define COVERAGE 	99.032

//#define NPVALUES	100		// number of pvalues
//#define COVERAGE 	99.5

#include "TGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/SpecFuncMathCore.h"
#include <iostream>

double vpvalue[NTRIALS][NPVALUES];

double InverseProbabilityBH(int N, int k, double p);
TGraph* DrawProbabilityBH(int N, double p, Color_t color, int style=9, TString mode="ALP");

void PvalueLocal2GlobalBelts(double coverage=COVERAGE, bool draw_pvalues_plot=true) {
//
// input coverage is the local coverage
// NOTE: NPVALUES must be setted 
// this macro printout the global coverage
//

  if(coverage<0 || coverage>100) {cout << "Error: wrong input coverage, must be [0,1]" << endl;exit(1);}

  for(int i=0;i<NTRIALS;i++) {						// loop over NTRIALS 
    double pvalue[NPVALUES];
    int index[NPVALUES];
    for(int j=0;j<NPVALUES;j++) pvalue[j] = gRandom->Uniform(0,1);	// fill pvalue array with NPVALUES uniform random numbers [0,1] 	
    TMath::Sort(NPVALUES,pvalue,index,kFALSE);				// sort NPVALUES with increasing order
    for(int j=0;j<NPVALUES;j++) vpvalue[i][j]=pvalue[index[j]];		// save sorted pvalues into vpvalue
  }

  // compute lower/upper belts for coverage

  double median_belt[NPVALUES];						// median belt array
  double lower_belt[NPVALUES];						// lower belt array
  double upper_belt[NPVALUES];						// upper belt array

  double median      = (50.)/100.;					// median
  double lower_tail  = ((100.-coverage)/2.)/100.;			// lower tail
  double upper_tail  = 1.-lower_tail;					// upper tail

  int imedian = int(NTRIALS*median);        if(imedian>=NTRIALS) imedian=NTRIALS-1;
  int ilower  = int(NTRIALS*lower_tail);    if(ilower>=NTRIALS) ilower=NTRIALS-1;
  int iupper  = int(NTRIALS*upper_tail);    if(iupper>=NTRIALS) iupper=NTRIALS-1;

  for(int j=0;j<NPVALUES;j++) {
    int index[NTRIALS];
    double pvalue[NTRIALS];
    for(int i=0;i<NTRIALS;i++) pvalue[i]=vpvalue[i][j];
    TMath::Sort(NTRIALS,pvalue,index,false);
    median_belt[j] = pvalue[index[imedian]];
    lower_belt[j]  = pvalue[index[ilower]];
    upper_belt[j]  = pvalue[index[iupper]];
  }

  // compute the fraction of trials inside coverage

  int nInside=0;
  for(int i=0;i<NTRIALS;i++) {
    bool iSinside=true;
    for(int j=0;j<NPVALUES;j++) if(vpvalue[i][j]<lower_belt[j] || vpvalue[i][j]>upper_belt[j]) iSinside=false;
    if(iSinside) nInside++;
  }
  double global_coverage = 100.*nInside/double(NTRIALS);

  cout << endl;
  cout << "percentage of trials inside coverage " << coverage << "% = " << global_coverage << endl;
  cout << endl;

  if(!draw_pvalues_plot) exit(0);

  // draw numerical median_belt, lower_belt, upper_belt

  gStyle->SetFrameBorderMode(0);     // remove the red box around canvas
  gROOT->ForceStyle();
  gStyle->SetTitleFont(12,"D");
  TCanvas* canvas = new TCanvas("pvalue", "pvalue", 300,40, 1000, 600);
  canvas->Clear();
  canvas->SetLogx();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);

  double x[NPVALUES]; for(int i=0;i<NPVALUES;i++) x[i]=i+1;

  TGraph* gmedian_belt = new TGraph(NPVALUES,x,median_belt);
  TGraph* glower_belt  = new TGraph(NPVALUES,x,lower_belt);
  TGraph* gupper_belt  = new TGraph(NPVALUES,x,upper_belt);

  glower_belt->Draw("ALP");
  gmedian_belt->Draw("SAME");
  gupper_belt->Draw("SAME");

  // draw analytic belts (NOTE: numerical and analytic belts are compatible !!!)

  TGraph* gmedian_belt2 = DrawProbabilityBH(NPVALUES, median, kRed, 10,"same");
  TGraph* glower_belt2  = DrawProbabilityBH(NPVALUES, lower_tail, kBlue,10,"same");
  TGraph* gupper_belt2  = DrawProbabilityBH(NPVALUES, upper_tail, kBlue,10,"same");

  // draw analytic local 90% coverage belts 

  double lower_tail3  = ((100.-90.)/2.)/100.;	// lower tail
  double upper_tail3  = 1.-lower_tail3;		// upper tail

  TGraph* glower_belt3  = DrawProbabilityBH(NPVALUES, lower_tail3, kGreen,10,"same");
  TGraph* gupper_belt3  = DrawProbabilityBH(NPVALUES, upper_tail3, kGreen,10,"same");

  // draw labels
  glower_belt->SetTitle(TString::Format("p-values distribution (N = %d)",NPVALUES));
  glower_belt->GetXaxis()->SetTitle("Cumulative number of events");
  glower_belt->GetYaxis()->SetTitle("p-value");

  glower_belt->GetYaxis()->SetRangeUser(0,1);
  glower_belt->GetXaxis()->SetTitleOffset(0.80);
  glower_belt->GetYaxis()->SetTitleOffset(0.80);
  glower_belt->GetXaxis()->CenterTitle(kTRUE);
  glower_belt->GetYaxis()->CenterTitle(kTRUE);
  glower_belt->GetXaxis()->SetTitleFont(132);
  glower_belt->GetXaxis()->SetLabelFont(132);
  glower_belt->GetYaxis()->SetTitleFont(132);
  glower_belt->GetYaxis()->SetLabelFont(132);
  glower_belt->GetXaxis()->SetTitleSize(0.045);
  glower_belt->GetXaxis()->SetLabelSize(0.045);
  glower_belt->GetYaxis()->SetTitleSize(0.045);
  glower_belt->GetYaxis()->SetLabelSize(0.045);
  glower_belt->GetYaxis()->SetLabelOffset(0.01);
  glower_belt->GetYaxis()->SetNdivisions(3);

  // draw legend
  TLegend *leg = new TLegend(0.4148297,0.1184669,0.8897796,0.3501742,NULL,"brNDC");
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
  leg->AddEntry(gmedian_belt,"numerical median","lp");
  leg->AddEntry(glower_belt,TString::Format("numerical %.4g%% confidence region",coverage));
  leg->AddEntry(gmedian_belt2,"analytical median","lp");
  leg->AddEntry(glower_belt2,TString::Format("analytical %.4g%% confidence region",coverage));
  leg->AddEntry(glower_belt3,TString::Format("analytical %.4g%% confidence region",90.));
  leg->Draw();

  return;
}

double InverseProbabilityBH(int N, int k, double p) {

#define LMAX            100
#define PRECISION       1.e-6

  int L=100;

  double x;
  double xmin=0;
  double xmax=1;

  double a=0.;
  for(int i=0;i<LMAX;i++) {
    x = gRandom->Uniform(xmin,xmax);
    //cout << i << "\txmin " << xmin << "\txmax " << xmax << "\tx " << x << endl; 
    a=ROOT::Math::inc_beta(x, k, N-k+1);
    if(a<p) xmin=x; else xmax=x;
    if(fabs(a-p)<PRECISION) break;
  }
  return x;
}

TGraph* DrawProbabilityBH(int N, double p, Color_t color, int style, TString mode) {
// BENJAMINI-HOCHBERG

  double* x = new double[N];
  for(int i=0;i<N;i++) x[i]=i+1;

  double* prob = new double[N];
  for(int k=1;k<=N;k++) prob[k-1]=InverseProbabilityBH(N, k, p);

  TGraph* gprob = new TGraph(N,x,prob);
  gprob->SetLineColor(color);
  gprob->SetLineWidth(3);
  gprob->SetLineStyle(style);
  gprob->Draw(mode);
  gprob->GetXaxis()->SetRangeUser(1,N);
  gprob->GetYaxis()->SetRangeUser(0,1);

  return gprob;
}

