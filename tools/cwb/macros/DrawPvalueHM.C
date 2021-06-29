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

// this macro produces the PE HM plot
//
// How to use: root 'DrawPvalueHM("roc_config.txt")'
//
// roc_config.txt format
//
// signal_statistic_1.txt noise_statistic_1.txt roc#1
// signal_statistic_2.txt noise_statistic_2.txt roc#2
//
// the files *_statistic_*.txt are produced in the final output directory by the cwb_pereport command 
// - */distributions/pestat_ResidualEnergy_statistic.txt
//
// roc#* are the labels displayed in the output plot legend 
// 


#define nHM_MAX  10

//#define SHOW_RESIDUALS
#define DUMP_PVALUES

int readPvalue(TString fname, vector<float>& valpha, vector<float>& vpvalue, vector<float>& vre);
int readConfig(TString fname, vector<TString>& vfpvalue, vector<TString>& vname);
void PlotPvalueHM(int nroc, TString ofname="");

int gCOLOR[nHM_MAX];   

wavearray<double> wALPHA[nHM_MAX];
wavearray<double> wPVALUE[nHM_MAX];
wavearray<double> wRE[nHM_MAX];
int wSIZE[nHM_MAX];

TString gNAME[nHM_MAX];

void DrawPvalueHM(TString ifname, TString ofname="") {

  // init blind colors 
  gCOLOR[0] = CWB::Toolbox::getTableau10BlindColor("DeepSkyBlue4");
  gCOLOR[1] = CWB::Toolbox::getTableau10BlindColor("DarkOrange1");
  gCOLOR[2] = CWB::Toolbox::getTableau10BlindColor("DarkGray");
  gCOLOR[3] = CWB::Toolbox::getTableau10BlindColor("DimGray");
  gCOLOR[4] = CWB::Toolbox::getTableau10BlindColor("SkyBlue3");
  gCOLOR[5] = CWB::Toolbox::getTableau10BlindColor("Chocolate3");
  gCOLOR[6] = CWB::Toolbox::getTableau10BlindColor("Gray");
  gCOLOR[7] = CWB::Toolbox::getTableau10BlindColor("SlateGray1");
  gCOLOR[8] = CWB::Toolbox::getTableau10BlindColor("SandyBrown");
  gCOLOR[9] = CWB::Toolbox::getTableau10BlindColor("LightGray");

  vector<TString> vfpvalue;
  vector<TString> vname;
  int nHM = readConfig(ifname, vfpvalue, vname);

  for(int n=0;n<nHM;n++) {

    gNAME[n] = vname[n];
    gNAME[n].ReplaceAll("*"," ");

    vector<float> valpha,vpvalue,vre;
     
    int size = readPvalue(vfpvalue[n], valpha, vpvalue, vre);
    wSIZE[n]=size;

    wALPHA[n].resize(size);
    wPVALUE[n].resize(size);
    wRE[n].resize(size);
    int j=0; 
    for(int i=0;i<size;i++) {
#ifdef DUMP_PVALUES
      if(valpha[i]>=0.2 && valpha[i]<3.5) {
        cout << valpha[i] << "\t" << vpvalue[i] << endl;
      }
#endif
      wALPHA[n][j] = valpha[i];
      wPVALUE[n][j] = vpvalue[i];
      wRE[n][j] = vre[i];
      j++;
    }
    wSIZE[n]=j;
    wALPHA[n].resize(j);
    wPVALUE[n].resize(j);
    wRE[n].resize(j);
  }

  PlotPvalueHM(nHM, ofname);
}

void PlotPvalueHM(int nHM, TString ofname) {

  // create plots
  gStyle->SetFrameBorderMode(0);     // remove the red box around canvas
  gROOT->ForceStyle();               

  gStyle->SetMarkerColor(50);
  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleW(0.98);     
  gStyle->SetTitleH(0.08);     
  gStyle->SetTitleY(0.98);     
  gStyle->SetFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleFont(12,"D");
//  gStyle->SetTitleSize(0.18);

  TCanvas *canvas = new TCanvas("Higher Modes", "Higher Modes", 300, 40, 800, 500);
  canvas->Clear();                                                       
  canvas->ToggleEventStatus();                                           
//  canvas->SetLogx();                                                     
  canvas->SetLogy();                                                     
  canvas->SetGridx();                                                    
  canvas->SetGridy();                                                    
  canvas->SetFillColor(kWhite);                                          

  double xmin=2e20;
  double xmax=0;
  double ymin=2e20;
  double ymax=0;
  TGraph* gr[nHM_MAX];
  for(int n=0;n<nHM;n++) {
#ifdef SHOW_RESIDUALS
    gr[n] = new TGraph(wSIZE[n],wALPHA[n].data,wRE[n].data);
#else
    gr[n] = new TGraph(wSIZE[n],wALPHA[n].data,wPVALUE[n].data);
#endif
    for(int i=0;i<wSIZE[n];i++) {
      if(wALPHA[n][i]<xmin)  xmin=wALPHA[n][i];
      if(wALPHA[n][i]<xmax)  xmax=wALPHA[n][i];
      if(wPVALUE[n][i]<ymin) ymin=wPVALUE[n][i];
      if(wPVALUE[n][i]<ymax) ymax=wPVALUE[n][i];
    }
  }

  for(int n=0;n<nHM;n++) {                                            
    gr[n]->SetLineWidth(2);                                            
    gr[n]->SetMarkerColor(gCOLOR[n]);                                  
    gr[n]->SetMarkerStyle(20);                                  
    gr[n]->SetMarkerSize(0.8);                                  
    gr[n]->SetLineColor(gCOLOR[n]);                                    
    gr[n]->SetLineStyle(1);                                    
  }                                                                    

  TMultiGraph* mg = new TMultiGraph();
  char gTitle[256]; 
  sprintf(gTitle,"Higher Modes");
  mg->SetTitle(gTitle); 
  for(int n=0;n<nHM;n++) mg->Add(gr[n]);  
  //mg->Paint("APL");                        
  mg->Paint("AP");                        

  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);  
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);  
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5); 

  mg->GetHistogram()->GetXaxis()->SetRangeUser(0.001,7.5);
#ifndef SHOW_RESIDUALS
  mg->GetHistogram()->GetYaxis()->SetRangeUser(0.001,1);
#endif

  mg->GetHistogram()->GetXaxis()->SetNdivisions(509);
  mg->GetHistogram()->GetYaxis()->SetNdivisions(509);

  mg->GetXaxis()->SetTitle(gr[0]->GetXaxis()->GetTitle());
  mg->GetXaxis()->SetLabelFont(42);                       
  mg->GetYaxis()->SetLabelFont(42);                       
  mg->GetXaxis()->SetTitleFont(42);                       
  mg->GetYaxis()->SetTitleFont(42);                       
  mg->GetXaxis()->SetTitleOffset(0.80);                   
  mg->GetYaxis()->SetTitleOffset(0.90);                   
  mg->GetXaxis()->SetTitleSize(0.05);                     
  mg->GetYaxis()->SetTitleSize(0.05);                     
  mg->GetXaxis()->SetTitle("alpha");
#ifdef SHOW_RESIDUALS
  mg->GetYaxis()->SetTitle("residuals");
#else
  mg->GetYaxis()->SetTitle("p-value");
#endif
  mg->GetXaxis()->CenterTitle(true);
  mg->GetYaxis()->CenterTitle(true);

  mg->Draw("ALP");

  // draw the legend

  TLegend* leg;                
  double hleg = 0.18+nHM*0.06; 
  leg = new TLegend(0.593985,0.128692,0.8834586,hleg,NULL,"brNDC");

  leg->SetBorderSize(1);
  leg->SetTextAlign(22);
  leg->SetTextFont(12); 
  leg->SetLineColor(1); 
  leg->SetLineStyle(1); 
  leg->SetLineWidth(1); 
  leg->SetFillColor(0); 
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.05); 
  leg->SetLineColor(kBlack);
  leg->SetFillColor(kWhite);

  for(int n=0;n<nHM;n++) {
    char legLabel[256];
    sprintf(legLabel,"%s",gNAME[n].Data());
    leg->AddEntry(gr[n],legLabel,"lp");
  }
  leg->Draw();

  // save plot
  if(ofname!="") canvas->Print(ofname);

  return;
}

int readPvalue(TString fname, vector<float>& valpha, vector<float>& vpvalue, vector<float>& vre) {

  float a1,a2,pvalue,re;     

  ifstream in;
  in.open(fname.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fname.Data() << endl;exit(1);}

// S190814bv : ResidualEnergy 2.0572 pvalue 0.0537 ccut_wdm_fres 16 ccut_bchirp 2.2500 ccut_uchirp 2.5500 ccut_ltime 0.5000 ccut_rtime 0.0000

  char dy[32];

  while (1) {
    in >> dy >> dy >> dy >> re >> dy >> pvalue >> dy >> dy >> dy >> a1 >> dy >> a2 >> dy >> dy >> dy >> dy;
    if (!in.good()) break;
    valpha.push_back((a1+a2)/2);
    vpvalue.push_back(pvalue);
    vre.push_back(re);
  }
  in.close();

  return vpvalue.size();
}

int readConfig(TString fname, vector<TString>& vfpvalue, vector<TString>& vname) {

  char fpvalue[1024];
  char name[1024];

  ifstream in;
  in.open(fname.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fname.Data() << endl;exit(1);}

  while (1) {
    in >> fpvalue >> name;
    if (!in.good()) break;
    if(TString(fpvalue[0])=="#") continue;
    vfpvalue.push_back(fpvalue);
    vname.push_back(name);
  }
  in.close();

  return vname.size();
}

