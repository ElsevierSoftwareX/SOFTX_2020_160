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

// this macro produces the PE ROC plot
//
// How to use: root 'DrawROCPE("roc_config.txt")'
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


#define nROC_MAX  10
#define wSIZE	  100

int readStatistic(TString fname, vector<float>& vstat) ;
int readConfig(TString fname, vector<TString>& vfsignal, vector<TString>& vfnoise, vector<TString>& vname);
void PlotROC(int nroc, TString ofname="");

int gCOLOR[nROC_MAX];   

wavearray<double> wSIGNAL[nROC_MAX];
wavearray<double> wNOISE[nROC_MAX];

TString gNAME[nROC_MAX];

void DrawROCxPE(TString ifname, TString ofname="") {

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

  vector<TString> vfsignal;
  vector<TString> vfnoise; 
  vector<TString> vname;
  int nROC = readConfig(ifname, vfsignal, vfnoise, vname);

  for(int n=0;n<nROC;n++) {

    gNAME[n] = vname[n];
    gNAME[n].ReplaceAll("*"," ");

    vector<float> vsignal;
    vector<float> vnoise;
     
    int nsignal = readStatistic(vfsignal[n], vsignal);
    int nnoise  = readStatistic(vfnoise[n], vnoise);

    wSIGNAL[n].resize(wSIZE);
    wNOISE[n].resize(wSIZE);
 
    wSIGNAL[n]=1.0;
 
    float dthr = vnoise[vnoise.size()-1]/wSIZE;

    for(int i=0;i<wSIZE;i++) {

      float thr = i*dthr;

      int isignal=0;
      for(int j=vsignal.size();j>=0;j--) if(vsignal[j]>thr) isignal++;
      int inoise=0;
      for(int j=vnoise.size();j>=0;j--)  if(vnoise[j]>thr)  inoise++;

      wSIGNAL[n][i] = (float)isignal/(float)vsignal.size();
      wNOISE[n][i]  = (float)inoise/(float)vnoise.size();
    }
  }

  PlotROC(nROC, ofname);
}

void PlotROC(int nROC, TString ofname) {

  // create plots
  gStyle->SetFrameBorderMode(0);     // remove the red box around canvas
  gROOT->ForceStyle();               

  gStyle->SetMarkerColor(50);
  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleW(0.98);     
  gStyle->SetTitleH(0.05);     
  gStyle->SetTitleY(0.98);     
  gStyle->SetFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleFont(12,"D");
//  gStyle->SetTitleSize(0.18);

  TCanvas *canvas = new TCanvas("ROC", "ROC", 300, 40, 800, 500);
  canvas->Clear();                                                       
  canvas->ToggleEventStatus();                                           
  canvas->SetLogx();                                                     
  canvas->SetLogy();                                                     
  canvas->SetGridx();                                                    
  canvas->SetGridy();                                                    
  canvas->SetFillColor(kWhite);                                          

  double xmin=2e20;
  double xmax=0;
  double ymin=2e20;
  double ymax=0;
  TGraph* gr[nROC_MAX];
  for(int n=0;n<nROC;n++) {
    gr[n] = new TGraph(wSIZE,wNOISE[n].data,wSIGNAL[n].data);
    for(int i=0;i<wSIZE;i++) {
      if(wNOISE[n][i]<xmin)  xmin=wNOISE[n][i];
      if(wNOISE[n][i]<xmax)  xmax=wNOISE[n][i];
      if(wSIGNAL[n][i]<ymin) ymin=wSIGNAL[n][i];
      if(wSIGNAL[n][i]<ymax) ymax=wSIGNAL[n][i];
    }
  }

  for(int n=0;n<nROC;n++) {                                            
    gr[n]->SetLineWidth(2);                                            
    gr[n]->SetMarkerColor(gCOLOR[n]);                                  
    gr[n]->SetMarkerStyle(20);                                  
    gr[n]->SetMarkerSize(0.0);                                  
    gr[n]->SetLineColor(gCOLOR[n]);                                    
    gr[n]->SetLineStyle(1);                                    
  }                                                                    

  TMultiGraph* mg = new TMultiGraph();
  char gTitle[256]; 
  sprintf(gTitle,"ROC");
  mg->SetTitle(gTitle); 
  for(int n=0;n<nROC;n++) mg->Add(gr[n]);  
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

  mg->GetHistogram()->GetXaxis()->SetRangeUser(xmin,xmax);
  mg->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);

  mg->GetHistogram()->GetYaxis()->SetNdivisions(509);

  mg->GetXaxis()->SetTitle(gr[0]->GetXaxis()->GetTitle());
  mg->GetXaxis()->SetLabelFont(42);                       
  mg->GetYaxis()->SetLabelFont(42);                       
  mg->GetXaxis()->SetTitleFont(42);                       
  mg->GetYaxis()->SetTitleFont(42);                       
  mg->GetXaxis()->SetTitleOffset(0.70);                   
  mg->GetYaxis()->SetTitleOffset(0.80);                   
  mg->GetXaxis()->SetTitleSize(0.06);                     
  mg->GetYaxis()->SetTitleSize(0.06);                     
  mg->GetXaxis()->SetTitle("p-value");
  mg->GetYaxis()->SetTitle("efficiency");
  mg->GetXaxis()->CenterTitle(true);
  mg->GetYaxis()->CenterTitle(true);

  mg->Draw("ALP");

  // draw the legend

  TLegend* leg;                
  double hleg = 0.18+nROC*0.06; 
  leg = new TLegend(0.4636591,0.1624473,0.8759398,hleg,NULL,"brNDC");

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

  for(int n=0;n<nROC;n++) {
    char legLabel[256];
    sprintf(legLabel,"%s",gNAME[n].Data());
    leg->AddEntry(gr[n],legLabel,"lp");
  }
  leg->Draw();

  // save plot
  if(ofname!="") canvas->Print(ofname);

  return;
}

int readStatistic(TString fname, vector<float>& vstat) {

  float stat;     

  ifstream in;
  in.open(fname.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fname.Data() << endl;exit(1);}

  while (1) {
    in >> stat;
    if (!in.good()) break;
    vstat.push_back(stat);
  }
  in.close();

  return vstat.size();
}

int readConfig(TString fname, vector<TString>& vfsignal, vector<TString>& vfnoise, vector<TString>& vname) {

  char fsignal[1024];
  char fnoise[1024];
  char name[1024];

  ifstream in;
  in.open(fname.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fname.Data() << endl;exit(1);}

  while (1) {
    in >> fsignal >> fnoise >> name;
    if (!in.good()) break;
    if(TString(fsignal[0])=="#") continue;
    vfsignal.push_back(fsignal);
    vfnoise.push_back(fnoise);
    vname.push_back(name);
  }
  in.close();

  return vname.size();
}

