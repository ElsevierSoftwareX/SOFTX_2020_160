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


// This macro is used to computed the false alarm probability of combined searches
// Uses as input parameters the files produced by the report sim when pp_fad is enabled
// 
// fapList :
// the fapList is the name of the file which contains the list of searches
// for each search the format is :
//
// $SEARCH TAG COLOR STYPE 
// - TAG   : label to be used in the legend
// - COLOR : line color
// - STYPE : line style
//
// $FAD FADvsRHO PRODvsFAD  
// - FADvsRHO  : file name of FAD vs RHO
// - PRODvsFAD : file name of PROD vs FAD
//
// $FAR FARvsRHO OBSTYPEvsFAR  
// - FARvsRHO     : file name of FAR vs RHO
// - OBSTIMEvsFAR : file name of OBSTIME vs FAR
//
// Example :
//
// $SEARCH LH_S6B_LF_LF    1       1
// $FAR    FARvsRHO_ALL.txt  OBSTIMEvsFAR_ALL.txt
// $FAD    FADvsRHO_ALL.txt  PRODvsFAD_ALL.txt
//
// fapMode :
// the fapMode is the mode used to produce the combined searches
// the fapMode can contains a combination of the following options : FAD/FAR/RHO
// the combined curves are obtained from the minimum FAP_FAD, FAP_FAR, FAP_RHO
// where FAP_XXX is the False Alarm Probability computed using the ranked 
// events according to XXX=(FAD/FAR/RHO)
//
// Examples :
//
// fadMode = "FAD"
// fadMode = "FAD FAR"
// fadMode = "FAD FAR RHO"
//
// Auxiliary parameters :
// rho_min : minimum rho showed in the plot (def=0, the value is computed from input files)
// rho_max : maximum rho showed in the plot (def=0, the value is computed from input files) 
// refline : vertical line showed as reference in the plot (def=0, not showed)
// reflabel: refernce label used in the legend (def="reference")
// pfile   : saved file name. If !="" -> save plot with name pfile (def="") 
//           if == "batch" the macro is executed in batch mode
//
// How to run (example) :
// root 'cwb_mkcfad.C("FAP.lst","FAD FAR RHO",5,40,7.35,"BigDog")'
//

#define MAX_FAP 20
#define SHOW_3SIGMA			// plot 3 sigma probability line
#define PROB_3SIGMA 	(1./370.398)	// percentage outside 3 sigma gaussian probability

// NOTE : FAX = FAD or FAR

//#define XAXIS_IS_IFAX	// uncomment to use as xaxis IFAX instead of RHO

double GetFAP(double rho, int n, int nFAP, TString* FAXvsRHO, TString* PRODvsFAX);

void GetFAPvsRHO(int n, vector<double>& x, vector<double>& y, 
                 int nFAP, TString* FAXvsRHO, TString* PRODvsFAX);

void GetFAPvsRHO(vector<double>& x, vector<double>& y, 
                 int nFAP, TString* FARvsRHO, TString* OBSTIMEvsFAR);

TGraph* GetGraph(vector<double> x, vector<double> y, 
                 TString ptitle, TString xtitle, TString ytitle);

double gRHO_MIN =  0;
double gRHO_MAX =  0;

double gFAP_MIN =  1.79769e+308;
double gFAP_MAX = -1.79769e+308;

void cwb_mkcfad(TString fapList, TString fapMode, double rho_min=0, double rho_max=0, 
                double refline=0, TString reflabel="reference", TString pfile="") {

  bool batch=false;
  if(pfile=="batch") {gROOT->SetBatch(true);pfile="";batch=true;}

  if(pfile!="" && !pfile.EndsWith(".png")) {
    cout << endl;
    cout << "cwb_mkcfad : Error in pfile name : " << pfile << endl;
    cout << "             file name must ends with .png" << endl << endl;
    exit(1);
  }

  gRHO_MIN = rho_min;
  gRHO_MAX = rho_max;

  TString FARvsRHO[MAX_FAP];
  TString OBSTIMEvsFAR[MAX_FAP];
  TString FADvsRHO[MAX_FAP];
  TString PRODvsFAD[MAX_FAP];
  TString FAP_LABEL[MAX_FAP];
  int     FAP_COLOR[MAX_FAP];
  int     FAP_STYLE[MAX_FAP];
  unsigned int FAP_TYPE[MAX_FAP]; 
  for(int i=0;i<MAX_FAP;i++) FAP_TYPE[i]=0;

  // Open list
  ifstream in;
  in.open(fapList.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fapList.Data() << endl;exit(1);}

  int size=0;
  char str[1024];
  int fpos=0;
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') size++;
  }
  cout << "size " << size << endl;
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);
  if (size==0) {cout << "Error : File " << fapList.Data() << " is empty" << endl;exit(1);}

  char iFARvsRHO[1024];
  char iOBSTIMEvsFAR[1024];
  char iFADvsRHO[1024];
  char iPRODvsFAD[1024];
  char iFAP_LABEL[1024];
  char iFAP_TAG[1024];
  int  iFAP_COLOR;
  int  iFAP_STYLE;
  int  nFAP=0;
  char line[1024];
  while(true) {
    in.getline(line,1024);
    if (in.eof()) break;
    std::stringstream linestream(line);
    if(TString(line).BeginsWith("#")) continue;
    if(TString(line)=="") continue;
    if(TString(line).BeginsWith("$SEARCH")) {
      if(FAP_TYPE[nFAP]>1) nFAP++;
      linestream >> iFAP_TAG >> iFAP_LABEL >> iFAP_COLOR >> iFAP_STYLE;
      FAP_LABEL[nFAP]=iFAP_LABEL;
      FAP_COLOR[nFAP]=iFAP_COLOR;
      FAP_STYLE[nFAP]=iFAP_STYLE;
      FAP_TYPE[nFAP]=1; 
    }
    if(TString(line).BeginsWith("$FAR")) {
      linestream >> iFAP_TAG >> iFARvsRHO >> iOBSTIMEvsFAR;
      FARvsRHO[nFAP]=iFARvsRHO;
      OBSTIMEvsFAR[nFAP]=iOBSTIMEvsFAR;
      FAP_TYPE[nFAP]|=2; 
    }
    if(TString(line).BeginsWith("$FAD")) {
      linestream >> iFAP_TAG >> iFADvsRHO >> iPRODvsFAD;
      FADvsRHO[nFAP]=iFADvsRHO;
      PRODvsFAD[nFAP]=iPRODvsFAD;
      FAP_TYPE[nFAP]|=4; 
    }
  }
  if(FAP_TYPE[nFAP]>1) nFAP++;
  in.close();

  for(int i=0;i<nFAP;i++) { 
    cout << "SEARCH : " << FAP_LABEL[i] << " " << FAP_COLOR[i] << " " << FAP_STYLE[i] << endl;
    if(FAP_TYPE[i]&2) cout << FARvsRHO[i] << " " << OBSTIMEvsFAR[i] << endl;
    if(FAP_TYPE[i]&4) cout << FADvsRHO[i] << " " << PRODvsFAD[i] << endl;
    cout << endl;
    for(int j=i;j<nFAP;j++) if(FAP_TYPE[i]!=FAP_TYPE[j]) {
      cout << "cwb_mkcfap : Error - searches not consistent !!!" << endl;
      cout << "             Missing FAR or FAD declaration" <<  endl << endl;
      exit(1);
    }
    if((!fapMode.Contains("RHO"))&&(!fapMode.Contains("FAR"))&&(!fapMode.Contains("FAD"))) {
      cout << "cwb_mkcfap : Error - fapMode not contains valid declarations (FAD/FAR/RHO)" << endl;
      exit(1);
    }
    if((!FAP_TYPE[i]&2)&&(fapMode.Contains("RHO"))) {
      cout << "cwb_mkcfap : Error - fapMode=RHO needs FAR declarations" << endl;
      exit(1);
    }
    if((!FAP_TYPE[i]&2)&&(fapMode.Contains("FAR"))) {
      cout << "cwb_mkcfap : Error - fapMode=FAR needs FAR declarations" << endl;
      exit(1);
    }
    if((!FAP_TYPE[i]&2)&&(fapMode.Contains("FAD"))) {
      cout << "cwb_mkcfap : Error - fapMode=FAD needs FAD declarations" << endl;
      exit(1);
    }
  }

  if((fapMode.Contains("RHO"))&&(!fapMode.Contains("FAR"))&&(!fapMode.Contains("FAD"))) {
    for(int i=0;i<nFAP;i++) FAP_COLOR[i]=2;	// red
  }

  // get FAP for RHO=refline
  if(refline>0) {
    cout << endl << "---------------------------------------------------------------------" << endl;
    for(int i=0;i<nFAP;i++) {
      double FAP;
      if(FAP_TYPE[i]&2) FAP = GetFAP(refline, i, nFAP, FARvsRHO, OBSTIMEvsFAR);
      if(FAP_TYPE[i]&4) FAP = GetFAP(refline, i, nFAP, FADvsRHO, PRODvsFAD);
      cout << "Search : " << FAP_LABEL[i] << "\t -> FAP at RHO = " << refline << " is " << FAP << endl;
    }
    cout << "---------------------------------------------------------------------" << endl << endl;
    if(batch) exit(0);
  }
 
  // create canvas
  TCanvas* gCANVAS;
  gCANVAS = new TCanvas("canvas", "canvas",16,30,825,546);
  gCANVAS->Range(-19.4801,-9.25,-17.4775,83.25);
  gCANVAS->SetBorderSize(2);
  gCANVAS->SetFrameFillColor(0);
  gCANVAS->SetGridx();
  gCANVAS->SetGridy();
//#ifdef XAXIS_IS_IFAD
  gCANVAS->SetLogx();
//#endif
  gCANVAS->SetLogy();
  gStyle->SetOptFit(kTRUE);

  // create graphs
  double rho_max=0;
  TGraph* gr[MAX_FAP];
  vector<double> xFAP[MAX_FAP],yFAP[MAX_FAP]; 
  vector<double> xFAD[MAX_FAP],yFAD[MAX_FAP]; 
  vector<double> xFAR[MAX_FAP],yFAR[MAX_FAP]; 
  vector<double> xRHO[MAX_FAP],yRHO[MAX_FAP]; 

  for(int i=0;i<nFAP;i++) {
    if(fapMode.Contains("FAD")) {
      // get FAP from events ranked by FAD
      GetFAPvsRHO(i, xFAD[i], yFAD[i], nFAP, FADvsRHO, PRODvsFAD);
      int xsize = xFAP[i].size();
      if(xsize) {for(int j=0;j<xsize;j++) if(yFAD[i][j]<yFAP[i][j]) yFAP[i][j]=yFAD[i][j];}
      else      {xFAP[i]=xFAD[i];yFAP[i]=yFAD[i];}
    }
    if(fapMode.Contains("FAR")) {
      // get FAP from events ranked by FAR
      GetFAPvsRHO(i, xFAR[i], yFAR[i], nFAP, FARvsRHO, OBSTIMEvsFAR);
      // select minimum between events ranked by FAR and by FAD
      int xsize = xFAP[i].size();
      if(xsize) {for(int j=0;j<xsize;j++) if(yFAR[i][j]<yFAP[i][j]) yFAP[i][j]=yFAR[i][j];}
      else      {xFAP[i]=xFAR[i];yFAP[i]=yFAR[i];}
    }
    if(fapMode.Contains("RHO")) {
      // get FAP from events ranked by RHO
      GetFAPvsRHO(xRHO[i], yRHO[i], nFAP, FARvsRHO, OBSTIMEvsFAR);
      // select minimum between events ranked by RHO and by FAR/FAD
      int xsize = xFAP[i].size();
      if(xsize) {for(int j=0;j<xsize;j++) if(yRHO[i][j]<yFAP[i][j]) yFAP[i][j]=yRHO[i][j];}
      else      {xFAP[i]=xRHO[i];yFAP[i]=yRHO[i];}
    }
  }

  for(int i=0;i<nFAP;i++) {
    int nMODE=0;
    TString subtitle=""; 
    if(fapMode.Contains("FAD")) {subtitle+=" FAD ";nMODE++;}
    if(fapMode.Contains("FAR")) {subtitle+=" FAR ";nMODE++;}
    if(fapMode.Contains("RHO")) {subtitle+=" RHO ";nMODE++;}
    TString title = "";
    if(nMODE==1) {
      title = TString::Format("False Alarm Probability - ranked by %s",subtitle.Data());
    } else {
      title = TString::Format("False Alarm Probability - minimum FAP ranked by (%s)",subtitle.Data());
    }
#ifdef XAXIS_IS_IFAD
    gr[i] = GetGraph(xFAP[i],yFAP[i],title,"IFAD","FAP");
#else
    gr[i] = GetGraph(xFAP[i],yFAP[i],title,"rho","FAP");
#endif
    gr[i]->SetLineColor(FAP_COLOR[i]);
    gr[i]->SetLineStyle(FAP_STYLE[i]);
    if(xFAP[i][xFAP[i].size()-1] > rho_max) rho_max=xFAP[i][xFAP[i].size()-1];
  }
  if((gRHO_MAX>0)&&(rho_max>gRHO_MAX)) rho_max = gRHO_MAX;

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(title);
  for(int i=0;i<nFAP;i++) mg->Add(gr[i]);

  // add reference line 
  if(refline>0) {
    gr[nFAP] = new TGraph;
    gr[nFAP]->SetPoint(0,refline,0);
    gr[nFAP]->SetPoint(1,refline,1);
    gr[nFAP]->SetLineColor(kGreen);
    gr[nFAP]->SetLineWidth(2);
    mg->Add(gr[nFAP]);
  }

  // add 3sigma reference line 
#ifdef SHOW_3SIGMA
  gr[nFAP+1] = new TGraph;
  gr[nFAP+1]->SetPoint(0,gRHO_MIN,PROB_3SIGMA);
  gr[nFAP+1]->SetPoint(1,rho_max,PROB_3SIGMA);
  gr[nFAP+1]->SetLineColor(kBlue);
  gr[nFAP+1]->SetLineWidth(2);
  gr[nFAP+1]->SetLineStyle(10);
  mg->Add(gr[nFAP+1]);
#endif  

  mg->Paint("APL");

  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetXaxis()->SetRangeUser(gRHO_MIN,gRHO_MAX);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  mg->GetHistogram()->GetYaxis()->SetRangeUser(0.5*gFAP_MIN,1.5*gFAP_MAX);

  mg->GetXaxis()->SetLabelFont(42);
  mg->GetYaxis()->SetLabelFont(42);
  mg->GetXaxis()->SetTitleFont(42);
  mg->GetYaxis()->SetTitleFont(42);
  mg->GetXaxis()->SetTitleOffset(1.20);
  mg->GetYaxis()->SetTitleOffset(1.05);
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleSize(0.04);
#ifdef XAXIS_IS_IFAD
  mg->GetXaxis()->SetTitle("IFAD");
#else
  mg->GetXaxis()->SetTitle("rho");
#endif
  mg->GetYaxis()->SetTitle("FAP");

  mg->GetXaxis()->SetMoreLogLabels();

  mg->Draw("ALP");

  // draw the legend
  TLegend *leg;
  double hleg = 0.84-nFAP*0.03;
  leg = new TLegend(0.7369062,hleg,0.9914738,0.9265385,NULL,"brNDC");
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

  for(int i=0;i<nFAP;i++) leg->AddEntry(gr[i],FAP_LABEL[i].Data(),"lp");
  if(refline>0) leg->AddEntry(gr[nFAP],reflabel.Data(),"lp");
#ifdef SHOW_3SIGMA
  leg->AddEntry(gr[nFAP+1],"3 sigma","lp");
#endif
  leg->Draw();

  // save plot
  if(pfile!="") {
    TString gfile=pfile;
    gfile.ReplaceAll(".png",".gif");
    gCANVAS->Print(gfile);
    char cmd[1024];
    sprintf(cmd,"convert %s %s",gfile.Data(),pfile.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",gfile.Data());
    gSystem->Exec(cmd);
    exit(0);
  }
}

double GetFAP(double rho, int n, int nFAP, TString* FAXvsRHO, TString* PRODvsFAX) {

  double RHO  = rho;
  double FAX  = CWB::Toolbox::GetStepFunction("Y",FAXvsRHO[n],RHO);
  double PROD = 0; 
  for(int j=0;j<nFAP;j++) PROD += CWB::Toolbox::GetStepFunction("Y",PRODvsFAX[j],FAX);
  double MU   = FAX*PROD;
  double FAP  = 1.-exp(-MU);

  return FAP;
}

void GetFAPvsRHO(int n, vector<double>& x, vector<double>& y, int nFAP, TString* FAXvsRHO, TString* PRODvsFAX) {

  x.clear();y.clear();

  double rho_min  = CWB::Toolbox::GetStepFunction("xmin",FAXvsRHO[n]);
  double rho_max  = CWB::Toolbox::GetStepFunction("xmax",FAXvsRHO[n]);
  if((gRHO_MIN>0)&&(rho_min<gRHO_MIN)) rho_min = gRHO_MIN;
  if((gRHO_MAX>0)&&(rho_max>gRHO_MAX)) rho_max = gRHO_MAX;

  double drho = 0.1;
  int    nrho = TMath::Nint((rho_max-rho_min)/drho)+2;

  for(int i=0;i<nrho;i++) {
    double RHO  = rho_min+i*drho;
    double FAX  = CWB::Toolbox::GetStepFunction("Y",FAXvsRHO[n],RHO);
    double PROD = 0; 
    for(int j=0;j<nFAP;j++) PROD += CWB::Toolbox::GetStepFunction("Y",PRODvsFAX[j],FAX);
    double MU   = FAX*PROD;
    double FAP  = 1.-exp(-MU);

#ifdef XAXIS_IS_IFAX
    x.push_back(1./FAX);
#else
    x.push_back(RHO);
#endif
    y.push_back(FAP);

    if((FAP!=0)&&(FAP<gFAP_MIN)) gFAP_MIN = FAP; // skip last y[i]=0 to avoid TGraph log issue
    if(FAP>gFAP_MAX) gFAP_MAX = FAP; 
  }
}

void GetFAPvsRHO(vector<double>& x, vector<double>& y, int nFAP, TString* FARvsRHO, TString* OBSTIMEvsFAR) {

  x.clear();y.clear();

  double rho_min = 1e10;
  double rho_max = 0;
  for(int j=0;j<nFAP;j++) {
    double _rho_min  = CWB::Toolbox::GetStepFunction("xmin",FARvsRHO[j]);
    double _rho_max  = CWB::Toolbox::GetStepFunction("xmax",FARvsRHO[j]);
    if(_rho_min<rho_min) rho_min=_rho_min;
    if(_rho_max>rho_max) rho_max=_rho_max;
  }
  if((gRHO_MIN>0)&&(rho_min<gRHO_MIN)) rho_min = gRHO_MIN;
  if((gRHO_MAX>0)&&(rho_max>gRHO_MAX)) rho_max = gRHO_MAX;

  double drho = 0.1;
  int    nrho = TMath::Nint((rho_max-rho_min)/drho)+2;

  double MIN_FAR[MAX_FAP];
  for(int j=0;j<MAX_FAP;j++) MIN_FAR[j]=1e10;
  for(int i=0;i<nrho;i++) {
    double RHO = rho_min+i*drho;
    double MU = 0;
    for(int j=0;j<nFAP;j++) {
      double FAR = CWB::Toolbox::GetStepFunction("Y",FARvsRHO[j],RHO);
      double OBSTIME = CWB::Toolbox::GetStepFunction("Y",OBSTIMEvsFAR[j],FAR);
      if(FAR>0) if(FAR<MIN_FAR[j]) MIN_FAR[j]=FAR;	// FAR min hold
      MU += MIN_FAR[j]*OBSTIME;
    }
    double FAP = 1.-exp(-MU);

    x.push_back(RHO);
    y.push_back(FAP);

    if((FAP!=0)&&(FAP<gFAP_MIN)) gFAP_MIN = FAP; // skip last y[i]=0 to avoid TGraph log issue
    if(FAP>gFAP_MAX) gFAP_MAX = FAP; 
  }
}

TGraph* 
GetGraph(vector<double> x, vector<double> y, TString ptitle, TString xtitle, TString ytitle) {

  int size = x.size();

  TGraph* gr = new TGraph;

  int cnt=0;

  for(int i=1;i<size;i++) {
    if(y[i]==0) continue;	// skip y[i]=0 to avoid TGraph log issue
    double dx = (x[i]-x[i-1])/10000.;
    gr->SetPoint(cnt++,x[i]-dx,y[i-1]);
    gr->SetPoint(cnt++,x[i]+dx,y[i]);
  }

  gr->GetHistogram()->GetXaxis()->SetTitle(xtitle.Data());
  gr->GetHistogram()->GetYaxis()->SetTitle(ytitle.Data());
  gr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
  // add more log labels when y range is small
  double max=0;double min=1e80;
  for(int i=0;i<size;i++) {
    if(y[i]==0) continue;	// skip y[i]=0 to avoid TGraph log issue
    if(y[i]>max) max=y[i]; if(y[i]<min && y[i]!=0) min=y[i];
  }
  if(max/min<10) {
    gr->GetHistogram()->GetYaxis()->SetMoreLogLabels();
  }
  gr->SetTitle(ptitle);
  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);

  return gr;
}
