//
// Draw Cumulative fraction of the sky as a function of the perc of error region
// Author : Gabriele Vedovato
// Ex : broot 'DrawLarsHistogramPRC.C("skyloc.lst","MEDIAN50","ADV_SIM_BRST_LF_L1H1V1_2G.M1.png")'
// ilist_file format : skymap.root label     
// skymap.root is a root file produced by the macro DrawSkyMapPRC.C
//

#define APPLY_FIT

#define MAX_FILES 20


void DrawLarsHistogramPRC(TString ilist_file, TString pType="", TString ofName="") {

  if((ofName!="")&&(!ofName.Contains(".png"))) {
    cout << "Error Outpu File : " << ilist_file.Data() << endl;
    cout << "Must have .png extension" << endl;
    exit(1);
  } else {
    ofName.ReplaceAll(".png",".gif");
  }

  TString ptitle="Cumulative fraction of the sky as a function of the ";

  if((pType!="MEDIAN50")&&(pType!="MEDIAN90")&&(pType!="WRC50")&&(pType!="EFFICIENCY")&&(pType!="")) {
    ptitle=pType;
  } else {
    if(pType=="") pType="MEDIAN50";
    if(pType=="MEDIAN50")   ptitle=ptitle+"50% error region";
    if(pType=="MEDIAN90")   ptitle=ptitle+"90% error region";
    if(pType=="WRC50")      ptitle=ptitle+"50% normalized residual energy";
    if(pType=="EFFICIENCY") ptitle=ptitle+"efficiency";
  }

  TString ifile[MAX_FILES];
  TString title[MAX_FILES];

  // Open list
  ifstream in;
  in.open(ilist_file.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << ilist_file.Data() << endl;exit(1);}

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
  if (size==0) {cout << "Error : File " << ilist_file.Data() << " is empty" << endl;exit(1);}

  char sfile[256];
  char stitle[256];
  int NFILE=0;
  while(true) {
    in >> sfile >> stitle;
    if(!in.good()) break;
    if(sfile[0]=='#') continue;
    ifile[NFILE]=sfile;
    title[NFILE]=stitle;
    cout << NFILE+1 << " " << ifile[NFILE] << " " << title[NFILE] << endl;
    NFILE++;
  }
  in.close();

#ifdef APPLY_FIT
  // sigmoid fit function
  TF1 *f1;
  f1=new TF1("logNfit",logNfit,2,1024,5);
  f1->SetParameters(16,0.7,1.,1.,0);
  f1->SetParNames("ptype","sigma","betam","betap");
  f1->SetParLimits(0,log10(1.),log10(1024.));
  f1->SetParLimits(1,0.1,10.);
  f1->SetParLimits(2,0.1,4.);
  f1->SetParLimits(3,0.1,4.);
#endif

  TCanvas* canvas;
  canvas = new TCanvas("lars histogram", "canvas", 35, 46, 817, 472);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetGridx(false);
  canvas->SetGridx(false);
  canvas->SetFillColor(kWhite);
  canvas->SetLeftMargin(0.08);
  canvas->SetBottomMargin(0.13);
  canvas->SetBorderMode(0);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gROOT->SetStyle("Plain");
  gPad->UseCurrentStyle();
  gROOT->ForceStyle();

  gStyle->SetTitleFont(12,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  //gStyle->SetTitleAlign(11);
  gStyle->SetTitleW(1.0);
  gStyle->SetTitleH(0.050);
  gStyle->SetTitleY(0.98);

  gnetwork gNET[MAX_FILES];
  gskymap*  gSM[MAX_FILES];
  for(int n=0;n<NFILE;n++) {
    cout << "Process File : " << ifile[n].Data() << endl;
    gNET[n].LoadObject((char*)ifile[n].Data());
    gSM[n] = gNET[n].GetGskymap();
  }

  int L = gSM[0]->size();

  double pi = TMath::Pi();
  double sr=0;
  sr = cos(gSM[0]->theta_1*pi/180.)-cos(gSM[0]->theta_2*pi/180.);
  sr*= (gSM[0]->phi_2-gSM[0]->phi_1)*180/pi/L;
  sr = double(int(1000*sr))/1000.;
  cout << "sky resolution (deg^2) : " << sr << endl; 

  double mean=0;
  for(int l=0;l<L;l++) mean+=gSM[0]->get(l);
  mean/=L;
  cout << mean << endl;
  double integral = 0;

  Color_t colors[32] = { 6, 3, 2, 8,43, 7, 8, 4, 797, 2,43, 1, 3, 1, 6, 7,
                         6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7 };

  int lstyle[32] = {2, 0, 0, 0, 9, 9, 9, 2, 2, 2, 2,43, 1, 3, 1, 6, 7,
                         6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6 };

  double ptype[MAX_FILES];
  TH1F* h[MAX_FILES];
  for(int n=0;n<NFILE;n++) {
    cout << "Process : " << title[n].Data() << endl;
    if(pType=="EFFICIENCY") {
      h[n] = new TH1F(title[n].Data(),title[n].Data(),20,0,1.0001);  // NOTE! 1.0001 fix hist when filled with value=1 
    } else if(pType=="WRC50") {
      h[n] = new TH1F(title[n].Data(),title[n].Data(),20,0,0.3);
    } else {
      h[n] = new TH1F(title[n].Data(),title[n].Data(),20,0,log(1024)/log(2));
    }
    //h[n]->SetTitle(title[n].Data());
    h[n]->SetMarkerStyle(20);
    h[n]->SetMarkerSize(1);
    h[n]->SetMarkerColor(colors[n]);
    h[n]->SetLineColor(colors[n]);
    // h[n]->SetLineStyle(lstyle[n]);
    h[n]->SetLineWidth(2);
    h[n]->SetStats(kFALSE);
    h[n]->GetXaxis()->SetLabelFont(42);
    h[n]->GetYaxis()->SetLabelFont(42);
    h[n]->GetYaxis()->SetRangeUser(0,1.0);
    h[n]->GetXaxis()->SetTitleFont(42);
    h[n]->GetYaxis()->SetTitleFont(42);
    if(pType=="EFFICIENCY") {
      h[n]->GetXaxis()->SetTitle("efficiency           ");
    } else if(pType=="WRC50") {
      h[n]->GetXaxis()->SetTitle("normalized residual energy");
    } else {
      h[n]->GetXaxis()->SetTitle("#sigma^{2} (deg^{2})           ");
    }
    h[n]->GetYaxis()->SetTitle("#frac{#Omega}{4#pi}");
    if(pType=="EFFICIENCY") {
      h[n]->GetXaxis()->SetRangeUser(0.5,1.0);
      for(int l=0;l<L;l++) {
        double sigma2 = gSM[n]->get(l);
        h[n]->Fill(sigma2);
      }
    } else if(pType=="WRC50") {
      h[n]->GetXaxis()->SetRangeUser(0,0.3);
      for(int l=0;l<L;l++) {
        double sigma2 = gSM[n]->get(l);
        h[n]->Fill(sigma2);
      }
    } else {
      h[n]->GetXaxis()->SetRangeUser(0,log(1024)/log(2));
      for(int l=0;l<L;l++) {
        double sigma2 = log(pow(gSM[n]->get(l),2))/log(2);
        if(sigma2<0) sigma2=0.25; 
        if(gSM[n]->get(l)==0) sigma2=-0.25; 
        h[n]->Fill(sigma2);
      }
    }

    // Normalization
    integral = 0;
    for (int i=0;i<=h[n]->GetNbinsX();i++) integral+=h[n]->GetBinContent(i);
    for (int i=0;i<=h[n]->GetNbinsX();i++) h[n]->SetBinContent(i,h[n]->GetBinContent(i)/integral);
    for (int i=2;i<=h[n]->GetNbinsX();i++) h[n]->SetBinContent(i,h[n]->GetBinContent(i)+h[n]->GetBinContent(i-1));

    if(n==0) {
      if(pType=="EFFICIENCY" || pType=="WRC50") {
        h[0]->GetXaxis()->SetLabelFont(42);
        h[0]->GetXaxis()->SetLabelOffset(0.01);
        h[0]->GetXaxis()->LabelsOption("h");
      } else {
        char xlabel[16];
        for (int i=0;i<=h[0]->GetNbinsX();i++) {
          double xvalue = h[0]->GetXaxis()->GetBinCenter(i)+h[0]->GetXaxis()->GetBinWidth(i)/2.;
          if(i%4==0) {
            sprintf(xlabel,"%2.0f",pow(2,xvalue));
          } else {
            sprintf(xlabel,"");
          }
          h[0]->GetXaxis()->SetBinLabel(i,xlabel);
          h[0]->GetXaxis()->SetLabelSize(0.07);
          h[0]->GetXaxis()->SetLabelFont(42);
          h[0]->GetXaxis()->SetLabelOffset(0.01);
          h[0]->GetXaxis()->LabelsOption("h");
        }
      }
      h[0]->SetTitle(ptitle);
      h[0]->Draw("lp");
    } else h[n]->Draw("lpsame"); 

#ifdef APPLY_FIT
    int size = h[n]->GetNbinsX();
    double* x = new double[size];
    double* y = new double[size];
    if(pType=="EFFICIENCY" || pType=="WRC50") {
      for (int i=1;i<=size;i++) {
        x[i-1] = h[n]->GetXaxis()->GetBinCenter(i);
        y[i-1] = h[n]->GetBinContent(i);
      }
    } else {
      for (int i=1;i<=size;i++) {
        x[i-1] = h[n]->GetXaxis()->GetBinCenter(i)+h[n]->GetXaxis()->GetBinWidth(i)/2.;
        x[i-1] = pow(2,x[i-1]);
        y[i-1] = h[n]->GetBinContent(i);
        //cout << i-1 << " : " << x[i-1] << " " << y[i-1] << endl; 
      }
    }
    TGraph* gr = new TGraph(size,x,y);
    // Fit with sigmoid
    if(pType=="EFFICIENCY" || pType=="WRC50") {
      TGraphSmooth* gs = new TGraphSmooth("normal");
      wavearray<double> xs(1000);for(int i=0;i<xs.size();i++) xs[i]=double(i)/double(xs.size());
      TGraph* grs = gs->Approx(gr,"linear", xs.size(), xs.data);
      // find ptype[n] for ys=0.5
      for(int i=0;i<xs.size();i++) {
        double ys;
        grs->GetPoint(i,ptype[n],ys);
        if(ys>0.5) break;
      }
    } else {
      f1->SetLineColor(colors[n]);
      f1->SetNpx(10000);
      gr->Fit("logNfit");
      ptype[n]=pow(10.,f1->GetParameter(0));
      double chi2=f1->GetChisquare();
      double sigma=f1->GetParameter(1);
      double betam=f1->GetParameter(2);
      double betap=f1->GetParameter(3);
      cout << "Fit ptype% " << ptype[n] << endl;
      //f1->Draw("same");
    }
    delete gr; 
#endif
  }

  // draw the legend
  double hleg = 0.20+NFILE*0.05;
  TLegend* leg;  
  if(pType=="EFFICIENCY" || pType=="MEDIAN90") {
    double hleg = 0.8-NFILE*0.05;
    leg = new TLegend(0.1291513,hleg,0.6482165,0.8738739,NULL,"brNDC");
  } else {
    double hleg = 0.20+NFILE*0.05;
    leg = new TLegend(0.5202952,0.1390135,0.9901599,hleg,NULL,"brNDC");
  }

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

  for(int n=0;n<NFILE;n++) {
    char legLabel[256];
#ifdef APPLY_FIT
    sprintf(legLabel,"%s [ %2.2f ]",title[n].Data(),ptype[n]);
#else
    sprintf(legLabel,"%s",title[n].Data());
#endif
    leg->AddEntry(h[n],legLabel,"lp");
  }
  leg->Draw();

  canvas->SetGridy(true);

  if(ofName!="") {
    TString gfileName=ofName;
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
  }
  
  gSystem->Exit(0);
}

