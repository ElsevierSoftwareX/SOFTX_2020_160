#define L1_SHH_FILE_NAME "report/dump/sensitivity_L1_841104609_S5_SIM_BRST_L1H1H2_SGW_SMDC_run5_tst1_job1.txt"
#define H1_SHH_FILE_NAME "report/dump/sensitivity_H1_841104609_S5_SIM_BRST_L1H1H2_SGW_SMDC_run5_tst1_job1.txt"
#define H2_SHH_FILE_NAME "report/dump/sensitivity_H2_841104609_S5_SIM_BRST_L1H1H2_SGW_SMDC_run5_tst1_job1.txt"

#define GPS_TIME 841104609

//#define GREYSCALE

#define SIZE 1000000

#define nIFO 3

void DrawSensitivitiesS5(TString ofname="") {

  ifstream in[nIFO];

  in[0].open(L1_SHH_FILE_NAME,ios::in);
  if (!in[0].good()) {cout << "Error Opening File : " << L1_SHH_FILE_NAME << endl;exit(1);}

  in[1].open(H1_SHH_FILE_NAME,ios::in);
  if (!in[1].good()) {cout << "Error Opening File : " << H1_SHH_FILE_NAME << endl;exit(1);}

  in[2].open(H2_SHH_FILE_NAME,ios::in);
  if (!in[2].good()) {cout << "Error Opening File : " << H2_SHH_FILE_NAME << endl;exit(1);}

  double freq[nIFO][SIZE]; 
  double shh[nIFO][SIZE]; 
  int size[nIFO]={0,0,0};

  for(int i=0;i<nIFO;i++) {
    while (1) {
      double F,S;
      in[i] >> F >> S;
      if (!in[i].good()) break;
      if(S!=0) {
        freq[i][size[i]]=F; 
        shh[i][size[i]]=S;
        size[i]++;
      }
    }
    in[i].close();
  }

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gStyle->SetTitleFont(72);
  gStyle->SetMarkerColor(50);
  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleW(0.98);
  gStyle->SetTitleH(0.05);
  gStyle->SetTitleY(0.98);
  gStyle->SetFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleFont(12,"D");

  TGraph* gr[nIFO];
  for(int n=0;n<nIFO;n++) gr[n] = new TGraph(size[n],freq[n],shh[n]);
  for(int n=0;n<nIFO;n++) gr[n]->SetLineWidth(3);

#ifdef GREYSCALE
  gr[0]->SetLineColor(kBlack);
  gr[0]->SetMarkerColor(kBlack);
  gr[0]->SetLineStyle(1);

  gr[1]->SetLineColor(kBlack);
  gr[1]->SetMarkerColor(kBlack);
  gr[1]->SetLineStyle(2);

  gr[2]->SetLineColor(kBlack);
  gr[2]->SetMarkerColor(kBlack);
  gr[2]->SetLineStyle(9);
#else
  gr[0]->SetLineColor(kBlue);
  gr[0]->SetMarkerColor(kBlue);

  gr[1]->SetLineColor(kBlack);
  gr[1]->SetMarkerColor(kBlack);

  gr[2]->SetLineColor(kRed);
  gr[2]->SetMarkerColor(kRed);
#endif

  TCanvas *canvas = new TCanvas("Sensitivity", "Sh", 300,40, 1200, 800);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetLogx();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);
#ifdef GREYSCALE
  canvas->SetGrayscale(true);
#endif

  TMultiGraph* mg = new TMultiGraph();
  char title[256];sprintf(title,"LIGO Detectors : S5 Sensitivity Curves - GPS : %d",GPS_TIME);
  mg->SetTitle(title);
  //mg->SetTitle("");

  mg->Add(gr[0]);
  mg->Add(gr[1]);
  mg->Add(gr[2]);

  mg->Paint("APL");

  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);

  mg->GetHistogram()->GetXaxis()->SetRangeUser(30,8192);
  mg->GetHistogram()->GetYaxis()->SetRangeUser(1e-23,1e-19);

  mg->GetXaxis()->SetTitle(gr[0]->GetXaxis()->GetTitle());
  mg->GetXaxis()->SetLabelFont(42);
  mg->GetYaxis()->SetLabelFont(42);
  mg->GetXaxis()->SetTitleFont(42);
  mg->GetYaxis()->SetTitleFont(42);
  mg->GetXaxis()->SetTitleOffset(1.20);
  mg->GetYaxis()->SetTitleOffset(1.05);
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitle("Frequency (Hz)      ");
  mg->GetYaxis()->SetTitle("#frac{1}{#sqrt{Hz}}          ");

  mg->Draw("APL");

  // draw the legend
  TLegendEntry* entry;
  //TLegend *legend = new TLegend(0.6580773,0.1684699,0.8810704,0.3261206,NULL,"brNDC");
  TLegend *legend = new TLegend(0.6705686,0.7269985,0.8929766,0.8853695,NULL,"brNDC");

  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(kBlack);
  legend->SetFillColor(kWhite);
  //legend->SetLineWidth(3);
  legend->AddEntry(gr[0],"      L1","lp");
  legend->AddEntry(gr[1],"      H1","lp");
  legend->AddEntry(gr[2],"      H2","lp");
  legend->Draw();

  if(ofname!="") {
    cout << ofname.Data() << endl;
    ofname.ReplaceAll(".png",".gif");
    canvas->Print(ofname.Data());
    char cmd[256];
    TString pfname(ofname);
    pfname.ReplaceAll(".gif",".png");
    sprintf(cmd,"convert %s %s",ofname.Data(),pfname.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",ofname.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    exit(0);
  }

}
