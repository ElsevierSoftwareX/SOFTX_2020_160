// Plot Error Regions Percentage vs Coverage 
// In this examples we use 4 mdc types at SNR = 10*sqrt(nIFO)

#define NMDC 4

#define NETSNR (10*sqrt(3.))

//#define nIFO  2
#define nIFO  3


double factors[3] = {10.*sqrt(nIFO), 20.*sqrt(nIFO), 30.*sqrt(nIFO)};
TString mdc[4] = {"WNB250_100_0d100","SG235Q3","SG235Q8d9","SGC235Q9"};

double GetPercentage(int idmdc, int iderA, int idfactor, TString fname);

void DrawCoverageVsPercentagePRC(TString ifile, int idfactor, bool save=false)
{

  TObjArray* token = TString(ifile).Tokenize(TString('/'));
  TObjString* sfile = (TObjString*)token->At(token->GetEntries()-1);
  TString TITLE = sfile->GetString();
  TString ofile = sfile->GetString();
  ofile.ReplaceAll(".root",".txt");
  TITLE.ReplaceAll(".root","");

  ofstream out;
  out.open(ofile.Data(),ios::out);
  if (!out.good()) {cout << "Error Opening File : " << ofile.Data() << endl;exit(1);}
  cout << "Create file : " << ofile.Data() << endl;

  // if coverage > percentage then erA is too large, the true erA is smaller (erA is conservative) 
  // if coverage < percentage then erA is too small, the true erA is greater (erA is optimistic) 

  double coverage[NMDC][11];
  for (int i=0;i<NMDC;i++) {
    coverage[i][0]=0.;
    coverage[i][10]=100.;
    for (int j=1;j<10;j++) {
      coverage[i][j]=GetPercentage(i+1,j,idfactor,ifile);
      out << mdc[i].Data() << " " << j*10 << " " << coverage[i][j] << endl;
      //cout << mdc[i].Data() << " " << j*10 << " " << coverage[i][j] << endl;
    }
  }
  out.close();

  TCanvas* canvas = new TCanvas("fom", "PRC", 300,40, 600, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetFillColor(0);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLogy(false);
  canvas->SetLogx(false);

  gStyle->SetTitleH(0.032);
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
  sprintf(title,"%s  -  netSNR = %3.2f",TITLE.Data(),NETSNR);
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
    sprintf(label,"%s ",mdc[i].Data());
    leg->AddEntry(gr[i],label,"lp");
  }
  leg->Draw();

  if(save) {
    char label[64];sprintf(label,"_%g_CovVsPerc.gif",factors[idfactor]);
    TString gfileName=ofile;
    gfileName.ReplaceAll(".txt",label);
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
    exit(0);
  }
}

double GetPercentage(int idmdc, int iderA, int idfactor, TString fname) {

  TFile *ifile = TFile::Open(fname.Data());
  TTree* itree = (TTree *) gROOT->FindObject("waveburst");
  itree->SetEstimate(itree->GetEntries());

  char selection[1024];
  char tree_cut[1024];
  sprintf(tree_cut,"abs(time[0]-time[3])<0.1 && type[1]==%d  && abs(factor-%f)<0.1",idmdc,factors[idfactor]);
  //cout << tree_cut << endl;

  sprintf(selection,"erA[0]-erA[%d]>>hist_cumulative_erA1(2000,-100,100)",iderA);
  itree->Draw(selection,tree_cut,"goff");
  int size = itree->GetSelectedRows();
  TH1D* hist_cumulative_erA1 = (TH1D*)gDirectory->Get("hist_cumulative_erA1");
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

