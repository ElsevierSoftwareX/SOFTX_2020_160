
#define nFAR_MAX        32

int nFAR;
TString FAR_PLOT_NAME;
TString FAR_TITLE;
TString FAR_NAME[nFAR_MAX];
int LINE_COLOR[nFAR_MAX];
int LINE_STYLE[nFAR_MAX];
int LINE_MARKER[nFAR_MAX];
double OBS_TIME[nFAR_MAX];
TString FAR_PATH[nFAR_MAX];
double ZERO_LAG_TIME[nFAR_MAX];

int readParameters(TString fname, wavearray<double>& RHO, wavearray<double>& eRHO, wavearray<double>& FAR, wavearray<double>& eFAR);
int writeParameters(TString fname, wavearray<double>& RHO, wavearray<double>& FAR);
void PlotFAR(int nfar, TString far_plot_name, double rho_min=-1, double rho_max=-1);

TString PLOT_PATH = "/home/vedovato/O1/ER8b_12Sep20Oct_C0101/ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run0a/report/dump/";

TString FARvsLowRHO_PATH = "/home/vedovato/O1/ER8b_12Sep20Oct_C0101/ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run0a/report/dump/";
TString FARvsRHO_PATH = "/home/vedovato/O1/ER8b_12Sep20Oct_C0101/ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a/report/dump/";

void 
Draw_FARvsRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_LTDQ() {


  Style_t markers[32]= {20,21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,
                        21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,20 };

  Color_t colors[32] = { 8, 0, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7,
                         6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7 };

  for(int i=0;i<nFAR_MAX;i++) {
    LINE_COLOR[i]  = colors[i];
    LINE_STYLE[i]  = 1;
    LINE_MARKER[i] = markers[i];
    ZERO_LAG_TIME[i] = -1; 
    OBS_TIME[i]  = 1;
  }

  FAR_TITLE="";
  FAR_PLOT_NAME="";

// ===========================================================================
// C01+C01+LTDQ
// ===========================================================================

  #define MODE "inclusive"
  nFAR = 6;
  FAR_PATH[0]=FARvsRHO_PATH+"FARvsRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_unmodeled_LTDQ.txt";
  FAR_PATH[1]=FARvsRHO_PATH+"FARvsRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_constrained_LTDQ.txt";
  FAR_PATH[2]=FARvsRHO_PATH+"FARvsRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_chirp_LTDQ.txt";
  FAR_PATH[3]=FARvsLowRHO_PATH+"FARvsLowRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_unmodeled_LTDQ.txt";
  FAR_PATH[4]=FARvsLowRHO_PATH+"FARvsLowRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_constrained_LTDQ.txt";
  FAR_PATH[5]=FARvsLowRHO_PATH+"FARvsLowRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_chirp_LTDQ.txt";

/*
  #define MODE "exclusive"
  nFAR = 6;
  FAR_PATH[0]=FARvsRHO_PATH+"FARvsRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_eunmodeled_LTDQ.txt";
  FAR_PATH[1]=FARvsRHO_PATH+"FARvsRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_econstrained_LTDQ.txt";
  FAR_PATH[2]=FARvsRHO_PATH+"FARvsRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_chirp_LTDQ.txt";
  FAR_PATH[3]=FARvsLowRHO_PATH+"FARvsLowRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_eunmodeled_LTDQ.txt";
  FAR_PATH[4]=FARvsLowRHO_PATH+"FARvsLowRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_econstrained_LTDQ.txt";
  FAR_PATH[5]=FARvsLowRHO_PATH+"FARvsLowRHO_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a_chirp_LTDQ.txt";
*/

  PlotFAR(nFAR, FAR_PLOT_NAME, -1, -1);

}

void PlotFAR(int nfar, TString far_plot_name, double rho_min, double rho_max) {

  // create plots
  gStyle->SetFrameBorderMode(0);     // remove the red box around canvas
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

  TCanvas *canvas = new TCanvas("roc", "roc", 300,40, 800, 500);
  canvas->Clear();                                                       
  canvas->ToggleEventStatus();                                           
  canvas->SetLogy();                                                     
  canvas->SetGridx();                                                    
  canvas->SetGridy();                                                    
  canvas->SetFillColor(kWhite);                                          

  double far_min=1e10;
  double rmin=1e20;
  double rmax=0;
  char fname[1024];
  int nGR=0;
  int line_style[3*nFAR_MAX];
  int line_marker[3*nFAR_MAX];                                    
  int line_color[3*nFAR_MAX];                                    
  TString FAR_name[3*nFAR_MAX];
  bool sigma_lines[3*nFAR_MAX];
  TGraphErrors* gr[3*nFAR_MAX];
  for(int n=0;n<nFAR;n++) {
    cout << n << " OBS_TIME : " << OBS_TIME[n] << endl;
    if(OBS_TIME[n]>0) if(far_min>1./OBS_TIME[n]) far_min=1./OBS_TIME[n];
    vector<double> rho_far;
    vector<double> far;
    // read FAR
    sprintf(fname,"%s",FAR_PATH[n].Data());
    cout << "rate_threshold : " << fname << endl;
    wavearray<double> RHO;
    wavearray<double> FAR;
    wavearray<double> eRHO;
    wavearray<double> eFAR;
    int far_size = readParameters(fname, RHO, eRHO, FAR, eFAR);
#if ( MODE == "exclusive")
    TString ofname = fname;
    ofname.ReplaceAll(".txt","_comp_paper.txt");
    writeParameters(ofname, RHO, FAR); 
#endif
    if(rmin<RHO[0]) rmin=RHO[0];
    if(rmax>RHO[far_size-1]) rmax=RHO[far_size-1];
    sigma_lines[nGR] = false;
    line_style[nGR]=LINE_STYLE[n];
    line_marker[nGR]=LINE_MARKER[n];                                    
    line_color[nGR]=LINE_COLOR[n];                                    
    FAR_name[nGR]=FAR_NAME[n];
    gr[nGR++] = new TGraphErrors(far_size,RHO.data,FAR.data,eRHO.data,eFAR.data);
    gr[nGR-1]->SetMarkerStyle(20);
    gr[nGR-1]->SetMarkerSize(0.35);
    if(n==0) gr[nGR-1]->SetMarkerColor(8);
    if(n==1) gr[nGR-1]->SetMarkerColor(kBlack);
    if(n==2) gr[nGR-1]->SetMarkerColor(kRed);
    if(n==3) gr[nGR-1]->SetMarkerColor(8);
    if(n==4) gr[nGR-1]->SetMarkerColor(kBlack);
    if(n==5) gr[nGR-1]->SetMarkerColor(kRed);
    gr[nGR-1]->SetName(TString::Format("gr%d",nGR-1));
    
  }

  // GW150914
  RHO.data[0]=14.1*sqrt(2); FAR.data[0]=1./(1369200.00/(365*24*3600));
  eRHO.data[0]=0; eFAR.data[0]=0;
  gr[nGR] = new TGraphErrors(1,RHO.data,FAR.data,eRHO.data,eFAR.data);
  gr[nGR]->SetMarkerColor(4);
  gr[nGR]->SetMarkerStyle(22);
  gr[nGR]->SetMarkerSize(2);
  FAR_name[nGR]="GW150914";
  nGR++;


  TMultiGraph* mg = new TMultiGraph();
  char gTitle[256]; 
  if(FAR_TITLE!="") sprintf(gTitle,FAR_TITLE.Data());
  else              sprintf(gTitle,"FAR Comparison");
// PAPER
#if ( MODE == "inclusive")
  sprintf(gTitle,"ER8b/O1 (C0101) background of unmodeled/constrained/chirp data sets (Sep 12  - Oct 20)");
#endif
#if ( MODE == "exclusive")
  sprintf(gTitle,"ER8b/O1 (C0101) background of edge/bulk/chirp data sets (Sep 12  - Oct 20)");
#endif
  sprintf(gTitle,"");
  mg->SetName("mg");
  mg->SetTitle(gTitle); 
  for(int n=0;n<nGR;n++) mg->Add(gr[n]);  
  //mg->Paint("APL");                        
  mg->Paint("AL");                        

  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);  
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);  
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5); 
//  mg->GetHistogram()->SetMinimum(far_min/2.); 

  if(rho_min>=0) rmin=rho_min*sqrt(2);
  if(rho_max>=0) rmax=rho_max*sqrt(2);
  rmax=16*sqrt(2);
  mg->GetHistogram()->GetXaxis()->SetRangeUser(rmin,rmax);

  mg->GetXaxis()->SetTitle(gr[0]->GetXaxis()->GetTitle());
  mg->GetXaxis()->SetLabelFont(42);                       
  mg->GetYaxis()->SetLabelFont(42);                       
  mg->GetXaxis()->SetTitleFont(42);                       
  mg->GetYaxis()->SetTitleFont(42);                       
  mg->GetXaxis()->SetTitleOffset(1.20);                   
  mg->GetYaxis()->SetTitleOffset(1.20);                   
  mg->GetXaxis()->SetTitleSize(0.04);                     
  mg->GetYaxis()->SetTitleSize(0.04);                     
  mg->GetXaxis()->CenterTitle(true);            //PAPER
  mg->GetYaxis()->CenterTitle(true);            //PAPER
  mg->GetXaxis()->SetTitle("Coherent Network SNR ( #eta_{c} )");
  //mg->GetXaxis()->SetTitle("#eta_{c}");
  mg->GetYaxis()->SetTitle("FAR ( yr^{-1} )");

  mg->Draw("ALP");

  // draw label GW150914
  //TLatex *pS1 = new TLatex(RHO[0], FAR[0], "GW150914");       //PAPER
  //TLatex *pS1 = new TLatex(18.30,2.14, "GW150914");       //PAPER
  TLatex *pS1 = new TLatex(18.93,2.61, "GW150914");       //PAPER
  pS1->SetTextFont(52);
  pS1->SetTextSize(0.040);      //PAPER
  pS1->SetLineWidth(2);
  pS1->SetTextColor(1);
//  pS1->Draw();

  TPad *pad = canvas->GetPad(0);;
  pad->cd();  
  TFile *froot = new TFile("chirp_set_low_high_rho_paper.root", "RECREATE");
  pad->Write("pad");
  froot->Close();

  // draw the legend

  TLegend* leg;                
  //leg = new TLegend(0.65,0.24,0.87,0.40,NULL,"brNDC");
  //leg = new TLegend(0.55,0.14,0.87,0.50,NULL,"brNDC");
//  leg = new TLegend(0.55,0.14,0.87,0.40,NULL,"brNDC");
//  leg = new TLegend(0.404,0.495,0.704,0.855,NULL,"brNDC");	// 3 sets
//  leg = new TLegend(0.5791457,0.6751055,0.879397,0.8586498,NULL,"brNDC");
//  leg = new TLegend(0.3278894,0.6603376,0.6281407,0.8438819,NULL,"brNDC");
  leg = new TLegend(0.4120603,0.6265823,0.6306533,0.8607595,NULL,"brNDC");


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

cout << nGR << endl;
  for(int n=0;n<nGR;n++) {
    if(n==3) continue;
    if(n==4) continue;
    if(n==5) continue;
    char legLabel[256];
#if ( MODE == "inclusive")
    //if(n==0) strcpy(legLabel,"inclusive UnModeled");
    //if(n==1) strcpy(legLabel,"inclusive Constrained");
    if(n==0) strcpy(legLabel,"UnModeled");
    if(n==1) strcpy(legLabel,"Constrained");
#endif
#if ( MODE == "exclusive")
    //if(n==0) strcpy(legLabel,"exclusive UnModeled");
    //if(n==1) strcpy(legLabel,"exclusive Constrained");
    //if(n==0) strcpy(legLabel,"Edge");
    //if(n==1) strcpy(legLabel,"Bulk");
    if(n==0) strcpy(legLabel,"B1");
    if(n==1) strcpy(legLabel,"B2");
#endif
    if(n==2) strcpy(legLabel,"Chirp");
    //if(n==2) strcpy(legLabel,"B3");
    if(n==6) strcpy(legLabel,"GW150914");
    //sprintf(legLabel,"%s",FAR_name[n].Data());
    TLegendEntry* eleg = leg->AddEntry(gr[n],legLabel,"lp");
    //eleg->SetMarkerSize(5);
    //eleg->SetMarkerColor(kBlack);
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
  leg->Draw();

   // draw markers in the legend (the defaults are too small)

  TMarker *mP1 = new TMarker(13.2658,2512.81, 20);
  mP1->SetMarkerStyle(20);
  mP1->SetMarkerSize(1);
  mP1->SetMarkerColor(kGreen);
  mP1->Draw();

  TMarker *mP2 = new TMarker(13.2658,479.617, 20);
  mP2->SetMarkerStyle(20);
  mP2->SetMarkerSize(1);
  mP2->SetMarkerColor(kBlack);
  mP2->Draw();

  TMarker *mP3 = new TMarker(13.2658,91.7268, 20);
  mP3->SetMarkerStyle(20);
  mP3->SetMarkerSize(1);
  mP3->SetMarkerColor(kRed);
  mP3->Draw();

  // save plot
#if ( MODE == "inclusive")
  far_plot_name = PLOT_PATH+"CompareFAR_C0101_UnModeled_Constrained_Chirp_run1ato10a_inclusive_LTDQ";
#endif
#if ( MODE == "exclusive")
  far_plot_name = PLOT_PATH+"CompareFAR_C0101_eUnModeled_eConstrained_Chirp_run1ato10a_exclusive_LTDQ";
#endif
  if(far_plot_name!="") {
    char gfileName[1024];
    sprintf(gfileName,"%s.gif",far_plot_name.Data());
    canvas->Print(gfileName);
    TString pfileName=gfileName;
    pfileName.ReplaceAll(".gif",".png");
    char cmd[1024];
    sprintf(cmd,"convert %s %s",gfileName,pfileName.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",gfileName);
    cout << cmd << endl;
    gSystem->Exec(cmd);
//    gSystem->Exit(0);
  }

  return;
}

int readParameters(TString fname, wavearray<double>& RHO, wavearray<double>& eRHO, wavearray<double>& FAR, wavearray<double>& eFAR) {

  double rho;     
  double erho;     
  double far;     
  double efar;     

  double year = (24.*3600.*365.);

  RHO.resize(1000000);
  eRHO.resize(1000000);
  FAR.resize(1000000);
  eFAR.resize(1000000);

  ifstream in;
  in.open(fname.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fname.Data() << endl;exit(1);}

  int n=0;
  while (1) {
    //in >> rho >> far >> erho >> efar;
    in >> rho >> far;
    if (!in.good()) break;
rho*=sqrt(2.);
if(rho>22.) continue;
far*=year;
    RHO[n] = rho;
    eRHO[n] = erho;
    FAR[n] = far;
//    eFAR[n] = efar;
    eFAR[n] = 0;

    n++;
  }
  in.close();

  RHO.resize(n); 
  eRHO.resize(n); 
  FAR.resize(n); 
  eFAR.resize(n); 

  return n;
}

int writeParameters(TString fname, wavearray<double>& RHO, wavearray<double>& FAR) {

  ofstream out;
  out.open(fname.Data(),ios::out);
  if (!out.good()) {cout << "Error Opening File : " << fname.Data() << endl;exit(1);}

  for(int i=0;i<RHO.size();i++) {
    out << RHO[i] << "\t" << FAR[i] << endl;	
  }

  out.close();
}


