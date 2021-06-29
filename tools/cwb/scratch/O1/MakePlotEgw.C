#define nINJ_MAX	40
#define nSET_MAX	40

// -------------------------------------------------------------------------------------
// IFAR=8
// -------------------------------------------------------------------------------------

#define IFAR	8

//#define MODE "Inclusive"
//#define IFILE	"ifit_parameters_ifar8_ER8b_12Sep20Oct_C01_SIM1_BRST1_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST1_LF_L1H1.inj"
//#define IFILE	"ifit_parameters_ifar8_ER8b_12Sep20Oct_C01_SIM1_BRST2_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST2_LF_L1H1.inj"
//#define IFILE	"ifit_parameters_ifar8_ER8b_12Sep20Oct_C01_SIM1_BRST3_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST3_LF_L1H1.inj"

#define MODE "Exclusive"
//#define IFILE	"efit_parameters_ifar8_ER8b_12Sep20Oct_C01_SIM1_BRST1_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST1_LF_L1H1.inj"
//#define IFILE	"efit_parameters_ifar8_ER8b_12Sep20Oct_C01_SIM1_BRST2_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST2_LF_L1H1.inj"
//#define IFILE	"efit_parameters_ifar8_ER8b_12Sep20Oct_C01_SIM1_BRST3_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST3_LF_L1H1.inj"
#define IFILE	"xfit_parameters_ifar8_ER8b_12Sep20Oct_C0101_SIM1_LALBRST2_LF_rMRA_run0a.txt"
#define MDC_INJ_FILE "O1_BRST2_LF_L1H1.inj"


// -------------------------------------------------------------------------------------
// IFAR=100
// -------------------------------------------------------------------------------------
/*
#define IFAR	100

//#define MODE "Inclusive"
//#define IFILE	"ifit_parameters_ifar100_ER8b_12Sep20Oct_C01_SIM1_BRST1_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST1_LF_L1H1.inj"
//#define IFILE	"ifit_parameters_ifar100_ER8b_12Sep20Oct_C01_SIM1_BRST2_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST2_LF_L1H1.inj"
//#define IFILE	"ifit_parameters_ifar100_ER8b_12Sep20Oct_C01_SIM1_BRST3_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST3_LF_L1H1.inj"

#define MODE "Exclusive"
//#define IFILE	"efit_parameters_ifar100_ER8b_12Sep20Oct_C01_SIM1_BRST1_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST1_LF_L1H1.inj"
//#define IFILE	"efit_parameters_ifar100_ER8b_12Sep20Oct_C01_SIM1_BRST2_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST2_LF_L1H1.inj"
//#define IFILE	"efit_parameters_ifar100_ER8b_12Sep20Oct_C01_SIM1_BRST3_LF_rMRA_run0a.txt"
//#define MDC_INJ_FILE "O1_BRST3_LF_L1H1.inj"
#define IFILE	"xfit_parameters_ifar100_ER8b_12Sep20Oct_C0101_SIM1_LALBRST2_LF_rMRA_run0a.txt"
#define MDC_INJ_FILE "O1_BRST2_LF_L1H1.inj"
*/

void PlotEgw(int ninj, double* freq, double* egw, int nset, size_t* iset, TString* set);

void MakePlotEgw() {

  double C  = watconstants::SpeedOfLightInVacuo();      // m s^-1
  double G  = watconstants::GravitationalConstant();    // N m^2 kg^-2 (N=Kg m s^-2)
  double Mo = watconstants::SolarMass();                // Kg
  double Pc = watconstants::Parsec();                   // m
  double Pi = TMath::Pi();


  // read injection file types
  char    imdc_set[nINJ_MAX][128];    // injection set
  size_t  imdc_type[nINJ_MAX];        // injection type
  char    imdc_name[nINJ_MAX][128];   // injection name
  double  imdc_fcentral[nINJ_MAX];    // injection central frequencies
  double  imdc_fbandwidth[nINJ_MAX];  // injection bandwidth frequencies
  size_t  imdc_index[nINJ_MAX];       // type reference array
  size_t  imdc_iset[nINJ_MAX];        // injection set index
  double  imdc_egw[nINJ_MAX];         // injection Egw

  int ninj=ReadInjType(MDC_INJ_FILE,nINJ_MAX,imdc_set,imdc_type,
                       imdc_name,imdc_fcentral,imdc_fbandwidth);

  TString imdc_sset[nINJ_MAX];
  for(int i=0;i<ninj;i++) imdc_sset[i] = imdc_set[i];

  TString* imdc_set_name = new TString[ninj];
  int nset=0;
  for(int i=0;i<ninj;i++) {
    bool bnew=true;
    for(int j=0;j<nset;j++) if(imdc_set[i]==imdc_set_name[j]) bnew=false;
    if(bnew) imdc_set_name[nset++]=imdc_set[i];
  }
  cout << "nset : " << nset << endl;
  for(int i=0;i<nset;i++) {
    for(int j=0;j<ninj;j++) if(imdc_set[j]==imdc_set_name[i]) imdc_iset[j]=i;
  }

  int ecount[nINJ_MAX];
  TString piumeno[nINJ_MAX];
  float chi2[nINJ_MAX], err[nINJ_MAX], par1[nINJ_MAX], par2[nINJ_MAX], par3[nINJ_MAX];
  double hrss50[nINJ_MAX];
  TString ewaveform[nINJ_MAX];

  ifstream in;
  in.open(IFILE,ios::in);
  if (!in.good()) {cout << "Error Opening File : " << IFILE << endl;exit(1);}

  for(int k=0; k<ninj; k++) {

    in >> ecount[k] >> chi2[k] >> hrss50[k] >> piumeno[k]
       >> err[k] >> par1[k] >> par2[k] >> par3[k] >> ewaveform[k];
    if (!in.good()) break;

//    cout << k << "\t" << ewaveform[k] << "\t" << hrss50[k] << "\t" << imdc_fcentral[k] << "\t" << imdc_name[k] <<  endl;

    double Ro = 10000*Pc;                                 // m
    double Fo = imdc_fcentral[k];                         // s^-1
    double Ho = hrss50[k] ;                               // s^(1/2)
    double Egw = pow(Pi,2) * pow(C,3)/G * pow(Ro*Fo*Ho,2);
    cout << imdc_name[k] << "\tEgw(Mo) : " << Egw/(Mo*C*C) 
         << "\thrss50 : " << hrss50[k] << "\t" << imdc_fcentral[k] << endl;

    imdc_egw[k] = Egw/(Mo*C*C);
  }

  PlotEgw(ninj, imdc_fcentral, imdc_egw, nset, imdc_iset, imdc_sset);

}

void PlotEgw(int ninj, double* freq, double* egw, int nset, size_t* iset, TString* set) {

  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetCanvasColor(kWhite);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  TCanvas* canvas;

  TGraphErrors* gr50[nINJ_MAX];
  TMultiGraph*  mg;
  TLegend*      legend;

  TString set_name[nINJ_MAX];

  mg     = new TMultiGraph();
  //legend = new TLegend(0.6676955,0.1163636,0.8919753,0.4436364,NULL,"brNDC");
  legend = new TLegend(0.1111111,0.5563636,0.3353909,0.8836364,NULL,"brNDC");

  for(int k=0;k<nset;k++) {

    gr50[k] = new TGraphErrors();
    gr50[k]->SetLineColor(2);
    gr50[k]->SetLineWidth(1);
    gr50[k]->SetLineStyle(7);

    if(k==0) gr50[k]->SetMarkerColor(kBlack);
    if(k==1) gr50[k]->SetMarkerColor(kRed);
    if(k==2) gr50[k]->SetMarkerColor(kGreen);
    if(k==3) gr50[k]->SetMarkerColor(kBlue);
    if(k==4) gr50[k]->SetMarkerColor(kMagenta);

    if(k==0) gr50[k]->SetMarkerStyle(20);
    if(k==1) gr50[k]->SetMarkerStyle(21);
    if(k==2) gr50[k]->SetMarkerStyle(22);
    if(k==3) gr50[k]->SetMarkerStyle(23);
    if(k==4) gr50[k]->SetMarkerStyle(28);

    if(k==0) gr50[k]->SetMarkerSize(1.4);
    if(k==1) gr50[k]->SetMarkerSize(1.2);
    if(k==2) gr50[k]->SetMarkerSize(1.5);
    if(k==3) gr50[k]->SetMarkerSize(1.5);
    if(k==4) gr50[k]->SetMarkerSize(1.5);

    for(int i=0; i<ninj; i++) {
      if(iset[i]==k) gr50[k]->SetPoint(i,freq[i],egw[i]);
      if(iset[i]==k) set_name[k]=set[i];
    }
    mg->Add(gr50[k]);
  }

  canvas = new TCanvas("Egw","Egw",125,82,976,576);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetLogy();
  canvas->SetLogx(true);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);
  canvas->cd();

  mg->SetTitle(TString::Format("Egw vs Frequency : %s - IFAR = %d",MODE,IFAR));
  mg->Paint("alp");
  mg->GetHistogram()->GetXaxis()->SetTitle("Frequency (Hz)");
  mg->GetHistogram()->GetXaxis()->CenterTitle(true);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetXaxis()->SetTitleFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetTitleFont(42);
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.20);
  mg->GetYaxis()->SetNdivisions(10);
  mg->GetHistogram()->GetYaxis()->SetTitle("Egw (Mo)");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(32,1024+256);
#if (IFAR == 8) 
  mg->GetHistogram()->GetYaxis()->SetRangeUser(2e-10,2e-6);
#endif
#if (IFAR == 100) 
  mg->GetHistogram()->GetYaxis()->SetRangeUser(5e-10,5e-6);
#endif

  mg->Draw("ap");

  legend->SetBorderSize(1);
  legend->SetTextAlign(22);
  legend->SetTextFont(12);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(1001);
  legend->SetTextSize(0.04);
  legend->SetLineColor(kBlack);
  legend->SetFillColor(kWhite);

  for(int k=0;k<nset;k++) {
    legend->AddEntry(gr50[k],set_name[k].Data(),"lp");
  }

  legend->Draw();

  
  TString ofname = IFILE;
  ofname.ReplaceAll(".txt",".gif");
  cout << ofname << endl;
  canvas->Print(ofname);
  char cmd[256];
  TString pfname(ofname);
  pfname.ReplaceAll(".gif",".png");
  sprintf(cmd,"convert %s %s",ofname.Data(),pfname.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"rm %s",ofname.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);
}
