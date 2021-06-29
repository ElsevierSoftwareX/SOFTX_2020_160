// This macro is used to computed the false alarm probability of combined searches
// Uses as input parameters the files produced by the report sim when pp_fad is enabled
// 
// the fadList is the name of the file which contains the list of searches
// for each line is a search and it is define as :
// 1 - file name of FAD vs RHO
// 2 - file name of PROD vs FAD	
// 3 - label to be used in the legend
// 4 - line color
// 5 - line style
//
// NOTE : if PROD is constant for every rho then it is possible to input the value 
//        of PROD instead of file name of PROD vs FAD
//
//
// Example :
// S6A_FADvsRHO_ALL_lt.txt   S6A_PRODvsFAD_ALL_lt.txt  S6A_INET_LT_0d5     1       9
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
// root 'CombineSearchesWithFAD.C("FAD.lst",5,40,7.35,"BigDog")'
//

#define MAX_FAP 20
#define SHOW_3SIGMA			// plot 3 sigma probability line
#define PROB_3SIGMA 	0.002699796063	// percentage outside 3 sigma gaussian probability
#define SHOW_4SIGMA			// plot 4 sigma probability line
#define PROB_4SIGMA 	0.000063342484	// percentage outside 4 sigma gaussian probability
//#define SHOW_5SIGMA			// plot 5 sigma probability line
#define PROB_5SIGMA 	0.000000573303	// percentage outside 5 sigma gaussian probability

//#define XAXIS_IS_IFAD	// uncomment to use as xaxis IFAD instead of RHO

double GetFAP(double rho, int n, int nFAP, TString* FADvsRHO, TString* PRODvsFAD);

void GetFAPvsRHO(int n, vector<double>& x, vector<double>& y, 
                 int nFAP, TString* FADvsRHO, TString* PRODvsFAD);

TGraph* GetGraph(vector<double> x, vector<double> y, 
                 TString ptitle, TString xtitle, TString ytitle);

double gRHO_MIN =  0;
double gRHO_MAX =  0;

double gFAP_MIN =  1.79769e+308;
double gFAP_MAX = -1.79769e+308;

void CombineSearchesWithFAD(TString fadList, double rho_min=0, double rho_max=0, 
                            double refline=0, TString reflabel="reference", TString pfile="") {

  bool batch=false;
  if(pfile=="batch") {gROOT->SetBatch(true);pfile="";batch=true;}

  if(pfile!="" && !pfile.EndsWith(".png")) {
    cout << endl;
    cout << "CombineSearchesWithFAD : Error in pfile name : " << pfile << endl;
    cout << "                         file name must ends with .png" << endl << endl;
    exit(1);
  }

  gRHO_MIN = rho_min;
  gRHO_MAX = rho_max;

  TString FADvsRHO[MAX_FAP];
  TString PRODvsFAD[MAX_FAP];
  TString FAP_LABEL[MAX_FAP];
  int     FAP_COLOR[MAX_FAP];
  int     FAP_STYLE[MAX_FAP];

  // Open list
  ifstream in;
  in.open(fadList.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fadList.Data() << endl;exit(1);}

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
  if (size==0) {cout << "Error : File " << fadList.Data() << " is empty" << endl;exit(1);}

  char iFADvsRHO[1024];
  char iPRODvsFAD[1024];
  char iFAP_LABEL[1024];
  int  iFAP_COLOR;
  int  iFAP_STYLE;
  int nFAP=0;
  char line[1024];
  while(true) {
    in.getline(line,1024);
    if (in.eof()) break;
    std::stringstream linestream(line);
    if((line[0]=='#')||(TString(line)=="")) continue;
    linestream >> iFADvsRHO >> iPRODvsFAD >> iFAP_LABEL >> iFAP_COLOR >> iFAP_STYLE;
    FADvsRHO[nFAP]=iFADvsRHO;
    PRODvsFAD[nFAP]=iPRODvsFAD;
    FAP_LABEL[nFAP]=iFAP_LABEL;
    FAP_COLOR[nFAP]=iFAP_COLOR;
    FAP_STYLE[nFAP]=iFAP_STYLE;
    cout << nFAP+1 << " " << FADvsRHO[nFAP] << " " << PRODvsFAD[nFAP] 
         << " " << FAP_LABEL[nFAP] << " " << FAP_COLOR[nFAP] << " " << FAP_STYLE[nFAP] << endl;
    nFAP++;
  }
  in.close();

  // get FAP for RHO=refline
  if(refline>0) {
    cout << endl << "---------------------------------------------------------------------" << endl;
    for(int i=0;i<nFAP;i++) {
      double FAP = GetFAP(refline, i, nFAP, FADvsRHO, PRODvsFAD);
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
  rho_max=0;
  TGraph* gr[MAX_FAP];
  vector<double> x[MAX_FAP]; 
  vector<double> y[MAX_FAP]; 
  for(int i=0;i<nFAP;i++) {
    GetFAPvsRHO(i, x[i], y[i], nFAP, FADvsRHO, PRODvsFAD);
#ifdef XAXIS_IS_IFAD
    gr[i] = GetGraph(x[i],y[i],"False Alarm Probability","IFAD","FAP");
#else
    gr[i] = GetGraph(x[i],y[i],"False Alarm Probability","rho","FAP");
#endif
    gr[i]->SetLineColor(FAP_COLOR[i]);
    gr[i]->SetLineStyle(FAP_STYLE[i]);
    if(x[i][x[i].size()-1] > rho_max) rho_max=x[i][x[i].size()-1];
  }
  if((gRHO_MAX>0)&&(rho_max>gRHO_MAX)) rho_max = gRHO_MAX;

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle("False Alarm Probability");
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

  // add 4sigma reference line 
#ifdef SHOW_4SIGMA
  gr[nFAP+2] = new TGraph;
  gr[nFAP+2]->SetPoint(0,gRHO_MIN,PROB_4SIGMA);
  gr[nFAP+2]->SetPoint(1,rho_max,PROB_4SIGMA);
  gr[nFAP+2]->SetLineColor(kRed);
  gr[nFAP+2]->SetLineWidth(2);
  gr[nFAP+2]->SetLineStyle(10);
  mg->Add(gr[nFAP+2]);
#endif  

  // add 5sigma reference line 
#ifdef SHOW_5SIGMA
  gr[nFAP+3] = new TGraph;
  gr[nFAP+3]->SetPoint(0,gRHO_MIN,PROB_5SIGMA);
  gr[nFAP+3]->SetPoint(1,rho_max,PROB_5SIGMA);
  gr[nFAP+3]->SetLineColor(kGreen);
  gr[nFAP+3]->SetLineWidth(2);
  gr[nFAP+3]->SetLineStyle(10);
  mg->Add(gr[nFAP+3]);
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
  double hleg = 0.9-nFAP*0.03;
  leg = new TLegend(0.7369062,hleg,0.9914738,0.9846154,NULL,"brNDC");
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
#ifdef SHOW_4SIGMA
  leg->AddEntry(gr[nFAP+2],"4 sigma","lp");
#endif
#ifdef SHOW_5SIGMA
  leg->AddEntry(gr[nFAP+3],"5 sigma","lp");
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

double GetFAP(double rho, int n, int nFAP, TString* FADvsRHO, TString* PRODvsFAD) {

  double RHO  = rho;
  double FAD  = CWB::Toolbox::GetStepFunction("Y",FADvsRHO[n],RHO);
  double PROD = 0; 
  for(int j=0;j<nFAP;j++) {
    if(PRODvsFAD[j].IsFloat()) {	// read PROD from file name
      PROD += PRODvsFAD[j].Atof();
    } else {				// read PROD from file
      PROD += CWB::Toolbox::GetStepFunction("Y",PRODvsFAD[j],FAD);
    }
  }
  double MU   = FAD*PROD;
  double FAP  = 1.-exp(-MU);

  return FAP;
}

void GetFAPvsRHO(int n, vector<double>& x, vector<double>& y, int nFAP, TString* FADvsRHO, TString* PRODvsFAD) {

  double rho_min  = CWB::Toolbox::GetStepFunction("xmin",FADvsRHO[n]);
  double rho_max  = CWB::Toolbox::GetStepFunction("xmax",FADvsRHO[n]);
  if((gRHO_MIN>0)&&(rho_min<gRHO_MIN)) rho_min = gRHO_MIN;
  if((gRHO_MAX>0)&&(rho_max>gRHO_MAX)) rho_max = gRHO_MAX;

  double drho = 0.1;
  int    nrho = TMath::Nint((rho_max-rho_min)/drho)+2;

  for(int i=0;i<nrho;i++) {
    double RHO  = rho_min+i*drho;
    double FAD  = CWB::Toolbox::GetStepFunction("Y",FADvsRHO[n],RHO);
    double PROD = 0; 
    for(int j=0;j<nFAP;j++) {
      if(PRODvsFAD[j].IsFloat()) {	// read PROD from file name
        PROD += PRODvsFAD[j].Atof();
      } else {				// read PROD from file
        PROD += CWB::Toolbox::GetStepFunction("Y",PRODvsFAD[j],FAD);
      }
    }
    double MU   = FAD*PROD;
    double FAP  = 1.-exp(-MU);

#ifdef XAXIS_IS_IFAD
    x.push_back(1./FAD);
#else
    x.push_back(RHO);
#endif
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
