// this macro produces the ROC report page
// the ROC is reported for each MDC @ fixed injection factor
// Author : Gabriele Vedovato
//
// run ex : root 'MakeROC.C("definitions.C")'
// definitions.C contains all definitions
//
// --------------------------------------------------------------
// ROC definitions
// --------------------------------------------------------------
// #define nROC N               // number of ROC compared
// #TString ROC_TITLE = "...";  // report title
// #TString ROC_DIR   = "...";  // output report directory (Ex: report/dump/roc)
// #TString ROC_WWW   = "...";  // www link of the ROC plots
//
// --------------------------------------------------------------
// reports SIM/BKG dir/www-links
// --------------------------------------------------------------
// TString ROC_NAME[nROC] = {"name1",...};
// TString SIM_WWW[nROC]  = {"https:sim_link1",...};
// TString BKG_WWW[nROC]  = {"https:bkg_link1",...};
// TString SIM_DIR[nROC]  = {"sim_dir1",...};
// TString BKG_DIR[nROC]  = {"bkg_dir1",...};
//
// --------------------------------------------------------------
// plot line definitions
// --------------------------------------------------------------
// int LINE_COLOR[nROC] = {4,2,3};	// Ex: {blue.red,green}
// int LINE_STYLE[nROC] = {1,1,1};	// Ex: {solid,solid,solid}
//

#define MDC_TYPE_NULL	100000
#define nROC_MAX        20

int nROC;
TString REP_WWW;
TString ROC_TITLE;
TString ROC_DIR;
TString ROC_WWW;
TString ROC_NAME[nROC_MAX];
int LINE_COLOR[nROC_MAX];
int LINE_STYLE[nROC_MAX];
TString SIM_DIR[nROC_MAX];
TString BKG_DIR[nROC_MAX];
TString SIM_WWW[nROC_MAX];
TString BKG_WWW[nROC_MAX];

int readParameters(TString fname, vector<TString>& mdc_name, 
                   vector<int>& mdc_type, vector<float>& factor);
int readParameters(TString fname, vector<double>& rho, vector<double>& par);
void PlotROC(int nroc, int ieff, TString mdc_name, float factor, TString ofDir);
void MakeHtml(TString roc_dir, int ieff, TString ROC_WWW, 
              TString* BKG_WWW, TString* SIM_WWW, TString* ROC_NAME);

void 
MakeROC(TString cfgFile, int ieff=50, int imdc_type=MDC_TYPE_NULL) {

  gROOT->Macro(cfgFile);

  char cmd[1024];

  char roc_dir[1024];
  sprintf(roc_dir,"%s_%d",ROC_DIR.Data(),ieff);

  if(ROC_DIR!="") CWB::Toolbox::mkDir(roc_dir,true);
  if(ROC_DIR!="") gROOT->SetBatch(true);

  // copy input list file to output plot dir
  sprintf(cmd,"cp %s %s/.",cfgFile.Data(),roc_dir);
  gSystem->Exec(cmd);

  vector<TString> mdc_name; 
  vector<int> mdc_type;
  vector<float> factor;
  char fname[1024];
  sprintf(fname,"%s/eff_%d_threshold_factors.txt",SIM_DIR[0].Data(),ieff);
  int nmdc = readParameters(fname, mdc_name, mdc_type, factor);
  // copy eff_threshold_factors.txt file to output plot dir
  sprintf(cmd,"cp %s %s/.",fname,roc_dir);
  gSystem->Exec(cmd);
  //cout << "nmdc : " << nmdc << endl;
  // check if it is consisten with the others parameters
  for(int n=1;n<nROC;n++) {
    vector<TString> xmdc_name; 
    vector<int> xmdc_type;
    vector<float> xfactor;
    char xfname[1024];
    sprintf(xfname,"%s/eff_%d_threshold_factors.txt",SIM_DIR[n].Data(),ieff);
    cout << "Process File : " << xfname << endl;
    int xnmdc = readParameters(xfname, xmdc_name, xmdc_type, xfactor);
    if(nmdc!=xnmdc) {
      cout << "Num parameters files inconsistents : " << xfname << endl;
      exit(1);
    }
    wavearray<int> check(xnmdc);
    check=-1;
    for(int i=0;i<nmdc;i++) {
      for(int j=0;j<nmdc;j++) {
        if(TString(xmdc_name[i]).Contains(mdc_name[j])) check[j]=i;
        if(TString(mdc_name[j]).Contains(xmdc_name[i])) check[j]=i; 
      }
    }
    for(int i=0;i<nmdc;i++) if(check[i]) if(factor[i]!=xfactor[check[i]]) check[i]=-1;
    for(int i=0;i<nmdc;i++) {
      if(check[i]==-1) {
        cout << "Parameters files inconsistents : " << endl;
        cout << "xmdc_name : " << xmdc_name[i] << " mdc_name : " << mdc_name[check[i]] << endl; 
        cout << "xfactor : " << xfactor[i] << " factor : " << factor[check[i]] << endl; 
        exit(1);
      }
    }
  }

  // list available mdc names
  if(imdc_type<0) {
    cout << "-------------------------------------------------" << endl;
    cout << "list of availables mdc types" << endl;
    cout << "-------------------------------------------------" << endl;
    for(int j=0;j<nmdc;j++) {
      cout << "mdc_type : " << mdc_type[j] << "\tmdc_name : " << mdc_name[j] << endl;
    }
    cout << "-------------------------------------------------" << endl;
    exit(1);
  }

  if((imdc_type!=MDC_TYPE_NULL)&&((imdc_type<0)||(imdc_type>=nmdc))) {
    cout << "mdc_type not allowed" << endl;
    exit(1);
  }

  for(int j=0;j<nmdc;j++) {
    if((imdc_type!=MDC_TYPE_NULL)&&(imdc_type!=mdc_type[j])) continue;
    PlotROC(nROC, ieff, mdc_name[j], factor[j], roc_dir);
  }

  if((ROC_DIR!="")&&(imdc_type==MDC_TYPE_NULL)) {
    MakeHtml(roc_dir, ieff, ROC_WWW, BKG_WWW, SIM_WWW, ROC_NAME);
    cout << endl;
    cout << "ROC report page : " << endl;
    cout << roc_dir << endl;
    cout << ROC_WWW << "_" << ieff << endl;
    cout << endl;
    gSystem->Exit(0);
  }
}

void PlotROC(int nroc, int ieff, TString mdc_name, float factor, TString ofDir) {

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
  canvas->SetLogx();                                                     
  canvas->SetGridx();                                                    
  canvas->SetGridy();                                                    
  canvas->SetFillColor(kWhite);                                          

  double eff_min=1;
  double eff_max=0;
  char fname[1024];
  TGraph* gr[nROC_MAX];
  for(int n=0;n<nROC;n++) {
    vector<double> rho_far;
    vector<double> rho_eff;
    vector<double> far;
    vector<double> eff;
    // read FAR
    sprintf(fname,"%s/rate_threshold_veto.txt",BKG_DIR[n].Data());
    cout << "rate_threshold : " << fname << endl;
    int far_size = readParameters(fname, rho_far, far);
    // read EFF
    sprintf(fname,"%s/eff_%d_threshold_%s.txt",SIM_DIR[n].Data(),ieff,mdc_name.Data());
    cout << "eff_threshold : " << fname << endl;
    int sim_size = readParameters(fname, rho_eff, eff);
    if(far_size!=sim_size) {
      cout << "rate_threshold not compatible with eff_threshold" << endl;
      cout << "far_size : " << far_size << " sim_size : " << sim_size << endl;
      gSystem->Exit(1);
    }
    for(int j=0;j<sim_size;j++) {
      if(rho_far[j]!=rho_eff[j]) {
        cout << "rate_threshold not compatible with eff_threshold" << endl;
        cout << "thresholds values are different" << endl;
        cout << j << " rho_far : " << rho_far[j] << " =rho_eff : " << rho_eff[j] << endl;
        gSystem->Exit(1);
      } 
    }
    int wsize=0;
    wavearray<double> wfar(far_size);
    wavearray<double> weff(sim_size);
    for(int j=0;j<sim_size;j++) {
      wfar[wsize]=far[j];
      weff[wsize]=eff[j];
      if((wfar[wsize]>0)&&(weff[wsize]>0)) {
         if(weff[wsize]<eff_min) eff_min=weff[wsize]; 
         if(weff[wsize]>eff_max) eff_max=weff[wsize]; 
         //cout << wsize << " " << wfar[wsize] << " " << weff[wsize] << endl;
         wsize++;
      }
    }
    gr[n] = new TGraph(wsize,wfar.data,weff.data);
  }
  for(int n=0;n<nROC;n++) {                                            
    gr[n]->SetLineWidth(2);                                            
    gr[n]->SetMarkerColor(LINE_COLOR[n]);                                  
    gr[n]->SetMarkerStyle(20);                                  
    gr[n]->SetLineColor(LINE_COLOR[n]);                                    
    gr[n]->SetLineStyle(LINE_STYLE[n]);                                    
  }                                                                    

  TMultiGraph* mg = new TMultiGraph();
  char gTitle[256]; 
  sprintf(gTitle,"ROC : mdc = %s : factor = %g",mdc_name.Data(),factor);
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

//  mg->GetHistogram()->GetXaxis()->SetRangeUser(1e-9,1e-3);
//  mg->GetHistogram()->GetYaxis()->SetRangeUser(eff_min,eff_max);
  mg->GetHistogram()->GetYaxis()->SetRangeUser(0,1);

  mg->GetXaxis()->SetTitle(gr[0]->GetXaxis()->GetTitle());
  mg->GetXaxis()->SetLabelFont(42);                       
  mg->GetYaxis()->SetLabelFont(42);                       
  mg->GetXaxis()->SetTitleFont(42);                       
  mg->GetYaxis()->SetTitleFont(42);                       
  mg->GetXaxis()->SetTitleOffset(1.20);                   
  mg->GetYaxis()->SetTitleOffset(1.20);                   
  mg->GetXaxis()->SetTitleSize(0.04);                     
  mg->GetYaxis()->SetTitleSize(0.04);                     
  mg->GetXaxis()->SetTitle("FAR");
  mg->GetYaxis()->SetTitle("Efficiency");

  mg->Draw("ALP");

  // draw the legend

  TLegend* leg;                
  double hleg = 0.15+nROC*0.05; 
  leg = new TLegend(0.61,hleg,0.99,0.15,NULL,"brNDC");

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

  for(int n=0;n<nROC;n++) {
    char legLabel[256];
    sprintf(legLabel,"%s",ROC_NAME[n].Data());
    leg->AddEntry(gr[n],legLabel,"lp");
  }
  leg->Draw();

  // save plot
  if(ofDir!="") {
    char gfileName[1024];
    sprintf(gfileName,"%s/roc_%s.gif",ofDir.Data(),mdc_name.Data());
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
  }

  return;
}

int readParameters(TString fname, vector<TString>& mdc_name, 
    vector<int>& mdc_type, vector<float>& factor) {
  // read factor50 from custom input file             

  char xmdc_name[256];
  int  xmdc_type;     
  float xfactors;     

  ifstream in;
  in.open(fname.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fname.Data() << endl;exit(1);}

  while (1) {
    in >> xmdc_name >> xmdc_type >> xfactors;
    if (!in.good()) break;
    //cout << xmdc_name <<" "<<  xmdc_type <<" "<< xfactors << endl;
    mdc_name.push_back(xmdc_name);
    mdc_type.push_back(xmdc_type);
    factor.push_back(xfactors);
  }
  in.close();

  return mdc_name.size();
}

int readParameters(TString fname, vector<double>& rho, vector<double>& par) {
  // read rho,par input file             

  rho.clear();
  par.clear();

  double xrho;     
  double xpar;     

  ifstream in;
  in.open(fname.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fname.Data() << endl;exit(1);}

  while (1) {
    in >> xrho >> xpar;
    if (!in.good()) break;
    //cout << xrho <<" "<<  xpar << endl;
    rho.push_back(xrho);
    par.push_back(xpar);
  }
  in.close();

  return rho.size();
}

void MakeHtml(TString roc_dir, int ieff, TString ROC_WWW, 
              TString* BKG_WWW, TString* SIM_WWW, TString* ROC_NAME) {

  vector<TString> mdc_name;
  vector<int> mdc_type;
  vector<float> factor;
  char fname[1024];
  sprintf(fname,"%s/eff_%d_threshold_factors.txt",roc_dir.Data(),ieff);
  int nmdc = readParameters(fname, mdc_name, mdc_type, factor);

  char roc_html_file[1024];sprintf(roc_html_file,"%s/index.html",roc_dir.Data());
  ofstream out;
  out.open(roc_html_file,ios::out);             // create index.html
  char oline[1024];
  out << "<html>" << endl;
  out << "<font color=\"black\" style=\"font-weight:bold;\"><center><p><h2>"
      << ROC_TITLE << "</h2><p><center></font>" << endl;

  for(int n=0;n<nROC;n++) {
    out << "<h4>" << ROC_NAME[n] << " : " << endl;
    sprintf(oline,"<a target=\"_parent\" href=\"%s\">bkg report page</a>",BKG_WWW[n].Data());
    out << oline << endl;
    out << "---" << endl;
    sprintf(oline,"<a target=\"_parent\" href=\"%s\">sim report page</a>",SIM_WWW[n].Data());
    out << oline << endl;
    out << "</h4>" << endl;
  }

  out << "<hr><br>" << endl;
  for(int i=0;i<mdc_name.size();i++) {
    out << "<table border=0 cellpadding=2 align=\"center\">" << endl;
    out << "<tr align=\"center\">" << endl;
    out << "<td><font style=\"font-weight:bold;\"><center><p><h2>"
        << mdc_name[i] << "</h2><p><center></font></td>" << endl;
    out << "</tr>" << endl;
    out << "<tr align=\"center\">" << endl;
    sprintf(oline,"<td><a href=\"%s_%d/roc_%s.png\"><img src=\"%s_%d/roc_%s.png\" width=650></a></td>",
            ROC_WWW.Data(),ieff,mdc_name[i].Data(),ROC_WWW.Data(),ieff,mdc_name[i].Data());
    out << oline << endl;
    out << "</tr>" << endl;
    out << "</table>" << endl;
  }
  out << "</html>" << endl;
  out.close();

  return;
}

