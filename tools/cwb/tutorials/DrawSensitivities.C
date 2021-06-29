//
// Draw Detector Sensitivity Curves
// Author : Gabriele Vedovato
//
// how to use it
// root -l 'DrawSensitivities.C("strain_psd.lst")'
//
// strain_psd.lst format
//
// strain psd file name							title
// ../plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt                advLIGO_NSNS_Opt
//../plugins/strains/AdV_May2010.txt                                    AdV_May2010
//

#define MAX_PSD 10
#define MAX_SIZE 10000

void DrawSensitivities(TString ilist_file="", TString ofName="", TString title="", 
                       double min_freq=-1, double max_freq=-1, double min_strain=-1, double max_strain=-1) {

  if(ilist_file=="") {
    cout << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << "HOW TO USE : " << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << "root 'DrawSensitivities(\"list_of_strain_psd_file_names.txt\",\"otput_psd_file_name.png\", \\"<< endl;
    cout << "                        \"psd_title\", min_freq, max_freq, min_strain, max_strain)" << endl;
    cout << endl;
    cout << "only list_of_strain_psd_file_names.txt is mandatory" << endl;
    cout << "if list_of_strain_psd_file_names.txt=\"\" then print this help" << endl; 
    cout << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << "this is the strain_psd.lst format" << endl;
    cout << "strain psd file name (list of freq & strain)             title" << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << "../plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt    advLIGO_NSNS_Opt" << endl;
    cout << "plugins/strains/AdV_May2010.txt                          AdV_May2010" << endl;
    cout << endl;
    exit(1);
  }

  if((ofName!="")&&(!ofName.Contains(".png"))) {
    cout << "Error Outpu File : " << ilist_file.Data() << endl;
    cout << "Must have .png extension" << endl;
    exit(1);
  } else {
    ofName.ReplaceAll(".png",".gif");
  }


  TString ifile[MAX_PSD];
  TString ltitle[MAX_PSD];

  // Read PSD file list
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
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);
  if (size==0) {cout << "Error : File " << ilist_file.Data() << " is empty" << endl;exit(1);}

  char sfile[256];
  char stitle[256];
  int NPSD=0;
  while(true) {
    fpos=in.tellg();
    in.getline(str,1024);
    if(str[0] == '#') continue;
    in.seekg(fpos, ios::beg);

    in >> sfile >> stitle;
    if(!in.good()) break;

    if(NPSD>=MAX_PSD) {cout << "Error - Input Files exceed MAX_PSD " << MAX_PSD << endl;exit(1);}
    ifile[NPSD]=sfile;
    ltitle[NPSD]=stitle;
    cout << NPSD+1 << " " << ifile[NPSD] << " " << ltitle[NPSD] << endl;
    NPSD++;

    fpos=in.tellg();
    in.seekg(fpos+1, ios::beg);
  }
  in.close();

  // Read PSD files

  double freq[MAX_PSD][MAX_SIZE];
  double psd[MAX_PSD][MAX_SIZE];
  int lenght[MAX_PSD];

  for(int i=0;i<NPSD;i++) {
    ifstream in;
    in.open(ifile[i].Data(),ios::in);
    if (!in.good()) {cout << "Error Opening File : " << ifile[i].Data() << endl;exit(1);}
    lenght[i]=0;
    while (1) {
      double F,S;
      in >> F >> S;
      if (!in.good()) break;
      if(S!=0) {
        freq[i][lenght[i]]=F;
        psd[i][lenght[i]]=S;
        lenght[i]++;
      }
      if(lenght[i]>=MAX_SIZE) {cout << "Error size " << ifile[i].Data() << " > MAX_SIZE = " << MAX_SIZE << endl;exit(1);}
    }
    in.close();
  }

  // create plots

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

  Color_t colors[MAX_PSD] = { kBlue, kBlack, kRed, 6, 3, 8, 43, 7, 8, 4};
  int lstyle[MAX_PSD] = {2, 0, 0, 0, 9, 9, 9, 2, 2, 2};


  TCanvas *canvas = new TCanvas("Sensitivity", "psd", 300,40, 1000, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetLogx();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);

  TGraph* gr[MAX_PSD];
  for(int n=0;n<NPSD;n++) gr[n] = new TGraph(lenght[n],freq[n],psd[n]);
  for(int n=0;n<NPSD;n++) {
    gr[n]->SetLineWidth(3);
    gr[n]->SetMarkerColor(colors[n]);
    gr[n]->SetLineColor(colors[n]);
  }

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle("Sensitivity Curves");
  if(title!="") mg->SetTitle(title.Data());
  for(int n=0;n<NPSD;n++) mg->Add(gr[n]);
  mg->Paint("APL");

  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);

  if(min_freq>=0 && max_freq>min_freq) 
    mg->GetHistogram()->GetXaxis()->SetRangeUser(min_freq,max_freq);
  if(min_strain>=0 && max_freq>min_strain) 
    mg->GetHistogram()->GetYaxis()->SetRangeUser(min_strain,max_strain);

  mg->GetXaxis()->SetTitle(gr[0]->GetXaxis()->GetTitle());
  mg->GetXaxis()->SetLabelFont(42);
  mg->GetYaxis()->SetLabelFont(42);
  mg->GetXaxis()->SetTitleFont(42);
  mg->GetYaxis()->SetTitleFont(42);
  mg->GetXaxis()->SetTitleOffset(1.20);
  mg->GetYaxis()->SetTitleOffset(1.22);
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitle("Frequency (Hz)      ");
  mg->GetYaxis()->SetTitle("#frac{1}{#sqrt{Hz}}          ");

  mg->Draw("APL");

  // draw the legend

  TLegend* leg;  
  double hleg = 0.8-NPSD*0.05;
  leg = new TLegend(0.6120401,hleg,0.9615385,0.8721805,NULL,"brNDC");

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

  for(int n=0;n<NPSD;n++) {
    char legLabel[256];
    sprintf(legLabel,"%s",ltitle[n].Data());
    leg->AddEntry(gr[n],legLabel,"lp");
  }
  leg->Draw();

  // save plot

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
}

