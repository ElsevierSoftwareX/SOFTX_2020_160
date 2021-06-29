//
// read & display diagnostic histogram from job supercluster file 
// Author : Gabriele Vedovato

{
  #define JOB_SUPERCLUSTER_FILE "data/supercluster_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_MSA_30_job1.root"
  #define APPLY_CUT

  TString jname = JOB_SUPERCLUSTER_FILE;

  // open input job file
  TFile* jfile = new TFile(jname);
  if(jfile==NULL||!jfile->IsOpen())
    {cout << "Error : file " << jname << " not found" <<  endl;exit(1);}

  // read diagnostic histogram from supercluster job file
  TH2F* ph = (TH2F*)jfile->Get("hdng");
  if(ph==NULL) {cout<<"Error : diagnostic histogram not present"<<endl;exit(1);}

  // draw diagnostic histogram
  TCanvas* canvas;
  canvas = new TCanvas("hdng", "hdng", 35, 46, 600, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetGridx(true);
  canvas->SetGridy(true);
  canvas->SetLogz(true);
  canvas->SetFillColor(kWhite);
  canvas->SetLeftMargin(0.08);
  canvas->SetBottomMargin(0.13);
  canvas->SetBorderMode(0);

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");

  ph->SetTitle("subNetCut Diagnostic Histogram");
  ph->SetStats(kFALSE);

  ph->GetXaxis()->SetTitle("(Lm-Nm*Es)/(Em-Nm*Es)");
  ph->GetXaxis()->SetLabelFont(42);
  ph->GetXaxis()->SetLabelSize(0.03);
  ph->GetXaxis()->SetTitleFont(42);
  ph->GetXaxis()->SetTitleSize(0.03);
  ph->GetXaxis()->SetTitleOffset(1.5);
  ph->GetXaxis()->CenterTitle(true);
  ph->GetXaxis()->SetRangeUser(0,1);

  ph->GetYaxis()->SetTitle("STAT/Em");
  ph->GetYaxis()->SetLabelFont(42);
  ph->GetYaxis()->SetLabelSize(0.03);
  ph->GetYaxis()->SetTitleFont(42);
  ph->GetYaxis()->SetTitleSize(0.03);
  ph->GetYaxis()->SetTitleOffset(1.2);
  ph->GetYaxis()->CenterTitle(false);
  ph->GetYaxis()->SetRangeUser(0,1);

  ph->GetZaxis()->SetLabelFont(42);
  ph->GetZaxis()->SetLabelSize(0.025);
 
  ph->Draw("colz");

  // compute entries
  int NPIX=0;
  int npix=0;
  for (int i=0;i<=ph->GetNbinsX();i++) {
    for (int j=0;j<=ph->GetNbinsY();j++) {
      double X = ph->GetXaxis()->GetBinCenter(i)+ph->GetXaxis()->GetBinWidth(i)/2.;
      double Y = ph->GetYaxis()->GetBinCenter(i)+ph->GetYaxis()->GetBinWidth(i)/2.;
      double Z = ph->GetBinContent(i,j);
      if(X+Y>1) npix+=Z;
      NPIX+=Z;
    }
  }
  cout << "npix/NPIX : " << npix << " / " << NPIX << endl;
}
