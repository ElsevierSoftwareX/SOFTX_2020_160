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


// used by the cwb_mkeff command
{

  #define NINJ_MAX 50
  #define NTYPE_MAX 20 
  #define NMDC_MAX 64

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.2,"Y");
//  gStyle->SetLabelOffset(0.014,"X");
//  gStyle->SetLabelOffset(0.010,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
//  gStyle->SetLabelSize(0.03,"X");
//  gStyle->SetLabelSize(0.03,"Y");

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
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(kFALSE);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  // read injection file types
  char    imdc_set[NMDC_MAX][128];    // injection set
  size_t  imdc_type[NMDC_MAX];        // injection type
  char    imdc_name[NMDC_MAX][128];   // injection name
  double  imdc_fcentral[NMDC_MAX];    // injection central frequencies
  double  imdc_fbandwidth[NMDC_MAX];  // injection bandwidth frequencies
  size_t  imdc_index[NMDC_MAX];       // type reference array
  size_t  imdc_iset[NMDC_MAX];        // injection set index

  int ninj=ReadInjType(mdc_inj_file,NMDC_MAX,imdc_set,imdc_type,
                       imdc_name,imdc_fcentral,imdc_fbandwidth);
  if(ninj==0) {cout << "Error - no injection - terminated" << endl;exit(1);}

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
  for (int iset=0;iset<nset;iset++) cout << iset << " " << imdc_set_name[iset].Data() << endl;

  char etitle[256];
  char ofile[256];
  TCanvas* canvas[NTYPE_MAX];
  int ecount[NINJ_MAX];
  TString piumeno[NINJ_MAX];
  float chi2[NINJ_MAX], err[NINJ_MAX], par1[NINJ_MAX], par2[NINJ_MAX], par3[NINJ_MAX];
  double ehrss10[NINJ_MAX], ehrss50[NINJ_MAX], ehrss90[NINJ_MAX];
  double hrss50_bis[NINJ_MAX];
  TString ewaveform[NINJ_MAX];
   
  TF1* fFit[NINJ_MAX];
  double hrss50[NTYPE_MAX][NINJ_MAX], hrss90[NTYPE_MAX][NINJ_MAX], hrss10[NTYPE_MAX][NINJ_MAX];

  double inf = simulation==2 ? log10(factors[0]) : -25;
  double sup = simulation==2 ? log10(factors[nfactor-1]) : -18.5;

  if(simulation==1 && pp_factor2distance) {
    inf = log10(pp_factor2distance/factors[nfactor-1]);
    sup = log10(pp_factor2distance/factors[0]);
  }

  int k=0;
  for (int iset=0;iset<nset;iset++) {

    char file[256];
    sprintf(file,"%s/fit_parameters_%s.txt",netdir, imdc_set_name[iset].Data());
    cout << file << endl;

    ifstream in2;
    in2.open(file,ios::in);
    if (!in2.good()) {cout << "Error Opening File : " << file << endl;exit(1);}

    for (int j=0; j<NINJ_MAX; j++) {
      hrss50_bis[j]=0;
      hrss10[iset][j]=0;
      hrss50[iset][j]=0;
      hrss90[iset][j]=0;
      ecount[j]=0;
      ewaveform[j]="";
    }

    for (int l=0; l<NINJ_MAX; l++) {

      in2 >> ecount[k] >> chi2[k] >> hrss50[iset][k] >> piumeno[k] 
          >> err[k] >> par1[k] >> par2[k] >> par3[k] >> ewaveform[k];
      if (!in2.good()) break;
      cout << ewaveform[k].Data() << endl;
      double par0=TMath::Log10(hrss50[iset][k]);
      fFit[k] = new TF1("logNfit",logNfit,pow(10.0,inf),pow(10.0,sup),5);
      fFit[k]->SetNpx(100000);
      fFit[k]->SetParameters(par0,par1[k],par2[k],par3[k],pp_factor2distance);
      hrss10[iset][k]=fFit[k]->GetX(.1,pow(10.0,inf),pow(10.0,sup));
      hrss50_bis[k]  =fFit[k]->GetX(.5,pow(10.0,inf),pow(10.0,sup));
      hrss90[iset][k]=fFit[k]->GetX(.9,pow(10.0,inf),pow(10.0,sup));
      if(fFit[k]->Eval(hrss90[iset][k])<0.89) hrss90[iset][k]=1e-10;

      par0=TMath::Log10(hrss50[iset][k]+err[k]);
      fFit[k]->SetParameters(par0,par1[k],par2[k],par3[k],pp_factor2distance);
      ehrss10[k]=fabs(fFit[k]->GetX(.1,pow(10.0,inf),pow(10.0,sup))-hrss10[iset][k]);
      ehrss50[k]=fabs(fFit[k]->GetX(.5,pow(10.0,inf),pow(10.0,sup))-hrss50[iset][k]);
      ehrss90[k]=fabs(fFit[k]->GetX(.9,pow(10.0,inf),pow(10.0,sup))-hrss90[iset][k]);

      cout << hrss10[iset][k] << " " << hrss50[iset][k] << " " << hrss90[iset][k] << endl;

      k++;
    }
  }

  TGraphErrors* gr10[NTYPE_MAX];
  TGraphErrors* gr50[NTYPE_MAX];
  TGraphErrors* gr90[NTYPE_MAX];
  TMultiGraph*  mg[NTYPE_MAX];
  TLegend*      legend[NTYPE_MAX];

  for(int iset=0;iset<nset;iset++) {

    sprintf(etitle,"%s",imdc_set_name[iset].Data());

    mg[iset]     = new TMultiGraph();
    legend[iset] = new TLegend(0.732,0.125,0.956,0.354,NULL,"brNDC");


    gr10[iset] = new TGraphErrors();
    gr10[iset]->SetLineColor(1);
    gr10[iset]->SetLineWidth(1);
    gr10[iset]->SetMarkerColor(1);
    gr10[iset]->SetMarkerStyle(20);
    gr10[iset]->SetLineStyle(7);

    gr50[iset] = new TGraphErrors();
    gr50[iset]->SetLineColor(2);
    gr50[iset]->SetLineWidth(1);
    gr50[iset]->SetMarkerColor(2);
    gr50[iset]->SetMarkerStyle(20);
    gr50[iset]->SetLineStyle(7);

    gr90[iset] = new TGraphErrors();
    gr90[iset]->SetLineColor(kBlue);
    gr90[iset]->SetLineWidth(1);
    gr90[iset]->SetMarkerColor(kBlue);
    gr90[iset]->SetMarkerStyle(20);
    gr90[iset]->SetLineStyle(7);

    int npoint=0; 
    for(int i=0; i<ninj; i++) {
      if(imdc_iset[i]!=iset) continue;  // skip unwanted injection types
      if(pp_show_eff_fit_curve) { 	// set hrss@10/50/90 from efficiency fits
        gr10[iset]->SetPoint(npoint,imdc_fcentral[i],hrss10[iset][i]);
        gr50[iset]->SetPoint(npoint,imdc_fcentral[i],hrss50[iset][i]);
        gr90[iset]->SetPoint(npoint,imdc_fcentral[i],hrss90[iset][i]);
      } else {
        gr10[iset]->SetPoint(npoint,imdc_fcentral[i],0.);
        gr50[iset]->SetPoint(npoint,imdc_fcentral[i],0.);
        gr90[iset]->SetPoint(npoint,imdc_fcentral[i],0.);
      }
      npoint++;   
    }

    mg[iset]->Add(gr10[iset]);
    mg[iset]->Add(gr50[iset]);
    mg[iset]->Add(gr90[iset]);

    char namecanvas[8];
    sprintf(namecanvas,"c%i",iset);
    canvas[iset] = new TCanvas(namecanvas,namecanvas,125+iset*10,82,976,576);
    canvas[iset]->Clear();
    canvas[iset]->ToggleEventStatus();
    canvas[iset]->SetLogy();
#ifdef CWB_MKEFF_LINX 
    canvas[iset]->SetLogx(false);
#else
    canvas[iset]->SetLogx(true);
#endif
    canvas[iset]->SetGridx();
    canvas[iset]->SetGridy();
    canvas[iset]->SetFillColor(kWhite); 
    canvas[iset]->cd();

    mg[iset]->SetTitle(etitle);
    mg[iset]->Paint("alp");
    mg[iset]->GetHistogram()->GetXaxis()->SetTitle("Frequency (Hz)");
    mg[iset]->GetHistogram()->GetXaxis()->CenterTitle(true);
    mg[iset]->GetHistogram()->GetXaxis()->SetLabelFont(42);
    mg[iset]->GetHistogram()->GetXaxis()->SetTitleFont(42);
    mg[iset]->GetHistogram()->GetYaxis()->SetLabelFont(42);
    mg[iset]->GetHistogram()->GetYaxis()->SetTitleFont(42);
    if(simulation==1 && pp_factor2distance) {
      mg[iset]->GetHistogram()->GetYaxis()->SetTitle("distance (Kpc)");
      mg[iset]->GetHistogram()->GetYaxis()->SetRangeUser(pp_factor2distance/factors[nfactor-1],pp_factor2distance/factors[0]);
    } else if(simulation==2) {
      mg[iset]->GetHistogram()->GetYaxis()->SetTitle("snr");
      mg[iset]->GetHistogram()->GetYaxis()->SetRangeUser(factors[0],factors[nfactor-1]);
    } else {
      mg[iset]->GetHistogram()->GetYaxis()->SetTitle("hrss (1/Hz^{-1/2})");
      mg[iset]->GetHistogram()->GetYaxis()->SetRangeUser(pp_hrss_min,pp_hrss_max);
    }

    mg[iset]->Draw("alp");

    legend[iset]->SetBorderSize(1);
    legend[iset]->SetTextAlign(22);
    legend[iset]->SetTextFont(12);
    legend[iset]->SetLineColor(1);
    legend[iset]->SetLineStyle(1);
    legend[iset]->SetLineWidth(1);
    legend[iset]->SetFillColor(0);
    legend[iset]->SetFillStyle(1001);
    legend[iset]->SetTextSize(0.04);
    legend[iset]->SetLineColor(kBlack);
    legend[iset]->SetFillColor(kWhite);

    if(simulation==1 && pp_factor2distance) {
      legend[iset]->AddEntry(gr90[iset],"dstance (Kpc) 90%","lp");
      legend[iset]->AddEntry(gr50[iset],"dstance (Kpc) 50%","lp");
      legend[iset]->AddEntry(gr10[iset],"dstance (Kpc) 10%","lp");
    } else if(simulation==2) {
      legend[iset]->AddEntry(gr90[iset],"snr 90%","lp");
      legend[iset]->AddEntry(gr50[iset],"snr 50%","lp");
      legend[iset]->AddEntry(gr10[iset],"snr 10%","lp");
    } else {
      legend[iset]->AddEntry(gr90[iset],"Hrss 90%","lp");
      legend[iset]->AddEntry(gr50[iset],"Hrss 50%","lp");
      legend[iset]->AddEntry(gr10[iset],"Hrss 10%","lp");
    }
    legend[iset]->Draw();

    sprintf(ofile,"%s/eff_freq_%s.gif",netdir,imdc_set_name[iset].Data());
    cout << ofile << endl;
    canvas[iset]->SaveAs(ofile);
  }

  exit(0);
}
