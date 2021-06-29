
{

#define NBINS 3000
#define BINMIN 5.
#define BINMAX 305.
#define REBIN_MIN_EVTS_PER_BIN 10.0
  CWB::Toolbox TB;
  cout<<"Checking ENV "<<endl;
  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

 
#ifdef myBATCH 


 TCanvas *c0 = new TCanvas("c0","c0",0,0,800,600);
// c0->Range(-1.216392,-477.6306,508.8988,2814.609);
 c0->SetFillColor(0);
 c0->SetBorderMode(0);
 c0->SetBorderSize(2);
 c0->SetGridx();
 c0->SetGridy();
 c0->SetRightMargin(0.1154618);
 c0->SetTopMargin(0.07642487);
 c0->SetBottomMargin(0.1450777);
#else
cout << "Setting up the canvas and the pads..."<<endl;
TCanvas *canvas1 = new TCanvas("canvas1","C1",0,0,1200,800);
canvas1->Divide(2,3);
//TCanvas c1;
canvas1->GetPad(1)->SetGridx();
canvas1->GetPad(1)->SetGridy();
canvas1->GetPad(1)->SetLogy(kTRUE);
canvas1->GetPad(1)->SetLogx(kTRUE);
canvas1->GetPad(2)->SetLogx(kTRUE);
canvas1->GetPad(2)->SetLogy(kTRUE);
canvas1->GetPad(2)->SetGridx();
canvas1->GetPad(2)->SetGridy();
canvas1->GetPad(3)->SetGridx();
canvas1->GetPad(3)->SetGridy();
canvas1->GetPad(3)->SetLogx(kTRUE);
canvas1->GetPad(4)->SetGridx();
canvas1->GetPad(4)->SetGridy();
canvas1->GetPad(4)->SetLogx(kTRUE);
canvas1->GetPad(5)->SetLogy(kTRUE);
canvas1->GetPad(5)->SetGridx();
canvas1->GetPad(5)->SetGridy();
canvas1->GetPad(6)->SetLogy(kTRUE);
canvas1->GetPad(6)->SetGridx();
canvas1->GetPad(6)->SetGridy();

TCanvas *canvas2 = new TCanvas("canvas2","C2",0,0,800,600);

canvas2->SetGridx();
canvas2->SetGridy();
//canvas2->SetLogy(kTRUE);
//canvas2->SetLogx(kTRUE);

TCanvas *canvas3 = new TCanvas("canvas3","C3",0,0,800,600);

canvas3->SetGridx();
canvas3->SetGridy();

#endif


  //If not vetoes are defined, then the char string veto_not_vetoed is forced to be equal to ch2 
 if (strlen(veto_not_vetoed) == 0){sprintf(veto_not_vetoed,"%s",ch2);}
 
  gStyle->SetTitleFillColor(kWhite);
 // gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);
 // gStyle->SetOptStat(kFALSE);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

 //If not vetoes are defined, then the char string veto_not_vetoed is forced to be equal to ch2 
 if (strlen(veto_not_vetoed) == 0){sprintf(veto_not_vetoed,"%s",ch2);}
 

  TChain wave("waveburst");
  TChain live("liveTime");

  wave.Add(net_file_name);
  live.Add(liv_file_name);

 char cnet_file_name[1024];
 char cliv_file_name[1024];
  
  TString Compare_tree;
  if(gSystem->Getenv("CWB_COMPARE_TREE")==NULL) {
    cout << "Error : environment CWB_COMPARE_TREE is not defined!!!" << endl;exit(1);
  } else {
    Compare_tree=TString(gSystem->Getenv("CWB_COMPARE_TREE"));
  }
#ifdef myBATCH
 // char cdata_label[1024];
   TObjArray* token = TString(Compare_tree).Tokenize(TString("/"));
   TString comp_tag =((TObjString *)(token->At(token->GetEntries()-1)))->String();
   comp_tag.ReplaceAll("wave_","");
   comp_tag.ReplaceAll(".root","");
  //sprintf(cdata_label,((TObjString*)token->At(token->GetEntries()-1))->GetString().Data());
    sprintf(netdir,"%s/bkgcomp_vs_%s",netdir,comp_tag.Data());
   TB.mkDir(netdir,true);

#endif

	TString Compare_live = Compare_tree;
	Compare_live.ReplaceAll("wave_","live_");
    sprintf(cnet_file_name,"%s",Compare_tree.Data());
    sprintf(cnet_file_name,"%s",Compare_tree.Data());
    sprintf(cliv_file_name,"%s",Compare_live.Data());
    cout << cnet_file_name << endl;
    cout << cliv_file_name << endl;

    // Check if files exist
    TB.checkFile(cnet_file_name);
    TB.checkFile(cliv_file_name);
    TChain cnet("waveburst");
    TChain cliv("liveTime");
    cnet.Add(cnet_file_name);
    cliv.Add(cliv_file_name);

// TCut c2 ="netcc[0]>0.700000&&!(lag[3]==0&&slag[3]==0)&&neted[0]/ecor<0.400000&&penalty>0.600000&&(!veto_hveto_L1&&!veto_hveto_H1&&!veto_hveto_V1&&!veto_cat3_L1&&!veto_cat3_H1&&!veto_cat3_V1)"; 
   Int_t bkg_size = (Int_t)wave.GetEntries();
   cout << "bkg_size current dir : " << bkg_size << endl;
   Int_t bkg_size2 = (Int_t)cnet.GetEntries();
   cout << "bkg_size compare dir : " << bkg_size2 << endl;
 


 // float xlag;
 // Int_t xrun, rUn;
 // double xlive, sTARt, sTOp, gps;
  double liveTot = 0.;
  double cliveTot = 0.;
 // double Live;
/* 
  wavearray<double> Trun(500000); Trun = 0.;
  wavearray<double>* Wlag = new wavearray<double>[nIFO+1];
  wavearray<double>* Wslag = new wavearray<double>[nIFO+1];
  wavearray<double> Tdlag;
  wavearray<double> Tlag;
 // wavearray<double> Rlag;
 // wavearray<double> Elag;
 // wavearray<double> Olag;

  cout << "Start CWB::Toolbox::getLiveTime of current ntuple : " << endl;
  liveTot=TB.getLiveTime(nIFO,live,Trun,Wlag,Wslag,Tlag,Tdlag,cwb_lag_number,cwb_slag_number);cout<<liveTot<<" s"<<endl;
   cout << "Start CWB::Toolbox::getLiveTime of current ntuple : " << endl;
  cliveTot=TB.getLiveTime(nIFO,cliv,Trun,Wlag,Wslag,Tlag,Tdlag,cwb_lag_number,cwb_slag_number);cout<<cliveTot<<" s"<<endl;
*/

 CWB::CBCTool cbcTool;
 cout << "Start CWB::CBCTool::getLiveTime2 of current ntuple : " << endl;
 liveTot=cbcTool.getLiveTime2(live);  //BEWARE!! getLiveTime2 returns bkg+zerolag time (usually a small error)
  cout << "Start CWB::CBCTool::getLiveTime2 of current ntuple : " << endl;
 cliveTot=cbcTool.getLiveTime2(cliv); //BEWARE!! getLiveTime2 returns bkg+zerolag time (usually a small error)
 
  double K = liveTot/cliveTot;
  cout << "K (ratio of live times) : "<< K <<endl;
  //int run;
   //UChar_t v;
   //float rho[2];
   //trbkg->SetBranchAddress("rho",rho);
   
   TH1D* curh = new TH1D("curh","",NBINS,BINMIN,BINMAX);
   TH1D* comph = new TH1D("comph","",NBINS,BINMIN,BINMAX);
   TH1D* Sh = new TH1D("Sh","Significance vs rho bin",NBINS,BINMIN,BINMAX);
   TH1D* Sh2 = new TH1D("Sh2","Distribution of significance",200,-10.,10.);

//cout <<"1"<<endl;
  char st[256];
  TString mycut = veto_not_vetoed;
  sprintf(st,"rho[%d]>>hcomp(%d,%f,%f)",pp_irho,NBINS,BINMIN,BINMAX);
  cout<<"Rho min: "<<T_cut<<endl;
  mycut += "&& rho["+mycut.Itoa(pp_irho,10)+"]>";
  mycut += T_cut;
  
 //  cnet.Draw(st,veto_not_vetoed,"goff");
 //TCut cut=c_veto_not_vetoed;

  // cnet.Draw(st,cut,"goff");
  // cout<<"2bis"<<endl;
 // Long64_t entriesx = (Long64_t)cnet.GetSelectedRows(); 
 //  cout<< entriesx<<endl;
 // Long64_t entriesx = cnet.Draw(st,veto_not_vetoed,"goff");
  Long64_t entriesx = cnet.Draw(st,mycut,"goff");
 cout<<"Comparison bkg events above threshold: "<< entriesx <<endl; 
 double* rhox = cnet.GetV1();
//cout<<"3"<<endl; 
   comph = (TH1D*)gDirectory->Get("hcomp");
  //  comph = (TH1D*) trbkg2->GetHistogram();
  sprintf(st,"rho[%d]>>hcur(%d,%f,%f)",pp_irho,NBINS,BINMIN,BINMAX);
  //Long64_t entriesy = wave.Draw(st,veto_not_vetoed,"goff");
  Long64_t entriesy = wave.Draw(st,mycut,"goff");
cout<<"Current bkg events above threshold: "<< entriesy <<endl;
// Long64_t entriesy = wave.Draw(st,cut,"goff");
  double* rhoy = wave.GetV1();
  double RHOMAX=TMath::Max(TMath::MaxElement(entriesy,rhoy),TMath::MaxElement(entriesx,rhox)); 
  cout<<"Rho max: "<<RHOMAX<<endl;
   // curh =(TH1D*) trbkg->GetHistogram(); 
//  TH1D *comph = (TH1D*)gDirectory->Get("hnew");
  curh = (TH1D*)gDirectory->Get("hcur");
//  curh->Scale(1/OLIVE);

//for(int i=0;i<300;i++){comph->SetBinContent(i+1,TMath::Nint(nrate[i]*NLIVE));}
//comph->Scale(1/NLIVE);
comph->GetXaxis()->SetRangeUser(BINMIN,BINMAX);

cout<<"Un-binned KS test"<<endl;
TMath::KolmogorovTest(entriesx,rhox,entriesy,rhoy,"D");
cout<<"Binned KS test (no rescaling)"<<endl;
curh->KolmogorovTest(comph,"D");
curh->AndersonDarlingTest(comph,"D");


curh->GetXaxis()->SetTitle("#rho");
curh->GetXaxis()->CenterTitle(true);
curh->GetYaxis()->SetTitle("#events");
curh->GetYaxis()->CenterTitle(true);
curh->SetFillColor(kBlue);
curh->SetTitle("Distribution of rho");

TH1D* comphr=(TH1D*)comph->Clone("comphr");
comphr->Scale(K);
double maxy = TMath::Max(comphr->GetMaximum(),curh->GetMaximum());
curh->GetYaxis()->SetRangeUser(0.1,1.2*maxy);
curh->GetXaxis()->SetRangeUser(T_cut,1.2*RHOMAX);
//curhr->GetXaxis()->SetTitle("RHO");
//curhr->GetYaxis()->SetTitle("#events");
comphr->SetLineColor(kRed);
comphr->SetFillColor(kRed);
comphr->SetFillStyle(3001);


TLegend* leg = new TLegend(0.65,0.7,0.99,0.96);

  //leg->AddEntry("");
leg->SetTextSize(0.03);
char entry[256];
sprintf(entry,"Rescaling factor K = %f", K);
 leg->SetHeader(entry);
sprintf(entry,"Current BKG");
leg->AddEntry(curh,entry,"f");
sprintf(entry,"Rescaled Comparison BKG");
leg->AddEntry(comphr,entry,"f");


#ifdef myBATCH
c0->cd();
c0->Clear();
c0->SetLogy(kTRUE);
c0->SetLogx(kTRUE);
curh->Draw();
comphr->Draw("same");
leg->Draw();
char fname[1024];
sprintf(fname,"%s/Rho_Distribution.png",netdir);
c0->Update();
c0->SaveAs(fname);
#else

cout <<"Plotting pad 1"<<endl;
canvas1->cd(1);
curh->Draw();
comphr->Draw("same");
    leg->Draw();

#endif

const Int_t nq = entriesy/REBIN_MIN_EVTS_PER_BIN;cout<<"nq ="<<nq<<endl;
Double_t* xq = new Double_t[nq];  // position where to compute the quantiles in [0,1]
Double_t* yq = new Double_t[nq]; 
Double_t* yq2 = new Double_t[nq+1]; 
for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
curh->GetQuantiles(nq,yq,xq);

for(int i=1;i<nq+1;i++){yq2[i]=yq[i-1];}
yq2[0]= BINMIN;
yq2[nq]= RHOMAX*1.2;
TH1D *hnew2 = new TH1D("hnew","rebinned",nq,yq2);
TH1D *hnew3 = new TH1D("hnew3","rebinned",nq,yq2);
TH1D* Sh3 = new TH1D("Sh3","Significance vs rho bin",nq,yq2);
TH1D* Sh4 = new TH1D("Sh4","Distribution of significance",200,-10.,10.);
TH1D* Shn = new TH1D("Shn","Normal Distribution",200,-10.,10.);
TRandom3 P;
P.SetSeed(1234.);
for(int i=0;i<entriesy;i++){hnew2->Fill(rhoy[i]);}
for(int i=0;i<entriesx;i++){hnew3->Fill(rhox[i]);}

hnew2->GetXaxis()->SetTitle("#rho");
hnew2->GetXaxis()->CenterTitle(true);
hnew2->GetYaxis()->SetTitle("#events");
hnew2->GetYaxis()->CenterTitle(true);
hnew2->SetFillColor(kBlue);
hnew2->SetTitle("Distribution of rho");
hnew3->Scale(liveTot/cliveTot);
hnew3->SetLineColor(kRed);
hnew3->SetFillColor(kRed);
hnew3->SetFillStyle(3001);


#ifdef myBATCH
c0->cd();
c0->Clear();
c0->SetLogy(kTRUE);
c0->SetLogx(kTRUE);
hnew2->Draw();
hnew3->Draw("SAME");
leg->Draw();
sprintf(fname,"%s/Rho_Distribution_rebin.png",netdir);
c0->Update();
c0->SaveAs(fname);
#else
cout<<"Plotting pad 2"<<endl;
canvas1->cd(2);

hnew2->Draw();
hnew3->Draw("SAME");
#endif
cout<<"After Rebinning"<<endl;
cout <<"Anderson-Darling test (no rescaling)"<<endl;
hnew2->AndersonDarlingTest(hnew3,"D");




for(int i=0;i<NBINS+1;i++){
//    cout << i <<" "<<curh->GetBinContent(i+1)<<" "<<comph->GetBinContent(i+1)<<" ";
	if((curh->GetBinContent(i+1)>0) || (comph->GetBinContent(i+1)>0)){
        
		Sh->SetBinContent(i+1,(curh->GetBinContent(i+1)-comph->GetBinContent(i+1)*K)/TMath::Sqrt(pow(curh->GetBinError(i+1),2)+pow(comph->GetBinError(i+1)*K,2)));
	        Sh2->Fill(Sh->GetBinContent(i+1)); 
	} 
	else{Sh->SetBinContent(i+1,0.0);}
  //   cout<<Sh->GetBinContent(i+1)<<" - ";
}
Sh->GetXaxis()->SetTitle("#rho");
Sh->GetXaxis()->CenterTitle(true);
Sh->GetYaxis()->SetTitle("Observed normalized significances");
Sh->GetYaxis()->CenterTitle(true);


#ifdef myBATCH
c0->cd();
c0->Clear();
c0->SetLogy(kFALSE);
c0->SetLogx(kTRUE);
Sh->GetXaxis()->SetRangeUser(T_cut,1.2*RHOMAX);
Sh->Draw();
sprintf(fname,"%s/Significance_rho.png",netdir);
c0->Update();
c0->SaveAs(fname);
#else
cout<<"Plotting pad 3"<<endl;
canvas1->cd(3);
Sh->Draw();
#endif


double w_tmp;
Double_t* w = new Double_t[nq+2];
int count =0;
for(int i=0;i<nq+2;i++){
//     cout << i <<" "<<hnew2->GetBinContent(i+1)<<" "<<hnew3->GetBinContent(i+1)<<" ";
	if((hnew2->GetBinContent(i+1)>0) ||(hnew3->GetBinContent(i+1)>0)){
		w_tmp = (hnew2->GetBinContent(i+1)-hnew3->GetBinContent(i+1))/TMath::Sqrt(pow(hnew2->GetBinError(i+1),2)+pow(hnew3->GetBinError(i+1),2));
		Sh3->SetBinContent(i+1,w_tmp);
	    Sh4->Fill(Sh3->GetBinContent(i+1));
		w[i]=w_tmp;
		count++;  	
	} 
	else{Sh3->SetBinContent(i+1,0.0);w[i]=0.0;}
  //  cout<<w[i]<<" - ";
}
cout<<endl;
Sh3->GetXaxis()->SetTitle("#rho");
Sh3->GetXaxis()->CenterTitle(true);
Sh3->GetYaxis()->SetTitle("Observed normalized significances");
Sh3->GetYaxis()->CenterTitle(true);
#ifdef myBATCH
c0->cd();
c0->Clear();
c0->SetLogy(kFALSE);
c0->SetLogx(kTRUE);
Sh3->GetXaxis()->SetRangeUser(T_cut,1.2*RHOMAX);
Sh3->Draw();
sprintf(fname,"%s/Significance_rebinrho.png",netdir);
c0->Update();
c0->SaveAs(fname);
#else
cout<<"Plotting pad 4"<<endl;
canvas1->cd(4);
Sh3->Draw();
#endif

Sh2->GetXaxis()->SetTitle("Observed normalized significances");
Sh2->GetXaxis()->CenterTitle(true);
Sh2->GetYaxis()->SetTitle("#bins");
Sh2->GetYaxis()->CenterTitle(true);
Sh2->SetFillColor(kBlue);
Sh2->SetFillStyle(3001);

#ifdef myBATCH
c0->cd();
c0->Clear();
c0->SetLogy(kTRUE);
c0->SetLogx(kFALSE);
Sh2->Draw();
sprintf(fname,"%s/Significance.png",netdir);
c0->Update();
c0->SaveAs(fname);
#else

cout<<"Plotting pad 5"<<endl;
canvas1->cd(5);
Sh2->Draw();
#endif

Sh4->GetXaxis()->SetTitle("Observed normalized significances");
Sh4->GetXaxis()->CenterTitle(true);
Sh4->GetYaxis()->SetTitle("#bins");
Sh4->GetYaxis()->CenterTitle(true);
Sh4->SetFillColor(kBlue);
Sh4->SetFillStyle(3001);

#ifdef myBATCH
c0->cd();
c0->Clear();
c0->SetLogy(kTRUE);
c0->SetLogx(kFALSE);
Sh4->Draw();
sprintf(fname,"%s/Significance_rebin.png",netdir);
c0->Update();
c0->SaveAs(fname);
#else
cout<<"Plotting pad 6"<<endl;
canvas1->cd(6);
Sh4->Draw();
#endif

for(int i=0;i<(int)Sh4->GetEntries();i++){Shn->Fill(P.Gaus(0,1));}
cout <<"Anderson-Darling 2-sample test significance distribution vs Normal (no rescaling)"<<endl;
Sh4->AndersonDarlingTest(Shn,"D");

ROOT::Math::GoFTest* goftest_1 = new ROOT::Math::GoFTest(count, w, ROOT::Math::GoFTest::kGaussian); 
Double_t pvalueAD_1 = goftest_1-> AndersonDarlingTest();
cout <<"Anderson-Darling 1-sample test significance distribution vs Normal (no rescaling) : "<<pvalueAD_1<<endl;

cout<<"Starting Q-Q plot ..."<<endl;
TGraphQQ* gr;
if (entriesy>=entriesx){

	gr = new TGraphQQ(entriesy,rhoy,entriesx,rhox);
	gr->GetHistogram()->GetXaxis()->SetTitle("Quantiles current histogram");
	gr->GetHistogram()->GetYaxis()->SetTitle("Quantiles compare histogram");
}
else{
        gr = new TGraphQQ(entriesx,rhox,entriesy,rhoy);
        gr->GetHistogram()->GetXaxis()->SetTitle("Quantiles compare histogram");
        gr->GetHistogram()->GetYaxis()->SetTitle("Quantiles current histogram");
}


 gr->SetMarkerColor(kRed);
 gr->GetHistogram()->GetYaxis()->CenterTitle(true);
 gr->GetHistogram()->GetXaxis()->CenterTitle(true);
 // gr->SetLineColor(kRed);
  gr->SetMarkerSize(1.0);
 gr->SetMarkerStyle(20);
 gr->SetTitle("Quantile-Quantile Plot");
 TLine *line = new TLine(BINMIN, BINMIN,BINMAX,BINMAX);
line->SetLineColor(kGreen);
line->SetLineWidth(2);
 line->SetLineStyle(2);

#ifdef myBATCH
c0->cd();
c0->SetLogy(kFALSE);
c0->SetLogx(kFALSE);
cout<<"Drawing Q-Q plot ..."<<endl;
gr->Draw("ALP");
cout<<"Finished Q-Q plot ..."<<endl;

line->Draw();
sprintf(fname,"%s/Q-Q_Plot.png",netdir);
c0->Update();
c0->SaveAs(fname);
#else
canvas2->cd();
gr->Draw("ALP");
line->Draw();
#endif

TF1* g1 = new TF1("g1","gaus",-10.,10.);
g1->SetParameter(0,count/sqrt(TMath::TwoPi()));
g1->SetParameter(1,0.);
g1->SetParameter(2,1.);

TGraphQQ* gr2 = new TGraphQQ(nq,w,g1);
 gr2->SetMarkerColor(kRed);
 gr2->GetHistogram()->GetYaxis()->CenterTitle(true);
 gr2->GetHistogram()->GetXaxis()->CenterTitle(true);
 // gr->SetLineColor(kRed);
  gr2->SetMarkerSize(1.0);
 gr2->SetMarkerStyle(20);
 gr2->SetTitle("Quantile-Quantile Plot");
 double binmin = TMath::Max(gr2->GetX()[0],gr2->GetY()[0]);
double binmax = TMath::Min(gr2->GetX()[nq-1],gr2->GetY()[nq-1]);
 TLine *line2 = new TLine(binmin,binmin,binmax,binmax);
line2->SetLineColor(kGreen);
line2->SetLineWidth(2);
 line2->SetLineStyle(2);
 cout << "Normalized significance distribution vs Gauss(0,1)"<<endl;
 cout << "Correlation Factor : "<<gr2->GetCorrelationFactor()<<endl;
 cout  << "Covariance : " <<gr2->GetCovariance()<<endl;;

#ifdef myBATCH
c0->cd();
c0->Clear();
c0->SetLogy(kFALSE);
c0->SetLogx(kFALSE);
gr2->Draw("ALP");
line2->Draw();
sprintf(fname,"%s/Q-Q_Plot_significance.png",netdir);
c0->Update();
c0->SaveAs(fname);

char exec[1024];
sprintf(exec,"ln -s %s/%s/index.html %s/index.html",work_dir, pp_dir, netdir);
gSystem->Exec(exec);
sprintf(exec,"ln -s %s/%s/header.html %s/header.html",work_dir, pp_dir, netdir);
gSystem->Exec(exec);
sprintf(exec,"cp $CWB_MACROS/../postproduction/burst/compare_bkg.html %s/body.html",netdir);
gSystem->Exec(exec);



exit(0);
#else
canvas3->cd();
gr2->Draw("ALP");
line2->Draw();
#endif

 

//exit(0);
}
