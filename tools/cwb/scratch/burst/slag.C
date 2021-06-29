{

  cout<<"slag starts..."<<endl;

  if(nIFO>3) {cout << "slag.C : Skipping plot production - wave tree should come from a 2- or 3-fold network : " << net_file_name << endl;exit(0);} 
  
  //CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(kFALSE);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  

  // ----------------------------------------------------------
  // canvas
  // ----------------------------------------------------------
  TCanvas *c1 = new TCanvas("Slag", "Slag",32,55,750,502);
  c1->Range(-0.75,-1.23433,6.75,8.17182);
  c1->SetBorderSize(1);
  c1->SetFillColor(0);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetFrameFillColor(0);
  gStyle->SetPalette(1,0);

  TChain wave("waveburst");
  TChain live("liveTime");

  wave.Add(net_file_name);
  live.Add(liv_file_name);

  // check if slag is presents in liv tree
  TBranch* branch;
  bool slagFound=false;
  TIter mynext(wave.GetListOfBranches());
  while ((branch=(TBranch*)mynext())) {
    if (TString("slag").CompareTo(branch->GetName())==0) slagFound=true;
  }
  mynext.Reset();
  if(!slagFound) {cout << "slag.C : Error - wave tree does not contain slag leaf : " << net_file_name << endl;exit(1);}
  
  char fname[1024];
char cut[512];
char cut2[256];

sprintf(cut,"rho[1]>%f",T_cut);
//test if vetoes are applied
TString test_veto = veto_not_vetoed;
if(!test_veto.IsNull())sprintf(cut,"%s && %s",cut,veto_not_vetoed);
sprintf(cut2,"!(lag[%d]==0&&slag[%d]==0)",nIFO,nIFO);
   long int num = wave.GetEntries();
   wave.SetEstimate(num);
 //  TString sel("slag[");
   TString sel("slag[1]:slag[2]");
  // sel+= TString::Itoa(nIFO-2,10);
  // sel+="]:slag[";
  // sel+=Itoa(nIFO-1,10);
  // sel+="]";
  if(nIFO==2) {sel.ReplaceAll("1","0");sel.ReplaceAll("2","1");}
  // TString sel="slag[1]:slag[2]";
  num=wave.Draw(sel,cut,"goff");
  const int slag_entries = num;
  cout<<"Number of Entries: "<<num<<endl;
  double* slag1 = new double[slag_entries];
  slag1 = wave.GetV1();
  double* slag2 = new double[slag_entries]; 
  slag2 = wave.GetV2(); 
 char mytitle[256];
 double SlagMax = wave.GetMaximum("slag")+segLen/2.;
 double SlagMin = wave.GetMinimum("slag")-segLen/2.;
 int NSlag = TMath::FloorNint((SlagMax-SlagMin)/segLen); 
 cout << "SLAG MAX : "<< wave.GetMaximum("slag")<< " s  SLAG MIN : "<< wave.GetMinimum("slag")<< " s  #SLAGS : "<<NSlag-1<<endl;

 if(NSlag==1){cout <<"Just one slag....Skipping further execution!"<<endl; exit(0);}
sprintf(mytitle,"FAR distribution over slags (post cat3 & rho>%f)",T_cut);
	TH2F* Slag = new TH2F("SLAG",mytitle,NSlag,SlagMin/86400.,SlagMax/86400.,NSlag,SlagMin/86400.,SlagMax/86400.);
	Slag->GetXaxis()->SetTitle("slag[1] shift [day]");
// 	Slag->GetXaxis()->SetNdivisions(10,kFALSE);
	Slag->GetYaxis()->SetTitle("slag[2] shift [day]");
//	Slag->GetYaxis()->SetNdivisions(10,kFALSE);
	Slag->SetStats(kFALSE);

	TH2F* lSlag = new TH2F("LSLAG","Live time distribution over slags",NSlag,SlagMin/86400.,SlagMax/86400.,NSlag,SlagMin/86400.,SlagMax/86400.);
	for(long int l=0; l<num; l++){
		Slag->Fill(slag1[l]/86400.,slag2[l]/86400.);
//	if(l%1000==0) {cout << scientific << (double)l <<" - ";}

	}

    float  xlag[NIFO_MAX+1];
    float  xslag[NIFO_MAX+1];
    double xlive;
    live.SetBranchAddress("slag",xslag);
    live.SetBranchAddress("lag",xlag);
    live.SetBranchAddress("live",&xlive);
    live.SetBranchStatus("gps",false);

    int ntrg = live.GetEntries();
    for(int i=0; i<ntrg; i++) {
					live.GetEntry(i);
                    lSlag->Fill(xslag[0]/86400.,xslag[1]/86400.,xlive);

        }

	Slag->Divide(lSlag);
	c1->Clear();
	Slag->Draw("colz");
 	sprintf(fname,"%s/slag.png",netdir);
  	c1->Update(); c1->SaveAs(fname);
    
    TH1D* lSlag2 = lSlag->ProjectionY();
    lSlag2->GetXaxis()->SetTitle("slag shift [day]");
    //lSlag2->GetXaxis()->SetNdivisions(10,kFALSE);
    lSlag2->GetYaxis()->SetTitle("Live time [day]");
    //lSlag2->GetYaxis()->SetNdivisions(10,kFALSE);
    lSlag2->SetMarkerStyle(20);
    lSlag2->Scale(1./86400);  //Scaling to get the live time in days
    sprintf(mytitle,"Live time (%d lags) over slag shift",lagSize);
    lSlag2->SetTitle(mytitle);
    lSlag2->Draw("P");
    sprintf(fname,"%s/Live_vs_slagtime.png",netdir);
    c1->Update(); c1->SaveAs(fname);
    

	TAxis *xaxis = Slag->GetXaxis();
	TAxis *yaxis = Slag->GetYaxis();
	int mycount = 0;
        double FAR[50000];
        double LSlag[50000];
        double eFAR[50000];
        double D[50000];
        double Dx[50000];
        double Dy[50000];
        for(int i =1 ;i<NSlag+1;i++){
                for(int j =1 ;j<NSlag+1;j++){
                        if(Slag->GetBinContent(i,j)>=1.E-15){
				FAR[mycount] = Slag->GetBinContent(i,j)*86400.*365;
				LSlag[mycount] = lSlag->GetBinContent(i,j);	
				LSlag[mycount] = lSlag->GetBinContent(i,j);	
				Dx[mycount] = xaxis->GetBinCenter(i);
				Dy[mycount] = yaxis->GetBinCenter(j);
				D[mycount] = TMath::Sqrt(pow(Dx[mycount],2)+pow(Dy[mycount],2));
				//D[mycount] = yaxis->GetBinCenter(j);
				eFAR[mycount] = Slag->GetBinError(i,j)/lSlag->GetBinError(i,j)*86400.*365;
				mycount++;
				if(i*j%100==0){cout << "slag[1] :"<<xaxis->GetBinCenter(i)<<" slag[2] :"<<yaxis->GetBinCenter(j) <<" FAR :"<<FAR[mycount-1] << "+/-" <<eFAR[mycount-1] << " " <<endl;}
			}

		}		
        }

	
 	cout << "Number of slags : " << mycount << endl;
	double liveTot = lSlag->GetSumOfWeights();
	double liveMax = lSlag->GetMaximum();
	cout << "Total BKG live time : " << liveTot << endl;
	//const int mycount2=mycount;
	double* FAR2 = new double[mycount];
        double* eFAR2 = new double[mycount];
        double* N = new double[mycount];
        double* Ny = new double[mycount];
        double* eN = new double[mycount];
        for(int i =0 ;i<mycount;i++){N[i]=D[i];}
	Int_t  *myindex = new Int_t[mycount];
        TMath::Sort(mycount,Ny,myindex,false);
	for(int i =0 ;i<mycount;i++){FAR2[i]=FAR[myindex[i]];eFAR2[i]=eFAR[myindex[i]];N[i]=D[myindex[i]];Ny[i]=Dy[myindex[i]];eN[i]=0.0;}
	TGraphErrors* FarPlot= new TGraphErrors(mycount, N, FAR2, eN, eFAR2);
	FarPlot->SetLineColor(kBlue);
	FarPlot->SetMarkerColor(kBlue);
	//FarPlot->GetHistogram()->GetYaxis()->SetRangeUser(1.e-9,1.e-7);
	FarPlot->SetTitle("FAR vs slag distance");
	FarPlot->SetMarkerStyle(20);
  	FarPlot->SetMarkerSize(1);
    FarPlot->SetLineStyle(3);
  	//FarPlot->SetMinimum(0.5/liveMax);
  	FarPlot->GetHistogram()->SetXTitle("slag distance [day]");
  	FarPlot->GetHistogram()->SetYTitle("FAR, [1/year]");
	c1->Clear();
	c1->SetLogy(kFALSE);
    gStyle->SetOptStat();
    gStyle->SetOptFit();
    
	FarPlot->Draw("AP");
    FarPlot->Fit("pol1");
 	sprintf(fname,"%s/FAR_vs_slag.png",netdir);
        c1->Update(); c1->SaveAs(fname);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(kFALSE);

   // CWB::CBCTool cbcTool;
    // cbcTool.domyPoissonPlot(nIFO, Wlag, Tlag, Rlag, netdir);

    TCanvas *c2 = new TCanvas("Slag2", "Slag2",32,55,1450,502);
    //c1->Range(-0.75,-1.23433,6.75,8.17182);
    c2->SetBorderSize(1);
    c2->SetFillColor(0);
    c2->SetGridx();
    c2->SetGridy();
    c2->SetFrameFillColor(0);

  // if(write_ascii){
				  sprintf(fname,"%s/FAR_vs_slag2.txt",netdir);
				  FILE* frate = fopen(fname,"w");
				  fprintf(frate,"#FAR[years^-1] eFAR Tshift[days]\n");
				  for(int i=0; i<mycount;i++){fprintf(frate,"%.4f %.4f %.4f\n",FAR[i],eFAR[i],Dy[i]);}
				  fclose(frate);

  // }


    TGraphErrors* FarPlot2= new TGraphErrors(mycount, Dy, FAR, eN, eFAR);	
    FarPlot2->SetLineColor(kBlue);
    FarPlot2->SetMarkerColor(kBlue);
    //FarPlot2->SetTitle("False Alarm rate vs time shift");
    FarPlot2->SetTitle("");
    FarPlot2->SetMarkerStyle(20);
    FarPlot2->SetMarkerSize(1);
    FarPlot2->SetLineStyle(3);
    //FarPlot2->SetMinimum(0.5/liveMax);
    FarPlot2->GetHistogram()->GetXaxis()->SetRangeUser(TMath::FloorNint(SlagMin/86400.),TMath::CeilNint(SlagMax/86400.));
    FarPlot2->GetHistogram()->SetXTitle("time shift [day]");
    FarPlot2->GetHistogram()->SetYTitle("FAR [1/year]");
    FarPlot2->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
    FarPlot2->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    FarPlot2->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    FarPlot2->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    c2->Clear();
    c2->SetLogy(kFALSE);
    
    FarPlot2->Draw("AP");
    TF1 *g1    = new TF1("g1","pol1",0,SlagMax/86400.);
    g1->SetLineColor(2);
    g1->SetLineWidth(5);
    TF1 *g2    = new TF1("g2","pol1",-SlagMax/86400.,0);
    g2->SetLineColor(3);
    g2->SetLineWidth(5);
//    TF1 *total = new TF1("total","pol1(0)+pol1(2)",-SlagMax/86400.,SlagMax/86400.);
    FarPlot2->Fit(g1,"R");
    FarPlot2->Fit(g2,"R+","SAME");
//    Double_t par[4];
//   g1->GetParameters(&par[0]);
 //  g2->GetParameters(&par[2]);
   legr = new TLegend(0.6,0.8,0.885,0.98,"","brNDC");
   legl = new TLegend(0.15,0.8,0.4,0.98,"","brNDC");
   char lab[256];
   sprintf(lab,"#chi^{2}/ndf(right) : %3.3g / %d",g1->GetChisquare(),g1->GetNDF());
   legr->AddEntry(g1, lab, "l");
   sprintf(lab,"P0(right) : %.2f +/- %.2f",g1->GetParameters()[0],g1->GetParErrors()[0]);
   legr->AddEntry("", lab, "a");
   sprintf(lab,"P1(right) : %.2g +/- %.2g",g1->GetParameters()[1],g1->GetParErrors()[1]);
   legr->AddEntry("", lab, "a");
   legr->SetLineColor(1);

   sprintf(lab,"#chi^{2}/ndf(left) : %3.3g / %d",g2->GetChisquare(),g2->GetNDF());
   legl->AddEntry(g2, lab, "l");
   sprintf(lab,"P0(left) : %.2f +/- %.2f",g2->GetParameters()[0],g2->GetParErrors()[0]);
   legl->AddEntry("", lab, "a");
   sprintf(lab,"P1(left) : %.2g +/- %.2g",g2->GetParameters()[1],g2->GetParErrors()[1]);
   legl->AddEntry("", lab, "a");
   legl->SetLineColor(1);

   legr->Draw(); 
   legl->Draw(); 
//   total->SetParameters(par);
 //   FarPlot2->Fit(total,"R+");
    g1->Draw("same");
    g2->Draw("same");
    sprintf(fname,"%s/FAR_vs_slag2.png",netdir);
        c2->Update(); c2->SaveAs(fname);
//    gStyle->SetOptStat(kFALSE);
 //   gStyle->SetOptFit(kFALSE);

	exit(0);
	
}
