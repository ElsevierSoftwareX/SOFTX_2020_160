{
#define NBINS 1100
#define BINMIN 4.50
#define BINMAX 114.5
#define REBIN_MIN_EVTS_PER_BIN 30.0
#define NSLAG_ORDER 5
  //cout<<"slag starts..."<<endl;

  if(nIFO>3) {cout << "slag.C : Skipping plot production - wave tree should come from a 2- or 3-fold network : " << net_file_name << endl;exit(0);} 
  
  int NSLAG = pow(2,NSLAG_ORDER);
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

   int num = (int)wave.GetEntries(cut);
   wave.SetEstimate(num);
   TString sel("slag[1]:rho[1]");
   //TString sel("slag[1]:slag[2]");
  // sel+= TString::Itoa(nIFO-2,10);
  // sel+="]:slag[";
  // sel+=Itoa(nIFO-1,10);
  // sel+="]";
 // if(nIFO==2) {sel.ReplaceAll("1","0");sel.ReplaceAll("2","1");}
  // TString sel="slag[1]:slag[2]";
   wave.Draw(sel,cut,"goff");
  cout<<"Number of Entries: "<<num<<endl;
  double* slag1 = wave.GetV1(); 
  double* rho = wave.GetV2(); 


  int num2 = (int)live.GetEntries(cut2);
  live.SetEstimate(num2);
  //sel+=":live";
  live.Draw("slag[1]:live",cut2,"goff");
  double* lslag1 = live.GetV1();
  //double* lslag2 = live.GetV2();
  double* Live = live.GetV2();

 char mytitle[256];
 double SlagMax = wave.GetMaximum("slag")+segLen/2.;
 double SlagMin = wave.GetMinimum("slag")-segLen/2.;
 int NSlag = TMath::FloorNint((SlagMax-SlagMin)/segLen); 
 cout << "SLAG MAX : "<< wave.GetMaximum("slag")<< " s  SLAG MIN : "<< wave.GetMinimum("slag")<< " s  #SLAGS : "<<NSlag-1<<endl;
sprintf(mytitle,"FAR distribution over slags (post cat3 & rho>%f)",T_cut);
	TH2F* Slag = new TH2F("SLAG",mytitle,NBINS,BINMIN,BINMAX,NSLAG,SlagMin,SlagMax);
	Slag->GetYaxis()->SetTitle("slag[1] shift [s]");
	//Slag->GetXaxis()->SetNdivisions(10,kFALSE);
	Slag->GetXaxis()->SetTitle("rho");
	//Slag->GetYaxis()->SetNdivisions(10,kFALSE);
	Slag->SetStats(kFALSE);
	TH1F* lSlag = new TH1F("LSLAG","FAR distribution over slags",NSLAG,SlagMin,SlagMax);
	for(int i =0 ;i<num;i++){
		Slag->Fill(rho[i],slag1[i]);

	}
 	for(int i =0 ;i<num2;i++){
                lSlag->Fill(lslag1[i],Live[i]);

        }

TH1D* P0 = new TH1D("P0","P0",100,0.0,1.0);
//TH1D* P1 = new TH1D("P1","P1",NBINS,BINMIN,BINMAX);
double Pval[NSLAG-1];
int test = 0;

for(int j =0 ;j<NSLAG/2;j++){
  for(int i =1 ;i<NSLAG;i++){
	  TH1D* P1 = (TH1D*) Slag->ProjectionX("",i,i+j)->Clone();
	//P[i-1]= (TH1D*)Px0->Clone();
	//P1->Scale(1./lSlag->GetBinContent(i));	
	//TH1D* P2 = (TH1D*) Slag->ProjectionX("",i+1,i+1)->Clone();	
        //P[i]= (TH1D*)Px0->Clone();
	//P2->Scale(1./lSlag->GetBinContent(i+1));
	  P1->AndersonDarlingTest(Slag->ProjectionX("",i+j+1,i+2*j+1));	
//	Pval= Slag->ProjectionX("",i,i)->AndersonDarlingTest(Slag->ProjectionX("",i+1,i+1));
		cout <<j<<" "<<i<<" bins ["<<i<<","<<i+j<<"] bins ["<<i+j+1<<","<<i+2*j+1<<"]"<<endl;
      //  cout <<" Anderson-Darling test (no rescaling) "<<Pval<<endl;
	 // P0->Fill(Pval);
	  const Int_t nq = P1->GetEntries()/REBIN_MIN_EVTS_PER_BIN; cout<<"Rebinning to "<<nq<<" bins"<<endl;
	  Double_t* xq = new Double_t[nq];  // position where to compute the quantiles in [0,1]
	  Double_t* yq = new Double_t[nq];
	  Double_t* yq2 = new Double_t[nq+3];
	  for (Int_t k=0;k<nq;k++) xq[k] = Float_t(k+1)/nq;
	  P1->GetQuantiles(nq,yq,xq);

	  for(int l=1;l<nq+2;l++){yq2[l]=yq[l-2];}
	  yq2[0]= BINMIN;
	  yq2[1]= T_cut;
	  yq2[nq+2]= BINMAX;

	  TH1D *hnew2 = new TH1D("hnew2","rebinned",nq+2,yq2);
	  TH1D *hnew3 = new TH1D("hnew3","rebinned",nq+2,yq2);
      TH1D* Sh3 = new TH1D("Sh3","Significance vs rho bin",nq+2,yq2);
      TH1D* Sh4 = new TH1D("Sh4","Distribution of significance",200,-10.,10.);
      TH1D* Shn = new TH1D("Shn","Normal Distribution",200,-10.,10.);
      TRandom3 P;
      P.SetSeed(1234.);
      for(int m=0;m<num;m++){
				if((slag1[m]>=Slag->GetYaxis()->GetBinLowEdge(i))&&(slag1[m]<Slag->GetYaxis()->GetBinUpEdge(i+j))) {hnew2->Fill(rho[m]);}
				if((slag1[m]>=Slag->GetYaxis()->GetBinLowEdge(i+j+1))&&(slag1[m]<Slag->GetYaxis()->GetBinUpEdge(i+2*j+1))) {hnew3->Fill(rho[m]);}
 
	   }

      cout<<"After Rebinning"<<endl;
	 // cout <<"Anderson-Darling test (no rescaling)"<<endl;
	  Pval[test]= hnew2->AndersonDarlingTest(hnew3,"D");
	  double w_tmp;
	  //double w[nq];
	  int count =0;
	  for(int n=0;n<nq;n++){
		  if((hnew2->GetBinContent(n+1)>0.0) ||(hnew3->GetBinContent(n+1)>0.0)){
				w_tmp = (hnew2->GetBinContent(n+1)-hnew3->GetBinContent(n+1))/TMath::Sqrt(pow(hnew2->GetBinError(n+1),2)+pow(hnew3->GetBinError(n+1),2));
				Sh3->SetBinContent(n+1,w_tmp);
				Sh4->Fill(Sh3->GetBinContent(n+1));
			//	w[n]=w_tmp;
				count++;
		  } 
		//  else{Sh3->SetBinContent(n+1,0.0);w[n]=0.0;}
	  }
	  for(int n=0;n<count;n++){Shn->Fill(P.Gaus(0,1));}
	  cout <<"Anderson-Darling 2-sample test significance distribution vs Normal (no rescaling)"<<endl;
	//  Pval[test]=Sh4->AndersonDarlingTest(Shn,"D");
      Sh4->AndersonDarlingTest(Shn,"D");

	  //ROOT::Math::GoFTest* goftest_1 = new ROOT::Math::GoFTest(count, w, ROOT::Math::GoFTest::kGaussian);
	  //Double_t pvalueAD_1 = goftest_1-> AndersonDarlingTest();
	  //cout <<"Anderson-Darling 1-sample test significance distribution vs Normal (no rescaling) : "<<pvalueAD_1<<endl;
      /*hnew2->Delete();
      hnew3->Delete();
      Sh3->Delete();
      Sh4->Delete();
      Shn->Delete();*/
      //delete w;
         
	   P0->Fill(Pval[test]);
	   test++;
      i+=2*j+1;	
	//Slag->ProjectionX("",i,i)->AndersonDarlingTest(Slag->ProjectionX("",i+1,i+1),"D");
	//P0->AndersonDarlingTest(P1,"D");
  }
  j+=j;
}

Int_t *index = new Int_t[test];
TMath::Sort(test,Pval,index,false);

float *pvals = new float[test];
float *rank = new float[test];
for(int i=0;i<test;i++){pvals[i]=Pval[index[i]];rank[i]=(i+1.)/test;cout<<pvals[i]<<" "<<rank[i]<<endl;}
TGraph Gr(test,rank,pvals);
Gr.Draw("alp");
//P0->Draw();

	//Slag->Divide(lSlag);
	//c1->Clear();
	//Slag->Draw("colz");
 	//sprintf(fname,"%s/slag.png",netdir);
  	//c1->Update(); c1->SaveAs(fname);
/*
	TAxis *xaxis = Slag->GetXaxis();
	TAxis *yaxis = Slag->GetYaxis();
	int mycount = 0;
        double FAR[1000];
        double eFAR[1000];
        double D[1000];
        for(int i =1 ;i<NSlag+1;i++){
                for(int j =1 ;j<NSlag+1;j++){
                        if(Slag->GetBinContent(i,j)>=1.E-15){
				FAR[mycount] = Slag->GetBinContent(i,j);
				D[mycount] = TMath::Sqrt(pow(xaxis->GetBinCenter(i),2)+pow(yaxis->GetBinCenter(j),2));
				//D[mycount] = yaxis->GetBinCenter(j);
				eFAR[mycount] = Slag->GetBinError(i,j)/lSlag->GetBinError(i,j);
				mycount++;
				cout << "slag[1] :"<<xaxis->GetBinCenter(i)<<" slag[2] :"<<yaxis->GetBinCenter(j) <<" FAR :"<<FAR[mycount-1] << "+/-" <<eFAR[mycount-1] << " " <<endl;
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
        double* eN = new double[mycount];
        for(int i =0 ;i<mycount;i++){N[i]=D[i];}
	Int_t  *myindex = new Int_t[mycount];
        TMath::Sort(mycount,N,myindex,false);
	for(int i =0 ;i<mycount;i++){FAR2[i]=FAR[myindex[i]];eFAR2[i]=eFAR[myindex[i]];N[i]=D[myindex[i]];eN[i]=0.0;}
	TGraphErrors* FarPlot= new TGraphErrors(mycount, N, FAR2, eN, eFAR2);
	FarPlot->SetLineColor(kBlue);
	FarPlot->SetMarkerColor(kBlue);
	//FarPlot->GetHistogram()->GetYaxis()->SetRangeUser(1.e-9,1.e-7);
	FarPlot->SetTitle("FAR vs slag distance");
	FarPlot->SetMarkerStyle(20);
  	FarPlot->SetMarkerSize(1);
  	//FarPlot->SetMinimum(0.5/liveMax);
  	FarPlot->GetHistogram()->SetXTitle("#slag distance [s]");
  	FarPlot->GetHistogram()->SetYTitle("rate, Hz");
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
	exit(0);
*/	
}
