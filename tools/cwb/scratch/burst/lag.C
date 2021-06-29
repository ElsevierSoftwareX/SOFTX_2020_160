{

  cout<<"Lag starts..."<<endl;

  if(nIFO!=3) {cout << "lag.C : Skipping plot production - wave tree should come from a 3-fold network : " << net_file_name << endl;exit(0);}
 
  CWB::Toolbox TB;

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
  TCanvas *c1 = new TCanvas("lag", "lag",32,55,750,502);
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
   wave.Draw("lag[0]-lag[1]:lag[1]-lag[2]",cut,"goff");
  cout<<"Number of Entries: "<<num<<endl;
  double* lag01 = wave.GetV1(); 
  double* lag12 = wave.GetV2(); 


  int num2 = (int)live.GetEntries(cut2);
  live.SetEstimate(num2);
  live.Draw("lag[0]-lag[1]:lag[1]-lag[2]:live",cut2,"goff");
  double* llag01 = live.GetV1();
  double* llag12 = live.GetV2();
  double* Live = live.GetV3();

 char title[256];
 double xplotlagMax = TMath::MaxElement(num,lag01);
 double xplotlagMin = TMath::MinElement(num,lag01);
double yplotlagMax = TMath::MaxElement(num,lag12);
 double yplotlagMin = TMath::MinElement(num,lag12);

 cout<<xplotlagMin<<" < lag[o]-lag[1] < "<<xplotlagMax<<endl;
 cout<<yplotlagMin<<" < lag[1]-lag[2] < "<<yplotlagMax<<endl;
// int Nlag = TMath::FloorNint((plotlagMax-plotlagMin)/segLen); 
 
//double plotlagMax = 300;
//double plotlagMin = -300;
int Nlag = 100;

sprintf(title,"FAR distribution over lags (post cat3 & rho>%f)",T_cut);
	TH2F* lag = new TH2F("SLAG",title,Nlag,xplotlagMin,xplotlagMax,Nlag,yplotlagMin,yplotlagMax);
	lag->GetXaxis()->SetTitle("lag[0]-lag[1] shift [s]");
	lag->GetXaxis()->SetNdivisions(10,kFALSE);
	lag->GetYaxis()->SetTitle("lag[1]-lag[2] shift [s]");
	lag->GetYaxis()->SetNdivisions(10,kFALSE);
	lag->SetStats(kFALSE);
	TH2F* llag = new TH2F("LSLAG","FAR distribution over lags",Nlag,xplotlagMin,xplotlagMax,Nlag,yplotlagMin,yplotlagMax);
	for(int i =0 ;i<num;i++){
		lag->Fill(lag01[i],lag12[i]);

	}
 	for(int i =0 ;i<num2;i++){
                llag->Fill(llag01[i],llag12[i],Live[i]);

        }

	lag->Divide(llag);
	c1->Clear();
	lag->Draw("colz");
 	sprintf(fname,"%s/lag.png",netdir);
  	c1->Update(); c1->SaveAs(fname);

	TAxis *xaxis = lag->GetXaxis();
	TAxis *yaxis = lag->GetYaxis();
	int count = 0;
        double FAR[100000];
        double eFAR[100000];
        double D[100000];
        for(int i =1 ;i<Nlag+1;i++){
                for(int j =1 ;j<Nlag+1;j++){
                        if(lag->GetBinContent(i,j)>=1.E-15){
				FAR[count] = lag->GetBinContent(i,j);
				D[count] = TMath::Sqrt(pow(xaxis->GetBinCenter(i),2)+pow(yaxis->GetBinCenter(j),2));
				eFAR[count] = lag->GetBinError(i,j)/llag->GetBinError(i,j);
				count++;
				cout << "i :"<<i<<" j :"<<j<<" FAR :"<<FAR[count-1] << "+/-" <<eFAR[count-1] << " " <<endl;
			}

		}		
        }

	
 	cout << "Number of lags : " << count << endl;
	double liveTot = llag->GetSumOfWeights();
	double liveMax = llag->GetMaximum();
	cout << "Total BKG live time : " << liveTot << endl;
	double FAR2[count];
        double eFAR2[count];
        double N[count];
        double eN[count];
        for(int i =0 ;i<count;i++){N[i]=D[i];}
	Int_t  *index = new Int_t[count];
        TMath::Sort(count,N,index,false);
	for(int i =0 ;i<count;i++){FAR2[i]=FAR[index[i]];eFAR2[i]=eFAR[index[i]];N[i]=D[index[i]];eN[i]=0.0;}
	TGraphErrors* FarPlot= new TGraphErrors(count, N, FAR2, eN, eFAR2);
	FarPlot->SetLineColor(kBlue);
	FarPlot->SetMarkerColor(kBlue);
	//FarPlot->GetHistogram()->GetYaxis()->SetRangeUser(1.e-9,1.e-7);
	FarPlot->SetTitle("FAR vs lag distance");
	FarPlot->SetMarkerStyle(20);
  	FarPlot->SetMarkerSize(1);
  	//FarPlot->SetMinimum(0.5/liveMax);
  	FarPlot->GetHistogram()->SetXTitle("#lag distance [s]");
  	FarPlot->GetHistogram()->SetYTitle("rate, Hz");
	c1->Clear();
	//c1->SetLogy(kTRUE);
	FarPlot->Draw("AP");
 	sprintf(fname,"%s/FAR_vs_lag.png",netdir);
        c1->Update(); c1->SaveAs(fname);	

	exit(0);

}
