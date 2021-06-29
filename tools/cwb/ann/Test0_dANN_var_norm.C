/*
# Copyright (C) 2019 Serena Vinciguerra
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


#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPaletteAxis.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TMath.h"
#include "mdc.hh"
#include "watplot.hh"
#include <vector>


//#define nINP 64
#define nIFO 3
#define nRHO 3
#define nANN 10
#define deltaANN 0.1
#define nCC 4
#define NPIX 20
#define RHOth 5.5
#define nth 10
using namespace std;
void graph(TString ifile);
void Annth(TString ifile);
void Plots(TString ifile);
//void Test0(TString NN_FILE,TString TEST_FILE,TString ofile, int TS=0, int TB=0, int s=0, int b=0, int uf=0){
void Test0_dANN_var_norm(TString NN_FILE,TString TEST_FILE, int TS=0, int TB=0, int s=0, int b=0, int uf=0){

   
   int p=0;
   char NNi[1024];
   char NNi2[1024];
   sprintf(NNi2,"%s",NN_FILE.Data());
   while (NNi2[p]){
	//cout<<NNi2[p]<<endl;
	if (NNi2[p]=='N'){
//	cout<<p<<endl;
                        if((NNi2[p+1]=='N')&&(NNi2[p+2]=='\/')) {
//					cout<<" dentro if"<<endl;
					int hh=p+3;
					while (NNi2[hh]!='\0'){NNi[hh-p-3]=NNi2[hh];hh=hh+1;}
					break;
					}
	}
	p=p+1;
   }

   int q=0;
   char Filei[1024];
//   for (int i=0;i<1024;i++)Filei[i]='';
   char Filei2[1024];
   sprintf(Filei2,"%s",TEST_FILE.Data());
 while (Filei2[q]){
// cout<<Filei2[q]<<endl;
 if(Filei2[q]=='n'&&Filei2[q+1]=='n'&&Filei2[q+2]=='T'&&Filei2[q+3]=='R'&&Filei2[q+4]=='E'&&(Filei2[q+5]=='E')&&(Filei2[q+6]=='\/')) {
               int hh=q+7;
               //while (Filei2[hh]!='\0') {Filei[hh-p-8]=Filei2[hh];hh=hh+1;cout<<Filei2[hh]<<endl;}
               while (Filei2[hh]!='\0') {Filei[hh-p-7]=Filei2[hh];hh=hh+1;
					//cout<<Filei2[hh]<<endl;
					}
	       for (int h0=hh-p-7;h0<1024;h0++) Filei[h0]='\0';
               break;
               }
        q=q+1;
   }

   cout<<NNi<<" original String: "<<NNi2<<endl;
   cout<<Filei<<" original String: "<<Filei2<<endl;
   
   TString NNi0(NNi);
   TString Filei0(Filei);

   NNi0.ReplaceAll(".root","");
   Filei0.ReplaceAll(".root","");
   

   TFile* fnet =TFile::Open(NN_FILE.Data());
   TMultiLayerPerceptron *mlp = (TMultiLayerPerceptron*)fnet->Get("TMultiLayerPerceptron");
   if(mlp==NULL) {cout << "Error getting mlp" << endl;exit(1);}
   TTree* infot=(TTree*)fnet->Get("info");
   int NNs;
   int NNb;
   int NNnTS;
   int NNnTB;
   infot->SetBranchAddress("Rand_start_Sig",&NNs);
   infot->SetBranchAddress("Rand_start_Bg",&NNb);
   infot->SetBranchAddress("#trainSig",&NNnTS);
   infot->SetBranchAddress("#trainBg",&NNnTB);
   infot->GetEntry(0);
   
   TFile* fTEST =TFile::Open(TEST_FILE.Data());
   TTree* NNTree=(TTree*)fTEST->Get("nnTree");
   int entries=NNTree->GetEntries();
   cout<<"entries: "<<entries<<endl;

//estraggo le info che servono anche al tree di def dell'mlp
   int ndim;
   int ninp;
   int y;
   int entriesTot;
   int bg_entries;
   int sig_entries;
   
   NNTree->SetBranchAddress("#Entries_type",&entriesTot);
   NNTree->SetBranchAddress("Matrix_dim",&ndim);
   NNTree->SetBranchAddress("#inputs",&ninp);
   NNTree->SetBranchAddress("amplitude_mode",&y);

   NNTree->GetEntry(0);
   sig_entries=entriesTot;
   NNTree->GetEntry(entries-1);
   bg_entries=entriesTot;

   int const NDIM=ndim;
   int const nINP=ninp;
   
   cout<<"NDIM "<<NDIM<<endl;
   cout<<"nINP"<<nINP<<endl;
   cout<<"sig e "<<sig_entries<<endl;
   cout<<"bg e "<<bg_entries<<endl;
   int minevents=0;
   if (sig_entries>bg_entries) minevents=bg_entries;
   else minevents=sig_entries;

   if(b==0) b=sig_entries;
   if(TB==0 && TS==0){
     TS=minevents;
     TB=minevents;
   }
   if (b<sig_entries){
	cout<<"Error: Bg index<sig_entries"<<endl;
	exit(0);
	} 
   if((TS>sig_entries||TB>bg_entries)&&(TS==TB)) {TS=minevents-s;TB=minevents-b+sig_entries;}
   if((TS>sig_entries||TB>bg_entries)&&(TS!=TB)) {TS=sig_entries-s;TB=bg_entries-b+sig_entries;}
   

   char NOMEtot[1024];
   sprintf(NOMEtot,"NN_norm_%s_FILE_%s_TS%i_TB%i_st%i_bt%i_dANN%1.2f_uf%i",NNi0.Data(),Filei0.Data(),TS,TB,s,b,deltaANN,uf);
   cout<<"nome: "<<NOMEtot<<endl;

   TString SIG_FILE;
   TString BG_FILE;
   char FILE_NAME[516];
   NNTree->SetBranchAddress("Files_name",&FILE_NAME);
   NNTree->GetEntry(0);
   SIG_FILE=FILE_NAME;
   NNTree->GetEntry(entries-1);
   BG_FILE=FILE_NAME;
   cout<<"fine ifdef RHO_CC"<<endl;
  
   TChain sigTree("waveburst");//cerca il Tree "waveburst nei file
   sigTree.Add(SIG_FILE.Data());//determina i file
   netevent signal(&sigTree,nIFO);
   int sig_entries2 = signal.GetEntries();
   cout << "sig entries2 : " << sig_entries2 << endl;

   TChain bgTree("waveburst");
   bgTree.Add(BG_FILE.Data());
   netevent background(&bgTree,nIFO);
   int bg_entries2 = background.GetEntries();
   cout << "bg entries2 : " << bg_entries2 << endl;
  
   cout<<"b: "<<b<<endl;
   cout<<"s: "<<s<<endl;

   // add leaf
   Float_t x[nINP];
   Float_t xo[nINP];
   for (int jj=0; jj<nINP;jj++) x[jj]=0.;
   char ilabel[nINP][16];
  
      //costruzione una foglia  per ogni input
   for(int i=0;i<nINP;i++) {
       sprintf(ilabel[i],"x%i",i+1);
       NNTree->SetBranchAddress(ilabel[i], &xo[i]);
   }
  
char ofile[1024];
sprintf(ofile,"outfile/%s.root",NOMEtot);
//TFile*f=new TFile(ofile.Data(),"RECREATE");
TFile*f=new TFile(ofile,"RECREATE");
   TTree* NNTree2=new TTree("Parameters","Parameters");
NNTree2->SetDirectory(f);
   double out=0.;
   double NNcc=0.;
   double NNrho=0.;
  NNTree2->Branch("ANNout",&out,"ANNout/D");
  NNTree2->Branch("cc",&NNcc,"cc/D");
  NNTree2->Branch("rho",&NNrho,"rho/D");
  char NNfilename[1024];
  NNTree2->Branch("NNname",&NNfilename,"NNname/C");
  sprintf(NNfilename,"%s",NN_FILE.Data());
  char Testf[1024];
  NNTree2->Branch("TestFile",&Testf,"TestFile/C");
  sprintf(Testf,"%s",TEST_FILE.Data());
  int nTestS;
  NNTree2->Branch("#TestSig",&nTestS,"#TestSig/I");
//  nTestS=TS;
  cout<<"nTestS: "<<nTestS<<" TS: "<<TS<<endl; 
  cout<<"dopo def tree"<<endl; 
/*   TChain TreeS("waveburst");//cerca il Tree "waveburst nei file
   TreeS.Add(TEST_FILE.Data());//determina i file
   netevent entryS(&TreeS,nIFO);
   int entriesTot=0;
   entriesTot = entryS.GetEntries();
   cout << "Sig entries : " << entriesTot << endl;
   std::vector<double>* frameS = new vector<double>;
   entryS.fChain->SetBranchAddress("nnframe", &frameS);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc

   //cout<<"cc"<<(double)entryS.netcc[1]<<endl;
*/   //cout<<"rho"<<(double)entryS.rho[0]<<endl;

int scount=0;
  if(uf==0) nTestS=TS;
  else {
for(int n=s;n<s+TS;n++) {
     if (uf!=0&&n>=NNs&&n<=(NNs+NNnTS)) continue;
        scount=scount+1;
        }
nTestS=scount;
}
   double params[nINP];
   int sig_05=0;
   for (int i=0; i<nINP; i++) params[i]=0.;
     //for(int n=0;n<entryS.GetEntries();n++) {
   for(int n=s;n<s+TS;n++) {
     if (uf!=0&&n>=NNs&&n<=(NNs+NNnTS)) continue;	
     NNTree->GetEntry(n);
     signal.GetEntry(n);
     NNcc=(double)signal.netcc[1];
     NNrho=(double)signal.rho[0];
     float xo_max=0.;
     xo_max=xo[0];
    for (int i=0; i<nINP; i++) {if(xo[i]>xo_max) xo_max=xo[i];}
     for (int i=0; i<nINP; i++){
	//x[i]=(*frameS)[i];
        x[i]=xo[i]/xo_max;
	params[i]=x[i];	
	}
     double output=mlp->Evaluate(0,params);
     out=output;
    if (out<0.6) sig_05=sig_05+1;
    cout<<"rho: "<<signal.rho[0]<<" cc "<<signal.netcc[1]<<" out: "<<out<<endl; 
     NNTree2->Fill();
    }
       
cout<<"riempito sig"<<endl;
int bg_05=0;
   //for(int n=sig_entries+b;n<sig_entries+b+TB;n++) {
   for(int n=b;n<b+TB;n++) {
     if (uf!=0&&n>=NNb&&n<=(NNb+NNnTB)) continue;	
     NNTree->GetEntry(n);
     cout<<"n: "<<n<<"Bg index"<<(n-sig_entries)<<endl;
     background.GetEntry(n-sig_entries);
     NNcc=(double)background.netcc[1];
     NNrho=(double)background.rho[0];
     float xo_max=0.;
     xo_max=xo[0];
     for (int i=0; i<nINP; i++){if(xo[i]>xo_max) xo_max=xo[i];}
     for (int i=0; i<nINP; i++){
        //x[i]=(*frameS)[i];
        x[i]=xo[i]/xo_max;
        params[i]=x[i];
        }

     double output=mlp->Evaluate(0,params);
     out=output;
	if(out>0.6) bg_05=bg_05+1;
    cout<<"rho: "<<background.rho[0]<<" cc "<<background.netcc[1]<<" out: "<<out<<endl; 
     NNTree2->Fill();
     }
cout<<"s "<<s<<" TS "<<TS<<" b "<<b<<" TB "<< TB<<endl;
cout<<"riempito bg"<<endl;
int Realentries=0;
Realentries=NNTree2->GetEntries();
cout<<" Realentries "<<Realentries<<endl;
NNTree2->Write();
f->Close();
cout<<" Realentries "<<Realentries<<endl;
cout<<"chiuso file"<<endl;
//graph(ofile.Data());
graph(ofile);
cout<<"dopo richiamo funzione"<<endl;
//Annth(ofile);
Plots(ofile);
cout<<ofile<<endl;
cout<<" Sig with out<06: "<<sig_05<<endl;
cout<<" Bg with out>06: "<<bg_05<<endl;
cout<<" Realentries "<<Realentries<<endl;

}

void graph(TString ifile){
	TString name(ifile);
	name.ReplaceAll("outfile/","");
	TFile* fTEST =TFile::Open(ifile.Data());
	TTree* NNTree2=(TTree*)fTEST->Get("Parameters");	           
	double out;
        double cc;	
	double rho;
	int nSi;
	NNTree2->SetBranchAddress("ANNout",&out);
	NNTree2->SetBranchAddress("cc",&cc);
	NNTree2->SetBranchAddress("rho",&rho);
	NNTree2->SetBranchAddress("#TestSig",&nSi);
	NNTree2->GetEntry(0);
	int const nSig=nSi;
	cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
        int const ncurve=nANN*nCC;
	cout<<"dentro funzione dopodef"<<endl;
//LOG(NBg)vs RHO------------------------------------------
        double* rhoSig[ncurve];
        for (int i=0;i<ncurve;i++) rhoSig[i]=new double[nSig];
	cout<<"dopo def rhoSig"<<endl;
        //double rhoSig[20][10];
        int NSig[ncurve];
//	int nBg=0;
	int const nBg=NNTree2->GetEntries()-nSig;
	for (int i=0;i<ncurve;i++) {
		NSig[i]=0;
		for (int j=0;j<nSig;j++) rhoSig[i][j]=0.;	
	}
	cout<<"dopo def rhoSig"<<endl;
        double* rhoBg[ncurve];
        for (int i=0;i<ncurve;i++) rhoBg[i]=new double[nBg];
        //double rhoBg[20][10];
        int NBg[ncurve];
	for (int i=0;i<ncurve;i++) {
		NBg[i]=0;
		for (int j=0;j<nBg;j++) rhoBg[i][j]=0.;	
	}
	cout<<"dopo def rhoBg"<<endl;
	double ccTh[nCC];
        for (int i=0;i<nCC;i++) ccTh[i]=0.;
        double NNTh[nANN];
        for (int i=0;i<nANN;i++) NNTh[i]=0.;
	double deltacc=0.;	
//	double deltaANN=0.;
	deltacc=0.2/nCC;
	//deltaANN=0.6/(nANN-1);

	cout<<NNTree2->GetEntries()<<endl;	
	for(int n=0;n<NNTree2->GetEntries();n++){
		cout<<n<<endl;
		NNTree2->GetEntry(n);
		cout<<"rho "<<rho<<" cc "<<cc<<" out "<<out<<endl;
        	for(int i=0; i<nCC;i++){
       			ccTh[i]=0.5+i*deltacc;
			if(cc<ccTh[i]) continue;
			for(int j=0; j<nANN;j++){
				if(j==0) NNTh[j]=-1000.;
				//else NNTh[j]=0.5+(j-1)*deltaANN;
				//else NNTh[j]=0.1+(j-1)*deltaANN;
				else NNTh[j]=0.+(j)*deltaANN;
				//else NNTh[j]=0.+(j-1)*deltaANN/1000.;
				//else NNTh[j]=1.;
				if(out<NNTh[j]) continue;
					int ni=0;
					if(n>nSig) {
						NBg[i*nANN+j]= NBg[i*nANN+j]+1;
						while(rhoBg[i*nANN+j][ni]!=0)ni=ni+1;
						rhoBg[i*nANN+j][ni]=rho;
					//	cout<<"rho: "<<rho<<" colonna "<<i*nANN+j<<" riga "<<ni<<endl;
//						cout<<" soglia_cc "<<ccTh[i]<<" soglia_ANN "<<NNTh[j]<<endl;
						}
					else {
						NSig[i*nANN+j]= NSig[i*nANN+j]+1;
						while(rhoSig[i*nANN+j][ni]!=0)ni=ni+1;
	                                        rhoSig[i*nANN+j][ni]=rho;
					//	cout<<"rho: "<<rho<<" colonna "<<i*nANN+j<<" riga "<<ni<<endl;
//						cout<<" soglia_cc "<<ccTh[i]<<" soglia_ANN "<<NNTh[j]<<endl;
					}
			//	}
			}
		}
	}		
	cout<<"dopo riempimento variabili"<<endl;
	int* indexS[ncurve];
        for (int i=0;i<ncurve;i++) indexS[i]=new int[nSig];

        TGraph * gS[ncurve];
	for (int y=0;y<ncurve;y++) {
		int igS=0;
		int igS_p=0;
	/*	int indexS[nSig];
		double rhoS[nSig];
		for(int i=0;i<nSig;i++) {
			rhoS[i]=0.;
			indexS[i]=0;
		}*/
	//	ig=1;
	//	for(int i=0;i<nSig;i++) rhoS[i]=rhoSig[y][i];
		gS[y]=new TGraph();
//		gS[y]->SetMarkerStyle(7);
		TMath::Sort(nSig,rhoSig[y],indexS[y],false);
//		cout<<"dopo Sort "<<y<<endl;
		for (int k=0;k<nSig;k++) {
			int ii=indexS[y][k];
			int yy=0;
			if (k>0){
//				cout<<"k "<<k<<endl;
				int ij=indexS[y][k-1];
				//if(rhoSig[y][ii]!=0&&rhoSig[y][ii]!=rhoSig[y][ij]) {
				if(rhoSig[y][ii]!=0) {
							yy=NSig[y]-igS;
							//gS[y]->SetPoint(igS,rhoSig[y][ii],yy);
							if(rhoSig[y][ii]!=rhoSig[y][ij]) gS[y]->SetPoint(igS_p++,rhoSig[y][ii],yy);
							cout<<"igS"<<igS<<" x "<<rhoSig[y][ii]<<" y: "<<yy<<endl;
							igS=igS+1;
							}
			}
			else {
				if(rhoSig[y][ii]!=0){
					yy=NSig[y]-igS;
					gS[y]->SetPoint(0,rhoSig[y][ii],yy);
					igS=igS+1;
//					cout<<" x "<<rhoSig[y][ii]<<" y: "<<yy<<endl;
				}
	
			}
		}
	}
 //       cout<<"dopo inserimento puntiiS"<<endl;
//	cout<<ncurve<<endl;
        int* indexB[ncurve];
        for (int i=0;i<ncurve;i++) indexB[i]=new int[nBg];

	 TGraph * gB[ncurve];
        for (int y=0;y<ncurve;y++) {
		int igB=0;
		int igB_p=0;
                gB[y]=new TGraph();
                TMath::Sort(nBg,rhoBg[y],indexB[y],false);
//		cout<<"dopo Sort "<<y<<endl;
                for (int k=0;k<nBg;k++) {
                        int ii=indexB[y][k];
                        int yy=0;
			if (k>0){
				int ij=indexB[y][k-1];
                                //if(rhoBg[y][ii]!=0&&rhoBg[y][ii]!=rhoBg[y][ij]) {
                                if(rhoBg[y][ii]!=0) {
						yy=NBg[y]-igB;
						//gB[y]->SetPoint(igB,rhoBg[y][ii],yy);
						if(rhoBg[y][ii]!=rhoBg[y][ij]) gB[y]->SetPoint(igB_p++,rhoBg[y][ii],yy);
						igB=igB+1;
//						cout<<"igB"<<igB<<" x "<<rhoBg[y][ii]<<" y: "<<yy<<endl;
						}
			}
                      else {
				if(rhoBg[y][ii]!=0)	{
					yy=NBg[y]-igB;
					gB[y]->SetPoint(0,rhoBg[y][ii],yy);
					igB=igB+1;
//					cout<<"igB"<<igB<<" x "<<rhoBg[y][ii]<<" y: "<<yy<<endl;
				}
               		}
		 }
        }
        cout<<"dopo inserimento puntiB"<<endl;
	//gB[1]->SetMarkerStyle(7);
	//gB[1]->SetPoint(1,1,2);
	TCanvas* cS=new TCanvas("Efficiency_vs_rho","Efficiency_vs_rho",0,0,1200,700);
	cS->Divide(2,2);
	cS->cd(1)->SetLogy();
	TMultiGraph* mg1=new TMultiGraph();
//	gS[0]->SetMarkerStyle(7);
	gS[0]->SetMarkerColor(2);
	gS[0]->SetLineColor(2);
	mg1->SetTitle("cc=0.5;rho;#Events");
	if(gS[0]->GetN()!=0) mg1->Add(gS[0]);
//	gS[0]->Draw("apl");
	for (int h=1;h<nANN;h++){
		gS[h]->SetMarkerColor(3);	
		gS[h]->SetLineColor(3);
//		gS[h]->SetMarkerStyle(h+1);	
		if(gS[h]->GetN()!=0) mg1->Add(gS[h]);
	//	gS[h]->Draw("apl,same");	
		}
	mg1->Draw("apl");
	cS->cd(2)->SetLogy();
	TMultiGraph* mg2=new TMultiGraph();
//	gS[nANN]->SetMarkerStyle(7);
	gS[nANN]->SetMarkerColor(2);
	gS[nANN]->SetLineColor(2);
        mg2->SetTitle("cc=0.55;rho;#Events");
	if(gS[nANN]->GetN()!=0) mg2->Add(gS[nANN]);
	for (int h=1;h<nANN;h++){
		gS[nANN+h]->SetMarkerColor(3);	
		gS[nANN+h]->SetLineColor(3);
//		gS[nANN+h]->SetMarkerStyle(h+1);	
		if(gS[nANN+h]->GetN()!=0) mg2->Add(gS[nANN+h]);
          //     gS[nANN+h]->Draw("apl,same");
                }
	mg2->Draw("apl");
	cS->cd(3)->SetLogy();
	TMultiGraph* mg3=new TMultiGraph();
//	gS[nANN*2]->SetMarkerStyle(7);
	gS[nANN*2]->SetMarkerColor(2);
	gS[nANN*2]->SetLineColor(2);
	mg3->SetTitle("cc=0.6;rho;#Events");
	if(gS[nANN*2]->GetN()!=0) mg3->Add(gS[nANN*2]);
        for (int h=1;h<nANN;h++){
		gS[2*nANN+h]->SetMarkerColor(3);	
		gS[2*nANN+h]->SetLineColor(3);
//		gS[2*nANN+h]->SetMarkerStyle(h+1);	
		if(gS[2*nANN+h]->GetN()!=0) mg3->Add(gS[2*nANN+h]);
             //   gS[2*nANN+h]->Draw("apl,same");
                }
	mg3->Draw("apl");
	cS->cd(4)->SetLogy();
	TMultiGraph* mg4=new TMultiGraph();
//	gS[nANN*3]->SetMarkerStyle(7);
	gS[nANN*3]->SetMarkerColor(2);
	gS[nANN*3]->SetLineColor(2);
	mg4->SetTitle("cc=0.65;rho;#Events");
	if(gS[nANN*3]->GetN()!=0) mg4->Add(gS[nANN*3]);
   //     gS[nANN*3]->Draw("apl");
        for (int h=1;h<nANN;h++){
		gS[3*nANN+h]->SetMarkerColor(3);	
		gS[3*nANN+h]->SetLineColor(3);
//		gS[3*nANN+h]->SetMarkerStyle(h+1);	
		if(gS[3*nANN+h]->GetN()!=0) mg4->Add(gS[3*nANN+h]);
              //  gS[3*nANN+h]->Draw("apl,same");
                }
	mg4->Draw("apl");
	cout<<"nuovo canv"<<endl;
	TCanvas* cB=new TCanvas("Number_vs_rho","Number_vs_rho",0,0,1200,700);
	cB->Divide(2,2);
	cB->cd(1)->SetLogy();
	TMultiGraph* mg1B=new TMultiGraph();
//	gB[0]->SetMarkerStyle(7);
	gB[0]->SetMarkerColor(2);
	gB[0]->SetLineColor(2);
	mg1B->SetTitle("cc=0.5;rho;#Events");
	if(gB[0]->GetN()!=0) mg1B->Add(gB[0]);
	for (int h=1;h<nANN;h++){
		gB[h]->SetMarkerColor(3);	
		gB[h]->SetLineColor(3);
	//	gB[h]>SetMarkerStyle(h+1);	
		if(gB[h]->GetN()!=0) mg1B->Add(gB[h]);
		//gB[h]->Draw("apl,same");	
		}
	mg1B->Draw("apl");	
	cB->cd(2)->SetLogy();
	TMultiGraph* mg2B=new TMultiGraph();
	//gB[nANN]->SetMarkerStyle(7);
	gB[nANN]->SetMarkerColor(2);
	gB[nANN]->SetLineColor(2);
	mg2B->SetTitle("cc=0.55;rho;#Events");
	if(gB[nANN]->GetN()!=0) mg2B->Add(gB[nANN]);
	for (int h=1;h<nANN;h++){
		gB[nANN+h]->SetMarkerColor(3);	
		gB[nANN+h]->SetLineColor(3);
		//gB[nANN+h]>SetMarkerStyle(h+1);	
		if(gB[nANN+h]->GetN()!=0) mg2B->Add(gB[nANN+h]);
		//gB[h]->Draw("apl,same");	
		}
	mg2B->Draw("apl");	
	cB->cd(3)->SetLogy();
	TMultiGraph* mg3B=new TMultiGraph();
//	gB[2*nANN]->SetMarkerStyle(7);
	gB[2*nANN]->SetMarkerColor(2);
	gB[2*nANN]->SetLineColor(2);
	mg3B->SetTitle("cc=0.6;rho;#Events");
	if(gB[2*nANN]->GetN()!=0) mg3B->Add(gB[2*nANN]);
	for (int h=1;h<nANN;h++){
		gB[2*nANN+h]->SetMarkerColor(3);	
		gB[2*nANN+h]->SetLineColor(3);
//		gB[2*nANN+h]>SetMarkerStyle(h+1);	
		if(gB[2*nANN+h]->GetN()!=0) mg3B->Add(gB[2*nANN+h]);
		//gB[h]->Draw("apl,same");	
		}
	mg3B->Draw("apl");	
	cB->cd(4)->SetLogy();
	TMultiGraph* mg4B=new TMultiGraph();
//	gB[3*nANN]->SetMarkerStyle(7);
	gB[3*nANN]->SetMarkerColor(2);
	gB[3*nANN]->SetLineColor(2);
	mg4B->SetTitle("cc=0.65;rho;#Events");
	if(gB[3*nANN]->GetN()!=0) mg4B->Add(gB[3*nANN]);
	for (int h=1;h<nANN;h++){
		gB[3*nANN+h]->SetMarkerColor(3);	
		gB[3*nANN+h]->SetLineColor(3);
//		gB[3*nANN+h]>SetMarkerStyle(h+1);	
		if(gB[3*nANN+h]->GetN()!=0) mg4B->Add(gB[3*nANN+h]);
		//gB[h]->Draw("apl,same");	
		}
	mg4B->Draw("apl");	

//	cout<<"dopo Draw()"<<endl;
	TString CnameS(name);	
	TString CnameB(name);	
	TString CnameSroot(name);	
	TString CnameBroot(name);	
	char CnameS2[1024];
	char CnameB2[1024];
	char CnameS2root[1024];
	char CnameB2root[1024];
	CnameS.ReplaceAll(".root",".png");
	CnameB.ReplaceAll(".root",".png");
	sprintf(CnameS2,"logN_rho/logN_rho_S_dANN%1.2f_%s",deltaANN,CnameS.Data());
	sprintf(CnameB2,"logN_rho/logN_rho_B_dANN%1.2f_%s",deltaANN,CnameB.Data());
	sprintf(CnameS2root,"logN_rho/logN_rho_S_dANN%1.2f_%s",deltaANN,CnameSroot.Data());
	sprintf(CnameB2root,"logN_rho/logN_rho_B_dANN%1.2f_%s",deltaANN,CnameBroot.Data());
	cS->SaveAs(CnameS2);
	cB->SaveAs(CnameB2);
	cS->Print(CnameS2root);
	cB->Print(CnameB2root);
//	cout<<"fine"<<endl;
	
}


void Annth(TString ifile){
	TString name(ifile);
	name.ReplaceAll("outfile/","");
	TFile* fTEST =TFile::Open(ifile.Data());
	TTree* NNTree2=(TTree*)fTEST->Get("Parameters");	           
	double out;
        double cc;	
	double rho;
	int nSi;
	NNTree2->SetBranchAddress("ANNout",&out);
	NNTree2->SetBranchAddress("cc",&cc);
	NNTree2->SetBranchAddress("rho",&rho);
	NNTree2->SetBranchAddress("#TestSig",&nSi);
	NNTree2->GetEntry(0);
	int const nSig=nSi;
//	cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
        int const ncurve2=nRHO*nCC;
	

 double* ANNSig[ncurve2];
       for (int i=0;i<ncurve2;i++) ANNSig[i]=new double[nSig];
	 //double rhoSig[20][10];
        int NSig[ncurve2];
//      int nBg=0;
        int const nBg=NNTree2->GetEntries()-nSig;
        for (int i=0;i<ncurve2;i++) {
                NSig[i]=0;
                for (int j=0;j<nSig;j++) ANNSig[i][j]=0.;
        }
        double* ANNBg[ncurve2];
	for (int i=0;i<ncurve2;i++) ANNBg[i]=new double[nBg];

        //double rhoBg[20][10];
        int NBg[ncurve2];
        for (int i=0;i<ncurve2;i++) {
                NBg[i]=0;
                for (int j=0;j<nBg;j++) ANNBg[i][j]=0.;
        }
        double ccTh[nCC];
        for (int i=0;i<nCC;i++) ccTh[i]=0.;
        double rhoTh[nRHO];
        for (int i=0;i<nRHO;i++) rhoTh[i]=0.;
        double deltacc=0.;
        double deltarho=0.;
        deltacc=0.2/nCC;
	deltarho=1./(nRHO);


        for(int n=0;n<NNTree2->GetEntries();n++){
                NNTree2->GetEntry(n);
                cout<<"rho "<<rho<<" cc "<<cc<<" out "<<out<<endl;
                for(int i=0; i<nCC;i++){
                        ccTh[i]=0.5+i*deltacc;
                        if(cc<ccTh[i]) continue;
                        for(int j=0; j<nRHO;j++){
                                rhoTh[j]=5+j*deltarho;
                                if(rho<rhoTh[j]) continue;
                                        int ni=0;
                                        if(n>nSig) {
                                                NBg[i*nRHO+j]= NBg[i*nRHO+j]+1;
                                                while(ANNBg[i*nRHO+j][ni]!=0)ni=ni+1;
                                                ANNBg[i*nRHO+j][ni]=out;
                             //                   cout<<" soglia_cc "<<ccTh[i]<<" soglia_RHO "<<rhoTh[j]<<endl;
                                                }
                                        else {
                                                NSig[i*nRHO+j]= NSig[i*nRHO+j]+1;
                                                while(ANNSig[i*nRHO+j][ni]!=0)ni=ni+1;
                                                ANNSig[i*nRHO+j][ni]=out;
                               //                 cout<<" soglia_cc "<<ccTh[i]<<" soglia_RHO "<<rhoTh[j]<<endl;
                                        }
                        }
                }
        }

  int* indexS[ncurve2];
         for (int i=0;i<ncurve2;i++) indexS[i]=new int[nSig];

        TGraph * gS[ncurve2];
        for (int y=0;y<ncurve2;y++) {
                int igS=0;
                int igS_p=0;
                gS[y]=new TGraph();
                TMath::Sort(nSig,ANNSig[y],indexS[y],false);
                for (int k=0;k<nSig;k++) {
                        int ii=indexS[y][k];
                        int yy=0;
                        if (k>0){
                                //cout<<"k "<<k<<endl;
                                int ij=indexS[y][k-1];
                                //if(ANNSig[y][ii]!=0&&ANNSig[y][ii]!=ANNSig[y][ij]) {
                                if(ANNSig[y][ii]!=0) {
                                                        yy=NSig[y]-igS;
                                                        //gS[y]->SetPoint(igS,ANNSig[y][ii],yy);
                                                        if(ANNSig[y][ii]!=ANNSig[y][ij]) gS[y]->SetPoint(igS_p++,ANNSig[y][ii],yy);
                   //                                     cout<<"igS"<<igS<<" x "<<ANNSig[y][ii]<<" y: "<<yy<<endl;
                                                        igS=igS+1;

                                                        }
                        }
                        else {
                                if(ANNSig[y][ii]!=0){
                                        yy=NSig[y]-igS;
                                        gS[y]->SetPoint(0,ANNSig[y][ii],yy);
                                        igS=igS+1;
                 //                       cout<<" x "<<ANNSig[y][ii]<<" y: "<<yy<<endl;
                                }

                        }
                }
        }



        int* indexB[ncurve2];
	for (int i=0;i<ncurve2;i++) indexB[i]=new int[nBg];

        TGraph * gB[ncurve2];
        for (int y=0;y<ncurve2;y++) {
                int igB=0;
		int igB_p=0;
                gB[y]=new TGraph();
                TMath::Sort(nBg,ANNBg[y],indexB[y],false);
               // cout<<"dopo Sort "<<y<<endl;
                for (int k=0;k<nBg;k++) {
                        int ii=indexB[y][k];
                        int yy=0;
                        if (k>0){
                                int ij=indexB[y][k-1];
                                //if(ANNBg[y][ii]!=0&&ANNBg[y][ii]!=ANNBg[y][ij]) {
                                if(ANNBg[y][ii]!=0) {
                                                yy=NBg[y]-igB;
                                                //gB[y]->SetPoint(igB,ANNBg[y][ii],yy);
                                                if(ANNBg[y][ii]!=ANNBg[y][ij]) gB[y]->SetPoint(igB_p++,ANNBg[y][ii],yy);
                                                igB=igB+1;
             //                                   cout<<"igB"<<igB<<" x "<<ANNBg[y][ii]<<" y: "<<yy<<endl;
                                                }
                        }
                      else {
                                if(ANNBg[y][ii]!=0)     {
                                        yy=NBg[y]-igB;
                                        gB[y]->SetPoint(0,ANNBg[y][ii],yy);
                                        igB=igB+1;
           //                             cout<<"igB"<<igB<<" x "<<ANNBg[y][ii]<<" y: "<<yy<<endl;
                                }
                        }
                 }
        }


  TCanvas* cS=new TCanvas("Efficiency_vs_ANN","Efficiency_vs_ANN",0,0,1200,700);
        cS->Divide(2,2);
        //cS->cd(1)->SetLogy();
        cS->cd(1);
        TMultiGraph* mg1=new TMultiGraph();
        mg1->SetTitle("cc=0.5;ANN;#Events");
        for (int h=0;h<nRHO;h++){
                gS[h]->SetLineColor(4);
                if(gS[h]->GetN()!=0) mg1->Add(gS[h]);
                }
        mg1->Draw("al");
        cS->cd(2);
     //   cS->cd(2)->SetLogy();
        TMultiGraph* mg2=new TMultiGraph();
        mg2->SetTitle("cc=0.55;ANN;#Events");
        for (int h=0;h<nRHO;h++){
                gS[nRHO+h]->SetLineColor(4);
                if(gS[nRHO+h]->GetN()!=0) mg2->Add(gS[nRHO+h]);
                }
        mg2->Draw("al");

	//cS->cd(3)->SetLogy();
	cS->cd(3);
        TMultiGraph* mg3=new TMultiGraph();
        mg3->SetTitle("cc=0.6;ANN;#Events");
        for (int h=0;h<nRHO;h++){
                gS[2*nRHO+h]->SetLineColor(4);
                if(gS[2*nRHO+h]->GetN()!=0) mg3->Add(gS[2*nRHO+h]);
                }
        mg3->Draw("al");
        //cS->cd(4)->SetLogy();
        cS->cd(4);
        TMultiGraph* mg4=new TMultiGraph();
        mg4->SetTitle("cc=0.65;ANN;#Events");
        for (int h=0;h<nRHO;h++){
                gS[3*nRHO+h]->SetLineColor(4);
                if(gS[3*nRHO+h]->GetN()!=0) mg4->Add(gS[3*nRHO+h]);
                }
        mg4->Draw("al");

	TCanvas* cB=new TCanvas("Number_vs_ANN","Number_vs_ANN",0,0,1200,700);
        cB->Divide(2,2);
        cB->cd(1)->SetLogy();
        TMultiGraph* mg1B=new TMultiGraph();
        mg1B->SetTitle("cc=0.5;ANN;#Events");
        for (int h=0;h<nRHO;h++){
                gB[h]->SetLineColor(4);
                if(gB[h]->GetN()!=0) mg1B->Add(gB[h]);
                }
        mg1B->Draw("al");
        cB->cd(2)->SetLogy();
        TMultiGraph* mg2B=new TMultiGraph();
        mg2B->SetTitle("cc=0.55;ANN;#Events");
        for (int h=0;h<nRHO;h++){
                gB[nRHO+h]->SetLineColor(4);
                if(gB[nRHO+h]->GetN()!=0) mg2B->Add(gB[nRHO+h]);
                }
        mg2B->Draw("al");
        cB->cd(3)->SetLogy();
        TMultiGraph* mg3B=new TMultiGraph();
        mg3B->SetTitle("cc=0.6;ANN;#Events");
        for (int h=0;h<nRHO;h++){
                gB[2*nRHO+h]->SetLineColor(4);
                if(gB[2*nRHO+h]->GetN()!=0) mg3B->Add(gB[2*nRHO+h]);
                }
        mg3B->Draw("al");
        cB->cd(4)->SetLogy();
        TMultiGraph* mg4B=new TMultiGraph();
        mg4B->SetTitle("cc=0.6;ANN;#Events");
        for (int h=0;h<nRHO;h++){
                gB[3*nRHO+h]->SetLineColor(4);
                if(gB[3*nRHO+h]->GetN()!=0) mg4B->Add(gB[3*nRHO+h]);
                }
        mg4B->Draw("al");


        //cout<<"dopo Draw()"<<endl;
        TString CnameS(name);
        TString CnameB(name);
        TString CnameSroot(name);
        TString CnameBroot(name);
        char CnameS2[1024];
        char CnameB2[1024];
        char CnameS2root[1024];
        char CnameB2root[1024];
        CnameS.ReplaceAll(".root",".png");
        CnameB.ReplaceAll(".root",".png");
        sprintf(CnameS2,"ANNthres/N_ANN_S_%s",CnameS.Data());
        sprintf(CnameB2,"ANNthres/N_ANN_B_%s",CnameB.Data());
        sprintf(CnameS2root,"ANNthres/N_ANN_S_%s",CnameSroot.Data());
        sprintf(CnameB2root,"ANNthres/N_ANN_B_%s",CnameBroot.Data());
        cS->SaveAs(CnameS2);
        cB->SaveAs(CnameB2);
        cS->Print(CnameS2root);
        cB->Print(CnameB2root);
        //cout<<"fine"<<endl;

//CARTELLA Annth



}


void Plots(TString ifile){
	TString name(ifile);
	name.ReplaceAll("outfile/","");
	TFile* fTEST =TFile::Open(ifile.Data());
	TTree* NNTree2=(TTree*)fTEST->Get("Parameters");	           
	double out;
        double cc;	
	double rho;
	int nSi;
	NNTree2->SetBranchAddress("ANNout",&out);
	NNTree2->SetBranchAddress("cc",&cc);
	NNTree2->SetBranchAddress("rho",&rho);
	NNTree2->SetBranchAddress("#TestSig",&nSi);
	NNTree2->GetEntry(0);
	int const nSig=nSi;
	cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
	int const nBg=NNTree2->GetEntries()-nSig;

/*	TH2D* hglitch=new TH2D("Purezza","Purezza",NCC,minCC,maxCC,NRHO,minRHO,maxRHO);
	TH2D* h2BgO_R=new TH2D("rho_cc_out","rho_cc_out",50,0.,50,50,-0.5,1.5);
	TH2D* h2SigO_R=new TH2D("rho_cc_out","rho_cc_out",50,0.,50,50,-0.5,1.5);
	TH2D* h2BgR_C=new TH2D("rho_cc_out","rho_cc_out",50,0.,1.2,50,0.,50);
	TH2D* h2SigR_C=new TH2D("rho_cc_out","rho_cc_out",50,0.,1.2,50,0.,50);
	TH2D* h2Bg=new TH2D("rho_cc_out","rho_cc_out",50,0.,1.2,50,-0.5,1.5);
	TH2D* h2Sig=new TH2D("rho_cc_out","rho_cc_out",50,0.,1.2,50,-0.5,1.5);
	h2BgO_R->SetDirectory(0);
	h2SigO_R->SetDirectory(0);
	h2BgR_C->SetDirectory(0);
	h2SigR_C->SetDirectory(0);
	h2Bg->SetDirectory(0);
	h2Sig->SetDirectory(0);
	hglitch->SetDirectory(0);
*/

	TGraph * gS[3];
	TGraph * gB[3];
	gB[0]=new TGraph();
	gS[0]=new TGraph();
	gB[1]=new TGraph();
	gS[1]=new TGraph();
	gB[2]=new TGraph();
	gS[2]=new TGraph();
	double aANNS[nSig];
	double aRHOS[nSig];
	double aCCS[nSig];
	for (int i=0;i<nSig;i++){
		aCCS[i]=0.;
		aRHOS[i]=0.;
		aANNS[i]=0.;
	}
	double aANNB[nBg];
	double aRHOB[nBg];
	double aCCB[nBg];
	for (int i=0;i<nBg;i++){
		aCCB[i]=0.;
		aRHOB[i]=0.;
		aANNB[i]=0.;
	}
//	cout<<"nSig: "<<nSig<<endl;
	for (int n = 0; n <NNTree2->GetEntries(); n++){
		 NNTree2->GetEntry(n);
		 if(n<nSig) {
			aRHOS[n]=rho;
			aANNS[n]=out;
			aCCS[n]=cc;	 
/*			gS[0]->SetPoint(n,cc,rho);
	                gS[1]->SetPoint(n,cc,out);
	                gS[2]->SetPoint(n,rho,out);
	                cout<<"Sig_graph1: x="<<cc<<" y: "<<rho<<endl;
	                cout<<"Sig_graph2: x="<<cc<<" y: "<<out<<endl;
	                cout<<"Sig_graph3: x="<<rho<<" y: "<<out<<endl;
*/
		}
		else {
			aRHOB[n-nSig]=rho;
                        aANNB[n-nSig]=out;
                        aCCB[n-nSig]=cc;
/*			gB[0]->SetPoint(n-nSig,cc,rho);
	                gB[1]->SetPoint(n-nSig,cc,out);
	                gB[2]->SetPoint(n-nSig,rho,out);
	                cout<<"Bg_graph1: x="<<cc<<" y: "<<rho<<endl;
	                cout<<"Bg_graph2: x="<<cc<<" y: "<<out<<endl;
	                cout<<"Bg_graph3: x="<<rho<<" y: "<<out<<endl;
*/
		}
	}

	int indexB[nBg];
	for (int k=0;k<nBg;k++)indexB[k]=0;
	TMath::Sort(nBg,aRHOB,indexB,false);
	int Z[(nth+1)*nth];
	for (int k=0;k<(nth+1)*nth;k++)Z[k]=0;
	double cc1th[nth];
	double ANN1th[nth+1]; 
	ANN1th[0]=0.;
	for (int k=0;k<(nth);k++){
		ANN1th[k+1]=0.;
		cc1th[k]=0.;
	}
	double dANNth=0.;
	double dccth=0.;
	dANNth=1./nth;
	dccth=0.5/nth;
	ANN1th[0]=-100.;
	for (int k=0;k<(nth);k++){
		ANN1th[k+1]=0.+dANNth*(k+1);
		cc1th[k]=0.5+dccth*k;
	 }

	/*for (int a=0;a<nBg;a++){
	  if(aRHOB[a]>=RHOth){
		for (int b=0;b<nth;b++){
		  if(aCCB[a]>=cc1th[b]){
			for(int c=0;c<nth+1;c++){
				if(aANNB[a]>=ANN1th[c]) Z[b*nth+c]=Z[b*nth+c]+1;
				if(c=10) cout<<" ANN "<<ANN1th[c]<<" Zcount: "<<Z[b*nth+c]<<" aNNB "<<aANNB[a]<<" a "<<a<<" b "<<b<<endl;
			}
		  }
		}
	  }
	}*/
	for (int a=0;a<nBg;a++){
	  cout<<"tutti: a "<<a<<" rho "<<aRHOB[a]<<" aNNB "<<aANNB[a]<<" cc "<<aCCB[a]<<endl;
	  if(aRHOB[a]<RHOth)continue;
	  cout<<"dopo rho: a "<<a<<" rho "<<aRHOB[a]<<" aNNB "<<aANNB[a]<<" cc "<<aCCB[a]<<endl;
	  for (int b=0;b<nth;b++){
	     if(aCCB[a]<cc1th[b])continue;
	  cout<<"dopo cc: "<<cc1th[b]<<"a "<<a<<" rho "<<aRHOB[a]<<" aNNB "<<aANNB[a]<<" cc "<<aCCB[a]<<endl;
		for(int c=0;c<nth+1;c++){
			if(aANNB[a]<ANN1th[c]) continue;
	  cout<<"dopo ANN:"<<ANN1th[c]<<" a "<<a<<" rho "<<aRHOB[a]<<" aNNB "<<aANNB[a]<<" cc "<<aCCB[a]<<endl;
			Z[b*(nth+1)+c]=Z[b*(nth+1)+c]+1;
	  cout<<"Zindex: "<<b*(nth+1)+c<<" Z: "<<Z[b*(nth+1)+c]<<endl;

			//if(c=10) cout<<" ANN "<<ANN1th[c]<<" Zcount: "<<Z[b*(nth+1)+c]<<" aNNB "<<aANNB[a]<<" a "<<a<<" b "<<b<<endl;
			}
		  
		}
	  
	}
			
	//TH2D* hglitch=new TH2D("Survived_glitches(RHO_cut:5.5)","Survived_glitches(RHO_cut:5.5)",nth,0.5,1,nth+1,-1,1);
	double bins[nth+1];
	for (int d=0; d<nth+1;d++) bins[d]=0.;
	//bins[0]=0.;
	for (int d=0; d<nth+1;d++) bins[d+1]=0.+(d+1)*dANNth;
	TH2D* hglitch=new TH2D("Survived_glitches(RHO_cut:5.5)","Survived_glitches(RHO_cut:5.5)",nth,0.5,1,nth+1,bins);
//	TH2D* hglitch=new TH2D("Survived_glitches(RHO_cut:5.5)","Survived_glitches(RHO_cut:5.5)",nth,0.5,1,nth+1,0.,1.);
	hglitch->SetDirectory(0);
	//if(tot[NRHO*ii+ji]>0) hglitch->Fill(cc1[ji]+deltacc/2,rho1[ii]+deltarho/2,pur[NRHO*ii+ji]);
	for (int b=0;b<nth;b++) for(int c=0;c<nth+1;c++) {
			//hglitch->Fill(cc1th[b]+dccth/2.,ANN1th[c]+dANNth/2.,Z[(nth+1)*b+c]);
			hglitch->Fill(cc1th[b]+dccth/2.,bins[c]+dANNth/2.,Z[(nth+1)*b+c]);
			cout<<" Zindex: "<<(nth+1)*b+c<< " x: "<<cc1th[b]<<" y: "<<ANN1th[c]<<" z: "<<Z[(nth+1)*b+c]<<endl;
			}
	
	gB[0]=new TGraph(nBg,aCCB,aRHOB);
	gS[0]=new TGraph(nSig,aCCS,aRHOS);
	gB[1]=new TGraph(nBg,aCCB,aANNB);
	gS[1]=new TGraph(nSig,aCCS,aANNS);
	gB[2]=new TGraph(nBg,aRHOB,aANNB);
	gS[2]=new TGraph(nSig,aRHOS,aANNS);

	for(int i=0;i<3;i++){
	gS[i]->SetMarkerColor(2);
	gB[i]->SetMarkerColor(4);
	gS[i]->SetMarkerStyle(6);
	gB[i]->SetMarkerStyle(7);
	}

        TCanvas* c=new TCanvas("Plots","Plots",0,0,1200,700);
	c->Divide(2,2);
	c->cd(1);
        TMultiGraph* mg1=new TMultiGraph();
	mg1->SetTitle("cc_rho;cc;rho");
	if(gB[0]->GetN()!=0) mg1->Add(gB[0]);
	if(gS[0]->GetN()!=0) mg1->Add(gS[0]);
	mg1->Draw("ap");
	c->cd(2);
        TMultiGraph* mg2=new TMultiGraph();
	mg2->SetTitle("cc_ANNoutput;cc;ANNoutput");
	if(gB[1]->GetN()!=0) mg2->Add(gB[1]);
	if(gS[1]->GetN()!=0) mg2->Add(gS[1]);
	mg2->Draw("ap");
	c->cd(3);
        TMultiGraph* mg3=new TMultiGraph();
	mg3->SetTitle("rho_ANNoutput;rho;ANNoutput");
	if(gB[2]->GetN()!=0) mg3->Add(gB[2]);
	if(gS[2]->GetN()!=0) mg3->Add(gS[2]);
	mg3->Draw("ap");
	c->cd(4);
	TText* text=new TText(0.37,0.0,"no cuts on ANN");
	hglitch->SetStats(0);
	hglitch->GetXaxis()->SetTitle("cc");
	hglitch->GetYaxis()->SetTitle("ANNoutput");
	hglitch->GetZaxis()->SetTitle("count");
	hglitch->Draw("colz");
	text->Draw();
        gPad->SetLogz();

   TString Cname(name);
   TString Cnameroot(name);
   char Cname2[1024];
   char Cname2root[1024];
   Cname.ReplaceAll(".root",".png");
   sprintf(Cname2,"CC_RHO_ANN_Plots/CC_RHO_ANN_Plots_%s",Cname.Data());
   sprintf(Cname2root,"CC_RHO_ANN_Plots/CC_RHO_ANN_Plots_%s",Cnameroot.Data());
   c->SaveAs(Cname2);
   c->Print(Cname2root);

        TCanvas* c2=new TCanvas("Plots_Bkg_on_Sig","Plots_Bkg_on_Sig",0,0,1200,700);
	c2->Divide(2,2);
	c2->cd(1);
	gS[0]->Draw("ap");
	gB[0]->Draw("p,same");
	
	c2->cd(2);
	gS[1]->Draw("ap");
	gB[1]->Draw("p,same");
	c2->cd(3);
	gS[2]->Draw("ap");
	gB[2]->Draw("p,same");
	c2->cd(4);
	TText* text2=new TText(0.37,0.0,"no cuts on ANN");
	hglitch->SetStats(0);
	hglitch->GetXaxis()->SetTitle("cc");
	hglitch->GetYaxis()->SetTitle("ANNoutput");
	hglitch->GetZaxis()->SetTitle("count");
	hglitch->Draw("colz");
	text2->Draw();
        gPad->SetLogz();

   TString Cname_2(name);
   TString Cname_2root(name);
   char Cname2_2[1024];
   char Cname2_2root[1024];
   Cname_2.ReplaceAll(".root",".png");
   sprintf(Cname2_2,"CC_RHO_ANN_Plots/CC_RHO_ANN_Plots_BoS_%s",Cname_2.Data());
   sprintf(Cname2_2root,"CC_RHO_ANN_Plots/CC_RHO_ANN_Plots_BoS_%s",Cname_2root.Data());
   c2->SaveAs(Cname2_2);
   c2->Print(Cname2_2root);

//CARTELLA CCi_RHO_ANN_Plots 
}

