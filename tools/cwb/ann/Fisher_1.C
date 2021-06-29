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
//#define deltaANN 0.1
#define nCC 4
#define NPIX 20
#define RHOth 5.5
#define nth 10
#define rhomin 6
#define ccmin 0.6
#define iMc 3
using namespace std;
void Fisher_ex(TString ifile, int nFs, int nFb, int Fs, int Fb);
void graph_Fish(TString ifile);
//void Annth(TString ifile);
void PlotsAv_cc(TString ifile);
void PlotsAv_Mc(TString ifile);


//uf !=0 means the Test file of events is the same used for train the network, ni==0 select the event used for the training, in this way we test the events alway for NNn-1 network, if ni!=0 also the events not used for the training are considered and there the average in made on all the nNN network
void Fisher_1(TString NN_FILE,TString NN_FILE2, TString NN_FILE3, TString NN_FILE4, TString NN_FILE5, TString NN_FILE6, TString NN_FILE7, TString NN_FILE8, TString TEST_FILE, TString NOMEtot_S, int nFs, int nFb, int Fs=0, int Fb=0,int TS=0, int TB=0, int s=0, int b=0, int uf=1, int ni=0){
if(uf==0) ni=1; 
int n_NN=0;   
if (NN_FILE.CompareTo("")){
        	       n_NN=n_NN+1; 
        }
if (NN_FILE2.CompareTo("")){
        	       n_NN=n_NN+1; 
        }
if (NN_FILE3.CompareTo("")){
        	       n_NN=n_NN+1; 
        }
if (NN_FILE4.CompareTo("")){
        	       n_NN=n_NN+1; 
        }
if (NN_FILE5.CompareTo("")){
        	       n_NN=n_NN+1; 
        }
if (NN_FILE6.CompareTo("")){
        	       n_NN=n_NN+1; 
        }
if (NN_FILE7.CompareTo("")){
        	       n_NN=n_NN+1; 
        }
if (NN_FILE8.CompareTo("")){
        	       n_NN=n_NN+1; 
        }
int const nNN=n_NN;

//   int p=0;
   char NNi[nNN][1024];
   char NNi2[nNN][1024];
   sprintf(NNi2[0],"%s",NN_FILE.Data());
   if(nNN>1) sprintf(NNi2[1],"%s",NN_FILE2.Data());
   if(nNN>2) sprintf(NNi2[2],"%s",NN_FILE3.Data());
   if(nNN>3) sprintf(NNi2[3],"%s",NN_FILE4.Data());
   if(nNN>4) sprintf(NNi2[4],"%s",NN_FILE5.Data());
   if(nNN>5) sprintf(NNi2[5],"%s",NN_FILE6.Data());
   if(nNN>6) sprintf(NNi2[6],"%s",NN_FILE7.Data());
   if(nNN>7) sprintf(NNi2[7],"%s",NN_FILE8.Data());
   for (int u=0;u<nNN;u++){
   int p=0;
   while (NNi2[u][p]){
	//cout<<NNi2[p]<<endl;
	if (NNi2[u][p]=='N'){
//	cout<<p<<endl;
                        if((NNi2[u][p+1]=='N')&&(NNi2[u][p+2]=='\/')) {
//					cout<<" dentro if"<<endl;
					int hh=p+3;
					while (NNi2[u][hh]!='\0'){NNi[u][hh-p-3]=NNi2[u][hh];hh=hh+1;}
					break;
					}
	}
	p=p+1;
   }
}
   int q=0;
   char Filei[1024];
   char Filei2[1024];
   sprintf(Filei2,"%s",TEST_FILE.Data());
 while (Filei2[q]){
 if(Filei2[q]=='n'&&Filei2[q+1]=='n'&&Filei2[q+2]=='T'&&Filei2[q+3]=='R'&&Filei2[q+4]=='E'&&(Filei2[q+5]=='E')&&(Filei2[q+6]=='\/')) {
               int hh=q+7;
               while (Filei2[hh]!='\0') {Filei[hh-q-7]=Filei2[hh];hh=hh+1;
					}
	       for (int h0=hh-q-7;h0<1024;h0++) Filei[h0]='\0';
               break;
               }
        q=q+1;
   }

   cout<<Filei<<" original String: "<<Filei2<<endl;
   

   TString NNi0[nNN];
   NNi0[0]=NNi[0];
   if(nNN>1) NNi0[1]=NNi[1];
   if(nNN>2) NNi0[2]=NNi[2];
   if(nNN>3) NNi0[3]=NNi[3];
   if(nNN>4) NNi0[4]=NNi[4];
   if(nNN>5) NNi0[5]=NNi[5];
   if(nNN>6) NNi0[6]=NNi[6];
   if(nNN>7) NNi0[7]=NNi[7];
   TString Filei0(Filei);
cout<<"NN"<<endl;
int min_NNb=5000000;
   for (int u=0;u<nNN;u++){
   	NNi0[u].ReplaceAll(".root","");
   }
   Filei0.ReplaceAll(".root","");
    
   TFile *fnet[nNN];
   TMultiLayerPerceptron* mlp[nNN];
   TTree* infot[nNN];
   int NNs[nNN];
   int NNb[nNN];
   int NNnTS[nNN];
   int NNnTB[nNN];
   for (int u=0;u<nNN;u++){
    fnet[u]=TFile::Open(NNi2[u]);
    cout<<NNi2[u]<<endl;
    cout<<"dopo file"<<endl;
    mlp[u]=(TMultiLayerPerceptron*)fnet[u]->Get("TMultiLayerPerceptron");
    
    cout<<"dopo TMLP"<<endl;
    if(mlp[u]==NULL) {cout << "Error getting mlp" << endl;exit(1);}
    infot[u]=(TTree*)fnet[u]->Get("info");
   infot[u]->SetBranchAddress("Rand_start_Sig",&NNs[u]);
   infot[u]->SetBranchAddress("Rand_start_Bg",&NNb[u]);
   infot[u]->SetBranchAddress("#trainSig",&NNnTS[u]);
   infot[u]->SetBranchAddress("#trainBg",&NNnTB[u]);
   infot[u]->GetEntry(0);
   if(NNb[u]<min_NNb) min_NNb=NNb[u];
   cout<<"n "<<u<<" b "<<NNb[u]<<" min_NNb "<<min_NNb<<endl;

   cout<<"da Tree"<<endl;
   cout<<"n "<<u<<" s  "<<NNs[u]<<" b "<<NNb[u]<<" s fin  "<<NNs[u]+NNnTS[u]<<" b fin "<<NNb[u]+NNnTB[u]<<" nTS "<<NNnTS[u]<<" nTB "<<NNnTB[u]<<endl;
 }
if (uf!=0&&ni==0&&(b<min_NNb)&&(TS!=0||TB!=0)){ b=min_NNb;cout<<" change in b value to have BKG events: "<<b<<endl;}
//exit(0);
cout<<"Def. tree e mlp"<<endl;
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
   sprintf(NOMEtot,"%s",NOMEtot_S.Data());
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
   for (int jj=0; jj<nINP;jj++) x[jj]=0.;
   char ilabel[nINP][16];
  
      //costruzione una foglia  per ogni input
   for(int i=0;i<nINP;i++) {
       sprintf(ilabel[i],"x%i",i+1);
       NNTree->SetBranchAddress(ilabel[i], &x[i]);
   }
  
char ofile[1024];
sprintf(ofile,"average_file/%s.root",NOMEtot);
TFile*f =new TFile(ofile,"RECREATE");
   TTree* NNTree2=new TTree("Parameters","Parameters");
  NNTree2->SetDirectory(f);
   double average=0.;
   double out[nNN];
   for(int y=0; y<nNN; y++) out[y]=0.;
   double Mc=0.;
   double NNcc=0.;
   double NNrho=0.;
   double std=0.;
  NNTree2->Branch("Average",&average,"Average/D");
  NNTree2->Branch("StandardDevaition",&std,"StandardDeviation/D");
  NNTree2->Branch("cc",&NNcc,"cc/D");
  NNTree2->Branch("Mchirp",&Mc,"Mchirp/D");
  NNTree2->Branch("rho",&NNrho,"rho/D");
  for(int u=0;u<nNN;u++){
    char NNoutl[512];
    char NNoutl2[512];
    char NNnamel[512];
    char NNnamel2[512];
    sprintf(NNoutl,"NNout%i",u);
    sprintf(NNoutl2,"NNout%i/D",u);
    sprintf(NNnamel,"NNname%i",u);
    sprintf(NNnamel2,"NNname%i/C",u);
    NNTree2->Branch(NNoutl,&out[u],NNoutl2);
    NNTree2->Branch(NNnamel,&NNi[u],NNnamel2);}
  char Testf[1024];
  NNTree2->Branch("TestFile",&Testf,"TestFile/C");
  sprintf(Testf,"%s",TEST_FILE.Data());
  int nTestS=0;
  NNTree2->Branch("#TestSig",&nTestS,"#TestSig/I");
  cout<<"nTestS: "<<nTestS<<" TS: "<<TS<<endl; 
  cout<<"dopo def tree"<<endl; 
int scount=0;
  if(uf==0) nTestS=TS;
  else {
for(int n=s;n<s+TS;n++) {
    int countNN=0;
    for(int u=0;u<nNN;u++){
     if (uf!=0&&n>=NNs[u]&&n<(NNs[u]+NNnTS[u])) countNN=countNN+1;
	}
        //if(nNN==1){ni=1;cout<<"out==Averege:only 1NN is considered"<<endl;}
	if(countNN==1&&ni==0) scount=scount+1;
	if(countNN<2&&ni!=0) scount=scount+1;
        if(countNN>1){cout<<"Error: training non independent"; exit(0);}
	}
nTestS=scount;
}

cout<<"s test "<<s<< " s+TS "<<s+TS<<" b "<<b<<" b+ TB "<<b+TB<<" nTestS "<<nTestS<<endl;
//exit(0);
   double params[nINP];
   int sig_05=0;
   for (int i=0; i<nINP; i++) params[i]=0.;
     //for(int n=0;n<entryS.GetEntries();n++) {

  int S_eff=0; 
  int TD=0; 
  for(int n=s;n<s+TS;n++) {
     average=0.;
     int indNN=-1;
     int countNN=0;
     for(int u=0;u<nNN;u++){
	     if (uf!=0&&n>=NNs[u]&&n<(NNs[u]+NNnTS[u])) { indNN=u;countNN=countNN+1;}	
	}
     if (countNN>1){cout<<"Error: training non independent"; exit(0);}
    if(nNN==1){cout<<"out==Averege:only 1NN is considered"<<endl;}
    if(ni==0&&countNN==0)continue;
     NNTree->GetEntry(n);
     signal.GetEntry(n);
     NNcc=(double)signal.netcc[1];
     NNrho=(double)signal.rho[0];
     Mc=(double)signal.chirp[iMc];
     for (int i=0; i<nINP; i++){
        params[i]=x[i];	
	}
    for (int u=0;u<nNN;u++) {
     double output=mlp[u]->Evaluate(0,params);
     cout<<output<<endl;
     out[u]=output;
    }
    for (int u=0;u<nNN;u++){
        if((u!=indNN&&ni==0)||ni!=0) {average=out[u]+average;std+=out[u]*out[u];}
      }
    //if(nNN!=1&&ni==0) {average=average/(nNN-countNN); if((nNN-1-countNN)!=0) std=pow((std/(nNN-1-countNN)-average*average),0.5); else std=pow((std/(nNN-countNN)-average*average)0,5);}
      if(nNN!=1&&ni==0) {average=average/(nNN-countNN); if((nNN-1-countNN)!=0)std=pow((std/(nNN-1-countNN)-average*average),0.5);else std=pow((std/(nNN-countNN)-average*average),0.5);}
    if(nNN!=1&&ni!=0) {average=average/(nNN); if((nNN-1)!=0)std=pow((std/(nNN-1)-average*average),0.5);else std=pow((std/nNN-average*average),0.5);}
    cout<<"average: "<<average<<endl;
  //  cout<<"Mc "<<Mc<<endl;
    if (average<0.6) TD=TD+1;
    //if (outbis<0.5) sig_05=sig_05+1;
    //if (outbis>0.5&&outbis<0.6) if(out<0.1) sig_05=sig_05+1;
   // if (outbis<0.5 || out<0.5) sig_05=sig_05+1;
    //if (outbis<0.6) sig_05=sig_05+1;
   // cout<<"rho: "<<signal.rho[0]<<" cc "<<signal.netcc[1]<<" out: "<<out<<" outbis "<<outbis<<endl; 
     NNTree2->Fill();
     S_eff=S_eff+1;
    }
int B_eff=0; 
int FA=0;      
cout<<"riempito sig"<<endl;
int bg_05=0;
cout<<"b "<<b<<" TB "<<TB<<endl;
//for(int n=b;n<b+TB;n++) {
//for(int u=0;u<nNN;u++){
  //                      cout<<" u "<<u<<" NN b[u] "<<NNb[u]<<" NNnTB[u] "<<NNnTB[u]<<" NNb[u]+NNnTB[u]"<<NNb[u]+NNnTB[u]<<endl;
   //              }
//}
//exit(0);
cout<<" prima ciclo fro bg"<<endl;
   //for(int n=sig_entries+b;n<sig_entries+b+TB;n++) {
   for(int n=b;n<b+TB;n++) {
      cout<<n<<endl;
      average=0.;
      int indNN=-1;
      int countNN=0;
      for(int u=0;u<nNN;u++){
      cout<<" uf "<<uf<<" n "<<n<<" u "<<u<<" NN b[u] "<<NNb[u]<<" (primo non preso) NNb[u]+NNnTB[u]"<<NNb[u]+NNnTB[u]<<endl;
             if (uf!=0&&n>=NNb[u]&&n<(NNb[u]+NNnTB[u])) { indNN=u;countNN=countNN+1;
	//		cout<<"n "<<n<<" u "<<u<<" NN b[u] "<<NNb[u]<<" (primo non preso)NNb[u]+NNnTB[u]"<<NNb[u]+NNnTB[u]<<endl;
		 }
        }
     if (countNN>1){cout<<"Error: training non independent"; exit(0);}
     if(nNN==1){cout<<"out==Averege:only 1NN is considered"<<endl;}
     if(ni==0&&countNN==0){
        cout<<"countNN"<<countNN<<endl;
	continue;}
     //if (uf!=0&&(n>=NNb&&n<=(NNb+NNnTB) || n>=NNbbis&&n<=(NNbbis+NNnTBbis))) continue;	
     NNTree->GetEntry(n);
     cout<<"n: "<<n<<"Bg index"<<(n-sig_entries)<<endl;
     background.GetEntry(n-sig_entries);
     NNcc=(double)background.netcc[1];
     NNrho=(double)background.rho[0];
     Mc=(double)background.chirp[iMc];
     for (int i=0; i<nINP; i++){
        //x[i]=(*frameS)[i];
        params[i]=x[i];
        }

	for(int u=0;u<nNN;u++) {
	     double output=mlp[u]->Evaluate(0,params);
    	     cout<<output<<endl;
	     out[u]=output;
	    }
    for (int u=0;u<nNN;u++){
           if((u!=indNN&&ni==0)||ni!=0) { average=out[u]+average;std=out[u]*out[u];}
      }
    if(nNN!=1&&ni==0) {cout<<average<<endl;average=average/(nNN-countNN);
//cout<<"ni==0"<<average<<endl; if((nNN-1-countNN)!=0)std=pow((std/(nNN-1-countNN)-average*average),0.5);else std=pow((std/(nNN-countNN)-average*average)0,5);}
cout<<"ni==0"<<average<<endl; if((nNN-1-countNN)!=0)std=pow((std/(nNN-1-countNN)-average*average),0.5);else std=pow((std/(nNN-countNN)-average*average),0.5);}

    //if(nNN!=1&&ni!=0) {cout<<average<<endl;average=average/(nNN);if((nNN-1)!=0)std=pow((std/(nNN-1)-average*average),0.5);else std=pow((std/(nNN)-average*average)0,5);cout<<"ni!=0"<<average<<endl;}
    if(nNN!=1&&ni!=0) {cout<<average<<endl;average=average/(nNN);if((nNN-1)!=0)std=pow((std/(nNN-1)-average*average),0.5);else std=pow((std/nNN-average*average),0.5); cout<<"ni!=0"<<average<<endl;}
    cout<<"average: "<<average<<" nNN-1-countNN "<<nNN-1-countNN<<endl;
    if(average>0.6) FA=FA+1;
/*     d uble output=mlp->Evaluate(0,params);
     out=output;
     double outputbis=mlpbis->Evaluate(0,params);
     outbis=outputbis;*/
//	if(outbis>0.6) bg_05=bg_05+1;
//	if(outbis<=0.6&&outbis>0.5) if(out>0.1)bg_05=bg_05+1;

//	if(outbis>0.5 && out>0.5) bg_05=bg_05+1;
	//if(outbis>0.6) bg_05=bg_05+1;
    
//    cout<<"rho: "<<background.rho[0]<<" cc "<<background.netcc[1]<<" out: "<<out<<endl; 
     NNTree2->Fill();
      B_eff=B_eff+1;
     }
cout<<"riempito bg"<<endl;
//TFile*f0 =new TFile(ofile,"RECREATE");
//NNTree2->SetDirectory(f0);
NNTree2->Write();
f->Close();
//f0->Close();
cout<<"chiuso file"<<endl;


//graph(ofile.Data());
//graph(ofile);
cout<<ofile<<endl;
cout<<"dopo richiamo funzione"<<endl;
Fisher_ex(ofile,nFs,nFb,Fs,Fb);
//Annth(ofile);
//exit(0);
graph_Fish(ofile);
//PlotsAv_cc(ofile);
//PlotsAv_Mc(ofile);
cout<<" S_eff "<<S_eff<<" B_eff "<<B_eff<<endl;
cout<<" FA>0.6 "<<FA<<" TD<0.6 "<<TD<<" FA/tot_BG "<<FA/B_eff<<" TD/tot_S "<<TD/S_eff<<endl;
//cout<<" Sig with out<06: "<<sig_05<<endl;
//cout<<" Bg with double threshold "<<bg_05<<endl;
}


void Fisher_ex(TString ifile, int nFs, int nFb, int Fs, int Fb){


TString name(ifile);
        name.ReplaceAll("average_file/","");
        TFile* fTEST =TFile::Open(ifile.Data());
        TTree* NNTree2=(TTree*)fTEST->Get("Parameters");
        double av;
        double cc;
        int nSi;
        NNTree2->SetBranchAddress("Average",&av);
        NNTree2->SetBranchAddress("cc",&cc);
        NNTree2->SetBranchAddress("#TestSig",&nSi);
        NNTree2->GetEntry(0);

//INIZIO DISCRIMINAZIONE FISHER

double qS_o=0.;
double qS_c=0.;
double sS_co=0.;
double mS_o=0.;
double mS_c=0.;
double mB_o=0.;
double mB_c=0.;
double qB_c=0.;
double qB_o=0.;
double sB_co=0.;

if (nFs+Fs>nSi){ Fs=0; if(nFs>nSi) nFs=nSi; cout<<" change: Fs "<<Fs<<" nFs "<<nFs<<endl;}
if(Fb==0)Fb=nSi;
if (Fb<nSi) {Fb=nSi+Fb; cout<<" change: Fb "<<Fb<<endl;}

for(int n=Fs; n<Fs+nFs;n++){
	 NNTree2->GetEntry(n);
	 mS_o+=av;
	 mS_c+=cc;
	 qS_o+=av*av;
	 qS_c+=cc*cc;
	 sS_co=av*cc;

}
mS_o=mS_o/nFs;
mS_c=mS_c/nFs;

int const TOTent=NNTree2->GetEntries();
int const nBg=NNTree2->GetEntries()-nSi;
cout<<"nSi "<<nSi<<" nBg "<<nBg<<" entries "<<TOTent<<endl;

if (nFb+Fb>TOTent){ Fb=nSi; if(nFs>nBg) nFb=nBg; cout<<" change: Fb "<<Fb<<" nFb "<<nFb<<endl;}

for(int n=Fb; n<Fb+nFb;n++){
         NNTree2->GetEntry(n);
         mB_o+=av;
         mB_c+=cc;
         qB_o+=av*av;
         qB_c+=cc*cc;
         sB_co=av*cc;

}
mB_o=mB_o/nFb;
mB_c=mB_c/nFb;

//ELEMENTS OF WITHIN MATRIX 
double a_s=0.;
double b_s=0.;
double d_s=0.;
double a_b=0.;
double b_b=0.;
double d_b=0.;

a_b=qB_c-nFb*(mB_c*mB_c);
d_b=qB_o-nFb*(mB_o*mB_o);
b_b=sB_co-nFb*(mB_c*mB_o);

cout<<"a_b "<<a_b<<endl;
cout<<"d_b "<<d_b<<endl;
cout<<"a_s "<<a_s<<endl;
cout<<"d_s "<<d_s<<endl;
a_s=qS_c-nFs*(mS_c*mS_c);
d_s=qS_o-nFs*(mS_o*mS_o);
b_s=sS_co-nFs*(mS_c*mS_o);


double a=0.;
double b=0.;
double d=0.;

a=a_b+a_s;
b=b_b+b_s;
d=d_b+d_s;

//ELEMENT OF THE INVERSE WITHIN MATRIX

double det=0.;
double a_i=0.;
double b_i=0.;
double d_i=0.;

det=a*d-b*b;
cout<<"ca standard deviation cc  bg "<<pow(a_b/nFb,0.5)<<" ca standard deviation out  bg "<<pow(d_b/nFb,0.5)<<" det "<<det<<" media cc  bg "<<mB_c<<" media out  bg "<<mB_o<<endl;
cout<<"ca standard deviation cc sig "<<pow(a_s/nFs,0.5)<<" ca standard deviation out sig "<<pow(d_s/nFs,0.5)<<" det "<<det<<" media cc sig "<<mS_c<<" media out sig "<<mS_o<<endl;
if(det==0){cout<<"Problem: det==0"<<endl;exit(0);}
a_i=d/det;
b_i=-b/det;
d_i=a/det;

//DIFFERENCES OF THE AVERAGES
double dm_o=0.;
double dm_c=0.;
dm_c=mS_c-mB_c;
dm_o=mS_o-mB_o;
cout<<" a_i "<<a_i<<" d_i "<<d_i<<" b_i  "<<b_i<<" dm_c "<<dm_c<<" dm_o "<<dm_o<<endl;
//w1 and w2
double w_c=0.;
double w_o=0.;
w_c=a_i*dm_c+b_i*dm_o;
w_o=b_i*dm_c+d_i*dm_o;

//FINDING THE THRESHOLD
TGraph *ms=new TGraph();
TGraph *mb=new TGraph();
TGraph *pg=new TGraph();
TGraph *g=new TGraph();
TGraph *g0=new TGraph();
TGraph *g1=new TGraph();
int count=0;
double w0=0.;
for(int n=Fs; n<Fs+nFs;n++){
	NNTree2->GetEntry(n);
//	w0*=-(w_c*cc+w_o*av);
	cout<<" sig"<<n<<" "<<w_c*cc+w_o*av<<endl;
	g0->SetPoint(count,cc,av);
	count=count+1;
	
} 
cout<<"count sig "<<count<<endl;
count=0;
/*for(int n=Fb; n<Fb+nFb;n++){
        NNTree2->GetEntry(n);
  //      w0+=-(w_c*cc+w_o*av);
count=0;*/
for(int n=Fb; n<Fb+nFb;n++){
        NNTree2->GetEntry(n);
  //      w0+=-(w_c*cc+w_o*av);
	cout<<" bg"<<n<<" "<<w_c*cc+w_o*av<<endl;
	g1->SetPoint(count,cc,av);
	count=count+1;
}
double mb0=w_c*(mB_c)+w_o*(mB_o);
mb0=mb0/2;

double ms0=w_c*(mS_c)+w_o*(mS_o);
ms0=ms0/2;

w0=w_c*(mS_c+mB_c)+w_o*(mS_o+mB_o);
cout<<" mS_c+mB_c "<<mS_c+mB_c<<" mS_o+mB_o "<<mS_o+mB_o<<endl;
w0=w0/2.;
cout<<" count bg "<<count<<endl;
double x_ms=0.;
double y_ms=0.;
double x_mb=0.;
double y_mb=0.;
double x_thres=0.;
double y_thres=0.;
/*x_mb=pow(mb0*mb0/(1+(w_c*w_c/(w_o*w_o))),0.5);
y_mb=-w_c/w_o*x_mb;
x_ms=pow(ms0*ms0/(1+(w_c*w_c/(w_o*w_o))),0.5);
y_ms=-w_c/w_o*x_ms;*/
x_thres=pow(w0*w0/(1+(w_c*w_c/(w_o*w_o))),0.5);
y_thres=-w_c/w_o*x_thres;
pg->SetPoint(0,x_thres,y_thres);
TF1 retta_ms("retta_p","[1]*x+[0]",-0.5,1.);
TF1 retta_mb("retta","[1]*x+[0]",-0.5,1.);
TF1 retta_p("retta_p","[1]*x+[0]",-0.5,1.);
TF1 retta("retta","[1]*x+[0]",-0.5,1.);
retta.SetParameter(0,0.);
double c_ang=0.;
c_ang=-w_c/w_o;
retta.SetParameter(1,c_ang);
double c_ang_p=-1./c_ang;
retta_p.SetParameter(1,c_ang_p);
double q_p=0.;
q_p=(mS_o+mB_o)/2.-c_ang_p*(mS_c+mB_c)/2.;
//q_p=y_thres-c_ang_p*x_thres;
retta_p.SetParameter(0,q_p);
retta_ms.SetParameter(1,c_ang_p);
double q_ms=0.;
//q_ms=y_ms-c_ang_p*x_ms;
q_ms=mS_o-c_ang_p*mS_c;
retta_ms.SetParameter(0,q_ms);
retta_mb.SetParameter(1,c_ang_p);
double q_mb=0.;
//q_mb=y_mb-c_ang_p*x_mb;
q_mb=mB_o-c_ang_p*mB_c;
retta_mb.SetParameter(0,q_mb);
cout<<"q_mb "<<q_mb<<" q_ms "<<q_ms<<" mS_o "<<mS_o<<" mS_c"<<mS_c<<" mB_o "<<mB_o<<" mB_c"<<mB_c<<endl;
cout<<" wc "<<w_c<<" w_o "<<w_o<<" coef ang "<<c_ang<<" x_thres "<<x_thres<<" y_thres "<<y_thres<<" w0 "<<w0<<" q_p "<<q_p<<" c_ang_p"<<c_ang_p<<endl;
TCanvas *canv=new TCanvas("reg_wErr_woErr_media","reg_wErr_woErr_media",0,0,600,400);
g1->SetMarkerStyle(6);
g1->SetMarkerColor(4);
g0->SetMarkerStyle(7);
g0->SetMarkerColor(2);
g->SetPoint(0,1,1.3);
g->SetPoint(1,-0.3,-0.3);
g->SetMarkerColor(10);
g->Draw("ap");
g0->Draw("p,same");
g1->Draw("p,same");
pg->SetMarkerStyle(29);
pg->SetMarkerColor(1);
retta_ms.SetLineColor(3);
retta_mb.SetLineColor(3);
retta_p.SetLineColor(6);
retta_ms.Draw("same");
retta_mb.Draw("same");
retta_p.Draw("same");
retta.SetLineColor(1);
retta.Draw("same");
//retta.Draw("al");
char nomec[1024];
char nomecroot[1024];
sprintf(nomec,"Fisher_png/%s_.png",name.Data());
sprintf(nomecroot,"Fisher_png/%s_.root",name.Data());
canv->SaveAs(nomec);
canv->SaveAs(nomecroot);
canv->Close();
fTEST->Close();

TString NOMEtot(ifile);
NOMEtot.ReplaceAll("average_file/","");
NOMEtot.ReplaceAll(".root","");
char namefF[1024];
sprintf(namefF,"Fisher_file/%s_Fisher.root",NOMEtot.Data());
cout<<" namefF "<<namefF<<" NOMEtot "<<NOMEtot.Data()<<endl;
TFile* Fishf =new TFile(namefF,"RECREATE");
  TTree* TreeFis=new TTree("Fisher_Discriminant","Fisher_Discriminant");
  TreeFis->SetDirectory(Fishf);
	  double wc=0.;
	  double wo=0.;
	  TreeFis->Branch("Fish_coef_cc",&wc,"Fish_coef_cc/D");
	  TreeFis->Branch("Fish_coef_ANNout",&wo,"Fish_coef_ANNout/D");
	  TreeFis->Branch("nFs",&nFs,"nFs/I");
	  TreeFis->Branch("nFb",&nFb,"nFb/I");
	  TreeFis->Branch("Fb_start",&Fb,"Fb_start/I");
	  TreeFis->Branch("Fs_start",&Fs,"Fs_start/I");
	  char nomeif[2048];
          sprintf(nomeif,"%s",ifile.Data());
	  TreeFis->Branch("File_name",&nomeif,"File_name/C");

	wo=w_o;
	wc=w_c;  
	TreeFis->Fill();
	cout<<" wc "<<wc<<" wo "<<wo<<endl;
	cout<<" entries "<<TreeFis->GetEntries()<<endl;
	TreeFis->Write();

Fishf->Close();

}


void graph_Fish(TString ifile){
	TString name(ifile);
	name.ReplaceAll("average_file/","Fisher_file/");
	TFile* fTEST =TFile::Open(ifile.Data());
	TTree* NNTree2=(TTree*)fTEST->Get("Parameters");	       
        
	double w_o;
	double w_c;
	int nFs;
	int nFb;
	int Fs;
	int Fb;

        TString namefF(name);
        namefF.ReplaceAll(".root","_Fisher.root");
        TFile* fF =TFile::Open(namefF.Data());
	TTree* FTree= (TTree*)fF->Get("Fisher_Discriminant");
        FTree->SetBranchAddress("Fs_start",&Fs);
        FTree->SetBranchAddress("Fb_start",&Fb);
        FTree->SetBranchAddress("nFs",&nFs);
        FTree->SetBranchAddress("nFb",&nFb);
        FTree->SetBranchAddress("Fish_coef_cc",&w_c);
        FTree->SetBranchAddress("Fish_coef_ANNout",&w_o);
	FTree->GetEntry(0);
	double w_i=0.;
	int ind=0;
	FTree->Branch("w_i",&w_i,"w_i/D");
	FTree->Branch("index",&ind,"index/I");

	//fF->Close();
	name.ReplaceAll("Fisher_file/","");

	double av;
        double cc;	
	double rho;
	int nSi;
	NNTree2->SetBranchAddress("Average",&av);
	NNTree2->SetBranchAddress("cc",&cc);
	NNTree2->SetBranchAddress("rho",&rho);
	NNTree2->SetBranchAddress("#TestSig",&nSi);
	NNTree2->GetEntry(0);
        
	int const nSig=nSi;
	cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
        cout<<" BG: "<<NNTree2->GetEntries()-nSi<<endl;
        int nBgi=NNTree2->GetEntries()-nSi; 
	double minw_i=10000000;
	double maxw_i=0.;
	
	for(int u=0;u<nSi;u++ ){
		if (u>=Fs&&u<(Fs+nFs))continue;
		ind=u;
		NNTree2->GetEntry(u);
		cout<<" u "<<u<<endl;
		w_i=w_c*cc+w_o*av;
		FTree->Fill();
		if(w_i<minw_i)minw_i=w_i;
		if(w_i>maxw_i)maxw_i=w_i;
	}
	int const nSigF=FTree->GetEntries();
	cout<<" nSigF "<<nSigF<<endl;
	for(int u=nSi;u<nSi+nBgi;u++ ){ 
		if (u>=Fb&&u<(Fb+nFb))continue;
		ind=u;
		NNTree2->GetEntry(u);
		cout<<" u "<<u<<endl;
		w_i=w_c*cc+w_o*av;
		FTree->Fill();
		if(w_i<minw_i)minw_i=w_i;
		if(w_i>maxw_i)maxw_i=w_i;
	}

	cout<<"F Entries "<<FTree->GetEntries()<<endl;
	
	double deltaw=0.;	
	deltaw=(maxw_i-minw_i)/50.;
	TH1D *hbg=new TH1D("Hbg", "projection[0]", 50, minw_i-deltaw, maxw_i+ deltaw);
	TH1D *hsig=new TH1D("Hsig", "projection[0]", 50, minw_i-deltaw, maxw_i+deltaw);
	hsig->SetDirectory(0);
	hbg->SetDirectory(0);


	
	for( int u=0; u<FTree->GetEntries()-1;u++){ // //-1 perchèxk hai tot c'è in npiù l'uno da Fisher_ex: i nuovi branches ripartono da zero
		FTree->GetEntry(u);
	cout<<" u "<<u<<" ind "<<ind<<" w_i "<<w_i<<" nSigF-1 "<<nSigF-1<<endl;
		if(u<nSigF-1) hsig->Fill(w_i);
		else hbg->Fill(w_i);
	}

//exit(0);
//Efficiency_vs_rho","Efficiency_vs_rho",0,0,1200,700
   TCanvas *c_histo= new TCanvas("Histogram_Fisher", "Histogram_Fisher", 0,0,900,600);
   hbg->SetLineColor(kBlue);
   hbg->SetFillStyle(3008);   hbg->SetFillColor(kBlue);
   hsig->SetLineColor(kRed);
   hsig->SetFillStyle(3003); hsig->SetFillColor(kRed);
   hsig->SetStats(0);
   hbg->SetStats(0);
//   isig->GetYaxis()->SetRangeUser(0,maxTrain/2);
//   ibg->GetYaxis()->SetRangeUser(0,maxTrain/2);
   hsig->Draw();
   hbg->Draw("same");
   TLegend *hlegend = new TLegend(.75, .80, .95, .95);
   hlegend->AddEntry(hbg, "Background");
   hlegend->AddEntry(hsig, "Signal");
   hlegend->Draw();
   char namech[2048];
   char namechroot[2048];
   TString namch(name);
   namch.ReplaceAll(".root",".png");
   sprintf(namech,"Fisher_png/Histog_%s",namch.Data());
   sprintf(namechroot,"Fisher_png/Histog_%s",name.Data());
   c_histo->SaveAs(namech); 
   c_histo->SaveAs(namechroot); 

/*	int const ncurve=nANN*nCC;
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
	double deltaANN=0.;
	deltacc=0.2/nCC;
	//deltaANN=0.6/(nANN-1);
	deltaANN=(maxw_i-minw_i)/nANN;

	cout<<NNTree2->GetEntries()<<endl;	
	for(int n=0;n<NNTree2->GetEntries();n++){
		cout<<n<<endl;
		NNTree2->GetEntry(n);
		w_i=w_c*cc+w_o*av;
		cout<<"rho "<<rho<<" cc "<<cc<<" w_i  "<<w_i<<endl;
        	for(int i=0; i<nCC;i++){
       			ccTh[i]=0.5+i*deltacc;
			if(cc<ccTh[i]) continue;
			for(int j=0; j<nANN;j++){
				if(j==0) NNTh[j]=-1000.;
				//else NNTh[j]=0.5+(j-1)*deltaANN;
				//else NNTh[j]=0.1+(j-1)*deltaANN;
				//else NNTh[j]=0.+(j)*deltaANN;
				else NNTh[j]=minw_i+(j)*deltaANN;
				//else NNTh[j]=0.+(j-1)*deltaANN/1000.;
				//else NNTh[j]=1.;
				if(w_i<NNTh[j]) continue;
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
//		int indexS[nSig];
//		double rhoS[nSig];
//		for(int i=0;i<nSig;i++) {
//			rhoS[i]=0.;
//			indexS[i]=0;
//		}
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
	fTEST->Close();
	fF->Close();
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
	sprintf(CnameS2,"logN_rho_Fish/logN_rho_S_dANN%1.2f_%s",deltaANN,CnameS.Data());
	sprintf(CnameB2,"logN_rho_Fish/logN_rho_B_dANN%1.2f_%s",deltaANN,CnameB.Data());
	sprintf(CnameS2root,"logN_rho_Fish/logN_rho_S_dANN%1.2f_%s",deltaANN,CnameSroot.Data());
	sprintf(CnameB2root,"logN_rho_Fish/logN_rho_B_dANN%1.2f_%s",deltaANN,CnameBroot.Data());
	cS->SaveAs(CnameS2);
	cB->SaveAs(CnameB2);
	cS->Print(CnameS2root);
	cB->Print(CnameB2root);
//	cout<<"fine"<<endl;
	*/
}

/*
void Annth(TString ifile){
        TString name(ifile);
        name.ReplaceAll("average_file/","Fisher_file/");
        TFile* fTEST =TFile::Open(ifile.Data());
        TTree* NNTree2=(TTree*)fTEST->Get("Parameters");

        double w_o;
        double w_c;
        int nFs;
        int nFb;
        int Fs;
        int Fb;

        TString namefF(name);
        namefF.ReplaceAll(".root","_Fisher.root");
        TFile* fF =TFile::Open(namefF.Data());
        TTree* FTree= (TTree*)fF->Get("Fisher_Discriminant");
        FTree->SetBranchAddress("Fs_start",&Fs);
        FTree->SetBranchAddress("Fb_start",&Fb);
        FTree->SetBranchAddress("nFs",&nFs);
        FTree->SetBranchAddress("nFb",&nFb);
        FTree->SetBranchAddress("Fish_coef_cc",&w_c);
        FTree->SetBranchAddress("Fish_coef_ANNout",&w_o);
        FTree->GetEntry(0);
        double w_i;
        int ind;
        FTree->SetBranchAddress("w_i",&w_i,"w_i/D");
        FTree->SetBranchAddress("index",&ind,"index/I");

	//TString name(ifile);
	//name.ReplaceAll("average_file/","");
	//TFile* fTEST =TFile::Open(ifile.Data());
	//TTree* NNTree2=(TTree*)fTEST->Get("Parameters");	           
	double av;
        double cc;	
	double rho;
	int nSi;
	NNTree2->SetBranchAddress("Average",&av);
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


        for(int n=0;n<FTree->GetEntries();n++){
		FTree->GetEntry(n)
                NNTree2->GetEntry(ind);
                cout<<"rho "<<rho<<" cc "<<cc<<" w_i  "<<w_i<<endl;
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
                                                ANNBg[i*nRHO+j][ni]=w_i;
                             //                   cout<<" soglia_cc "<<ccTh[i]<<" soglia_RHO "<<rhoTh[j]<<endl;
                                                }
                                        else {
                                                NSig[i*nRHO+j]= NSig[i*nRHO+j]+1;
                                                while(ANNSig[i*nRHO+j][ni]!=0)ni=ni+1;
                                                ANNSig[i*nRHO+j][ni]=w_i;
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
        sprintf(CnameS2,"ANNthres_Fish/N_ANN_S_%s",CnameS.Data());
        sprintf(CnameB2,"ANNthres_Fish/N_ANN_B_%s",CnameB.Data());
        sprintf(CnameS2root,"ANNthres_Fish/N_ANN_S_%s",CnameSroot.Data());
        sprintf(CnameB2root,"ANNthres_Fish/N_ANN_B_%s",CnameBroot.Data());
        cS->SaveAs(CnameS2);
        cB->SaveAs(CnameB2);
        cS->Print(CnameS2root);
        cB->Print(CnameB2root);
        //cout<<"fine"<<endl;

//CARTELLA Annth



}

*/
void PlotsAv_cc(TString ifile){
	TString name(ifile);
	name.ReplaceAll("outfile/","");
	name.ReplaceAll("average_file/","");
	TFile* fTEST =TFile::Open(ifile.Data());
	TTree* NNTree2=(TTree*)fTEST->Get("Parameters");	           
	//double outbis;
	double av;
	double out;
        double cc;	
	double rho;
	int nSi;
//	NNTree2->SetBranchAddress("ANNout",&out);
	NNTree2->SetBranchAddress("Average",&av);
	NNTree2->SetBranchAddress("cc",&cc);
	NNTree2->SetBranchAddress("rho",&rho);
	NNTree2->SetBranchAddress("#TestSig",&nSi);
	NNTree2->GetEntry(0);
	int const nSig=nSi;
	cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
	int const nBg=NNTree2->GetEntries()-nSig;
	cout<<"nBg: "<<nBg<<" Entries: "<<NNTree2->GetEntries()<<endl;

	TGraph * gS[3];
	TGraph * gB[3];
	gB[0]=new TGraph();
	gS[0]=new TGraph();
//	cout<<"nSig: "<<nSig<<endl;
	for (int n = 0; n <NNTree2->GetEntries(); n++){
		 NNTree2->GetEntry(n);
		 if(n<nSig) {
			gS[0]->SetPoint(n,cc,av);
	                cout<<"Sig_graph1: x="<<cc<<" y: "<<av<<endl;
		}
		else {
			gB[0]->SetPoint(n-nSig,cc,av);
	                cout<<"Bg_graph1: x="<<cc<<" y: "<<av<<endl;
		}
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
        gS[0]->SetMarkerColor(2);
        gB[0]->SetMarkerColor(4);
        gS[0]->SetMarkerStyle(6);
        gB[0]->SetMarkerStyle(7);

        TCanvas* c=new TCanvas("Plots","Plots",0,0,1200,700);
	c->Divide(1,2);
	c->cd(1);
        TMultiGraph* mg1=new TMultiGraph();
	mg1->SetTitle("Av_cc");
	if(gB[0]->GetN()!=0) mg1->Add(gB[0]);
	if(gS[0]->GetN()!=0) mg1->Add(gS[0]);
	mg1->Draw("ap");
	c->cd(2);
        TMultiGraph* mg2=new TMultiGraph();
	mg2->SetTitle("Av_cc");
	if(gS[0]->GetN()!=0) mg2->Add(gS[0]);
	if(gB[0]->GetN()!=0) mg2->Add(gB[0]);
	mg2->Draw("ap");

   cout<<" name "<<name<<endl;
   TString Cname(name);
   TString Cnameroot(name);
   char Cname2[1024];
   char Cname2root[1024];
   Cname.ReplaceAll(".root",".png");
   sprintf(Cname2,"average_png/out_Plots_%s",Cname.Data());
   sprintf(Cname2root,"average_png/out_Plots_%s",Cnameroot.Data());
   c->SaveAs(Cname2);
   c->Print(Cname2root);
/*
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

//CARTELLA CCi_RHO_ANN_Plots */
}
void PlotsAv_Mc(TString ifile){
        TString name(ifile);
        name.ReplaceAll("outfile/","");
        name.ReplaceAll("average_file/","");
        TFile* fTEST =TFile::Open(ifile.Data());
        TTree* NNTree2=(TTree*)fTEST->Get("Parameters");
        //double outbis;
        double av;
        double out;
        double cc;
        double Mc;
        double rho;
        int nSi;
//      NNTree2->SetBranchAddress("ANNout",&out);
        NNTree2->SetBranchAddress("Mchirp",&Mc);
        NNTree2->SetBranchAddress("Average",&av);
        NNTree2->SetBranchAddress("cc",&cc);
        NNTree2->SetBranchAddress("rho",&rho);
        NNTree2->SetBranchAddress("#TestSig",&nSi);
        NNTree2->GetEntry(0);
        int const nSig=nSi;
        cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
        int const nBg=NNTree2->GetEntries()-nSig;

        TGraph * gS[3];
        TGraph * gB[3];
        gB[0]=new TGraph();
        gS[0]=new TGraph();
//      cout<<"nSig: "<<nSig<<endl;
        for (int n = 0; n <NNTree2->GetEntries(); n++){
                 NNTree2->GetEntry(n);
	  if(n<nSig) {
                        gS[0]->SetPoint(n,Mc,av);
                        cout<<"Sig_graph1: x="<<cc<<" y: "<<av<<endl;
                }
                else {
                        gB[0]->SetPoint(n-nSig,Mc,av);
                        cout<<"Bg_graph1: x="<<cc<<" y: "<<av<<endl;
                }
        }


        /*for (int a=0;a<nBg;a++){
          if(aRHOB[a]>=RHOth){
                                if(c=10) cout<<" ANN "<<ANN1th[c]<<" Zcount: "<<Z[b*nth+c]<<" aNNB "<<aANNB[a]<<" a "<<a<<" b "<<b<<endl;
                        }
                  }
                }
          }
        }*/
        gS[0]->SetMarkerColor(2);
        gB[0]->SetMarkerColor(4);
        gS[0]->SetMarkerStyle(6);
        gB[0]->SetMarkerStyle(7);

        TCanvas* c=new TCanvas("Plots","Plots",0,0,1200,700);
        c->Divide(1,2);
        c->cd(1);
        TMultiGraph* mg1=new TMultiGraph();
        mg1->SetTitle("Av_Mc");
        if(gB[0]->GetN()!=0) mg1->Add(gB[0]);
        if(gS[0]->GetN()!=0) mg1->Add(gS[0]);
        mg1->Draw("ap");
        c->cd(2);
        TMultiGraph* mg2=new TMultiGraph();
        mg2->SetTitle("Av_Mc");
        if(gS[0]->GetN()!=0) mg2->Add(gS[0]);
        if(gB[0]->GetN()!=0) mg2->Add(gB[0]);
        mg2->Draw("ap");

   cout<<" name "<<name<<endl;
   TString Cname(name);
   TString Cnameroot(name);
   char Cname2[1024];
   char Cname2root[1024];
   Cname.ReplaceAll(".root",".png");
   sprintf(Cname2,"average_png/Mc_Plots_%s",Cname.Data());
   sprintf(Cname2root,"average_png/Mc_Plots_%s",Cnameroot.Data());
   c->SaveAs(Cname2);
   c->Print(Cname2root);
}
