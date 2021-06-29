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
#define Ncc3D 50
#define maxCC 1
#define minCC 0
#define nIFO 3
#define nRHO 3
#define nANN 10
#define deltaANN 0.1
#define deltacc 0.05
#define deltarho 1
#define nCC 4
#define NPIX 20
#define RHOth 5.5
#define nth 10
#define rhomin 6
#define ccmin 0.6
#define iMc 1
#define NMc_inj 50
#define minMc_inj 0
#define maxMc_inj 100
#define NMc_r 50
#define maxMc_r 40
#define minMc_r -5
#define NANN_av 50
#define minANN_av -0.2
#define maxANN_av 1.5

//definition parameters for SNR_plots
#define N_ANNth 50
#define N_RHOth 50
#define ANNth_max 1.5
#define ANNth_min -0.5
#define RHOth_max 10.0
#define RHOth_min 4.5

#define x_asym 0.00009999999999999999999999999999
#define nANNth 6
#define min_yasym -0.1
#define max_yasym  0.6
#define max_Mcth 2.5 
#define min_Mcth 1.5
#define npoints 50
#define max_k 10
#define min_k 0.01

#define cutMc1 0.5
#define cutMc2 1.5
#define cutANN1 0.15
#define cutANN2 0.25
#define maxrho 20.
#define minrho 0.


using namespace std;
void graph(TString ifile);
void Annth(TString ifile, int &FAth_n, int &FDth_n, double &thresFA);
void Mcth(TString ifile, int &McFAth_n, int &McFDth_n, double &McthresFA);
void PlotsAv_cc(TString ifile);
void PlotsAv_Mc(TString ifile);
void SNR_plots(TString ifile);
void CoG_plots(TString ifile);
void CutMcANN(TString ifile);
//void Test0(TString NN_FILE,TString TEST_FILE,TString ofile, int TS=0, int TB=0, int s=0, int b=0, int uf=0){

//NN_FILE-NN_FILE8: name of the files where the interested ANNs are saved (!!they must become from a directory whose name ends in "NN").
//TEST_FILE: name of the file which contains the events you are going to analyse (!!they must become from a directory whose name ends in "nnTREE").
//TS, TB: numbers of events to test for each class.
//s, b: index values of the first events considered for both the classes. If b < #SIG_events, b is newly set to b=b+#SIG_events
//uf !=0 means the Test file of events is the same used to train the network, in this case ni==0. It means that the events used for the training are analysed only by the other ANNs, in this way we test the events with at least NNn-1 networks. If uf==0 the Test and Train files are different, in this case ni!=0 and all the events can beanalysed by all the nNN networks.
//consider_all=1(0): the events not used for the training procedure are (aren't) considered for the test.
//ni: =0 excludes from the test the training set of events
void JoinCutANN_Mchirp_ROCcurves_rho2(TString NN_FILE,TString NN_FILE2, TString NN_FILE3, TString NN_FILE4, TString NN_FILE5, TString NN_FILE6, TString NN_FILE7, TString NN_FILE8, TString TEST_FILE, TString NOMEtot_S, int TS=0, int TB=0, int s=0, int b=0, int uf=1, int consider_all=0, int av_on_nNN=0){
int ni=0;
double FAdes=0.01;
if (uf==0&&av_on_nNN!=0) cout<<"all events were not used for the training procedure"<<endl;
if(uf==0) {ni=1; consider_all=1; av_on_nNN=0;}
else ni=0; //events used for the training not considered for the ANN trained on them

int TS0=0;
TS0=TS;
int TB0=0;
TB0=TB;
//if(uf!=0&&consider_all==0) ni=1;

//COUNT THE ANNs
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

if(nNN==1) {av_on_nNN=1; consider_all=0; ni=0;}
//saving the ANN names
   for (int u=0;u<nNN;u++){
   int p=0;
   while (NNi2[u][p]){
	if (NNi2[u][p]=='N'){
                        if((NNi2[u][p+1]=='N')&&(NNi2[u][p+2]=='\/')) {
					int hh=p+3;
					while (NNi2[u][hh]!='\0'){NNi[u][hh-p-3]=NNi2[u][hh];hh=hh+1;}
					break;
					}
	}
	p=p+1;
   }
}
   
//saving the Test-file name
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

//removing the extensions   
for (int u=0;u<nNN;u++){
   	NNi0[u].ReplaceAll(".root","");
   }
   Filei0.ReplaceAll(".root","");

//Definitions of some parameters 
   TFile *fnet[nNN];
   TMultiLayerPerceptron* mlp[nNN];
   TTree* infot[nNN];
   int TotBg[nNN];
   int FD0[nNN]; //FD=False Dismissal
   int FA0[nNN]; //FA=False Alarm

   for (int yy=0; yy<nNN; yy++){TotBg[yy]=0;FA0[yy]=0;FD0[yy]=0;}

   int NNs[nNN]; //first SIG-event considered for ANN trainings
   int NNb[nNN]; //first BKG-event considered for ANN trainings
   int NNnTS[nNN]; //number of SIG-events considered for ANN trainings
   int NNnTB[nNN]; //number of BKG-events considered for ANN trainings

int min_NNb=5000000;
int min_NNs=5000000;
int max_NNb=0;
int max_NNs=0;

//Extracting ANNs and information from the files
//if (uf!=0&&ni==0&&(TS!=0||TB!=0)){
   for (int u=0;u<nNN;u++){
    fnet[u]=TFile::Open(NNi2[u]);
    mlp[u]=(TMultiLayerPerceptron*)fnet[u]->Get("TMultiLayerPerceptron");
    
    cout<<"dopo TMLP"<<endl;
    if(mlp[u]==NULL) {cout << "Error getting mlp" << endl;exit(1);}
    infot[u]=(TTree*)fnet[u]->Get("info");
   
   infot[u]->SetBranchAddress("Rand_start_Sig",&NNs[u]);
   infot[u]->SetBranchAddress("Rand_start_Bg",&NNb[u]);
   infot[u]->SetBranchAddress("#trainSig",&NNnTS[u]);
   infot[u]->SetBranchAddress("#trainBg",&NNnTB[u]);
   infot[u]->GetEntry(0);
   cout<<"n "<<u<<" b "<<NNb[u]<<" min_NNb "<<min_NNb<<endl;

   if(NNs[u]<min_NNs) min_NNs=NNs[u];
   if(NNb[u]<min_NNb) min_NNb=NNb[u];
   if(NNs[u]>max_NNs) max_NNs=NNs[u];
   if(NNb[u]>max_NNb) max_NNb=NNb[u];
    cout<<"n "<<u<<" s  "<<NNs[u]<<" b "<<NNb[u]<<" s fin  "<<NNs[u]+NNnTS[u]<<" b fin "<<NNb[u]+NNnTB[u]<<" nTS "<<NNnTS[u]<<" nTB "<<NNnTB[u]<<endl;
   }

//Changing b value to its sum to the first event considered for the training
//if (uf!=0&&ni==0&&(b<min_NNb)&&(TS!=0||TB!=0)){ b+=min_NNb;cout<<" change in b value to have BKG events: "<<b<<endl;}
/*if (uf!=0&&ni==0&&(b<min_NNb)){ b+=min_NNb;cout<<" change in b value to have BKG events: "<<b<<endl;}
if (uf!=0&&ni==0&&(s<min_NNs)){ s+=min_NNs;cout<<" change in s value to have SIG events: "<<s<<endl;}*/
//if (uf!=0&&consider_all==0&&av_on_nNN==0&&(b<min_NNb)){ b+=min_NNb;cout<<" change in b value to have BKG events: "<<b<<endl;}
//if (uf!=0&&consider_all==0&&av_on_nNN==0&&(s<min_NNs)){ s+=min_NNs;cout<<" change in s value to have SIG events: "<<s<<endl;}

//if (uf!=0&&av_on_nNN!=0&&(b>min_NNb)&&(b<))
cout<<" min_NNb  "<<min_NNb<<" b "<<b<<endl;
cout<<"Def. tree e mlp"<<endl;
 
//opening the Test-file
   TFile* fTEST =TFile::Open(TEST_FILE.Data());
cout<<"Def. tree e mlp"<<endl;
   TTree* NNTree=(TTree*)fTEST->Get("nnTree");
cout<<"Def. tree e mlp"<<endl;
  int entries=NNTree->GetEntries();
//  int entries=0;

   cout<<"entries: "<<entries<<endl;
   cout<<" NNTree->GetEntries() "<<NNTree->GetEntries()<<endl;
//exit(0);
//Extract information necessary for the tree which defines the mlp
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
   cout<<"sig "<<sig_entries<<" bg "<<bg_entries<<endl;
   int const NDIM=ndim;
   int const nINP=ninp;
   cout<<"NDIM: "<<NDIM<<endl;
   cout<<"nINP: "<<nINP<<endl;
   cout<<"sig e: "<<sig_entries<<endl;
   cout<<"bg e: "<<bg_entries<<endl;

//defining the sample which has less events
   int minevents=0;
   if (sig_entries>bg_entries) minevents=bg_entries;
   else minevents=sig_entries;

   //if(b==0) b=sig_entries;
   if(TB==0 && TS==0){
	   TS=minevents;
	   TB=minevents;
   }
cout<<"TS "<<TS<<" TB "<<TB<<" minevents "<<minevents<<endl;
//defining compatible parameters
   if (b<sig_entries){
	int a =b;
   	b+=sig_entries;
	cout<<"Error: Bg index<sig_entries: new set of b parameter-> b="<<b<<" instead of b="<<a<<endl;
	//exit(0);
} 
 if(s>sig_entries){
  int a=s;
  s-=sig_entries;
 cout<<"Error: Sig index>sig_entries: new set of s parameter-> s="<<b<<" instead of s="<<a<<endl;
}
if((TS>sig_entries-s||TB>(bg_entries-(b-sig_entries)))&&(TS==TB)) {TS=minevents-s;TB=minevents-(b-sig_entries);
	if (TS<TB) TB=TS;
	else TS=TB;
	 cout<<"Error:S>sig_entries or TB>bg_entries, to maintain ugual number of analysed events TS=TB="<<TB<<endl;
	}
   if((TS>sig_entries-s||TB>bg_entries-(b-sig_entries))&&(TS!=TB)) {
	if(TS>sig_entries) TS=sig_entries-s;
	if (TB>bg_entries) TB=bg_entries-(b-sig_entries);
	cout<<"Error:S>sig_entries or TB>bg_entries, new TS and TB values are thus define-> TS="<<TS<<" TB="<<TB<<endl;
	}
   

   char NOMEtot[1024];
   sprintf(NOMEtot,"%s",NOMEtot_S.Data());
   cout<<"output name: "<<NOMEtot<<endl;

//extracting original files
   TString SIG_FILE;
   TString BG_FILE;
   char FILE_NAME[516];
   NNTree->SetBranchAddress("Files_name",&FILE_NAME);
   NNTree->GetEntry(0);
   SIG_FILE=FILE_NAME;
   NNTree->GetEntry(entries-1);
   BG_FILE=FILE_NAME;
   cout<<"fine ifdef RHO_CC"<<endl;
  
//extracting other information
   TChain sigTree("waveburst");//search the Tree "waveburst" in files
   sigTree.Add(SIG_FILE.Data());
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
  
   // define a branch for each input
   for(int i=0;i<nINP;i++) {
       sprintf(ilabel[i],"x%i",i+1);
       NNTree->SetBranchAddress(ilabel[i], &x[i]);
   }
  
char ofile[1024];
sprintf(ofile,"average_file/%s.root",NOMEtot);
TFile*f =new TFile(ofile,"RECREATE");
  TTree* NNTree2=new TTree("Parameters","Parameters");
  NNTree2->SetDirectory(f);

// defining useful information
   double average=0.;
   double out[nNN];
   for(int y=0; y<nNN; y++) out[y]=0.;
   double Mc=0.;
   double NNcc=0.;
   double NNrho=0.;
   double Mc_inj=0.;
   double std=0.;
   int ev_ind=0;
   int CoG=0;
  // int ev_ind;
   int index_ev=0;
   char TestFile[1024];
//adding new information to the original tree
  NNTree2->Branch("Center_of_Gravity",&CoG,"Center_of_Gravity/I");
  NNTree2->Branch("Average",&average,"Average/D");
  NNTree2->Branch("StandardDevaition",&std,"StandardDeviation/D");
  NNTree2->Branch("cc",&NNcc,"cc/D");
  NNTree2->Branch("Mchirp",&Mc,"Mchirp/D");
  NNTree2->Branch("Mchirp_injected",&Mc_inj,"Mchirp_injected/D");
  NNTree2->Branch("rho",&NNrho,"rho/D");
  NNTree2->Branch("index",&index_ev,"index/I");

//  NNTree2->Branch("event_index",&ev_ind,"event_index/I");
  NNTree2->Branch("Test_file",&TestFile,"Test_file/C");

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
    	NNTree2->Branch(NNnamel,&NNi[u],NNnamel2);
	}
  
  char Testf[1024];
  NNTree2->Branch("TestFile",&Testf,"TestFile/C");
  sprintf(Testf,"%s",TEST_FILE.Data());
  int nTestS=0;


//defining how many examples to use
  int scount=0;
  if(uf==0) nTestS=TS;
  else {
 	 //for(int n=s;n<s+TS;n++) {
 	 for(int n=s;n<sig_entries;n++) {
 		 int countNN=0;
 	   	 for(int u=0;u<nNN;u++){
 	    		 if (uf!=0&&n>NNs[u]&&n<(NNs[u]+NNnTS[u])) countNN=countNN+1;
			}

  //  if(consider_all==0&&av_on_nNN==0&&countNN==0)continue; //skip if no network is trained on the considered event
  //  if(consider_all==0&&av_on_nNN!=0&&countNN!=0) {sigskip=sigskip+1; continue; }

         	if(countNN==0&&av_on_nNN!=0) scount=scount+1;
		else{
         	//if(nNN==1&&countNN==1&&av_on_nNN!=0) scount=scount+1;
         	if(countNN==1&&consider_all==0&&av_on_nNN==0) scount=scount+1;
	 	if(countNN<2&&consider_all!=0) scount=scount+1; //if you are indifferently intereseted in average on nNN or nNN-1
         	if(countNN>1){cout<<"Error: training non independent"<<n<<endl; exit(0);} //exit if the ANNs are training on dependent set of events
		cout<<"n "<<n<<"; countNN "<<countNN<<"; scount "<< scount<<"; consider_all "<<consider_all<<"; av_on_nNN "<<av_on_nNN<<endl;
		if(scount==TS) break;}
	}
	nTestS=scount;
}
cout<<" nTestS "<< nTestS<<" TS "<<TS<<endl;
//exit(0);
int B_eff=0;
int FA=0;
int bg_05=0;
int cont_su=0;
int cont_su5=0;
int cont_su4=0;
int cont_su3=0;

int nTestB=0;
  int bcount=0;
  if(uf==0) nTestB=TB;
  else {
         //for(int n=b;n<b+TB;n++) {
         for(int n=b;n<sig_entries+bg_entries;n++) {
                 int countNN=0;
                 for(int u=0;u<nNN;u++){
                         if (uf!=0&&n>NNb[u]&&n<(NNb[u]+NNnTB[u])) countNN=countNN+1;
                        }
                if(countNN==0&&av_on_nNN!=0) bcount=bcount+1;
                //if(nNN==1&&countNN==1&&av_on_nNN!=0) bcount=bcount+1;
                if(countNN==1&&consider_all==0&&av_on_nNN==0) bcount=bcount+1;
                if(countNN<2&&consider_all!=0) bcount=bcount+1; //if you are indifferently intereseted in average on nNN or nNN-1
                if(countNN>1){cout<<"Error: training non independent"<<n<<endl; exit(0);} //exit if the ANNs are training on dependent set of events
                if(bcount==TB) break;
        }
        nTestB=bcount;
}

if(uf!=0) nTestS=nTestS-1;
if(uf!=0) nTestB=nTestB-1;
//TB=TB-1;
//TS=TS-1;
cout<<"TS "<<TS<<" TB "<<TB<<" nTestS "<<nTestS<<" nTestB "<<nTestB<<endl;
if(TS0==0&&TB0==0) {
			if(nTestB>nTestS) {nTestB=nTestS; TB=nTestS; TS=nTestS;}
			else {nTestS=nTestB; TB=nTestB; TS=nTestB;}
		}

if(nTestS<TS) TS= nTestS;
else nTestS=TS;

//  NNTree2->Branch("#TestSig",&TS,"#TestSig/I");
  NNTree2->Branch("#TestSig",&TS,"#TestSig/I");
  //NNTree2->Branch("#TestSig",&nTestS,"#TestSig/I");
  cout<<"nTestS: "<<nTestS<<" TS: "<<TS<<endl; 
  cout<<"dopo def tree"<<endl; 
cout<<"s test "<<s<< " s+TS "<<s+TS<<" b "<<b<<" b+ TB "<<b+TB<<" nTestS "<<nTestS<<" TB "<<TB<<" nTestB "<<nTestB<<endl;
double params[nINP];
   int sig_05=0;
   for (int i=0; i<nINP; i++)  params[i]=0.; 
 
  int S_eff=0; 
  int FD=0; 
int sigskip=0;

        TString txt_originaldata(NOMEtot_S);
        txt_originaldata.ReplaceAll(".root","_originalData.txt");
        //txt_originaldata.ReplaceAll(".root","_originalData.txt");
        char txt_originaldata_c[1024];
        sprintf(txt_originaldata_c,"txt_files/%s.txt",txt_originaldata.Data());
        ofstream od_txt(txt_originaldata_c);


  //for(int n=s;n<s+TS;n++) {
  for(int n=s;n<sig_entries;n++) {
double ti_0=0;
double fi_0=0;
     
     index_ev=n;
     average=0.;
     int indNN=-1;
     int countNN=0;
     for(int u=0;u<nNN;u++){
	     if (uf!=0&&n>=NNs[u]&&n<(NNs[u]+NNnTS[u])) { indNN=u;countNN=countNN+1;}	
	}
    if (countNN>1){cout<<"Error: training non independent"; exit(0);}
    if(nNN==1){cout<<"out==Averege:only 1NN is considered"<<endl;}
    if(consider_all==0&&av_on_nNN==0&&countNN==0)continue; //skip if no network is trained on the considered event
    if(consider_all==0&&av_on_nNN!=0&&countNN!=0) {sigskip=sigskip+1; continue; }
     NNTree->GetEntry(n);
     signal.GetEntry(n);
     char line_od[1024];

     // exracting parameters
     NNcc=(double)signal.netcc[0];
     NNrho=(double)signal.rho[1];
     Mc=(double)signal.chirp[iMc];
     Mc_inj=(double)signal.chirp[0];

     sprintf(TestFile,"%s",TEST_FILE.Data());

     sprintf(line_od,"%i \n",n);
     od_txt<<line_od;
     char line_od2[1024];
     // evaluating the event
     for (int i=0; i<nINP; i++){
        params[i]=x[i];	
     
     int ti_i=0;
     int fi_i=0;
     ti_i=i/NDIM;
     fi_i=i-ti_i*NDIM;
     ti_0+=(double)x[i]*ti_i; 
     fi_0+=(double)x[i]*fi_i;
        sprintf(line_od2,"%s %1.5f",line_od2, x[i]);
	}
        sprintf(line_od2,"%s \n",line_od2);

    if (ti_0<fi_0) CoG=1;
    else CoG=0;
 
//cout<<" ti_0 "<<ti_0<<" fi_0 "<<fi_0<<" CoG "<<CoG<<endl;
     od_txt<<line_od2;
     for (int u=0;u<nNN;u++) {
//        cout<<" u "<<u<<" nNN "<<nNN<<endl;
        double output=mlp[u]->Evaluate(0,params);
//        cout<<output<<endl;
        out[u]=output;
   	}
    
//    cout<<"Dopo param"<<endl;

//calculating the average and the standard deviation

  int c_nNN0=0;
    for (int u=0;u<nNN;u++){
        if((u!=indNN&&ni==0)||uf==0) {average=out[u]+average;std+=out[u]*out[u];c_nNN0=c_nNN0+1; if(out[u]<0.6) FD0[u]=FD0[u]+1;}
      }

      average=average/(c_nNN0); 
      if((c_nNN0)!=1)std=pow((std/(nNN-1-countNN)-average*average),0.5);
      else std=-1;

     char line_od3[1024];
      sprintf(line_od3,"ev_ind %i average %1.5f  Mc %1.5f  Mc_inj %1.5f n %i \n",index_ev, average,Mc,Mc_inj,n);
     od_txt<<line_od3;

    if (average<0.6) FD=FD+1;
     NNTree2->Fill();
     S_eff=S_eff+1;
    if(S_eff==1) s=n;
    if(S_eff==TS) break;
   cout<<" cc "<<NNcc<<" NNrho "<<NNrho<<"  Mc_inj "<< Mc_inj<<endl;
}

od_txt.close();
//exit(0);
cout<<txt_originaldata_c<<endl;
//exit(0);
//closing the cycle on SIG events
cout<<" S_eff "<<S_eff<<" s "<<s<<" b "<<b<<" sigskip "<<sigskip<<endl;
//exit(0);
//defining BKG parameters

int bgskip=0;
//opening the cycle on BKG events
   for(int n=b;n<sig_entries+bg_entries;n++) {
double ti_0=0;
double fi_0=0;
      index_ev=0;
      index_ev=n;
      average=0.;
      int indNN=-1;
      int countNN=0;
      for(int u=0;u<nNN;u++){
      cout<<" uf "<<uf<<" n "<<n<<" u "<<u<<" NN b[u] "<<NNb[u]<<" NNb[u]+NNnTB[u]"<<NNb[u]+NNnTB[u]<<endl;
             if (uf!=0&&n>=NNb[u]&&n<(NNb[u]+NNnTB[u])) { indNN=u;countNN=countNN+1; }
        }
    if (countNN>1){cout<<"Error: training non independent"; exit(0);}
    if(nNN==1){cout<<"out==Averege:only 1NN is considered"<<endl;}
    cout<<"consider_all "<<consider_all<<" av_on_nNN "<<av_on_nNN<<" countNN "<<countNN<<" indNN "<<indNN<<endl;
    if(consider_all==0&&av_on_nNN==0&&countNN==0)continue; //skip if no network is trained on the considered event
    if(consider_all==0&&av_on_nNN!=0&&countNN!=0) { cout<<" n skipped "<<n<<endl; bgskip=bgskip+1; continue;}
     NNTree->GetEntry(n);
     cout<<"BKG->n: "<<n<<"Bg index"<<(n-sig_entries)<<endl;
     background.GetEntry(n-sig_entries);
     //NNcc=(double)background.netcc[1];
     NNcc=(double)background.netcc[0];
     //NNrho=(double)background.rho[0];
     NNrho=(double)background.rho[1];
     Mc=(double)background.chirp[iMc];
     Mc_inj=(double)background.chirp[0];
     /*int hn=0;
     ev_ind=0;
     hn=n;
     ev_ind=hn;*/
     sprintf(TestFile,"%s",TEST_FILE.Data());
  
for (int i=0; i<nINP; i++){
        params[i]=x[i];
     int ti_i=0;
     int fi_i=0;
     ti_i=i/NDIM;
     fi_i=i-ti_i*NDIM;
     ti_0+=(double)x[i]*ti_i;
     fi_0+=(double)x[i]*fi_i;

        }
    if (ti_0<fi_0) CoG=1;
    else CoG=0;

	for(int u=0;u<nNN;u++) {
	     double output=mlp[u]->Evaluate(0,params);
    	     cout<<output<<endl;
	     out[u]=output;
	    }

     int c_nNN0=0;
    int nnsu=0;
    for (int u=0;u<nNN;u++){
        if((u!=indNN&&ni==0)||uf==0) {average=out[u]+average;std+=out[u]*out[u];c_nNN0=c_nNN0+1; TotBg[u]=TotBg[u]+1;if(out[u]>0.6) {FA0[u]=FA0[u]+1;nnsu=nnsu+1;}}
      }


      average=average/(c_nNN0);
      if((c_nNN0)!=1)std=pow((std/(nNN-1-countNN)-average*average),0.5);
      else std=-1;

    if(nnsu==6) cont_su=cont_su+1;
    if(nnsu>4) cont_su5=cont_su5+1;
    if(nnsu>3) cont_su4=cont_su4+1;
    if(nnsu>2) cont_su3=cont_su3+1;

//calculating the final average value


//saving false alarms
    if(average>0.6) FA=FA+1;
     
    NNTree2->Fill();
    B_eff=B_eff+1;
    if(B_eff==1) b=n;
    if(B_eff==TB) break;
     }
//cout<<"bkg filled"<<endl;

NNTree2->Write();
f->Close();
cout<<"B_eff "<<B_eff<<" TB "<<TB<<endl;
//cout<<"closed file"<<endl;
cout<<ofile<<endl;

int cont=0;

double freq_c[nNN];
double meanf=0.;
double dev=0.;

for (int yy=0; yy<nNN; yy++){
	freq_c[yy]=0.;
//	cout<<" FA0[yy] "<<FA0[yy]<<" FD0[yy] "<<FD0[yy]<<" FD "<<FD<<" TotBg[yy] "<<TotBg[yy]<<endl;
	if(TotBg[yy]!=0){freq_c[yy]=(double)FA0[yy]/TotBg[yy];
			//cout<<" freq_c[yy] "<<freq_c[yy]<<endl;
			}
	if(freq_c[yy]!=0) {meanf=freq_c[yy]+meanf;cont=cont+1;
			//cout<<" cont "<<cont<<endl;
			}
	dev=freq_c[yy]*freq_c[yy]+dev;
}

if(cont==0) cout<<"Error cont==0"<<endl;
if(cont!=0) meanf=meanf/cont;
dev=0.;
for (int yy=0; yy<nNN; yy++){
//cout<<" dev "<<dev<<endl;
//cout<<"freq_c[yy] "<<freq_c[yy]<<" meanf "<<meanf<<endl;
if(freq_c[yy]!=0)dev+=(freq_c[yy]-meanf)*(freq_c[yy]-meanf);
}
double McFAdes=FAdes;
int FAth_n=FAdes*B_eff;
int McFAth_n=FAdes*B_eff;
cout<<"FAth_n "<<FAth_n<<" FAdes*B_eff "<<FAdes*B_eff<<endl;
int FDth_n=-1;
int McFDth_n=-1;
double thresFA=-1;
double McthresFA=-1;
//Annth(ofile, FAth_n, FDth_n, thresFA);
//Mcth(ofile, McFAth_n, McFDth_n, McthresFA);
//graph(ofile);
//PlotsAv_cc(ofile);
//PlotsAv_Mc(ofile);
//SNR_plots(ofile);
//CoG_plots(ofile);
CutMcANN(ofile);
if(cont!=0) dev=pow(dev/(cont-1),0.5);
cout<<"meanf "<<meanf<<" dev "<<dev<<endl;
cout<<" FA average "<<FA<<endl;
cout<<" cont_su "<<cont_su<<" cont_su5 "<<cont_su5<<" cont_su4 "<<cont_su4<<" cont_su3 "<<cont_su3<<endl;

cout<<" S_eff "<<S_eff<<" B_eff "<<B_eff<<" s "<<s<<" b "<<b<<" bgskip "<<bgskip<<" sigskip "<<sigskip<<endl;
cout<<" FA0[1] "<<FA0[1]<<" TotBg[1] "<<TotBg[1]<<endl;
cout<<" FD0[1] "<<FD0[1]<<" TotBg[1] "<<TotBg[1]<<endl;
//double S_eff2=(double
double FDdes=(double)FDth_n/S_eff;
double McFDdes=(double)McFDth_n/S_eff;
FAdes=(double)FAth_n/B_eff;
McFAdes=(double)McFAth_n/B_eff;
cout<<"FAth_n "<<FAth_n<<" FAdes: "<<FAdes<<" ANN_av_FDdes "<<FDdes<<" ANN_av_FDth_n "<<FDth_n<<" ANN_av_thresFA "<<thresFA<<endl;
cout<<"McFAth_n "<<McFAth_n<<" McFAdes: "<<McFAdes<<" Mc_FDdes "<<McFDdes<<" McFDth_n "<<McFDth_n<<" McthresFA "<<McthresFA<<endl;
}

void graph(TString ifile){
	TString name(ifile);
	name.ReplaceAll("average_file/","");
	TFile* fTEST =TFile::Open(ifile.Data());
	TTree* NNTree2=(TTree*)fTEST->Get("Parameters");	           
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
	//double deltacc=0.;	
//	double deltaANN=0.;
//	deltacc=0.2/nCC;
	//deltaANN=0.6/(nANN-1);

	cout<<NNTree2->GetEntries()<<endl;	
	for(int n=0;n<NNTree2->GetEntries();n++){
		cout<<n<<endl;
		NNTree2->GetEntry(n);
		cout<<"rho "<<rho<<" cc "<<cc<<" av  "<<av<<endl;
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
				if(av<NNTh[j]) continue;
					int ni=0;
					if(n>=nSig) {
		
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
	//cS->cd(1)->SetLogy();
	cS->cd(1);
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
	//cS->cd(2)->SetLogy();
	cS->cd(2);
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
	//cS->cd(3)->SetLogy();
	cS->cd(3);
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
	//cS->cd(4)->SetLogy();
	cS->cd(4);
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
	sprintf(CnameS2,"logN_rho_av/logN_rho_S_dANN%1.2f_%s",deltaANN,CnameS.Data());
	sprintf(CnameB2,"logN_rho_av/logN_rho_B_dANN%1.2f_%s",deltaANN,CnameB.Data());
	sprintf(CnameS2root,"logN_rho_av/logN_rho_S_dANN%1.2f_%s",deltaANN,CnameSroot.Data());
	sprintf(CnameB2root,"logN_rho_av/logN_rho_B_dANN%1.2f_%s",deltaANN,CnameBroot.Data());
	cS->SaveAs(CnameS2);
	cB->SaveAs(CnameB2);
	cS->Print(CnameS2root);
	cB->Print(CnameB2root);
//	cout<<"fine"<<endl;
	
}


void Annth(TString ifile, int &FAth_n, int &FDth_n, double &thresFA){
	TString name(ifile);
	name.ReplaceAll("average_file/","");
	TFile* fTEST =TFile::Open(ifile.Data());
	TTree* NNTree2=(TTree*)fTEST->Get("Parameters");	           
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
       // double deltacc=0.;
       // double deltarho=0.;
       // deltacc=0.2/nCC;
//	deltarho=1./(nRHO);

//int thresFA=-1;
        for(int n=0;n<NNTree2->GetEntries();n++){
                NNTree2->GetEntry(n);
                cout<<"rho "<<rho<<" cc "<<cc<<" av "<<av<<endl;
                for(int i=0; i<nCC;i++){
                        ccTh[i]=0.5+i*deltacc;
			ccTh[0]=0.;
                        if(cc<ccTh[i]) continue;
                        for(int j=0; j<nRHO;j++){
                                rhoTh[j]=5+j*deltarho;
				rhoTh[0]=0.;
                                if(rho<rhoTh[j]) continue;
                                        int ni=0;
                                        if(n>=nSig) {
                                                NBg[i*nRHO+j]= NBg[i*nRHO+j]+1;
						if(i*nRHO+j==0) cout<<" NBg[i*nRHO+j] "<<NBg[i*nRHO+j]<<endl;
                                                while(ANNBg[i*nRHO+j][ni]!=0)ni=ni+1;
                                                ANNBg[i*nRHO+j][ni]=av;
				//		if (NBg[0]==FAth_n) thresFA=av;
                             //                   cout<<" soglia_cc "<<ccTh[i]<<" soglia_RHO "<<rhoTh[j]<<endl;
                                                }
                                        else {
                                                NSig[i*nRHO+j]= NSig[i*nRHO+j]+1;
						if(i*nRHO+j==0) cout<<" NSig[i*nRHO+j] "<<NSig[i*nRHO+j]<<endl;
                                                while(ANNSig[i*nRHO+j][ni]!=0)ni=ni+1;
                                                ANNSig[i*nRHO+j][ni]=av;
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


int thres_ind=0;
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
                        int ik=indexB[y][k-1];
                        int yy=0;
                        if (k>0){
                                int ij=indexB[y][k-1];
                                //if(ANNBg[y][ii]!=0&&ANNBg[y][ii]!=ANNBg[y][ij]) {
                                if(ANNBg[y][ii]!=0) {
                                                yy=NBg[y]-igB;
							if(y==0 && yy==FAth_n+1) {thresFA=ANNBg[y][ii];}
							if (y==0 && yy<FAth_n+1) if(ANNBg[y][ii]==ANNBg[y][ik]){FAth_n=FAth_n-1;}
						if(y==0) cout<<"yy "<<yy<<" y "<<y<<" FAth_n "<<FAth_n<<" thresFA "<<thresFA<<" ANNBg[y][ii] "<<ANNBg[y][ii]<<endl;
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
//int FDth_n=0;
//for(int i=0; i<nCC;i++){
//	for(int j=0; j<nRHO;j++){
FDth_n=0;
		for(int ni=1;ni<nSig;ni++){
			cout<<"FAth_n "<<FAth_n<<" thresFA "<<thresFA<<" ANNSig[0][indexS[0][ni]] "<<ANNSig[0][indexS[0][ni]]<<" FDth_n "<<FDth_n<<" NSig[0]-(ni-1 )"<<NSig[0]-(ni-1)<<" NSig[0] "<<NSig[0]<<endl;
//			if (ANNSig[0][indexS[0][ni]]==thresFA || (ANNSig[0][indexS[0][ni-1]]<thresFA && ANNSig[0][indexS[0][ni]]>thresFA)) FDth_n=NSig[0]-(ni-1);
			if (ANNSig[0][indexS[0][ni]]<=thresFA) FDth_n=FDth_n+1;
		}
//	}
//}
//exit(0);
  TCanvas* cS=new TCanvas("Efficiency_vs_ANN","Efficiency_vs_ANN",0,0,1200,700);
        cS->Divide(2,2);
        //cS->cd(1)->SetLogy();
        cS->cd(1);
        TMultiGraph* mg1=new TMultiGraph();
        mg1->SetTitle("cc: no_cut;ANN;#Events");
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
        mg1B->SetTitle("cc: no_cut;ANN;#Events");
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
        mg4B->SetTitle("cc=0.65;ANN;#Events");
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
        sprintf(CnameS2,"ANNthres_av/N_ANN_S_%s",CnameS.Data());
        sprintf(CnameB2,"ANNthres_av/N_ANN_B_%s",CnameB.Data());
        sprintf(CnameS2root,"ANNthres_av/N_ANN_S_%s",CnameSroot.Data());
        sprintf(CnameB2root,"ANNthres_av/N_ANN_B_%s",CnameBroot.Data());
        cS->SaveAs(CnameS2);
        cB->SaveAs(CnameB2);
        cS->Print(CnameS2root);
        cB->Print(CnameB2root);
        //cout<<"fine"<<endl;

//CARTELLA Annth



}

void Mcth(TString ifile, int &McFAth_n, int &McFDth_n, double &McthresFA){
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
        double Mc_inj;
        double rho;
        int nSi;
//      NNTree2->SetBranchAddress("ANNout",&out);
        NNTree2->SetBranchAddress("Mchirp",&Mc);
        NNTree2->SetBranchAddress("Average",&av);
        NNTree2->SetBranchAddress("cc",&cc);
        NNTree2->SetBranchAddress("rho",&rho);
        NNTree2->SetBranchAddress("Mchirp_injected",&Mc_inj);
        NNTree2->SetBranchAddress("#TestSig",&nSi);
        NNTree2->GetEntry(0);
        int const nSig=nSi;
  cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
        int const ncurve2=nRHO*nCC;

 double* McSig[ncurve2];
       for (int i=0;i<ncurve2;i++) McSig[i]=new double[nSig];
         //double rhoSig[20][10];
        int NSig[ncurve2];
//      int nBg=0;
        int const nBg=NNTree2->GetEntries()-nSig;
        for (int i=0;i<ncurve2;i++) {
                NSig[i]=0;
                for (int j=0;j<nSig;j++) McSig[i][j]=0.;
        }
        double* McBg[ncurve2];
        for (int i=0;i<ncurve2;i++) McBg[i]=new double[nBg];

        //double rhoBg[20][10];
        int NBg[ncurve2];
        for (int i=0;i<ncurve2;i++) {
                NBg[i]=0;
                for (int j=0;j<nBg;j++) McBg[i][j]=0.;
        }
        double ccTh[nCC];
        for (int i=0;i<nCC;i++) ccTh[i]=0.;
        double rhoTh[nRHO];
        for (int i=0;i<nRHO;i++) rhoTh[i]=0.;
//        double deltacc=0.;
//        double deltarho=0.;
//        deltacc=0.2/nCC;
//        deltarho=1./(nRHO);

        for(int n=0;n<NNTree2->GetEntries();n++){
                NNTree2->GetEntry(n);
                cout<<"rho "<<rho<<" cc "<<cc<<" av "<<av<<endl;
                for(int i=0; i<nCC;i++){
                        ccTh[i]=0.5+i*deltacc;
                        ccTh[0]=0.;
                        if(cc<ccTh[i]) continue;
                        for(int j=0; j<nRHO;j++){
                                rhoTh[j]=5+j*deltarho;
                                rhoTh[0]=0.;
                                if(rho<rhoTh[j]) continue;
                                        int ni=0;
                                        if(n>=nSig) {
                                                NBg[i*nRHO+j]= NBg[i*nRHO+j]+1;
                                                if(i*nRHO+j==0) cout<<" NBg[i*nRHO+j] "<<NBg[i*nRHO+j]<<endl;
                                                while(McBg[i*nRHO+j][ni]!=0)ni=ni+1;
                                                McBg[i*nRHO+j][ni]=Mc;
                                //              if (NBg[0]==FAth_n) thresFA=av;
                             //                   cout<<" soglia_cc "<<ccTh[i]<<" soglia_RHO "<<rhoTh[j]<<endl;
                                                }
                                        else {
                                                NSig[i*nRHO+j]= NSig[i*nRHO+j]+1;
                                                if(i*nRHO+j==0) cout<<" NSig[i*nRHO+j] "<<NSig[i*nRHO+j]<<endl;
                                                while(McSig[i*nRHO+j][ni]!=0)ni=ni+1;
                                                McSig[i*nRHO+j][ni]=Mc;
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
                TMath::Sort(nSig,McSig[y],indexS[y],false);
                for (int k=0;k<nSig;k++) {
                        int ii=indexS[y][k];
                        int yy=0;
                        if (k>0){
                                //cout<<"k "<<k<<endl;
                                int ij=indexS[y][k-1];
                                //if(ANNSig[y][ii]!=0&&ANNSig[y][ii]!=ANNSig[y][ij]) {
                                if(McSig[y][ii]!=0) {
                                                        yy=NSig[y]-igS;
                                                        //gS[y]->SetPoint(igS,ANNSig[y][ii],yy);
                                                        if(McSig[y][ii]!=McSig[y][ij]) gS[y]->SetPoint(igS_p++,McSig[y][ii],yy);
                   //                                     cout<<"igS"<<igS<<" x "<<ANNSig[y][ii]<<" y: "<<yy<<endl;
                                                        igS=igS+1;

                                                        }
                        }
			else {
                                if(McSig[y][ii]!=0){
                                        yy=NSig[y]-igS;
                                        gS[y]->SetPoint(0,McSig[y][ii],yy);
                                        igS=igS+1;
                 //                       cout<<" x "<<ANNSig[y][ii]<<" y: "<<yy<<endl;
                                }

                        }
                }
        }

int thres_ind=0;
        int* indexB[ncurve2];
        for (int i=0;i<ncurve2;i++) indexB[i]=new int[nBg];

        TGraph * gB[ncurve2];
        for (int y=0;y<ncurve2;y++) {
                int igB=0;
                int igB_p=0;
                gB[y]=new TGraph();
                TMath::Sort(nBg,McBg[y],indexB[y],false);
               // cout<<"dopo Sort "<<y<<endl;
                for (int k=0;k<nBg;k++) {
                        int ii=indexB[y][k];
                        int ik=indexB[y][k-1];
                        int yy=0;
                        if (k>0){
                                int ij=indexB[y][k-1];
                                //if(ANNBg[y][ii]!=0&&ANNBg[y][ii]!=ANNBg[y][ij]) {
                                if(McBg[y][ii]!=0) {
                                                yy=NBg[y]-igB;
                                                        if(y==0 && yy==McFAth_n+1) {McthresFA=McBg[y][ii];
                                                                }
							if (y==0 && yy<McFAth_n+1) { cout<<"dentro primo if"<<endl;
								 //if(McBg[y][ii]==McBg[y][ik]){McFAth_n=McFAth_n-1; cout<<"dentro secondo if"<<endl;}
								 if(McBg[y][ii]>McBg[y][ik]){McFAth_n=McFAth_n;cout<<"maggiore ii "<<" diff"<<(McBg[y][ii]-McBg[y][ik])<<endl;}
								 if(McBg[y][ii]<McBg[y][ik]){McFAth_n=McFAth_n;cout<<"minore ii"<<endl;}
								 if(McBg[y][ii]==McBg[y][ik]){McFAth_n=McFAth_n-1;cout<<"uguale"<<endl;}
								 //else { McFAth_n=McFAth_n-1;cout<<"dentro secondo if"<<endl;}
}
                                               if(y==0) cout<<"ii "<<ii<<" yy "<<yy<<" y "<<y<<" McFAth_n "<<McFAth_n<<"McthresFA "<<McthresFA<<" McBg[y][ii] "<<McBg[y][ii]<<" McBg[y][indexB[y][k-1] "<<McBg[y][indexB[y][k-1]]<<endl;
                                                //gB[y]->SetPoint(igB,ANNBg[y][ii],yy);
                                                if(McBg[y][ii]!=McBg[y][ij]) gB[y]->SetPoint(igB_p++,McBg[y][ii],yy);
                                                igB=igB+1;
             //                                   cout<<"igB"<<igB<<" x "<<ANNBg[y][ii]<<" y: "<<yy<<endl;
                                                }
                        }
                     else {
                                if(McBg[y][ii]!=0)     {
                                        yy=NBg[y]-igB;
                                        gB[y]->SetPoint(0,McBg[y][ii],yy);
                                        igB=igB+1;
           //                             cout<<"igB"<<igB<<" x "<<ANNBg[y][ii]<<" y: "<<yy<<endl;
                                }
                        }
                 }
        }

McFDth_n=0;
                for(int ni=1;ni<nSig;ni++){
                        cout<<"McFAth_n "<<McFAth_n<<" McthresFA "<<McthresFA<<" McSig[0][indexS[0][ni]] "<<McSig[0][indexS[0][ni]]<<" McFDth_n "<<McFDth_n<<" NSig[0]-(ni-1 )"<<NSig[0]-(ni-1)<<" NSig[0] "<<NSig[0]<<endl;
//                      if (ANNSig[0][indexS[0][ni]]==thresFA || (ANNSig[0][indexS[0][ni-1]]<thresFA && ANNSig[0][indexS[0][ni]]>thresFA)) FDth_n=NSig[0]-(ni-1);
                        if (McSig[0][indexS[0][ni]]<=McthresFA) McFDth_n=McFDth_n+1;
                }


  TCanvas* cS=new TCanvas("Efficiency_vs_Mchirp","Efficiency_vs_Mchirp",0,0,1200,700);
        cS->Divide(2,2);
        //cS->cd(1)->SetLogy();
        cS->cd(1);
        TMultiGraph* mg1=new TMultiGraph();
        mg1->SetTitle("cc: no_cut;Mchirp;#Events");
        for (int h=0;h<nRHO;h++){
                gS[h]->SetLineColor(4);
                if(gS[h]->GetN()!=0) mg1->Add(gS[h]);
                }
        mg1->Draw("al");
        cS->cd(2);
     //   cS->cd(2)->SetLogy();
        TMultiGraph* mg2=new TMultiGraph();
        mg2->SetTitle("cc=0.55;Mchirp;#Events");
        for (int h=0;h<nRHO;h++){
                gS[nRHO+h]->SetLineColor(4);
                if(gS[nRHO+h]->GetN()!=0) mg2->Add(gS[nRHO+h]);
                }
        mg2->Draw("al");

        //cS->cd(3)->SetLogy();
        cS->cd(3);
        TMultiGraph* mg3=new TMultiGraph();
        mg3->SetTitle("cc=0.6;Mchirp;#Events");
        for (int h=0;h<nRHO;h++){
                gS[2*nRHO+h]->SetLineColor(4);
                if(gS[2*nRHO+h]->GetN()!=0) mg3->Add(gS[2*nRHO+h]);
                }
      mg3->Draw("al");
        //cS->cd(4)->SetLogy();
        cS->cd(4);
        TMultiGraph* mg4=new TMultiGraph();
        mg4->SetTitle("cc=0.65;Mchirp;#Events");
        for (int h=0;h<nRHO;h++){
                gS[3*nRHO+h]->SetLineColor(4);
                if(gS[3*nRHO+h]->GetN()!=0) mg4->Add(gS[3*nRHO+h]);
                }
        mg4->Draw("al");

       TCanvas* cB=new TCanvas("Number_vs_Mchirp","Number_vs_Mchirp",0,0,1200,700);
        cB->Divide(2,2);
        cB->cd(1)->SetLogy();
        TMultiGraph* mg1B=new TMultiGraph();
        mg1B->SetTitle("cc: no_cut;Mchirp;#Events");
        for (int h=0;h<nRHO;h++){
                gB[h]->SetLineColor(4);
                if(gB[h]->GetN()!=0) mg1B->Add(gB[h]);
                }
        mg1B->Draw("al");
        cB->cd(2)->SetLogy();
        TMultiGraph* mg2B=new TMultiGraph();
        mg2B->SetTitle("cc=0.55;Mchirp;#Events");
        for (int h=0;h<nRHO;h++){
                gB[nRHO+h]->SetLineColor(4);
                if(gB[nRHO+h]->GetN()!=0) mg2B->Add(gB[nRHO+h]);
                }
        mg2B->Draw("al");
        cB->cd(3)->SetLogy();
        TMultiGraph* mg3B=new TMultiGraph();
        mg3B->SetTitle("cc=0.6;Mchirp;#Events");
        for (int h=0;h<nRHO;h++){
                gB[2*nRHO+h]->SetLineColor(4);
                if(gB[2*nRHO+h]->GetN()!=0) mg3B->Add(gB[2*nRHO+h]);
                }
        mg3B->Draw("al");
        cB->cd(4)->SetLogy();
        TMultiGraph* mg4B=new TMultiGraph();
        mg4B->SetTitle("cc=0.6;Mchirp;#Events");
        for (int h=0;h<nRHO;h++){
                gB[3*nRHO+h]->SetLineColor(4);
    if(gB[3*nRHO+h]->GetN()!=0) mg4B->Add(gB[3*nRHO+h]);
                }
        mg4B->Draw("al");


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
        sprintf(CnameS2,"Mcthres/N_Mc_S_%s",CnameS.Data());
        sprintf(CnameB2,"Mcthres/N_Mc_B_%s",CnameB.Data());
        sprintf(CnameS2root,"Mcthres/N_Mc_S_%s",CnameSroot.Data());
        sprintf(CnameB2root,"Mcthres/N_Mc_B_%s",CnameBroot.Data());
        cS->SaveAs(CnameS2);
        cB->SaveAs(CnameB2);
        cS->Print(CnameS2root);
        cB->Print(CnameB2root);


}

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
        //gS[0]->GetHistogram()->GetXaxis()->SetTitle("Mchirp");
        //gB[0]->GetHistogram()->GetYaxis()->SetTitle("Average on ANN ouputs");
        //gS[0]->GetHistogram()->GetXaxis()->SetTitle("Mchirp");
        //gB[0]->GetHistogram()->GetYaxis()->SetTitle("Average on ANN ouputs");
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
	mg1->GetHistogram()->GetXaxis()->SetTitle("cc");
	mg1->GetHistogram()->GetYaxis()->SetTitle("Average on ANN ouputs");
//	mg1->Draw("ap");
	c->cd(2);
        TMultiGraph* mg2=new TMultiGraph();
	mg2->SetTitle("Av_cc");
	if(gS[0]->GetN()!=0) mg2->Add(gS[0]);
	if(gB[0]->GetN()!=0) mg2->Add(gB[0]);
	mg2->Draw("ap");
	mg2->GetHistogram()->GetXaxis()->SetTitle("cc");
	mg2->GetHistogram()->GetYaxis()->SetTitle("Average on ANN ouputs");

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
        double Mc_inj;
        double rho;
        int nSi;
//      NNTree2->SetBranchAddress("ANNout",&out);
        NNTree2->SetBranchAddress("Mchirp",&Mc);
        NNTree2->SetBranchAddress("Average",&av);
        NNTree2->SetBranchAddress("cc",&cc);
        NNTree2->SetBranchAddress("rho",&rho);
        NNTree2->SetBranchAddress("Mchirp_injected",&Mc_inj);
        NNTree2->SetBranchAddress("#TestSig",&nSi);
        NNTree2->GetEntry(0);
        int const nSig=nSi;
        cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
	//exit(0);
        int const nBg=NNTree2->GetEntries()-nSig;

        TGraph * gS[3];
        TGraph * gB[3];
        gB[0]=new TGraph();
        gS[0]=new TGraph();


        TCanvas* ANN_Mc_c3D=new TCanvas("ANN_average_vs_Mc_reconstructed","ANN_average_vs_Mc_reconstructed",0,0,1200,700);
        ANN_Mc_c3D->Divide(1,2);
        TH2D* ANN_Mc_SIG3D=new TH2D("ANN_average_vs_Mc_recostructed_sig","ANN_average_vs_Mc_recostructed_sig",NMc_r,minMc_r,maxMc_r,NANN_av,minANN_av,maxANN_av);
        TH2D* ANN_Mc_BKG3D=new TH2D("ANN_average_vs_Mc_recostructed_bkg","ANN_average_vs_Mc_recostructed_bkg",NMc_r,minMc_r,maxMc_r,NANN_av,minANN_av,maxANN_av);
       // Mc_inj_g->SetPoint(0,rhoSig[y][ii],yy);
        ANN_Mc_BKG3D->SetDirectory(0);
        ANN_Mc_BKG3D->SetStats(0);
        ANN_Mc_SIG3D->SetDirectory(0);
        ANN_Mc_SIG3D->SetStats(0);

	  TCanvas* ANN_cc_c3D=new TCanvas("ANN_average_vs_cc","ANN_average_vs_cc",0,0,1200,700);
        ANN_cc_c3D->Divide(1,2);
        TH2D* ANN_cc_SIG3D=new TH2D("ANN_average_vs_cc_SIG","ANN_average_vs_cc_SIG",Ncc3D, minCC, maxCC,NANN_av,minANN_av,maxANN_av);
        TH2D* ANN_cc_BKG3D=new TH2D("ANN_average_vs_cc_BKG","ANN_average_vs_cc_BKG",Ncc3D, minCC, maxCC,NANN_av,minANN_av,maxANN_av);
        ANN_cc_SIG3D->SetDirectory(0);
        ANN_cc_SIG3D->SetStats(0);
        ANN_cc_BKG3D->SetDirectory(0);
        ANN_cc_BKG3D->SetStats(0);
        double deltacc3D=0;
        deltacc3D=(maxCC-minCC)/(Ncc3D);


        TCanvas* Mc_inj_c=new TCanvas("ANN_vs_Mc_injected","ANN_vs_Mc_injected",0,0,1200,700);
	Mc_inj_c->Divide(1,2);
	TH2D* Mc_inj_gMc=new TH2D("Mc_injected_Mc_recostructed","Mc_injected_recostructed",NMc_inj,minMc_inj,maxMc_inj,NMc_r,minMc_r,maxMc_r);
	TH2D* Mc_inj_g=new TH2D("Mc_injected_comparison","Mc_injected_comparison",NMc_inj,minMc_inj,maxMc_inj,NANN_av,minANN_av,maxANN_av);
       // Mc_inj_g->SetPoint(0,rhoSig[y][ii],yy);
	Mc_inj_gMc->SetDirectory(0);
	Mc_inj_gMc->SetStats(0);
	Mc_inj_g->SetDirectory(0);
	Mc_inj_g->SetStats(0);
	double deltaMc_inj=0;
	deltaMc_inj=(maxMc_inj-minMc_inj)/(NMc_inj);
	double deltaANN_av_M=0;
	deltaANN_av_M=(maxANN_av-minANN_av)/(NANN_av);
	double deltaMc_r=0;
	deltaMc_r=(maxMc_r-minMc_r)/(NMc_r);


        TString txt_name0(name);
	txt_name0.ReplaceAll(".root",".txt");
	char txt_name[1024];
	sprintf(txt_name,"txt_files/%s",txt_name0.Data());
	ofstream file_txt_out(txt_name);
  

//      cout<<"nSig: "<<nSig<<endl;
        int aa1=0;
        int bb1=0;
        for (int n = 0; n <NNTree2->GetEntries(); n++){
                 NNTree2->GetEntry(n);
			char line[2048];
	  if(n<nSig&&n!=nSig) {
			Mc_inj_gMc->Fill(Mc_inj,Mc);
			Mc_inj_g->Fill(Mc_inj,av);
                        gS[0]->SetPoint(n,Mc,av);
			cout<<" n "<<n<<" Mc "<<Mc<<" av "<<av<<" M_inj "<<Mc_inj<<endl;
                        cout<<"Sig_graph1: x="<<cc<<" y: "<<av<<endl;
			aa1=aa1+1;
			ANN_Mc_SIG3D->Fill(Mc,av);
                        ANN_cc_SIG3D->Fill(cc,av);
		        sprintf(line,"SIG:Average_%1.5f; Mc_rec_%1.5f; Mc_inj_%1.5f; event %i; count_%i; nSig_%i\n",av,Mc,Mc_inj,n,aa1,nSig);
			file_txt_out<<line;		
	  		cout<<" n "<<n<<" nSig "<<nSig<<endl;
                }
                else {
			bb1=bb1+1;
                        gB[0]->SetPoint(n-nSig,Mc,av);
                        ANN_Mc_BKG3D->Fill(Mc,av);
                        ANN_cc_BKG3D->Fill(cc,av);
	//		cout<<"Bg_graph1: x="<<cc<<" y: "<<av<<endl;
	 		sprintf(line,"BKG:Average %1.5f; Mc_rec %1.5f; Mc_inj_%1.5f; event %i; count_%i\n",av,Mc,Mc_inj,n,bb1);
			file_txt_out<<line;		
			
                }
        }
	  //exit(0);


        gS[0]->SetMarkerColor(2);
        gB[0]->SetMarkerColor(4);
        gS[0]->SetMarkerStyle(6);
        gB[0]->SetMarkerStyle(7);

        ANN_Mc_c3D->cd(1);
   ANN_Mc_SIG3D->GetXaxis()->SetTitle("SIG:Mchirp_reconstructed");
   ANN_Mc_SIG3D->GetYaxis()->SetTitle("SIG:ANN_average");
   ANN_Mc_SIG3D->GetZaxis()->SetTitle("count");
   ANN_Mc_SIG3D->Draw("colz");
        ANN_Mc_c3D->cd(2);
   ANN_Mc_BKG3D->GetXaxis()->SetTitle("BKG:Mchirp_reconstructed");
   ANN_Mc_BKG3D->GetYaxis()->SetTitle("BKG:ANN_average");
   ANN_Mc_BKG3D->GetZaxis()->SetTitle("count");
   ANN_Mc_BKG3D->Draw("colz");


   TString ANN_Mc_name3D(name);
   TString ANN_Mc_rootname3D(name);
   char ANN_Mc_cname23D[1024];
   char ANN_Mc_crootname23D[1024];
   ANN_Mc_name3D.ReplaceAll(".root",".png");
   sprintf(ANN_Mc_cname23D,"average_png/nSig%i_nBg%i_3DMc_ANNav_Plots_%s",nSig,nBg,ANN_Mc_name3D.Data());
   sprintf(ANN_Mc_crootname23D,"average_png/nSig%i_nBg%i_3DMc_ANNav_Plots_%s",nSig,nBg,ANN_Mc_rootname3D.Data());
   ANN_Mc_c3D->SaveAs(ANN_Mc_cname23D);
   ANN_Mc_c3D->Print(ANN_Mc_crootname23D);


        ANN_cc_c3D->cd(1);
   ANN_cc_SIG3D->GetXaxis()->SetTitle("SIG:cc");
   ANN_cc_SIG3D->GetYaxis()->SetTitle("SIG:ANN_average");
   ANN_cc_SIG3D->GetZaxis()->SetTitle("count");
   ANN_cc_SIG3D->Draw("colz");
        ANN_cc_c3D->cd(2);
   ANN_cc_BKG3D->GetXaxis()->SetTitle("BKG:cc");
   ANN_cc_BKG3D->GetYaxis()->SetTitle("BKG:ANN_average");
   ANN_cc_BKG3D->GetZaxis()->SetTitle("count");
   ANN_cc_BKG3D->Draw("colz");


   TString ANN_cc_cname3D(name);
   TString ANN_cc_crootname3D(name);
   char ANN_cc_cname23D[1024];
   char ANN_cc_crootname23D[1024];
   ANN_cc_cname3D.ReplaceAll(".root",".png");
   sprintf(ANN_cc_cname23D,"average_png/nSig%i_nBg%i_3DANNav_cc_Plots_%s",nSig,nBg,ANN_cc_cname3D.Data());
   sprintf(ANN_cc_crootname23D,"average_png/nSig%i_nBg%i_3DANNav_cc_Plots_%s",nSig,nBg,ANN_cc_crootname3D.Data());
   ANN_cc_c3D->SaveAs(ANN_cc_cname23D);
   ANN_cc_c3D->Print(ANN_cc_crootname23D);

	Mc_inj_c->cd(1);
   Mc_inj_g->GetXaxis()->SetTitle("Mchirp_injected");
   Mc_inj_g->GetYaxis()->SetTitle("ANN_average");
   Mc_inj_g->GetZaxis()->SetTitle("count");
   Mc_inj_g->Draw("colz");
	Mc_inj_c->cd(2);
   Mc_inj_gMc->GetXaxis()->SetTitle("Mchirp_injected");
   Mc_inj_gMc->GetYaxis()->SetTitle("Mchirp_reconstructed");
   Mc_inj_gMc->GetZaxis()->SetTitle("count");
   Mc_inj_gMc->Draw("colz");


   TString Mc_inj_cname(name);
   TString Mc_inj_crootname(name);
   char Mc_inj_cname2[1024];
   char Mc_inj_crootname2[1024];
   Mc_inj_cname.ReplaceAll(".root",".png");
   sprintf(Mc_inj_cname2,"average_png/nSig%i_nBg%i_Mc_inj_Plots_%s",nSig,nBg,Mc_inj_cname.Data());
   sprintf(Mc_inj_crootname2,"average_png/nSig%i_nBg%i_Mc_inj_Plots_%s",nSig,nBg,Mc_inj_crootname.Data());
  Mc_inj_c->SaveAs(Mc_inj_cname2);
   Mc_inj_c->Print(Mc_inj_crootname2);


	TCanvas* c=new TCanvas("Plots","Plots",0,0,1200,700);
        c->Divide(1,2);
        c->cd(1);
        TMultiGraph* mg1=new TMultiGraph();
        mg1->SetTitle("Av_Mc");
        if(gB[0]->GetN()!=0) mg1->Add(gB[0]);
        if(gS[0]->GetN()!=0) mg1->Add(gS[0]);
        mg1->Draw("ap");
	mg1->GetHistogram()->GetXaxis()->SetTitle("Mchirp");
	mg1->GetHistogram()->GetYaxis()->SetTitle("Average on ANN ouputs");
        c->cd(2);
        TMultiGraph* mg2=new TMultiGraph();
        mg2->SetTitle("Av_Mc");
        if(gS[0]->GetN()!=0) mg2->Add(gS[0]);
        if(gB[0]->GetN()!=0) mg2->Add(gB[0]);
        mg2->Draw("ap");
	mg2->GetHistogram()->GetXaxis()->SetTitle("Mchirp");
	mg2->GetHistogram()->GetYaxis()->SetTitle("Average on ANN ouputs");

   cout<<" name "<<name<<endl;
   TString Cname(name);
   TString Cnameroot(name);
   char Cname2[1024];
   char Cname2root[1024];
   Cname.ReplaceAll(".root",".png");
   sprintf(Cname2,"average_png/nSig%i_nBg%i_Mc_Plots_%s",nSig,nBg,Cname.Data());
   sprintf(Cname2root,"average_png/nSig%i_nBg%i_Mc_Plots_%s",nSig,nBg,Cnameroot.Data());
   cout<<"Cname2 "<<Cname2<<endl;
   c->SaveAs(Cname2);
   c->Print(Cname2root);
   cout<<" gS[0]->GetN() "<<gS[0]->GetN()<<endl;
        //cout<<"entries "<<mg1->GetHistogram()->GetEntries()<<endl;
   cout<<" nSig "<<nSig<<" nBg "<<nBg<<" aa1 "<<aa1<<endl;
}


void SNR_plots(TString ifile){

        TString name(ifile);
        name.ReplaceAll("outfile/","");
        name.ReplaceAll("average_file/","");
        TFile* fTEST =TFile::Open(ifile.Data());
        TTree* NNTree2=(TTree*)fTEST->Get("Parameters");
        double av;
        double out;
        double cc;
        double Mc;
        double Mc_inj;
        double rho;
        int nSi;
        NNTree2->SetBranchAddress("Mchirp",&Mc);
        NNTree2->SetBranchAddress("Average",&av);
        NNTree2->SetBranchAddress("cc",&cc);
        NNTree2->SetBranchAddress("rho",&rho);
        NNTree2->SetBranchAddress("Mchirp_injected",&Mc_inj);
        NNTree2->SetBranchAddress("#TestSig",&nSi);
        NNTree2->GetEntry(0);
        int const nSig=nSi;
        cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
        int const nBg=NNTree2->GetEntries()-nSig;
        int const nTot=NNTree2->GetEntries();
	int const n_points=N_ANNth*N_RHOth;
	int nBgSur[3][n_points];
	int nSigSur[3][n_points];
	double zax[3][n_points];
	
	for (int i=0; i<3; i++) for (int u=0; u<n_points; u++) {nBgSur[i][u]=0; nSigSur[i][u]=0;zax[i][u]=0.;}
	
	double pANNth[N_ANNth];
	double pRHOth[N_RHOth];

	double pCCth[3];

	pCCth[0]=0.;
	pCCth[1]=0.;
	pCCth[2]=0.;

	pCCth[0]=0.5;
	pCCth[1]=0.6;
	pCCth[2]=0.7;

	double deltaANNth=0.;
	double deltaRHOth=0.;

	deltaANNth=(double)(ANNth_max-ANNth_min)/N_ANNth;
	deltaRHOth=(double)(RHOth_max-RHOth_min)/N_RHOth;


        TCanvas* SNR_c=new TCanvas("SIG-BKG_reduction","SIG-BKG_reduction",0,0,1200,700);
        SNR_c->Divide(1,3);
	
        TH2D* SNR3Dcc0=new TH2D("cc=0.5","cc=0.5",N_ANNth,ANNth_min,ANNth_max,N_RHOth,RHOth_min,RHOth_max);
        TH2D* SNR3Dcc1=new TH2D("cc=0.6","cc=0.6",N_ANNth,ANNth_min,ANNth_max,N_RHOth,RHOth_min,RHOth_max);
        TH2D* SNR3Dcc2=new TH2D("cc=0.7","cc=0.7",N_ANNth,ANNth_min,ANNth_max,N_RHOth,RHOth_min,RHOth_max);
        SNR3Dcc0->SetDirectory(0);
        SNR3Dcc0->SetStats(0);
        SNR3Dcc1->SetDirectory(0);
        SNR3Dcc1->SetStats(0);
        SNR3Dcc2->SetDirectory(0);
        SNR3Dcc2->SetStats(0);

	for (int u=0; u<N_ANNth; u++) { pANNth[u]=0.; pANNth[u]=ANNth_min+u*deltaANNth; cout<<" pANNth[u] "<<pANNth[u]<<" u "<<u<<endl;}
	for (int u=0; u<N_RHOth; u++) { pRHOth[u]=0.;  pRHOth[u]=RHOth_min+u*deltaRHOth; cout<<" pRHOth[u] "<<pRHOth[u]<<" u "<<u<<endl;}

cout<<"fino def cose"<<endl;
/*
	for (int u=0; u<nTot; u++){
		NNTree2->GetEntry(u);
		cout<<" cc "<<cc<<" av "<<av<<" rho "<< rho<<endl;
		for (int y=2; y>-1; y--){
		//	cout<<" cc "<<cc<<" pCCth "<<pCCth[y]<<" y "<<y<<endl;
			if(cc>pCCth[y]) {		
				for (int i=N_ANNth-1; i>-1; i--){
		//		cout<<" av "<<av<<" pANNth[i] "<<pANNth[i]<<" i "<<i<<endl;
					if(av>pANNth[i]){
						for (int j=N_RHOth-1; j>-1; j++){
							//cout<<" rho "<<rho<<" pRHOth[j] "<<pRHOth[j]<<" j "<<j<<" nSig "<<nSig<<" u "<<u<<endl;
							 if(rho>pRHOth[j]){
								for(int yy=y; yy>-1;yy--) for (int ii=i; ii>-1; ii--) for (int jj=j; jj>-1;jj--){
						 			if(u<nSig) nSigSur[yy][jj+N_ANNth*ii]=nSigSur[yy][jj+N_ANNth*ii]+1;
						 			else nBgSur[yy][jj+N_ANNth*ii]=nBgSur[yy][jj+N_ANNth*ii]+1;
								//	cout<<" yy "<<yy<<" ii "<<ii<<" jj "<<jj<<" nBgSur[yy][jj+N_ANNth*ii] "<<nBgSur[yy][jj+N_ANNth*ii]<<" nSigSur[yy][jj+N_ANNth*ii] "<<nSigSur[yy][jj+N_ANNth*ii]<<endl;
								}
							 break;
							 }
							 else continue;	
						}
					break;
					}
					else continue;
				}
			break;
			}
			else continue;	
		}
	}	
//exit(0);
cout<<"dopo ciclo"<<endl;
*/
for (int u=0; u<nTot; u++){
                NNTree2->GetEntry(u);
		for (int y=0; y<3; y++){
			if(cc>pCCth[y]){
				for (int i=0; i<N_ANNth; i++){
					if(av>pANNth[i]) {
						for (int j=0; j<N_RHOth; j++){
							if(rho>pRHOth[j]){
								 if(u<nSig) nSigSur[y][j+N_ANNth*i]=nSigSur[y][j+N_ANNth*i]+1;
                                                                 else nBgSur[y][j+N_ANNth*i]=nBgSur[y][j+N_ANNth*i]+1;
							}
							else continue;
						}
					}
					else continue;
				}			
			}
			else continue;
		}	
	}


//exit(0);
int nBgSur0=0;
nBgSur0=nBgSur[0][0];
int nSigSur0=0;
nSigSur0=nSigSur[0][0];
if(nBgSur0==0){nBgSur0=nBg;}
if(nSigSur0==0){nSigSur0=nSig;}
int count_zax=0;
cout<<"dopo ciclo"<<" nSigSur0 "<<nSigSur0<<" nBgSur0 "<<nBgSur0<<endl;
for (int u=0; u<3; u++) for (int j=0; j<N_RHOth; j++) for (int i=0; i<N_ANNth; i++) {
	//cout<<" u "<<u<<" j "<<j<<" i "<<i<<" nSigSur[u][i+j*N_ANNth] "<<nSigSur[u][i+j*N_ANNth]<<endl;
	//if(nSigSur[u][i+j*N_ANNth]!=0){zax[u][i+j*N_ANNth]=(nBgSur[u][i+j*N_ANNth]/nBgSur0)/(nSigSur[u][i+j*N_ANNth]/nSigSur0);};
	if(nSigSur[u][i+j*N_ANNth]==0) { 
		if (nSigSur[u][i+j*N_ANNth]==0 && nBgSur[u][i+j*N_ANNth]==0) zax[u][i+j*N_ANNth]=-1; 
		else zax[u][i+j*N_ANNth]=-2;
	}
	else {
		double num=0.;
		double den=0.;
		num=(double)nBgSur[u][i+j*N_ANNth]/nBgSur0;
		den=(double)nSigSur0/nSigSur[u][i+j*N_ANNth];
		zax[u][i+j*N_ANNth]=(double)num/den; count_zax=count_zax+1;
		cout<<(double)(nBgSur[u][i+j*N_ANNth]/nBgSur0)*(nSigSur0/nSigSur[u][i+j*N_ANNth])<<endl;
		}
//		cout<<" nBgSur[yy][jj+N_ANNth*ii] "<<nBgSur[u][j+N_ANNth*i]<<" nSigSur[y][j+N_ANNth*i] "<<nSigSur[u][j+N_ANNth*i]<<" zax[u][i+j*N_ANNth] "<<zax[u][i+j*N_ANNth]<<endl;
	//else zax[u][i+j*N_ANNth]=(nBgSur[u][i+j*N_ANNth]/nBgSur0)/(nSigSur[u][i+j*N_ANNth]/nSigSur0);
	if(u==2){
	/*cout<<" zax[0][i+j*N_ANNth] "<<zax[0][i+j*N_ANNth]<<endl;
	cout<<" zax[1][i+j*N_ANNth] "<<zax[1][i+j*N_ANNth]<<endl;
	cout<<" zax[2][i+j*N_ANNth] "<<zax[2][i+j*N_ANNth]<<endl;*/
	SNR3Dcc0->Fill(pANNth[i]+ deltaANNth/2.,pRHOth[j]+ deltaRHOth/2.,zax[0][i+j*N_ANNth]);
	SNR3Dcc1->Fill(pANNth[i]+ deltaANNth/2.,pRHOth[j]+ deltaRHOth/2.,zax[1][i+j*N_ANNth]);
	SNR3Dcc2->Fill(pANNth[i]+ deltaANNth/2.,pRHOth[j]+ deltaRHOth/2.,zax[2][i+j*N_ANNth]);
	}
}
for(int yy=0; yy<3;yy++) {
	cout<<" pCCth "<<pCCth[yy]<<endl;
	for (int ii=N_ANNth-1; ii>-1; ii--) {
		cout<<"		pANNth[i] "<<pANNth[ii]<<endl;
		for (int jj=N_RHOth-1; jj>-1;jj--){
			cout<<" 			pRHOth[jj] "<<pRHOth[jj]<<endl;
			cout<<" 				 nBgSur[yy][ii+N_ANNth*jj] "<<nBgSur[yy][ii+N_ANNth*jj]<<endl;
			cout<<" 				 nSigSur[yy][ii+N_ANNth*jj] "<<nSigSur[yy][ii+N_ANNth*jj]<<endl;
			cout<<" 				 zax[0][i+j*N_ANNth]"<<zax[yy][ii+jj*N_ANNth]<<endl;
			cout<<" (double)(nBgSur[u][i+j*N_ANNth]/nBgSur0)"<<(double)nBgSur[yy][ii+jj*N_ANNth]/nBgSur0<<endl;
   //cout<<" yy "<<yy<<" ii "<<ii<<" jj "<<jj<<" nBgSur[yy][jj+N_ANNth*ii] "<<nBgSur[yy][jj+N_ANNth*ii]<<" nSigSur[yy][jj+N_ANNth*ii] "<<nSigSur[yy][jj+N_ANNth*ii]<<endl;
                                                                }
						}
			}
cout<<"coutn_zax "<<count_zax<<endl;
   SNR_c->cd(1);
   SNR3Dcc0->GetXaxis()->SetTitle("ANN_average_thresholds");
   SNR3Dcc0->GetYaxis()->SetTitle("RHO_thresholds");
   SNR3Dcc0->GetZaxis()->SetTitle("normBKG/normSIG");
   SNR3Dcc0->Draw("colz");

   SNR_c->cd(2);
   SNR3Dcc1->GetXaxis()->SetTitle("ANN_average_thresholds");
   SNR3Dcc1->GetYaxis()->SetTitle("RHO_thresholds");
   SNR3Dcc1->GetZaxis()->SetTitle("normBKG/normSIG");
   SNR3Dcc1->Draw("colz");

   SNR_c->cd(3);
   SNR3Dcc2->GetXaxis()->SetTitle("ANN_average_thresholds");
   SNR3Dcc2->GetYaxis()->SetTitle("RHO_thresholds");
   SNR3Dcc2->GetZaxis()->SetTitle("normBKG/normSIG");
   SNR3Dcc2->Draw("colz");

   TString SNR_n(name);
   TString SNR_nroot(name);
   char SNR_c_n[1024];
   char SNR_c_nroot[1024];
   SNR_n.ReplaceAll(".root",".png");
   sprintf(SNR_c_n,"SNR_plots/nSig%i_nBg%i_invSNR_Plots_%s",nSig,nBg,SNR_n.Data());
   sprintf(SNR_c_nroot,"SNR_plots/nSig%i_nBg%i_invSNR_Plots_%s",nSig,nBg,SNR_nroot.Data());
   SNR_c->SaveAs(SNR_c_n);
   SNR_c->Print(SNR_c_nroot);


}

void CoG_plots(TString ifile){

        TString name(ifile);
        name.ReplaceAll("outfile/","");
        name.ReplaceAll("average_file/","");
        TFile* fTEST =TFile::Open(ifile.Data());
        TTree* NNTree2=(TTree*)fTEST->Get("Parameters");
        double av;
        double out;
        double cc;
        double Mc;
        double Mc_inj;
        double rho;
        int nSi;
        int CoG;
        NNTree2->SetBranchAddress("Center_of_Gravity",&CoG);
        NNTree2->SetBranchAddress("Mchirp",&Mc);
        NNTree2->SetBranchAddress("Average",&av);
        NNTree2->SetBranchAddress("cc",&cc);
        NNTree2->SetBranchAddress("rho",&rho);
        NNTree2->SetBranchAddress("Mchirp_injected",&Mc_inj);
        NNTree2->SetBranchAddress("#TestSig",&nSi);
        NNTree2->GetEntry(0);
        int const nSig=nSi;
        cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
	int const nTot=NNTree2->GetEntries();
	int const nBg=nTot-nSig;
	
	TCanvas* CoG_c=new TCanvas("SIG-BKG_reduction","SIG-BKG_reduction",0,0,1200,700);

        TH1D* imp=new TH1D("Center_of_Gravity_VS_ANN","Center_of_Gravity_VS_ANN",N_ANNth,ANNth_min,ANNth_max);
        TH1D* CoG1Sig=new TH1D("Center_of_Gravity_VS_ANN","Center_of_Gravity_VS_ANN",N_ANNth,ANNth_min,ANNth_max);
        TH1D* CoG1Bg=new TH1D("Center_of_Gravity_VS_ANN","Center_of_Gravity_VS_ANN",N_ANNth,ANNth_min,ANNth_max);
        TH1D* CoG0Sig=new TH1D("Center_of_Gravity_VS_ANN","Center_of_Gravity_VS_ANN",N_ANNth,ANNth_min,ANNth_max);
        TH1D* CoG0Bg=new TH1D("Center_of_Gravity_VS_ANN","Center_of_Gravity_VS_ANN",N_ANNth,ANNth_min,ANNth_max);

       CoG1Sig->SetDirectory(0);
       CoG1Sig->SetStats(0);
       CoG0Sig->SetDirectory(0);
       CoG0Sig->SetStats(0);
       CoG1Bg->SetDirectory(0);
       CoG1Bg->SetStats(0);
       CoG0Bg->SetDirectory(0);
       CoG0Bg->SetStats(0);

	for (int n=0; n<nTot; n++){
	 NNTree2->GetEntry(n);
	 if (n<nSig) {
		if(CoG==0) {CoG0Sig->Fill(av);}
		else {CoG1Sig->Fill(av);}
	 }
	 else {
		 if(CoG==0) {CoG0Bg->Fill(av);}
                else {CoG1Bg->Fill(av);}
	 }
	 cout<<" CoG "<<CoG<<endl;
	}
	
   imp->SetDirectory(0);
   imp->SetStats(0);
   int max_n=0;
 //  if(nBg>nSig) max_n=nBg;
 //  else max_n=nSig;

   int max[4];
   for(int i=0; i<4; i++) max[i]=0;
 
   max[0]=CoG1Sig->GetMaximum();
   max[1]=CoG0Sig->GetMaximum();
   max[2]=CoG1Bg->GetMaximum();
   max[3]=CoG0Bg->GetMaximum();
   
   max_n=max[0];
   for(int i=0; i<4; i++) if(max_n<max[i]) max_n=max[i] ;

   imp->Fill(ANNth_min,max_n);
   imp->SetLineColor(10);
   imp->SetBins(N_ANNth,ANNth_min,ANNth_max);
   imp->GetXaxis()->SetTitle("ANN");
   imp->GetYaxis()->SetTitle("Count");
   imp->Draw();  

   CoG1Sig->SetFillStyle(3018);
   CoG0Sig->SetFillStyle(3018);
   CoG0Bg->SetFillStyle(3004);
   CoG1Bg->SetFillStyle(3004);
   CoG1Sig->SetFillColor(2);
   CoG0Sig->SetFillColor(6);
   CoG0Bg->SetFillColor(4);
   CoG1Bg->SetFillColor(9);
   CoG1Sig->SetLineColor(2);
   CoG0Sig->SetLineColor(6);
   CoG0Bg->SetLineColor(4);
   CoG1Bg->SetLineColor(9);
 //  CoG1Sig->SetBins(N_ANNth,ANNth_min,ANNth_max);
 //  CoG1Sig->GetXaxis()->SetTitle("ANN");
 //  CoG1Sig->GetYaxis()->SetTitle("Count");
   CoG1Sig->Draw("same");
   CoG0Sig->Draw("same");
   CoG0Bg->Draw("same");
   CoG1Bg->Draw("same");

   TString CoG_n(name);
   TString CoG_nroot(name);
   char CoG_c_n[1024];
   char CoG_c_nroot[1024];
   CoG_n.ReplaceAll(".root",".png");
   sprintf(CoG_c_n,"CoG_plots/nSig%i_nBg%i_invSNR_Plots_%s",nSig,nBg,CoG_n.Data());
   sprintf(CoG_c_nroot,"CoG_plots/nSig%i_nBg%i_invSNR_Plots_%s",nSig,nBg,CoG_nroot.Data());
   CoG_c->SaveAs(CoG_c_n);
   CoG_c->Print(CoG_c_nroot);

}

void CutMcANN(TString ifile){

        TString name(ifile);
        name.ReplaceAll("outfile/","");
        name.ReplaceAll("average_file/","");
        TFile* fTEST =TFile::Open(ifile.Data());
        TTree* NNTree2=(TTree*)fTEST->Get("Parameters");
        double av;
        double out;
        double cc;
        double Mc;
        double Mc_inj;
        double rho;
        int nSi;
        int CoG;
        NNTree2->SetBranchAddress("Center_of_Gravity",&CoG);
        NNTree2->SetBranchAddress("Mchirp",&Mc);
        NNTree2->SetBranchAddress("Average",&av);
        NNTree2->SetBranchAddress("cc",&cc);
        NNTree2->SetBranchAddress("rho",&rho);
        NNTree2->SetBranchAddress("Mchirp_injected",&Mc_inj);
        NNTree2->SetBranchAddress("#TestSig",&nSi);
        NNTree2->GetEntry(0);
        int const nSig=nSi;
        cout<<"nSig: "<<nSig<<" nSi: "<<nSi<<endl;
        int const nTot=NNTree2->GetEntries();
        int const nBg=nTot-nSig;


TGraph *Mc_g[nANNth];
TGraph *ANNMc_g[nANNth];
TGraph * rho_g[7];

double Mccuts[7];
double ANNcuts[7];
double TA_rho[7][npoints];
double FA_rho[7][npoints];
double FA_Mc[nANNth+1][npoints];
double FA_ANNMc[nANNth][npoints];
double TA_Mc[nANNth+1][npoints];
double TA_ANNMc[nANNth][npoints];
double th_Mc[npoints];
double th_k[npoints];
double th_rho[npoints];

for(int i=0; i<7; i++){
Mccuts[i]=0.0;
ANNcuts[i]=0.0;
}
ANNcuts[0]=-1.0;
ANNcuts[1]=cutANN2;
ANNcuts[2]=-1.0;
Mccuts[2]=cutMc2;
ANNcuts[3]=cutANN1;
Mccuts[3]=cutMc1;
Mccuts[4]=cutMc2;
ANNcuts[4]=cutANN2;
ANNcuts[5]=cutANN2;
Mccuts[5]=cutMc1;
Mccuts[6]=cutMc2;
ANNcuts[6]=cutANN1;

for(int i=0; i<7; i++){
	for (int j=0; j<npoints; j++){
		TA_rho[i][j]=0.;
		FA_rho[i][j]=0.;
	}
}
double delta_rho=0.;
delta_rho=(double)(maxrho-minrho)/npoints;

for(int i=0; i<npoints; i++){
th_rho[i]=0.;
th_rho[i]=minrho+i*delta_rho;
}
double yasym[nANNth];
//Mc_g[0]=new TGraph();
for (int i=0; i<nANNth; i++){
//	Mc_g[i+1]=new TGraph();
//	ANNMc_g[i]=new TGraph();
	for (int j=0; j<npoints; j++){
		FA_Mc[i+1][j]=0.;
		FA_ANNMc[i][j]=0.;
		TA_Mc[i+1][j]=0.;
		TA_ANNMc[i][j]=0.;
	}
}

//ANN thresholds
double delta_yasym=0.;
delta_yasym=(double)(max_yasym-min_yasym)/nANNth;

cout<<" delta_yasym "<<delta_yasym<<endl;
for (int i=0; i<nANNth; i++){
	yasym[i]=0.;
	yasym[i]=min_yasym+i*delta_yasym;
	cout<<" yasym[i] "<<yasym[i]<<" i "<<i<<endl;
}


//thresholds on the parameters necessary to built the TA-FA curves
double delta_Mc=0.;
//delta_Mc=(max_Mcth-min_Mcth)/nANNth;
delta_Mc=(double)(max_Mcth-min_Mcth)/npoints;

double delta_logk=0.;
double min_logk=0.;
double max_logk=0.;

min_logk=log10(min_k);
max_logk=log10(max_k);
delta_logk=(log10(max_k)-log10(min_k))/npoints;

for (int i=0; i<npoints; i++){
	FA_Mc[0][i]=0.;
	TA_Mc[0][i]=0.;
	th_Mc[i]=min_Mcth+i*delta_Mc;
	th_k[i]=pow(10,min_logk+i*delta_logk);
	cout<<" th_k[i] "<<th_k[i]<<" th_Mc[i] "<<th_Mc[i]<<endl;
	}

//FA and TA definitions in function of rho
int countS=0;
int countB=0;
for(int n=0; n<nTot; n++){
NNTree2->GetEntry(n);
        if(cc<0.5) continue;
        if(n<nSig) countS=countS+1;
        else countB=countB+1;
        for(int i=0; i<7; i++){
                if(Mc>Mccuts[i]&&av>ANNcuts[i]) {
                        for (int j=0; j<npoints;j++) {
                                if(rho>th_rho[j]&&n<nSig) TA_rho[i][j]=TA_rho[i][j]+1;
                                if(rho>th_rho[j]&&(n>nSig||n==nSig)) FA_rho[i][j]=FA_rho[i][j]+1;
                                }
                        }
                }
}


for(int i=0; i<7; i++){
        for (int j=0; j<npoints;j++) {
                if(i==6) cout<<"TA_rho[i][j] "<<TA_rho[i][j]<<" countS "<<countS<<" FA_rho[i][j] "<<FA_rho[i][j]<<" countB "<<countB<<endl;
                FA_rho[i][j]=(double)FA_rho[i][j]/countB;
                TA_rho[i][j]=(double)TA_rho[i][j]/countS;
		if(FA_rho[i][j]==0) FA_rho[i][j]=0.00001;
                if(i==6) cout<<"TA_rho[i][j] "<<TA_rho[i][j]<<" FA_rho[i][j] "<<FA_rho[i][j]<<endl;
                }
        }


// FA and TA definitions
for (int n=0; n<nTot; n++){
	NNTree2->GetEntry(n);
	if(n<nSig){
		for(int i=0; i<nANNth; i++){
			for(int j=0; j<npoints; j++){
				if(Mc>th_Mc[j]&&i==0)TA_Mc[0][j]=TA_Mc[0][j]+1;
				if(Mc>th_Mc[j]&&av>yasym[i]) TA_Mc[i+1][j]=TA_Mc[i+1][j]+1;
				if(av>((th_k[j]/(Mc-x_asym))+yasym[i]))TA_ANNMc[i][j]=TA_Mc[i+1][j]+1;
//			if(i==0&&j==0)cout<<" av "<<av<<" Mc "<<Mc<<" th_k[j] "<<th_k[j]<<" x_asym "<<x_asym<<" asym[i] "<<yasym[i]<<" ((th_k[j]/(Mc-x_asym))+yasym[i])) "<<(th_k[j]/(Mc-x_asym))+yasym[i]<<" TA_ANNMc[i][j] "<<TA_ANNMc[i][j]<<endl;
			}
		}
	}
		
	else {
		 for(int i=0; i<nANNth; i++){
			for(int j=0; j<npoints; j++){
				if(Mc>th_Mc[j]&&i==0)FA_Mc[0][j]=FA_Mc[0][j]+1;
				if(Mc>th_Mc[j]&&av>yasym[i])  FA_Mc[i+1][j]=FA_Mc[i+1][j]+1;
				if(av>(th_k[j]/(Mc-x_asym)+yasym[i]))FA_ANNMc[i][j]=FA_Mc[i+1][j]+1;
			}
		}
	}

}
//exit(0);

for(int j=0; j<npoints; j++){
	TA_Mc[0][j]=TA_Mc[0][j]/nSig;
        FA_Mc[0][j]=FA_Mc[0][j]/nBg;
	for(int i=0; i<nANNth; i++){
		TA_Mc[i+1][j]=(double)TA_Mc[i+1][j]/nSig;
		FA_Mc[i+1][j]=(double)FA_Mc[i+1][j]/nBg;
		TA_ANNMc[i][j]=(double)TA_ANNMc[i][j]/nSig;
		FA_ANNMc[i][j]=(double)FA_ANNMc[i][j]/nBg;
	}
}

Mc_g[0]=new TGraph(npoints,FA_Mc[0],TA_Mc[0]);
TGraph* g0=new TGraph();
g0->SetPoint(0.00001,1,0);
for (int i=0; i<7; i++){
        rho_g[i]=new TGraph(npoints,FA_rho[i],TA_rho[i]);
        }
 TCanvas* canrho=new TCanvas("Efficiency_vs_FalseAlarm","Efficiency_vs_Alarm",0,0,1200,700);
       canrho->SetLogx();
        TLegend *legrho = new TLegend(.75, .2, .99, .4);
        TMultiGraph* mg0=new TMultiGraph();
        if(g0->GetN()!=0) mg0->Add(g0);
        g0->SetMarkerColor(0);

        mg0->SetTitle("Comparison;FalseAlarms;Efficiency");
 for (int i=0;i<7;i++){
                rho_g[i]->SetLineColor(2+i);
                rho_g[i]->SetMarkerColor(2+i);
                rho_g[i]->SetMarkerStyle(i+20);
                //if((i+2)!=6&&(i+2)!=7&&(i+2)!=1) ANNMc_g[i]->SetMarkerSize(0.5);
                rho_g[i]->SetMarkerSize(0.5);
               // char def[1024];
                char defANNMc_rho[1024];
                sprintf(defANNMc_rho,"Mc_%1.2f & ANNcut_%1.2f",Mccuts[i],ANNcuts[i]);
                if(rho_g[i]->GetN()!=0) legrho->AddEntry(rho_g[i],defANNMc_rho,"pl");
                if(rho_g[i]->GetN()!=0) mg0->Add(rho_g[i]);
              //  if(ANNMc_g[0]->GetN()!=0) mg1->Add(ANNMc_g[0]);
                }
        mg0->Draw("apl");
        legrho->Draw();
 TString rhoCname(name);
 TString rhoCnamer(name);
 rhoCname.ReplaceAll(".root",".png");
 char rhoCname2[1024];
 char rhoCname2r[1024];
 sprintf(rhoCname2,"Cuts_caROC/TAvsFA_changingRHO%s",rhoCname.Data());
 sprintf(rhoCname2r,"Cuts_caROC/TAvsFA_changingRHO%s",rhoCnamer.Data());
 canrho->SaveAs(rhoCname2);
 canrho->Print(rhoCname2r);


for(int j=0; j<npoints; j++) cout<<" th_Mc[j] "<<th_Mc[j]<<" j "<<j<<" TA_Mc[0][j] "<<TA_Mc[0][j]<<" FA_Mc[0][j] "<<FA_Mc[0][j]<<endl;


for (int i=0; i<nANNth; i++){
	Mc_g[i+1]=new TGraph(npoints,FA_Mc[i+1],TA_Mc[i+1]);
	ANNMc_g[i]=new TGraph(npoints,FA_ANNMc[i],TA_ANNMc[i]);
}

 TCanvas* can=new TCanvas("Efficiency_vs_FalseAlarm","Efficiency_vs_Alarm",0,0,1200,700);
       can->SetLogx();
	TLegend *leg = new TLegend(.75, .4, .99, .65);
	//TLegend *leg =new TLengend(.75, .80, .95, .95);
        TMultiGraph* mg1=new TMultiGraph();
        mg1->SetTitle("Comparison;FalseAlarms;Efficiency");
        Mc_g[0]->SetLineColor(3);
        Mc_g[0]->SetMarkerStyle(1);
        for (int i=0;i<nANNth;i++){
                Mc_g[i+1]->SetLineColor(4);
                //if((i+2)!=6&&(i+1)!=7&&(i+2)!=1)Mc_g[i+1]->SetMarkerSize(0.5);
                Mc_g[i+1]->SetMarkerSize(0.5);
		Mc_g[i+1]->SetMarkerStyle(i+20);
                Mc_g[i+1]->SetMarkerColor(4);
		ANNMc_g[i]->SetLineColor(2);
		ANNMc_g[i]->SetMarkerColor(2);
		ANNMc_g[i]->SetMarkerStyle(i+20);
		//if((i+2)!=6&&(i+2)!=7&&(i+2)!=1) ANNMc_g[i]->SetMarkerSize(0.5);
		ANNMc_g[i]->SetMarkerSize(0.5);
		char defMc[1024];
		char defANNMc[1024];
		sprintf(defMc,"Mc & ANNcut_%1.2f",yasym[i]);
		sprintf(defANNMc,"Hyperbola_cut,y_asymtoto_%1.2f","yasym[i]");
                if(Mc_g[i+1]->GetN()!=0) leg->AddEntry(Mc_g[i+1],defMc,"pl");
                if(Mc_g[i+1]->GetN()!=0) mg1->Add(Mc_g[i+1]);
                if(ANNMc_g[i]->GetN()!=0) leg->AddEntry(ANNMc_g[i],defANNMc,"pl");
                if(ANNMc_g[i]->GetN()!=0) mg1->Add(ANNMc_g[i]);
              //  if(ANNMc_g[0]->GetN()!=0) mg1->Add(ANNMc_g[0]);
                }
	if(Mc_g[0]->GetN()!=0) mg1->Add(Mc_g[0]);
	if(Mc_g[0]->GetN()!=0) leg->AddEntry(Mc_g[0],"Mc, no ANNcut","pl");
        mg1->Draw("apl");
	leg->Draw();
// ANNMc_g[0]= new TGraph(npoints, )
// cout<<" TA_ANNMc[0][0] "<<TA_ANNMc[0][0]<<" FA_ANNMc[0][0] "<<FA_ANNMc[0][0]<<endl;
 TString Cname(name);
 TString Cnamer(name);
 Cname.ReplaceAll(".root",".png");
 char Cname2[1024];
 char Cname2r[1024];
 sprintf(Cname2,"Cuts_caROC/TAvsFA_%s",Cname.Data());
 sprintf(Cname2r,"Cuts_caROC/TAvsFA_%s",Cnamer.Data());
 can->SaveAs(Cname2);
 can->Print(Cname2r);
 
}


