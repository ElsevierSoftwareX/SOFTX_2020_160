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


//
// Neural Network Analysis
// Author : Serena Vinciguerra
// Note : run with ACLiC -> root -l macro/cwb_mlp.C+
//

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


//#define SIG_FILE "/home/serena.vinciguerra/test0/S6D_R11_SIM_EOBNRv2_advL1H1V1_2G_run1_2/merge/NN_wave_S6D_R11_SIM_EOBNRv2_advL1H1V1_2G_run1_2.M1.root"
//#define BG_FILE  "/home/serena.vinciguerra/test0/S6D_R11_BKG_L1H1V1_2G_MP_run1/merge/NN_wave_S6D_R11_BKG_L1H1V1_2G_MP_run1.M1.root"
#define NRHO 50
#define NCC 50
#define nIFO  3		// number of detectors
//#define differentW
#define RHO_CC
#define CREATE_NETWORK
#define trainTEST
TString select7Train_non_es(TString ORIGINAL_FILE,std::vector<int>* choice,TString IDname);

//NOTA_BENE:OGNI VOLTA CHE VIENE RICHIAMATA LA MACRO LA NN PRECEDENTE VIENE SOVRASCRITTA
//NOTA_BENE: OGNI VOLTA CHE FACCIO UN GRAFICO ROC SI SOVRASCRIVE SE CAMBIO SOLO COSE DIVERSE DA NTRAINBG E NTRAINS E Y
//NOTA_BENE: DOVE SI USA QST MACRO DEVONO ESISTERE LA CARRTELLA png(se def Graph), tree, NN 

void cwbANN_purity_SoB(int nTrain,int nTrainBg, TString option,TString TREE_FILE,int q,int p,double rm=0., double cm=0., int lm=5, Int_t nepoch=700) {
   
   if (!gROOT->GetClass("TMultiLayerPerceptron")) {
      gSystem->Load("libMLP");
   }
   
  //gSystem->Load("select7_C.so");  
//  gROOT->LoadMacro("select7Train_non_es.C");
// cout<<soglia<<endl;
   double error;
   //int y=0;
//   int q=0;//sig
//   int p=0;//bg
   double maxTrain=0;
   // cout<<"y:"<<y<<endl;
   if(nTrain>nTrainBg) maxTrain=nTrain;
   else maxTrain=nTrainBg;

   TFile* ifile=TFile::Open(TREE_FILE.Data());  
   TTree* NNTree=(TTree*)ifile->Get("nnTree");
   //NNTree->SetDirectory(0);
   //ifile->Close();

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

//contolli
    if(sig_entries==bg_entries) cout<<"Error: training with a single type of event"<<endl;

    if(bg_entries<nTrainBg){cout<<"Error: Bg: num train>num events"<<endl;
			exit(0);}
    if(sig_entries<nTrain) {cout<<"Error: Sig: num train>num events"<<endl;
			exit(0);}

    if(!nTrainBg) nTrainBg=bg_entries;
    if(!nTrain) nTrain=sig_entries;
//fine controlli

//definizione start
  // q=gRandom->Integer(sig_entries-nTrain);
  // p=gRandom->Integer(bg_entries-nTrainBg)+sig_entries;
    //q=1;
    //p=1;

   cout<<"p"<<p<<endl;
   cout<<"q"<<q<<endl;
//fine defnizione start
   TString TreeF(TREE_FILE);
   TreeF.ReplaceAll("/home/serena.vinciguerra/","");
   TreeF.ReplaceAll("/","_");
   TreeF.ReplaceAll(".root","");
TString structure(option);
   structure.ReplaceAll(":","_");
   char nomeNN[512];
   sprintf(nomeNN,"NN/mlpNetwork_colored_nTS%i_nTB%i_structure%s_lm%i_epochs%i_s%i_b%i_ROC_TREEfile_%s.root",nTrain,nTrainBg,structure.Data(),lm,nepoch,q,p,TreeF.Data());
cout<<nomeNN<<endl;
   char nameTrain[512];
   sprintf(nameTrain,"PropertyGraphs/mlpNetwork_colored_nTS%i_nTB%i_s%i_b%i_TREEfile_%s.root",nTrain,nTrainBg,q,p,TreeF.Data());
#ifdef RHO_CC
TString SIG_FILE;
TString BG_FILE;
char FILE_NAME[516];
NNTree->SetBranchAddress("Files_name",&FILE_NAME);
NNTree->GetEntry(0);
SIG_FILE=FILE_NAME;
NNTree->GetEntry(entries-1);
BG_FILE=FILE_NAME;
cout<<"fine ifdef RHO_CC"<<endl;
#endif

TString FileSelectName(nameTrain);
TString nomefile;
FileSelectName.ReplaceAll("PropertyGraphs/mlpNetwork","selectTREE/TreeTrain");
TFile *f0=TFile::Open(FileSelectName.Data(),"read");

if (!f0) {

//RICOSTRUZIONE TREE per Train
std::vector<int>* choice=new vector<int>;
cout<<"vector definition"<<endl;
//choice->push_back(q);
for (int z=0;z<nTrain;z++) choice->push_back(q+z);
//(*choice)[nTrain]=p;
for (int z=0;z<nTrainBg;z++) choice->push_back(z+p);
//   gSystem->Load("select7_C.so");  

nomefile=select7Train_non_es(TREE_FILE.Data(),choice,nameTrain);
cout<<"new file created"<<endl;
}
else {
nomefile=FileSelectName.Data();
cout<<"open file-Train pre-existing: "<<nomefile.Data()<<endl;
f0->Close();}

//char nomefile[516];
//sprintf(nomefile,"selectTREE/TreeTrain_size%i_s%i_b%i.root",choice->size(),q,p);
TFile* fsel=TFile::Open(nomefile.Data());
TTree* tTrain1=(TTree*)fsel->Get("nnTree");
cout<<"def Tree per Train"<<endl;
 // add leaf
   Int_t type;
   Float_t x[nINP];
   char ilabel[nINP][16];   
  // char llabel[nINP][16];
   for(int i=0;i<nINP;i++) {
     sprintf(ilabel[i],"x%d",i+1); 
     // cout<<"ilabel "<<ilabel[i]<<endl;
    // sprintf(llabel[i],"%s/F",ilabel[i]);
    // tTrain->Branch(ilabel[i], &x[i],llabel[i]);
     tTrain1->SetBranchAddress(ilabel[i], &x[i]);
   }
   tTrain1->SetBranchAddress("type",&type);
   #ifdef differentW
   Float_t w[nINP];
   char wlabel[nINP][16];   
  // char llabel2[nINP][16];

   for(int i=0;i<nINP;i++) {
     sprintf(wlabel[i],"w%d",i+1); 
  //   sprintf(llabel2[i],"%s/F",wlabel[i]);
     tTrain1->SetBranchAddress(wlabel[i], &w[i]);
   }
   #else
   float w;
   tTrain1->SetBranchAddress("w",&w);
   #endif
 ifile->Close();

#ifdef RHO_CC
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

   TH1F *ibg = new TH1F("bgh", "rho[0]", 50, 0, 50);
   TH1F *isig = new TH1F("sigh", "rho[0]", 50, 0, 18);
   TH1F *ibg2 = new TH1F("bgh2", "netcc[1]", 50, 0, 1.2);
   TH1F *isig2 = new TH1F("sigh2", "netcc[1]", 50, 0, 1.1);
//   TH1F *ibg3 = new TH1F("bgh3", "test secondo giro", 50, -1, 2);
//   TH1F *isig3 = new TH1F("sigh3", "test secondo giro", 50, -1, 2);

   ibg->SetDirectory(0);
   isig->SetDirectory(0);
   ibg2->SetDirectory(0);
   isig2->SetDirectory(0);
//   ibg3->SetDirectory(0);
//   isig3->SetDirectory(0);
   Int_t i;
   double minRHO=0.;
   double minCC=0.;
   double maxRHO=0.;
   double maxCC=0.;

TTree* t_NNout=new TTree("rho_cc_NNout","rho_cc_NNout");
double rho_NN=0.;
double cc_NN=0.;
double NNout=0.;
int NNtype;
t_NNout->Branch("rho",&rho_NN,"rho/D");
t_NNout->Branch("cc",&cc_NN,"cc/D");
t_NNout->Branch("type",&NNtype,"type/I");
t_NNout->Branch("NNout",&NNout,"NNout/D");

for(int h=q;h<q+nTrain;h++){
        int resto0=0;
        if(q%2==0) resto0=(h)%2;
        else resto0=(h+1)%2;
	sigTree.GetEntry(h);
        if (resto0==0&&h<q+2) minCC=(double)signal.netcc[1];
        if (resto0==0&&h<q+2) minRHO=(double)signal.rho[0];
        if (signal.rho[0]<minRHO&&resto0==0)minRHO=(double)signal.rho[0];
        if (signal.netcc[1]<minCC&&resto0==0)minCC=(double)signal.netcc[1];
        if (signal.rho[0]>maxRHO&&resto0==0)maxRHO=(double)signal.rho[0];
        if (signal.netcc[1]>maxCC&&resto0==0)maxCC=(double)signal.netcc[1];
	isig->Fill(signal.rho[0]);
	isig2->Fill(signal.netcc[1]);
	}

for(int h=(p-sig_entries);h<(p-sig_entries+nTrainBg);h++){
        int resto0=0;
        if((p-sig_entries+nTrainBg)%2==0) resto0=(h)%2;
        else resto0=(h+1)%2;
	bgTree.GetEntry(h);
        if (resto0==0&&h<q+2) minCC=(double)signal.netcc[1];
        if (resto0==0&&h<q+2) minRHO=(double)signal.rho[0];
        if (background.rho[0]<minRHO&&resto0==0)minRHO=(double)background.rho[0];
        if (background.netcc[1]<minCC&&resto0==0)minCC=(double)background.netcc[1];
        if (background.rho[0]>maxRHO&&resto0==0)maxRHO=(double)background.rho[0];
        if (background.netcc[1]>maxCC&&resto0==0)maxCC=(double)background.netcc[1];
	ibg->Fill(background.rho[0]);
	ibg2->Fill(background.netcc[1]);
	}
cout<<"max RHO"<<maxRHO<<" maxCC "<<maxCC<<" minRHO "<<minRHO<<" minCC "<<minCC<<endl;
#endif
if (cm!=0) maxCC=cm;
if (rm!=0) maxRHO=rm;


TString layout;
for (int i=0;i<nINP-1;i++) {
char add[512];
sprintf(add,"@%s,",ilabel[i]);
layout+=add;

}
char add1[512];
sprintf(add1,"@%s:%s:type",ilabel[nINP-1],option.Data());
layout+=add1;
cout<<layout.Data()<<endl;

//  char layout[5096]="";
//   for (int i=0;i<nINP-1;i++) sprintf(layout,"%s@%s,",layout,ilabel[i]);
//   sprintf(layout,"%s@%s:%s:type",layout,ilabel[nINP-1],option.Data());
  // cout << "MultiLayer Perceptron layout -> " << layout << endl;

#ifdef CREATE_NETWORK
   TCanvas *errors=new TCanvas("errors","errors",0,0,600,600);
   TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron(layout,"w",tTrain1,"Entry$%2","(Entry$+1)%2");
	if (lm==1)mlp->SetLearningMethod(TMultiLayerPerceptron::kStochastic);
	if(lm==2)mlp->SetLearningMethod(TMultiLayerPerceptron::kBatch);
	if(lm==3)mlp->SetLearningMethod(TMultiLayerPerceptron::kSteepestDescent);
	if(lm==4)mlp->SetLearningMethod(TMultiLayerPerceptron::kRibierePolak);
	if(lm==5)mlp->SetLearningMethod(TMultiLayerPerceptron::kFletcherReeves);
	if(lm==6)mlp->SetLearningMethod(TMultiLayerPerceptron::kBFGS);
	if(lm>6||lm<1){
	cout<<"Invalid Learning Method:1<=lm<=6"<<endl;
	exit(0);
	}

   mlp->Train(nepoch, "text,current,graph,update=10");
  // char nomePNG[128];
   TString PNGname(nomeNN);
   TString PNGnameroot(nomeNN);
   PNGname.ReplaceAll("NN/","ErrorGraphs/");
   PNGnameroot.ReplaceAll("NN/","ErrorGraphs/");
   PNGname.ReplaceAll(".root",".png");
//   sprintf(nomePNG,"ErrorGraphs/mlpNetwork_nTS%i_nTB%i_lm%i_epochs%i_s%i_b%i.png",nTrain,nTrainBg,q,p);
   errors->SaveAs(PNGname.Data());
   errors->Print(PNGnameroot.Data());
  // char nomeNN[128];
  // sprintf(nomeNN,"NN/mlpNetwork_nTS%i_nTB%i__lm%i_epochs%i_s%i_b%i.root",nTrain,nTrainBg,q,p);
  // TFile fnet(nomeNN,"RECREATE");
   TFile fnet(nomeNN,"RECREATE");

   mlp->Write();
//   fnet.Close(); 

#else
cout<<"ci sono"<<endl;
   TFile fnet(nomeNN);
cout<<"non ci sono più"<<endl;
//   TFile fnet = TFile::Open(nomeNN);
   TMultiLayerPerceptron *mlp = (TMultiLayerPerceptron*)fnet.Get("TMultiLayerPerceptron");
   if(mlp==NULL) {cout << "Error getting mlp" << endl;exit(1);}
   //fnet->Close();
#endif

   // Use TMLPAnalyzer to see what it looks for
   TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis");
#ifdef RHO_CC
   mlpa_canvas->Divide(2,2);
  
   //DISEGNA L'ISTROGRMMA DI RHO
   mlpa_canvas->cd(1);
   ibg->SetLineColor(kBlue);
   ibg->SetFillStyle(3008);   ibg->SetFillColor(kBlue);
   isig->SetLineColor(kRed);
   isig->SetFillStyle(3003); isig->SetFillColor(kRed);
   ibg->SetStats(0);
   isig->SetStats(0);
//   isig->GetYaxis()->SetRangeUser(0,maxTrain/2);
//   ibg->GetYaxis()->SetRangeUser(0,maxTrain/2);
   ibg->Draw();
   isig->Draw("same");
   TLegend *ilegend = new TLegend(.75, .80, .95, .95);
   ilegend->AddEntry(ibg, "Background");
   ilegend->AddEntry(isig, "Signal");
   ilegend->Draw();
   cout << "Integral ibg : " << ibg->GetEntries() << endl;
   cout << "Integral isig : " << isig->GetEntries() << endl;


   //DISEGNA L'ISTOGRMMA CC
   mlpa_canvas->cd(2);
  // ibg2->GetYaxis()->SetRangeUser(0,maxTrain/2);
  // isig2->GetYaxis()->SetRangeUser(0,maxTrain/2);
   ibg2->SetLineColor(kBlue);
   ibg2->SetFillStyle(3008);   ibg2->SetFillColor(kBlue);
   isig2->SetLineColor(kRed);
   isig2->SetFillStyle(3003); isig2->SetFillColor(kRed);
   ibg2->SetStats(0);
   isig2->SetStats(0);
   ibg2->Draw();
   isig2->Draw("same");
   TLegend *ilegend2 = new TLegend(.75, .80, .95, .95);
   ilegend2->AddEntry(ibg2, "Background");
   ilegend2->AddEntry(isig2, "Signal");
   ilegend2->Draw();
   cout << "Integral ibg : " << ibg2->GetEntries() << endl;
   cout << "Integral isig : " << isig2->GetEntries() << endl;
#endif

   // Build and train the NN ptsumf is used as a weight since we are primarly 

   TMLPAnalyzer ana(mlp);
   // Initialisation
   ana.GatherInformations();
   // output to the console
   ana.CheckNetwork();
#ifdef CREATE_NETWORK
#ifdef RHO_CC
   mlpa_canvas->cd(4);
   // shows how each variable influences the network
   ana.DrawDInputs();
  // mlpa_canvas->cd(2);
   // shows the network structure
  // mlp->Draw();
   //ana.DrawDInputs();
   mlpa_canvas->cd(3);
   // draws the resulting network
  gPad->SetLogy();
   ana.DrawNetwork(0,"type==1","type==0"); //disegna solo il test sample
#else
mlpa_canvas->Divide(1,2);
mlpa_canvas->cd(1);
   ana.DrawDInputs();
mlpa_canvas->cd(2);
  gPad->SetLogy();
   ana.DrawNetwork(0,"type==1","type==0"); //disegna solo il test sample
#endif

//DEFINISCO IL TREE DOVE SALVO TUTTE LE INFO SUL SET TEST

TTree* t= new TTree("info","info");
t->Branch("amplitude_mode",&y,"amplitude_mode/I");
// int nDim=NDIM;
t->Branch("Matrix_dim",&ndim,"Matrix_dim/I");
// int nInp=nINP;
t->Branch("#inputs",&ninp,"#inputs/I");
 char NNname[516];
 sprintf(NNname,"%s",option.Data());
t->Branch("NNname",&NNname,"NNname/C");
t->Branch("#trainBg",&nTrainBg,"#trainBg/I");
t->Branch("#trainSig",&nTrain,"#trainSig/I");
//double threshold;
//t->Branch("threshold",&soglia,"threshold/D");
t->Branch("Rand_start_Bg",&p,"Rand_start_Bg/I");
t->Branch("Rand_start_Sig",&q,"Rand_start_Sig/I");
t->Branch("epochs",&nepoch,"epocs/I");
t->Branch("learning_method",&lm,"learning_method/I");
t->Fill();

t->Write();
//t->SetDirectory(fnet);

double errTot=mlp->GetError(TMultiLayerPerceptron::kTest);                        
	

mlpa_canvas->cd(0);
mlpa_canvas->Write();
TString CanvNN(nomeNN);
TString CanvNNroot(nomeNN);
CanvNN.ReplaceAll("NN/","mlpa_canv/");
CanvNNroot.ReplaceAll("NN/","mlpa_canv/");
CanvNN.ReplaceAll(".root",".png");
mlpa_canvas->SaveAs(CanvNN.Data());
mlpa_canvas->Print(CanvNNroot.Data());
#endif
fnet.Close();	
//fnet.Close();	

//Testin the Train sample
   int const npunti=40;
   double params[nINP];
   float TA[npunti+1];
   float FA[npunti+1];
   float thres0[npunti+1];
   float TA1[npunti+1];
   float FA1[npunti+1];
//   float thres01[npunti+1];
   
   for (int i=0;i<=npunti;i++){
     TA[i]=0;
     TA1[i]=0;
     FA[i]=0;
     FA1[i]=0;
  //   thres01[i]=0;
     thres0[i]=0;
     }
   float passo=0;
   passo=(log(1.1)-log(0.5))/(npunti/2);
   cout<<"passo "<<passo<<endl;
   thres0[npunti]=0.5;
   cout<<"soglia0 "<<thres0[0]<<endl;
   for(int i=npunti;i>npunti/2;i--) {
        thres0[i-1]=thres0[i]*exp(passo);}
   for(int i=0;i<(npunti/2);i++) {
        thres0[npunti/2+i+1]=1.6 -thres0[npunti/2+i+1];
        thres0[npunti/2-i-1]=1.-thres0[npunti/2+i+1];
//        cout<<"i: "<<i+11<<" "<<thres0[11+i]<<endl;
//        cout<<"i: "<<11-i<<" "<<thres0[11-i]<<endl;
//        thres0[i+1]=thres0[i]+passo;
        }
thres0[npunti/2]=0.5;
   for(int i=0;i<(npunti);i++) {
         cout<<thres0[i]<<endl;
        }
   
/*for(int i=0;i<(npunti/2);i++) {
        thres0[npunti/2+i+1]=thres0[npunti/2+i]*exp(passo);}
   for(int i=0;i<(npunti/2);i++) {
        thres0[npunti/2+i+1]=1.5 -thres0[npunti/2+i+1];
        thres0[npunti/2-i-1]=1.-thres0[npunti/2+i+1];
//        thres0[i+1]=thres0[i]+passo;
        }
*/
//TH2D* PUREZZA_S=new TH2D("Purezza","Purezza",NCC,0.,1.2,NRHO,0.,50);
TH2D* PUREZZA_S=new TH2D("Purezza","Purezza",NCC,minCC,maxCC,NRHO,minRHO,maxRHO);
TH2D* h2BgO_R=new TH2D("rho_cc_out","rho_cc_out",50,0.,50,50,-0.5,1.5);
TH2D* h2SigO_R=new TH2D("rho_cc_out","rho_cc_out",50,0.,50,50,-0.5,1.5);
TH2D* h2BgR_C=new TH2D("rho_cc_out","rho_cc_out",50,0.,1.2,50,0.,50);
TH2D* h2SigR_C=new TH2D("rho_cc_out","rho_cc_out",50,0.,1.2,50,0.,50);
TH2D* h2Bg=new TH2D("rho_cc_out","rho_cc_out",50,0.,1.2,50,-0.5,1.5);
TH2D* h2Sig=new TH2D("rho_cc_out","rho_cc_out",50,0.,1.2,50,-0.5,1.5);
h2Bg->SetDirectory(0);
h2Sig->SetDirectory(0);
PUREZZA_S->SetDirectory(0);
    double tt1=0;
    double yy1=0;
    double ty1=0;
    double tmedio1=0;
    double ymedio1=0;
    double tt=0;
    double yy=0;
    double ty=0;
    double tmedio=0;
    double ymedio=0;
   double rho1[NRHO];
   double cc1[NCC];
int const CCRHO=NRHO*NCC;
   double pur[CCRHO];
   double tot[CCRHO];
   for (int ii=0;ii<NRHO; ii++) {
		rho1[ii]=0.;
		for (int ji=0;ji<NCC; ji++) {
					if (ii==0) cc1[ji]=0.;
					pur[NRHO*ii+ji]=0.;
					tot[NRHO*ii+ji]=0.;
					}
		}
    double deltacc=0.;
    //deltacc=1.2/NRHO;
    deltacc=(maxCC-minCC)/NCC;
    double deltarho=0.;
    //deltarho=50/NRHO;
    deltarho=(maxRHO-minRHO)/NRHO;
    //for (int ii=0;ii<NRHO;ii++) { rho1[ii]=50.0-ii*deltarho;  
    for (int ii=0;ii<NRHO;ii++) { 
	rho1[ii]=maxRHO-(ii+1)*deltarho; 
/*	cout<<" rho_i "<<rho1[ii];
	rho1[ii]=100*rho1[ii];
        int rho0=0;
        rho0=(int)(rho1[ii]+0.5); 
	rho1[ii]=0.;
        rho1[ii]=(double)rho0/100;
	cout<<" rho_f "<<rho1[ii]<<endl;
*/	}
    for (int ji=0;ji<NCC;ji++) { 
	cc1[ji]=(ji)*deltacc+minCC; 
	/*cout<<" cc_i "<<cc1[ji]<<endl;
        cc1[ji]=100*cc1[ji];
        int cc0=0;
        cc0=(int)(cc1[ji]+0.5);
        cc1[ji]=0.;
        cc1[ji]=(double)cc0/100;
        cout<<" cc_f "<<cc1[ji]<<endl;
*/	}
    for (int i = 0; i <(double)tTrain1->GetEntries(); i++) //o va +1 o +0,5?
         {
           tTrain1->GetEntry(i);
           for (int j=0;j<nINP;j++)
             {
                params[j]=0;
                params[j]=(double)x[j];
             }
           int j=npunti;
           double output=mlp->Evaluate(0,params);
           cout<<"output "<<output;
		//<<endl;
    int resto=i%2;
           if (resto==0){
		NNout=output;
                if(type==1){
			NNtype=1;
			rho_NN=(double)signal.rho[0];
			cc_NN=(double)signal.netcc[1];
			}
		if(type==0){
			NNtype=0;
			rho_NN=(double)background.rho[0];
			cc_NN=(double)background.netcc[1];
			}
	   t_NNout->Fill();
		}
#ifdef RHO_CC
           if (type==1&&resto==0){
		sigTree.GetEntry(i+q);
		h2SigO_R->Fill(signal.rho[0],output);
		for (int ii=0; ii<NRHO;ii++){
				if (signal.rho[0]>=rho1[ii]) {
			for (int ji=0;ji<NCC;ji++){
				  if (signal.netcc[1]>=cc1[ji]) { 
					if (output>0.5){ pur[ii*NRHO+ji]=pur[ii*NRHO+ji]+1;
							 tot[ii*NRHO+ji]=tot[ii*NRHO+ji]+1;
						}
					}
				}
			}
		}
  		if (output<0.1) cout<<"output<0.1, while type=1 event_index: "<<i+q<<" of file: "<<SIG_FILE.Data()<<endl;
		h2SigR_C->Fill(signal.netcc[1],signal.rho[0]);
		h2Sig->Fill(signal.netcc[1],output);
		}
          if (type==0&&resto==0){
			bgTree.GetEntry(p-sig_entries+i-nTrain);
			h2BgO_R->Fill(background.rho[0],output);
		for (int ii=0; ii<NRHO;ii++){
			for (int ji=0;ji<NCC;ji++){
				if (background.rho[0]>rho1[ii]) {
				  if (background.netcc[1]>cc1[ji]) { 
					if (output>0.5) { tot[ii*NRHO+ji]=tot[ii*NRHO+ji]+1;
						}
					}
				}
			}
		}
  			if (output>0.9) cout<<"output>0.9, while type=0 event_index: "<<p-sig_entries+i-nTrain<<" of file: "<<BG_FILE.Data()<<endl;
			h2BgR_C->Fill(background.netcc[1],background.rho[0]);
			h2Bg->Fill(background.netcc[1],output);
			}
#endif
       // cout<<"output "<<output<<" type "<<type;
       // cout<<" #entry "<<i<<endl;
    int entries=nTrain; 
    if(resto==0){
       while(output<thres0[j]&&(j>0)) j=j-1;
        for (int k=npunti;k>=0;k--){
                if(type==0) {
                        if(k<=j) FA[k]=FA[k]+1;}
                //      cout<<"FA di "<<k<<" : "<<FA[k]<<" relativa a soglia "<<thres0[k]<<endl;}
                if(type==1) {
                        if(k<=j) TA[k]=TA[k]+1; }
        }
         tmedio=tmedio+type;
         ymedio=ymedio+output;
         yy=output*output+yy;
         tt=type*type+tt;
         ty=type*output+ty;
     }
   if(resto!=0){
       while(output<thres0[j]&&(j>0)) j=j-1;
        for (int k=npunti;k>=0;k--){
                if(type==0) {
                        if(k<=j) FA1[k]=FA1[k]+1;}
                //      cout<<"FA di "<<k<<" : "<<FA[k]<<" relativa a soglia "<<thres0[k]<<endl;}
                if(type==1) {
                        if(k<=j) TA1[k]=TA1[k]+1; }
        }
         tmedio1=tmedio1+type;
         ymedio1=ymedio1+output;
         yy1=output*output+yy1;
         tt1=type*type+tt1;
         ty1=type*output+ty1;
     }

    }
/*   h2Bg->SetLineColor(kBlue);
   h2Bg->SetFillStyle(3008);   ibg->SetFillColor(kBlue);
   h2Sig->SetLineColor(kRed);
   h2Sig->SetFillStyle(3003); isig->SetFillColor(kRed);
*/

   TString NNoutName(nomeNN);
   NNoutName.ReplaceAll("NN/","NNout/NNout_");
    TFile* f_NNout=new TFile(NNoutName,"RECREATE");
   t_NNout->Write();
   f_NNout->Close();
   
cout<<""<<endl;
   for (int ii=0;ii<NRHO; ii++) for (int ji=0;ji<NCC; ji++) {
        cout<<" pur "<<pur[NRHO*ii+ji]<<" tot "<<tot[NRHO*ii+ji]<<" cc "<<cc1[ji]<< " rho "<<rho1[ii];
//<<endl;
	if(tot[NRHO*ii+ji]>0)  pur[NRHO*ii+ji]=pur[NRHO*ii+ji]/tot[NRHO*ii+ji];
        cout<<" pur_fin "<<pur[NRHO*ii+ji];
        //cout<<" index "<<NRHO*ii+ji;
        cout<<""<<endl;
	//if(tot[NRHO*ii+ji]>0) PUREZZA_S->Fill(cc1[ji],rho1[ii],pur[NRHO*ii+ji]);
	if(tot[NRHO*ii+ji]>0) PUREZZA_S->Fill(cc1[ji]+deltacc/2,rho1[ii]+deltarho/2,pur[NRHO*ii+ji]);
	} 
  h2Bg->SetMarkerStyle(7);
   h2Sig->SetMarkerStyle(6);
   h2Bg->SetMarkerColor(4);
   h2Sig->SetMarkerColor(2);
   h2BgO_R->SetMarkerStyle(7);
   h2SigO_R->SetMarkerStyle(6);
   h2BgO_R->SetMarkerColor(4);
   h2SigO_R->SetMarkerColor(2);
   h2BgR_C->SetMarkerStyle(7);
   h2SigR_C->SetMarkerStyle(6);
   h2BgR_C->SetMarkerColor(4);
   h2SigR_C->SetMarkerColor(2);
  // h2Bg->SetStats(0);
  // h2Bg->SetStats(0);
   
   TCanvas *cRHO_CC_OUT=new TCanvas("RHO_CC_NNOUT","RHO_CC_NNOUT",1200,700);
   //cRHO_CC_OUT->Divide(3,1);
   cRHO_CC_OUT->Divide(2,2);
   cRHO_CC_OUT->cd(1)->SetTitle("OUT_CC");
   h2Bg->GetXaxis()->SetTitle("cc");
   h2Bg->GetYaxis()->SetTitle("NNoutput");
   h2Bg->GetZaxis()->SetTitle("count");
   h2Bg->Draw();
   h2Sig->Draw("same");
   cRHO_CC_OUT->cd(2)->SetTitle("OUT_RHO");
   h2BgO_R->GetXaxis()->SetTitle("rho");
   h2BgO_R->GetYaxis()->SetTitle("NNoutput");
   h2BgO_R->GetZaxis()->SetTitle("count");
   h2BgO_R->Draw();
   h2SigO_R->Draw("same");
   cRHO_CC_OUT->cd(3)->SetTitle("RHO_CC");
   h2BgR_C->GetXaxis()->SetTitle("cc");
   h2BgR_C->GetYaxis()->SetTitle("rho");
   h2BgR_C->GetZaxis()->SetTitle("count");
   h2BgR_C->Draw();
   h2SigR_C->Draw("same");
   cRHO_CC_OUT->cd(4)->SetTitle("PUREZZA");
   PUREZZA_S->SetStats(0);
   PUREZZA_S->GetXaxis()->SetTitle("cc");
   PUREZZA_S->GetYaxis()->SetTitle("rho");
   PUREZZA_S->GetZaxis()->SetTitle("count");
   PUREZZA_S->Draw("colz");
//   PUREZZA_B->Draw("same");
   TString cRHO_CC_OUTname(nomeNN);
   TString cRHO_CC_OUTnameroot(nomeNN);
   cRHO_CC_OUTname.ReplaceAll(".root",".png");
   cRHO_CC_OUTnameroot.ReplaceAll("NN/mlpNetwork","CC_OUT/RHO_CC_OUTgraph");
   cRHO_CC_OUTname.ReplaceAll("NN/mlpNetwork","CC_OUT/RHO_CC_OUTgraph");
   cRHO_CC_OUT->SaveAs(cRHO_CC_OUTname.Data());
   cRHO_CC_OUT->Print(cRHO_CC_OUTnameroot.Data());
//        }
   tmedio1=tmedio1/nTrain;
    ymedio1=ymedio1/nTrain;
//    tt1=tt1/nTrain;
//    yy1=yy1/nTrain;
//    ty1=ty1/nTrain;
    tmedio=tmedio/nTrain;
    ymedio=ymedio/nTrain;
//    tt=tt/nTrain;
//    yy=yy/nTrain;
//    ty=ty/nTrain;
    cout<<"tmedio1: "<<tmedio1<<" ymedio1 "<<ymedio1<<" tt1 "<<tt1<<" yy1 "<<yy1<<" ty1 "<<ty1<<" entries "<<entries<<endl;
    double corr=0.;
    double corr1=0.;
    corr1=(ty1-nTrain*ymedio1*tmedio1)/sqrt((tt1-nTrain*tmedio1*tmedio1)*(yy1-nTrain*ymedio1*ymedio1));
    corr=(ty-nTrain*ymedio*tmedio)/sqrt((tt-nTrain*tmedio*tmedio)*(yy-nTrain*ymedio*ymedio));
    cout<<"corr "<<corr<<endl;
    cout<<"corr1 "<<corr1<<endl;  

  
  for (int k=0;k<=npunti;k++) {
      //  cout<<"su "<<nTrain/2<<" TA: "<<TA[k]<<endl;
      //  cout<<"su "<<nTrainBg/2<<" FA: "<<FA[k]<<endl;
        TA1[k]=(TA1[k]/nTrain)*2;
        TA[k]=(TA[k]/nTrain)*2;
        //cout<<"Eff: "<<TA[k]<<endl;
        FA1[k]=(FA1[k]/nTrainBg)*2;
        FA[k]=(FA[k]/nTrainBg)*2;
         }
    TString nf(nomeNN);
    cout<<"ng "<<nf.Data()<<endl;
    char nomefileTest[516];
    int npunti2=npunti;
    //sprintf(nomefile,"thresFiles/thres_np%i_nTS%i_nTB%i.root",npunti,nTrainS,nTrainB);
    sprintf(nomefileTest,"thresFiles/thres_np%i_Test_onTrainSet_",npunti2);
    nf.ReplaceAll("NN/mlpNetwork",nomefileTest);
    cout<<"nf "<<nf.Data()<<endl;
    TFile *gfile=new TFile(nf.Data(),"RECREATE");

//#ifdef ROC_GRAPH
    TGraph* gROC1=new TGraph(npunti+1,FA1,TA1);
    gROC1->SetTitle("ROC_train(threshold)");
    gROC1->SetMarkerStyle(7);
    TGraph* gFA1=new TGraph(npunti+1,thres0,FA1);
    gFA1->SetTitle("False Alarm _train(threshold)");
    gFA1->SetMarkerStyle(7);
    TGraph* gTA1=new TGraph(npunti+1,thres0,TA1);
    gTA1->SetTitle("Efficiency _train(threshold)");
    gTA1->SetMarkerStyle(7);

    TGraph* gROC=new TGraph(npunti+1,FA,TA);
    gROC->SetTitle("ROC_test(threshold)");
    gROC->SetMarkerStyle(7);
    TGraph* gFA=new TGraph(npunti+1,thres0,FA);
    gFA->SetTitle("False Alarm _test(threshold)");
    gFA->SetMarkerStyle(7);
    TGraph* gTA=new TGraph(npunti+1,thres0,TA);
    gTA->SetTitle("Efficiency _test(threshold)");
    gTA->SetMarkerStyle(7);

    TCanvas *cROC=new TCanvas("ROC_graph","ROC_graph");
    cROC->Divide(2,3);
    cROC->cd(2)->SetLogx();
   // cROC->cd(2)->SetLogy();
    gROC->Draw("apl");
    cROC->cd(6);
    cROC->cd(6)->SetLogy();
    gFA->Draw("apl");
    cROC->cd(4);
    cROC->cd(4)->SetLogy();
    gTA->Draw("apl");
    cROC->cd(1)->SetLogx();
   // cROC->cd(1)->SetLogy();
    gROC1->Draw("apl");
    cROC->cd(5);
    cROC->cd(5)->SetLogy();
    gFA1->Draw("apl");
    cROC->cd(3);
    cROC->cd(3)->SetLogy();
    gTA1->Draw("apl");

    cROC->Write();
    TString CANV_NAME(nomeNN);
    TString CANV_NAMEroot(nomeNN);
    CANV_NAME.ReplaceAll("NN/mlpNetwork","TRAIN_GRAPHS/ROC");
    CANV_NAMEroot.ReplaceAll("NN/mlpNetwork","TRAIN_GRAPHS/ROC");
    CANV_NAME.ReplaceAll(".root",".png");
    cROC->SaveAs(CANV_NAME.Data());
    cROC->Print(CANV_NAMEroot.Data());
//#endif
//#ifdef TH_TREE
    TTree *ThTree=new TTree("ROC","ROC");
char FAlabel1[128];
char TAlabel1[128];
char FAlabel[128];
char TAlabel[128];
char thlabel[128];
sprintf(FAlabel1,"FalseAlarm_train[%i]/F",npunti);
sprintf(TAlabel1,"TrueAlarm_train[%i]/F",npunti);
sprintf(FAlabel,"FalseAlarm_test[%i]/F",npunti);
sprintf(TAlabel,"TrueAlarm_test[%i]/F",npunti);
sprintf(thlabel,"thresholds[%i]/F",npunti);
//    float FAl=0;
    ThTree->Branch("FalseAlarm_test",&FA,FAlabel);
    ThTree->Branch("FalseAlarm_train",&FA1,FAlabel1);
//    float TAl=0;
    ThTree->Branch("Efficiency_train",&TA1,TAlabel1);
    ThTree->Branch("Efficiency_test",&TA,TAlabel);
//    double thresholds=0;
    ThTree->Branch("thresholds",&thres0,thlabel);
    ThTree->Branch("correlation_train",&corr1,"correlation_train/D");
    ThTree->Branch("correlation_test",&corr,"correlation_test/D");
//    ThTree->Branch("#TrainNN_Sig",&nTrainSNN,"#TrainNN_Sig/I");
    int u=npunti;
    ThTree->Branch("#points",&u,"#points/I");
//    ThTree->Branch("#TrainNN_Bg",&nTrainBNN,"#TrainNN_Bg/I");
    ThTree->Fill();
    ThTree->Write();
//#endif

    gfile->Close();
}




TString select7Train_non_es(TString ORIGINAL_FILE,std::vector<int>* choice,TString nomeNN) {
TString IDname(nomeNN);

//IDname.ReplaceAll("NN/","PropertyGraphs/");

TFile* fOriginal=TFile::Open(ORIGINAL_FILE.Data());
TTree* tOriginal=(TTree*)fOriginal->Get("nnTree");

if((choice->size())>(tOriginal->GetEntries())){cout<<"Error:startB+nTrainB>Entries"<<endl;
                                     // exit(0);
                                        }

//estrazione nINP
int ninp;
tOriginal->SetBranchAddress("#inputs",&ninp);
tOriginal->GetEntry(0);

int const nINP=ninp;
int const NDIM=sqrt(nINP);

//definizione nuovo tree
TTree* t=new TTree("nnTree","nnTree");

int type;
tOriginal->SetBranchAddress("type",&type);
t->Branch("type",&type,"type/I");

float x[nINP];
char xlabel[nINP][16];
char xllabel[nINP][16];
#ifdef defferentW
float w[nINP];
char wlabel[nINP][16];
char wllabel[nINP][16];
#else
float w;
tOriginal->SetBranchAddress("w",&w);
t->Branch("w",&w,"w/F");
#endif

#ifdef trainTEST
float SDs[nINP];
char SDlabels[nINP][16];
char SDllabels[nINP][16];
float RMSs[nINP];
char RMSlabels[nINP][16];
char RMSllabels[nINP][16];
float averages[nINP];
char averagelabels[nINP][16];
char averagellabels[nINP][16];
float SDb[nINP];
char SDlabelb[nINP][16];
char SDllabelb[nINP][16];
float RMSb[nINP];
char RMSlabelb[nINP][16];
char RMSllabelb[nINP][16];
float averageb[nINP];
char averagelabelb[nINP][16];
char averagellabelb[nINP][16];
float SD[nINP];
char SDlabel[nINP][16];
char SDllabel[nINP][16];
float RMS[nINP];
char RMSlabel[nINP][16];
char RMSllabel[nINP][16];
float average[nINP];
char averagelabel[nINP][16];
char averagellabel[nINP][16];

TTree *tTrainTest=new TTree("trainTEST","trainTEST");
#endif

for (int i=0; i<nINP;i++){
   #ifdef trainTEST
    RMS[i]=0;
    average[i]=0;
    SD[i]=0;
    RMSs[i]=0;
    averages[i]=0;
    SDs[i]=0;
    RMSb[i]=0;
    averageb[i]=0;
    SDb[i]=0;
  #endif
    sprintf(xlabel[i],"x%i",i+1);
    sprintf(xllabel[i],"%s/F",xlabel[i]);
    tOriginal->SetBranchAddress(xlabel[i],&x[i]);
    t->Branch(xlabel[i],&x[i],xllabel[i]);
    #ifdef differentW
    sprintf(wlabel[i],"w%i",i+1);
    sprintf(wllabel[i],"%s/F",wlabel[i]);
    tOriginal->SetBranchAddress(wlabel[i],&w[i]);
    t->Branch(wlabel[i],w[i],wllabel[i]);
    #endif
    }
int size=0;
size=choice->size();
int nTSB=0;
nTSB=size/2;
for (int k=0;k<choice->size();k++){
    //cout<<"choice di "<<k<<" : "<<(*choice)[k]<<endl;
    tOriginal->GetEntry((*choice)[k]);
    #ifdef trainTEST
      for (int i=0;i<nINP;i++){
        RMS[i]=RMS[i]+x[i]*x[i];
        average[i]=average[i]+x[i];
	if(k<nTSB){
        RMSs[i]=RMSs[i]+x[i]*x[i];
        averages[i]=averages[i]+x[i];
       	}
	if(k>=nTSB){
        RMSb[i]=RMSb[i]+x[i]*x[i];
        averageb[i]=averageb[i]+x[i];
       	}
        }
    #endif
    for (int u=0;u<k;u++) if( (*choice)[u]==(*choice)[k]) {cout<<"Error:same event taken 2 times"<<"u "<<u<<" choise "<<(*choice)[u]<<endl; exit(0);}
//  h[k]=(*choice)[k];
    t->Fill();
    }
/*
for (int k=b;k<b+nTrainB;k++){
    tOriginal->GetEntry(k);
    t->Fill();

*/
//cout<<t->GetEntries();
int Ntest=choice->size();
#ifdef trainTEST
/*TGraph *RMSg=new TGraph();
TGraph *AVg=new TGraph();
TGraph *SDg=new TGraph();
*/
TH2D *RMShs=new TH2D("RMShs","RMShs",NDIM,0,NDIM,NDIM,0,NDIM);
TH2D *AVhs=new TH2D("AVhs","Avhs",NDIM,0,NDIM,NDIM,0,NDIM);
TH2D *SDhs=new TH2D("SDhs","SDhs",NDIM,0,NDIM,NDIM,0,NDIM);
TH2D *RMShb=new TH2D("RMShb","RMShb",NDIM,0,NDIM,NDIM,0,NDIM);
TH2D *AVhb=new TH2D("AVhb","Avhb",NDIM,0,NDIM,NDIM,0,NDIM);
TH2D *SDhb=new TH2D("SDhb","SDhb",NDIM,0,NDIM,NDIM,0,NDIM);
TH2D *RMSh=new TH2D("RMSh","RMSh",NDIM,0,NDIM,NDIM,0,NDIM);
TH2D *AVh=new TH2D("AVh","Avh",NDIM,0,NDIM,NDIM,0,NDIM);
TH2D *SDh=new TH2D("SDh","SDh",NDIM,0,NDIM,NDIM,0,NDIM);
//TH1D *RMSh=new TH1D("RMSh","RMSh",nINP,0,nINP-1);
//TH1D *AVh=new TH1D("AVh","Avh",nINP,0,nINP-1);
//TH1D *SDh=new TH1D("SDh","SDh",nINP,0,nINP-1);
RMSh->SetStats(0);
SDh->SetStats(0);
AVh->SetStats(0);
SDhb->SetStats(0);
AVhb->SetStats(0);
RMShb->SetStats(0);
SDhs->SetStats(0);
AVhs->SetStats(0);
RMShs->SetStats(0);



for (int i=0;i<nINP;i++){
    RMSs[i]=RMSs[i]/Ntest*2;
    averages[i]=averages[i]/Ntest*2;
    SDs[i]=RMSs[i]-pow(averages[i],2);
    SDs[i]=pow(SDs[i],0.5);
    RMSs[i]=pow(RMSs[i],0.5);

    cout << "SIGNAL: " << averages[i] << " " << RMSs[i] << " " << SDs[i] << endl;

    RMSb[i]=RMSb[i]/Ntest*2;
    averageb[i]=averageb[i]/Ntest*2;
    SDb[i]=RMSb[i]-pow(averageb[i],2);
    SDb[i]=pow(SDb[i],0.5);
    RMSb[i]=pow(RMSb[i],0.5);

    cout << "BKG: " << averageb[i] << " " << RMSb[i] << " " << SDb[i] << endl;

    RMS[i]=RMS[i]/Ntest;
    average[i]=average[i]/Ntest;
    SD[i]=RMS[i]-pow(average[i],2);
    SD[i]=pow(SD[i],0.5);
    RMS[i]=pow(RMS[i],0.5);

    cout << "ALL: " << average[i] << " " << RMS[i] << " " << SD[i] << endl;

    sprintf(SDlabels[i],"SDs%i",i+1);
    sprintf(SDllabels[i],"%s/F",SDlabels[i]);
    tTrainTest->Branch(SDlabels[i],&SDs[i],SDllabels[i]);
    sprintf(RMSlabels[i],"RMSs%i",i+1);
    sprintf(RMSllabels[i],"%s/F",RMSlabels[i]);
    tTrainTest->Branch(RMSlabels[i],&RMSs[i],RMSllabels[i]);
    sprintf(averagelabels[i],"averages%i",i+1);
    sprintf(averagellabels[i],"%s/F",averagelabels[i]);
    tTrainTest->Branch(averagelabels[i],&averages[i],averagellabels[i]);
    sprintf(SDlabelb[i],"SDb%i",i+1);
    sprintf(SDllabelb[i],"%s/F",SDlabelb[i]);
    tTrainTest->Branch(SDlabelb[i],&SDb[i],SDllabelb[i]);
    sprintf(RMSlabelb[i],"RMSb%i",i+1);
    sprintf(RMSllabelb[i],"%s/F",RMSlabelb[i]);
    tTrainTest->Branch(RMSlabelb[i],&RMSb[i],RMSllabelb[i]);
    sprintf(averagelabelb[i],"averageb%i",i+1);
    sprintf(averagellabelb[i],"%s/F",averagelabelb[i]);
    tTrainTest->Branch(averagelabelb[i],&averageb[i],averagellabelb[i]);
    sprintf(SDlabel[i],"SD%i",i+1);
    sprintf(SDllabel[i],"%s/F",SDlabel[i]);
    tTrainTest->Branch(SDlabel[i],&SD[i],SDllabel[i]);
    sprintf(RMSlabel[i],"RMS%i",i+1);
    sprintf(RMSllabel[i],"%s/F",RMSlabel[i]);
    tTrainTest->Branch(RMSlabel[i],&RMS[i],RMSllabel[i]);
    sprintf(averagelabel[i],"average%i",i+1);
    sprintf(averagellabel[i],"%s/F",averagelabel[i]);
    tTrainTest->Branch(averagelabel[i],&average[i],averagellabel[i]);
   /* RMSg->SetPoint(i,i,RMS[i]);
    AVg->SetPoint(i,i,average[i]);
    SDg->SetPoint(i,i,SD[i]);
    */
    int fi=0;
    int ti=0;
    ti=i/NDIM;
    fi=i-ti*NDIM;
//    cout<<"ti: "<<ti<<" fi: "<<fi<<endl;
    RMShs->Fill(ti,fi,RMSs[i]);
    AVhs->Fill(ti,fi,averages[i]);
    SDhs->Fill(ti,fi,SDs[i]);
    RMShb->Fill(ti,fi,RMSb[i]);
    AVhb->Fill(ti,fi,averageb[i]);
    SDhb->Fill(ti,fi,SDb[i]);
    RMSh->Fill(ti,fi,RMS[i]);
    AVh->Fill(ti,fi,average[i]);
    SDh->Fill(ti,fi,SD[i]);
//    RMSh->Fill(i,RMS[i]);
//    AVh->Fill(i,average[i]);
//    SDh->Fill(i,SD[i]);
    }
TCanvas *cProp=new TCanvas ("TrainSet_Properties","TrainSet_Properties",0,0,900,600);
cProp->Divide(3,3);
cProp->cd(1);
RMSh->SetTitle("RMS");
RMSh->Draw("COLZ");
cProp->cd(2);
AVh->SetTitle("Average");
AVh->Draw("COLZ");
cProp->cd(3);
SDh->SetTitle("SD");
SDh->Draw("COLZ");
cProp->cd(4);
RMShs->SetTitle("RMS_s");
RMShs->Draw("COLZ");
cProp->cd(5);
AVhs->SetTitle("Average_s");
AVhs->Draw("COLZ");
cProp->cd(6);
SDhs->SetTitle("SD_s");
SDhs->Draw("COLZ");
cProp->cd(7);
RMShb->SetTitle("RMS_b");
RMShb->Draw("COLZ");
cProp->cd(8);
AVhb->SetTitle("Average_b");
AVhb->Draw("COLZ");
cProp->cd(9);
SDhb->SetTitle("SD_b");
SDhb->Draw("COLZ");

TString Prop_name(IDname);
TString Prop_nameroot(IDname);
Prop_name.ReplaceAll(".root","Prop.png");
Prop_nameroot.ReplaceAll(".root","Prop.root");
cProp->SaveAs(Prop_name.Data());
cProp->Print(Prop_nameroot.Data());
cProp->Close();
#endif
int ent=t->GetEntries();
if (choice->size()!=ent) cout<<"Error in Tree dimension: Entries new tree="<<t->GetEntries()<<" vector->size"<<choice->size()<<endl;

//char nomefile[516];
//sprintf(nomefile,"selectTREE/TreeTrain_size%i_s%i_b%i.root",choice->size(),s,b);
TString NOME_FILE(IDname);
NOME_FILE.ReplaceAll("PropertyGraphs/mlpNetwork","selectTREE/TreeTrain");
TFile* ofile=new TFile(NOME_FILE.Data(),"RECREATE");
t->SetDirectory(ofile);
t->Write();
ofile->Close();
fOriginal->Close();

#ifdef trainTEST
/*TCanvas *matrixIndex_RMS=new TCanvas("matrixIndex_RMS","matrixIndex_RMS",0,0,600,600);
tTrainTest->Draw("RMS");
matrixIndex*/
TFile* TrainTest= new TFile(IDname.Data(),"RECREATE");
tTrainTest->Write();
TrainTest->Close();
#endif
cout<<NOME_FILE.Data()<<endl;
return NOME_FILE;
}

                                                  
