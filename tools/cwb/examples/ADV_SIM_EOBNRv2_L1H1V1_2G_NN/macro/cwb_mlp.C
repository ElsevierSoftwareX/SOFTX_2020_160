//
// Neural Network Analysis
// Author : Gabriele Vedovato
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


#define SIG_FILE "data/nn_931158208_60_ADV_SIM_EOBNRv2_L1H1V1_2G_NN_17.3205_job1.root"
#define BG_FILE  "insert bckground file !!!"

#define nIFO  3		// number of detectors

#define NDIM 6		
#define nINP 36

#define CREATE_NETWORK

void cwb_mlp(Int_t ntrain=100) {
   
   if (!gROOT->GetClass("TMultiLayerPerceptron")) {
      gSystem->Load("libMLP");
   }

   // SIGNAL
   TChain sigTree("waveburst");
   sigTree.Add(SIG_FILE);
   netevent signal(&sigTree,nIFO);
   int sig_entries = signal.GetEntries();
   cout << "sig entries : " << sig_entries << endl;
   std::vector<double>* sig_frame = new vector<double>;
   signal.fChain->SetBranchAddress("nnframe", &sig_frame);

   // BACKGROUND
   TChain bgTree("waveburst");
   bgTree.Add(BG_FILE);
   netevent background(&bgTree,nIFO);
   int bg_entries = background.GetEntries();
   cout << "bg entries : " << bg_entries << endl;
   std::vector<double>* bg_frame = new vector<double>;
   background.fChain->SetBranchAddress("nnframe", &bg_frame);

   // MPL ROOT FILE
   // The 2 trees are merged into one, and a "type" branch, 
   // equal to 1 for the signal and 0 for the background is added.
   TTree* nnTree = (TTree*)sigTree.CloneTree(0);
   netevent nn(nnTree,nIFO);
   nnTree->SetDirectory(0);
   // add leaf
   Int_t type;
   float s[nINP];
   Float_t x[nINP],w;
   char ilabel[nINP][16];   
   for(int i=0;i<nINP;i++) {
     sprintf(ilabel[i],"x%d",i+1); 
     char tlabel[16]; sprintf(tlabel,"x%d/F",i+1); 
     nnTree->Branch(ilabel[i], &x[i], tlabel);
   }
   nnTree->Branch("w",&w,"w/F");
   nnTree->Branch("type",&type,"type/I");

   w = 1;
   TH1F *ibg = new TH1F("bgh", "NN output", 50, 0, 18);
   TH1F *isig = new TH1F("sigh", "NN output", 50, 0, 18);
   ibg->SetDirectory(0);
   isig->SetDirectory(0);
   Int_t i;

   // ADD SIGNAL
   type = 1;
   for(int n=0;n<sig_entries;n++) {
     sigTree.GetEntry(n);
     for (int j=0;j<nINP;j++) x[j]=(*sig_frame)[j];
     isig->Fill(signal.rho[0]);
     nnTree->Fill();
   }

   // ADD BACKGROUND
   type = 0;
   for(int n=0;n<bg_entries;n++) {
     bgTree.GetEntry(n);
     for (int j=0;j<nINP;j++) x[j]=(*bg_frame)[j];
     isig->Fill(background.rho[0]);
     nnTree->Fill();
   }

   // Define MultiLayer Perceptron layout
   char layout[256]="";
   for (int i=0;i<nINP-1;i++) sprintf(layout,"%s%s,",layout,ilabel[i]);
   sprintf(layout,"%s%s:10:40:10:type",layout,ilabel[nINP-1]);
   cout << "MultiLayer Perceptron layout -> " << layout << endl;

#ifdef CREATE_NETWORK
   TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron(layout,"w",nnTree,"Entry$%2","(Entry$+1)%2");
   //mlp->SetLearningMethod(TMultiLayerPerceptron::kStochastic);
   //mlp->SetLearningMethod(TMultiLayerPerceptron::kBatch);
   mlp->SetLearningMethod(TMultiLayerPerceptron::kSteepestDescent);
   //mlp->SetLearningMethod(TMultiLayerPerceptron::kRibierePolak);
   //mlp->SetLearningMethod(TMultiLayerPerceptron::kFletcherReeves);
   //mlp->SetLearningMethod(TMultiLayerPerceptron::kBFGS);
   mlp->Train(ntrain, "text,graph,update=10");

   TFile fnet("mplNetwork.root","RECREATE");
   mlp->Write();
   fnet.Close(); 
#else
   TFile* fnet = new TFile("mplNetwork.root");
   TMultiLayerPerceptron *mlp = (TMultiLayerPerceptron*)fnet->Get("TMultiLayerPerceptron");
   if(mlp==NULL) {cout << "Error getting mlp" << endl;exit(1);}
   //fnet->Close(); 
#endif

   // Use TMLPAnalyzer to see what it looks for
   TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis");
   mlpa_canvas->Divide(2,2);

   mlpa_canvas->cd(1);
   ibg->SetLineColor(kBlue);
   ibg->SetFillStyle(3008);   ibg->SetFillColor(kBlue);
   isig->SetLineColor(kRed);
   isig->SetFillStyle(3003); isig->SetFillColor(kRed);
   ibg->SetStats(0);
   isig->SetStats(0);
   ibg->Draw();
   isig->Draw("same");
   TLegend *ilegend = new TLegend(.75, .80, .95, .95);
   ilegend->AddEntry(ibg, "Background");
   ilegend->AddEntry(isig, "Signal");
   ilegend->Draw();
   cout << "Integral ibg : " << ibg->GetEntries() << endl;
   cout << "Integral isig : " << isig->GetEntries() << endl;

#ifdef CREATE_NETWORK
   // Build and train the NN ptsumf is used as a weight since we are primarly 

   TMLPAnalyzer ana(mlp);
   // Initialisation
   ana.GatherInformations();
   // output to the console
   ana.CheckNetwork();
   mlpa_canvas->cd(1);
   // shows how each variable influences the network
   //ana.DrawDInputs();
   mlpa_canvas->cd(2);
   // shows the network structure
   mlp->Draw();
   mlpa_canvas->cd(3);
   // draws the resulting network
   gPad->SetLogy();
   ana.DrawNetwork(0,"type==1","type==0");
   mlpa_canvas->cd(4);

#endif

   TH1F *bg = new TH1F("bgh", "NN output", 50, -0.5, 1.5);
   TH1F *sig = new TH1F("sigh", "NN output", 50, -0.5, 2.5);
   bg->SetDirectory(0);                                     
   sig->SetDirectory(0);                                    
   sig->Reset();                                            
   bg->Reset();                                             
   sig->Clear();                                            
   bg->Clear();                                             
   Double_t params[nINP];                                   
   cout << "#Entries : " << nnTree->GetEntries() << endl;     
   for (i = 0; i < nnTree->GetEntries(); i++) {               
     nnTree->GetEntry(i);
     for (int j=0;j<nINP;j++) params[j]=x[j];             
     float output = mlp->Evaluate(0,params);
     //if(type==0) for (int j=0;j<nINP;j++) cout << i << " " << x[j] << endl;             
     if(nn.netcc[1]>0.7) { 
       if(type==1) sig->Fill(output);
       else         bg->Fill(output);
     }
     if(type==0) {
       if((output>0.2)&&(nn.netcc[1]>0.7)) 
         cout << i << " " << mlp->Evaluate(0,params) << " " << nn.rho[0] << " " << nn.netcc[1] << endl;             
     }
   }

#ifndef CREATE_NETWORK
   mlpa_canvas->Divide(1,1);
   mlpa_canvas->cd(0);
#endif

   gPad->SetLogy();
   bg->SetLineColor(kBlue);
   bg->SetFillStyle(3008);   bg->SetFillColor(kBlue);
   sig->SetLineColor(kRed);
   sig->SetFillStyle(3003); sig->SetFillColor(kRed);
   bg->SetStats(0);
   sig->SetStats(0);
   sig->Draw();
   bg->Draw("same");
   TLegend *legend = new TLegend(.75, .80, .95, .95);
   legend->AddEntry(bg, "Background");
   legend->AddEntry(sig, "Signal");
   legend->Draw();
   mlpa_canvas->cd(0);
}
