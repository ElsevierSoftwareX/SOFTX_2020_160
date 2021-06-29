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


/*
   nn_otree->Branch("dt_corepixels",&dt,NPIX,0);
   //char deltat_l[1024];
   //sprintf(deltat_l,"dt_corepixels[%i]",NPIX);
   //nn_otree->Branch("dt_corepixels",&deltat,deltat_l);
   nn_otree->Branch("index",&ind,"ind/I");
   nn_otree->Branch("t_duration",&dT,"duration/D");
   nn_otree->Branch("f_duration",&dF,"f_duration/D");
   nn_otree->Branch("central_f",&Fc,"central_f/D");

*/



//#define differentW

//#define SIG_FILE "/home/serena.vinciguerra/test0/S6D_R11_SIM_EOBNRv2_advL1H1V1_2G_run1_2/merge/NN_wave_S6D_R11_SIM_EOBNRv2_advL1H1V1_2G_run1_2.M1.root"
//#define BG_FILE  "/home/serena.vinciguerra/test0/S6D_R11_BKG_L1H1V1_2G_MP_run1/merge/NN_wave_S6D_R11_BKG_L1H1V1_2G_MP_run1.M1.root"
#define nIFO 3
#define NPIX 20
using namespace std;

void TreeTest_CTFmod(TString ifileS,TString ifileB,TString ofile,int const y=1, int const NDIM=8,int const nINP=64){


/*if(ifileS!=0)
//{
   // SIGNAL
   TChain TreeS("waveburst");//cerca il Tree "waveburst nei file
   TreeS.Add(ifileS.Data());//determina i file
   netevent entryS(&TreeS,nIFO);
   int entriesS = entryS.GetEntries();
   cout << "Sig entries : " << entriesS << endl;
   std::vector<double>* frameS = new vector<double>;
   entryS.fChain->SetBranchAddress("nnframe", &frameS);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
//}

//if(ifileB!=0)
//{

   TChain TreeB("waveburst");//cerca il Tree "waveburst nei file
   TreeB.Add(ifileB.Data());//determina i file
   netevent entryB(&TreeB,nIFO);
   int entriesB = entryB.GetEntries();
   cout << "Bg entries : " << entriesB << endl;
   std::vector<double>* frameB = new vector<double>;
   entryB.fChain->SetBranchAddress("nnframe", &frameB);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
// }
*/

   //cout<<"frame "<<(*frame)[0]<<endl;
   TTree* nnTree2 = new TTree ("nnTree","nnTree");
//   netevent nn(Tree,nIFO);
  nnTree2->SetDirectory(0);

   // add leaf
   Int_t type;
   float s[nINP];
   Float_t x[nINP];
   for (int jj=0; jj<nINP;jj++) x[jj]=0.;
   char ilabel[nINP][16];
   double dt[NPIX];
   for (int jj=0; jj<NPIX;jj++) dt[jj]=0.;
   char ilabeldt[NPIX][16];
   double dT=0.;
   double Fc=0.;
   double dF=0.;

   //costruzione una foglia  per ogni input
   for(int i=0;i<nINP;i++) {
     sprintf(ilabel[i],"x%i",i+1);
     char tlabel[16]; sprintf(tlabel,"x%i/F",i+1);
     nnTree2->Branch(ilabel[i], &x[i], tlabel);
   }

//costruzione una foglia per ogni dt
   for(int i=0;i<NPIX;i++) {
     sprintf(ilabeldt[i],"dt%i",i+1);
     char tlabeldt[16]; sprintf(tlabeldt,"dt%i/D",i+1);
     nnTree2->Branch(ilabeldt[i], &dt[i], tlabeldt);
   }
//add info time-frequncy plane
     nnTree2->Branch("duration", &dT, "duration/D");
     nnTree2->Branch("frequency_width", &dF, "frequency_width/D");
     nnTree2->Branch("central_frequency", &Fc, "central_frequency/D");
   
	
   //costruzione foglie peso e tipo
   nnTree2->Branch("type",&type,"type/I");
   int yi=0;
   yi=y;
   int nDim=0;
  nDim=NDIM;
   int nInp=0;
  nInp=nINP;
int index=0;
int entriesTot=0;
   nnTree2->Branch("amplitude_mode",&yi,"amplitude_mode/I");
  nnTree2->Branch("Matrix_dim",&nDim,"Matrix_dim/I");
   nnTree2->Branch("#inputs",&nInp,"#inputs/I");
   nnTree2->Branch("index",&index,"index/I");
   nnTree2->Branch("#Entries_type",&entriesTot,"#Entries_type/I");
char fname[516];
   nnTree2->Branch("Files_name",&fname,"Files_name/C");

#ifdef differentW
 char wlabel[nINP][16];
 float w[nINP];
 for(int i=0;i<nINP;i++) {
     sprintf(wlabel[i],"w%d",i+1);
     char llabel[16]; sprintf(llabel,"w%d/F",i+1);
     nnTree2->Branch(ilabel[i], &w[i], tlabel);
     w[i]=1;
     }

#else
 float w;
   nnTree2->Branch("w",&w,"w/F");
   w=1;
#endif

cout<<"fine definizione foglie"<<endl;
if(ifileS!=0)
{
   TChain TreeS("waveburst");//cerca il Tree "waveburst nei file
   TreeS.Add(ifileS.Data());//determina i file
   netevent entryS(&TreeS,nIFO);
   entriesTot = entryS.GetEntries();
   cout << "Sig entries : " << entriesTot << endl;
   std::vector<double>* frameS = new vector<double>;
   entryS.fChain->SetBranchAddress("nnframe", &frameS);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
   std::vector<double>* dtvS = new vector<double>;
   entryS.fChain->SetBranchAddress("dt_corepixels", &dtvS);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
   
   double dTS;
   double dFS;
   double FcS;

   entryS.fChain->SetBranchAddress("t_duration", &dTS);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
   entryS.fChain->SetBranchAddress("f_duration", &dFS);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
   entryS.fChain->SetBranchAddress("central_f", &FcS);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
 
  type = 1;
  sprintf(fname,"%s",ifileS.Data());
  cout<<fname<<endl;  
 for(int n=0;n<entryS.GetEntries();n++) {
     TreeS.GetEntry(n);
     dT=dTS;
     dF=dFS;
     Fc=FcS;
     for (int j=0;j<NPIX;j++) dt[j]=(*dtvS)[j];
     for (int j=0;j<nINP;j++)
     //for (int j=0;j<nINP;j++)
        {  if(yi==0){ 
               if((*frameS)[j]>0.001) x[j]=1;
                             else x[j]=0;
                  }
          //else x[j]=frameS->at(j);
          else x[j]=(*frameS)[j];
        }
     index=index+1;
     nnTree2->Fill();
	 }
}
cout<<"Fine ciclo sul segnale"<<endl;
if(ifileB!=0)
   {
   TChain TreeB("waveburst");//cerca il Tree "waveburst nei file
   TreeB.Add(ifileB.Data());//determina i file
   netevent entryB(&TreeB,nIFO);
   entriesTot = entryB.GetEntries();
   cout << "Bg entries : " << entriesTot << endl;
   std::vector<double>* frameB = new vector<double>;
   entryB.fChain->SetBranchAddress("nnframe", &frameB);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
   std::vector<double>* dtvB = new vector<double>;
   entryB.fChain->SetBranchAddress("dt_corepixels", &dtvB);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
  
   double dTB;
   double dFB;
   double FcB;

   entryB.fChain->SetBranchAddress("t_duration", &dTB);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
   entryB.fChain->SetBranchAddress("f_duration", &dFB);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
   entryB.fChain->SetBranchAddress("central_f", &FcB);//fa esattamente come ha fatto per le altre variabili nella macro netevent.cc
 
   type = 0;
 sprintf(fname,"%s",ifileB.Data());
 cout<<fname<<endl;
   for(int n=0;n<entryB.GetEntries();n++) {
     TreeB.GetEntry(n);
     dT=dTB;
     dF=dFB;
     Fc=FcB;
     for (int j=0;j<NPIX;j++) dt[j]=(*dtvB)[j];
     for (int j=0;j<nINP;j++)
     //for (int j=0;j<nINP;j++)
        {  if(yi==0){ 
               if((*frameB)[j]>0.001) x[j]=1;
                             else x[j]=0;
                  }
          //else x[j]=frameB->at(j);
          else x[j]=(*frameB)[j];
        }
     index=index+1;
     nnTree2->Fill();
	 }
}
cout<<"fine ciclo sul Bg"<<endl;
//SAVE THE TREE
/*
TString ofileS(ifileS.Data());
ofileS.ReplaceAll(".root","_nnTree");
TString ofileB(ifileB.Data());
ofileB.ReplaceAll(".root","_nnTree");
*/
char nomefile[516];
sprintf(nomefile,"nnTREE/nnTree_%s.root",ofile.Data());
//char nomefile[516]="nnTREE/nnTree_S6D_and_SIM0run3.root";
TFile* f=new TFile(nomefile,"RECREATE");
nnTree2->Write();
f->Close();

}
