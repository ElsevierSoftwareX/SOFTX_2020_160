/*
# Copyright (C) 2019 Gabriele Vedovato, Sergey Klimenko
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


#include "livetime.hh"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"

ClassImp(livetime)	// used by THtml doc

livetime& livetime::operator=(livetime& a)
{
  this->fChain=                NULL;             
  this->fCurrent=              a.fCurrent;                                                        
  this->run=                   a.run;             
  this->gps=                   a.gps;             
  this->live=                  a.live;              
  this->lag[NIFO_MAX+1]=       a.lag[NIFO_MAX+1];                 
  this->slag[NIFO_MAX+1]=      a.slag[NIFO_MAX+1];                 
  for(size_t n=0; n<(NIFO_MAX+1); n++) {
    this->lag[n]=              a.lag[n];                 
    this->slag[n]=             a.slag[n];                 
  }
  for(size_t n=0; n<NIFO_MAX; n++) {
    this->start[n]=            a.start[n];                 
    this->stop[n]=             a.stop[n];                 
  }
  return *this;
}

//   Set branch addresses
void livetime::Init(TTree *tree)
{
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",&(this->run));
   fChain->SetBranchAddress("gps",&(this->gps));
   fChain->SetBranchAddress("live",&(this->live));
   fChain->SetBranchAddress("lag",(this->lag));
   fChain->SetBranchAddress("slag",(this->slag));
   fChain->SetBranchAddress("start",(this->start));
   fChain->SetBranchAddress("stop",(this->stop));

   Notify();
}

void livetime::allocate()
{
   if (!lag)         lag=      (Float_t*)malloc((NIFO_MAX+1)*sizeof(Int_t));
   else              lag=      (Float_t*)realloc(lag,(NIFO_MAX+1)*sizeof(Int_t));
   if (!slag)       slag=      (Float_t*)malloc((NIFO_MAX+1)*sizeof(Int_t));
   else             slag=      (Float_t*)realloc(slag,(NIFO_MAX+1)*sizeof(Int_t));
   if (!start)     start=      (Double_t*)malloc((NIFO_MAX)*sizeof(Double_t));
   else            start=      (Double_t*)realloc(start,(NIFO_MAX)*sizeof(Double_t));
   if (!stop)       stop=      (Double_t*)malloc((NIFO_MAX)*sizeof(Double_t));
   else             stop=      (Double_t*)realloc(start,(NIFO_MAX)*sizeof(Double_t));
   return;
}

//++++++++++++++++++++++++++++++++++++++++++++++
// set noise livetime tree
//++++++++++++++++++++++++++++++++++++++++++++++
TTree* livetime::setTree()
{
   TTree* waveTree = new TTree("liveTime","liveTime");

 //==================================
 // Define livetime tree
 //==================================

   char clag[16];sprintf(clag,"lag[%d]/F",NIFO_MAX+1);
   char cslag[16];sprintf(cslag,"slag[%d]/F",NIFO_MAX+1);
   char cstart[16];sprintf(cstart,"start[%d]/D",NIFO_MAX);
   char cstop[16];sprintf(cstop,"stop[%d]/D",NIFO_MAX);

   waveTree->Branch("run",       &(this->run),        "run/I");
   waveTree->Branch("gps",       &(this->gps),        "gps/D");
   waveTree->Branch("live",      &(this->live),       "live/D");
   waveTree->Branch("lag",        (this->lag),        clag);
   waveTree->Branch("slag",       (this->slag),       cslag);
   waveTree->Branch("start",      (this->start),      cstart);
   waveTree->Branch("stop",       (this->stop),       cstop);
   
   return waveTree;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++
// dump livetime into tree
//++++++++++++++++++++++++++++++++++++++++++++++++++++
void livetime::output(TTree* waveTree, network* net, float* slag, vector<waveSegment> detSegs)
{
  size_t n,m;
  size_t M = net->ifoListSize();
  if(!M) return;

  // add detectors to tree user info
//  if(waveTree!=NULL) for(m=0;m<M;m++) waveTree->GetUserInfo()->Add(net->getifo(m));

//Fill tree
  for(m=0; m<(NIFO_MAX+1); m++) {
    this->lag[m]  = -1;
    this->slag[m] = -1;
  }                     
  for(m=0; m<NIFO_MAX; m++) {
    this->start[m] = 0;
    this->stop[m]  = 0;
  }                     
  for(m=0; m<detSegs.size(); m++) {
    this->start[m] = detSegs[m].start;
    this->stop[m]  = detSegs[m].stop;
  }
  this->run = net->nRun;
  this->gps = net->getifo(0)->getTFmap()->start();

  if(slag!=NULL) for(m=0; m<M+1; m++) this->slag[m] = slag[m];

  for(n=0; n<net->nLag; n++){               // loop on time lags

    this->live = net->getliveTime(n);

    for(m=0; m<NIFO_MAX; m++){                     // loop on detectors
      if(m<M) {
	this->lag[m] = net->getifo(m)->lagShift.data[n];
      }
    }
    this->lag[M] = net->getwc(n)->shift;
    waveTree->Fill();
  }
}



Bool_t livetime::Notify()
{
   // Called when loading a new file.
   // Get branch pointers.
   b_run  = fChain->GetBranch("run");
   b_gps  = fChain->GetBranch("gps");
   b_live = fChain->GetBranch("live");
   b_lag  = fChain->GetBranch("lag");
   b_slag = fChain->GetBranch("slag");
   b_start= fChain->GetBranch("start");
   b_stop = fChain->GetBranch("stop");
   return kTRUE;
}

Int_t livetime::GetEntry(Int_t entry) 
{ 
  if (!fChain) return 0; 
  return fChain->GetEntry(entry); 
};

void livetime::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}


/*
Int_t livetime::LoadTree(Int_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

Int_t livetime::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void livetime::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L livetime.C
//      Root > livetime t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show rmss of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//    This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Int_t nentries = Int_t(fChain->GetEntriesFast());

   Int_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
*/
