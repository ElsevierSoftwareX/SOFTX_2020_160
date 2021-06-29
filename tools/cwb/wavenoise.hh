/*
# Copyright (C) 2019 Gabriele Vedovato
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


//////////////////////////////////////////////////////////
//   noise variance class 
//   Sergey Klimenko, University of Florida
//////////////////////////////////////////////////////////


#ifndef wavenoise_h
#define wavenoise_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "wseries.hh"

/* Structure of WaveBurst single IFO event */

class wavenoise {
   public :

  TTree          *fChain;        //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent;      //!current Tree number in a TChain

//Declaration of leaves types
// for arrays: ifo1 - first index, ifo2 - second index

  Int_t           layer;           // event count                                                  
  Int_t           ifo;             // ifo ID: 1/2/3 - L1/H1/H2                                     
  Double_t        rms;             // noise standard deviation                                   
  Float_t         frequency;       // frequency band                                   
  Double_t        time;            // measurement gps time                                   
  Double_t        gps;             // segment gps time                                   

//List of branches

   TBranch        *b_layer;   //!
   TBranch        *b_ifo;   //!
   TBranch        *b_rms;   //!
   TBranch        *b_frequency;   //!
   TBranch        *b_time;   //!
   TBranch        *b_gps;   //!

   wavenoise() {fChain=NULL;};
   wavenoise(TTree *tree) {if(tree) Init(tree);};
   virtual ~wavenoise() { if (!fChain) return; delete fChain->GetCurrentFile(); };

   virtual wavenoise& operator<<(wavenoise &);

   Int_t  GetEntry(Int_t);
   void   Init(TTree *);
   Bool_t Notify();
   TTree* setTree();
   void output(TTree*, WSeries<double>*, int=0, double=1.);

//   void   Loop();
//   Int_t  Cut(Int_t entry);
//   Int_t  LoadTree(Int_t entry);

   void   Show(Int_t entry = -1);

   // used by THtml doc
   ClassDef(wavenoise,1)	
};

#endif




