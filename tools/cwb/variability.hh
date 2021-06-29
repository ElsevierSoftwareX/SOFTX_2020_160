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


//////////////////////////////////////////////////////////
//   noise variability class 
//   Sergey Klimenko, University of Florida
//////////////////////////////////////////////////////////


#ifndef variability_h
#define variability_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "wavearray.hh"

/* Structure of WaveBurst single IFO event */

class variability {
   public :

  TTree          *fChain;        //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent;      //!current Tree number in a TChain

//Declaration of leaves types
// for arrays: ifo1 - first index, ifo2 - second index

  Int_t           nevent;            // event count                                                  
  Int_t           ifo;               // ifo ID: 1/2/3 - L1/H1/H2                                     
  Float_t         value;             // noise variability                                   
  Double_t        time;              // average center_of_snr time                                   
  Double_t        gps;               // segment gps time                                   

//List of branches

   TBranch        *b_nevent;   //!
   TBranch        *b_ifo;   //!
   TBranch        *b_value;   //!
   TBranch        *b_time;   //!
   TBranch        *b_gps;   //!

   variability() {fChain=NULL;};
   variability(TTree *tree) {if(tree) Init(tree);};
   virtual ~variability() { if (!fChain) return; delete fChain->GetCurrentFile(); };

   virtual variability& operator<<(variability &);

   Int_t  GetEntry(Int_t);
   void   Init(TTree *);
   Bool_t Notify();
   TTree* setTree();
   void output(TTree*, wavearray<float>*, int=0, double=0.);

//   void   Loop();
//   Int_t  Cut(Int_t entry);
//   Int_t  LoadTree(Int_t entry);

   void   Show(Int_t entry = -1);

   // used by THtml doc
   ClassDef(variability,1)	
};

#endif




