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
//   root tree class to store live time   
//   Sergey Klimenko, University of Florida
//////////////////////////////////////////////////////////


#ifndef livetime_h
#define livetime_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <stdlib.h>
#include "wat.hh"
#include "network.hh"

/* Structure of WaveBurst single IFO event */

static vector<waveSegment> DEFAULT_WAVESEGMENT;

class livetime {
   public :

  TTree          *fChain;        //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent;      //!current Tree number in a TChain

//Declaration of leaves types
// for arrays: ifo1 - first index, ifo2 - second index

  Int_t           run;            // run                                                  
  Double_t        gps;            // gps start time segment of the first detector
  Double_t        live;           // live time                                   
  Float_t*        lag;            //! time lag [sec]                                                   
  Float_t*        slag;           //! time slag [sec]                                                   
  Double_t*       start;          //! array of gps start time segment for each detector                                  
  Double_t*       stop;           //! array of gps stop time segment for each detector        

//List of branches

   TBranch        *b_run;   	//!
   TBranch        *b_gps;   	//!
   TBranch        *b_live;  	//!
   TBranch        *b_lag;   	//!
   TBranch        *b_slag;   	//!
   TBranch        *b_start;   	//!
   TBranch        *b_stop;   	//!

   livetime() : lag(NULL), slag(NULL), start(NULL), stop(NULL) {allocate(); fChain=NULL;};
   livetime(TTree *tree) : lag(NULL), slag(NULL), start(NULL), stop(NULL) {allocate(); if(tree) Init(tree);};
   virtual ~livetime() { 
     if (lag)         free(lag);         // lag array                            
     if (slag)        free(slag);        // slag array                            
     if (start)       free(start);       // start array                            
     if (stop)        free(stop);        // stop array                            
     if (fChain)      delete fChain->GetCurrentFile();
     return;
   };

   virtual livetime& operator=(livetime &);

   Int_t  GetEntry(Int_t);
   void   Init(TTree *);
   void   allocate();
   Bool_t Notify();
   TTree* setTree();
   void output(TTree* waveTree, network* net, float* slag=NULL, 
               vector<waveSegment> detSegs = DEFAULT_WAVESEGMENT);

//   void   Loop();
//   Int_t  Cut(Int_t entry);
//   Int_t  LoadTree(Int_t entry);

   void   Show(Int_t entry = -1);

   // used by THtml doc
   ClassDef(livetime,2)		
};

#endif




