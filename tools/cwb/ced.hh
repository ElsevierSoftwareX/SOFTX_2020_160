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


/**********************************************************
 * Package:      ced Class Library
 * File name:    ced.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef CED_HH
#define CED_HH

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "watfun.hh"
#include "injection.hh"
#include "detector.hh"
#include "netcluster.hh"
#include "network.hh"
#include "netevent.hh"


using namespace std;

namespace CWB {

class ced {

public:

  // rbasedirCED : plots are stored in the sbasedirCED directory 
  ced(network* NET, netevent* EVT, char* sbasedirCED) {
    if(NET==NULL) {cout << "ced::ced Error : NET is NULL" << endl;exit(1);}
    if(EVT==NULL) {cout << "ced::CED Error : EVT is NULL" << endl;exit(1);}
    if(sbasedirCED==NULL) {cout << "ced::ced Error : sbasedirCED is NULL" << endl;exit(1);}
    this->NET = NET;
    this->EVT = EVT;
    this->sbasedirCED = sbasedirCED;
    this->rbasedirCED = NULL;
    this->simulation=0;
    this->rho=0;
    this->inRate=0;
    this->spectrogram_zmax=0;
    strcpy(this->gtype,"png");
    this->paletteId=0;
    strcpy(this->chName,"");
    this->useSparse=false;
  };

  // rbasedirCED : plots are stored in root file (opened externally) 
  ced(network* NET, netevent* EVT, TDirectory* rbasedirCED) {
    if(NET==NULL) {cout << "ced::ced Error : NET is NULL" << endl;exit(1);}
    if(EVT==NULL) {cout << "ced::ced Error : EVT is NULL" << endl;exit(1);}
    if(rbasedirCED==NULL) {cout << "ced::ced Error : rbasedirCED is NULL" << endl;exit(1);}
    this->NET = NET;
    this->EVT = EVT;
    this->sbasedirCED = NULL;
    this->rbasedirCED = rbasedirCED;
    this->simulation=0;
    this->rho=0;
    this->inRate=0;
    this->spectrogram_zmax=0;
    strcpy(this->gtype,"png");
    this->paletteId=0;
    strcpy(this->chName,"");
    this->useSparse=false;
  };

  ~ced() {};

  void SetOptions(int simulation, double rho, double inRate, bool useSparse=false, 
                  char* gtype=const_cast<char*>("png"), int paletteId=0) {
    this->simulation=simulation; // 1/2 for simulation(hrss/snr), 0 for production
    this->rho=rho;	         // ced is produced for events with rho >= this->rho threshold
    this->inRate=inRate;	 // original data rate, used to rescale data when data are resampled
    if(inRate<=0) {cout << "ced::SetOptions : Error - inrate must be > 0" << endl;exit(1);}
    strcpy(this->gtype,gtype);
    this->paletteId=paletteId;
    this->useSparse=useSparse;   // if true sparse map are used instead of TF
  }

  void SetOptions(TString options);

  void SetChannelName(char* chName) {strcpy(this->chName,chName);}
  int Write(double factor, bool fullCED=true); 

private:

  void Write(double factor, size_t iID, int LAG, char* dirCED);

  // x     : input array
  // P     : energy percentage
  // bT,eT : begin,end time of the range which contains percentage P of the total energy E
  double GetBoundaries(wavearray<double> x, double P, double& bT, double& eT);

  network*  NET;		//!
  netevent* EVT;		//!

  TDirectory* rbasedirCED;	//!
  char*  sbasedirCED;		//!
  int    simulation;
  double rho;
  double inRate;	
  char   gtype[32];
  int    paletteId; 
  char   chName[64];         
  bool   useSparse;
  double spectrogram_zmax;

  // used by THtml doc
  ClassDef(ced,2)
};

} // end namespace

#endif

