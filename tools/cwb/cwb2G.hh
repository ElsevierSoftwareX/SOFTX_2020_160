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
 * Package:      cwb2G Class Library
 * File name:    cwb2G.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef CWB2G_HH
#define CWB2G_HH

#include "cwb.hh"
#include "WDM.hh"


class cwb2G : public cwb {

public:
 
  // Constructor 
  cwb2G(CWB_STAGE jstage=CWB_STAGE_FULL) : cwb(jstage) {
    for(int n=0;n<NRES_MAX;n++) pwdm[n]=NULL; 
    for(int n=0;n<NIFO_MAX;n++)  hot[n]=NULL; 
  }
 
  // Constructor 
  cwb2G(TString fName, TString xName="", CWB_STAGE jstage=CWB_STAGE_FULL) : cwb(fName,xName,jstage) {Init();}
 
  // Constructor 
  cwb2G(CWB::config cfg, CWB_STAGE jstage=CWB_STAGE_FULL) : cwb(cfg,jstage) {Init();}

  ~cwb2G();  

//  void run(int runID=0);

  void   WriteSparseTFmap(TFile* jfile, int ifactor, TString tdir, TString tname);
  void   FillSparseTFmap(TFile* jfile, int ifactor, TString tname);

  CWB::ced* GetCED() {return ced;}	// return pointer to CED object

//private:	// ROOT6 fix, commented out otherwise not visibles from cWB plugins

  double ReadData(double mdcShift, int ifactor);
  void   DataConditioning(int ifactor);
  void   DataConditioning(TString fName, int ifactor);
  void   Coherence(int ifactor);
  void   SuperCluster(int ifactor);
  bool   Likelihood(int ifactor, char* ced_dir, netevent* output = NULL, 
                    TTree* net_tree = NULL, char* outDump = NULL);

  void   Init();

  void   SaveWaveforms(TFile* jfile, detector* pD, int ifactor, int wfSAVE=1);
  void   LoadWaveforms(TFile* ifile, detector* pD, int ifactor, int wfSAVE=1);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// declarations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  int nRES;     		     // number of frequency resolution levels
  double TDRate;   		     // time-delay filter rate

  WDM<double>* pwdm[NRES_MAX];       //! wavelet pointers: pwdm[0] - l_high, wdm[nRES-1] l_low
  wavearray<double>* hot[NIFO_MAX];  //! temporary time series

  CWB::ced *ced;                     //! CED pointer

  ClassDef(cwb2G,3)
};  

#endif
