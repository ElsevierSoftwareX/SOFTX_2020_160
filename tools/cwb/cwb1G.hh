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
 * Package:      cwb1G Class Library
 * File name:    cwb1G.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef CWB1G_HH
#define CWB1G_HH

#include "cwb.hh"

class cwb1G : public cwb {

public:

  // Contructor
  cwb1G(CWB_STAGE jstage=CWB_STAGE_FULL) : cwb(jstage) {}

  // Contructor
  cwb1G(TString fName,TString xName="", CWB_STAGE jstage=CWB_STAGE_FULL) : cwb(fName,xName) {Init();}

  // Contructor
  cwb1G(CWB::config cfg, CWB_STAGE jstage=CWB_STAGE_FULL) : cwb(cfg) {Init();}

  virtual ~cwb1G();

protected:

  double ReadData(double mdcShift, int ifactor);
  void   DataConditioning(int ifactor);
  void   Coherence(int ifactor);
  void   SuperCluster(int ifactor);
  bool   Likelihood(int ifactor, char* ced_dir, netevent* output = NULL,
                    TTree* net_tree = NULL, char* outDump = NULL);

  void   Init();

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// declarations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ClassDef(cwb1G,1)
};  

#endif
