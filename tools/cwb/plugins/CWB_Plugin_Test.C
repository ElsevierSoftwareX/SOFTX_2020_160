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


#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
/*
#define IMPORT_TST(TYPE,VAR,SIZE1,SIZE2) { 					\
  TGlobal* global = gROOT->GetGlobal(#VAR,true); 				\
  if(global!=NULL) memcpy(&VAR,global->GetAddress(),SIZE1*SIZE2*sizeof(TYPE)); }
*/
void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// Test

  cout << "CWB PLUGIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
 
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;

/*
  int check=0;
  int check2=0;
  cfg->configPlugin.Exec(); 
  //TGlobal* global = (TGlobal*)gROOT->GetGlobal("check",true);
  //check = *(int*)global->GetAddress();
  //check = *(int*)((TGlobal*)gROOT->GetGlobal("check",true))->GetAddress();
  //check = *(int*)gROOT->GetGlobal("check",true)->GetAddress();
  IMPORT(int,check,1,1);
  //memcpy((void*)&check,(void*)&check2,sizeof(check));

//  TGlobal* global = gROOT->GetGlobal("check",true);
//  if(global!=NULL) memcpy((void*)&check,global->GetAddress(),sizeof(check)); 

  cout << "check " << check << endl;

  gSystem->Exit(0);
*/
  return;
}
