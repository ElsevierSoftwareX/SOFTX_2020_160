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
#include "cwb2G.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "gwavearray.hh"
#include "TString.h"    
#include "TObjArray.h"  
#include "TObjString.h" 
#include "TRandom.h"    
#include "TComplex.h"   
#include "TMath.h"      
#include "mdc.hh"       
#include <vector>       

//!DATA_CONDITIONING

// This plugin is used to save the nRMS wseries into the job root file
// nRMS contains the estimated noise rms used to whitening the strain data
// Must be used with lagSize=1, slagSize=1, simulation=0 and input stage = CSTRAIN
//
// To use in interactive mode do :
//   cwb_inet2G config/user_parameters.C CSTRAIN #job 
// the cstrain job file is saved under the data directory 
//
// To use in batch mode do :
//   cwb_condor create CSTRAIN
//   cwb_condor submit
// the cstrain job files are saved under the output directory 
//

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  if(type==CWB_PLUGIN_CONFIG) {

    // import istage,jstage
    int gISTAGE=-1; IMPORT(int,gISTAGE)
    int gJSTAGE=-1; IMPORT(int,gJSTAGE)

    if(gJSTAGE!=CWB_STAGE_CSTRAIN) 
      {cout << "CWB_Plugin_nRMS - Error : input stage must be CSTRAIN" << endl;gSystem->Exit(1);}

    if(cfg->simulation!=0)
      {cout << "CWB_Plugin_nRMS - Error : simulation must be 0" << endl;gSystem->Exit(1);}

    if(cfg->lagSize!=1)
      {cout << "CWB_Plugin_nRMS - Error : lagSize must be 1" << endl;gSystem->Exit(1);}

    if(cfg->slagSize!=1)
      {cout << "CWB_Plugin_nRMS - Error : slagSize must be 1" << endl;gSystem->Exit(1);}
  }

  if(type==CWB_PLUGIN_OREADDATA) {

    // import pointer to object cwb2G
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)

    // disable the builtin save cstrain & rms
    gCWB2G->jobfOptions&=(CWB_JOBF_SAVE_ALL-CWB_JOBF_SAVE_CSTRAIN);	
  }

  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {

    // get ifo ID
    int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {ifoID=n;break;}
    // get detector pointer
    detector* pD = net->getifo(ifoID);  
    // import factor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR) 

    // create the rms directory in root file
    jfile->cd();
    TDirectory* drms = (TDirectory*)jfile->Get("rms");
    if(drms==NULL) drms=jfile->mkdir("rms");
    drms->cd(); 
    char cdrms_name[32];sprintf(cdrms_name,"rms-f%d",gIFACTOR);
    TDirectory* cdrms = (TDirectory*)drms->Get(cdrms_name);
    if(cdrms==NULL) cdrms=drms->mkdir(cdrms_name);

    // save nRMS
    cdrms->cd();pD->nRMS.Write(pD->Name);
  }

}
