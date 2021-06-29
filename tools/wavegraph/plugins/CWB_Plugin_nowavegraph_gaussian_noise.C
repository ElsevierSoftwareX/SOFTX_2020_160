/*
# Copyright (C) 2019 Eric Chassande-Mottin, Philippe Bacon, Gayathri V, Archana Pai
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
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TMath.h"
#include "TPolyMarker.h"
#include "mdc.hh"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iterator>

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// This plugin implements the wavegraph in coherence stage (only 2G)

  cout << endl;
  cout << "-----> CWB_Plugin_nowavegraph.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->dataPlugin=true;       // disable read data from frames
    cfg->mdcPlugin=true;        // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_DATA) {

    CWB::Toolbox TB;

    int seed;
    if(ifo.CompareTo("L1")==0) seed=1001;
    if(ifo.CompareTo("L2")==0) seed=1002;
    if(ifo.CompareTo("L3")==0) seed=1003;

    if(ifo.CompareTo("H1")==0) seed=2001;
    if(ifo.CompareTo("H2")==0) seed=2002;
    if(ifo.CompareTo("H3")==0) seed=2003;

    if(ifo.CompareTo("V1")==0) seed=3001;
    if(ifo.CompareTo("V2")==0) seed=3002;
    if(ifo.CompareTo("V3")==0) seed=3003;

    if(ifo.CompareTo("I1")==0) seed=4001;
    if(ifo.CompareTo("I2")==0) seed=4002;
    if(ifo.CompareTo("I3")==0) seed=4003;

    if(ifo.CompareTo("J1")==0) seed=5001;
    if(ifo.CompareTo("J2")==0) seed=5002;
    if(ifo.CompareTo("J3")==0) seed=5003;

    TString fName;
    if(ifo.CompareTo("L1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("L2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("L3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";

    if(ifo.CompareTo("H1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";

    if(ifo.CompareTo("V1")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";
    if(ifo.CompareTo("V2")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";
    if(ifo.CompareTo("V3")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";

    if(ifo.CompareTo("I1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("I2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("I3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";

    if(ifo.CompareTo("J1")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";
    if(ifo.CompareTo("J2")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";
    if(ifo.CompareTo("J3")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";

    int size=x->size();
    double start=x->start();
    TB.getSimNoise(*x, fName, seed, net->nRun);
    x->resize(size);
    x->start(start);

    /* 
    // dump spectrum                
    char file[1024]; 
    sprintf(file,"%s/sensitivity_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun); 
    cout << endl << "Dump Sensitivity : " << file << endl << endl;            
    TB.makeSpectrum(file, *x);  
    if(TString(ifo).CompareTo("V1")==0) exit(0);
    */
  }

  if(type==CWB_PLUGIN_MDC) {

    char cmd[128];
    sprintf(cmd,"network* net = (network*)%p;",net);
    gROOT->ProcessLine(cmd);
    sprintf(cmd,"CWB::config* cfg = (CWB::config*)%p;",cfg);
    gROOT->ProcessLine(cmd);

    CWB::mdc MDC(net);

    // ---------------------------------
    // read plugin config
    // ---------------------------------

    cfg->configPlugin.Exec();

    // ---------------------------------
    // set list of mdc waveforms
    // ---------------------------------

    IMPORT(CWB::mdc,MDC)
    MDC.Print();

    // ---------------------------------
    // get mdc data
    // ---------------------------------

    MDC.Get(*x,ifo);

    // ---------------------------------
    // set mdc list in the network class
    // ---------------------------------

    if(ifo.CompareTo(net->ifoName[0])==0) {
      net->mdcList.clear();
      net->mdcType.clear();
      net->mdcTime.clear();
      net->mdcList=MDC.mdcList;
      net->mdcType=MDC.mdcType;
      net->mdcTime=MDC.mdcTime;
    }

    cout.precision(14);
    for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;

  }

}
