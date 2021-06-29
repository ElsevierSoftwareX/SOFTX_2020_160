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


//!NOISE_MDC_SIMULATION
// Plugin to injected MDC 'on the fly'  

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
#include "TComplex.h"
#include "TMath.h"
#include "mdc.hh"
#include <vector>

//#define DUMP_LOG	 // uncomment to enable dump log

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->mdcPlugin=true;         // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_MDC) {  

    cout << "Execute CWB_Plugin_MDC_OTF.C : Inject On The Fly MDC ..." << endl;

    // ---------------------------------
    // Declare mdc class 
    // On The Fly MDC Injections
    // ---------------------------------

    CWB::mdc MDC(net);
    CWB_PLUGIN_EXPORT(MDC)

    // export global variables to the config plugin
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    int xstart = (int)x->start();
    int xstop  = (int)x->stop();
    CWB_PLUGIN_EXPORT(net)
    CWB_PLUGIN_EXPORT(cfg)
    CWB_PLUGIN_EXPORT(xstart)
    CWB_PLUGIN_EXPORT(xstop)
    CWB_PLUGIN_EXPORT(gIFACTOR)
    int gIFOID=0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {gIFOID=n;break;}
    CWB_PLUGIN_EXPORT(gIFOID)

    // ---------------------------------
    // read plugin config 
    // CWB_Plugin_MDC_OTF_Config.C
    // ---------------------------------

    int error=0;
    // execute config plugin
    cfg->configPlugin.Exec(NULL,&error);
    if(error) {
      cout << "Error executing macro : " << cfg->configPlugin.GetTitle() << endl;
      cfg->configPlugin.Print(); 
      gSystem->Exit(1); 
    }
    // print list waveforms declared in the config plugin 
    MDC.Print();

    // ---------------------------------
    // get mdc data
    // fill x array with MDC injections
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

    // ---------------------------------
    // write MDC log file  
    // if enabled (uncomment #define DUMP_LOG) then the txt log file is created under the output dir
    // the cwb_merge collect all log job files on a single file under the merge directory
    // ---------------------------------

#ifdef DUMP_LOG
    char logFile[512];
    int runID = net->nRun;
    int Tb=x->start()+cfg->segEdge;
    int dT=x->size()/x->rate()-2*cfg->segEdge;
    sprintf(logFile,"%s/log_%d_%d_%s_job%d.txt",cfg->output_dir,int(Tb),int(dT),cfg->data_label,runID);
    cout << "Dump : " << logFile << endl;
    if(ifo==cfg->ifo[0]) MDC.DumpLog(logFile);
#endif

    // ---------------------------------
    // print MDC injections list 
    // ---------------------------------

    cout.precision(14);
    if(ifo.CompareTo(net->ifoName[0])==0) {
      for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
      for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
      for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;
    }
  }

  return;
}
