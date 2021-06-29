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

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// Plugin to generate simulated gaussian noise and injected 'on the fly' EBBH

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_EBBH.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
//    cfg->dataPlugin=true; // disable read data from frames
    cfg->mdcPlugin=true;  // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_MDC || type==CWB_PLUGIN_INIT_JOB) {  

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

    if(type==CWB_PLUGIN_INIT_JOB) {  
      if(ifo.CompareTo(net->ifoName[0])==0) {
        net->mdcList.clear();
        net->mdcType.clear();
        net->mdcTime.clear();
        net->mdcList=MDC.mdcList;
        net->mdcType=MDC.mdcType;
        net->mdcTime=MDC.mdcTime;
        double Tb = x->start();
        double dT = x->size()/x->rate();
        char log_label[512];
        char tmpFile[1024];
        sprintf(log_label,"%d_%d_%s_job%d",int(Tb),int(dT),cfg->data_label,net->nRun);
        sprintf(tmpFile,"%s/%s-LogTMP.txt",cfg->log_dir,log_label);
        cout << "Write MDC Log : " << tmpFile << endl;
        MDC.DumpLog(tmpFile);  
      }
    }

    cout.precision(14);
    for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;
  }

  if(type==CWB_PLUGIN_CLOSE_JOB) {  
    if(ifo.CompareTo(net->ifoName[0])==0) {
      double Tb = x->start();
      double dT = x->size()/x->rate();
      char log_label[512];
      char tmpFile[1024];
      char outFile[1024];
      sprintf(log_label,"%d_%d_%s_job%d",int(Tb),int(dT),cfg->data_label,net->nRun);
      sprintf(tmpFile,"%s/%s-LogTMP.txt",cfg->log_dir,log_label);
      sprintf(outFile,"%s/%s-Log.txt",cfg->log_dir,log_label);
      char command[1024];
      sprintf(command,"/bin/mv %s %s", tmpFile, outFile);
      gSystem->Exec(command);
    }
  }

  return;
}
