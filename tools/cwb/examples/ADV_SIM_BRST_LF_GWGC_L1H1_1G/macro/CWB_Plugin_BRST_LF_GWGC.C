#define XIFO 5

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

// Plugin to generate simulated gaussian noise and injected 'on the fly' CBC MDC

  cout << endl;
  cout << "-----> CWB_Plugin_BRST_LF_GWGC.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->dataPlugin=true; // disable read data from frames
    cfg->mdcPlugin=true;  // disable read mdc from frames
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
  }

  if(type==CWB_PLUGIN_MDC) {  

    char cmd[128];
    sprintf(cmd,"network* net = (network*)%p;",net);
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

    // ---------------------------------
    // write MDC log file
    // ---------------------------------

    char logFile[512];
    int runID = net->nRun;
    int Tb=x->start()+cfg->segEdge;
    int dT=x->size()/x->rate()-2*cfg->segEdge;
    sprintf(logFile,"%s/log_%d_%d_%s_job%d.txt",cfg->output_dir,int(Tb),int(dT),cfg->data_label,runID);
    cout << "Dump : " << logFile << endl;
    if(ifo==cfg->ifo[0]) MDC.DumpLog(logFile);

    // ---------------------------------
    // print MDC injections list
    // ---------------------------------

    cout.precision(14);
    for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;
  }

  return;
}
