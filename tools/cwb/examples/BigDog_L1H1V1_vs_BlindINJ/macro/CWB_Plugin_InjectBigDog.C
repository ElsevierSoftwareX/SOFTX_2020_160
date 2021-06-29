//!NOISE_MDC_SIMULATION
// Plugin to injected 'on the fly' BigDog MDC

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
#include "frame.hh"
#include <vector>

//#define DUMP_LOG	 // uncomment to enable dump log

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->mdcPlugin=true;         // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_DATA) {

    // correct wrong HW INJ sign for L1,H1
    if(ifo=="L1") (*x)*=-1.;
    if(ifo=="H1") (*x)*=-1.;
  }

  if(type==CWB_PLUGIN_MDC) {  

    cout << "Execute CWB_Plugin_MDC_OTF.C : Inject On The Fly MDC ..." << endl;

    // ---------------------------------
    // Declare mdc class 
    // On The Fly MDC Injections
    // ---------------------------------

    CWB::mdc MDC(net);
    CWB_PLUGIN_EXPORT(MDC)

    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    int xstart = (int)x->start();
    int xstop  = (int)x->stop();
    CWB_PLUGIN_EXPORT(net)
    CWB_PLUGIN_EXPORT(cfg)
    CWB_PLUGIN_EXPORT(xstart)
    CWB_PLUGIN_EXPORT(xstop)
    CWB_PLUGIN_EXPORT(gIFACTOR)

    // ---------------------------------
    // read plugin config 
    // CWB_Plugin_MDC_OTF_Config.C
    // ---------------------------------

#ifndef _USE_ROOT6
    char cmd[128]; 
    // export gIFOID to the config plugin
    int gIFOID=0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {gIFOID=n;break;}
    sprintf(cmd,"int gIFOID = %d;",gIFOID);
    gROOT->ProcessLine(cmd);
    // export network object to the config plugin
    sprintf(cmd,"network* net = (network*)%p;",net);
    gROOT->ProcessLine(cmd);
    // export config object to the config plugin
    sprintf(cmd,"CWB::config* cfg = (CWB::config*)%p;",cfg);
    gROOT->ProcessLine(cmd);
    // export job start time to the config plugin
    sprintf(cmd,"int xstart = %d;",(int)x->start());
    gROOT->ProcessLine(cmd);
    // export job stop time to the config plugin
    sprintf(cmd,"int xstop = %d;",(int)x->stop());
    gROOT->ProcessLine(cmd);
#endif
    int error=0;
    // execute config plugin
    cfg->configPlugin.Exec(NULL,&error);
    if(error) {
      cout << "Error executing macro : " << cfg->configPlugin.GetTitle() << endl;
      cfg->configPlugin.Print(); 
      gSystem->Exit(1); 
    }
#ifndef _USE_ROOT6
    // import mdc object initialized in the config plugin
    IMPORT(CWB::mdc,MDC) 
#endif
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

  if(type==CWB_PLUGIN_STRAIN_AND_MDC) {

    // get ifo index
    int xIFO =0;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==cfg->ifo[n]) {xIFO=n;break;}

    int runID = net->nRun;
    char flabel[512];
    int Tb=x->start();
    int dT=x->size()/x->rate();
    sprintf(flabel,"%d_%d_%s_job%d",int(Tb),int(dT),cfg->data_label,runID);

    int size=0;

    // in channel
    cout << "x->rate() : " << x->rate() << endl;
    wavearray<double> xx(x->size());
    xx.rate(x->rate());
    xx.start(x->start());
    xx.stop(x->stop());
    WSeries<double> w;

    // read strain again and write to x array
    CWB::frame frl(cfg->frFiles[xIFO],"READ");
    frl.setChName(cfg->channelNamesRaw[xIFO]);
    frl >> xx;
    cout << "xx.rate() : " << xx.rate() << endl;
    Meyer<double> B(1024);           // set wavelet for resampling
    w.Forward(xx,B,cfg->levelR);
    cout << "w.rate() : " << w.rate() << endl;
    w.getLayer(xx,0);
    if(TString(cfg->analysis)=="1G") {
      Meyer<double> S(1024,2);
      x->Forward(xx,S,cfg->levelD);
    }
    if(TString(cfg->analysis)=="2G") {
      for(int i=0;i<(int)x->size();i++) x->data[i] = xx[i];
    }
    cout << "x->rate() : " << x->rate() << endl;

    // correct wron HW INJ sign for L1,H1
    if(ifo=="L1") (*x)*=-1.;
    if(ifo=="H1") (*x)*=-1.;
  }

  return;
}
