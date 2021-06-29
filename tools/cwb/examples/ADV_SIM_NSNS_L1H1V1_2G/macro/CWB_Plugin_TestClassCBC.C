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

// Plugin to generate simulated gaussian noise and injected 'on the fly' Burst MDC

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_TestClassMDC.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
//    cfg->dataPlugin=true; // disable read data from frames
    cfg->mdcPlugin=true;  // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_DATA) {  
/*
    CWB::Toolbox TB;

    int seed;
    if(ifo.CompareTo("L1")==0) seed=1000;
    if(ifo.CompareTo("H1")==0) seed=2000;
    if(ifo.CompareTo("V1")==0) seed=3000;
    if(ifo.CompareTo("J1")==0) seed=4000;
    if(ifo.CompareTo("A2")==0) seed=5000;
    if(ifo.CompareTo("Y2")==0) seed=6000;
    if(ifo.CompareTo("Y3")==0) seed=7000;

    TString fName;
    if(ifo.CompareTo("L1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("V1")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";
    if(ifo.CompareTo("J1")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";
    if(ifo.CompareTo("A2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("Y2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("Y3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";

    int size=x->size();
    double start=x->start();
    TB.getSimNoise(*x, fName, seed, net->nRun);
    x->resize(size);
    x->start(start);
*/
/*
    // dump spectrum
    char file[1024];
    sprintf(file,"%s/sensitivity_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 8, cfg->segEdge);
    if(TString(ifo).CompareTo("V1")==0) exit(0);
*/
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

    // Dump MDC
/*
    char mdcFile[1024];
    double Tb = x->start();
    double dT = x->size()/x->rate();
    sprintf(mdcFile,"%s/MDC_%s_%d_%d_%s_job%d.root",
            cfg->output_dir,ifo.Data(),int(Tb),int(dT),cfg->data_label,net->nRun);
    TFile* ifile = new TFile(mdcFile, "RECREATE");
    if(ifile==NULL||!ifile->IsOpen())
      {cout << "cwb::cwb - Error opening root file : " << mdcFile << endl;exit(1);}
    cout << "Write MDC data : " << mdcFile << endl;
    x->Write(ifo.Data());   
    ifile->Close();
    //if(TString(ifo).CompareTo("V1")==0) exit(0);
*/
  }

  if(type==CWB_PLUGIN_WHITE) {  
/*
    CWB::Toolbox TB;
    //int level=x->getLevel(); 
    //x->Inverse(-1);
    // dump spectrum
    char file[1024];
    sprintf(file,"%s/sensitivity_white_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 8, cfg->segEdge);
    if(TString(ifo).CompareTo("V1")==0) exit(0);
    //x->Forward(level);
*/
/*
    bool isWavearray=false;
    if(TString(cfg->analysis)=="1G") isWavearray=false;
    if(TString(cfg->analysis)=="2G") isWavearray=true;

    WSeries<double>* w = NULL;
    if(isWavearray) { // x is a wavearray
      // input x is a wavearray and must be converted in a WSeries<double> type and apply Forward
      cout << "Input data is a wavearray" << endl;
      Meyer<double> S(1024,2);
      w = new WSeries<double>();
      w->Forward(*x,S,cfg->levelD);
    } else { // x is wseries
      cout << "Input data is a wseries" << endl;
      w = x;
    }

    cout << "Decomposition Level : " << w->getLevel() << endl;

    watplot WTS(const_cast<char*>("wtswrc"));
    //scalogram maps
    double start = w->start()+cfg->segEdge;
    double stop  = w->start()+w->size()/w->rate()-cfg->segEdge;
    double flow  = 64;
    double fhigh = 2048;
    WTS.plot(w, 2, start, stop,const_cast<char*>("COLZ"));
    WTS.hist2D->GetYaxis()->SetRangeUser(flow, fhigh);
    // dump spectrum
    char fname[1024];
    sprintf(fname,"%s/scalogram_white_%s_%d_%s_lev%d_job%lu.root",cfg->dump_dir,ifo.Data(),
            int(w->start()),cfg->data_label,w->getLevel(),net->nRun);
    cout << endl << "Dump Scalogram : " << fname << endl << endl;
    WTS.canvas->Print(fname);
    int nIFO=net->ifoListSize();
    if(isWavearray) delete w;
    if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) exit(0);  // last ifo
*/
  }

  if(type==CWB_PLUGIN_LIKELIHOOD) {  
/*
    cout << "CWB_PLUGIN_LIKELIHOOD :  store net into net file" << endl;
    // store net into net file
    TString ename = jfile->GetPath();
    ename.ReplaceAll(":/","");
    ename.ReplaceAll("job","net");
    ename.ReplaceAll(TString(cfg->nodedir)+"/",TString(cfg->output_dir)+"/");
    cout << ename.Data() << endl;
    TFile* efile = new TFile(ename,"CREATE");
    if(TString(ifo).CompareTo(net->ifoName[0])==0) {
      efile = new TFile(ename,"CREATE");
    } else {
      efile = new TFile(ename,"UPDATE");
    }
    if(efile==NULL||!efile->IsOpen())
      {cout << "CWB_Plugin_TestClassCBC - Error : file " << ename << " not found" <<  endl;exit(1);}
    cfg->Write("config");
    net->Write("network");
    efile->ls();
    efile->Write("",TObject::kOverwrite);
    efile->Close();
*/
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
