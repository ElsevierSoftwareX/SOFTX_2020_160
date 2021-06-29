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
#include "TSystem.h"
#include "mdc.hh"
#include "WDM.hh"
#include "regression.hh"
#include "frame.hh"
#include "watplot.hh"
#include "Biorthogonal.hh"
#include <fstream>
#include <vector>


using namespace CWB;

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  cout << endl;
  cout << "-----> macro/CWB_Plugin_Vaibhav_Regression.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << endl;


  //---------------------------------------------------------------------
  // MAIN
  //---------------------------------------------------------------------

  bool APPLY_REGRESSION                 = true;
  bool SAVE_FRAME                       = true;
  bool SAVE_PLOT                        = true;
  bool SAVE_PLOT_INDATA                 = false;
  bool CUT_0_32_HZ                      = true;

  TString OPLOT_DIR                     = "oplots";
  TString OFRAME_DIR                    = "oframes";

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS FOR H1
  //---------------------------------------------------------------------

  TString H1_FRLIST_WITNESS             = "input/H1_RDS_R_L1_94240000_942500000.frl";
  TString H1_WITNESS_CHANNEL_LIST       = "input/WitnessChannels_H1.lst";

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS
  //---------------------------------------------------------------------

  int    nFILTER                        = 5;
  double LAYER_WIDTH                    = 1.0;
  int    nPL                            = 5;       // number of processed power lines
  double FWIDTH                         = 5;       // frequency width used by regression

  // parameters to clean power lines

  double APPLY_THRESHOLD                = 0.2;
  double SOLVE_THRESHOLD                = 0.0;
  double SOLVE_NEIGEN_PER_LAYER         = 0;
  char   SOLVE_REGULATOR                = 'h';

  // ---------------------------------
  // read plugin config
  // ---------------------------------

  cfg->configPlugin.Exec();

  IMPORT(bool,APPLY_REGRESSION)
  IMPORT(bool,SAVE_FRAME)
  IMPORT(bool,SAVE_PLOT)
  IMPORT(bool,SAVE_PLOT_INDATA)
  IMPORT(bool,CUT_0_32_HZ)

  IMPORT(TString,OPLOT_DIR)
  IMPORT(TString,OFRAME_DIR)
  IMPORT(TString,H1_FRLIST_WITNESS)
  IMPORT(TString,H1_WITNESS_CHANNEL_LIST)

  IMPORT(int,nFILTER)
  IMPORT(double,LAYER_WIDTH)
  IMPORT(int,nPL)
  IMPORT(double,FWIDTH)

  IMPORT(double,APPLY_THRESHOLD)
  IMPORT(double,SOLVE_THRESHOLD)
  IMPORT(double,SOLVE_NEIGEN_PER_LAYER)
  IMPORT(char,SOLVE_REGULATOR)

  if(type==CWB_PLUGIN_DATA) {  

    // create output regression report dirs
    gSystem->Exec(TString("mkdir -p ")+TString(OPLOT_DIR));
    gSystem->Exec(TString("mkdir -p ")+TString(OFRAME_DIR));

    if(SAVE_PLOT_INDATA) {
      cout << "CWB_PLUGIN_DATA : ..." << endl;
      gROOT->SetBatch(true);

      char ofName[256];
      char gtitle[256];
      TString gfile;
      watplot plot(const_cast<char*>("plot"),200,20,800,500);

      // save psd cleaned/dirty data
      sprintf(gtitle,"Dirty Data %s",cfg->channelNamesRaw[0]);
      plot.goptions("alp logy", 1, 0., 0., true, 40, 240, true, 32);
      plot.gtitle(gtitle,"frequency (Hz)","strain/#sqrt{Hz}");
      sprintf(ofName,"oframes/PSD-H-H1_LDAS_C02_L2_ORIG-%lu-%lu.png",(int)x->start(),(int)x->stop()-(int)x->start());
      gfile=ofName;

      *x >> plot;
      plot >> gfile;
    }
  }

  if(type==CWB_PLUGIN_DATA_MDC) {  

    // get ifo index
    int xIFO =0;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==cfg->ifo[n]) {xIFO=n;break;}

    if(!APPLY_REGRESSION) { 
      if(ifo==cfg->ifo[xIFO]) {
        cout << "remove temporary file ..." << endl;
        TString jname = jfile->GetPath();
        jname.ReplaceAll(":/","");
        cout << jname.Data() << endl;
        gSystem->Exec(TString("rm "+jname).Data());
        cout << "end job" << endl;
        exit(0);
      }
      return;
    }

    cout << "CWB_PLUGIN_DATA_MDC : Apply Regression ..." << endl;

    CWB::Toolbox TB;

    Biorthogonal<double> Bio(512);
    WSeries<double> wB(Bio);

    TString frlist_witness;
    TString wit_channel_file_list; 
    if(ifo=="H1") {
      frlist_witness		= H1_FRLIST_WITNESS;
      wit_channel_file_list	= H1_WITNESS_CHANNEL_LIST; 
    } else {
      cout << "Error : ifo " << ifo.Data() << " not supported" << endl;  
      exit(1);
    }

    // target channel
    int level=x->getLevel();
    x->Inverse(-1);   
    wavearray<double> xx(x->size());
    xx.start(x->start());
    xx.rate(x->rate());
    for(int i=0;i<x->size();i++) xx.data[i]=x->data[i];

    if(CUT_0_32_HZ) {
      // set range [0-32]Hz to 0 (remove large dynamics at low frequency)
      int SR = int(xx.rate());
      int mm = 0;
      while (((SR % 2) == 0) && SR > 64) {SR /= 2;mm++;}
      wB.Forward(xx,mm);
      for(int i=0;i<1;i++) {wB.getLayer(xx,i);xx=0;wB.putLayer(xx,i);}
      wB.Inverse();
      wB.getLayer(xx,0);
    }

    // resample target to 2048Hz
    int sr = int(xx.rate());
    int nn = 0;
    while (((sr % 2) == 0) && sr > 2048) {sr /= 2;nn++;}
    wB.Forward(xx,nn);
    wB.getLayer(xx,0);

    int rrlevel=xx.rate()/(2*LAYER_WIDTH);
    WDM<double> WD(rrlevel, 2*rrlevel, 6, 12);
    double scratch = WD.m_H/xx.rate();
    if(scratch>cfg->segEdge+0.001) {
       cout << endl;
       cout << "Regression Plugin : Error - filter scratch must be <= cwb scratch!!!" << endl; 
       cout << "filter scratch : " << scratch << " sec" << endl;
       cout << "cwb    scratch : " << cfg->segEdge << " sec" << endl;
       exit(1);
    }

    WSeries<double> wT;

    // clean power lines
    frame frw(frlist_witness,"","README",true);
    wavearray<double> wit;
    wit.start(x->start()); wit.stop(x->stop());
    wT.Forward(xx,WD);
    regression rr(wT,"target");
    rr.mask(0,0.,xx.rate()/2.);
    for(int n=1;n<=nPL;n++) {double f=n*60.;rr.unmask(0,f-FWIDTH,f+FWIDTH);}
    // read witness channels
    ifstream iwit;
    iwit.open(wit_channel_file_list,ios::in);
    if (!iwit.good()) {cout << "Error Opening File : " << wit_channel_file_list.Data() << endl;exit(1);}
    char witness[1024];
    cout << endl;
    while(true) {
       iwit >> witness;
       if (!iwit.good()) break;
       if(witness[0]=='#') continue;
       cout << "read witness channel : \t" << witness << endl;
       frw.setChName(witness);
       frw >> wit;
       rr.add(wit,witness);
    }
    cout << endl;
    iwit.close();
    rr.setFilter(nFILTER);
    rr.setMatrix(cfg->segEdge,1.);
    rr.solve(SOLVE_THRESHOLD,SOLVE_NEIGEN_PER_LAYER,SOLVE_REGULATOR);
    rr.apply(APPLY_THRESHOLD);

    double OS=0;
    char ofName[512];
    char frName[512];
    char chName[512];

    // cget cleaned data
    wavearray<double> cc = rr.getClean();

    // save cleaned data into frame 
/*
    // restore original rate
    wB.putLayer(cc,0);
    wB.Inverse();
    wB.getLayer(cc,0);
*/
    // remove scratch 
    wavearray<double> CC = cc;
    OS = cfg->segEdge*cc.rate();
    CC.start(cc.start()+cfg->segEdge);
    CC.stop(cc.stop()-cfg->segEdge);
    CC.resize(cc.size()-2*OS);
    for(int i=0;i<CC.size()-2*OS;i++) CC[i]=CC[i+OS];

    // save noisy data into frame 
/*
    // restore original rate
    wB.Forward(cc,nn);
    wB.putLayer(xx,0);
    wB.Inverse();
    wB.getLayer(xx,0);
*/
    // remove scratch 
    wavearray<double> XX = xx;
    OS = cfg->segEdge*cc.rate();
    XX.start(xx.start()+cfg->segEdge);
    XX.stop(xx.stop()-cfg->segEdge);
    XX.resize(xx.size()-2*OS);
    for(int i=0;i<XX.size()-2*OS;i++) XX[i]=XX[i+OS];

    if(SAVE_PLOT) {
      gROOT->SetBatch(true);

      char gtitle[256];
      TString gfile;
      watplot plot(const_cast<char*>("plot"),200,20,800,500);

      double flow[4] ={40,55,115,175};
      double fhigh[4]={240,65,125,185};

      // save psd cleaned/dirty data
      for(int n=0;n<4;n++) {
        sprintf(gtitle,"Dirty/Cleaned Data %s - %3.0f-%3.0f Hz",cfg->channelNamesRaw[xIFO],flow[n],fhigh[n]);
        plot.gtitle(gtitle,"frequency (Hz)","strain/#sqrt{Hz}");
        plot.goptions("alp logy", 1, 0., 0., true, flow[n], fhigh[n], true, 32);
        sprintf(ofName,"%s/PSD-H-H1_LDAS_C02_L2-%lu-%lu-F-%3.0f-%3.0f.png",
                OPLOT_DIR.Data(),(int)CC.start(),(int)CC.stop()-(int)CC.start(),flow[n],fhigh[n]);
        gfile=ofName;
        XX >> plot; CC >> plot; plot >> gfile;
      }
    }

    if (SAVE_FRAME) {
      sprintf(chName,cfg->channelNamesRaw[xIFO]);
      sprintf(frName,"LHO_4k");
      sprintf(ofName,"%s/H-H1_LDAS_C02_L2_REG-%lu-%lu.gwf",
                     OFRAME_DIR.Data(),(int)CC.start(),(int)CC.stop()-(int)CC.start());
      cout << ofName << " " << CC.start() << " " << CC.stop()-CC.start() << endl;
      frame cfr(ofName,chName,"WRITE");
      cfr.setFrName(frName);
      CC >> cfr;
      cfr.close();

      sprintf(chName,cfg->channelNamesRaw[xIFO]);
      sprintf(frName,"LHO_4k");
      sprintf(ofName,"%s/H-H1_LDAS_C02_L2-%lu-%lu.gwf",
                     OFRAME_DIR.Data(),(int)XX.start(),(int)XX.stop()-(int)XX.start());
      cout << ofName << " " << XX.start() << " " << XX.stop()-XX.start() << endl;
      frame xfr(ofName,chName,"WRITE");
      xfr.setFrName(frName);
      XX >> xfr;
      xfr.close();
    }  
/*
    // restore original rate
    wB.putLayer(cc,0);
    wB.Inverse();
    wB.getLayer(cc,0);
    // return in x the cleaned data
    for(int i=0;i<x->size();i++) x->data[i]=cc.data[i];
    x->Forward(level);
*/
    // end job - clean temporary files
    if(ifo==cfg->ifo[xIFO]) {
      cout << "remove temporary file ..." << endl;
      TString jname = jfile->GetPath();
      jname.ReplaceAll(":/","");
      cout << jname.Data() << endl;
      gSystem->Exec(TString("rm "+jname).Data());
      cout << "end job" << endl;
      exit(0);
    }
  }

  return;
}
