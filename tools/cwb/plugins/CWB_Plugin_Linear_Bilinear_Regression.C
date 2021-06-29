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
//!REGRESSION
// Plugin for linear/bilinear regression analysis for one detector

  cout << endl;
  cout << "-----> CWB_Plugin_Linear_Bilinear_Regression.C : " << ifo.Data() << endl;
  cout << endl;

  char cmd[128];
  sprintf(cmd,"int type = %d;",type);
  gROOT->ProcessLine(cmd);
  sprintf(cmd,"network* net = (network*)%p;",net);
  gROOT->ProcessLine(cmd);

  //---------------------------------------------------------------------
  // MAIN REGRESSION
  //---------------------------------------------------------------------

  bool EXIT_AFTER_REGRESSION            = true;         // if true cwb terminate after the regression analysis
  bool APPLY_LINEAR_REGRESSION          = true;         // if true then linear regression is applied
  bool APPLY_BILINEAR_REGRESSION        = true;         // if true then blinear regression is applied;
  bool SAVE_INFRAME                     = false;        // if true the input data are saved to frame files
  bool SAVE_OUTFRAME                    = true;         // if true the cleaned data are saved to frame files
  bool SAVE_PSD_PLOT                    = false;        // if true psd cleaned vs noisy data are saved to png
  bool SAVE_EIGEN_PLOT                  = false;        // if true the regression filter eigenvaules are saved to png
  int  CUT_LOW_FREQ                     = 32;           // if > 0 the data in [0:CUT_LOW_FREQ] are set to 0

  TString OFRDIR                        = "oframes";    // output frame files directory
  TString FRNAME                        = "";           // frame name
  TString FRLABEL                       = "";           // name used to label the output frames

  bool DISABLE_MDC_FROM_FRAMES		= true;		// if true then disable read mdc from frames
  bool DISABLE_MDC_FROM_PLUGIN		= false;	// if true then disable read mdc from plugin

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS
  //---------------------------------------------------------------------

  TString FRLIST_WITNESS                = "";           // list of witness frame files
  TString CHLIST_LINEAR                 = "";           // list of channel names used for linear regression
  TString CHLIST_BILINEAR               = "";           // list of channel names used for bilinear regression

  TString CHNAME_LINEAR                 = "";           // channel name used to made the bilinear channels

  int	  RESAMPLING_INDEX		= 11;		// resample target/witness channels to pow(2,RESAMPLING_INDEX)

  double  WFLOW                         = 0;            // lower frequency used to made the bilinear channels
  double  WFHIGH                        = 10;           // high  frequency used to made the bilinear channels

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS
  //---------------------------------------------------------------------

  double LAYER_WIDTH                    = 1.0;          // frequency layer resolution used in the regression analysis
  double fPOWERLINE                     = 60.0;         // powerline frequency (Hz)
  int    lPOWERLINE                     = 1;            // low power line harmonic (lPOWERLINE*fPOWERLINE)
  int    hPOWERLINE                     = 3;            // high power line harmonic (hPOWERLINE*fPOWERLINE)
  double FWIDTH                         = 5;            // frequency width used by regression

  // PARAMETERS FOR LINEAR REGRESSION

  int    L_NFILTER                      = 5;            // half-size length of a unit filter (setFilter)
  double L_APPLY_THRESHOLD              = 0.2;          // threshold used in apply
  double L_SOLVE_THRESHOLD              = 0.0;          // eigenvalue threshold (solve)
  double L_SOLVE_NEIGEN_PER_LAYER       = 0;            // number of selected eigenvalues (solve)
  char   L_SOLVE_REGULATOR              = 'h';          // regulator (solve)

  // PARAMETERS FOR BILINEAR REGRESSION

  int    B_NFILTER                      = 5;            // half-size length of a unit filter (setFilter)
  double B_APPLY_THRESHOLD              = 0.2;          // threshold used in apply
  double B_SOLVE_THRESHOLD              = 0.0;          // eigenvalue threshold (solve)
  double B_SOLVE_NEIGEN_PER_LAYER       = 0;            // number of selected eigenvalues (solve)
  char   B_SOLVE_REGULATOR              = 'h';          // regulator (solve)

  // ---------------------------------
  // read plugin config
  // ---------------------------------

  cfg->configPlugin.Exec();

  IMPORT(bool,EXIT_AFTER_REGRESSION)
  IMPORT(bool,APPLY_LINEAR_REGRESSION)
  IMPORT(bool,APPLY_BILINEAR_REGRESSION)
  IMPORT(bool,SAVE_INFRAME)
  IMPORT(bool,SAVE_OUTFRAME)
  IMPORT(bool,SAVE_PSD_PLOT)
  IMPORT(bool,SAVE_EIGEN_PLOT)
  IMPORT(int,CUT_LOW_FREQ)
  IMPORT(bool,DISABLE_MDC_FROM_FRAMES)
  IMPORT(bool,DISABLE_MDC_FROM_PLUGIN)
  IMPORT(int,RESAMPLING_INDEX)
  IMPORT(double,WFLOW)
  IMPORT(double,WFHIGH)

  IMPORT(TString,OFRDIR)
  IMPORT(TString,FRNAME)
  IMPORT(TString,FRLABEL)
  IMPORT(TString,FRLIST_WITNESS)
  IMPORT(TString,CHLIST_LINEAR)
  IMPORT(TString,CHLIST_BILINEAR)
  IMPORT(TString,CHNAME_LINEAR)

  if(FRLIST_WITNESS=="")  {cout << "Error : witness frames list not defined" << endl;gSystem->Exit(1);}
  if(CHLIST_LINEAR=="")   {cout << "Error : linear witness channels list not defined" << endl;gSystem->Exit(1);}
  if(CHLIST_BILINEAR=="") {cout << "Error : bilinear witness channels list not defined" << endl;gSystem->Exit(1);}
  if(CHNAME_LINEAR=="")   {cout << "Error : main linear witness channel list not defined" << endl;gSystem->Exit(1);}
  if(FRNAME=="")          {cout << "Error : frame name not defined" << endl;gSystem->Exit(1);}
  if(FRLABEL=="")         {cout << "Error : frame label not defined" << endl;gSystem->Exit(1);}

  IMPORT(double,LAYER_WIDTH)
  IMPORT(double,fPOWERLINE)
  IMPORT(int,lPOWERLINE)
  IMPORT(int,hPOWERLINE)
  IMPORT(double,FWIDTH)

  IMPORT(int,L_NFILTER)
  IMPORT(double,L_APPLY_THRESHOLD)
  IMPORT(double,L_SOLVE_THRESHOLD)
  IMPORT(double,L_SOLVE_NEIGEN_PER_LAYER)
  IMPORT(char,L_SOLVE_REGULATOR)

  IMPORT(int,B_NFILTER)
  IMPORT(double,B_APPLY_THRESHOLD)
  IMPORT(double,B_SOLVE_THRESHOLD)
  IMPORT(double,B_SOLVE_NEIGEN_PER_LAYER)
  IMPORT(char,B_SOLVE_REGULATOR)

  if(!APPLY_LINEAR_REGRESSION && !APPLY_BILINEAR_REGRESSION) {
    cout << "Error : regression type [linear/bilinear] is not defined" << endl;
    gSystem->Exit(1);
  } 

  if(type==CWB_PLUGIN_CONFIG) {
    cfg->dataPlugin=false; 			// disable read data from frames
    cfg->mdcPlugin=DISABLE_MDC_FROM_FRAMES;   	// disable read mdc from frames
  }

  if(type==CWB_PLUGIN_MDC) {

    if(DISABLE_MDC_FROM_PLUGIN) {*x=0.;return;}

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


  if(type==CWB_PLUGIN_DATA_MDC) {  

    // create output regression report dirs
    gSystem->Exec(TString("mkdir -p ")+TString(OFRDIR));

    // get ifo index
    int xIFO =0;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==cfg->ifo[n]) {xIFO=n;break;}

    int runID = net->nRun;
    char flabel[512];
    int Tb=x->start();
    int dT=x->size()/x->rate();
    sprintf(flabel,"%d_%d_%s_job%d",int(Tb),int(dT),cfg->data_label,runID);

    Biorthogonal<double> Bio(512);
    WSeries<double> wB(Bio);
    WSeries<double> wT;
    int size=0;

    // target channel
    int level=x->getLevel();    // save the input decomposition level
    x->Inverse(-1);   
    wavearray<double> xx(x->size());
    xx.start(x->start());
    xx.rate(x->rate());
    for(int i=0;i<x->size();i++) xx.data[i]=x->data[i];

    if(CUT_LOW_FREQ>0) {
      // set range [0-CUT_LOW_FREQ]Hz to 0 (remove large dynamics at low frequency)
      int SR = int(xx.rate());
      int mm = 0;
      while (((SR % 2) == 0) && SR > 2*CUT_LOW_FREQ) {SR /= 2;mm++;}
      wB.Forward(xx,mm);
      for(int i=0;i<1;i++) {wB.getLayer(xx,i);xx=0;wB.putLayer(xx,i);}
      wB.Inverse();
      wB.getLayer(xx,0);
      cout << "Set Frequency Range [0:" << SR/2 << "] = 0" << endl;
    }

    // resample target to pow(2,RESAMPLING_INDEX) Hz
    int sr = int(xx.rate());
    int nn = 0;
    while (((sr % 2) == 0) && sr > (1<<RESAMPLING_INDEX)) {sr /= 2;nn++;}
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
       gSystem->Exit(1);
    }

    regression rr;

    // --------------------------------------------------------------
    // LINEAR CLEANING
    // --------------------------------------------------------------

    if(APPLY_LINEAR_REGRESSION) { 

      cout << "CWB_PLUGIN_DATA_MDC : Apply Linear Regression ..." << endl;

      frame frl(FRLIST_WITNESS,"","README",true);
      wavearray<double> wl;
      wl.start(x->start()); wl.stop(x->stop());
      wT.Forward(xx,WD);
      size = rr.add(wT,"target");
      if(size==0) {cout << "Regression Plugin - empty target channel" << endl;gSystem->Exit(1);}
      rr.mask(0,0.,xx.rate()/2.);
      for(int n=lPOWERLINE;n<=hPOWERLINE;n++) {double f=n*fPOWERLINE;rr.unmask(0,f-FWIDTH,f+FWIDTH);}

      // open linear channel list
      ifstream ifl;
      ifl.open(CHLIST_LINEAR.Data(),ios::in);
      if (!ifl.good()) {cout << "Error Opening File : " << CHLIST_LINEAR.Data() << endl;gSystem->Exit(1);}
  
      char linear[1024];
      cout << endl;
      while(true) {
        ifl >> linear;
        if (!ifl.good()) break;
        if(linear[0]=='#') continue;
        cout << "read linear channel : \t" << linear << endl;
        frl.setChName(linear);
        frl.setSRIndex(RESAMPLING_INDEX);
        frl >> wl;
        rr.add(wl,linear);
      }
      cout << endl;
      ifl.close();
      rr.setFilter(L_NFILTER);
      rr.setMatrix(cfg->segEdge,1.);
      rr.solve(L_SOLVE_THRESHOLD,L_SOLVE_NEIGEN_PER_LAYER,L_SOLVE_REGULATOR);
      rr.apply(L_APPLY_THRESHOLD);
      wl.resize(0);
      if(SAVE_EIGEN_PLOT) {
        gROOT->SetBatch(true);
        char ofName[256];
        char gtitle[256];
        
        wavearray<double> eigen = rr.getVEIGEN(0);
        watplot eplot(const_cast<char*>("eplot"),200,20,800,500);
        sprintf(gtitle,"Linear Regression Filter eigenvalues");
        eplot.gtitle(gtitle,"filter index","eigenvalue");
        eplot.goptions("alp", 1, 0., 0.);

        sprintf(ofName,"%s/eigen_linear_regression_%s_%s.png",cfg->dump_dir,ifo.Data(),flabel);
        cout << "write results to " << ofName << endl;

        TString gfile=ofName;
        eigen >> eplot; eplot >> gfile;
      }
    }

    // --------------------------------------------------------------
    // BILINEAR CLEANING
    // --------------------------------------------------------------

    if(APPLY_BILINEAR_REGRESSION) { 

      cout << "CWB_PLUGIN_DATA_MDC : Apply Bilinear Regression ..." << endl;

      wavearray<double> yy = rr.getClean();
      wT.Forward(yy,WD);
      rr.clear();
      size=rr.add(wT,"target");
      if(size==0) {cout << "Regression Plugin - empty target channel" << endl;gSystem->Exit(1);}
      rr.mask(0,0.,yy.rate()/2.);
      for(int n=lPOWERLINE;n<=hPOWERLINE;n++) {double f=n*fPOWERLINE;rr.unmask(0,f-FWIDTH,f+FWIDTH);}

      // add main linear channel for bicoherence removal
      cout << endl;
      cout << "read main linear channel : \t" << CHNAME_LINEAR.Data() << endl << endl;
      wavearray<double> ml;
      ml.start(x->start()); ml.stop(x->stop());
      frame frb(FRLIST_WITNESS,"","README",true);
      frb.setChName(CHNAME_LINEAR);
      frb.setSRIndex(RESAMPLING_INDEX);
      frb >> ml;
      size=rr.add(ml,const_cast<char*>(CHNAME_LINEAR.Data()));
      if(size==0) {cout << "Regression Plugin - empty main linear channel" << endl;gSystem->Exit(1);}

      // open bilinear channel list
      ifstream iwb;
      iwb.open(CHLIST_BILINEAR.Data(),ios::in);
      if (!iwb.good()) {cout << "Error Opening File : " << CHLIST_BILINEAR.Data() << endl;gSystem->Exit(2);}

      char bilinear[1024];
      wavearray<double> wb;
      wb.start(x->start()); wb.stop(x->stop());
      int nBILINEAR=0;
      cout << endl;
      while(true) {
        iwb >> bilinear;
        if (!iwb.good()) break;
        if(bilinear[0]=='#') continue;
        cout << "read bilinear channel : \t" << bilinear << endl;
        frb.setChName(bilinear);
        frb.setSRIndex(RESAMPLING_INDEX);
        frb >> wb;
        if(!rr.add(wb,bilinear,WFLOW,WFHIGH)) continue;
        nBILINEAR++;
      }
      iwb.close();
      cout << endl;
      for(int n=2;n<=nBILINEAR+1;n++) rr.add(1,n,"bilinear");
      for(int n=1;n<=nBILINEAR+1;n++) rr.mask(n);

      rr.setFilter(B_NFILTER);
      rr.setMatrix(cfg->segEdge,1.);
      rr.solve(B_SOLVE_THRESHOLD,B_SOLVE_NEIGEN_PER_LAYER,B_SOLVE_REGULATOR);
      rr.apply(B_APPLY_THRESHOLD);
      ml.resize(0);
      wb.resize(0);
      if(SAVE_EIGEN_PLOT) {
        gROOT->SetBatch(true);
        char ofName[256];
        char gtitle[256];
        
        wavearray<double> eigen = rr.getVEIGEN(0);
        watplot eplot(const_cast<char*>("eplot"),200,20,800,500);
        sprintf(gtitle,"Linear Regression Filter eigenvalues");
        eplot.gtitle(gtitle,"filter index","eigenvalue");
        eplot.goptions("alp", 1, 0., 0.);

        sprintf(ofName,"%s/eigen_bilinear_regression_%s_%s.png",cfg->dump_dir,ifo.Data(),flabel);
        cout << "write results to " << ofName << endl;

        TString gfile=ofName;
        eigen >> eplot; eplot >> gfile;
      }
    }

    // get cleaned data
    wavearray<double> cc = rr.getClean();

    if(SAVE_PSD_PLOT) {
      gROOT->SetBatch(true);

      char ofName[512];
      char gtitle[256];
      TString gfile;
      watplot plot(const_cast<char*>("plot"),200,20,800,500);

      double tstart = xx.start()+cfg->segEdge;
      double tstop  = xx.stop()-cfg->segEdge;

      // save psd cleaned/dirty data
      for(int n=lPOWERLINE;n<=hPOWERLINE;n++) {
        double flow=n*fPOWERLINE-FWIDTH;
        double fhigh=n*fPOWERLINE+FWIDTH;
        sprintf(gtitle,"Dirty/Cleaned Data %s - %3.0f-%3.0f Hz",
                cfg->channelNamesRaw[xIFO],flow,fhigh);
        plot.gtitle(gtitle,"frequency (Hz)","strain/#sqrt{Hz}");
        plot.goptions("alp logy", 1, tstart, tstop, true, flow,fhigh, true, 32);

        sprintf(ofName,"%s/psd_regression_%s_%s_F%2.0f_F%2.0f.png",cfg->dump_dir,ifo.Data(),flabel,flow,fhigh);
        cout << "write results to " << ofName << endl;

        gfile=ofName;
        xx >> plot; cc >> plot; plot >> gfile;
      }
    }

    if (SAVE_INFRAME) {
      double OS=0;
      char frName[512];
      char chName[512];
      char ofName[512];

      // save noisy data into frame 
      wavearray<double> XX = xx;
/*
      // restore original rate
      wB.Forward(cc,nn);
      wB.putLayer(xx,0);
      wB.Inverse();
      wB.getLayer(XX,0);
*/
      // remove scratch 
      OS = cfg->segEdge*cc.rate();
      XX.start(xx.start()+cfg->segEdge);
      XX.stop(xx.stop()-cfg->segEdge);
      XX.resize(xx.size()-2*OS);
      for(int i=0;i<XX.size()-2*OS;i++) XX[i]=XX[i+OS];

      sprintf(chName,cfg->channelNamesRaw[xIFO]);
      sprintf(frName,FRNAME.Data());
      sprintf(ofName,"%s/I-%s-%lu-%lu.gwf", 
                     OFRDIR.Data(),FRLABEL.Data(),(int)XX.start(),(int)XX.stop()-(int)XX.start());
      cout << "write frame file " << ofName << endl;
      frame xfr(ofName,chName,"WRITE");
      xfr.setFrName(frName);
      XX >> xfr;
      xfr.close();
      XX.resize(0);
    }  

    if (SAVE_OUTFRAME) {
      double OS=0;
      char frName[512];
      char chName[512];
      char ofName[512];

      // save cleaned data into frame 
      wavearray<double> CC = cc;
/*
      // restore original rate
      wB.putLayer(cc,0);
      wB.Inverse();
      wB.getLayer(CC,0);
*/
      // remove scratch 
      OS = cfg->segEdge*cc.rate();
      CC.start(cc.start()+cfg->segEdge);
      CC.stop(cc.stop()-cfg->segEdge);
      CC.resize(cc.size()-2*OS);
      for(int i=0;i<CC.size()-2*OS;i++) CC[i]=CC[i+OS];

      sprintf(chName,cfg->channelNamesRaw[xIFO]);
      sprintf(frName,FRNAME.Data());
      sprintf(ofName,"%s/R-%s-%lu-%lu.gwf", 
                     OFRDIR.Data(),FRLABEL.Data(),(int)CC.start(),(int)CC.stop()-(int)CC.start());
      cout << "write frame file " << ofName << endl;
      frame cfr(ofName,chName,"WRITE");
      cfr.setFrName(frName);
      CC >> cfr;
      cfr.close();
      CC.resize(0);
    }

    if(EXIT_AFTER_REGRESSION) {
      // end job - clean temporary files
      if(ifo==cfg->ifo[xIFO]) {
        cout << "remove temporary file ..." << endl;
        TString jname = jfile->GetPath();
        jname.ReplaceAll(":/","");
        cout << jname.Data() << endl;
        gSystem->Exec(TString("rm "+jname).Data());
        cout << "end job" << endl;
        gSystem->Exit(0);
      }
    } else {
      // restore original rate
      wB.putLayer(cc,0);
      wB.Inverse();
      wB.getLayer(cc,0);
      // return the cleaned data
      for(int i=0;i<x->size();i++) x->data[i]=cc.data[i];
      // restore original decomposition level
      x->Forward(level);
    }
  }

  return;
}
