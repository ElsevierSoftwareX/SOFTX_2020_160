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

#define NCH 1000

using namespace CWB;

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!REGRESSION
// Plugin for linear/bilinear ranking analysis for one detector

  cout << endl;
  cout << "-----> CWB_Plugin_Linear_Bilinear_Regression_Rank.C : " << ifo.Data() << endl;
  cout << endl;

  //---------------------------------------------------------------------
  // MAIN REGRESSION RANK
  //---------------------------------------------------------------------

  bool SAVE_TREE			= false;		// if true then results are saved to tree
  bool SAVE_ASCII			= false;		// if true then results are saved to ascii file

  TString REGRESSION_RANK               = "LINEAR";		// regression rank [LINEAR/BILINEAR]
  int  CUT_LOW_FREQ                     = 32;			// if > 0 the data in [0:CUT_LOW_FREQ] are set to 0

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS
  //---------------------------------------------------------------------

  TString FRLIST_WITNESS                = "";           // list of witness frame files
  TString CHLIST_LINEAR                 = "";           // list of channel names used for linear regression
  TString CHLIST_WITNESS                = "";           // list of channel names used for regression rank

  TString CHNAME_LINEAR                 = "";           // channel name used to made the bilinear channels

  int     RESAMPLING_INDEX              = 11;           // resample target/witness channels to pow(2,RESAMPLING_INDEX)

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

  IMPORT(bool,SAVE_TREE)
  IMPORT(bool,SAVE_ASCII)

  IMPORT(TString,REGRESSION_RANK)
  IMPORT(int,CUT_LOW_FREQ)
  IMPORT(int,RESAMPLING_INDEX)
  IMPORT(double,WFLOW)
  IMPORT(double,WFHIGH)

  IMPORT(TString,FRLIST_WITNESS)
  IMPORT(TString,CHLIST_LINEAR)
  IMPORT(TString,CHLIST_WITNESS)
  IMPORT(TString,CHNAME_LINEAR)

  if(FRLIST_WITNESS=="")  {cout << "Error : witness frames list not defined" << endl;gSystem->Exit(1);}
  if(CHLIST_LINEAR=="")   {cout << "Error : linear witness channels list not defined" << endl;gSystem->Exit(1);}
  if(CHLIST_WITNESS=="")  {cout << "Error : witness channels list not defined" << endl;gSystem->Exit(1);}
  if(CHNAME_LINEAR=="")   {cout << "Error : main linear witness channel list not defined" << endl;gSystem->Exit(1);}

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

  REGRESSION_RANK.ToUpper();
  if(REGRESSION_RANK!="LINEAR" && REGRESSION_RANK!="BILINEAR")  {
    cout << "Error : REGRESSION_RANK not defined [LINEAR/BILINEAR]" << endl;
    gSystem->Exit(1);
  }

  if(type==CWB_PLUGIN_DATA_MDC) {  

    int rank[4] = {50,30,10,1};

    int nCH=0;
    std::vector<std::string> CHN;
    std::vector<std::string>* pCHN=&CHN;
    double CHC[4][NCH];

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
    wavearray <double> e[4];
    int size=0;

    // target channel
    int level=x->getLevel();
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
    // LINEAR RANKING
    // --------------------------------------------------------------

    cout << "CWB_PLUGIN_DATA_MDC : Apply Linear Regression ..." << endl;

    frame frl(FRLIST_WITNESS,"","README",true);
    wavearray<double> wl;
    wl.start(x->start()); wl.stop(x->stop());
    wT.Forward(xx,WD);

    if(REGRESSION_RANK=="LINEAR") {
      cout << "CWB_PLUGIN_DATA_MDC : Apply Linear Regression Rank ..." << endl;

      // open witness channel list
      ifstream ifl;
      ifl.open(CHLIST_WITNESS.Data(),ios::in);
      if (!ifl.good()) {cout << "Error Opening File : " << CHLIST_WITNESS.Data() << endl;gSystem->Exit(1);}
  
      char witness[1024];
      cout << endl;

      char ofname[256];
      sprintf(ofname,"%s/linear_regression_rank_%s_%s.txt",cfg->output_dir,ifo.Data(),flabel);
      cout << "write results to " << ofname << endl;
      ofstream out;
      if(SAVE_ASCII) {
        out.open(ofname);
        if (!out.good()) {cout << "Error Opening File : " << ofname << endl;gSystem->Exit(1);}
      }
      nCH=0;
      char fstring[512];
      while(true) {
        ifl >> witness;
        if (!ifl.good()) break;
        if(witness[0]=='#') continue;
        //cout << "read witness channel : \t" << witness << endl;
        frl.setChName(witness);
        frl.setSRIndex(RESAMPLING_INDEX);
        frl >> wl;
        size = rr.add(wT,"target");
        if(size==0) {cout << "Regression Plugin - empty target channel" << endl;gSystem->Exit(1);}
        rr.mask(0,0.,xx.rate()/2.);
        for(int n=lPOWERLINE;n<=hPOWERLINE;n++) {double f=n*fPOWERLINE;rr.unmask(0,f-FWIDTH,f+FWIDTH);}
        size = rr.add(wl,witness);
        if(size>0) {
          rr.setFilter(L_NFILTER);
          rr.setMatrix(cfg->segEdge,1.);
          rr.solve(L_SOLVE_THRESHOLD,L_SOLVE_NEIGEN_PER_LAYER,L_SOLVE_REGULATOR);
          rr.apply(L_APPLY_THRESHOLD);

          for(int n=0;n<4;n++) {
            e[n] = rr.rank(rank[n],0);
            CHC[n][nCH]=e[n][0];
          }
        } else {
          for(int n=0;n<4;n++) {
            CHC[n][nCH]=-1;
          }
        }
        CHN.push_back(witness);
        sprintf(fstring,"%3d %35s %10.5f %9.4f %8.3f %7.2f",
                nCH,witness,CHC[0][nCH],CHC[1][nCH],CHC[2][nCH],CHC[3][nCH]);
        cout << fstring << endl;
        if(SAVE_ASCII) out << fstring << endl;out.flush();
        nCH++;
      }
      cout << endl;
      ifl.close();
      if(SAVE_ASCII) out.close();
    } else {

      // open linear channel list
      ifstream ifl;
      ifl.open(CHLIST_LINEAR.Data(),ios::in);
      if (!ifl.good()) {cout << "Error Opening File : " << CHLIST_LINEAR.Data() << endl;gSystem->Exit(1);}
  
      char linear[1024];
      cout << endl;

      size = rr.add(wT,"target");
      if(size==0) {cout << "Regression Plugin - empty target channel" << endl;gSystem->Exit(1);}
      rr.mask(0,0.,xx.rate()/2.);
      for(int n=lPOWERLINE;n<=hPOWERLINE;n++) {double f=n*fPOWERLINE;rr.unmask(0,f-FWIDTH,f+FWIDTH);}
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
    }
    wl.resize(0);

    // --------------------------------------------------------------
    // BILINEAR RANKING
    // --------------------------------------------------------------

    if(REGRESSION_RANK=="BILINEAR") {

      cout << "CWB_PLUGIN_DATA_MDC : Apply Bilinear Regression Rank ..." << endl;

      char ofname[256];
      sprintf(ofname,"%s/bilinear_regression_rank_%s_%s.txt",cfg->output_dir,ifo.Data(),flabel);
      cout << "write results to " << ofname << endl;
      ofstream out;
      if(SAVE_ASCII) {
        out.open(ofname);
        if (!out.good()) {cout << "Error Opening File : " << ofname << endl;gSystem->Exit(1);}
      }
      wavearray<double> yy = rr.getClean();
      wT.Forward(yy,WD);

      // open main linear frame list
      wavearray<double> ml;
      ml.start(x->start()); ml.stop(x->stop());
      frame frb(FRLIST_WITNESS,"","README",true);

      // open witness channel list
      ifstream iwb;
      iwb.open(CHLIST_WITNESS.Data(),ios::in);
      if (!iwb.good()) {cout << "Error Opening File : " << CHLIST_WITNESS.Data() << endl;gSystem->Exit(2);}

      char fstring[512];
      char witness[1024];
      wavearray<double> wb;
      wb.start(x->start()); wb.stop(x->stop());
      cout << endl;
      nCH=0;
      while(true) {
        iwb >> witness;
        if (!iwb.good()) break;
        //cout << "read witness channel : \t" << witness << endl;
        if(witness[0]=='#') continue;
        frb.setChName(witness);
        frb.setSRIndex(RESAMPLING_INDEX);
        frb >> wb;

        rr.clear();
        size=rr.add(wT,"target");
        if(size==0) {cout << "Regression Plugin - empty target channel" << endl;gSystem->Exit(1);}
        rr.mask(0,0.,yy.rate()/2.);
        for(int n=lPOWERLINE;n<=hPOWERLINE;n++) {double f=n*fPOWERLINE;rr.unmask(0,f-FWIDTH,f+FWIDTH);}

        // add main linear channel for bicoherence removal
        //cout << endl;
        //cout << "read main linear channel : \t" << CHNAME_LINEAR.Data() << endl << endl;
        frb.setChName(CHNAME_LINEAR);
        frb.setSRIndex(RESAMPLING_INDEX);
        frb >> ml;
        size=rr.add(ml,const_cast<char*>(CHNAME_LINEAR.Data()));
        if(size==0) {cout << "Regression Plugin - empty main linear channel" << endl;gSystem->Exit(1);}

        if(rr.add(wb,witness,WFLOW,WFHIGH)) {
          rr.add(1,2,"biwitness");
          rr.mask(1); rr.mask(2);

          rr.setFilter(B_NFILTER);
          rr.setMatrix(cfg->segEdge,1.);
          rr.solve(B_SOLVE_THRESHOLD,B_SOLVE_NEIGEN_PER_LAYER,B_SOLVE_REGULATOR);
          rr.apply(B_APPLY_THRESHOLD);

          for(int n=0;n<4;n++) {
            e[n] = rr.rank(rank[n],0);
            CHC[n][nCH]=e[n][0];
          }
        } else {
          for(int n=0;n<4;n++) {
            CHC[n][nCH]=-1;
          }
        }
        CHN.push_back(witness);
        nCH++;
        sprintf(fstring,"%3d %35s %10.5f %9.4f %8.3f %7.2f",
                nCH,witness,CHC[0][nCH],CHC[1][nCH],CHC[2][nCH],CHC[3][nCH]);
        cout << fstring << endl;
        if(SAVE_ASCII) out << fstring << endl;out.flush();
      }
      iwb.close();
      if(SAVE_ASCII) out.close();
      cout << endl;

      ml.resize(0);
      wb.resize(0);
    }

    if(SAVE_TREE) {
      // save data to root
      char outFile[512];
      if(REGRESSION_RANK=="LINEAR") {
        sprintf(outFile,"%s/linear_regression_rank_%s_%s.root",cfg->output_dir,ifo.Data(),flabel);
      } else {
        sprintf(outFile,"%s/bilinear_regression_rank_%s_%s.root",cfg->output_dir,ifo.Data(),flabel);
      }
      cout << "write results to " << outFile << endl;

      TFile* ofile = new TFile(outFile,"RECREATE");
      cfg->Write("config"); 			// write config object
      //cfg->configPlugin.Write("configPlugin"); // write configPlugin object

      TTree rrtree("regression","regression");
      double igps=x->start();
      double start=x->start();
      double stop=x->stop();
      char   sifo[256];sprintf(sifo,"%s",ifo.Data());
      rrtree.Branch("run",&runID,"run/I");
      rrtree.Branch("gps",&igps,"gps/D");
      rrtree.Branch("start",&start,"start/D");
      rrtree.Branch("stop",&stop,"stop/D");
      rrtree.Branch("sifo",sifo,"sifo/C");
      rrtree.Branch("ifo",&xIFO,"ifo/I");
      rrtree.Branch("nch",&nCH,"nch/I");
      rrtree.Branch("chc50",CHC[0],"chc[nch]/D");
      rrtree.Branch("chc30",CHC[1],"chc[nch]/D");
      rrtree.Branch("chc10",CHC[2],"chc[nch]/D");
      rrtree.Branch("chc01",CHC[3],"chc[nch]/D");
      rrtree.Branch("chn",&pCHN);

      rrtree.Fill();
      rrtree.Write();
      ofile->Close();
    }

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
  }

  return;
}
