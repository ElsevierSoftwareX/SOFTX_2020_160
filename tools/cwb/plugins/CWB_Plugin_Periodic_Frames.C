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


//!STRAIN_DATA_MANIPULATION_
// Plugin used to replicate periodicaly the input strain data 

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
#include "cwb2G.hh"
#include <vector>


// ---------------------------------------------------------------------------------
// HOW TO SETUP PLUGIN STRAIN IN CWB USER CONFIGURATION (EXAMPLE)
// ---------------------------------------------------------------------------------

// Set veto in the range [strain_veto_time-strain_veto_window/2, strain_veto_time+strain_veto_window/2]

/*
  // Example for the GW150914 event

  TString optstrain = "";                       	// NOTE : add space at the end of each line

  optstrain += "strain_veto_time=1126259462.404 ";      // Veto GPS time
  optstrain += "strain_veto_window=10 ";                // Veto window time

  strcpy(parPlugin,optstrain.Data());           	// set STRAIN plugin parameters
  strcpy(comment,"STRAIN");		    		// optional

*/

#define STRAIN_VETO_TIME             0.0             
#define STRAIN_VETO_WINDOW           0.0             

// ---------------------------------------------------------------------------------
// USER STRAIN PLUGIN OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  double veto_time;
  double veto_window;
};

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void SetCAT2(TString ifo, CWB::config* cfg, double gps, double window=0.0);
void SetTimeShiftCAT2(TString ifo, CWB::config* cfg, double shift);

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions();


// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions gOPT;                           // global User Options
double strainShift;


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  // this plugin permits to replicate periodicaly in time the same set of data declared in input frame list files
  //
  // the GPS range P=[frRange.start, frRange.stop] defines the Start data period to be used 
  //
  //                    frRange.start   frRange.stop
  //                    ^               ^ 
  // ...................|xxxxxx P xxxxxx|............................
  //
  // the period P is replicated back and forth starting from startFramesGPS
  //
  //                    frRange.start                    
  //                    ^                                  
  // xxx|xxxxxx P xxxxxx|xxxxxx P xxxxxx|xxxxxx P xxxxxx|xxxxxx P xxx
  //                                  
  // = J ==|== J ==|== J ==|== J ==|== J ==|== J ==|== J ==|== J ==|=
  //
  //   N       Y       N       Y       N       Y       N       Y       
  //                                  
  // WARNING: The jobs (J) with range overlapping 2 periods do not work with this method 
  //          It is necessary to use a procedure which select only the working jobs
  //                                  


  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->dataPlugin=true;         			// disable read strain from main pipeline

    ResetUserOptions();                                 // set default config options
    ReadUserOptions(cfg->parPlugin);                    // user config options : read from parPlugin
    PrintUserOptions();                                 // print config options
  }

  if(type==CWB_PLUGIN_INIT_JOB) {

    CWB::Toolbox TB;

    int ifoID=0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {ifoID=n;break;}

    int segEdge  = int(cfg->segEdge);

    // get data range
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
    double xstart = gCWB2G->Tb-segEdge;
    double xstop  = gCWB2G->Te+segEdge;
    cout << "CWB_Plugin_Periodic_Frames.C - " << "xstart : " << int(xstart) << " xstop : " << int(xstop) << endl;

    if(TString(cfg->frFiles[ifoID])=="") {
      cout << "CWB_Plugin_Periodic_Frames.C - Error : strain frame files list name is not defined" << endl;
      exit(1);
    }

    // open frame file
    CWB::frame fr;
    fr.open(cfg->frFiles[ifoID]);
    fr.setVerbose();
    fr.setRetryTime(cfg->frRetryTime);
    fr.setSRIndex(TMath::Nint(TMath::Log2(cfg->inRate))); 	// force data to be resampled to inRate
    int nfrFiles = fr.getNfiles();
    cout << "CWB_Plugin_Periodic_Frames.C - Strain " << " -> nfrFiles : " << nfrFiles << endl;

    // read strain range from strain frl files
    waveSegment frRange  = fr.getFrRange();
    cout << "CWB_Plugin_Periodic_Frames.C - frRange : " << frRange.start << " " << frRange.stop << endl;

    // compute strain shift time 
    double mlength = frRange.stop-frRange.start;		// frame files time range
    double toffset = xstart-frRange.start;
    int nShift = toffset>=0 ? int(toffset/mlength) : int(toffset/mlength)-1;
    strainShift = mlength*nShift;
    cout << "strainShift : " << strainShift << endl;

    // ------------------------------
    // init cat2 segments
    // ------------------------------

    if(gOPT.veto_time>0) SetCAT2(ifo, cfg, gOPT.veto_time);	// set veto CAT2 defined by the user
    SetTimeShiftCAT2(ifo, cfg, strainShift);			// apply strainShift to CAT2

    // set veto CAT2 with 2 sec window to the edges of frames
    SetCAT2(ifo, cfg, frRange.start, 2);				
    SetCAT2(ifo, cfg, frRange.stop,  2);				

    // apply CAT2
    if(TString(ifo).CompareTo(net->ifoName[cfg->nIFO-1])==0) {   // last ifo

      // get shifted merged dq cat0+cat1+cat2 list
      vector<waveSegment> cat2List=TB.readSegList(cfg->nDQF, cfg->DQF, CWB_CAT2);
      // store cat2List into network class
      net->segList=cat2List;

      // check if seg+cat2 data length is >= segTHR
      vector<waveSegment> detSegs(1);
      detSegs[0].start = xstart;
      detSegs[0].stop  = xstop;
      vector<waveSegment> detSegs_dq2;
      detSegs_dq2.push_back(detSegs[0]);
      detSegs_dq2 = TB.mergeSegLists(detSegs_dq2,cat2List);
      for(int i=0;i<(int)detSegs_dq2.size();i++) {
        cout << "detSegs_dq2[" << i << "]  GPS range : "
             << detSegs_dq2[i].start << "-" << detSegs_dq2[i].stop << endl;
      }
      cout << endl;
      double detSegs_ctime = TB.getTimeSegList(detSegs_dq2);
      cout << "live time after cat 2 : " << detSegs_ctime << " sec" << endl;
      if(detSegs_ctime<cfg->segTHR) {
        cout << "CWB_Plugin_Periodic_Frames.C : job segment live time after cat2 < segTHR="
             << cfg->segTHR << " sec, job terminated !!!" << endl;
        cout << endl << "To remove this check set segTHR=0 in the parameter file" << endl << endl;
        exit(1);
      }
    }
  }

  if(type==CWB_PLUGIN_STRAIN) {  

    cout << "CWB_Plugin_Periodic_Frames.C : Read Strain Frames ..." << endl;

    int ifoID=0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {ifoID=n;break;}

    int xstart   = (int)x->start();
    int xstop    = (int)x->stop();
    int segEdge  = int(cfg->segEdge);

    // open frame file
    CWB::frame fr;
    fr.open(cfg->frFiles[ifoID]);
    fr.setVerbose();
    fr.setRetryTime(cfg->frRetryTime);
    fr.setSRIndex(TMath::Nint(TMath::Log2(cfg->inRate))); // force data to be resampled to inRate
    int nfrFiles = fr.getNfiles();
    cout << "CWB_Plugin_Periodic_Frames.C - Strain " << " -> nfrFiles : " << nfrFiles << endl;

    // read frame file 
    frfile FRF = fr.getFrList(xstart-strainShift+segEdge, xstop-strainShift-segEdge, segEdge);
    fr.readFrames(FRF,const_cast<char*>(cfg->channelNamesRaw[ifoID]),*x);
    x->start(x->start()+strainShift);
    if(x->rate()!=cfg->inRate)
      {cout << "CWB_Plugin_Periodic_Frames.C - input rate from frame " << x->rate()
            << " do not match the one defined in config : " << cfg->inRate << endl;gSystem->Exit(1);}
  }

  return;
}

void SetCAT2(TString ifo, CWB::config* cfg, double gps, double window) {

  if(gps<0) return;

  if(window<=0) window = cfg->iwindow/2.;

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}

  for(int n=0; n<cfg->nIFO; n++) {

    if(ifo == cfg->DQF[n].ifo) {

      strcpy(cfg->DQF[cfg->nDQF].ifo, cfg->ifo[n]);
      sprintf(cfg->DQF[cfg->nDQF].file, "%s/%s_%s.gps_%d",cfg->tmp_dir,cfg->ifo[n],cfg->data_label,int(gps));
      cfg->DQF[cfg->nDQF].cat    = CWB_CAT2;
      cfg->DQF[cfg->nDQF].shift  = 0.;
      cfg->DQF[cfg->nDQF].invert = true;
      cfg->DQF[cfg->nDQF].c4     = true;
      cfg->nDQF++;

      cout << cfg->DQF[cfg->nDQF-1].file << endl;

      ofstream out;
      out.open(cfg->DQF[cfg->nDQF-1].file,ios::out);
      cout << "Write file : " << cfg->DQF[cfg->nDQF-1].file << endl;
      if (!out.good()) {cout << "Error Opening File : " << cfg->DQF[cfg->nDQF-1].file << endl;exit(1);}
      out.precision(14);
      int istart = int(gps)-window;
      int istop  = int(gps)+window;
      out << "1 " << istart << " " << istop << " " << 2*window << endl;
      out.close();
    }
  }
}

void SetTimeShiftCAT2(TString ifo, CWB::config* cfg, double shift) {

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}

  for(int n=0; n<(cfg->nDQF); n++) {
    if(ifo == cfg->DQF[n].ifo) cfg->DQF[n].shift  += shift;
  }
}

void ReadUserOptions(TString options) {

  // get plugin options 

  if(TString(options)!="") {

    //cout << "STRAIN options : " << options << endl;
    TObjArray* token = TString(options).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++) {

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();

      if(stok.Contains("strain_veto_time=")) {
        TString strain_veto_time=stok;
        strain_veto_time.Remove(0,strain_veto_time.Last('=')+1);
        if(strain_veto_time.IsFloat()) gOPT.veto_time=strain_veto_time.Atof();
      }

      if(stok.Contains("strain_veto_window=")) {
        TString strain_veto_window=stok;
        strain_veto_window.Remove(0,strain_veto_window.Last('=')+1);
        if(strain_veto_window.IsFloat()) gOPT.veto_window=strain_veto_window.Atof();
      }
    }
  }
}

void ResetUserOptions() {

    gOPT.veto_time    = STRAIN_VETO_TIME;
    gOPT.veto_window  = STRAIN_VETO_WINDOW;
} 

void PrintUserOptions() {

    cout.precision(14);
    cout << "-----------------------------------------"     << endl;
    cout << "STRAIN config options                    "     << endl;
    cout << "-----------------------------------------"     << endl << endl;
    cout << "STRAIN_VETO_TIME     " << gOPT.veto_time       << endl;
    cout << "STRAIN_VETO_WINDOW   " << gOPT.veto_window     << endl;
    cout << endl;
}

