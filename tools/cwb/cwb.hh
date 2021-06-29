/*
# Copyright (C) 2019 Gabriele Vedovato, Sergey Klimenko
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


/**********************************************************
 * Package:      CWB Class Library
 * File name:    CWB.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef CWB_HH
#define CWB_HH

#include "TObjString.h"
#include "TObjArray.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TTree.h"
#include "TFile.h"
#include "TNamed.h"
#include "TGlobal.h"
#include "TFileMerger.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <vector>

#include "wavecomplex.hh"
#include "wavearray.hh"
#include "wat.hh"

#include "History.hh"
#include "Toolbox.hh"
#include "config.hh"
#include "time.hh"
#include "ced.hh"
#include "frame.hh"

#include "livetime.hh"
#include "variability.hh"
#include "wavenoise.hh"

#include "CWB_Plugin.h"

// macro used by the pipeline to share parameters between pipeline and Plugin

#define IMPORT(TYPE,VAR) {      	                                        \
  char cmdline[128];                                                            \
  TGlobal* _global = (TGlobal*)gROOT->GetGlobal("gPOINTER",true);               \
  if(_global==NULL) sprintf(cmdline,"void* gPOINTER = (void*)&%s;",#VAR);       \
  else              sprintf(cmdline,"gPOINTER = (void*)&%s;",#VAR);             \
  gROOT->ProcessLine(cmdline);                                                  \
  TGlobal* _gvar = (TGlobal*)gROOT->GetGlobal(#VAR,true);                       \
  if(_gvar!=NULL) {                                                             \
    void* gPOINTER=NULL;                                                        \
    _global = (TGlobal*)gROOT->GetGlobal("gPOINTER",true);                      \
    memcpy((void*)&gPOINTER,(void*)_global->GetAddress(),sizeof(void*));        \
    VAR = *(TYPE*)gPOINTER;                                                     \
  }                                                                             \
}

#define EXPORT(TYPE,VAR,CMD) {      	                                        \
  TGlobal* _global = (TGlobal*)gROOT->GetGlobal(#VAR,true);                     \
  char cmdline[128];                                                            \
  if(_global==NULL) sprintf(cmdline,"%s %s;",#TYPE,CMD);                        \
  else              sprintf(cmdline,"%s;",CMD);                                 \
  gROOT->ProcessLine(cmdline);                                                  \
}

enum CWB_PLUGIN {
  CWB_PLUGIN_CONFIG 		= 0,
  CWB_PLUGIN_NETWORK 		= 1,
  CWB_PLUGIN_STRAIN    		= 2,
  CWB_PLUGIN_DATA    		= 2,
  CWB_PLUGIN_MDC     		= 3,
  CWB_PLUGIN_STRAIN_AND_MDC    	= 4,
  CWB_PLUGIN_DATA_MDC           = 4,
  CWB_PLUGIN_REGRESSION    	= 4,
  CWB_PLUGIN_IDATA_CONDITIONING = 4,
  CWB_PLUGIN_WHITE              = 5,
  CWB_PLUGIN_ODATA_CONDITIONING = 5,
  CWB_PLUGIN_ICOHERENCE     	= 6,
  CWB_PLUGIN_ISUPERCLUSTER     	= 7,
  CWB_PLUGIN_OLIKELIHOOD     	= 8,
  CWB_PLUGIN_INIT_JOB 		= 9,
  CWB_PLUGIN_CLOSE_JOB 		= 10,
  CWB_PLUGIN_RMDC     		= 11,
  CWB_PLUGIN_ILIKELIHOOD     	= 12,
  CWB_PLUGIN_OCOHERENCE     	= 13,
  CWB_PLUGIN_OREADDATA          = 14,
  CWB_PLUGIN_DATA_CONDITIONING  = 15,
  CWB_PLUGIN_XCOHERENCE     	= 16,
  CWB_PLUGIN_OSUPERCLUSTER     	= 17,
  CWB_PLUGIN_ECOHERENCE         = 18,
  CWB_PLUGIN_XLIKELIHOOD        = 19,
  CWB_PLUGIN_CED                = 20,
  CWB_PLUGIN_BCOHERENCE     	= 21
};

enum CWB_STAGE {
  CWB_STAGE_FULL 		= 0,
  CWB_STAGE_INIT 		= 1,
  CWB_STAGE_STRAIN 		= 2,
  CWB_STAGE_CSTRAIN	   	= 3,
  CWB_STAGE_COHERENCE     	= 4,
  CWB_STAGE_SUPERCLUSTER     	= 5,
  CWB_STAGE_LIKELIHOOD     	= 6,
  CWB_STAGE_SAVE     		= 7,
  CWB_STAGE_FINISH     		= 8
};

extern void CWB_Plugin(TFile* jfile, CWB::config*, network*, WSeries<double>*, TString, int);

class cwb : public TNamed {

public:

  cwb(CWB_STAGE jstage=CWB_STAGE_FULL); 
  cwb(TString fName, TString xName="", CWB_STAGE jstage=CWB_STAGE_FULL); 
  cwb(CWB::config cfg, CWB_STAGE jstage=CWB_STAGE_FULL); 
  virtual ~cwb();  

  virtual void run(int runID=0);

  static size_t GetProcInfo(bool mvirtual=true);

  // return the pointer to the config object
  CWB::config* GetConfig() {return &cfg;}

  // return the pointer to the network object
  network* GetNetwork() {return &NET;}

  // return the pointer to the history object
  CWB::History* GetHistory() {return history;}

  // return the final stage to process by the analysis
  CWB_STAGE GetStage() {return jstage;}

  // return the maximun stage value
  static int GetStageSize() {return CWB_STAGE_LIKELIHOOD;}

  static TString GetStageString(CWB_STAGE jstage);
  void SetupStage(CWB_STAGE jstage);

  int SetSkyMask(network* net, CWB::config* cfg, char* options, char skycoord, double skyres=-1);
  static void MakeSkyMask(skymap& SkyMask, double theta, double phi, double radius);

  char GetLagMode() {return lagMode[0];}
  TArrayC GetLagBuffer() {return lagBuffer;}

  void print() {PrintAnalysis(false);}  	// *MENU*

  // return the FRF frame files list
  vector<frfile> GetFrList(int ifoID=-1);
  vector<frfile> GetFrList(TString ifo);

  // Method used by TBrowser when double click the cwb object
  virtual void Browse(TBrowser *b) {print();}		 

  double GetSegBegin() {return Tb;}
  double GetSegEnd()   {return Te;}

  bool IsSingleDetector() {return singleDetector;}

  void PrintAnalysisInfo(CWB_STAGE stage, TString comment, TString info, 
                         bool out=true, bool log=true);

//protected:	// ROOT6 fix, commented out otherwise not visibles from cWB plugins

  CWB::config cfg;

  virtual void   InitNetwork();
  virtual void   InitNetwork(TString fName);
  virtual void   InitHistory();
  virtual double InitJob();
  virtual double InitJob(TString fName);

  // virtual base class ReadData method 
  virtual double ReadData(double mdcShift, int ifactor) {return 0.;}
  virtual double ReadData(TString fName);

  // virtual base class DataConditioning method 
  virtual void   DataConditioning(int ifactor) {return;}
  virtual void   DataConditioning(TString fName, int ifactor);

  // virtual base class Coherence method 
  virtual void   Coherence(int ifactor) {return;}
  virtual void   Coherence(TString fName);

  // virtual base class SuperCluster method 
  virtual void   SuperCluster(int ifactor) {return;}
  virtual void   SuperCluster(TString fName);

  // virtual base class Likelihood method 
  virtual bool   Likelihood(int ifactor, char* ced_dir, netevent* netburst = NULL, 
                            TTree* net_tree = NULL, char* outDump = NULL) {return false;}

  void   LoadPlugin(TMacro& plugin, TMacro& configPlugin);

  void   Init();
  void   PrintAnalysis(bool stageInfos=true);
  void   PrintElapsedTime(int job_elapsed_time, TString info);
  void   PrintStageInfo(CWB_STAGE stage, TString comment, 
                        bool out=true, bool log=true, TString fname="");
  TString GetStageInfo(CWB_STAGE stage, TString comment, TString fname="");
  TString GetAnalysisInfo(CWB_STAGE stage, TString comment, TString info);

  void   Exec(char* command, int maxtry=3, bool verbose=true);
  void   FileGarbageCollector(TString ifName, TString ofName="", 
                              vector<TString> delObjList = vector<TString>());

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// declarations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TStopwatch watchJob;				//! job benchmark 
  TStopwatch watchStage;			//! stage benchmark 

  TString iname; 				// input root data file

  CWB::Toolbox TB;				//! Toolbox 
  int  nIFO;					// number of detectors
  char ifo[NIFO_MAX][8];			// detector names array
  int  runID;					// job ID
  CWB_STAGE istage;				// initial stage			
  CWB_STAGE jstage;				// final stage

  // the first nIFO positions contains the ifo frFiles, nIFO+1 contains the sim frFile
  CWB::frame   fr[2*NIFO_MAX];			// array of frame objects used to read frame files	
  int          nfrFiles[2*NIFO_MAX];		// frame files list names 
  frfile       FRF[2*NIFO_MAX];	        	// frame file structures 

  // output root file variables
  TFile* froot;					//! output root file

  // temporary job file variables
  TFile* jfile;					//! job file object
  char jname[1024];				//  job file name

  WSeries<float> v[NIFO_MAX];      		//! noise variability
  detector* pD[NIFO_MAX];          		//! pointers to detectors
  WSeries<double>* pTF[NIFO_MAX];  		//! pointers to WSeries

  network NET;                     		//! network object

  injection* mdc;				//! injection object
  livetime live;				//! livetime object
  netevent* netburst;				//! netburst object
  variability wavevar;				//! variability object
  wavenoise noiserms;				//! wavenoise object

  CWB::History* history;			//! history object

  unsigned int jobfOptions;			//! used for the stage stuff
  bool singleDetector;				// true -> single detector mode

  size_t rateANA;				// analysis data rate

  double Tb;					// start segment gps time (segEdge excluded)
  double Te;					// stop segment gps time (segEdge excluded)
  double dT;       				// WB segment duration (segEdge excluded)

  double mTau;     				//! maximum time delay
  double dTau;     				//! time delay difference

  vector<waveSegment> detSegs;			// detector job segments
  vector<waveSegment> cat1List;			//! category 1 data quality list
  vector<waveSegment> cat2List;			//! category 2 data quality list
  int jobID;int slagID;int segID[20];
  float slagShift[20];				// used to store the super lags 
  size_t lags;					// number of lags
  char lagMode[1];				// used to store the the user lagFile mode
  TArrayC lagBuffer;				// used to store the contents of the user lagFile

  bool bplugin;					//! false/true -> defined/not-defined

  ClassDef(cwb,4)
};  

#endif
