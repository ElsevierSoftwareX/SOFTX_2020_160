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


#include "cwb.hh"
#include "mdc.hh"
#include "time.hh"
#include "Math/Polar3D.h"
#include "TMD5.h"

#define EXIT(ERR) gSystem->Exit(ERR)   // better exit handling for ROOT stuff

using namespace ROOT::Math;

ClassImp(cwb)

cwb::cwb(CWB_STAGE jstage) {
//
// Default Constructor
//
// used only for cwb class streaming 
//

  this->istage=CWB_STAGE_FULL;
  this->jstage=jstage;
  this->runID=0;
  this->singleDetector=false;
  this->lagBuffer.Set(0);
  SetupStage(jstage);
  Init();
}

cwb::cwb(TString fName, TString xName, CWB_STAGE jstage) {
//
// Constructor
//
// In this method the cwb configuration is loaded
// The default configuration is defined in the $CWB_PARAMETERS_FILE file
// The user can provide a customized configuration using the fName or xName parameters
//
// fName  : if fName ends with *.C the parameter is used as cwb user configuration file name 
//             The custom config must contains only the parameters which differ from the defaults
//             and must not be declared with types (Ex : fLow=64; instead of double fLow=64;)
//          if fName="" only the default cwb configuration file is used ($CWB_PARAMETERS_FILE)
//          if fName ends with *.root the parameter is used as job file name 
//             the configuration is read from the job file (the previous processed stage)
//             istage is initialized with the previous processed stage
// xName  : this parameter is used as auxiliary cwb configuration file when fName ends with *.root    
//          the xName file is used to change the configuration stored in the job file 
//
// jstage : is the final stage to be processed 
//

  this->runID=0;
  this->singleDetector=false;
  this->lagBuffer.Set(0);
  this->istage=CWB_STAGE_FULL;
  if(fName.EndsWith(".root")) {
    TFile* ifile = new TFile(fName);
    if(!ifile->IsOpen()) {
      cout << "cwb::cwb - Error opening root file : " << fName.Data() << endl;
      EXIT(1);
    }
    // read config object
    if(ifile->Get("config")!=NULL) {
      // read config object
      CWB::config* pcfg = (CWB::config*)ifile->Get("config");
      cfg = *pcfg;
      // read cwb object
      cwb* CWB = (cwb*)ifile->Get("cwb");
      if(CWB==NULL) {
        cout << "cwb::cwb - Error : cwb is not contained in root file " << fName.Data() << endl;
        EXIT(1); 
      }
      // read runID
      this->runID = CWB->runID;
      // read singleDetector
      this->singleDetector = CWB->singleDetector;
      // read initial stage
      this->istage = CWB->jstage;
      if(jstage<=istage) {
        cout << "cwb::cwb - Error : the stage(job file) " << GetStageString(istage).Data() << " >= " 
             << " input stage " << GetStageString(jstage).Data() << endl;EXIT(1);
        EXIT(1);
      }
      // read frame structures 
      for(int i=0;i<2*NIFO_MAX;i++) {
        this->nfrFiles[i] = CWB->nfrFiles[i];
        this->fr[i]       = CWB->fr[i];
        this->FRF[i]      = CWB->FRF[i];
      }
      // restore lag setting 
      lagBuffer = CWB->lagBuffer;		
      lagBuffer.Set(lagBuffer.GetSize()+1);
      lagBuffer[lagBuffer.GetSize()-1]=0;   	// set end of string 
      lagMode[0]=CWB->lagMode[0];lagMode[1]=0;	// restore lagMode
      // restore detSegs
      detSegs = CWB->detSegs;

      delete CWB;
      // set local environment
      if(gSystem->Getenv("HOME_WAT_FILTERS")==NULL) {
        cout << "cwb::cwb - Error : environment HOME_WAT_FILTERS is not defined!!!" << endl;EXIT(1);
      } else {
        strcpy(cfg.filter_dir,TString(gSystem->Getenv("HOME_WAT_FILTERS")).Data());
      }
      // export config to cint 
      cfg.Export();
      // set auxiliary configuration
      if(xName!="") {
        CWB::config icfg = cfg;	// cfg from input job file
        cfg.Import(xName); 	// cfg from input job file + auxiliary configuration
        if(cfg.simulation==4) {               		// initialize array factors
          if(int(cfg.factors[0])<=0) cfg.factors[0]=1;  // set 1 if <=0
          int ioffset = int(cfg.factors[0]);            // extract offset
          for(int i=ioffset;i<ioffset+cfg.nfactor;i++) cfg.factors[i]=i; // assign integer values
        }
        // check if auxiliary parameters are compatible with the previous stage config
        if(cfg.simulation!=icfg.simulation) {	// simulation
          cout << "cwb::cwb - Error : aux simulation " << cfg.simulation 
               << " != in cfg simulation " << icfg.simulation << endl; EXIT(1);}
        if(cfg.nfactor!=icfg.nfactor) {		// nfactor
          cout << "cwb::cwb - Error : aux nfactor " << cfg.nfactor 
               << " != in cfg nfactor " << icfg.nfactor << endl; EXIT(1);
        } else {				// factors
          for(int i=0;i<=cfg.nfactor;i++) if(cfg.factors[i]!=icfg.factors[i]) {
            cout << "cwb::cwb - Error : aux factors["<<i<<"]="<<cfg.factors[i] 
                 << " != in cfg factors["<<i<<"]=="<<icfg.factors[i]<<endl;EXIT(1);}
        } 
        if(cfg.l_low!=icfg.l_low) {		// l_low
          cout << "cwb::cwb - Error : aux l_low " << cfg.l_low 
               << " != in cfg l_low " << icfg.l_low << endl; EXIT(1);}
        if(cfg.l_high!=icfg.l_high) {		// l_high
          cout << "cwb::cwb - Error : aux l_high " << cfg.l_high 
               << " != in cfg l_high " << icfg.l_high << endl; EXIT(1);}
        if(cfg.fLow!=icfg.fLow) {		// fLow
          cout << "cwb::cwb - Error : aux fLow " << cfg.fLow 
               << " != in cfg fLow " << icfg.fLow << endl; EXIT(1);}
        if(cfg.fHigh!=icfg.fHigh) {		// fHigh
          cout << "cwb::cwb - Error : aux fHigh " << cfg.fHigh 
               << " != in cfg fHigh " << icfg.fHigh << endl; EXIT(1);}
        if(cfg.healpix!=icfg.healpix) {		// healpix
          cout << "cwb::cwb - Error : aux healpix " << cfg.healpix 
               << " != in cfg healpix " << icfg.healpix << endl; EXIT(1);}
      }
    } else {
      cout << "cwb::cwb - Error : config is not contained in root file " << fName.Data() << endl;
      EXIT(1); 
    }
    ifile->Close();
    iname=fName;
  } else if(fName.EndsWith(".C")) {
    if(gSystem->Getenv("CWB_PARAMETERS_FILE")==NULL) {
      cout << "cwb::cwb - Error : environment CWB_PARAMETERS_FILE is not defined!!!" << endl;
      EXIT(1);
    } else {
      // load default configuration $CWB_PARAMETERS_FILE
      cfg.Import("$CWB_PARAMETERS_FILE");
    }
    cfg.Import(fName);
    iname="";
  } else if(fName=="") {
    if(gSystem->Getenv("CWB_PARAMETERS_FILE")==NULL) {
      cout << "cwb::cwb - Error : environment CWB_PARAMETERS_FILE is not defined!!!" << endl;
      EXIT(1);
    } else {
      // load default configuration $CWB_PARAMETERS_FILE
      cfg.Import("$CWB_PARAMETERS_FILE");
    }
    iname="";
  } else {
    cout << "cwb::cwb - Error : bad input file extension [.C, .root] " << fName.Data() << endl;
    EXIT(1);
  }
  this->jstage=jstage;
  SetupStage(jstage);
  Init();
  if(fName!="") cfg.Check();	// check parameter's consistency
}

cwb::cwb(CWB::config cfg, CWB_STAGE jstage) {
//
// Constructor
//
// use the config object as configuration 
// 
// jstage : is the final stage to be processed 
//  

  this->cfg = cfg; 
  this->istage=CWB_STAGE_FULL;
  this->jstage=jstage;
  this->runID=0;
  this->singleDetector=false;
  this->lagBuffer.Set(0);
  SetupStage(jstage);
  Init();
  cfg.Check();			// check parameter's consistency
}

cwb::~cwb() {
//
// Destructor
//

  if(mdc!=NULL) delete mdc;
  //if(netburst!=NULL) delete netburst;  // commented because is already deleted by "delete froot"
  if(history!=NULL) delete history;
  if(jfile!=NULL) delete jfile;
  if(froot!=NULL) delete froot;
  if(singleDetector && pD[1]!=NULL) {pD[1]->IWFP.clear();pD[1]->RWFP.clear();}
  for(int n=0;n<NIFO_MAX;n++) if(pD[n]!=NULL) delete pD[n];
}

void
cwb::Init() {
//
// Reset/Initialize the pipeline parameters
//

  // There is a bug in ROOT > 5.32.04 & < 5.34.00
  // GarbageCollector don't remove the deleted objects from root file
  // Tree saved to root in the Coherence stage could be corrupted
#ifdef _USE_ROOT6 
  //cout << "cwb::Init - Warning !!! : This is a testing version for ROOT6" << endl;
#else
  if(gROOT->GetVersionInt()>53204 && gROOT->GetVersionInt()<53400) {
    cout << "cwb::Init - Error : cWB analysis don't works with ROOT version > 5.32.04 && version < 5.34.00" << endl;
    cout << "You are running version : ROOT " << gROOT->GetVersion() << endl << endl;
    EXIT(1);
  } 
#endif
  if((cfg.simulation==4)&&(jstage==CWB_STAGE_STRAIN)&&(cfg.nfactor>1)) {
    // NOTE : in the ReadData it is necessary to save MDC data for each factor
    cout << "cwb::Init - Error : stage STRAIN not implemented with simulation=4" << endl; 
    EXIT(1);
  }
 
  if((jstage!=CWB_STAGE_FULL)&&(jobfOptions&CWB_JOBF_SAVE_TRGFILE)) {
    cout<<"cwb::Init - Error : jobfOptions=CWB_JOBF_SAVE_TRGFILE is allowed only with CWB_STAGE_FULL!!!"<<endl;
    EXIT(1);
  }

  if(cfg.fResample>0) {                         // RESAMPLING
    rateANA=cfg.fResample>>cfg.levelR;
  } else {
    rateANA=cfg.inRate>>cfg.levelR;
  }

  // check if rate is an integer at resolution level = l_high
  if(rateANA-(1<<cfg.l_high)*(rateANA/(1<<cfg.l_high))!=0) {
    cout << "cwb::Init - Error : data rate : " << rateANA 
         << " is not a multiple of 2^l_high : " << (1<<cfg.l_high) << endl;
    EXIT(1);
  }

  froot=NULL;
  jfile=NULL;
  mdc=NULL;
  netburst=NULL;
  history=NULL;
  for(int n=0;n<NIFO_MAX;n++) pD[n]=NULL;
  return;
}

void
cwb::run(int runID) {
//
// The method used to start the analysis
//
// runID : is the job ID number, this is used in InitJob method to identify the 
//         the time range to be analyzed
//
// These are the main actions performed by this method
//
// - LoadPlugin 
// - InitNetwork
// - InitHistory
// - InitJob
// - Loop over Factors
//   - ReadData
//   - DataConditioning
//   - Coherence
//   - SuperCluster 
//   - Likelihood
//   - Save Recontructed Parameters
// - Save Job File (only for multi stage analysis)
//

  lags = 0;
  double factor = 1.0;          // strain factor
  int ioffset = 0; 	        // ifactor offset
  watchJob.Start();	        // start job benchmark
  watchStage.Start();	        // start stage benchmark

  bplugin = (TString(cfg.plugin.GetName())!="") ? true : false;     // user plugin

  unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files 

  double DATA_RATE = 0.;

  char tmpFile[1024];
  char outFile[1024];
  char endFile[1024];
  char outDump[1024];
  char endDump[1024];
  char out_CED[1024];
  char end_CED[1024];
  char command[1024];
  int  ecommand=0;
  FileStat_t fstemp;
  char cmd[1024]; 

  if(cfg.nIFO==0) {cout << "cwb::cwb - Error : no detector is presents in the configuration" << endl;EXIT(1);}

  this->nIFO=cfg.nIFO;
  for(int n=0;n<nIFO;n++) {
    if(strlen(cfg.ifo[n])>0) strcpy(this->ifo[n],cfg.ifo[n]);           // built in detector
    else                     strcpy(this->ifo[n],cfg.detParms[n].name); // user define detector
  }

  if(this->runID==0) { 			    // this->runID>0 for multistage analysis 
    if(runID>0) {
      this->runID=runID;
    } else {
      TString srunID = TString(gSystem->Getenv("CWB_JOBID"));
      this->runID = srunID.Atoi();
    }
  }

  sprintf(cfg.work_dir,gSystem->WorkingDirectory());
  sprintf(cfg.data_label,gSystem->BaseName(cfg.work_dir));

  cout<<"job ID   : "<<this->runID<<endl;
  cout<<"output   : "<<cfg.output_dir<<endl;
  cout<<"label    : "<<cfg.data_label<<endl;
  cout<<"nodedir  : "<<cfg.nodedir<<endl;
  cout<<"Pid      : "<<Pid<<endl;

  // create log, nodedir & output directories
  // Note : this step is necessary in multi stage analysis when pipeline
  //        start from an intermediate stage from a non structured working dir
  sprintf(cmd,"mkdir -p %s",cfg.log_dir);    gSystem->Exec(cmd);
  sprintf(cmd,"mkdir -p %s",cfg.nodedir);    gSystem->Exec(cmd);
  sprintf(cmd,"mkdir -p %s",cfg.output_dir); gSystem->Exec(cmd);

  // export to CINT istage,jstage (used by plugins)
  sprintf(cmd,"gISTAGE = (CWB_STAGE)%d;",istage); EXPORT(CWB_STAGE,gISTAGE,cmd)
  sprintf(cmd,"gJSTAGE = (CWB_STAGE)%d;",jstage); EXPORT(CWB_STAGE,gJSTAGE,cmd)
  // export to CINT ifactor (used by plugins)
  sprintf(cmd,"gIFACTOR = -1;"); EXPORT(int,gIFACTOR,cmd)	// init to gIFACTOR=-1

  // compile & load user plugin
  if(bplugin) LoadPlugin(cfg.plugin,cfg.configPlugin);

  if(bplugin) {CWB_Plugin(NULL,&cfg,&NET,NULL,"",CWB_PLUGIN_CONFIG);SetupStage(jstage);}

  // ---------------------------------------
  // init network
  // ---------------------------------------

  iname=="" ? InitNetwork() : InitNetwork(iname);

  PrintAnalysis();

  // ---------------------------------------
  // init history
  // ---------------------------------------

  InitHistory();

  // ---------------------------------------
  // init job
  // ---------------------------------------

  double mdcShift = iname=="" ? InitJob() : InitJob(iname);

  // temporary job file
  sprintf(jname,"%s/job_%d_%s_%d_%d.root",cfg.nodedir,int(Tb),cfg.data_label,this->runID,Pid);
  cout<<"temporary job file : " << jname<<endl;

  if(jstage==CWB_STAGE_INIT) goto JOB_END;

  // ---------------------------------------
  // read data (simulation!=4)
  // ---------------------------------------

  if(cfg.simulation!=4) {
    DATA_RATE = ReadData(mdcShift,0);
    if(jstage==CWB_STAGE_STRAIN) goto JOB_END;
  }

  // ---------------------------------------
  // if simulation!=0, loop on the injection strain factors
  // ---------------------------------------

  if(!cfg.simulation) cfg.nfactor = 1;

  ioffset = (cfg.simulation==4) ? int(cfg.factors[0]) : 0; // ifactor offset

  for(int ifactor=ioffset; ifactor<ioffset+cfg.nfactor; ifactor++) {

    int ceddir = 0; // flag if ced directory exists

    // export to CINT ifactor (used by plugins)
    sprintf(cmd,"gIFACTOR = %d;",ifactor); EXPORT(int,gIFACTOR,cmd) 

    // ---------------------------------------
    // declaration of output files
    // ---------------------------------------

    if(cfg.simulation) {
      factor=cfg.factors[ifactor];
      cout<<endl<< "---> Start processing factor["<<ifactor<< "]="<<factor<<endl<< endl;
      char sfactor[32];
      if(cfg.simulation==3) {
        if(factor<0)  sprintf(sfactor,"n%g",fabs(factor));
        if(factor==0) sprintf(sfactor,"z%g",factor);
        if(factor>0)  sprintf(sfactor,"p%g",factor);
      } else          sprintf(sfactor,"%g",factor);
      char sim_label[512];
      sprintf(sim_label,"%d_%d_%s_%s_job%d",int(Tb),int(dT),cfg.data_label,sfactor,this->runID);

      sprintf(outFile,"%s/wave_%s_%d.root",cfg.nodedir,sim_label,Pid);
      sprintf(endFile,"%s/wave_%s.root",cfg.output_dir,sim_label);
      sprintf(tmpFile,"%s/wave_%s_%d.root.tmp",cfg.nodedir,sim_label,Pid);
      sprintf(outDump,"%s/wave_%s_%d.txt",cfg.nodedir,sim_label,Pid);
      sprintf(endDump,"%s/wave_%s.txt",cfg.output_dir,sim_label);
      sprintf(out_CED,"%s/ced_%s_%d",cfg.nodedir,sim_label,Pid);
      sprintf(end_CED,"%s/ced_%s",cfg.output_dir,sim_label);

      if(!gSystem->GetPathInfo(endFile,fstemp)) {
        printf("The file %s already exists - skip\n",endFile);
        fflush(stdout);
        TFile rf(endFile);
        if(!rf.IsZombie()) continue;
      }
    }
    else {
      char prod_label[512];
      sprintf(prod_label,"%d_%d_%s_slag%d_lag%lu_%lu_job%d",
              int(Tb),int(dT),cfg.data_label,slagID,cfg.lagOff,cfg.lagSize,this->runID);

      sprintf(outFile,"%s/wave_%s_%d.root",cfg.nodedir,prod_label,Pid);
      sprintf(endFile,"%s/wave_%s.root",cfg.output_dir,prod_label);
      sprintf(tmpFile,"%s/wave_%s_%d.root.tmp",cfg.nodedir,prod_label,Pid);
      sprintf(outDump,"%s/wave_%s_%d.txt",cfg.nodedir,prod_label,Pid);
      sprintf(endDump,"%s/wave_%s.txt",cfg.output_dir,prod_label);
      sprintf(out_CED,"%s/ced_%s_%d",cfg.nodedir,prod_label,Pid);
      sprintf(end_CED,"%s/ced_%s",cfg.output_dir,prod_label);
    }

    // check if out_CED & end_CED already exist
    if(cfg.cedDump) {
      Long_t id,size=0,flags,mt;
      if (!gSystem->GetPathInfo(out_CED,&id,&size,&flags,&mt)) {
        cout << "cwb::run - Warning !!! - Dir \"" << out_CED << "\" already exist" << endl;
        EXIT(0);
      }
      if (!gSystem->GetPathInfo(end_CED,&id,&size,&flags,&mt)) {
        cout << "cwb::run - Warning !!! - Dir \"" << end_CED << "\" already exist" << endl;
        EXIT(0);
      }
    }

    // remove outDump file
    Long_t xid,xsize,xflags,xmt;
    int xestat = gSystem->GetPathInfo(outDump,&xid,&xsize,&xflags,&xmt);
    if (xestat==0) {
      sprintf(command,"/bin/rm %s",outDump);
      ecommand=gSystem->Exec(command);
      if(ecommand) {cout << "cwb::cwb - Warning -> " << command << endl;}
    }

    cout<<"output file on the node : "<<outFile<<endl;
    cout<<"final output file name  : "<<endFile<<endl;
    cout<<"temporary output file   : "<<tmpFile<<endl;

    // ---------------------------------------
    // read mdc (simulation==4)
    // ---------------------------------------

    if(cfg.simulation==4) {
      DATA_RATE = ReadData(mdcShift,ifactor);
      if(TString(cfg.analysis)=="2G" && jstage==CWB_STAGE_STRAIN) continue;  
    }

    // ---------------------------------------
    // data conditioning
    // ---------------------------------------

    DataConditioning(ifactor);
    if(TString(cfg.analysis)=="2G" && jstage==CWB_STAGE_CSTRAIN) continue;  

    // ---------------------------------------
    // create output root file
    // initialization of tree structures
    // create prod/sim/lag directories
    // ---------------------------------------

    froot = new TFile(tmpFile, "RECREATE");
    if(froot==NULL) {cout << "cwb::cwb - Error opening root file : " << tmpFile << endl;EXIT(1);}
    TTree* net_tree = !cfg.outPlugin ? netburst->setTree() : NULL; // outPlugin=true -> out tree is disabled
    TTree* live_tree= live.setTree();
    TTree* mdc_tree=NULL;
    TTree* var_tree=NULL;
    TTree* noise_tree=NULL;

    if(cfg.simulation) {
      mdc_tree = mdc->setTree();
    } else {
      var_tree = wavevar.setTree();
      noise_tree = noiserms.setTree();
    }

    // ---------------------------------------
    // start of the coherent search
    // ---------------------------------------

    size_t mlagSize=cfg.lagOff+cfg.lagSize;
    size_t mlagOff=cfg.lagOff;
    size_t mlagStep=cfg.mlagStep;

    if(mlagStep==0) mlagStep=cfg.lagSize; // if mlagStep=0 -> standard lag analysis

    for(size_t mlag=mlagOff;mlag<mlagSize;mlag+=mlagStep) { // multilag loop

      cfg.lagOff  = mlag;
      cfg.lagSize = cfg.lagOff+mlagStep<=mlagSize ? mlagStep : mlagSize-cfg.lagOff;
      if(cfg.lagSize==0) continue;

      cout << "lagSize : " << cfg.lagSize << " lagOff : " << cfg.lagOff << endl;

      std::vector<double> livTime;	 
      if(iname!="") livTime=NET.livTime;  		// save livTime 
      if(!cfg.simulation) {                          	// setup lags
        lags = NET.setTimeShifts(cfg.lagSize,cfg.lagStep,cfg.lagOff,cfg.lagMax,
                                 lagBuffer.GetArray(),lagMode,cfg.lagSite);
        cout<<"lag step: "<<cfg.lagStep<<endl;
        cout<<"number of time lags: "<<lags<<endl;
      }
      else if(!lags) lags = NET.setTimeShifts();
      if(iname!="") NET.livTime=livTime;		// restore livTime

      // ---------------------------------------
      // coherence analysis
      // ---------------------------------------

      Coherence(ifactor);
      if(TString(cfg.analysis)=="2G" && jstage==CWB_STAGE_COHERENCE) continue;  

      // ---------------------------------------
      // supercluster analysis
      // ---------------------------------------

      SuperCluster(ifactor);
      cout<<endl;
      if(TString(cfg.analysis)=="1G") cout<<"events in the buffer: "<<NET.events()<<"\n";
      if(TString(cfg.analysis)=="2G" && jstage==CWB_STAGE_SUPERCLUSTER) continue;  

      // ---------------------------------------
      // likelihood
      // ---------------------------------------

      ceddir=Likelihood(ifactor, out_CED, netburst, net_tree, outDump);
      cout<<"\nSearch done\n";
      if(!cfg.cedDump) cout<<"reconstructed events: "<<NET.events()<<"\n";
      if(cfg.simulation) NET.printwc(0);

      // ---------------------------------------
      // save data into root file
      // ---------------------------------------

      PrintStageInfo(CWB_STAGE_SAVE,"Data Save");
      froot->cd();

      live.output(live_tree,&NET,slagShift,detSegs);

      if(cfg.simulation) {
        double ofactor=0;
        if(cfg.simulation==4)      ofactor=-factor;
        else if(cfg.simulation==3) ofactor=-ifactor; 
        else                       ofactor=factor;
        if(TString(cfg.analysis)=="1G") {
          if(cfg.dump) netburst->dopen(outDump,const_cast<char*>("w"));
          netburst->output(net_tree,&NET,ofactor);
          if(cfg.dump) netburst->dclose();
        }
        mdc->output(mdc_tree,&NET,ofactor);
      } else {
        if(TString(cfg.analysis)=="1G") {
          if(cfg.dump) netburst->dopen(outDump,const_cast<char*>("w"));
          netburst->output(net_tree,&NET);
          if(cfg.dump) netburst->dclose();
        }
        for(int i=0; i<nIFO; i++) {
          if(cfg.outfOptions&CWB_OUTF_SAVE_VAR)  	// write var tree
            wavevar.output(var_tree,&v[i],i+1,cfg.segEdge);
          if(cfg.outfOptions&CWB_OUTF_SAVE_NOISE)  	// write noise tree
            noiserms.output(noise_tree,&pD[i]->nRMS,i+1,DATA_RATE/2);
        }
      }

      // save log info to final output root file
      PrintStageInfo(CWB_STAGE_FINISH,"Job Finished",false,true,endFile);
      history->AddLog( const_cast<char*>("FULL"), const_cast<char*>("STOP JOB"));
      history->Write("history");

      froot->Write();

    } // end mlag loop

    froot->Close();

    // ---------------------------------
    // move output files to end dirs
    // ---------------------------------

    if((TString(cfg.analysis)=="2G") && (jstage!=CWB_STAGE_FULL && jstage!=CWB_STAGE_LIKELIHOOD)) { 
      // delete temporary output file
      sprintf(command,"/bin/rm %s", tmpFile);
      ecommand=gSystem->Exec(command);
      if(ecommand) {cout << "cwb::cwb - Warning -> " << command << endl;}
      continue;
    } 
    sprintf(command,"/bin/mv %s %s", tmpFile, outFile);
    if(!cfg.cedDump || (cfg.cedDump && cfg.online)) Exec(command);
    if(cfg.cedDump && !cfg.online) sprintf(command,"/bin/rm %s",tmpFile);
    else                           sprintf(command,"/bin/mv %s %s",outFile,endFile);
    Exec(command);
    if(cfg.cedDump && !cfg.online) sprintf(command,"/bin/rm %s",outDump);
    else                           sprintf(command,"/bin/mv %s %s",outDump,endDump);
    if(cfg.dump) {
      // if file outDump exists then it is moved to endDump 
      if(!gSystem->GetPathInfo(outDump,&xid,&xsize,&xflags,&xmt)) Exec(command);
    }
    xestat = gSystem->GetPathInfo(end_CED,&xid,&xsize,&xflags,&xmt);
    if (xestat==0) {
      sprintf(command,"/bin/mv %s/* %s/.",out_CED,end_CED);
    } else {
      sprintf(command,"/bin/mv %s %s",out_CED,end_CED);
    }
    if(cfg.cedDump && ceddir && !(jobfOptions&CWB_JOBF_SAVE_CED)) Exec(command);
  }

  JOB_END:

  sprintf(cmd,"gIFACTOR = -1;"); EXPORT(int,gIFACTOR,cmd)	// set gIFACTOR=-1

  // -------------------------------------------------------------------
  // write cwb, history, config, network objects into temporary job file
  // -------------------------------------------------------------------

  if((jstage==CWB_STAGE_FULL)&&(jobfOptions&CWB_JOBF_SAVE_TRGFILE))  {
    // set options to save trigger file when full stage is used
    this->cwb::jstage = CWB_STAGE_SUPERCLUSTER;
    //this->cwb::jobfOptions = CWB_JOBF_SAVE_DISABLE;
    cfg.jobfOptions = CWB_JOBF_SAVE_DISABLE;
    // clear network skymap filled in Likelihood stage (save space)
    NET. nSensitivity=0; NET. nAlignment=0;   NET. nCorrelation=0;
    NET. nLikelihood=0;  NET. nNullEnergy=0;  NET. nPenalty=0;
    NET. nCorrEnergy=0;  NET. nNetIndex=0;    NET. nDisbalance=0;
    NET. nSkyStat=0;     NET. nEllipticity=0; NET. nPolarisation=0;
    NET. nProbability=0;
  }

  jfile = new TFile(jname,"UPDATE");
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb::cwb - Error opening root file : " << jname << endl;EXIT(1);}
  jfile->cd();
  if(jobfOptions&CWB_JOBF_SAVE_CONFIG)  {		// write config object
    TString tempName  = cfg.configPlugin.GetName();	// save temporary macro name
    TString tempTitle = cfg.configPlugin.GetTitle();	// save temporary macro title
    // get original macro name	
    TObjString* objn = cfg.configPlugin.GetLineWith("//#?CONFIG_PLUGIN_NAME "); 
    TString origName = objn ? objn->GetString() : ""; 
    cfg.configPlugin.SetName(origName.ReplaceAll("//#?CONFIG_PLUGIN_NAME ",""));
    // get original macro title
    TObjString* objt = cfg.configPlugin.GetLineWith("//#?CONFIG_PLUGIN_TITLE "); 
    TString origTitle = objt ? objt->GetString() : ""; 
    cfg.configPlugin.SetTitle(origTitle.ReplaceAll("//#?CONFIG_PLUGIN_TITLE ",""));
    cfg.Write("config"); 
    cfg.configPlugin.SetName(tempName); 		// restore temporary macro name
    cfg.configPlugin.SetTitle(tempTitle);		// restore temporary macro title
  }
  if(NET.pixeLHood.pWavelet->m_WaveType==WDMT) {
    // in 2G analysis the network::getNetworkPixels do "pixeLHood = *pTF;"
    // the Forwad operation after WDM::setTDFilter include in pTF the TD structures
    // since WDM and SymmObjArray streamers are incomplete because of its pointers
    // the read/save to/from root jfile is not a safe operation 
    // a workaround is to set to 0 the TD structures
    ((WDM<double>*)NET.pixeLHood.pWavelet)->T0.Resize(0);
    ((WDM<double>*)NET.pixeLHood.pWavelet)->Tx.Resize(0);
  }
  if(jobfOptions&CWB_JOBF_SAVE_NETWORK) NET.Write("network");  // write network object
  if(jobfOptions&CWB_JOBF_SAVE_JNET_MACRO) {                   // write cwb_jnet macro
    // check if cwb_jnet.C macro exists
    TString cwb_jnet_name = gSystem->ExpandPathName("$CWB_MACROS/cwb_jnet.C"); 
    TB.checkFile(cwb_jnet_name);
    TMacro cwb_jnet(cwb_jnet_name);
    cwb_jnet.Write("cwb_jnet");
  }
  PrintStageInfo(CWB_STAGE_FINISH,"Job Finished",false,true);
  history->AddLog( const_cast<char*>("FULL"), const_cast<char*>("STOP JOB"));
  if(jobfOptions&CWB_JOBF_SAVE_HISTORY) history->Write("history");
  if(jobfOptions&CWB_JOBF_SAVE_CWB) this->cwb::Write("cwb");
  if(bplugin) {
    wavearray<double> x;         // temporary time series
    for(int i=0; i<nIFO; i++) {
      x.rate(cfg.inRate); x.start(FRF[i].start); x.resize(int(x.rate()*(FRF[i].stop-FRF[i].start)));
      CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)&x,ifo[i],CWB_PLUGIN_CLOSE_JOB);
    }
    x.resize(0);
  }
  jfile->Close();

  // ---------------------------------
  // clean-up temporary job file
  // ---------------------------------

  if(jobfOptions) {
    // get rid of duplicated trees created by TFileMergers
    TFile *f = TFile::Open(jname,"UPDATE");
    TDirectory* d;
    for(int ifactor=0; ifactor<cfg.nfactor; ifactor++) {
      for(int j=0; j<(int)NET.nLag; j++) {
        netcluster* pwc = NET.getwc(j);
        if(pwc==NULL) continue; 
        int cycle = cfg.simulation ? ifactor : Long_t(pwc->shift);
        char trName[64];sprintf(trName,"clusters-cycle:%d;2",cycle);
        d = (TDirectory*)f->Get("coherence;1");
        if(d!=NULL) if(d->Get(trName)!=NULL) d->Delete(trName); 
        d = (TDirectory*)f->Get("supercluster;1");
        if(d!=NULL) if(d->Get(trName)!=NULL) d->Delete(trName); 
        d = (TDirectory*)f->Get("likelihood;1");
        if(d!=NULL) if(d->Get(trName)!=NULL) d->Delete(trName); 
      }
    }
    f->Close();

    // assign final file name according to the stage name
    char jlabel[512]; sprintf(jlabel,"%d_%d_%s",int(Tb),int(dT),cfg.data_label);
    char wlabel[512]; sprintf(wlabel,"wave_%s",jlabel);
    TString jLABEL=jlabel;
    if((jstage!=CWB_STAGE_FULL)&&(jstage!=CWB_STAGE_LIKELIHOOD)) 
      sprintf(endFile,"%s/wave_%s_job%d.root",cfg.output_dir,jlabel,this->runID);
    TString ojname = endFile;
    if(TString(cfg.analysis)=="1G")      ojname.ReplaceAll(wlabel,TString("job_")+jLABEL);
    else {			         // 2G 
      if(jstage==CWB_STAGE_FULL)         ojname.ReplaceAll(wlabel,TString("job_")+jLABEL); 
      if(jstage==CWB_STAGE_INIT)         ojname.ReplaceAll(wlabel,TString("init_")+jLABEL); 
      if(jstage==CWB_STAGE_STRAIN)       ojname.ReplaceAll(wlabel,TString("strain_")+jLABEL); 
      if(jstage==CWB_STAGE_CSTRAIN)      ojname.ReplaceAll(wlabel,TString("cstrain_")+jLABEL); 
      if(jstage==CWB_STAGE_COHERENCE)    ojname.ReplaceAll(wlabel,TString("coherence_")+jLABEL); 
      if(jstage==CWB_STAGE_SUPERCLUSTER) ojname.ReplaceAll(wlabel,TString("supercluster_")+jLABEL); 
      if(jstage==CWB_STAGE_LIKELIHOOD)   ojname.ReplaceAll(wlabel,TString("job_")+jLABEL); 
    }

    // if CWB_JOBF_SAVE_NODE then 2G to store the intermediate stage job files to nodedir
    // and creates a symbolic link to the final working dir output file name
    if((jobfOptions&CWB_JOBF_SAVE_NODE)&&(jstage!=CWB_STAGE_FULL)&&(jstage!=CWB_STAGE_LIKELIHOOD)) {
      TString olname = ojname;				// symbolic link path destination
      ojname.ReplaceAll(TString(cfg.output_dir)+TString("/"),TString(cfg.nodedir)+TString("/"));    
                                                        // symbolic link path origin
      if(ojname.BeginsWith("/")) {
        sprintf(command,"/bin/ln -sf %s %s",ojname.Data(),olname.Data());
      } else {
        sprintf(command,"/bin/ln -sf ../%s %s",ojname.Data(),olname.Data());
      }
      cout << command << endl;
      ecommand=gSystem->Exec(command);
      if(ecommand) {cout << "cwb::cwb - Warning -> " << command << endl;}
    }

    sprintf(endFile,"%s",ojname.Data());	// used by PrintStageInfo

    // move temporary job file to final position
    sprintf(command,"/bin/mv %s %s",jname,ojname.Data());
    ecommand=gSystem->Exec(command);
    if(ecommand) {cout << "cwb::cwb - Warning -> " << command << endl;}

  } else {

    // remove temporary job file 
    sprintf(command,"/bin/rm %s",jname);
    ecommand=gSystem->Exec(command);
    if(ecommand) {cout << "cwb::cwb - Warning -> " << command << endl;}
  }

  // remove temporary configPlugin file
  if(TString(cfg.configPlugin.GetTitle())!="") 
     gSystem->Exec(TString("rm ")+cfg.configPlugin.GetTitle());

  // ---------------------------------
  // print final infos
  // ---------------------------------

  cout << endl << endl;
  PrintStageInfo(CWB_STAGE_FINISH,"Job Finished",true,false,endFile);
  int job_data_size_sec = dT;
  double job_speed_factor = double(job_data_size_sec)/double(watchJob.RealTime());
  cout << endl;
  printf("Job Speed Factor - %2.2fX\n",job_speed_factor);
  cout << endl;

  watchJob.Stop();watchJob.Reset();

  // ONLINE : create empty file "finished" in the output directory
  if(cfg.online) {
    sprintf(command,"touch %s/finished",cfg.output_dir);
    gSystem->Exec(command);
  }

  return;
}

void
cwb::InitNetwork(TString fName) {
//
// Init the Network object : set detectors and network parameters
//
// read the network object from the job file (the previous processed stage)
//

  PrintStageInfo(CWB_STAGE_INIT,"cwb::InitNetwork from file");

  // read network object from job file
  TFile* ifile = new TFile(fName);
  if(ifile==NULL) {cout << "cwb::InitNetwork - Error opening root file : " << fName.Data() << endl;EXIT(1);}
  network* pnet = (network*)ifile->Get("network");
  if(pnet!=NULL) {
    // copy network object into local NET object
    NET = *pnet;	
    delete pnet;
  } else {
    cout << "cwb::InitNetwork - Error : net is not contained in root file " << fName.Data() << endl;
    EXIT(1); 
  }
  ifile->Close();
  // save detector pointers into local structures pD
  for(int n=0; n<nIFO; n++) pD[n] = NET.getifo(n);
  // set the original sampling rate
  for(int i=0; i<nIFO; i++) pD[i]->rate = cfg.fResample>0 ? cfg.fResample : cfg.inRate;

  // restore skymaps (to save space the network skymaps was not saved into job file)
  if(cfg.healpix) NET.setSkyMaps(int(cfg.healpix));
  else            NET.setSkyMaps(cfg.angle,cfg.Theta1,cfg.Theta2,cfg.Phi1,cfg.Phi2);
  NET.setAntenna();

  // restore network parameters
  NET.constraint(cfg.delta,cfg.gamma);
  NET.setDelay(cfg.refIFO);
  NET.Edge = cfg.segEdge;
  NET.netCC = cfg.netCC;
  NET.netRHO = cfg.netRHO;
  NET.EFEC = cfg.EFEC;
  NET.precision = cfg.precision;
  NET.nSky = cfg.nSky;
  NET.eDisbalance = cfg.eDisbalance;
  NET.setRunID(runID);
  NET.setAcore(cfg.Acore);
  NET.optim=cfg.optim;
  NET.pattern=cfg.pattern;

  // restore sky mask
  if(TString(cfg.analysis)=="1G") {
    if(cfg.mask>0) NET.setSkyMask(cfg.mask,cfg.skyMaskFile);
    else if(cfg.mask<0 && strlen(cfg.skyMaskFile)>0) 
           SetSkyMask(&NET,&cfg,cfg.skyMaskFile,'e'); 
  } else {
    if(strlen(cfg.skyMaskFile)>0) SetSkyMask(&NET,&cfg,cfg.skyMaskFile,'e'); 
  }
  if(strlen(cfg.skyMaskCCFile)>0) SetSkyMask(&NET,&cfg,cfg.skyMaskCCFile,'c'); 

  // restore detector names defined in the network class
  // there is a issue in network class (no TClass for char*)
  // detectors names are not saved properly
  NET.ifoName.clear();
  for(int n=0; n<nIFO; n++) NET.ifoName.push_back(pD[n]->Name);

  // execute user plugin
  if(bplugin) CWB_Plugin(NULL,&cfg,&NET,NULL,"",CWB_PLUGIN_NETWORK);

  // create injection (injected events) and netevent (reconstructed events) objects
  mdc = new injection(nIFO);
  netburst = new netevent(nIFO,cfg.Psave);

  // print system infos
  gSystem->Exec("/bin/date");
  gSystem->Exec("/bin/hostname");

  return;
}

void
cwb::InitNetwork() {
//
// Init the Network object : 
//
// - create object & add detector objects to network
// - set network parameters : search type, regulators, thresholds, sky map grid resolution
// - set the user sky mask : enable/disable the sky map tiles (earth coordinates)
// - set the user celestial sky mask : enable/disable the sky map tiles (celestial coordinates)
//

  PrintStageInfo(CWB_STAGE_INIT,"cwb::InitNetwork");

  // ---------------------------------------
  // Check single detector mode
  // ---------------------------------------

  // if nIFO=1 the analysis is done as a network of 2 equal detectors
  if(nIFO==1) {
    // ------> cwb in Sigle Detector Mode !!!
    cfg.SetSingleDetectorMode();
    nIFO = cfg.nIFO;
    strcpy(ifo[1],ifo[0]);
    singleDetector=true;
  } else singleDetector=false;

  // ---------------------------------------
  // init network
  // ---------------------------------------

  // when ifo!="" check if user detectors defined in detParams has been initialized 
  for(int i=0; i<nIFO; i++) {
    if(strlen(cfg.ifo[i])==0 && strlen(cfg.detParms[i].name)==0) {
      cout << "cwb::InitNetwork - Error : user detector name at position " 
           << i << " is not defined (detParms)" << endl;
      EXIT(1); 
    }
  }

  // create detector objects
  for(int i=0; i<nIFO; i++) {
    if(strlen(cfg.ifo[i])>0) pD[i] = new detector(cfg.ifo[i]);        // built in detector
    else                     pD[i] = new detector(cfg.detParms[i]);   // user define detector
  }
  // set the original sampling rate
  for(int i=0; i<nIFO; i++) pD[i]->rate = cfg.fResample>0 ? cfg.fResample : cfg.inRate;
  // add detector object to network
  for(int i=0; i<nIFO; i++) NET.add(pD[i]);

  // set network skymaps 
  if(cfg.healpix) NET.setSkyMaps(int(cfg.healpix));
  else            NET.setSkyMaps(cfg.angle,cfg.Theta1,cfg.Theta2,cfg.Phi1,cfg.Phi2);
  NET.setAntenna();

  // restore network parameters
  NET.constraint(cfg.delta,cfg.gamma);
  NET.setDelay(cfg.refIFO);
  NET.Edge = cfg.segEdge;
  NET.netCC = cfg.netCC;
  NET.netRHO = cfg.netRHO;
  NET.EFEC = cfg.EFEC;
  NET.precision = cfg.precision;
  NET.nSky = cfg.nSky;
  NET.eDisbalance = cfg.eDisbalance;
  NET.setRunID(runID);
  NET.setAcore(cfg.Acore);
  NET.optim=cfg.optim;
  NET.pattern=cfg.pattern;

  // set sky mask
  if(TString(cfg.analysis)=="1G") {
    if(cfg.mask>0) NET.setSkyMask(cfg.mask,cfg.skyMaskFile);
    else if(cfg.mask<0 && strlen(cfg.skyMaskFile)>0) 
           SetSkyMask(&NET,&cfg,cfg.skyMaskFile,'e'); 
  } else {
    if(strlen(cfg.skyMaskFile)>0) SetSkyMask(&NET,&cfg,cfg.skyMaskFile,'e'); 
  }
  if(strlen(cfg.skyMaskCCFile)>0) SetSkyMask(&NET,&cfg,cfg.skyMaskCCFile,'c'); 

  // execute user plugin
  if(bplugin) CWB_Plugin(NULL,&cfg,&NET,NULL,"",CWB_PLUGIN_NETWORK);

  // create injection (injected events) and netevent (reconstructed events) objects
  mdc = new injection(nIFO);
  netburst = new netevent(nIFO,cfg.Psave);

  // print system infos
  gSystem->Exec("/bin/date");
  gSystem->Exec("/bin/hostname");

  return;
}

void cwb::InitHistory() {
//
// Init History : it is used for bookkeping
//
// for each of these stages : 
// (FULL,INIT,STRAIN,CSTRAIN,COHERENCE,SUPERCLUSTER,LIKELIHOOD)
// the following informations are saved
//
// CWB_ENV        : cwb enviromental variable
// CWB_ENV_MD5    : the MD5 of the CWB_ENV string 
// WATVERSION     : WAT version
// XIFO           : XIFO (max detectors's number)
// GITVERSION     : git version
// GITBRANCH      : git branch
// WORKDIR        : working directory path 
// FRLIBVERSION   : framelib version
// ROOTVERSION    : ROOT version
// LALVERSION     : LAL version
// DATALABEL      : label used to tag the analysis 
// CMDLINE        : command line used to start the job
// ROOTLOGON      : the cwb rootlogon.C 
// ROOTLOGON_MD5  : the MD5 of the ROOTLOGON_MD5 string
// PARAMETERS     : the user_parameters.C file
// PARAMETERS_MD5 : the MD5 of the PARAMETERS_MD5 string
// CWB_CONFIG_URL    : cWB/config url
// CWB_CONFIG_PATH   : cWB/config path
// CWB_CONFIG_BRANCH : cWB/config branch
// CWB_CONFIG_TAG    : cWB/config tag
// CWB_CONFIG_HASH   : cWB/config hash code
// CWB_CONFIG_DIFF   : cWB/config difference (modified)
//

  PrintStageInfo(CWB_STAGE_INIT,"cwb::InitHistory");

  // check if configuration files exist
  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));

  // declare stages
  const char* STAGE_NAMES[7] = {"FULL","INIT","STRAIN","CSTRAIN","COHERENCE","SUPERCLUSTER","LIKELIHOOD"};
  // declare types
  const char* TYPE_NAMES[22] = {"CWB_ENV","WATVERSION","XIFO","GITVERSION","GITBRANCH","WORKDIR",
                                "FRLIBVERSION","ROOTVERSION","LALVERSION",
                                "DATALABEL","CMDLINE","ROOTLOGON","PARAMETERS",
                                "CWB_ENV_MD5","ROOTLOGON_MD5","PARAMETERS_MD5",
                                "CWB_CONFIG_URL","CWB_CONFIG_PATH","CWB_CONFIG_BRANCH","CWB_CONFIG_TAG","CWB_CONFIG_HASH","CWB_CONFIG_DIFF"};
  // create history object
  history = new CWB::History(const_cast<char**>(STAGE_NAMES), 7, const_cast<char**>(TYPE_NAMES), 22);
  history->SetSortOrder(InsertionOrder);

  // If any, add history from previous processed stage 
  if(iname!="") {
    TFile* ifile = new TFile(iname);
    if(ifile==NULL) {cout << "cwb::InitHistory - Error opening root file : " << iname.Data() << endl;EXIT(1);}
    // read network object
    if(ifile->Get("history")!=NULL) {
      // read history object
      CWB::History* hst = (CWB::History*)ifile->Get("history");
      TList* stageList = hst->GetStageNames();
      TList* typeList  = hst->GetTypeNames();
      for(int i=0;i<stageList->GetSize();i++) {
        TObjString* stageString = (TObjString*)stageList->At(i);
        for(int j=0;j<typeList->GetSize();j++) {
          TObjString* typeString = (TObjString*)typeList->At(j);
          char* histData = hst->GetHistory(const_cast<char*>(stageString->GetString().Data()),
                                           const_cast<char*>(typeString->GetString().Data()));
          if(histData!=NULL) {	// add previous history stage to current history
            history->AddHistory(const_cast<char*>(stageString->GetString().Data()),
                                const_cast<char*>(typeString->GetString().Data()),histData);
            delete histData;
          }
        }
      }
      delete stageList;
      delete typeList;
      // read logs 
      int log_size = hst->GetLogSize(const_cast<char*>("FULL"));
      for(int i=0;i<log_size;i++) {
        TString log = hst->GetLog(const_cast<char*>("FULL"),i);
        // Replace the final stage (STG:8) of the input job file with (STG:jstage)
        char stg_label[16];sprintf(stg_label,"STG:%d",jstage);
        if(log.Contains("STG:8")) log.ReplaceAll("STG:8",stg_label);
        // Add previous history logs to current history
        history->AddLog(const_cast<char*>("FULL"), const_cast<char*>(log.Data()));
      }
      delete hst;
    } else {
      cout << "cwb::InitHistory - Error : history is not contained in root file " << iname.Data() << endl;
      EXIT(1); 
    }
    ifile->Close();
  }

  // add current history stage
  TString jStageString = GetStageString(jstage).ReplaceAll("CWB_STAGE_","");
  char jStage[256];sprintf(jStage,jStageString.Data());

  // get cwb enviromental variable and save to history
  char* cwbBuffer = TB.getEnvCWB();
  if(cwbBuffer!=NULL) {
    history->AddHistory(jStage, const_cast<char*>("CWB_ENV"), cwbBuffer);
    TMD5 md5;md5.Update((UChar_t*)cwbBuffer,strlen(cwbBuffer));md5.Final();
    history->AddHistory(jStage, const_cast<char*>("CWB_ENV_MD5"), const_cast<char*>(md5.AsString()));
    delete [] cwbBuffer;
  }

  // save library versions
  CWB::mdc MDC;
  char framelib_version[32]; sprintf(framelib_version,"%f",FRAMELIB_VERSION);
  history->AddHistory(jStage, const_cast<char*>("WATVERSION"), watversion('s'));
  history->AddHistory(jStage, const_cast<char*>("XIFO"), watversion('i'));
  history->AddHistory(jStage, const_cast<char*>("GITVERSION"), watversion('r'));
  history->AddHistory(jStage, const_cast<char*>("GITBRANCH"), watversion('b'));
  history->AddHistory(jStage, const_cast<char*>("FRLIBVERSION"), framelib_version);
  history->AddHistory(jStage, const_cast<char*>("ROOTVERSION"), const_cast<char*>(gROOT->GetVersion()));
  history->AddHistory(jStage, const_cast<char*>("LALVERSION"), const_cast<char*>(GetLALVersion().Data()));

  // save work directory and data label
  history->AddHistory(jStage, const_cast<char*>("WORKDIR"), cfg.work_dir);
  history->AddHistory(jStage, const_cast<char*>("DATALABEL"), cfg.data_label);

  // save command line used to start the application
  char cmd_line[2048]="";
  int cmd_line_len=0;
  for(int i=0;i<gApplication->Argc();i++) cmd_line_len+=strlen(gApplication->Argv(i));
  if(cmd_line_len>2047)  
    {cout << "cwb::InitHistory - command line too long : " << cmd_line_len << endl;EXIT(1);}
  for(int i=0;i<gApplication->Argc();i++) sprintf(cmd_line,"%s %s",cmd_line,gApplication->Argv(i));
  history->AddHistory(jStage, const_cast<char*>("CMDLINE"), cmd_line);

  // save rootlogon script (ROOT initialization script)
  char* rootlogonBuffer = TB.readFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  if(rootlogonBuffer!=NULL) {
    history->AddHistory(jStage, const_cast<char*>("ROOTLOGON"), rootlogonBuffer);
    TMD5 md5;md5.Update((UChar_t*)rootlogonBuffer,strlen(rootlogonBuffer));md5.Final();
    history->AddHistory(jStage, const_cast<char*>("ROOTLOGON_MD5"), const_cast<char*>(md5.AsString()));
    delete [] rootlogonBuffer;
  }

  // save git cWB/config infos
  history->AddHistory(jStage, const_cast<char*>("CWB_CONFIG_URL"),    const_cast<char*>(GetGitInfos("url","$CWB_CONFIG").Data()));
  history->AddHistory(jStage, const_cast<char*>("CWB_CONFIG_PATH"),   const_cast<char*>(GetGitInfos("path","$CWB_CONFIG").Data()));
  history->AddHistory(jStage, const_cast<char*>("CWB_CONFIG_BRANCH"), const_cast<char*>(GetGitInfos("branch","$CWB_CONFIG").Data()));
  history->AddHistory(jStage, const_cast<char*>("CWB_CONFIG_TAG"),    const_cast<char*>(GetGitInfos("tag","$CWB_CONFIG").Data()));
  history->AddHistory(jStage, const_cast<char*>("CWB_CONFIG_HASH"),   const_cast<char*>(GetGitInfos("hash","$CWB_CONFIG").Data()));
  history->AddHistory(jStage, const_cast<char*>("CWB_CONFIG_DIFF"),   const_cast<char*>(GetGitInfos("diff","$CWB_CONFIG")!="" ? "M" : ""));

  // save configuration file
  char tmpFile[1024];
  unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files 
  sprintf(tmpFile,"%s/CWB_Config_%s_%d_job%d.XXXXXX",cfg.nodedir,cfg.data_label,Pid,this->runID);
  if(mkstemp(tmpFile)==-1) {       				// create temporary file
    cout << "cwb::InitHistory - mkstemp error in creating tmp file : " << tmpFile << endl;
    EXIT(1); 
  }
  char nodedir[1024];strcpy(nodedir,cfg.nodedir);		// clean nodedir before dump to history
  strcpy(cfg.nodedir,"");                                       // this makes config the same for all jobs	 
  TString tempTitle = cfg.configPlugin.GetTitle();		// save temporary macro title
  								// get original macro title
  TObjString* objt = cfg.configPlugin.GetLineWith("//#?CONFIG_PLUGIN_TITLE "); 
  TString origTitle = objt ? objt->GetString() : ""; 
  cfg.configPlugin.SetTitle(origTitle.ReplaceAll("//#?CONFIG_PLUGIN_TITLE ",""));
  cfg.Print(tmpFile);						// write config to temporary file
  cfg.configPlugin.SetTitle(tempTitle);				// restore temporary macro title
  strcpy(cfg.nodedir,nodedir);					// restore nodedir
  char* parametersBuffer = TB.readFile(tmpFile);		// read config from tmp file
  gSystem->Exec((TString("/bin/rm ")+TString(tmpFile)).Data()); // delete temporary file
  if(parametersBuffer!=NULL) {
    history->AddHistory(jStage, const_cast<char*>("PARAMETERS"), parametersBuffer);
    TMD5 md5;md5.Update((UChar_t*)parametersBuffer,strlen(parametersBuffer));md5.Final();
    history->AddHistory(jStage, const_cast<char*>("PARAMETERS_MD5"), const_cast<char*>(md5.AsString()));
    delete [] parametersBuffer;
  }

  // save START JOB to Log history
  history->AddLog(const_cast<char*>("FULL"), const_cast<char*>("START JOB"));

  return;
}

double cwb::InitJob(TString fName) {
//
// Init the Job 
//
// read the Job configuration parameters from the job file (the previous processed stage)
//

  PrintStageInfo(CWB_STAGE_INIT,"cwb::InitJob from file");

  double mdcShift=0.;

  // read cwb object from job file
  TFile* ifile = new TFile(fName);
  if(ifile==NULL||!ifile->IsOpen()) 
    {cout << "cwb::InitJob - Error opening root file : " << fName.Data() << endl;EXIT(1);}
  ifile->cd();
  cwb* CWB = (cwb*)ifile->Get("cwb");
  if(CWB==NULL) {
    cout << "cwb::InitJob - Error : cwb is not contained in root file " << fName.Data() << endl;
    EXIT(1); 
  }

  // restore cwb segment parameters from saved cwb object 
  jobID = CWB->jobID;
  Tb = CWB->Tb;
  dT = CWB->dT;
  Te = CWB->Te;
  slagID = CWB->slagID;
  for(int n=0;n<=(int)NET.ifoListSize();n++) segID[n]=CWB->segID[n];		// restore segID
  for(int n=0;n<=(int)NET.ifoListSize();n++) slagShift[n]=CWB->slagShift[n];	// restore slags
  netburst->setSLags(slagShift);	// set slags into netevent
  ifile->Close();

  // restore TFmap start,rate (size=0) , needed by likelihood to get the initial segment time
  for(int n=0;n<(int)NET.ifoListSize();n++) {
    NET.getifo(n)->getTFmap()->start(Tb-cfg.segEdge);
    NET.getifo(n)->getTFmap()->rate(rateANA);
  }

  // set a fake TFmap for the first detector, needed to restore setTimeShifts
  if(TString(cfg.analysis)=="2G" && jstage==CWB_STAGE_LIKELIHOOD) {
    NET.getifo(0)->getTFmap()->resize(rateANA*(dT+2*cfg.segEdge));
    for(int n=0;n<NET.ifoListSize();n++) {  // used for ced
      NET.getifo(n)->getTFmap()->setlow(cfg.fLow);
      NET.getifo(n)->getTFmap()->sethigh(cfg.fHigh);
    }
  }
   
  // init simulation structures for simulation=4 
  if(cfg.simulation==4) {						// initialize array factors
    // for simulation=4 factors are used as tags for the output root files 
    if(fmod(cfg.factors[0],1)) {					// must be integer
      cout<<"cwb::InitJob : when simulation=4 factors[0] is the offset and must be integer>=0"<<endl;
      EXIT(1);
    }
    if(cfg.factors[0]<0) {
      cout<<"cwb::InitJob : when simulation=4 factors[0] is the offset and must be integer>=0"<<endl;
      EXIT(1);
    }
    if(int(cfg.factors[0])<=0) cfg.factors[0]=1;			// set 1 if <=0
    int ioffset = int(cfg.factors[0]);					// extract offset 
    for(int i=ioffset;i<ioffset+cfg.nfactor;i++) cfg.factors[i]=i;   	// assign integer values
  } 

  return mdcShift;
}

double cwb::InitJob() {
//
// Init the Job structures
//
// Selection of the period to be analyzed by this job (job segment)
// - Read the data quality segment (DQ) files (CAT 0/1/2/4) : CWB::config::DFQ
//   CAT 0/1/4 are used to define where the data can be analyzed
//   CAT 2 are used to veto bad periods
//
// The job segment is computed according to the segment parameters : CWB::config::segXXX, CWB::config::slagXXX 
// Each job segment includes at beginning and at the end a scratch period which length is : CWB::config::segEdge
// There are 2 types of procedures used to select the job segment :
// - 1) the circular time-shifts (LAG) are performed within non shifted detector segments (same as 1G)
//      the length of the job segment is optimized to get the maximum livetime.
//      The minimum segment len is : CWB::config::segMLS , the maximum segment length is CWB::config::segLen
// - 2) the circular time-shifts (LAG) are performed within shifted detector segments (SLAG)
//      all job segments have the same length : CWB::config::segLen
//
// Initialization of frame file structures 
// - Read noise/MDC frame lists
//   noise/MDC can be read from files or from the plugin (generated 'On The Fly')
// - Read MDC log file (only for simulation)
//   Is the list of the injection parameters (inj time, hrss, source location, ...)
//

  PrintStageInfo(CWB_STAGE_INIT,"cwb::InitJob");

  double mdcShift=0.;

  vector<TString> ifos(nIFO);
  for(int n=0;n<nIFO;n++) ifos[n]=ifo[n];

  if((cfg.simulation)&&((cfg.slagSize>1)||(cfg.lagSize!=1))) {
    cout << "Error : slagSize<=1 & lagSize==1 in simulation mode !!!" << endl;EXIT(1);
  }

  // -----------------------------------------------------------------------
  // init cat1 segments
  //
  // job segments are computed using as segEdge the value segEdge+segOverlap
  // NB! it is used only to compute the job segments
  //
  // |x| = segEdge : |+| = segOverlap : |-| = segLen
  //                               
  // |xxx|-------------------|++++++|xxx|
  //                     |xxx|-------------------|++++++|xxx|
  // -----------------------------------------------------------------------

  if(cfg.slagSize>0) {   // SLAG Segments

    cout << endl << "START SLAG Init ..." << endl << endl;

    // Get zero lag merged dq cat0 & cat1 lists
    // The cat1List list is used only for slagSegs && slagJobList
    cat1List=TB.readSegList(cfg.nDQF, cfg.DQF, CWB_CAT1);

    // Get number/list of available super lag jobs
    // Compute the available segments with length=segLen contained between the interval [min,max] 
    // Where min,max are the minimum and macimum times in the cat1List list 
    // The start time of each segment is forced to be a multiple of segLen
    vector<waveSegment> slagJobList=TB.getSlagJobList(cat1List, int(cfg.segLen));
    int slagSegs=slagJobList.size();

    // Get super lag list : slagList
    // Is the list of available segment shifts according to the slag configuration parameters
    vector<slag> slagList=TB.getSlagList(nIFO, cfg.slagSize, slagSegs, cfg.slagOff,
                                         cfg.slagMin, cfg.slagMax, cfg.slagSite, cfg.slagFile);

    // Extract SLAG from the the slagList according the input runID
    slag SLAG=TB.getSlag(slagList,runID);
    if(SLAG.jobId!=runID) {cout << "jobID " << runID << " not found in the slag list !!!" << endl;EXIT(1);}
    slagID = SLAG.slagId[0];					// slag identification number
    jobID  = SLAG.jobId;					// for slag -> jobID=runID
    for(int n=0; n<nIFO; n++) segID[n]=SLAG.segId[n];		// segment shifts
    cout << "SuperLag=" << slagID << " jobID=" << jobID;
    for(int n=0; n<nIFO; n++) cout << " segID[" << ifo[n] << "]=" << segID[n];cout << endl;

    // Apply slag shifts to the DQF structures (data quality)
    for(int n=0;n<nIFO;n++) ifos[n]=ifo[n];
    TB.setSlagShifts(SLAG, ifos, cfg.segLen, cfg.nDQF, cfg.DQF);

    // Get merged dq cat0 & cat1 lists after slag shifts
    cat1List=TB.readSegList(cfg.nDQF, cfg.DQF, CWB_CAT1);

    // Extract ifo's slag segments ranges contained in the selected slag
    // Extract max length segment after the application of cat0+cat1 
    // Final segment length must be gt [cfg.segMLS+2*(cfg.segEdge+cfg.segOverlap)] & gt [cfg.segMLS]
    detSegs=TB.getSegList(SLAG, slagJobList, cfg.segLen, cfg.segMLS, cfg.segEdge+cfg.segOverlap, cat1List);
    cout << endl;
    cout << endl << "Segment type = SLAG" << endl;

    cout << endl << "END SLAG Init ..." << endl << endl;

  } else {     // Segments are not shifted (standard shifts)

    jobID  = runID;
    slagID = 0;
    for(int n=0;n<nIFO;n++) segID[n]=0;

    // Get merged dq cat0 & cat1 lists 
    cat1List=TB.readSegList(cfg.nDQF, cfg.DQF, CWB_CAT1);

    // Extract ifo's slag segments ranges
    // Extract max length segment after the application of cat0+cat1 
    // Final segment length must be gt [cfg.segMLS+2*(cfg.segEdge+cfg.segOverlap)] & gt [cfg.segMLS]
    detSegs=TB.getSegList(jobID, nIFO, cfg.segLen, cfg.segMLS, cfg.segEdge+cfg.segOverlap, cat1List);
    cout << endl;
    cout << "Segment type = LAG" << endl;
  }
  cout << "segLen       = " << cfg.segLen << " sec" << endl;
  cout << "segMLS       = " << cfg.segMLS << " sec" << endl;
  cout << "segOverlap   = " << cfg.segOverlap << " sec" << endl;
  cout << endl;

  // set slags in to the netevent class
  for(int n=0; n<nIFO; n++) slagShift[n]=(segID[n]-segID[0])*cfg.segLen;
  slagShift[nIFO]=slagID;
  netburst->setSLags(slagShift);

  // Check if seg length is compatible with WDM parity (only for 2G)
  // This condition is necessary to avoid mixing between odd
  // and even pixels when circular buffer is used for lag shift
  // The MRAcatalog distinguish odd and even pixels
  // If not compatible then the length is modified according the requirements 
  if(TString(cfg.analysis)=="2G") {
    int rate_min  = rateANA>>cfg.l_high;
    for(int i=0;i<(int)detSegs.size();i++) {
      double length = detSegs[i].stop - detSegs[i].start;
      if(int(length*rate_min+0.001)&1) detSegs[i].stop-=1;
    }
  }

  // add segOverlap to the detector's segments stop for this job
  for(int i=0;i<nIFO;i++) detSegs[i].stop+=cfg.segOverlap;
  // print detector's segments for this job
  cout.precision(14);
  for(int i=0;i<nIFO;i++) {
    cout << "detSegs_dq1[" << ifo[i] << "] GPS range : " 
         << detSegs[i].start << "-" << detSegs[i].stop << endl;
  }
  cout << endl;
  if(detSegs.size()==0) {cout << "no segments found for this job, job terminated !!!" << endl;EXIT(1);}

  // --------------------------------------------------------------------------
  // Fill Frame Files structures (FRF)
  // if (dataPlugin=false) noise data is read from frame files else from plugin 
  // if (mdcPlugin=false)    mdc data is read from frame files else from plugin 
  // --------------------------------------------------------------------------

  // init FRF structures for noise data 
  for(int n=0; n<nIFO; n++) {
    if(!cfg.dataPlugin) {		// noise data from frames
      if(TString(cfg.frFiles[n])=="") {
        cout << "cwb::InitJob - Error : noise frame files list name is not defined" << endl;
        EXIT(1); 
      }
      fr[n].setTimeRange(detSegs[n].start-cfg.dataShift[n]-2*cfg.segEdge,
                         detSegs[n].stop -cfg.dataShift[n]+2*cfg.segEdge);
      fr[n].open(cfg.frFiles[n]);	// read frame files list
      fr[n].setVerbose();
      fr[n].setRetryTime(cfg.frRetryTime); 
      fr[n].setSRIndex(TMath::Nint(TMath::Log2(cfg.inRate))); // force data to be resampled to inRate
      nfrFiles[n]=fr[n].getNfiles();	// number of frame files used for this job
      cout << ifo[n] << " -> nfrFiles : " << nfrFiles[n] << endl;
      FRF[n] = fr[n].getFrList(int(detSegs[n].start-cfg.dataShift[n]),
                               int(detSegs[n].stop-cfg.dataShift[n]), int(cfg.segEdge));
    } else {				// noise data from the user plugin
      FRF[n].start  = int(detSegs[n].start-cfg.dataShift[n]-cfg.segEdge);
      FRF[n].stop   = int(detSegs[n].stop -cfg.dataShift[n]+cfg.segEdge);
      FRF[n].length = FRF[n].stop-FRF[n].start;
    }
  }

  // init FRF structures for mdc data 
  if(cfg.simulation) { 
    for(int n=nIFO; n<2*nIFO; n++) {
      // mdc are declared in the array frFiles[nIFO,2*nIFO-1] 
      // or within a single frame file frFiles[nIFO]
      // in such case the frFiles[nIFO+1,2*nIFO-1] are a copy of frFiles[nIFO]
      if(TString(cfg.frFiles[n])=="" && n>nIFO) {FRF[n]=FRF[nIFO];continue;}
      if(!cfg.mdcPlugin) {		// mdc from frames
        // set frame file list
        if(TString(cfg.frFiles[n])=="") {
          cout << "cwb::InitJob - Error : MDC frame files list name is not defined" << endl;
          EXIT(1); 
        }
        if(cfg.mdc_shift.startMDC>=0) {	// cfg.mdc_shift.startMDC<0 -> automatically mdc time shift 	
          fr[n].setTimeRange(detSegs[n-nIFO].start-2*cfg.segEdge,detSegs[n-nIFO].stop+2*cfg.segEdge);
        }
        fr[n].open(cfg.frFiles[n]);
        fr[n].setVerbose();
        fr[n].setRetryTime(cfg.frRetryTime); 
        fr[n].setSRIndex(TMath::Nint(TMath::Log2(cfg.inRate))); // force data to be resampled to inRate
        nfrFiles[n]=fr[n].getNfiles();
        cout << "MDC " << " -> nfrFiles : " << nfrFiles[n] << endl;
        if(cfg.mdc_shift.startMDC<0) {  // the mdc shift is automatically selected
          // read mdc range from mdc frl files
          waveSegment mdc_range  = fr[n].getFrRange();
          cout << "mdc_range : " << mdc_range.start << " " << mdc_range.stop << endl;
          cfg.mdc_shift.startMDC=mdc_range.start;
          cfg.mdc_shift.stopMDC=mdc_range.stop;
        }
        mdcShift  = TB.getMDCShift(cfg.mdc_shift, detSegs[n-nIFO].start);
        cout << "mdcShift : " << mdcShift << endl;
        FRF[n] = fr[n].getFrList(int(detSegs[n-nIFO].start-mdcShift),
                                 int(detSegs[n-nIFO].stop-mdcShift), int(cfg.segEdge));
      } else {				// mdc from user plugin
        FRF[n]=FRF[n-nIFO];
      }
    }
   
    // init simulation structures for simulation=4 
    if(cfg.simulation==4) {						// initialize array factors
      // for simulation=4 factors are used as tags for the output root files 
      if(fmod(cfg.factors[0],1)) {					// must be integer
        cout<<"cwb::InitJob : when simulation=4 factors[0] is the offset and must be integer>=0"<<endl;
        EXIT(1);
      }
      if(cfg.factors[0]<0) {
        cout<<"cwb::InitJob : when simulation=4 factors[0] is the offset and must be integer>=0"<<endl;
        EXIT(1);
      }
      if(int(cfg.factors[0])<=0) cfg.factors[0]=1;			// set 1 if <=0
      int ioffset = int(cfg.factors[0]);				// extract offset 
      for(int i=ioffset;i<ioffset+cfg.nfactor;i++) cfg.factors[i]=i;   	// assign integer values
    } 
  }

  // ------------------------------
  // init cat2 segments
  // ------------------------------

  // get shifted merged dq cat0+cat1+cat2 list
  cat2List=TB.readSegList(cfg.nDQF, cfg.DQF, CWB_CAT2);
  // store cat2List into network class
  NET.segList=cat2List;

  // check if seg+cat2 data length is >= segTHR
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
  if(detSegs_ctime<cfg.segTHR) {
    cout << "cwb::InitJob : job segment live time after cat2 < segTHR=" 
         << cfg.segTHR << " sec, job terminated !!!" << endl;
    cout << endl << "To remove this check set segTHR=0 in the parameter file" << endl << endl;
    EXIT(1);
  }

  Tb=detSegs[0].start;			// start  seg time (excluded segEdge)
  Te=detSegs[0].stop;			// stop   seg time (excluded segEdge)
  dT = Te-Tb;                           // length seg time (excluded segEdge)

  // --------------------------------------------------------------------------------
  // read input mdc log file when data loaded from frame files (cfg.mdcPlugin==false)
  // --------------------------------------------------------------------------------
  if(cfg.simulation && !cfg.mdcPlugin) {             // read MDC log file
    int i=NET.readMDClog(cfg.injectionList,double(long(Tb))-mdcShift);
    printf("GPS: %16.6f saved,  injections: %d\n",double(long(Tb)),i);
    if(!cfg.mdcPlugin) TB.shiftBurstMDCLog(NET.mdcList, ifos, mdcShift);
    for(int i=0;i<(int)NET.mdcTime.size();i++) NET.mdcTime[i]+=mdcShift;

    // check if seg+cat2+inj data length is not zero
    vector<waveSegment> mdcSegs(NET.mdcTime.size());
    for(int k=0;k<(int)NET.mdcTime.size();k++) {
      mdcSegs[k].start=NET.mdcTime[k]-cfg.iwindow/2.;
      mdcSegs[k].stop=NET.mdcTime[k]+cfg.iwindow/2.;
    }
    vector<waveSegment> mdcSegs_dq2 = TB.mergeSegLists(detSegs_dq2,mdcSegs);
    double mdcSegs_ctime = TB.getTimeSegList(mdcSegs_dq2);
    cout << "live time in zero lag after cat2+inj : " << mdcSegs_ctime << endl;
    if(mdcSegs_ctime==0) {
      cout << "cwb::InitJob : job segment with zero cat2+inj "
           << "live time in zero lag, job terminated !!!" << endl;
      cout << "Warning : MDC data frames could be zero !!!" << endl;
      EXIT(1);
    }
  }

  // ------------------------------------------------------------
  // store the contents of the user lagFile to string (lagBuffer)
  // (used in multistages analysis to restore lags)
  // ------------------------------------------------------------
  lagBuffer.Set(0);
  lagMode[0]=cfg.lagMode[0];lagMode[1]=0;
  if(!cfg.simulation && cfg.lagFile!=NULL) {
    if(TString(cfg.lagMode)=="r") {	// read lags from file
      ifstream in;
      in.open(cfg.lagFile, ios::in);
      if(!in.good()) {
        cout << "cwb::InitJob : Error Opening File : " << cfg.lagFile << endl;
        EXIT(1);
      }
      // get length of file:
      in.seekg (0, in.end);
      int len = in.tellg();
      in.seekg (0, in.beg);
      lagBuffer.Set(len);
      in.read(lagBuffer.GetArray(),len);// store the contents of the file into the lagBuffer
      lagMode[0]='s';lagMode[1]=0;	// change lagMode : setTimeShifts read lags from string
      if(!in) {
        cout << "cwb::InitJob : Error Reading File : " << cfg.lagFile << endl;
        EXIT(1);
      }
      in.close();
      // set end of string
      lagBuffer.Set(lagBuffer.GetSize()+1);
      lagBuffer[lagBuffer.GetSize()-1]=0;   
    } else {				// write lags to file
      lagBuffer.Set(strlen(cfg.lagFile),cfg.lagFile);
      lagMode[0]='w';lagMode[1]=0;
    }
  } 

  // execute user plugin
  if(bplugin) {
    wavearray<double> x;         // temporary time series
    for(int i=0; i<nIFO; i++) {
      x.rate(cfg.inRate); x.start(FRF[i].start); x.resize(int(x.rate()*(FRF[i].stop-FRF[i].start)));
      CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)&x,ifo[i],CWB_PLUGIN_INIT_JOB);
    }
    x.resize(0);
  }

  return mdcShift;
}

double
cwb::ReadData(TString fName) {
//
// Read unconditioned Noise/MDC from input job file (the previous processed stage) 
// Data are saved into current job file
//

  PrintStageInfo(CWB_STAGE_STRAIN,"cwb::ReadData from file");

  // open in update mode the current job file (jfile)
  jfile = new TFile(jname,"UPDATE");
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb::ReadData - Error opening root file : " << jname << endl;EXIT(1);}
  TDirectory* cdstrain = (TDirectory*)jfile->Get("strain");
  if(cdstrain==NULL) cdstrain = jfile->mkdir("strain");

  // read  unconditioned Noise/MDC from input job file (ifile)
  // write unconditioned Noise/MDC to current job file (jfile)
  TFile* ifile = new TFile(fName);
  if(ifile==NULL) {cout << "cwb::ReadData - Error opening root file : " << fName << endl;EXIT(1);}

  WSeries<double>* pws=NULL;
  for(int i=0; i<nIFO; i++) {
    ifile->cd();
    // read strain from ifile
    pws = (WSeries<double>*)ifile->Get(TString("strain/")+pD[i]->Name);
    if(pws==NULL) 
      {cout << "cwb::ReadData - Error : data not present in file : " << fName << endl;EXIT(1);}
    // write strain to jfile
    pTF[i] = pD[i]->getTFmap();
    cdstrain->cd();pws->Write(ifo[i]);
    if(cfg.simulation) {
      ifile->cd();
      // read mdc from ifile
      pws = (WSeries<double>*)ifile->Get(TString("mdc/")+ifo[i]);
      if(pws==NULL) 
        {cout << "cwb::ReadData - Error : mdc not present in file : " << fName << endl;EXIT(1);}
      // write mdc to jfile
      TDirectory* cdmdc = (TDirectory*)jfile->Get("mdc");
      if(cdmdc==NULL) cdmdc = jfile->mkdir("mdc");
      cdmdc->cd();pws->Write(ifo[i]);
    }
    if(singleDetector) break;
  }

  // close ifile/jfile
  ifile->Close();
  jfile->Close();
  return pws->rate();
}

void
cwb::DataConditioning(TString fName, int ifactor) {
//
// Read conditioned/whitened Noise/MDC from input job file (the previous processed stage) 
// Data are saved into current job file
//

  PrintStageInfo(CWB_STAGE_CSTRAIN,"cwb::DataConditioning from file");

  char cdcstrain_name[32];sprintf(cdcstrain_name,"cstrain-f%d",ifactor);

  // open in update mode the current job file (jfile)
  jfile=NULL;
  TDirectory* jcdcstrain = NULL;
  if(jobfOptions&CWB_JOBF_SAVE_CSTRAIN) {
    jfile = new TFile(jname,"UPDATE");
    if(jfile==NULL||!jfile->IsOpen()) 
      {cout << "cwb::DataConditioning - Error opening root file : " << jname << endl;EXIT(1);}
    jcdcstrain=jfile->mkdir(cdcstrain_name);
  }

  // read  conditioned data from input job file (ifile)
  // write conditioned data to current job file (jfile)
  TFile* ifile = new TFile(fName);
  if(ifile==NULL) 
    {cout << "cwb::DataConditioning - Error opening root file : " << fName << endl;EXIT(1);}

  WSeries<double>* pws=NULL;
  for(int i=0; i<nIFO; i++) {
    ifile->cd();
    // read whitened data from ifile
    pws = (WSeries<double>*)ifile->Get(TString(cdcstrain_name)+TString("/")+pD[i]->Name);
    if(pws==NULL) 
      {cout << "cwb::DataConditioning - Error : data not present, job terminated!!!" << endl;EXIT(1);}
    // write whitened data to jfile
    if(jfile!=NULL) {jfile->cd();jcdcstrain->cd();pws->Write(ifo[i]);}
    pTF[i] = pD[i]->getTFmap();
    *pTF[i]=*pws;
    delete pws; 
    if(singleDetector) break;
  }

  // close ifile/jfile
  if(jfile!=NULL) jfile->Close();
  ifile->Close();
  return;
}

void
cwb::Coherence(TString fName) {
//
// Base class virtual Coherence method
// It is implemented in cwb1G/cwb2G
//

  PrintStageInfo(CWB_STAGE_COHERENCE,"cwb::Coherence from file");

  return;
}

void
cwb::SuperCluster(TString fName) {
//
// Base class virtual SuperCluster method
// It is implemented in cwb1G/cwb2G
//

  PrintStageInfo(CWB_STAGE_SUPERCLUSTER,"cwb::SuperCluster from file");

  return;
}

void
cwb::PrintAnalysis(bool stageInfos) {
//
// Print Analysis Setup
// input : stageInfos : if true stage infos are printed
//
// output format (Example) :
//
//                                jobID : 1
//
//                             Analysis : 2G
//                                stage : CSTRAIN
// 
//                            detectors : L1 H1 V1 
//                               search : un-modeled(r)-MRA-Standard
//                           simulation : 2                
// maximum time delay between detectors : 0.0264483        
//        maximum time delay difference : 0.0264483        
//                        HEALPix order : 7                
//         levelR, rateANA, fLow, fHigh : 4, 1024, 32, 512 
//                           bpp, Acore : 0.001, 1.5       
//                     netRHO and netCC : 5, 0.5           
//         clustering TFgap, Tgap, Fgap : 6.0, 0.0, 0.0
//                              pattern : 0
//             subnetcut subnet, subcut : 0.6, 0.33   
//               regulator delta, gamma : 0.05, 0.5   
//               precision csize, order : 50, 5
//

  if(stageInfos) PrintStageInfo(CWB_STAGE_INIT,"cwb::PrintAnalysis");

  cout << "                               jobID : " << runID << endl;
  cout << endl;
  cout << "                            Analysis : " << cfg.analysis << endl;
  cout << "                               stage : ";
  if(TString(cfg.analysis)=="1G") 	   cout << "FULL" << endl;
  else {				   // 2G
    if(jstage==CWB_STAGE_FULL)   	   cout << "FULL" << endl;
    if(jstage==CWB_STAGE_INIT) 	   	   cout << "INIT" << endl;
    if(jstage==CWB_STAGE_STRAIN) 	   cout << "STRAIN" << endl;
    if(jstage==CWB_STAGE_CSTRAIN)          cout << "CSTRAIN" << endl;
    if(jstage==CWB_STAGE_COHERENCE) 	   cout << "COHERENCE" << endl;
    if(jstage==CWB_STAGE_SUPERCLUSTER) 	   cout << "SUPERCLUSTER" << endl;
    if(jstage==CWB_STAGE_LIKELIHOOD) 	   cout << "LIKELIHOOD" << endl;
  }

  cout << "\n                           detectors : ";
  for(int i=0; i<cfg.nIFO; i++) cout<<ifo[i]<<" ";
  cout<<endl;

  cout << "                              search : ";
  if(TString(cfg.analysis)=="1G") {
    if(cfg.search=='E' || cfg.search=='E') cout<<"un-modeled (Energy)";
    if(cfg.search=='b' || cfg.search=='B') cout<<"un-modeled (single stream)";
    if(cfg.search=='r' || cfg.search=='R') cout<<"un-modeled (dual stream)";
    if(cfg.search=='i' || cfg.search=='I') cout<<"elliptical polarisation";
    if(cfg.search=='g' || cfg.search=='G') cout<<"circular polarisation";
    if(cfg.search=='s' || cfg.search=='S') cout<<"linear polarisation";
    cout << "(" << cfg.search << ")" << endl;
  } else if(TString(cfg.analysis)=="2G") {
    char _search = std::tolower(cfg.search);
    if(_search=='r')                       cout<<"un-modeled";
    if(_search=='i')                       cout<<"iota-wave";
    if(_search=='p')                       cout<<"psi-wave";
    if(_search=='l' || _search=='s')       cout<<"linear polarisation";
    if(_search=='c' || _search=='g')       cout<<"circular polarisation";
    if(_search=='e' || _search=='b')       cout<<"elliptical polarisation";
    cout << "(" << cfg.search << ")";
    if(cfg.optim) cout<<"-SRA";
    else          cout<<"-MRA";
    if(cfg.pattern==0)                     cout<<"-Standard";
    if((cfg.pattern!=0 && cfg.pattern<0))  cout<<"-Mixed";
    if((cfg.pattern!=0 && cfg.pattern>0))  cout<<"-Packet";
    cout << endl;
  } else {
    cout<<"\n unavailable analysis type !!!"<<endl;EXIT(1);
  }
  cout << "                          simulation : " << cfg.simulation << endl;;

  mTau=NET.getDelay(const_cast<char*>("MAX"));  // maximum time delay
  dTau=NET.getDelay(const_cast<char*>(""));     // time delay difference
    cout<<"maximum time delay between detectors : "<<mTau<<endl;
    cout<<"       maximum time delay difference : "<<dTau<<endl;
  if(cfg.healpix) {
    cout<<"                       HEALPix order : "<<cfg.healpix<<endl;
  } else {
    cout<<"           skymap angular resolution : "<<cfg.angle<<endl;
    cout<<"          skymap size in polar angle : "<<cfg.Theta1<<", "<<cfg.Theta2<<endl;
    cout<<"      skymap size in azimuthal angle : "<<cfg.Phi1<<", "<<cfg.Phi2<<endl;
  }
    cout<<"        levelR, rateANA, fLow, fHigh : "<<cfg.levelR<<", "<<rateANA<<", "
                                                   <<cfg.fLow<<", "<<cfg.fHigh<<endl;
    cout<<"                          bpp, Acore : "<<cfg.bpp<<", "<<cfg.Acore<<endl;
    cout<<"                    netRHO and netCC : "<<cfg.netRHO<<", "<<cfg.netCC<<endl;
  if(TString(cfg.analysis)=="1G") {
    cout<<"              regulator delta, gamma : "<<cfg.delta<<", "<<NET.gamma<<endl;
  } else {                                       // 2G
    cout<<"        clustering TFgap, Tgap, Fgap : "<<cfg.TFgap<<", "<<cfg.Tgap<<", "<<cfg.Fgap<<endl;

    cout<<"                             pattern : "<<cfg.pattern<<" -> ";
    int pattern = abs(cfg.pattern);
    if(pattern== 0) cout<<"('*' single pixel standard)"<<endl;
    if(pattern== 1) cout<<"('3|'  packet)"<<endl;
    if(pattern== 2) cout<<"('3-'  packet)"<<endl;
    if(pattern== 3) cout<<"('3/'  packet - chirp)"<<endl;
    if(pattern== 4) cout<<"('3\\'  packet - ringdown)"<<endl;
    if(pattern== 5) cout<<"('5/'  packet - chirp)"<<endl;
    if(pattern== 6) cout<<"('5\\'  packet - ringdown)"<<endl;
    if(pattern== 7) cout<<"('3+'  packet)"<<endl;
    if(pattern== 8) cout<<"('3x'  cross packet)"<<endl;
    if(pattern== 9) cout<<"('9p'  9-pixel square packet)"<<endl;
    if(pattern > 9) cout<<"('*' single pixel packet)"<<endl;

    cout<<"                      subnet, subcut : "<<cfg.subnet<<", "<<cfg.subcut<<endl;
    cout<<"              regulator delta, gamma : "<<cfg.delta<<", "<<NET.gamma<<endl;
  }
  if(cfg.precision!=0) {
    int precision = int(fabs(cfg.precision));
    int csize  = precision%65536;               // get number of pixels threshold per level
    int order = (precision-csize)/65536;        // get resampled order
    cout<<"              precision csize, order : "<<csize<<", "<<order<<endl;
  }
  if((TString(cfg.analysis)=="1G")&&(cfg.mask>0.)) {
    cout<<"              earth sky mask applied ; "<<cfg.skyMaskFile<<endl;
    cout<<"                                mask ; "<<cfg.mask<<endl;
  }
  if((TString(cfg.analysis)=="1G")&&(cfg.mask<0.)&&(strlen(cfg.skyMaskFile)>0)) 
    cout<<"              earth sky mask applied ; "<<cfg.skyMaskFile<<endl;
  if((TString(cfg.analysis)=="2G")&&(strlen(cfg.skyMaskFile)>0)) 
    cout<<"              earth sky mask applied ; "<<cfg.skyMaskFile<<endl;
  if(strlen(cfg.skyMaskCCFile)>0) 
    cout<<"          celestial sky mask applied ; "<<cfg.skyMaskCCFile<<endl;
  cout<<endl; 

  // print ifo infos
  //cout << " nIFO : " << nIFO << endl;
  //NET.print();
  //for(int n=0; n<nIFO; n++) pD[n]->print();

  return;
}

size_t
cwb::GetProcInfo(bool mvirtual) {
//
// Print to screen System Infos
// 
// mvirtual      : true/false -> {virtual memory)/(resident memory)
//
// output format (Example) :
//
// Mon Feb  3 17:34:12 UTC 2014              
// Memory             -  virtual : 596 (mb)  rss  : 315 (mb)
//
// 

   ProcInfo_t info;
   gSystem->GetProcInfo(&info); 

   cout << "Memory             -  virtual : " <<  info.fMemVirtual / 1024 
        << " (mb)  rss  : " <<  info.fMemResident / 1024 << " (mb)" << endl;
   return mvirtual ? info.fMemVirtual/1024 : info.fMemResident/1024;
}

void   
cwb::LoadPlugin(TMacro& plugin, TMacro& configPlugin) {
//
// Compile & load the plugin 
// Load the configPlugin 
//
// plugin.GetName()        : is the plugin macro name  
// plugin.GetTitle()       : is the macro path
//                           if ends with .C the plugin code is compiled
//                           if ends with .so then it is used as path for the pre-compiled lib
// configPlugin.GetName()  : is the configPlugin macro name  
//

  if(TString(plugin.GetName())=="") return;

  unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files 

  TString pluginPath = plugin.GetTitle();

  int error,check;

  // check if configPlugin exists and is readable  
  if(TString(configPlugin.GetName())!="") {
    // save original name/title into the macro file, is used in the InitHistory 
    configPlugin.AddLine(TString::Format("//#?CONFIG_PLUGIN_NAME %s",configPlugin.GetName()));
    configPlugin.AddLine(TString::Format("//#?CONFIG_PLUGIN_TITLE %s",configPlugin.GetTitle()));

    // set a unique name for macro, it is used by configPlugin.Exec() to create a temporary file 
    // extract file name
   
    char configPluginName[1024];                            
    if((gROOT->GetVersionInt()>=53400)&&(gROOT->GetVersionInt()<=53407)) {
      // starting from ROOT 5.34.00 to 5.34.07 TMacro::Exec() creates the temporary macro file under /tmp
      // the temporary configPlugin name must be defined without the directory name
      sprintf(configPluginName,"CWB_Plugin_Config_%s_%d_job%d.XXXXXX",
              cfg.data_label,Pid,this->runID);

    } else {
      sprintf(configPluginName,"%s/CWB_Plugin_Config_%s_%d_job%d.XXXXXX",
              cfg.nodedir,cfg.data_label,Pid,this->runID);
    }
    if(mkstemp(configPluginName)==-1) { // create temporary file, must be deleted only at the end of run
      cout << "cwb::LoadPlugin - Error : mkstemp error in creating tmp file : " << configPluginName << endl;
      EXIT(1); 
    }
    configPlugin.SetName("CWB_PluginConfig");
    configPlugin.SetTitle(configPluginName);
  }                                                

  // if the plugin path extension ends with 'so' the pre-compiled plugin is used
  check=0;
  if(pluginPath.EndsWith(".so")) {     	
    // check if the plugin compilation date is up to date  
    // the plugin compilation date must be > wat compilation date
    FileStat_t fStatus;
    int estat = gSystem->GetPathInfo(pluginPath.Data(),fStatus);
    if(!estat) {
      wat::Time plugin_comp_date(double(fStatus.fMtime)); // plugin compilation date (unix)
      plugin_comp_date.UnixToGps();		          // plugin compilation date (gps)
      wat::Time wat_comp_date(watversion('T'));	          // wat compilation date
      if(plugin_comp_date.GetGPS()<wat_comp_date.GetGPS()) {
        TString pluginSrc = pluginPath;
        pluginSrc.ReplaceAll("_C.so",".C");     
        cout << endl;        
        cout << "cwb::LoadPlugin : The plugin compilation date is not up to date !!! " << endl;
        cout << "                  plugin compilation date : " << plugin_comp_date.GetDateString() << endl;
        cout << "                  wat    compilation date : " << wat_comp_date.GetDateString() << endl;
        cout << endl;
        cout << "Recompile the plugin : " << endl;
        cout << "cwb_compile " << pluginSrc << endl << endl;
        EXIT(1); 
      }
    }                                                                             
    cout << endl << "cwb::LoadPlugin - Load pre-compiled plugin ..." << endl << endl; 
    check = gROOT->LoadMacro(pluginPath.Data(),&error);
    if(check==-1) {                                    
      cout << "cwb::LoadPlugin : Load pre-compiled Plugin Failed !!! " << endl;
      cout << "cwb::LoadPlugin : The plugin is compiled 'On-The-Fly'!!! " << endl;
    }                                                                             
  }

  // if the plugin path extension don't ends with 'so' or check!=0
  // then a copy is saved on tmp dir and the plugin is compiled 'On-The-Fly'
  if(!(pluginPath.EndsWith(".so"))||(check!=0)) { 
    cout << endl << "cwb::LoadPlugin - compile and load plugin ..." << endl << endl; 
    // create unique temporary file                     
    char tmpFile[1024]="";                            
    sprintf(tmpFile,"%s/CWB_Plugin_%s_%d_job%d.XXXXXX",cfg.nodedir,cfg.data_label,Pid,this->runID);
    if(mkstemp(tmpFile)==-1) {       	// create temporary file
      cout << "cwb::LoadPlugin - mkstemp error in creating tmp file : " << tmpFile << endl;
      EXIT(1); 
    }
    char pluginSrc[1024]="";                            
    sprintf(pluginSrc,"%s.C",tmpFile);	// must added .C to be compilable

    CWB::Toolbox::addCWBFlags(plugin,pluginSrc); // get macro + cWB compiler flags

    // compile & load plugin  
    int success = gSystem->CompileMacro(TString(pluginSrc).Data(),"f");
    gSystem->Exec(TString("/bin/rm ")+tmpFile);	// remove temporary file
    if(!success) {                                           
      cout << "cwb::LoadPlugin : Plugin Compilation Failed !!! " << endl;
      EXIT(1);                                                           
    }                                                                    

    // check if the compiled filename exists and is readable  
    TString pluginShr = pluginSrc;
    pluginShr.ReplaceAll(".C","_C.so");              
    check = gROOT->LoadMacro(pluginShr.Data(),&error,true);  
    if(check==-1) {                                           
      cout << "cwb::LoadPlugin : Plugin Compilation Failed !!! " << endl;
      EXIT(1);                                                           
    }                                                                    

    // removed temporary plugin files & shared object created by ACLiC
    TString pluginTmp = pluginSrc;                                  
    char tmpStr[1024];                                                 

    sprintf(tmpStr, "/bin/rm -f %s", pluginTmp.Data());                          
    system(tmpStr);                                                   
    pluginTmp.ReplaceAll(".C","_C.so");                      
    sprintf(tmpStr, "/bin/rm -f %s", pluginTmp.Data());           
    system(tmpStr);                                                   
    pluginTmp.ReplaceAll("_C.so","_C.d");                    
    sprintf(tmpStr, "/bin/rm -f %s", pluginTmp.Data());           
    system(tmpStr);                                                   
  }                                                                               

  return;
}        

void 
cwb::PrintElapsedTime(int job_elapsed_time, TString info) {
//
// convert job_elapsed_time to (hh:mm:ss) format and print it
//
// job_elapsed_time : time (seconds)
//
// info             : info string added to (hh:mm:ss)
//

  int job_elapsed_hour  = int(job_elapsed_time/3600);
  int job_elapsed_min   = int((job_elapsed_time-3600*job_elapsed_hour)/60);
  int job_elapsed_sec   = int(job_elapsed_time-3600*job_elapsed_hour-60*job_elapsed_min);
  char buf[1024];
  sprintf(buf,"%s %02d:%02d:%02d (hh:mm:ss)\n",info.Data(),job_elapsed_hour,job_elapsed_min,job_elapsed_sec);
  cout << buf;

  return;
}

void
cwb::PrintStageInfo(CWB_STAGE stage, TString comment, bool out, bool log, TString fname) {
// 
// Print stage infos to screen/history
//
// stage	: stage 
// comment 	: comment
// out		: true -> print to standard output
// log		: true -> save to history
// fname	: fname : job/output file name
//

  TString info = GetStageInfo(stage,comment,fname);
  // print to standard output
  if(out) cout << info.Data() << endl;
  // save to history 
  if(log&&history) history->AddLog(const_cast<char*>("FULL"), const_cast<char*>(info.Data()));
}

TString 
cwb::GetStageInfo(CWB_STAGE stage, TString comment, TString fname) {
// 
// return string filled using comment and internal analysis status parameters
//
// stage	: stage 
// comment 	: comment 
// fname	: fname : job/output file name
//
// The return string has this format (Example):
// 
// --------------------------------------------------------------------
// cwb::PrintAnalysis                                                  
// --------------------------------------------------------------------
// UTC                -  2014-02-03 21:17:17 UTC Mon                   
// Job   Elapsed Time -  00:00:09 (hh:mm:ss)                           
// Stage Elapsed Time -  00:00:01 (hh:mm:ss)                           
// Memory             -  virtual : 472 (mb)  rss  : 187 (mb)           
// Job File Size      -  0 (bytes) : 0 (kb) : 0 (mb)                   
// --------------------------------------------------------------------
// GPS:1075497453-JOB:1-STG:1-FCT:-1-JET:9-SET:1-MEM:472-JFS:0         
// --------------------------------------------------------------------
//
// where :
//
//   GPS : gps time (sec)
//   JOB : job number
//   STG : stage number (defined in cwb.hh)
//   FCT : factors index array
//   JET : Job Elapsed Time (sec)
//   SET : Stage Elapsed Time (sec)
//   MEM : Memory usage (MB)
//   JFS : Job File Size - temporary file (bytes)
//

  size_t GPS;		// GPS Time
  size_t JOB=runID;	// Job ID
  size_t STG=stage;	// stage ID
     int FCT;		// factor ID
  size_t JET;		// Job Elapsed Time	
  size_t SET;		// Stage Elapsed Time
  size_t MEM;		// Virtual Memory
  size_t JFS;		// Job File Size

  // get ifactor
  int gIFACTOR; IMPORT(int,gIFACTOR)
  FCT = gIFACTOR;

  // redirect cout to buffer and return contents
  std::stringstream buffer;
  std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());

  // get stage info
  cout << endl;
  cout << "--------------------------------------------------------------------" << endl;
  cout << comment.Data();  
  if(cfg.simulation&&(FCT>=0)) {
    cout << " - factor[" << FCT << "]=" << cfg.factors[FCT] << endl;
  } else cout<<endl;
  cout << "--------------------------------------------------------------------" << endl;
  cout << "UTC                -  "; cout.flush();
  wat::Time date("now"); cout << date.GetDateString() << endl; GPS=date.GetGPS();
  watchJob.Stop();
  JET = watchJob.RealTime();
  PrintElapsedTime(watchJob.RealTime(),"Job   Elapsed Time - ");
  watchJob.Continue();
  watchStage.Stop();
  SET = watchStage.RealTime();
  PrintElapsedTime(watchStage.RealTime(),"Stage Elapsed Time - ");
  watchStage.Reset();
  watchStage.Start();
  MEM=GetProcInfo();
  Long_t xid,xsize,xflags,xmt;
  TString xname = fname=="" ? jname : fname; 
  int xestat = gSystem->GetPathInfo(xname.Data(),&xid,&xsize,&xflags,&xmt);
  if(xestat==0) JFS=xsize; else JFS=0;
  cout << "Job File Size      -  " << JFS << " (bytes) : " 
       <<  int(JFS/1024.) << " (kb) : " <<  int(JFS/1024./1024) << " (mb)" << endl;
  cout << "--------------------------------------------------------------------" << endl;
  cout <<  "GPS:" << GPS << "-JOB:" << JOB << "-STG:" << STG << "-FCT:" << FCT  
       << "-JET:" << JET << "-SET:" << SET << "-MEM:" << MEM << "-JFS:" << JFS << endl;
  cout << "--------------------------------------------------------------------" << endl;
  cout << endl;

  // restore cout 
  std::cout.rdbuf(old);

  return buffer.str();
}

void
cwb::PrintAnalysisInfo(CWB_STAGE stage, TString comment, TString info, bool out, bool log) {
// 
// Print to screen/history the comment string
//
// stage	: stage 
// comment 	: comment to be showed/saved
// info         : extra user info 
// out		: true -> print to standard output
// log		: true -> save to history
//

  TString ainfo = GetAnalysisInfo(stage,comment,info);
  // print to standard output
  if(out) cout << ainfo.Data() << endl;
  // save to history 
  if(log&&history) history->AddLog(const_cast<char*>("FULL"), const_cast<char*>(ainfo.Data()));
}

TString 
cwb::GetAnalysisInfo(CWB_STAGE stage, TString comment, TString info) {
// 
// return string filled using comment/info and internal analysis status parameters
//
// stage	: stage 
// comment 	: comment 
// info         : extra user info 
//
// The return string has this format (Example):
// 
// --------------------------------------------------------------------
// GPS:1075497453-JOB:1-STG:1-FCT:-1-JET:9-SET:1-MEM:472-JFS:0         
// --------------------------------------------------------------------
//
// where :
//
//   GPS : gps time (sec)
//   JOB : job number
//   STG : stage number (defined in cwb.hh)
//   FCT : factors index array
//   JET : Job Elapsed Time (sec)
//   SET : Stage Elapsed Time (sec)
//   MEM : Memory usage (MB)
//   JFS : Job File Size - temporary file (bytes)
//

  size_t GPS;		// GPS Time
  size_t JOB=runID;	// Job ID
  size_t STG=stage;	// stage ID
     int FCT;		// factor ID

  // get ifactor
  TGlobal* global = gROOT->GetGlobal("gIFACTOR",true);
  if(global!=NULL) FCT = *(int*)global->GetAddress(); else FCT=-1;

  // redirect cout to buffer and return contents
  std::stringstream buffer;
  std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());

  // get stage info
  cout << endl;
  cout << "--------------------------------------------------------------------" << endl;
  cout << comment.Data();  
  if(cfg.simulation&&(FCT>=0)) {
    cout << " - factor[" << FCT << "]=" << cfg.factors[FCT] << endl;
  } else cout<<endl;
  cout << "--------------------------------------------------------------------" << endl;
  cout <<  "GPS:" << GPS << "-JOB:" << JOB << "-STG:" << STG << "-FCT:" << FCT << info << endl; 
  cout << "--------------------------------------------------------------------" << endl;
  cout << endl;

  // restore cout 
  std::cout.rdbuf(old);

  return buffer.str();
}

TString 
cwb::GetStageString(CWB_STAGE jstage) {
//
// return stage string format
//
// jstage : is the stage number defined the enum CWB_STAGE 
//

  switch(jstage) {
  case CWB_STAGE_FULL :
    return("CWB_STAGE_FULL"); 
    break;
  case CWB_STAGE_INIT :
    return("CWB_STAGE_INIT"); 
    break;
  case CWB_STAGE_STRAIN :
    return("CWB_STAGE_STRAIN"); 
    break;
  case CWB_STAGE_CSTRAIN :
    return("CWB_STAGE_CSTRAIN"); 
    break;
  case CWB_STAGE_COHERENCE :
    return("CWB_STAGE_COHERENCE"); 
    break;
  case CWB_STAGE_SUPERCLUSTER :
    return("CWB_STAGE_SUPERCLUSTER"); 
    break;
  case CWB_STAGE_LIKELIHOOD :
    return("CWB_STAGE_LIKELIHOOD"); 
    break;
  default :
    return(""); 
    break;
  }  
}

void 
cwb::SetupStage(CWB_STAGE jstage) {
//
// set save options for job file which is used to pass informations between stages 
//
// jstage : is the stage number defined the enum CWB_STAGE 
//

  jobfOptions = cfg.jobfOptions;

  switch(jstage) {
  case CWB_STAGE_INIT :
    jobfOptions |= CWB_JOBF_SAVE_CONFIG  |
                   CWB_JOBF_SAVE_NETWORK |
                   CWB_JOBF_SAVE_HISTORY |
                   CWB_JOBF_SAVE_CWB;
    break;
  case CWB_STAGE_STRAIN :
    jobfOptions |= CWB_JOBF_SAVE_CONFIG  |
                   CWB_JOBF_SAVE_NETWORK |
                   CWB_JOBF_SAVE_HISTORY |
                   CWB_JOBF_SAVE_CWB     |
                   CWB_JOBF_SAVE_STRAIN  |
                   CWB_JOBF_SAVE_MDC;
    break;
  case CWB_STAGE_CSTRAIN :
    jobfOptions |= CWB_JOBF_SAVE_CONFIG  |
                   CWB_JOBF_SAVE_NETWORK |
                   CWB_JOBF_SAVE_HISTORY |
                   CWB_JOBF_SAVE_CWB     |
                   CWB_JOBF_SAVE_CSTRAIN;
    break;
  case CWB_STAGE_COHERENCE :
    jobfOptions |= CWB_JOBF_SAVE_CONFIG  |
                   CWB_JOBF_SAVE_NETWORK |
                   CWB_JOBF_SAVE_HISTORY |
                   CWB_JOBF_SAVE_CWB     |
                   CWB_JOBF_SAVE_CSTRAIN |
                   CWB_JOBF_SAVE_COHERENCE;
    break;
  case CWB_STAGE_SUPERCLUSTER :
    jobfOptions |= CWB_JOBF_SAVE_CONFIG  |
                   CWB_JOBF_SAVE_NETWORK |
                   CWB_JOBF_SAVE_HISTORY |
                   CWB_JOBF_SAVE_CWB     |
                   CWB_JOBF_SAVE_SUPERCLUSTER;
    break;
  case CWB_STAGE_LIKELIHOOD :
    jobfOptions |= CWB_JOBF_SAVE_DISABLE;
    break;
  case CWB_STAGE_FULL :
    // if option CWB_JOBF_SAVE_TRGFILE and jstage=CWB_STAGE_FULL then 
    // all necessary objects are saved in the trigger file
    if(!(jobfOptions&CWB_JOBF_SAVE_TRGFILE)) break;
    jobfOptions |= CWB_JOBF_SAVE_CONFIG  |
                   CWB_JOBF_SAVE_NETWORK |
                   CWB_JOBF_SAVE_HISTORY |
                   CWB_JOBF_SAVE_CWB     |
                   CWB_JOBF_SAVE_SPARSE  |
                   CWB_JOBF_SAVE_SUPERCLUSTER;
    break;
  default :
    break;
  }  

  return;
}

void
cwb::Exec(char* command, int maxtry, bool verbose) {
//
// Execute System Command 
//
// command : system command 
// maxtry  : number of retry before to exit with error (default=3)
// verbose : show the retry infos (default=true)
//

  int ntry=0;
  int ecommand=1;
  while(ecommand&&(ntry<maxtry)) {
    ecommand=gSystem->Exec(command);
    if(ecommand) gSystem->Sleep(int(gRandom->Uniform(10000,30000)));  // wait [10,30] sec
    ntry++;
    if(verbose) {
      if(ecommand) { 
        cout << command << endl;
        cout << "cwb::Exec - NTRY " << ntry << endl;
        cout.flush(); 
      }
    }
  }
  if(ecommand) {cout << "cwb::Exec - Error -> " << command << endl;EXIT(1);}
  return;
}

void
cwb::FileGarbageCollector(TString ifName, TString ofName, vector<TString> delObjList) {
//                                                                                      
// for ROOT < 5.34.05                                                                   
// the objects marked as deleted do not free the disk space                             
// this method is used to free the deleted object space contained in the root file      
//                                                                                      
// for ROOT >= 5.34.05                                                                  
// All objected contained in ifName and are not in the delObjList are copied to ofName 
// This method is used to cleanup the ROOT file ifName                                  
//                                                                                      
// ifName      : input file name                                                         
// ofName      : output file name (cleaned) - if ofName="" -> ofName=ifName                                            
// delObjList  : list of object which must removed from the ifName (only for ROOT >= 5.34.05)
//                                                                                           

  bool replace = false;
  if(ofName=="") {     
    ofName = ifName;    
    ofName.ReplaceAll(".root","_tmp.root");
    replace = true;                        
  }                                        

  if(delObjList.size()==0) {
    if(!replace) {
      char command[1024];                                      
      sprintf(command,"/bin/mv %s %s",ofName.Data(),ifName.Data());
      gSystem->Exec(command);                                      
    } else return;
  }

  if(gROOT->GetVersionInt()<53405) {
    // delete the objects contained in the delObjList
    TFile* kfile = new TFile(ifName, "UPDATE");       
    for(int i=0;i<delObjList.size();i++) kfile->Delete(delObjList[i]+";1");
    kfile->Close();                                                   
  }                                                                   

  gErrorIgnoreLevel=kBreak;     // disable root kInfo & kWarning & kError messages
  TFileMerger M(false);                                                           
  M.AddFile(ifName,false);                                                        
  M.OutputFile(ofName);                                                           
  if(gROOT->GetVersionInt()<53400) {                                              
    // TFileMerger free the deleted object space                                  
    if(!M.Merge()) {                                                              
      cout << "cwb::FileGarbageCollector : Error - Merge failed !!!" << endl;EXIT(1);
    }                                                                                
  }                                                                                  
  if(gROOT->GetVersionInt()>=53400 && gROOT->GetVersionInt()<53405) {                
    // starting from ROOT 5.34.00 TFileMerger::Merge() merge also the deleted object 
    // The following is a workaround to fix this issue                               
    if(!M.PartialMerge(TFileMerger::kAllIncremental | TFileMerger::kRegular)) {      
       cout << "cwb::FileGarbageCollector : Error - Merge failed !!!" << endl;EXIT(1);
    }                                                                                 
  }                                                                                   
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,5)
  if(gROOT->GetVersionInt()>=53405) {                                                 
    // add files which must not be copied from ifName to ofName                       
    for(int i=0;i<delObjList.size();i++) M.AddObjectNames(delObjList[i]);           
    // Must add new merging flag on top of the the default ones                       
    Int_t default_mode = TFileMerger::kAll | TFileMerger::kIncremental;               
    Int_t mode = default_mode | TFileMerger::kSkipListed;                             
    if(!M.PartialMerge(mode)) {                                                       
      cout << "cwb::FileGarbageCollector : Error - Merge failed !!!" << endl;EXIT(1); 
    }                                                                                 
    M.Reset();                                                                        
  }                                                                                   
#endif
  gErrorIgnoreLevel=kUnset;     // re-enable root kInfo & kWarning messages           

  if(replace) {                 // replace ifName with ofName
    char command[1024];                                      
    sprintf(command,"/bin/mv %s %s",ofName.Data(),ifName.Data());
    gSystem->Exec(command);                                      
  }                                                              

  return;
}        

void
cwb::MakeSkyMask(skymap& SkyMask, double theta, double phi, double radius) {
//
// make a sky mask 
// is used to define what are the sky locations that are analyzed
// 
// theta,phi,radius : used to define the Celestial SkyMap
//                    define a circle centered in (theta,phi) and radius=radius 
//                    theta : [-90,90], phi : [0,360], radius : degrees
//
// SkyMask          : output sky celestial mask
//                    inside the circle is filled with 1 otherwise with 0
//

  int L = SkyMask.size();
  int healpix = SkyMask.getOrder();
  // check input parameters
  if(fabs(theta)>90 || (phi<0 || phi>360) || radius<=0 || L<=0) {
    cout << "cwb::MakeSkyMask : wrong input parameters !!! " << endl;
    if(fabs(theta)>90)   cout << theta << " theta must be in the range [-90,90]" << endl;      
    if(phi<0 || phi>360) cout << phi << " phi must be in the range [0,360]" << endl;
    if(radius<=0)        cout << radius << " radius must be > 0" << endl;
    if(L<=0)             cout << L << " SkyMask size must be > 0" << endl;
    EXIT(1);                                                                    
  }

  if (!gROOT->GetClass("Polar3DVector")) gSystem->Load("libMathCore");

  // compute the minimun available radius
  // must be greater than the side of a pixel  
  if(healpix) {
    int npix = 12*(int)pow(4.,(double)healpix);
    double sphere_solid_angle = 4*TMath::Pi()*pow(180./TMath::Pi(),2.);
    double skyres = sphere_solid_angle/npix;
    if(radius < sqrt(skyres)) radius = sqrt(skyres);
  } else {
    if(radius < SkyMask.sms) radius = SkyMask.sms;
  }

  double ph,th;
  GeographicToCwb(phi,theta,ph,th);	// Geographic to CWB coordinates

  Polar3DVector ov1(1, th*TMath::Pi()/180, ph*TMath::Pi()/180);

  int nset=0;
  for (int l=0;l<L;l++) {
    double phi = SkyMask.getPhi(l);
    double theta = SkyMask.getTheta(l);

    Polar3DVector ov2(1, theta*TMath::Pi()/180, phi*TMath::Pi()/180 );
    double Dot = ov1.Dot(ov2);
    double dOmega = 180.*TMath::ACos(Dot)/TMath::Pi();
    //cout << "dOmega : " << dOmega << endl;

    if(dOmega<=radius) {SkyMask.set(l,1);nset++;} else SkyMask.set(l,0);
  }

  if(!nset) {                                    
    cout << "cwb::MakeSkyMask : no sky positions setted !!! " << endl;
    cout << "check input mask parameters : theta = " 
         << theta << " phi = " << phi << " radius : " << radius << endl << endl; 
    EXIT(1);                                                                    
  }                                                                             

  return;
}

int
cwb::SetSkyMask(network* net, CWB::config* cfg, char* options, char skycoord, double skyres) {
//
// set earth/celestial sky mask 
// is used to define what are the sky locations that are analyzed
// by default all sky is used
// 
// net           : pointer to the network object
// cfg           : pointer to the cwb config object
// options       : used to define the earth/celestial SkyMap
//                 option A :
//                 skycoord='e'
//                   --theta THETA --phi PHI --radius RADIUS
//                   define a circle centered in (THETA,PHI) and radius=RADIUS 
//                   THETA : [-90,90], PHI : [0,360], RADIUS : degrees
//                 skycoord='c'
//                   --theta DEC --phi RA --radius RADIUS
//                   define a circle centered in (DEC,RA) and radius=RADIUS 
//                   DEC : [-90,90], RA : [0,360], RADIUS : degrees
//                 option B:
//                 file name
//                 format : two columns ascii file -> [sky_index	value]
//                 sky_index : is the sky grid index
//                 value     : if !=0 the index sky location is used for the analysis
// skycoord      : sky coordinates : 'e'=earth, 'c'=celestial
// skyres        : sky resolution : def=-1 -> use the value defined in parameters (angle,healpix)
// 
// return 0 if successful 

  if(skycoord!='e' && skycoord!='c') {
    cout << "cwb::SetSkyMask - Error : wrong input sky coordinates "
         << " must be 'e'/'c' earth/celestial" << endl;;
    EXIT(1);
  }

  if(strlen(options)>0) {              			  // set SkyMask
    if(!TString(options).Contains("--")) { // input parameter is the skyMask file
      if(skyres>=0) return 1;
      int ret = net->setSkyMask(options,skycoord);
      if(ret==0) {
        cout << endl;
        cout << "cwb::SetSkyMask - Error : skyMask file"  
             << " not exist or it has a wrong format" << endl; 
        cout << endl;
        cout << "  format : two columns ascii file -> [sky_index        value]" << endl;
        cout << "  sky_index : is the sky grid index" << endl;
        cout << "  value     : if !=0 the index sky location is used for the analysis" << endl;
        cout << endl;
        EXIT(1); 
      }
      if(skycoord=='e') net->setIndexMode(0);
      return 0;
    } 

    double theta=-1000;
    TString THETA = TB.getParameter(options,"--theta");   // get theta
    if(THETA.IsFloat()) theta=THETA.Atof();
    double phi=-1000;
    TString PHI = TB.getParameter(options,"--phi");       // get phi
    if(PHI.IsFloat()) phi=PHI.Atof();
    double radius=-1000;
    TString RADIUS = TB.getParameter(options,"--radius"); // get Radius
    if(RADIUS.IsFloat()) radius=RADIUS.Atof();

    if(theta==-1000 || phi==-1000 || radius==-1000) { // input parameter are the skyMask params
      cout << endl << "cwb::SetSkyMask - Error : wrong input skyMask params" << endl << endl; 
      cout << "wrong input options : " << options << endl;
      if(skycoord=='e')
        cout<<"options must be : --theta THETA --phi PHI --radius RADIUS"<<endl<<endl; 
      if(skycoord=='c')
        cout<<"options must be : --theta DEC --phi AR --radius RADIUS"<<endl<<endl; 
      if(fabs(theta)>90)   cout << theta << " theta must be in the range [-90,90]" << endl;      
      if(phi<0 || phi>360) cout << phi << " phi must be in the range [0,360]" << endl;
      if(radius<=0)        cout << radius << " radius must be > 0" << endl;
      cout << endl;
      EXIT(1); 
    } else {					   // create & set SkyMask
      skymap* SkyMask=NULL;
      if(skyres<0) skyres = cfg->healpix ? cfg->healpix : cfg->angle;
      if(cfg->healpix) SkyMask = new skymap(int(skyres));
      else             SkyMask = new skymap(skyres,cfg->Theta1,cfg->Theta2,cfg->Phi1,cfg->Phi2);
      MakeSkyMask(*SkyMask, theta, phi, radius);
      net->setSkyMask(*SkyMask,skycoord);
      if(skycoord=='e') net->setIndexMode(0);
      if(SkyMask!=NULL) delete SkyMask;
    }
  }
  return 0;
}

vector<frfile> 
cwb::GetFrList(int ifoID) {
//
// return the FRF frame files list
// input : ifoID : ifo iD , if ifoID=-1 then all ifo FRF are returned
//

  vector<frfile> frlist;

  // fill frlist
  for(int i=0;i<2*nIFO;i++) {
    if((ifoID==-1)||(i==ifoID)||(i==(ifoID+nIFO))) frlist.push_back(FRF[i]);
  }

  return frlist; 
}

vector<frfile> 
cwb::GetFrList(TString ifo) {
//
// return the FRF frame files list
// input : ifo : ifo name 
//

  // get ifo id
  int ifoID=-1;
  for(int n=0;n<nIFO;n++) if(ifo.CompareTo(this->ifo[n])==0) {ifoID=n;break;}
  
  if(ifoID==-1) {
    cout << "cwb::GetFrList - Error : requested ifo " << ifo
         << " not present !!!" << endl;;
    EXIT(1);
  }

  return GetFrList(ifoID); 
}

