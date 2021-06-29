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


//!NOISE_MDC_SIMULATION
// Plugin to injected MDC 'on the fly'  

#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "cwb2G.hh"
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

// ---------------------------------------------------------------------------------
// HOW TO SETUP PLUGIN MDC_PE IN CWB USER CONFIGURATION (EXAMPLE)
// ---------------------------------------------------------------------------------

/*

Example1:
  
  TString optmdc_pe = "";                   	   // NOTE : add space at the end of each line
  optmdc_pe += "mdc_pe_dummy=true ";       	   // enable dummy signal injections
  optmdc_pe += "mdc_pe_xml_file=xml_fname.xml ";   // xml_file used for dummy injections
  strcpy(parPlugin,optmdc_pe.Data());       	   // set MDC_PE plugin parameters

Example2:
  
  TString optmdc_pe = "";                   	   // NOTE : add space at the end of each line
  optmdc_pe += "mdc_pe_xml_file=xml_fname1.xml ";  // xml_file used for dummy injections
  optmdc_pe += "mdc_pe_xml_file=xml_fname2.xml ";  // xml_file used for real  injections
  strcpy(parPlugin,optmdc_pe.Data());       	   // set MDC_PE plugin parameters

*/

// ---------------------------------------------------------------------------------
// DEFINES
// ---------------------------------------------------------------------------------

#define MDC_PE_SIM	true 
#define MDC_PE_DUMMY	false 
#define MDC_PE_INVERT	true 

// ---------------------------------------------------------------------------------
// USER SGW PLUGIN OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  bool sim;         		// (def=true)  - if false the plugin is disabled
  bool dummy;         		// (def=false) - if enabled the plugin is used for dummy signal injections. 
                        	// the signal is whitened but not injected into the data stream.
                        	// This feature can be used in CED report for signal comparison or 
                        	// to report the whitened signal into the output root file the whitened signal
                                // if dummy=true then xml_file is read from xml_file[0] or from plugin config CWB_Plugin_MDC_PE_Config.C  
  bool invert;         		// (def=true) - if enabled & if two xml_file are defined then wf1 is wf2 are switched

  vector<TString> xml_file;     // xml file is provided in user plugin option
				// if not setup then xml file is read from plugin config CWB_Plugin_MDC_PE_Config.C
};

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions	gOPT;           // global User Options
int 		gXMLID;		// index used to select the xml file from the array gOPT.xml_file

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);

using namespace CWB;

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  if(type==CWB_PLUGIN_CONFIG) {  

    ResetUserOptions();                         // set default config options
    ReadUserOptions(cfg->parPlugin);            // user config options : read from parPlugin

    if(!gOPT.sim) {				// if gOPT.sim=false -> disaple plugin
      cfg->nfactor = 1;
      cfg->simulation = 0;
      cfg->mdcPlugin=false;
      gOPT.dummy=false; 
      gOPT.xml_file.resize(0);
      return;
    }

    // check if xml files exist
    for(int n=0;n<gOPT.xml_file.size();n++) CWB::Toolbox::checkFile(gOPT.xml_file[n]);

    cfg->mdcPlugin=true;         		// disable read mdc from frames

    gXMLID=0;					// first time -> read mdc using gOPT.xml_file[0]
  }

  if(type==CWB_PLUGIN_NETWORK) {
    PrintUserOptions(cfg);      // print config options
  }

  if(!gOPT.sim) return;		// skip plugin if gOPT.sim=false

  if(type==CWB_PLUGIN_MDC) {  

    cout << "Execute CWB_Plugin_MDC_PE.C : Inject On The Fly MDC ..." << endl;

    // ---------------------------------
    // Declare mdc class 
    // On The Fly MDC Injections
    // ---------------------------------

    CWB::mdc MDC(net);
    CWB_PLUGIN_EXPORT(MDC)

    // export global variables to the config plugin
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    int xstart = (int)x->start();
    int xstop  = (int)x->stop();
    CWB_PLUGIN_EXPORT(net)
    CWB_PLUGIN_EXPORT(cfg)
    CWB_PLUGIN_EXPORT(xstart)
    CWB_PLUGIN_EXPORT(xstop)
    CWB_PLUGIN_EXPORT(gIFACTOR)
    int gIFOID=0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {gIFOID=n;break;}
    CWB_PLUGIN_EXPORT(gIFOID)

    if(gOPT.xml_file.size()==0) { 	// read plugin config CWB_Plugin_MDC_PE_Config.C
      int error=0;
      // execute config plugin
      cfg->configPlugin.Exec(NULL,&error);
      if(error) {
        cout << "Error executing macro : " << cfg->configPlugin.GetTitle() << endl;
        cfg->configPlugin.Print(); 
        gSystem->Exit(1); 
      }
    } else { 				// xml file is provided in user plugin option 
      char symlink[1024];
      sprintf(symlink,gOPT.xml_file[gXMLID]);
      TString xmlFile = CWB::Toolbox::getFileName(symlink); // get path from symbolic link
      if(xmlFile!="") {                                     // it is a symbolic link
        if(xmlFile[0]=='.') {                               // if path is relative, start with '../' then '../' is removed 
          xmlFile = xmlFile(xmlFile.Index('/')+1, xmlFile.Sizeof()-xmlFile.Index('/')-2);
        }
        cout << endl;
        cout << "---------> CWB_Plugin_MDC_PE: injections.xml = " << xmlFile << endl;
      } else {
        xmlFile=symlink;                                    // not a symbolic link
      }

      vector<TString> fileList;
      void *dir = gSystem->OpenDirectory(xmlFile);
      if(dir) {     // is directory
        // get number of files
        fileList = CWB::Toolbox::getFileListFromDir(xmlFile, ".xml", "", "");
        if(fileList.size()<cfg->nfactor) cfg->nfactor=fileList.size();
        // get file with index gIFACTOR
        char containString[1024];
        sprintf(containString,"_%d.xml",gIFACTOR);
        fileList = CWB::Toolbox::getFileListFromDir(xmlFile, ".xml", "", containString);
        xmlFile=fileList[0];
      } else {      // is file
        cfg->nfactor=1;  // set nfactor=1
      }
      cout << "---------> CWB_Plugin_MDC_PE: xmlFile = " << xmlFile << " nfactor = " << cfg->nfactor << endl;
      cout << endl;
      CWB::Toolbox::checkFile(xmlFile);     // check if file exist
      TString clbFile = xmlFile;            // set clb file
      clbFile.ReplaceAll(".xml",".clb");

      TString inspOptions="";
      inspOptions+= "--dir "+TString(cfg->nodedir)+" ";
      inspOptions+= "--xml "+xmlFile;
      MDC.SetInspiral("PE",inspOptions);
    }

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

    if(!gOPT.dummy) return;

    // if dummy=true -> restore the original strain data -> read data again from frames

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
    wavearray<double> xx(x->size());
    xx.rate(x->rate());
    xx.start(x->start());
    xx.stop(x->stop());

    // read strain again and write to x array
    frame frl(cfg->frFiles[xIFO],"READ");
    frl.setChName(cfg->channelNamesRaw[xIFO]);
    frl >> xx;
    xx.Resample(xx.rate()/(1<<cfg->levelR));              // resampling
    xx*=sqrt(1<<cfg->levelR);                             // rescaling
    for(int i=0;i<(int)x->size();i++) x->data[i] = xx[i];
  }

  if(type==CWB_PLUGIN_BCOHERENCE) {

    if(gOPT.xml_file.size()<2) return;

    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
    
    vector<TString> delObjList;
    // supercluster clusters and parse maps are removed 
    delObjList.push_back("strain");
    delObjList.push_back("mdc");
    delObjList.push_back("cstrain");
    delObjList.push_back("rms");
    delObjList.push_back("waveform");
    delObjList.push_back("csparse");
    TString jname = jfile->GetPath();
    jname.ReplaceAll(":/","");
    jfile->Close();
    gCWB2G->FileGarbageCollector(jname,"",delObjList);

    gXMLID=1;	// read again strain data and mdc using gOPT.xml_file[1]

    // perform read data and data conditioning stages
    // the injected signals (second set) are stored after the previous in the pD->IWFP / pD->IWFID arrays
    // warning the ID are stored with the same values used for the first set, see below
    gCWB2G->ReadData(0, gIFACTOR);
    gCWB2G->DataConditioning(gIFACTOR);

    // the first  inj set is stored in the first pD->IWFP.size()/2 positions
    // the second inj set is stored from pD->IWFP.size()/2 to pD->IWFP.size()-1 positions
    // we use the second inj set to store the difference between the second and the first set
    int nIFO = net->ifoListSize();  // number of detectors
    for(int i=0; i<nIFO; i++) {
      detector* pD = net->getifo(i);
      int J = pD->IWFP.size()/2;
      for (int j=0;j<J;j++) {
        // the index id must be modified to be unique
        if(pD->IWFID[j+J]<0) pD->IWFID[j+J]-=J/2; else pD->IWFID[j+J]+=J/2;
        // compute the difference
        wavearray<double>* wf1 = (wavearray<double>*)pD->IWFP[j];
        wavearray<double>* wf2 = (wavearray<double>*)pD->IWFP[j+J];
        wavearray<double>  WF2 = *wf2;
        *wf2 = CWB::mdc::GetDiff(wf2,wf1);
        // move wf2 to wf1
        if(gOPT.invert) *wf1 = WF2;
      }
    }
  } 

  return;
}

void ReadUserOptions(TString options) {

  // get plugin options 

  if(TString(options)!="") {

    //cout << "WF options : " << options << endl;
    TObjArray* token = TString(options).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++) {

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();

      if(stok.Contains("mdc_pe_sim=")) {
        TString mdc_pe_sim=stok;
        mdc_pe_sim.Remove(0,mdc_pe_sim.Last('=')+1);
        if(mdc_pe_sim=="true")  gOPT.sim=true;
        if(mdc_pe_sim=="false") gOPT.sim=false;
      }

      if(stok.Contains("mdc_pe_dummy=")) {
        TString mdc_pe_dummy=stok;
        mdc_pe_dummy.Remove(0,mdc_pe_dummy.Last('=')+1);
        if(mdc_pe_dummy=="true")  gOPT.dummy=true;
        if(mdc_pe_dummy=="false") gOPT.dummy=false;
      }

      if(stok.Contains("mdc_pe_invert=")) {
        TString mdc_pe_invert=stok;
        mdc_pe_invert.Remove(0,mdc_pe_invert.Last('=')+1);
        if(mdc_pe_invert=="true")  gOPT.invert=true;
        if(mdc_pe_invert=="false") gOPT.invert=false;
      }

      if(stok.Contains("mdc_pe_xml_file=")) {
        TString mdc_pe_xml_file=stok;
        mdc_pe_xml_file.Remove(0,mdc_pe_xml_file.Last('=')+1);
        if(gOPT.xml_file.size()<2) gOPT.xml_file.push_back(mdc_pe_xml_file);
      }
    }
  }
}

void ResetUserOptions() {

    gOPT.sim             = MDC_PE_SIM;
    gOPT.dummy           = MDC_PE_DUMMY;
    gOPT.invert          = MDC_PE_INVERT;
}

void PrintUserOptions(CWB::config* cfg) {

    cout << "-----------------------------------------"     << endl;
    cout << "MDC_PE config options                    "     << endl;
    cout << "-----------------------------------------"     << endl << endl;

    cout << "MDC_PE_SIM           " << gOPT.sim             << endl;
    cout << "MDC_PE_DUMMY         " << gOPT.dummy           << endl;
    cout << "MDC_PE_INVERT        " << gOPT.invert          << endl;

    for(int n=0;n<gOPT.xml_file.size();n++) {
      cout << "MDC_PE_XML_FILE      " << n << "\t" << gOPT.xml_file[n] << endl;
    }

    cout << endl;
}

