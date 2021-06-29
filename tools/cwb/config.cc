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


#include "config.hh"
#include "xroot.hh"
#include "TGlobal.h"

TGlobal* global;
char cmdline[128];                                                      

#ifdef _USE_ROOT6 	// ------------------------------------------------------ ROOT6

#define SETVAR(MODE,VAR,TYPE,SIZE1,SIZE2) {                                           \
  if(MODE==0) {                                                                       \
    global = (TGlobal*)gROOT->GetGlobal("gPOINTER",true);                             \
    if(global==NULL) sprintf(cmdline,"void* gPOINTER = (void*)&%s;",#VAR);            \
    else             sprintf(cmdline,"gPOINTER = (void*)&%s;",#VAR);                  \
    gROOT->ProcessLine(cmdline);                                                      \
    global = (TGlobal*)gROOT->GetGlobal("gPOINTER",true);                             \
    if(global!=NULL) {                                                                \
      void* gPOINTER=NULL;                                                            \
      memcpy((void*)&gPOINTER,(void*)global->GetAddress(),sizeof(void*));             \
      memcpy((void*)&VAR,(void*)gPOINTER,SIZE1*SIZE2*sizeof(TYPE));                   \
    }                                                                                 \
  } else {                                                                            \
    global = (TGlobal*)gROOT->GetGlobal(#VAR,true);                                   \
    if(global==NULL) {                                                                \
      if(SIZE1==1 && SIZE2==1) strcpy(cmdline,#TYPE" "#VAR";");                       \
      if(SIZE1>1  && SIZE2==1) strcpy(cmdline,#TYPE" "#VAR"["#SIZE1"];");             \
      if(SIZE1>1  && SIZE2>1)  strcpy(cmdline,#TYPE" "#VAR"["#SIZE1"]["#SIZE2"];");   \
      gROOT->ProcessLine(cmdline);                                                    \
      global = (TGlobal*)gROOT->GetGlobal(#VAR,true);                                 \
    }                                                                                 \
    if(SIZE1==1 && SIZE2==1) {                                                        \
      void* pVAR = (void*)&VAR;                                                       \
      global = (TGlobal*)gROOT->GetGlobal(#VAR,true);                                 \
      sprintf(cmdline,#VAR" = *("#TYPE"*)%p;",pVAR);                                  \
    } else {                                                                          \
      void* pVAR = (void*)&VAR;                                                       \
      global = (TGlobal*)gROOT->GetGlobal(#VAR,true);                                 \
      sprintf(cmdline,"memcpy((void*)%p,(void*)%p,"#SIZE1"*"#SIZE2"*sizeof("#TYPE"));", \
              (void*)global->GetAddress(),pVAR);                                      \
    }                                                                                 \
    gROOT->ProcessLine(cmdline);                                                      \
 }} 

#else		// -------------------------------------------------------------- ROOT5

#define SETVAR(MODE,VAR,TYPE,SIZE1,SIZE2) {                                           \
  if(MODE==0) {                                                                       \
    global = (TGlobal*)gROOT->GetGlobal("gPOINTER",true);                             \
    if(global==NULL) sprintf(cmdline,"void* gPOINTER = (void*)&%s;",#VAR);            \
    else             sprintf(cmdline,"gPOINTER = (void*)&%s;",#VAR);                  \
    gROOT->ProcessLine(cmdline);                                                      \
    global = (TGlobal*)gROOT->GetGlobal("gPOINTER",true);                             \
    if(global!=NULL) {                                                                \
      void* gPOINTER=NULL;                                                            \
      memcpy((void*)&gPOINTER,(void*)global->GetAddress(),sizeof(void*));             \
      memcpy((void*)&VAR,(void*)gPOINTER,SIZE1*SIZE2*sizeof(TYPE));                   \
    }                                                                                 \
  } else {                                                                            \
    global = (TGlobal*)gROOT->GetListOfGlobals()->FindObject(#VAR);                   \
    if(global==NULL) {                                                                \
      if(SIZE1==1 && SIZE2==1) strcpy(cmdline,#TYPE" "#VAR";");                       \
      if(SIZE1>1  && SIZE2==1) strcpy(cmdline,#TYPE" "#VAR"["#SIZE1"];");             \
      if(SIZE1>1  && SIZE2>1)  strcpy(cmdline,#TYPE" "#VAR"["#SIZE1"]["#SIZE2"];");   \
      gROOT->ProcessLine(cmdline);                                                    \
      global = (TGlobal*)gROOT->GetGlobal(#VAR,true);                                 \
    }                                                                                 \
    sprintf(cmdline,"memcpy((void*)%p,(void*)%p,"#SIZE1"*"#SIZE2"*sizeof("#TYPE"));", \
            (void*)global->GetAddress(),(void*)&VAR);                                 \
    gROOT->ProcessLine(cmdline);                                                      \
 }} 

#endif		// ---------------------------------------------------------- END MACRO 

#define EXPORT(TYPE,VAR,CMD) {                                                  \
  TGlobal* global = (TGlobal*)gROOT->GetGlobal(#VAR,true);                      \
  char __cmdline[128];                                                          \
  if(global==NULL) sprintf(__cmdline,"%s %s;",#TYPE,CMD);                       \
  else             sprintf(__cmdline,"%s;",CMD);                                \
  gROOT->ProcessLine(__cmdline);                                                \
}


#define PRINT(ARG1,ARG2,ARG3,ARG4) { \
  printf("  "#ARG1" "#ARG2"\t= "#ARG3";\t\t// "#ARG4"\n",ARG2);  \
  }

////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
<p>The config class is designed to manage the cwb configuration

Overview:
<ol style="list-style-type: upper-roman;">
  <li><a href="#usage">Usage</a></li>
  <li><a href="#readwritetxt">Read/Write Configuration from/to text file</a></li>
  <li><a href="#readwriteobj">Read/Write Configuration Object from/to root file</a></li>
</ol>

<h3><a name="usage">I. Usage</a></h3>
<pre>
    root[] CWB::config config;           // define config object with the default parameters
    root[] cout << config.nIFO << endl;  // config.nIFO = 0 
    root[]
    root[] config.Export();              // export config parameters to CINT 
    root[]                               // all parameters are visible in CINT
    root[] cout << nIFO << endl;         // nIFO = 0 
    root[] nIFO = 3;              
    root[] cout << nIFO << endl;         // nIFO = 3 
    root[] config.Import();              // import CINT config parameters into config object
    root[] cout << config.nIFO << endl;  // config.nIFO = 3
</pre>

<h3><a name="writeframetxt">II. Read/Write Configuration from/to text file</a></h3>
<pre>
    root[] config.Import("config.C");    // import config parameters from file config.C
    root[]
    root[] config.Export("config.C");    // export config parameters from file config.C
</pre>

<h3><a name="writeframeobj">III. Read/Write Configuration Object from/to root file</a></h3>
<pre>
    root[] // --------------------------------------------------------
    root[] // Write object config to root file
    root[] // --------------------------------------------------------
    root[] TFile *ofroot = new TFile("config.root", "RECREATE");
    root[] config.Write("CC");
    root[] ofroot->Close();
    root[] 
    root[] // --------------------------------------------------------
    root[] // Read object config from root file
    root[] // --------------------------------------------------------
    root[] TFile *ifroot = new TFile("config.root");
    root[] ifroot.ls();
    root[] CWB::config *iconfig = (CWB::config*)f->Get("CC");
    root[] ifroot->Close();
</pre>

</p>

END_HTML */
////////////////////////////////////////////////////////////////////////////////

ClassImp(CWB::config)

//______________________________________________________________________________
CWB::config::config(TString umacro) {
//
// default constructor
//

  umacro=="" ? Init() : Import(umacro);
}

//______________________________________________________________________________
CWB::config::~config() {
//
// default destructor 
//

}

//______________________________________________________________________________
void
CWB::config::Browse(TBrowser *b) {
//
// Browse 
//

  View();
}

//______________________________________________________________________________
void                                                                            
CWB::config::Init() {                                             
//
// Reset parameters to 0/NULL/""/' '
//            

  strcpy(this->analysis,"");  
  online = false;   

  nIFO = 0;         
  search = ' ';                                 
  optim = false;   

  inRate= 0;
  bpp   = 0.;
  Tgap  = 0.;
  Fgap  = 0.;
  TFgap = 0.;
  fLow  = 0.;
  fHigh = 0.;
  fResample = 0;
  Acore = 0.;
  Tlpr  = 0.;

  x2or  = 0.;
  netRHO= 0.;

  levelR   = 0;
  levelF   = 0;
  levelD   = 0;
  l_low    = 0;
  l_high   = 0;

  segLen     = 0.;
  segMLS     = 0.;
  segTHR     = 0.;
  segEdge    = 0.; 
  segOverlap = 0.; 

  lagSize    = 0;
  lagStep    = 0.;
  lagOff     = 0;
  lagMax     = 0;
  lagFile    = NULL;
  strcpy(lagMode," ");
  lagSite    = NULL;
  for(int i=0;i<NIFO_MAX;i++) shift[i] = 0.;

  mlagStep   = 0;

  slagSize   = 0;
  slagMin    = 0;
  slagMax    = 0;                
  slagOff    = 0;                
  slagSite   = NULL;             
  slagFile   = NULL;             

  whiteWindow = 0.;
  whiteStride = 0.;     

  Psave = 0;                                                                               

  for(int i=0;i<NIFO_MAX;i++) dcCal[i] = 0.; 

  simulation = 0;
  iwindow    = 0.;
  nfactor    = 0;        

  for(int i=0;i<NIFO_MAX;i++) dataShift[i] = 0.; 

  mdc_shift.startMDC = 0;                                    
  mdc_shift.stopMDC  = 0;                                    
  mdc_shift.offset   = 0;                                    

  strcpy(wdmXTalk,"");    
  upTDF  = 0;
  TDSize = 0;                         
  strcpy(filter,"");                  

  pattern= 0;
  BATCH  = 0;
  LOUD   = 0;    
  subnet = 0.;  
  subcut = 0.;  

  delta = 0.;
  gamma = 0.;
  eDisbalance = false;                                                                                  

  EFEC   = false;
  mode   = 0;
  angle  = 0.;
  Theta1 = 0.;
  Theta2 = 0.;
  Phi1   = 0.; 
  Phi2   = 0.;
  mask   = 0.00;
  healpix= 0;

  precision = 0.;

  jobfOptions  = CWB_JOBF_SAVE_DISABLE;
  outfOptions  = CWB_OUTF_SAVE_DISABLE;

  dumpHistory = false;
  dump        = false;
  savemode    = false;
  cedDump     = false;
  cedRHO      = 0.;                                          
  nSky        = 0;

  for(int i=0;i<2*NIFO_MAX;i++) strcpy(frFiles[i],"");

  frRetryTime = 0;                                           

  strcpy(filter_dir,"");          

  strcpy(injectionList,"");
  strcpy(skyMaskFile,"");  
  strcpy(skyMaskCCFile,"");
  for(int i=0;i<NIFO_MAX;i++) strcpy(channelNamesRaw[i],"");
  for(int i=0;i<NIFO_MAX;i++) strcpy(channelNamesMDC[i],"");

  strcpy(work_dir   , "");
  strcpy(config_dir , "");
  strcpy(input_dir  , ""); 
  strcpy(output_dir , "");
  strcpy(merge_dir  , ""); 
  strcpy(condor_dir , "");
  strcpy(report_dir , "");
  strcpy(macro_dir  , ""); 
  strcpy(log_dir    , "");   
  strcpy(data_dir   , "");  
  strcpy(tmp_dir    , "");   
  strcpy(ced_dir    , "");
  strcpy(pp_dir     , "");
  strcpy(dump_dir   , "");    
  strcpy(www_dir    , "");    
  strcpy(data_label , "");    
  strcpy(condor_log , "");    
  strcpy(condor_tag , "");    
  strcpy(nodedir    , "");    

  nDQF = 0;                                                                                

  plugin.SetName("");
  configPlugin.SetName("");
  strcpy(parPlugin  , "");    
  dataPlugin=false;
  mdcPlugin=false;
  dcPlugin=false;
  cohPlugin=false;
  scPlugin=false;
  outPlugin=false;

  strcpy(comment    , "");    

  return;
}

//______________________________________________________________________________
void
CWB::config::Import(TString umacro) {
//
// Import from macro or CINT the configuration parameters
//
// Input: umacro  - unnamed root macro
//                  if umacro="" then parameters are imported from CINT
//
// NOTE: macro must be unnamed, parameters in macro must be declared with type
//
// WARNING: if umacro!="" all CINT global variables with the same name are overwritten !!!
//

#ifdef _USE_ROOT6
  // The interpreter of root6 if full c++ compliant
  // When unamed macros are loaded the symbols which have been redeclared with the same name are not allowed
  // The following code check if umacro has been already loaded from the root command line
  // if already loaded the macro umacro is not reloaded
  // This patch fix the job running with condor
  for(int i=0;i<gApplication->Argc();i++) {
    bool check=true;
    if(TString(gApplication->Argv(i)).EndsWith(".C")) {

      char* file1 = CWB::Toolbox::readFile(gApplication->Argv(i));
      if(file1==NULL) {check=false;continue;}
      char* file2 = CWB::Toolbox::readFile(umacro);
      if(file2==NULL) {delete [] file1;check=false;continue;}

      //cout << "file1 : " << gApplication->Argv(i) << " " << strlen(file2) << endl;
      //cout << "file2 : " << umacro << " " << strlen(file2) << endl;
      if(strlen(file1)==strlen(file2)) {
        for(int i=0;i<strlen(file1);i++) {if(file1[i]!=file2[i]) check=false;break;}
      } else check=false;
      delete [] file1; 
      delete [] file2; 
    } else check=false;
    if(check==true) return;
  } 
#endif

  int err=0;
  if(umacro!="") {
    gROOT->ProcessLine("#include \"xroot.hh\"");	// define macros SEARCH,GAMMA,XROOT
    gROOT->Macro(umacro,&err);
    if(err!=0) {
      cout << "CWB::config::Import : Error Loading Macro " << umacro.Data() << endl;
      exit(1);
    }
  }

  SetVar(0);
}

//______________________________________________________________________________
void
CWB::config::Export(TString fname) {
//
// Export to macro or CINT the configuration parameters
//
// Input: fname  - output unnamed macro file name
//                 if fname="" then parameters are exported to CINT  
//

  if(fname.Sizeof()>1) Print(fname);
  else                 SetVar(1);
}
   
//______________________________________________________________________________
void
CWB::config::SetVar(bool MODE) {
//
// Import/Export from/to CINT the configuration parameters
//
// Input: MODE  - 0/1  -> Import/Export 
//

// config parameters

  SETVAR(MODE,analysis,char,8,1);          
  SETVAR(MODE,online,bool,1,1);

  SETVAR(MODE,nIFO,int,1,1);                
#ifdef _USE_ROOT6
  char cfg_search; 
  if(MODE==0) { // import
    SETVAR(MODE,cfg_search,char,1,1); 
    search=cfg_search;
  } else {      // export
    cfg_search=search;
    SETVAR(MODE,cfg_search,char,1,1);             
  }
#else
  SETVAR(MODE,search,char,1,1);             
#endif
  SETVAR(MODE,optim,bool,1,1);

  SETVAR(MODE,ifo,char,NIFO_MAX,8);                
  SETVAR(MODE,refIFO,char,4,1);          
  SETVAR(MODE,detParms,detectorParams,NIFO_MAX,1);                

  // cWB settings

  SETVAR(MODE,inRate,size_t,1,1);             
  SETVAR(MODE,bpp,double,1,1);              
  SETVAR(MODE,Tgap,double,1,1);             
  SETVAR(MODE,Fgap,double,1,1);             
  SETVAR(MODE,TFgap,double,1,1);             
  SETVAR(MODE,fLow,double,1,1);             
  SETVAR(MODE,fHigh,double,1,1);            
  SETVAR(MODE,fResample,size_t,1,1);        

  SETVAR(MODE,Acore,double,1,1);            
  SETVAR(MODE,Tlpr,double,1,1);             

  SETVAR(MODE,x2or,double,1,1);             
  SETVAR(MODE,netRHO,double,1,1);           
  SETVAR(MODE,netCC,double,1,1);            

  // wavelet transformation settings

  SETVAR(MODE,levelR,int,1,1);              
  SETVAR(MODE,levelF,int,1,1);              
  SETVAR(MODE,levelD,int,1,1);              
  SETVAR(MODE,l_low,int,1,1);               
  SETVAR(MODE,l_high,int,1,1);              

  // segments
  SETVAR(MODE,segLen,double,1,1);           
  SETVAR(MODE,segMLS,double,1,1);           
  SETVAR(MODE,segTHR,double,1,1);           
  SETVAR(MODE,segEdge,double,1,1);          
  SETVAR(MODE,segOverlap,double,1,1);          

  // lags
  SETVAR(MODE,lagSize,size_t,1,1);          
  SETVAR(MODE,lagStep,double,1,1);          
  SETVAR(MODE,lagOff,size_t,1,1);           
  SETVAR(MODE,lagMax,size_t,1,1);           
  SETVAR(MODE,lagFile,char*,1,1);          
  SETVAR(MODE,lagMode,char,2,1);       
  SETVAR(MODE,lagSite,size_t*,1,1);         
  SETVAR(MODE,shift,double,NIFO_MAX,1);  

  // multi lags
  SETVAR(MODE,mlagStep,int,1,1);         

  // super lags
  SETVAR(MODE,slagSize,int,1,1);         
  SETVAR(MODE,slagMin,int,1,1);          
  SETVAR(MODE,slagMax,int,1,1);          
  SETVAR(MODE,slagOff,int,1,1);          
  SETVAR(MODE,slagSite,size_t*,1,1);        
  SETVAR(MODE,slagFile,char*,1,1);          

  // whitening parameters
  SETVAR(MODE,whiteWindow,double,1,1);  
  SETVAR(MODE,whiteStride,double,1,1);  

  // Skymap probability pixels to be saved in the final output root file
  SETVAR(MODE,Psave,int,1,1);

  // DC corrections
  SETVAR(MODE,dcCal,double,NIFO_MAX,1);

  // simulation parameters
  SETVAR(MODE,simulation,int,1,1);          
  SETVAR(MODE,iwindow,double,1,1);              
  SETVAR(MODE,nfactor,int,1,1);             
  SETVAR(MODE,factors,double,FACTORS_MAX,1);     

  // noise shift data
  SETVAR(MODE,dataShift,double,NIFO_MAX,1);

  // parameter to shift in time the injections (sec)
  SETVAR(MODE,mdc_shift,mdcshift,1,1);

  // delay filter
  SETVAR(MODE,wdmXTalk,char,1024,1);       
  SETVAR(MODE,upTDF,size_t,1,1);          
  SETVAR(MODE,TDSize,size_t,1,1);          
  SETVAR(MODE,filter,char,1024,1);       

  // coherence stage
  SETVAR(MODE,pattern,int,1,1);             

  // supercluster stage
  SETVAR(MODE,BATCH,int,1,1);             
  SETVAR(MODE,LOUD,int,1,1);             
  SETVAR(MODE,subnet,double,1,1);             
  SETVAR(MODE,subcut,double,1,1);             

  // regulator
  SETVAR(MODE,delta,double,1,1);            
#ifdef _USE_ROOT6
  double cfg_gamma;
  if(MODE==0) { // import
    SETVAR(MODE,cfg_gamma,double,1,1); 
    gamma=cfg_gamma;
  } else {      // export
    cfg_gamma=gamma;
    SETVAR(MODE,cfg_gamma,double,1,1); 
  }
#else
  SETVAR(MODE,gamma,double,1,1);            
#endif
                           
  SETVAR(MODE,eDisbalance,bool,1,1);

  // sky settings

  SETVAR(MODE,EFEC,bool,1,1);             
  SETVAR(MODE,mode,size_t,1,1);             
  SETVAR(MODE,angle,double,1,1);            
  SETVAR(MODE,Theta1,double,1,1);           
  SETVAR(MODE,Theta2,double,1,1);           
  SETVAR(MODE,Phi1,double,1,1);             
  SETVAR(MODE,Phi2,double,1,1);             
  SETVAR(MODE,mask,double,1,1);             
  SETVAR(MODE,healpix,size_t,1,1);          

  // error regions settings

  SETVAR(MODE,precision,double,1,1);        

  // file dump MODE

  SETVAR(MODE,jobfOptions,CWB_JOBF_OPTIONS,1,1);        
  SETVAR(MODE,outfOptions,CWB_OUTF_OPTIONS,1,1);        
  SETVAR(MODE,dumpHistory,bool,1,1);        
  SETVAR(MODE,dump,bool,1,1);               
  SETVAR(MODE,savemode,bool,1,1);           
  SETVAR(MODE,cedDump,bool,1,1);            
  SETVAR(MODE,cedRHO,double,1,1);
  SETVAR(MODE,nSky,long,1,1);             

  // directories, file names

  SETVAR(MODE,filter_dir,char,1024,1);

  SETVAR(MODE,injectionList,char,1024,1);
  SETVAR(MODE,skyMaskFile,char,1024,1);
  SETVAR(MODE,skyMaskCCFile,char,1024,1);

  SETVAR(MODE,channelNamesRaw,char,NIFO_MAX,50);
  SETVAR(MODE,channelNamesMDC,char,NIFO_MAX,50);

  // working dir
  SETVAR(MODE,work_dir,char,1024,1);

  SETVAR(MODE,config_dir,char,1024,1);
  SETVAR(MODE,input_dir,char,1024,1);
  SETVAR(MODE,output_dir,char,1024,1);
  SETVAR(MODE,merge_dir,char,1024,1);
  SETVAR(MODE,condor_dir,char,1024,1);
  SETVAR(MODE,report_dir,char,1024,1);
  SETVAR(MODE,macro_dir,char,1024,1);
  SETVAR(MODE,log_dir,char,1024,1);
  SETVAR(MODE,data_dir,char,1024,1);
  SETVAR(MODE,tmp_dir,char,1024,1);
  SETVAR(MODE,ced_dir,char,1024,1);
  SETVAR(MODE,pp_dir,char,1024,1);
  SETVAR(MODE,dump_dir,char,1024,1);
  SETVAR(MODE,www_dir,char,1024,1);

  // data label
  SETVAR(MODE,data_label,char,1024,1);

  // condor declarations
  SETVAR(MODE,condor_log,char,1024,1);

  // Define a Unique Tag for Condor Jobs
  SETVAR(MODE,condor_tag,char,1024,1);

  // frame files list : [0:nIFO-1]/[nIFO:2*nIFO-1] contains strain/mdc file names
  // If all mdc channels are in a single frame file -> mdc must be declared in the nIFO position
  SETVAR(MODE,frFiles,char,(2*NIFO_MAX),1024);
  // frame reading retry time (sec) : 0 -> disable
  SETVAR(MODE,frRetryTime,int,1,1);

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  SETVAR(MODE,nDQF,int,1,1);
  SETVAR(MODE,DQF,dqfile,DQF_MAX,1);

  // read and dump data on local disk (nodedir)
  SETVAR(MODE,nodedir,char,1024,1);

  // plugin  
  if(MODE==0) {     // import from CINT       
    global = (TGlobal*)gROOT->GetListOfGlobals()->FindObject("plugin"); 
    if(global!=NULL)
      plugin = *(TMacro*)global->GetAddress();
      if(TString(plugin.GetName())!="") {
        // check if code is empty
        TList* list = plugin.GetListOfLines();
        if(list->GetSize()==0) {
          cout << "CWB::config::SetVar - Error loading plugin : " 
               << plugin.GetTitle() << endl;
          exit(1); 
        } 
      }
  } else {	    // export to CINT
    global = (TGlobal*)gROOT->GetListOfGlobals()->FindObject("plugin");
    if(global==NULL) {
      sprintf(cmdline,"plugin;");
      EXPORT(TMacro,plugin,cmdline);
    }
    if(TString(plugin.GetName())!="") {
      // export to CINT macro plugin object 
      sprintf(cmdline,"pplugin = (TMacro*)%p;",&plugin);
      EXPORT(TMacro*,pplugin,cmdline);
      sprintf(cmdline,"plugin = TMacro(*pplugin);");	
      gROOT->ProcessLine(cmdline); 
      //sprintf(cmdline,"plugin.SetTitle("");");	// macro is only in memory (no disk)
      sprintf(cmdline,"plugin.SetTitle(\"%s\");",plugin.GetTitle());	
      gROOT->ProcessLine(cmdline); 
    }
  }
  // configPlugin  
  if(MODE==0) {     // import from CINT
    global = (TGlobal*)gROOT->GetListOfGlobals()->FindObject("configPlugin"); 
    if(global!=NULL) {
      configPlugin = *(TMacro*)global->GetAddress();
      if(TString(configPlugin.GetName())!="") {
        // check if code is empty
        TList* list = configPlugin.GetListOfLines();
        if(list->GetSize()==0) {
          cout << "CWB::config::SetVar - Error loading configPlugin : " 
               << configPlugin.GetTitle() << endl;
          exit(1); 
        } 
      }
    }
  } else {          // export to CINT
    global = (TGlobal*)gROOT->GetListOfGlobals()->FindObject("configPlugin");
    if(global==NULL) {
      sprintf(cmdline,"configPlugin;");
      EXPORT(TMacro,configPlugin,cmdline);
      global = (TGlobal*)gROOT->GetGlobal("configPlugin",true);
    }
    char tmpFile[1024]="";
    if(TString(configPlugin.GetName())!="") {
      // export to CINT macro configPlugin object 
      sprintf(cmdline,"pconfigPlugin = (TMacro*)%p;",&configPlugin);
      EXPORT(TMacro*,pconfigPlugin,cmdline);
      sprintf(cmdline,"configPlugin = TMacro(*pconfigPlugin);");	
      gROOT->ProcessLine(cmdline); 				
      //sprintf(cmdline,"configPlugin.SetTitle("");");	// macro is only in memory (no disk)	
      sprintf(cmdline,"configPlugin.SetTitle(\"%s\");",configPlugin.GetTitle());	
      gROOT->ProcessLine(cmdline); 
    }
  }
  SETVAR(MODE,parPlugin,char,1024,1);
  SETVAR(MODE,dataPlugin,bool,1,1);        
  SETVAR(MODE,mdcPlugin,bool,1,1);        
  SETVAR(MODE,dcPlugin,bool,1,1);        
  SETVAR(MODE,cohPlugin,bool,1,1);        
  SETVAR(MODE,scPlugin,bool,1,1);        
  SETVAR(MODE,outPlugin,bool,1,1);        

  // user defined comment
  SETVAR(MODE,comment,char,1024,1);

  return;
}
/*
void 
CWB::config::Streamer(TBuffer &R__b) {

   // Stream an object of class detector.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, config::IsA());
   } else {
      R__c = R__b.WriteVersion(config::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}
*/

//______________________________________________________________________________
void
CWB::config::Print(Option_t* option) {
//
// printf configuration parameters
//
// Input: option  - if option="" then parameters are displayed to the stdout
//                  otherwise to file name = option   
//

  FILE* stream = stdout;
  if(TString(option).Sizeof()>1) 
    if((stream = fopen(option, "w")) == NULL) {
      cout << "CWB::config::Print : Error opening file " << option << endl;
      exit(1);
    } 

  char comma[NIFO_MAX];
  if(nIFO==0) {for(int n=0;n<NIFO_MAX-1;n++) comma[n]=','; comma[NIFO_MAX-1]=' ';}
  else        {for(int n=0;n<nIFO-1;n++) comma[n]=','; comma[nIFO-1]=' ';}

  if(TString(option).Sizeof()>1) fprintf(stream,"{\n\n");

  fprintf(stream,"  #ifndef CWB_USER_PARAMETER_FILE\n");
  fprintf(stream,"  #define CWB_USER_PARAMETER_FILE\n\n");

  fprintf(stream,"  #include \"xroot.hh\"        // defines macro to manage ROOT5 vs ROOT6\n");

  fprintf(stream,"\n");
  fprintf(stream,"  char analysis[8] \t= \"%s\";\t\t// cWB analysis\n",analysis);
  fprintf(stream,"  bool online\t\t= %d;\t\t// true/false -> online/offline\n",online);
  fprintf(stream,"\n");
  fprintf(stream,"  int nIFO\t\t= %d;\t\t// size of network starting with first detector ifo[]\n",nIFO);
  fprintf(stream,"  SEARCH(char)\t\t= '%c';\t\t// see description below\n",search);
  fprintf(stream,"  bool optim\t\t= %d;\t\t// true -> optimal resolution likelihood analysis\n",optim);
  fprintf(stream,"\n");

  if(nIFO==0) {
    fprintf(stream,"  char ifo[NIFO_MAX][8];\n");
    fprintf(stream,"  char refIFO[4];\t\t\t// reference IFO\n");
    fprintf(stream,"  \n");
    fprintf(stream,"  // user define detectors list : is selected if detectorParams[n].name!=\"\"\n");
    fprintf(stream,"  // {name, latitude, longitude, elevation, AltX, AzX, AltY, AzY}\n");
    fprintf(stream,"  detectorParams detParms[NIFO_MAX];\n");
  } else {
    char ifos[32]="";for(int n=0;n<nIFO;n++) sprintf(ifos,"%s\"%s\"%c",ifos,ifo[n],comma[n]);
    fprintf(stream,"  char ifo[%d][8]\t= {%s};\n",nIFO,ifos);
    fprintf(stream,"  char refIFO[4]\t= \"%s\";\t\t// reference IFO\n",refIFO);
    fprintf(stream,"  \n");
    fprintf(stream,"  // user define detectors list : is selected if detectorParams[n].name!=\"\"\n");
    fprintf(stream,"  // {name, latitude, longitude, elevation, AltX, AzX, AltY, AzY}\n");
    fprintf(stream,"  detectorParams detParms[%d] = {\n",nIFO);
    for(int i=0;i<nIFO;i++) fprintf(stream,"                   {\"%s\", %g, %g, %g, %g, %g, %g, %g},\n",
                                    detParms[i].name,detParms[i].latitude,detParms[i].longitude,
                                    detParms[i].elevation,detParms[i].AltX,detParms[i].AzX,
                                    detParms[i].AltY,detParms[i].AzY); 
    fprintf(stream,"                 };\n");
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // cWB settings\n");
  fprintf(stream,"\n");

  fprintf(stream,"  size_t inRate\t\t= %lu;\t// input data rate\n",inRate);
  fprintf(stream,"  double bpp\t\t= %g;\t// probability for pixel selection\n",bpp);
  if(TString(analysis)=="1G") {
    fprintf(stream,"  double Tgap\t\t= %g;\t// time gap between clusters (sec)\n",Tgap);
    fprintf(stream,"  double Fgap\t\t= %g;\t// frequency gap between clusters (Hz)\n",Fgap);
  } else {
    fprintf(stream,"  double Tgap\t\t= %g;\t// defragmentation time gap between clusters (sec)\n",Tgap);
    fprintf(stream,"  double Fgap\t\t= %g;\t// defragmentation frequency gap between clusters (Hz)\n",Fgap);
    fprintf(stream,"  double TFgap\t\t= %g;\t// threshold on the time-frequency separation between two pixels\n",TFgap);
  }
  fprintf(stream,"  double fLow\t\t= %g;\t// low frequency of the search\n",fLow);
  fprintf(stream,"  double fHigh\t\t= %g;\t// high frequency of the search\n",fHigh);
  fprintf(stream,"  size_t fResample\t= %lu;\t\t// if>0 the inRate is resampled to fResample\n",fResample);
  fprintf(stream,"  double Acore\t\t= %g;\t// threshold for selection of core pixels\n",Acore);
  fprintf(stream,"  double Tlpr\t\t= %g;\t// training time for LPR filter\n",Tlpr);

  fprintf(stream,"\n");

  fprintf(stream,"  double x2or\t\t= %g;\t// 2 OR threshold\n",x2or);
  if(TString(analysis)=="1G") {
    fprintf(stream,"  double netRHO\t\t= %g;\t// threshold on rho\n",netRHO);
  } else {
    fprintf(stream,"  double netRHO\t\t= %g;\t// (4.5-5.5) - kills weak clusters to reduce output pixel rate (supercluster)\n",netRHO);
  }
  fprintf(stream,"  double netCC\t\t= %g;\t// threshold on network correlation\n",netCC);

  fprintf(stream,"\n");
  fprintf(stream,"  // wavelet transformation settings\n");
  fprintf(stream,"\n");

  fprintf(stream,"  int levelR\t\t= %d;\t\t// resampling level : inRate[fResample]/(2^levelR) Hz\n",levelR);
  fprintf(stream,"  int levelF\t\t= %d;\t\t// level where second LPR filter is applied\n",levelF);
  fprintf(stream,"  int levelD\t\t= %d;\t\t// decomposition level\n",levelD);
  fprintf(stream,"  int l_low\t\t= %d;\t\t// low  frequency resolution level (2^l_low Hz)\n",l_low);
  fprintf(stream,"  int l_high\t\t= %d;\t\t// high frequency resolution level (2^l_high Hz)\n",l_high);

  fprintf(stream,"\n");
  fprintf(stream,"  // time shift analysis settings\n");
  fprintf(stream,"\n");

  fprintf(stream,"  // segments\n");
  
  fprintf(stream,"  double segLen\t\t= %g;\t// Segment length [sec]\n",segLen);
  fprintf(stream,"  double segMLS\t\t= %g;\t// Minimum Segment Length after DQ_CAT1 [sec]\n",segMLS);
  fprintf(stream,"  double segTHR\t\t= %g;\t// Minimum Segment Length after DQ_CAT2 [sec]\n",segTHR);
  fprintf(stream,"  double segEdge\t= %g;\t// wavelet boundary offset [sec]\n",segEdge);
  fprintf(stream,"  double segOverlap\t= %g;\t// overlap between job segments [sec]\n",segOverlap);

  fprintf(stream,"\n");
  fprintf(stream,"  // lags\n");
  
  fprintf(stream,"  size_t lagSize\t= %lu;\t\t// number of lags (simulation=1)\n",lagSize);
  fprintf(stream,"  double lagStep\t= %g;\t// time interval between lags [sec]\n",lagStep);
  fprintf(stream,"  size_t lagOff\t\t= %lu;\t\t// first lag id (lagOff=0 - include zero lag )\n",lagOff);
  fprintf(stream,"  size_t lagMax\t\t= %lu;\t\t// 0/>0 -  standard/extended lags\n",lagMax);
  if(lagFile==NULL) fprintf(stream,"  char*  lagFile\t= NULL;\t\t// lag file list\n");
  else {
    fprintf(stream,"  char*  lagFile\t= new char[1024];\t\t// lag file list\n");
    fprintf(stream,"  sprintf(lagFile,\"%s\");\n",lagFile);
  }
  fprintf(stream,"  char   lagMode[2]\t= \"%s\";\t\t// w/r  -  write/read lag list\n",lagMode);
  if(lagSite==NULL) {
    fprintf(stream,"  size_t* lagSite\t= NULL;\t\t// site index starting with 0\n");
  } else { 
    char sites[32];for(int n=0;n<nIFO;n++) sprintf(sites,"%s%lu%c",sites,lagSite[n],comma[n]);
    fprintf(stream,"  size_t  lagSite[%d]\t= {%s};\t\t// site index starting with 0\n",nIFO,sites);
  }
  char shifts[64]="";
  if(nIFO==0) {for(int n=0;n<NIFO_MAX;n++) sprintf(shifts,"%s%1.0f%c",shifts,0.,comma[n]);
               fprintf(stream,"  double shift[NIFO_MAX] = {%s};\t// use for standard shifts\n",shifts);
  } else {     for(int n=0;n<nIFO;n++) sprintf(shifts,"%s%g%c",shifts,shift[n],comma[n]);
               fprintf(stream,"  double shift[%d]\t= {%s};\t// use for standard shifts\n",nIFO,shifts);
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // multi lags\n");
  fprintf(stream,"  int    mlagStep\t= %d;\t\t// if mlagStep=0 then 'standard lag mode'\n",mlagStep);
  fprintf(stream,"  \t\t\t\t\t// else cicle over lags with step mlagStep\n");

  fprintf(stream,"\n");
  fprintf(stream,"  // super lags\n");
  fprintf(stream,"  int    slagSize\t= %d;\t\t// number of super lags (simulation=1) - if slagSize=0 -> Standard Segments\n",slagSize);
  fprintf(stream,"  int    slagMin\t= %d;\t\t// select the minimum available slag distance : slagMin must be <= slagMax\n",slagMin);
  fprintf(stream,"  int    slagMax\t= %d;\t\t// select the maximum available slag distance\n",slagMax);
  fprintf(stream,"  int    slagOff\t= %d;\t\t// first slag id (slagOff=0 - include zero slag )\n",slagOff);
  if(slagSite==NULL) {
    fprintf(stream,"  size_t* slagSite\t= NULL;\t\t// site index starting with 0\n");
  } else { 
    char sites[32];
    for(int n=0;n<nIFO;n++) sprintf(sites,"%s%lu%c",sites,slagSite[n],comma[n]);
    fprintf(stream,"  size_t  slagSite[%d]\t= {%s};\t\t// site index starting with 0\n",nIFO,sites);
  }
  if(slagFile==NULL) fprintf(stream,"  char*  slagFile\t= NULL;\t\t// slag file list\n");
  else {
    fprintf(stream,"  char*  slagFile\t= new char[1024];\t\t// slag file list\n");
    fprintf(stream,"  sprintf(slagFile,\"%s\");\n",slagFile);
  }

  fprintf(stream,"\n");
  fprintf(stream," // whitening parameters\n");
  fprintf(stream,"  double whiteWindow\t= %g;\t// time window dT. if = 0 - dT=T, where T is wavearray duration\n",whiteWindow);
  fprintf(stream,"  double whiteStride\t= %g;\t// noise sampling interval (window stride)\n",whiteStride);

  fprintf(stream,"\n");
  fprintf(stream,"  int Psave\t\t= %d;\t\t// Skymap probability to be saved in the final output root file\n",Psave);
  fprintf(stream,"  \t\t\t\t\t// (saved if !=0 : see nSky)\n");

  fprintf(stream,"\n");
  fprintf(stream,"  // DC corrections\n");
  char dcCals[64]="";
  if(nIFO==0) {for(int n=0;n<NIFO_MAX;n++) sprintf(dcCals,"%s%1.0f%c",dcCals,0.,comma[n]);
               fprintf(stream,"  double dcCal[NIFO_MAX] = {%s};\t// use for standard dcCals\n",dcCals);
  } else {     for(int n=0;n<nIFO;n++) sprintf(dcCals,"%s%g%c",dcCals,dcCal[n],comma[n]);
               fprintf(stream,"  double dcCal[%d]\t= {%s};\t// use for standard dcCals\n",nIFO,dcCals);
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // simulation parameters\n");
  fprintf(stream,"  int simulation\t= %d;\t\t// 1 for simulation, 0 for production\n",simulation);
  fprintf(stream,"  double iwindow\t= %g;\t// analysis time window for injections (Range = Tinj +/- iwindow/2)\n",iwindow);
  fprintf(stream,"  int nfactor\t\t= %d;\t\t// number of strain factors\n",nfactor);
  if(nfactor==0) {
    fprintf(stream,"  double factors[FACTORS_MAX];\t\t\t// array of strain factors\n");
  } else {
    fprintf(stream,"  double factors[] = {\t\t\t// array of strain factors\n");
    for(int i=0;i<nfactor-1;i++) fprintf(stream,"                      %g,\n",factors[i]); 
    if(nfactor>0)                fprintf(stream,"                      %g\n",factors[nfactor-1]); 
    fprintf(stream,"                     };\n");
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // noise shift data\n");
  char dataShifts[64]="";
  if(nIFO==0) {for(int n=0;n<NIFO_MAX;n++) sprintf(dataShifts,"%s%1.0f%c",dataShifts,0.,comma[n]);
               fprintf(stream,"  double dataShift[NIFO_MAX] = {%s};\t// use for standard dataShifts\n",dataShifts);
  } else {     for(int n=0;n<nIFO;n++) sprintf(dataShifts,"%s%g%c",dataShifts,dataShift[n],comma[n]);
               fprintf(stream,"  double dataShift[%d]\t= {%s};\t// use for standard dataShifts\n",nIFO,dataShifts);
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // use this parameter to shift in time the injections (sec)\n");
  fprintf(stream,"  // use {0,0,0} to set mdc_shift to 0\n");
  fprintf(stream,"  // if {-1,0,0} the shift is automaticaly selected\n");
  fprintf(stream,"  // {startMDC, stopMDC}\n");
  fprintf(stream,"  mdcshift mdc_shift = {%g, %g, %g};\n",mdc_shift.startMDC,mdc_shift.stopMDC,mdc_shift.offset);

  fprintf(stream,"\n");
  fprintf(stream,"  // delay filter\n");
  fprintf(stream,"\n");

  if(TString(analysis)=="1G") {
    fprintf(stream,"  char   filter[1024] = \"%s\";\t\t// delay filter suffix: \"\", or \"up1\", or \"up2\" [1G]\n",filter);
  } else {
    fprintf(stream,"  \t\t\t\t\t// catalog of WDM cross-talk coefficients [2G]\n");
    fprintf(stream,"  char   wdmXTalk[1024] = \"%s\";\n",wdmXTalk);
    fprintf(stream,"  size_t upTDF\t\t= %lu;\t\t// upsample factor to obtain rate of TD filter\n",upTDF);
    fprintf(stream,"  \t\t\t\t\t// TDRate = (inRate>>levelR)*upTDF [2G]\n");
    fprintf(stream,"  size_t TDSize\t\t= %lu;\t\t// time-delay filter size (max 20) [2G]\n",TDSize);
  }

  if(TString(analysis)=="2G") {
    fprintf(stream,"\n");
    fprintf(stream,"  // coherence stage\n");
    fprintf(stream,"\n");

    fprintf(stream,"  int    pattern\t= %d;\t\t// select pixel pattern used to produce the energy max maps for pixel's selection [2G]\n",pattern);

    fprintf(stream,"\n");
    fprintf(stream,"  // supercluster stage\n");
    fprintf(stream,"\n");

    fprintf(stream,"  int    BATCH\t= %d;\t\t// max number of pixel to process in one loadTDamp batch [2G]\n",BATCH);
    fprintf(stream,"  int    LOUD\t= %d;\t\t\t// number of pixel per cluster to load TD amplitudes [2G]\n",LOUD);
    fprintf(stream,"  double subnet\t= %g;\t\t// sub network threshold (supercluster) [2G]\n",subnet);
    fprintf(stream,"  double subcut\t= %g;\t\t// sub network threshold in the skyloop (supercluster) [2G]\n",subcut);
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // regulator\n");
  fprintf(stream,"\n");

  if(TString(analysis)=="1G") {
    fprintf(stream,"  double delta\t= %g;\t\t//  [0/1] -> [weak/soft]\n",delta);
    fprintf(stream,"  GAMMA(double)\t= %g;\t\t//  set params in net5, [0/1]->net5=[nIFO/0],\n",gamma);
    fprintf(stream,"  \t\t\t\t\t//  if net5>[threshold=(nIFO-1)] weak/soft[according to delta] else hard\n");
    fprintf(stream,"  bool   eDisbalance\t= %d;\n",eDisbalance);
  } else {
    fprintf(stream,"  double delta\t= %g;\t\t// 2G: [-1:1] - regulate 2 Detector sky locations\n",delta);
    fprintf(stream,"  \t\t\t\t\t// 2G: delta=0 : regulator is disabled, delta<0 : select Lo as skystat instead of Lr\n");
    fprintf(stream,"  GAMMA(double)\t= %g;\t\t// 2G: [-1,1] - \n",gamma);
    fprintf(stream,"  \t\t\t\t\t// 2G: gamma=0 : regulator is disabled, gamma<0 : sky prior is applied\n");
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // sky settings\n");
  fprintf(stream,"\n");

  fprintf(stream,"  bool   EFEC\t\t= %d;\t\t// Earth Fixed / Selestial coordinates\n",EFEC);
  fprintf(stream,"  size_t mode\t\t= %lu;\t\t// sky search mode\n",mode);
  fprintf(stream,"  double angle\t\t= %g;\t// angular resolution\n",angle);
  fprintf(stream,"  double Theta1\t\t= %g;\t// start theta\n",Theta1);
  fprintf(stream,"  double Theta2\t\t= %g;\t// end theta\n",Theta2);
  fprintf(stream,"  double Phi1\t\t= %g;\t// start theta\n",Phi1);
  fprintf(stream,"  double Phi2\t\t= %g;\t// end theta\n",Phi2);
  fprintf(stream,"  double mask\t\t= %g;\t// sky mask fraction\n",mask);
  fprintf(stream,"  size_t healpix\t= %lu;\t\t// if not 0 use healpix sky map (healpix order)\n",healpix);

  fprintf(stream,"\n");
  fprintf(stream,"  // error regions settings\n");
  fprintf(stream,"\n");

  if(TString(analysis)=="1G") {
    fprintf(stream,"  double precision\t= %g;\t//  No = nIFO*(K+KZero)+precision*E\n",precision);
  } else {
    fprintf(stream,"  double precision\t= %g;\t//  set parameters for big clusters events management",precision);
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // file dump mode\n");
  fprintf(stream,"\n");

  fprintf(stream,"  // job file options\n");
  fprintf(stream,"  CWB_JOBF_OPTIONS jobfOptions;\n");
  fprintf(stream,"  jobfOptions = CWB_JOBF_SAVE_DISABLE;\n");
  if(jobfOptions&CWB_JOBF_SAVE_CONFIG)       fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_CONFIG;\n");
  if(jobfOptions&CWB_JOBF_SAVE_NETWORK)      fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_NETWORK;\n");
  if(jobfOptions&CWB_JOBF_SAVE_HISTORY)      fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_HISTORY;\n");
  if(jobfOptions&CWB_JOBF_SAVE_STRAIN)       fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_STRAIN;\n");
  if(jobfOptions&CWB_JOBF_SAVE_MDC)          fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_MDC;\n");
  if(jobfOptions&CWB_JOBF_SAVE_CSTRAIN)      fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_CSTRAIN;\n");
  if(jobfOptions&CWB_JOBF_SAVE_COHERENCE)    fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_COHERENCE;\n");
  if(jobfOptions&CWB_JOBF_SAVE_SUPERCLUSTER) fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_SUPERCLUSTER;\n");
  if(jobfOptions&CWB_JOBF_SAVE_LIKELIHOOD)   fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_LIKELIHOOD;\n");
  if(jobfOptions&CWB_JOBF_SAVE_CED)          fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_CED;\n");
  if(jobfOptions&CWB_JOBF_SAVE_JNET_MACRO)   fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_JNET_MACRO;\n");
  if(jobfOptions&CWB_JOBF_SAVE_CWB)          fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_CWB;\n");
  if(jobfOptions&CWB_JOBF_SAVE_SPARSE)       fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_SPARSE;\n");
  if(jobfOptions&CWB_JOBF_SAVE_WFINJ)        fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_WFINJ;\n");
  if(jobfOptions&CWB_JOBF_SAVE_WFREC)        fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_WFREC;\n");
  if(jobfOptions&CWB_JOBF_SAVE_NODE)         fprintf(stream,"  jobfOptions &= CWB_JOBF_SAVE_NODE;\n");
  fprintf(stream,"  // output root file options\n");
  fprintf(stream,"  CWB_OUTF_OPTIONS outfOptions;\n");
  fprintf(stream,"  outfOptions = CWB_OUTF_SAVE_DISABLE;\n");
  if(jobfOptions&CWB_OUTF_SAVE_VAR)          fprintf(stream,"  jobfOptions &= CWB_OUTF_SAVE_VAR;\n");
  if(jobfOptions&CWB_OUTF_SAVE_NOISE)        fprintf(stream,"  jobfOptions &= CWB_OUTF_SAVE_NOISE;\n");
  fprintf(stream,"\n");
  fprintf(stream,"  bool dumpHistory\t= %d;\t\t// dump history into output root file\n",dumpHistory);
  fprintf(stream,"  bool dump\t\t= %d;\t\t// dump triggers into ascii file\n",dump);
  fprintf(stream,"  bool savemode\t\t= %d;\t\t// temporary save clusters on disc\n",savemode);
  fprintf(stream,"  bool cedDump\t\t= %d;\t\t// dump ced plots with rho>cedRHO\n",cedDump);
  fprintf(stream,"  double cedRHO\t\t= %g;\n",cedRHO);
  fprintf(stream,"  long nSky\t\t= %lu;\t\t// # of skymap prob pixels dumped to ascii\n",nSky);
  fprintf(stream,"  \t\t\t\t\t//(nSky=0 -> (#pixels==1000 || cum prob > 0.99))\n");

  fprintf(stream,"\n");
  fprintf(stream,"  // directories, file names\n");
  fprintf(stream,"\n");

  fprintf(stream,"  char filter_dir[1024]\t\t= \"%s\";\n",filter_dir);

  fprintf(stream,"  char injectionList[1024]\t= \"%s\";\n",injectionList);
  fprintf(stream,"  char skyMaskFile[1024]\t= \"%s\";\n",skyMaskFile);
  fprintf(stream,"  char skyMaskCCFile[1024]\t= \"%s\";\n",skyMaskCCFile);
  if(nIFO==0) {
    fprintf(stream,"  char channelNamesRaw[NIFO_MAX][50];\n");
  } else {
    fprintf(stream,"  char channelNamesRaw[%d][50] = {\n",nIFO);
    for(int n=0;n<nIFO;n++) fprintf(stream,"                            \"%s\"%c\n",channelNamesRaw[n],comma[n]); 
    fprintf(stream,"                           };\n");
  }
  if(nIFO==0) {
    fprintf(stream,"  char channelNamesMDC[NIFO_MAX][50];\n");
  } else {
    fprintf(stream,"  char channelNamesMDC[%d][50] = {\n",nIFO);
    for(int n=0;n<nIFO;n++) fprintf(stream,"                            \"%s\"%c\n",channelNamesMDC[n],comma[n]); 
    fprintf(stream,"                           };\n");
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // working dir\n");
  fprintf(stream,"  char work_dir[1024]\t= \"%s\";\n",work_dir);

  fprintf(stream,"  char config_dir[1024]\t= \"%s\";\n",config_dir);
  fprintf(stream,"  char input_dir[1024]\t= \"%s\";\n",input_dir);
  fprintf(stream,"  char output_dir[1024]\t= \"%s\";\n",output_dir);
  fprintf(stream,"  char merge_dir[1024]\t= \"%s\";\n",merge_dir);
  fprintf(stream,"  char condor_dir[1024]\t= \"%s\";\n",condor_dir);
  fprintf(stream,"  char report_dir[1024]\t= \"%s\";\n",report_dir);
  fprintf(stream,"  char macro_dir[1024]\t= \"%s\";\n",macro_dir);
  fprintf(stream,"  char log_dir[1024]\t= \"%s\";\n",log_dir);
  fprintf(stream,"  char data_dir[1024]\t= \"%s\";\n",data_dir);
  fprintf(stream,"  char tmp_dir[1024]\t= \"%s\";\n",tmp_dir);
  fprintf(stream,"  char ced_dir[1024]\t= \"%s\";\n",ced_dir);
  fprintf(stream,"  char pp_dir[1024]\t= \"%s\";\n",pp_dir);
  fprintf(stream,"  char dump_dir[1024]\t= \"%s\";\n",dump_dir);
  fprintf(stream,"  char www_dir[1024]\t= \"%s\";\n",www_dir);

  fprintf(stream,"\n");
  fprintf(stream,"  // data label\n");
  fprintf(stream,"  char data_label[1024]\t= \"%s\";\n",data_label);

  fprintf(stream,"\n");
  fprintf(stream,"  // condor declarations\n");
  fprintf(stream,"  char condor_log[1024]\t= \"%s\";\n",condor_log);

  fprintf(stream,"\n");
  fprintf(stream,"  // Define a Unique Tag for Condor Jobs\n");
  fprintf(stream,"  char condor_tag[1024]\t= \"%s\";\n",condor_tag);

  fprintf(stream,"\n");
  fprintf(stream,"  // frame files list : [0:nIFO-1]/[nIFO:2*nIFO-1] contains strain/mdc file names\n");
  fprintf(stream,"  // If all mdc channels are in a single frame file -> mdc must be declared in the nIFO position\n");
  if(nIFO==0) {
    fprintf(stream,"  char frFiles[2*NIFO_MAX][1024];\n");
  } else {
    fprintf(stream,"  char frFiles[%d][1024] = {\n",2*nIFO);
    for(int n=0;n<2*nIFO-1;n++) fprintf(stream,"                    \"%s\",\n",frFiles[n]); 
                                fprintf(stream,"                    \"%s\" \n",frFiles[2*nIFO-1]); 
    fprintf(stream,"                   };\n");
    fprintf(stream,"\n");
  }
  fprintf(stream,"  // frame reading retry time (sec) : 0 -> disable\n");
  fprintf(stream,"  int frRetryTime = %d;\n",frRetryTime);

  fprintf(stream,"\n");
  fprintf(stream,"  // dq file list\n");
  fprintf(stream,"  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}\n");
  fprintf(stream,"  int nDQF = %d;\n",nDQF);
  if(nDQF==0) {
    fprintf(stream,"  dqfile DQF[DQF_MAX];\n");
  } else {
    fprintf(stream,"  dqfile DQF[%d] = {\n",nDQF);
    for(int i=0;i<nDQF;i++) fprintf(stream,"                   {\"%s\", \"%s\", (CWB_CAT)%d, %g, %d, %d},\n",
                                    DQF[i].ifo,DQF[i].file,DQF[i].cat,DQF[i].shift,DQF[i].invert,DQF[i].c4); 
    fprintf(stream,"                 };\n");
  }

  fprintf(stream,"\n");
  fprintf(stream,"  // read and dump data on local disk (nodedir)\n");
  fprintf(stream,"  char nodedir[1024]\t= \"%s\";\n",nodedir);
  fprintf(stream,"\n");

  fprintf(stream,"\n");
  fprintf(stream,"  // cwb config path\n");
  fprintf(stream,"  char cwb_config_env[1024]\t= \"%s\";\n",cwb_config_env);
  fprintf(stream,"  // cluster site name\n");
  fprintf(stream,"  char site_cluster_env[1024]\t= \"%s\";\n",site_cluster_env);
  fprintf(stream,"\n");

  fprintf(stream,"  // plugin name\n");
  fprintf(stream,"  TMacro plugin;\n");
  if(TString(plugin.GetTitle())!="") {
    TString pluginPath = plugin.GetTitle();
    pluginPath.ReplaceAll("_C.so",".C");
    fprintf(stream,"  plugin = TMacro(\"%s\");\n",pluginPath.Data());
    fprintf(stream,"  plugin.SetTitle(\"%s\");\n",plugin.GetTitle());
  }
  fprintf(stream,"  TMacro configPlugin;\n");
  if(TString(configPlugin.GetTitle())!="") {
    //fprintf(stream,"  configPlugin.SetName(\"%s\");\n",configPlugin.GetName()); // (removed to get unique md5)
    fprintf(stream,"  configPlugin = TMacro(\"%s\");\n",configPlugin.GetTitle());
  }
  fprintf(stream,"  char parPlugin[1024]\t= \"%s\";\n",parPlugin);
  fprintf(stream,"  bool dataPlugin\t= %d;\t\t// if dataPlugin=true disable read data from frames\n",dataPlugin);
  fprintf(stream,"  bool mdcPlugin\t= %d;\t\t// if mdcPlugin=true disable read data from frames\n",mdcPlugin);
  fprintf(stream,"  bool dcPlugin\t\t= %d;\t\t// if dcPlugin=true disable built-in data conditioning (only 2G)\n",dcPlugin);
  fprintf(stream,"  bool cohPlugin\t= %d;\t\t// if cohPlugin=true disable built-in coherence stage (only 2G)\n",cohPlugin);
  fprintf(stream,"  bool scPlugin\t\t= %d;\t\t// if scPlugin=true disable built-in supercluster function (only 2G)\n",scPlugin);
  fprintf(stream,"  bool outPlugin\t= %d;\t\t// if outPlugin=true disable built-in output wave file (only 2G)\n",outPlugin);
  fprintf(stream,"\n");

  fprintf(stream,"  char comment[1024]\t= \"%s\";\n",comment);

  if(TString(analysis)=="1G") {
    fprintf(stream,"  // statistics\n");
    fprintf(stream,"  // L  - likelihood\n");
    fprintf(stream,"  // c  - network correlation coefficient\n");
    fprintf(stream,"  // A  - energy disbalance asymmetry\n");
    fprintf(stream,"  // P  - penalty factor based on correlation coefficients <x,s>/sqrt(<x,x>*<s,s>)\n");
    fprintf(stream,"  // E  - total energy in the data streams\n");
    fprintf(stream,"\n");
    fprintf(stream,"  // search modes\n");
    fprintf(stream,"  // 'c' - un-modeled search, fast S5 cWB version, requires constraint settings\n");
    fprintf(stream,"  // 'h' - un-modeled search, S5 cWB version, requires constraint settings\n");
    fprintf(stream,"  // 'B' - un-modeled search, max(P*L*c/E)\n");
    fprintf(stream,"  // 'b' - un-modeled search, max(P*L*c*A/E)\n");
    fprintf(stream,"  // 'I' - elliptical plarisation, max(P*L*c/E)\n");
    fprintf(stream,"  // 'S' - linear plarisation, max(P*L*c/E)\n");
    fprintf(stream,"  // 'G' - circular plarisation, max(P*L*c/E)\n");
    fprintf(stream,"  // 'i' - elliptical plarisation, max(P*L*c*A/E)\n");
    fprintf(stream,"  // 's' - linear plarisation, max(P*L*c*A/E)\n");
    fprintf(stream,"  // 'g' - circular plarisation, max(P*L*c*A/E)\n");
  } else {
    fprintf(stream,"  //  search modes\n");
    fprintf(stream,"  //   r - un-modeled\n");
    fprintf(stream,"  //   i - iota - wave (no dispersion correction)\n");
    fprintf(stream,"  //   p - Psi - wave\n");
    fprintf(stream,"  // l,s - linear\n");
    fprintf(stream,"  // c,g - circular\n");
    fprintf(stream,"  // e,b - elliptical (no dispersion correction)\n");
  }
  fprintf(stream,"\n");

  fprintf(stream,"  #endif\n");

  if(TString(option).Sizeof()>1) fprintf(stream,"\n}\n");
  if(TString(option).Sizeof()>1) fclose(stream);
}

//______________________________________________________________________________
void 
CWB::config::DumpPlugin(const char* filename) {
//
// Save macro plugin source into filename
//
// Input: filename  -  filename where source macro is saved
//                     if filename="" then the source code is dumped to stdout
//

  if(sizeof(filename)==0) {
    cout << "CWB::config::DumpPlugin - Error : no input filename" << endl;
    cout.flush();
    return;
  }
  if(TString(filename)=="") {
    TList* list=plugin.GetListOfLines();
    for(int i=0;i<list->GetSize();i++) {
      TObjString* string = (TObjString*)list->At(i);
      cout << string->GetString() << endl;
    }
  } else {
    plugin.SaveSource(filename);
    cout << "plugin saved into : " << filename << endl;
    cout.flush();
  }
}

//______________________________________________________________________________
void 
CWB::config::DumpConfigPlugin(const char* filename) {
//
// Save macro configPlugin source into filename
//
// Input: filename  -  filename where source macro is saved
//                     if filename="" then the source code is dumped to stdout
//

  if(sizeof(filename)==0) {
    cout << "CWB::config::DumpConfigPlugin - Error : no input filename" << endl;
    cout.flush();
    return;
  }
  if(TString(filename)=="") {
    TList* list=configPlugin.GetListOfLines();
    for(int i=0;i<list->GetSize();i++) {
      TObjString* string = (TObjString*)list->At(i);
      cout << string->GetString() << endl;
    }
  } else {
    configPlugin.SaveSource(filename);
    cout << "configPlugin saved into : " << filename << endl;
    cout.flush();
  }
}

//______________________________________________________________________________
int  
CWB::config::Compare(CWB::config config) {
//
// compare config with this->config, display the differences
//
// Input: config  - config object to be compared with this object
//

  gRandom->SetSeed(0);
  int rnID = int(gRandom->Rndm(13)*1.e9);
  UserGroup_t* uinfo = gSystem->GetUserInfo();
  TString uname = uinfo->fUser;
  gSystem->Exec(TString("mkdir -p /dev/shm/")+uname);

  char ofile1[1024];
  sprintf(ofile1,"/dev/shm/%s/cwb_config_1_%d.txt",uname.Data(),rnID);
  char ofile2[1024];
  sprintf(ofile2,"/dev/shm/%s/cwb_config_2_%d.txt",uname.Data(),rnID);

  config.Print(ofile1);
  this->Print(ofile2);

  gSystem->Exec(TString("diff ")+TString(ofile1)+" "+TString(ofile2));
   
  gSystem->Exec(TString("rm ")+ofile1);
  gSystem->Exec(TString("rm ")+ofile2);

  return 0;
}

//______________________________________________________________________________
void  
CWB::config::View() {
//
// Show config with a editor defined in CWB_CONFIG_VIEWER environment
//

  gRandom->SetSeed(0);
  int rnID = int(gRandom->Rndm(13)*1.e9);
  UserGroup_t* uinfo = gSystem->GetUserInfo();
  TString uname = uinfo->fUser;
  gSystem->Exec(TString("mkdir -p /dev/shm/")+uname);

  char fName[1024];
  sprintf(fName,"/dev/shm/%s/cwb_config_%d.C",uname.Data(),rnID);
  Print(fName);

  int ret;
  char tmpStr[1024];

  if (getenv("CWB_CONFIG_VIEWER") != NULL) {
    sprintf(tmpStr, "%s %s", getenv("CWB_CONFIG_VIEWER"), fName);
    ret = system(tmpStr);
    if (ret != 0) {
      sprintf(tmpStr, "vim %s", fName);
      system(tmpStr);
    }
  }
  else {
    sprintf(tmpStr, "vim %s", fName);
    system(tmpStr);
  }

  sprintf(tmpStr, "rm -f %s", fName);
  system(tmpStr);

  return;
}

//______________________________________________________________________________
void
CWB::config::DumpConfig(const char* filename, Option_t* option) {
//
// Save configuration to file
//
// Input: filename - output file name
//

  if(sizeof(filename)==0) {
    cout << "CWB::config::DumpConfig - Error : no input filename" << endl;
    cout.flush(); 
    return;
  }
  Print(filename);
  return;
}

void 
CWB::config::SetSingleDetectorMode() {
//
// Set configuration in 1 detector mode
// The first detector (if declared) is selected
//

  if(TString(ifo[0])!="" || TString(detParms[0].name)!="") {
    cout << "------>  Set Sigle Detector Mode !!!" << endl;
    // all parameters are declared twice because cWB works only with nIFO>1
    // the analysis is done as a network of 2 equal detectors
    sprintf(ifo[1],"%s",ifo[0]);
    detParms[1] = detParms[0];
    dataShift[1] = dataShift[0];
    if(TString(ifo[0])!="") strcpy(refIFO,ifo[0]); 
    else                    strcpy(refIFO,detParms[0].name);  // user defined detector
    strcpy(channelNamesRaw[1],channelNamesRaw[0]);
    strcpy(channelNamesMDC[1],channelNamesMDC[0]);
    strcpy(frFiles[2],frFiles[nIFO]);  // mdc frFiles
    strcpy(frFiles[1],frFiles[0]);

    // select ifo[0] DQ
    int k=0;
    for(int i=0;i<nDQF;i++) if(TString(DQF[i].ifo)==refIFO) DQF[k++]=DQF[i];
    nDQF=k;
    for(int i=0;i<nDQF;i++) DQF[i+nDQF]=DQF[i];
    nDQF*=2;

    eDisbalance = false;   // disable energy disbalance
    if(slagSize!=0) {
      slagSize = 1;
      slagMin  = 0; 
      slagMax  = 0;
      slagOff  = 0;
      slagFile = NULL;
    } 
    lagSize = 1;
    lagStep = 1;
    lagOff  = 0;
    mode    = 1;           // 1 - exclude duplicate delay configurations
    delta   = 0;           // weak regulator
    if(TString(analysis)=="1G") {
      gamma = 0;           // force regulator not to be hard
    } else {
      gamma = 1.0; 
      precision = 0.0;
    }
    nSky    = 1;           // dump only 1 sky location (for single detector there is no sky loop)

    nIFO=2;

  } else {
    cout << "CWB::config::SetSingleDetectorMode - Error : ifo[0] not defined !!! " << endl;
    exit(1);
  }

  return;
}

void 
CWB::config::Check() {
//
// Check consistency of the parameters
//

  if(TString(analysis)!="1G" && TString(analysis)!="2G") 
    {cout<<"config::Check : analisys parameter non valid "<<analysis<<endl;exit(1);}

  // check consistency with the environment CWB_ANALYSIS
  if(gSystem->Getenv("CWB_ANALYSIS")!=NULL) {
    if(TString(analysis)!=TString(gSystem->Getenv("CWB_ANALYSIS"))) { 
      cout << "CWB::config::Check - Error : analysis=" << analysis; 
      cout << " is inconsistent with the environment CWB_ANALYSIS=" 
           << gSystem->Getenv("CWB_ANALYSIS") << endl;
      cout << "        check analysis parameter in user_parameters.C" << endl;
      cout << "        check CWB_ANALYSIS env in watenv setup" << endl;
      cout << "        use 'cwb_setpipe 1G/2G' command to switch analysis type" << endl;
      exit(1);
    }
  } else {
    cout << "" << endl;
    cout << "CWB_ANALYSIS env not defined in watenv setup" << endl;
    cout << "add in watenv setup the following statement" << endl;
    cout << "" << endl;
    cout << "setenv CWB_ANALYSIS        '1G'       # 1G  analysis " << endl;
    cout << "or" << endl;
    cout << "setenv CWB_ANALYSIS        '2G'       # 2G  analysis " << endl;
    cout << "" << endl;
    exit(1);
  }

  if(nIFO<=0 || nIFO>XIFO) {
    cout<<"config::Check : nIFO parameter non valid -> "<<nIFO<<endl;
    cout<<"                WAT is compiled with XIFO = "<<XIFO<<endl<<endl;
    exit(1);
  }

  if(nIFO==2) {
    TString detName1 = TString(ifo[0])!="" ? ifo[0] : detParms[0].name; 
    TString detName2 = TString(ifo[1])!="" ? ifo[1] : detParms[1].name; 
    if(detName1==detName2) {
      if(lagSize!=1) { 
        cout<<"config::Check : Error - when nIFO=2 & ifo[0]=ifo[1] -> lagSize must be 1"<<endl;
        exit(1);
      }
      if(slagSize==1) { 
        if((slagMin!=0)||(slagMax!=0)||(slagOff!=0)||(slagFile!=NULL)) { 
          cout<<"config::Check : Error - when nIFO=2 & ifo[0]=ifo[1] -> slagSize must be 0 or 1"<<endl;
          cout<<"                if slagSize=1 -> slagMin=slagMax=slagOff=0, slagFile=NULL "<<endl;
          exit(1);
        }
      }
    }
  }
  if(nIFO>2) {
    for(int n=0;n<nIFO;n++) {
      for(int m=n+1;m<nIFO;m++) {
        TString detName1 = TString(ifo[n])!="" ? ifo[n] : detParms[n].name; 
        TString detName2 = TString(ifo[m])!="" ? ifo[m] : detParms[m].name; 
        if(detName1==detName2) {
          cout<<"config::Check : Error - ifo["<<n<<"]=ifo["<<m<<"]"<<endl;
          cout<<"config::Check : when nIFO>2 detector names must be different"<<endl;
          exit(1);
        }
      }
    }
  }

  if(TString(analysis)=="1G") {
    if((search!='r')&&(search!='c')&&(search!='h')&&(search!='B')&&(search!='b')&&(search!='I')&& 
       (search!='S')&&(search!='G')&&(search!='i')&&(search!='s')&&(search!='g'))
      {cout<<"config::Check : 1G search parameter non valid "<<search<<endl;exit(1);}
  }
  if(TString(analysis)=="2G") {
    char _search = std::tolower(search);
    if((_search!='r')&&(_search!='i')&&(_search!='p')&& 
       (_search!='l')&&(_search!='s')&&(_search!='c')&&
       (_search!='g')&&(_search!='e')&&(_search!='b')) 
      {cout<<"config::Check : 2G search parameter non valid "<<search<<endl;exit(1);}
    if(fabs(delta)>1) 
      {cout<<"config::Check : 2G delta parameter non valid "<<delta<<" - must be [-1:1]"<<endl;exit(1);}
    if(fabs(gamma)>1) 
      {cout<<"config::Check : 2G gamma parameter non valid "<<gamma<<" - must be [-1:1]"<<endl;exit(1);}

    if(fmod(precision,1)!=0) {
      {cout<<"config::Check : precision must be integer : "<<precision<<endl;exit(1);}
    }
    if(precision!=0 && healpix==0) {
      {cout<<"config::Check : precision is enabled only for healpix>0 : "<<precision<<endl;exit(1);}
    }
    if(precision!=0 && healpix>0) {
      int iprecision = int(fabs(precision));
      int csize = iprecision%65536;               // get number of pixels threshold per level
      int order = (iprecision-csize)/65536;       // get resampled order
      if(csize==0)
        {cout<<"config::Check : precision must be defined with csize>0 && order>0 : "<<csize<<" " <<order<<endl;exit(1);}
      if(order==0 || order>healpix)
        {cout<<"config::Check : precision must be defined with order<=healpix : "<<order<<endl;exit(1);}
    }
  }

  if(TString(analysis)=="2G") {
    int x=upTDF; 	// must be power of 2
    if(!((x != 0) && ((x & (~x + 1)) == x)) || upTDF<=0)
      {cout<<"config::Check : upTDF  parameter non valid : must be power of 2 : "<<upTDF<<endl;exit(1);}
  }

  int x=inRate; 	// must be power of 2
  if(!((x != 0) && ((x & (~x + 1)) == x)) || inRate<=0)
    {cout<<"config::Check : inRate  parameter non valid : must be power of 2 : "<<inRate<<endl;exit(1);}

  if(!(inRate>>levelR))
    {cout<<"config::Check : levelR  parameter non valid "<<levelR<<endl;exit(1);}

  if((TString(analysis)=="1G")&&(levelF>levelD))
    {cout<<"config::Check : levelF must be <= levelD  "<<levelF<<" " <<levelD<<endl;exit(1);}

  if(l_high<l_low)
    {cout<<"config::Check : l_low must be <= l_high  "<<l_low<<" " <<l_high<<endl;exit(1);}

  if(l_low<0)
    {cout<<"config::Check : l_low must be >0 "<<l_low<<endl;exit(1);}

  if((TString(analysis)=="2G")&&(l_high-l_low+1>NRES_MAX)) {
    cout<<"config::Check : number of resolutions must le NRES_MAX="<<NRES_MAX<<endl;
    cout<<"                l_low          : "<<l_low <<endl;
    cout<<"                l_high         : "<<l_high<<endl;
    cout<<"                l_high-l_low+1 : "<<l_high-l_low+1<<endl;
    exit(1);
  }

  if((TString(analysis)=="1G")&&(l_high>levelD))
    {cout<<"config::Check : l_high must be <= levelD  "<<l_high<<" " <<levelD<<endl;exit(1);}

  if((TString(analysis)=="1G")&&((simulation<0)||(simulation>3)))
    {cout<<"config::Check : simulation parameter not valid [0:3]  "<<simulation<<endl;exit(1);}

  if((TString(analysis)=="2G")&&((simulation<0)||(simulation>4)))
    {cout<<"config::Check : simulation parameter not valid [0:4]  "<<simulation<<endl;exit(1);}

  if((simulation==0)&&(nfactor!=1))
    {cout<<"config::Check : nfactor must be 1 when simulation=0  "<<endl;exit(1);}

  if((simulation!=0)&&(lagSize!=1))
    {cout<<"config::Check : lagSize must be 1 when simulation!=0  "<<endl;exit(1);}

  if(segLen<whiteWindow)
    {cout<<"config::Check : segLen must be ge whiteWindow"<<endl;exit(1);}

  if((TString(analysis)=="1G")&&(dcPlugin || outPlugin))
    {cout<<"config::Check : dcPlugin ot outPlugin not enabled in 1G"<<endl;exit(1);}

  // segment parameters must be integers 
  if(fmod(segLen,1))       {cout<<"config::Check : segLen must be integer"<<endl;exit(1);}
  if(fmod(segMLS,1))       {cout<<"config::Check : segMLS must be integer"<<endl;exit(1);}
  if(fmod(segTHR,1))       {cout<<"config::Check : segTHR must be integer"<<endl;exit(1);}
  if(fmod(segEdge,1))      {cout<<"config::Check : segEdge must be integer"<<endl;exit(1);}
  if(fmod(segOverlap,1))   {cout<<"config::Check : segOverlap must be integer"<<endl;exit(1);}
  // segment parameters must be positive  
  if(segLen<=0)            {cout<<"config::Check : segLen must be >0"<<endl;exit(1);}
  if(segMLS<=0)            {cout<<"config::Check : segMLS must be >0"<<endl;exit(1);}
  // segment parameters must be greater than zero  
  if(segTHR<0)             {cout<<"config::Check : segTHR must be >=0"<<endl;exit(1);}
  if(segEdge<0)            {cout<<"config::Check : segEdge must be >=0"<<endl;exit(1);}
  if(segOverlap<0)         {cout<<"config::Check : segOverlap must be >=0"<<endl;exit(1);}
  // segment parameters consistent checks  
  if(segMLS>segLen)        {cout<<"config::Check : segMLS must be <= segLen"<<endl;exit(1);}

  if(simulation==0) {                                    
    if(lagStep<=0) {                                    
      cout<<"config::Check : when simulation=0 factors[0] lagStep must be >0"<<endl;
      exit(1);
    }
    if(lagMax>segMLS/lagStep) {                                    
      cout<<"config::Check : when simulation=0 lagMax must be <= segMLS/lagStep"<<endl;
      exit(1);
    }
  }

  // nfactor must be > 0
  if(nfactor<=0) {                                    
    cout<<"config::Check : nfactor must be > 0"<<endl;
    exit(1);
  }
  // if simulation==1 || simulation==2 factors must be > 0
  if(simulation==1 || simulation==2) {                                    
    for(int i=0;i<nfactor;i++) if(factors[i]<=0) {
      cout<<"config::Check : factors["<<i<<"]="<<factors[i]<<endl; 
      cout<<"config::Check : factors must be > 0"<<endl;
      exit(1);
    }
  }

  // jobfOptions=CWB_JOBF_SAVE_TRGFILE is not allowed in simulation mode
  if(simulation!=0) {                                    
    if(jobfOptions&CWB_JOBF_SAVE_TRGFILE) {
      cout<<"config::Check : jobfOptions=CWB_JOBF_SAVE_TRGFILE is not allowed in simulation mode!!!"<<endl;
      exit(1);
    } 
  } 

  // for simulation=4 factors are used as tags for the output root files
  if(simulation==4) {                                    
    if(fmod(factors[0],1)) {                                    
      cout<<"config::Check : when simulation=4 factors[0] is the offset and must be integer>0"<<endl;
      exit(1);
    }
    if(factors[0]<0) {                                    
      cout<<"config::Check : when simulation=4 factors[0] is the offset and must be integer>=0"<<endl;
      exit(1);
    }
    bool fcheck=true;
    for(int i=1;i<nfactor;i++) {
      if(factors[i]!=factors[0]+i-1) fcheck=false;
    }
    if(!fcheck) {
      for(int i=1;i<nfactor;i++) {
        if(factors[i]!=0) {                                    
          cout<<endl;
          cout<<"config::Check : when simulation=4 only factors[0] must be declared"<<endl;
          cout<<"                factors are generated automatically by the pipeline"<<endl;
          cout<<endl;
          cout<<"factors[1]=offset;"<<endl;  
          cout<<"factors[2]=offset+1;"<<endl;  
          cout<<"..."<<endl; 
          cout<<"factors[N]=offset+N-1;"<<endl;  
          cout<<endl;
          cout<<"where offset=factors[0] which must be integer>0"<<endl;
          cout<<"where N=nfactor"<<endl; 
          cout<<"if factors[0]==0 then factors[0] is set to 1"<<endl; 
          cout<<endl;
          exit(1);
        }
      }
    }
  }

  if(mdc_shift.startMDC>0 && mdc_shift.stopMDC<=mdc_shift.startMDC) {
    cout<<"config::Check : error : mdc_shift.stopMDC must be > =mdc_shift.startMDC"<<endl;
    exit(1);
  }

}

//______________________________________________________________________________
void 
CWB::config::Streamer(TBuffer &R__b) {

   // Stream an object of class CWB::config.

   //This works around a msvc bug and should be harmless on other platforms
   int size;                                                               
   typedef ::CWB::config thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      if(R__v > 1) R__b.ReadStaticArray((char*)analysis);
      if(R__v > 3) R__b >> online;
      R__b >> nIFO;
      R__b >> search;
      if(R__v > 13) R__b >> optim;
      R__b.ReadStaticArray((char*)ifo);
      R__b.ReadStaticArray((char*)refIFO);
      int R__i;
      for (R__i = 0; R__i < 9; R__i++)
         R__b.StreamObject(&(detParms[R__i]),typeid(detectorParams));
//         detParms[R__i].Streamer(R__b);
      if(R__v > 1) R__b >> inRate;
      R__b >> bpp;
      R__b >> Tgap;
      R__b >> Fgap;
      if(R__v > 11) R__b >> TFgap;
      R__b >> fLow;
      R__b >> fHigh;
      R__b >> fResample;
      R__b >> Acore;
      R__b >> Tlpr;
      R__b >> x2or;
      R__b >> netRHO;
      R__b >> netCC;
      R__b >> levelR;
      R__b >> levelF;
      R__b >> levelD;
      R__b >> l_low;
      R__b >> l_high;
      R__b >> segLen;
      R__b >> segMLS;
      R__b >> segTHR;
      R__b >> segEdge;
      if(R__v > 12) R__b >> segOverlap;
      R__b >> lagSize;
      R__b >> lagStep;
      R__b >> lagOff;
      R__b >> lagMax;
      if(R__v > 7) {                                                       
        R__b >> size; lagFile = size ? new char[size] : NULL;              
        if(lagFile) R__b.ReadStaticArray((char*)lagFile);                  
      }                                                                    
      R__b.ReadStaticArray((char*)lagMode);
      if(R__v > 7) {                                                       
        R__b >> size; lagSite = size ? new size_t[size] : NULL;            
        if(lagSite) R__b.ReadStaticArray((size_t*)lagSite);                
      }                                                                    
      R__b.ReadStaticArray((double*)shift);
      R__b >> mlagStep;
      R__b >> slagSize;
      R__b >> slagMin;
      R__b >> slagMax;
      R__b >> slagOff;
      if(R__v > 7) {                                                       
        R__b >> size; slagSite = size ? new size_t[size] : NULL;           
        if(slagSite) R__b.ReadStaticArray((size_t*)slagSite);              
      }                                                                    
      if(R__v > 7) {                                                       
        R__b >> size; slagFile = size ? new char[size] : NULL;             
        if(slagFile) R__b.ReadStaticArray((char*)slagFile);                
      }                                                                    
      R__b >> whiteWindow;
      R__b >> whiteStride;
      R__b >> Psave;
      R__b.ReadStaticArray((double*)dcCal);
      R__b >> simulation;
      R__b >> iwindow;
      R__b >> nfactor;
      R__b.ReadStaticArray((double*)factors);
      R__b.ReadStaticArray((double*)dataShift);
      R__b.StreamObject(&(mdc_shift),typeid(mdcshift));
      if(R__v > 3) R__b.ReadStaticArray((char*)wdmXTalk);
      if(R__v > 3) R__b >> upTDF;
      if(R__v > 4) R__b >> TDSize;
      R__b.ReadStaticArray((char*)filter);
      if(R__v > 20) R__b >> pattern;
      if(R__v > 4) R__b >> BATCH;
      if(R__v > 4) R__b >> LOUD;
      if(R__v > 4) R__b >> subnet;
      if(R__v > 18) R__b >> subcut;
      R__b >> delta;
      R__b >> gamma;
      R__b >> eDisbalance;
      R__b >> EFEC;
      R__b >> mode;
      R__b >> angle;
      R__b >> Theta1;
      R__b >> Theta2;
      R__b >> Phi1;
      R__b >> Phi2;
      R__b >> mask;
      if(R__v > 1) R__b >> healpix;
      if((R__v < 7)||(R__v > 19)) R__b >> precision;
      if(R__v > 2) {
         void *ptr_jobfOptions = (void*)&jobfOptions;
         R__b >> *reinterpret_cast<Int_t*>(ptr_jobfOptions);
      }
      if(R__v > 6) {
         void *ptr_outfOptions = (void*)&outfOptions;
         R__b >> *reinterpret_cast<Int_t*>(ptr_outfOptions);
      }
      bool bool_dummy;
      if(R__v < 3) R__b >> bool_dummy;  // previous was bool saveTemp
      R__b >> dumpHistory;
      R__b >> dump;
      R__b >> savemode;
      R__b >> cedDump;
      R__b >> cedRHO;
      R__b >> nSky;
      R__b.ReadStaticArray((char*)filter_dir);
      R__b.ReadStaticArray((char*)injectionList);
      R__b.ReadStaticArray((char*)skyMaskFile);
      R__b.ReadStaticArray((char*)skyMaskCCFile);
      R__b.ReadStaticArray((char*)channelNamesRaw);
      R__b.ReadStaticArray((char*)channelNamesMDC);
      R__b.ReadStaticArray((char*)work_dir);
      R__b.ReadStaticArray((char*)config_dir);
      R__b.ReadStaticArray((char*)input_dir);
      R__b.ReadStaticArray((char*)output_dir);
      R__b.ReadStaticArray((char*)merge_dir);
      R__b.ReadStaticArray((char*)condor_dir);
      R__b.ReadStaticArray((char*)report_dir);
      R__b.ReadStaticArray((char*)macro_dir);
      R__b.ReadStaticArray((char*)log_dir);
      R__b.ReadStaticArray((char*)data_dir);
      R__b.ReadStaticArray((char*)tmp_dir);
      R__b.ReadStaticArray((char*)ced_dir);
      R__b.ReadStaticArray((char*)pp_dir);
      R__b.ReadStaticArray((char*)dump_dir);
      R__b.ReadStaticArray((char*)www_dir);
      R__b.ReadStaticArray((char*)data_label);
      R__b.ReadStaticArray((char*)condor_log);
      if(R__v > 17) R__b.ReadStaticArray((char*)condor_tag);
      R__b.ReadStaticArray((char*)frFiles);
      if(R__v > 3) R__b >> frRetryTime;
      R__b >> nDQF;
      if(R__v > 23) 
        for (R__i = 0; R__i < DQF_MAX; R__i++) R__b.StreamObject(&(DQF[R__i]),typeid(dqfile));
      else
        for (R__i = 0; R__i < 20; R__i++) R__b.StreamObject(&(DQF[R__i]),typeid(dqfile));
//         DQF[R__i].Streamer(R__b);
      R__b.ReadStaticArray((char*)nodedir);
      if(R__v > 22) R__b.ReadStaticArray((char*)cwb_config_env);
      if(R__v > 22) R__b.ReadStaticArray((char*)site_cluster_env);
      plugin.Streamer(R__b);
      configPlugin.Streamer(R__b);
      if(R__v > 21) R__b.ReadStaticArray((char*)parPlugin);
      R__b >> dataPlugin;
      R__b >> mdcPlugin;
      if(R__v > 5) R__b >> dcPlugin;
      if(R__v > 14) R__b >> cohPlugin;
      if(R__v > 15) R__b >> scPlugin;
      if(R__v > 9) R__b >> outPlugin;
      if(R__v > 21) R__b.ReadStaticArray((char*)comment);
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b.WriteArray(analysis, 8);
      R__b << online;
      R__b << nIFO;
      R__b << search;
      R__b << optim;
      R__b.WriteArray((char*)ifo, 72);
      R__b.WriteArray(refIFO, 4);
      int R__i;
      for (R__i = 0; R__i < 9; R__i++)
         R__b.StreamObject(&(detParms[R__i]),typeid(detectorParams));
//         detParms[R__i].Streamer(R__b);
      R__b << inRate;
      R__b << bpp;
      R__b << Tgap;
      R__b << Fgap;
      R__b << TFgap;
      R__b << fLow;
      R__b << fHigh;
      R__b << fResample;
      R__b << Acore;
      R__b << Tlpr;
      R__b << x2or;
      R__b << netRHO;
      R__b << netCC;
      R__b << levelR;
      R__b << levelF;
      R__b << levelD;
      R__b << l_low;
      R__b << l_high;
      R__b << segLen;
      R__b << segMLS;
      R__b << segTHR;
      R__b << segEdge;
      R__b << segOverlap;
      R__b << lagSize;
      R__b << lagStep;
      R__b << lagOff;
      R__b << lagMax;
      if(lagFile!=NULL) {size=strlen(lagFile)+1;R__b << size;R__b.WriteArray(lagFile, size);}
      else R__b << 0;
      R__b.WriteArray(lagMode, 2);
      if(lagSite!=NULL) {size=NIFO_MAX;R__b << size;R__b.WriteArray(lagSite, size);}
      else R__b << 0;
      R__b.WriteArray(shift, 9);
      R__b << mlagStep;
      R__b << slagSize;
      R__b << slagMin;
      R__b << slagMax;
      R__b << slagOff;
      if(slagSite!=NULL) {size=NIFO_MAX;R__b << size;R__b.WriteArray(slagSite, size);}
      else R__b << 0;
      if(slagFile!=NULL) {size=strlen(slagFile)+1;R__b << size;R__b.WriteArray(slagFile, size);}
      else R__b << 0;
      R__b << whiteWindow;
      R__b << whiteStride;
      R__b << Psave;
      R__b.WriteArray(dcCal, 9);
      R__b << simulation;
      R__b << iwindow;
      R__b << nfactor;
      R__b.WriteArray(factors, FACTORS_MAX);
      R__b.WriteArray(dataShift, 9);
      R__b.StreamObject(&(mdc_shift),typeid(mdcshift));
      R__b.WriteArray(wdmXTalk, 1024);
      R__b << upTDF;
      R__b << TDSize;
      R__b.WriteArray(filter, 1024);
      R__b << pattern;
      R__b << BATCH;
      R__b << LOUD;
      R__b << subnet;
      R__b << subcut;
      R__b << delta;
      R__b << gamma;
      R__b << eDisbalance;
      R__b << EFEC;
      R__b << mode;
      R__b << angle;
      R__b << Theta1;
      R__b << Theta2;
      R__b << Phi1;
      R__b << Phi2;
      R__b << mask;
      R__b << healpix;
      R__b << precision;
      R__b << jobfOptions;
      R__b << outfOptions;
      R__b << dumpHistory;
      R__b << dump;
      R__b << savemode;
      R__b << cedDump;
      R__b << cedRHO;
      R__b << nSky;
      R__b.WriteArray(filter_dir, 1024);
      R__b.WriteArray(injectionList, 1024);
      R__b.WriteArray(skyMaskFile, 1024);
      R__b.WriteArray(skyMaskCCFile, 1024);
      R__b.WriteArray((char*)channelNamesRaw, NIFO_MAX*50);
      R__b.WriteArray((char*)channelNamesMDC, NIFO_MAX*50);
      R__b.WriteArray(work_dir, 1024);
      R__b.WriteArray(config_dir, 1024);
      R__b.WriteArray(input_dir, 1024);
      R__b.WriteArray(output_dir, 1024);
      R__b.WriteArray(merge_dir, 1024);
      R__b.WriteArray(condor_dir, 1024);
      R__b.WriteArray(report_dir, 1024);
      R__b.WriteArray(macro_dir, 1024);
      R__b.WriteArray(log_dir, 1024);
      R__b.WriteArray(data_dir, 1024);
      R__b.WriteArray(tmp_dir, 1024);
      R__b.WriteArray(ced_dir, 1024);
      R__b.WriteArray(pp_dir, 1024);
      R__b.WriteArray(dump_dir, 1024);
      R__b.WriteArray(www_dir, 1024);
      R__b.WriteArray(data_label, 1024);
      R__b.WriteArray(condor_log, 1024);
      R__b.WriteArray(condor_tag, 1024);
      R__b.WriteArray((char*)frFiles, 2*NIFO_MAX*1024);
      R__b << frRetryTime;
      R__b << nDQF;
      for (R__i = 0; R__i < DQF_MAX; R__i++)
         R__b.StreamObject(&(DQF[R__i]),typeid(dqfile));
//         DQF[R__i].Streamer(R__b);
      R__b.WriteArray(nodedir, 1024);
      R__b.WriteArray(cwb_config_env, 1024);
      R__b.WriteArray(site_cluster_env, 1024);
      plugin.Streamer(R__b);
      configPlugin.Streamer(R__b);
      R__b.WriteArray(parPlugin, 1024);
      R__b << dataPlugin;
      R__b << mdcPlugin;
      R__b << dcPlugin;
      R__b << cohPlugin;
      R__b << scPlugin;
      R__b << outPlugin;
      R__b.WriteArray(comment, 1024);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

