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


// multistage 2G interactive cwb pipeline

void 
cwb_inet2G(int runID, TString fName, CWB_STAGE jstage, TString uName="", bool eced=false) {

  cwb2G* CWB = NULL;

  if(eced) { 				// get fName from env option --cfg
    TString cwb_eced_opts=TString(gSystem->Getenv("CWB_ECED_OPTS"));
    if(cwb_eced_opts!="") {
      TString ECED_CFG = CWB::Toolbox::getParameter(cwb_eced_opts,"--cfg");
      if(ECED_CFG!="") {
        fName = gSystem->ExpandPathName(ECED_CFG.Data());
      }
    }
  }

  if(fName.EndsWith(".root")) {	        // job file
    CWB = new cwb2G(fName,uName,jstage);
    // Warning : define the file to be saved in CED, must be implemented
    // Must include the config in the root file and the aux user config file
    // Temporarily is declared as empty and it is not saved : see ced.cc
    gSystem->Setenv("CWB_UPARAMETERS_FILE","");	
  } else {				// user configuration file
    CWB::config icfg;
    icfg.Import("$CWB_PARAMETERS_FILE");
#ifdef _USE_ROOT6
    if(eced) { 				// skip config::check (used by eced) 
      int nIFO=1;int nfactor=1;   	
      EXPORT(int,nIFO,TString::Format("nIFO = %d",nIFO).Data())
      EXPORT(int,nfactor,TString::Format("nfactor = %d",nfactor).Data())
    }
#else
    if(eced) {nIFO=1;nfactor=1;}   	// skip config::check (used by eced) 
#endif
    icfg.Import(fName);
    strcpy(icfg.analysis,"2G");
    icfg.Export();
    CWB = new cwb2G(icfg,jstage);
    // updated standard user config (saved in CED)
    gSystem->Setenv("CWB_UPARAMETERS_FILE",fName);	
  }

  if(CWB) {
    CWB::config* cfg = CWB->GetConfig();
    TArrayC lagBuffer;  
    if(CWB->GetLagBuffer().GetSize()) {  
#ifdef _USE_ROOT6
      char* lagFile=NULL;
      char  lagMode[2];
#endif
      // read lags from job file (used in multistage analysis)
      lagBuffer = CWB->GetLagBuffer();
      lagFile = lagBuffer.GetArray();
      lagMode[0] = CWB->GetLagMode();
#ifdef _USE_ROOT6
      EXPORT(char*,lagFile,TString::Format("lagFile = (char*)%p",lagFile).Data())
      EXPORT(char*,lagMode,TString::Format("strcpy(lagMode,(char*)%p)",lagMode).Data())
      strcpy(cfg->lagMode,lagMode);
      cfg->Import();
#endif
    }
    if(eced) cfg->Import(gSystem->ExpandPathName("$CWB_MACROS/cwb_eced.C"));  // easy ced
    cfg->Import(gSystem->ExpandPathName("$CWB_MACROS/cwb_inet.C"));
    CWB->SetupStage(jstage);
#ifdef _USE_ROOT6
    char tmp_dir[1024]; strcpy(tmp_dir,cfg->tmp_dir);
#endif
    if(eced) cfg->Print(TString(tmp_dir)+"/eced_parameters.C");  // dump full parameters file
    CWB->run(runID);
    delete CWB;
  }

  gSystem->Exit(0);
}
