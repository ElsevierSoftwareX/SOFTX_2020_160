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


// merge output root job files : used by the cwb_merge command

{

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  TString cwb_merge_label;
  if(gSystem->Getenv("CWB_MERGE_LABEL")==NULL) {
    cout << "Error : environment CWB_MERGE_LABEL is not defined!!!" << endl;exit(1);
  } else {
    cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  }
  // check if label has the correct format (M#) & extract merge number
  int cwb_merge_number=0; 
  if(cwb_merge_label[0]!='M') {
    cout << "Error : label " << cwb_merge_label.Data() << " has bad format (M#)" << endl;exit(1);
  } else {
    TString lcheck=cwb_merge_label;
    lcheck.Remove(0,1);
    if(!lcheck.IsDigit()) {
      cout << "Error : label " << cwb_merge_label.Data() << " has bad format (M#)" << endl;exit(1);
    }
    cwb_merge_number=lcheck.Atoi();
  }

  // get merge options

  bool cwb_bpsm     = false;
  bool cwb_brms     = false;
  bool cwb_bvar     = false;
  int  cwb_nthreads = 0;
  TString cwb_utag  = "";

  TString cwb_merge_opts=TString(gSystem->Getenv("CWB_MERGE_OPTS"));
  if(cwb_merge_opts!="") {
    TString option="";
    // get bpsm
    option = TB.getParameter(cwb_merge_opts,"--psm"); 
    option.ToUpper();
    cwb_bpsm  = (option=="TRUE") ? true : false;
    // get brms
    option = TB.getParameter(cwb_merge_opts,"--rms");
    option.ToUpper();
    cwb_brms  = (option=="TRUE") ? true : false;
    // get bvar
    option = TB.getParameter(cwb_merge_opts,"--var");
    option.ToUpper();
    cwb_bvar  = (option=="TRUE") ? true : false;
    // get ntreads
    option = TB.getParameter(cwb_merge_opts,"--nthreads");
    if(option!="") {
      if(!option.IsDigit()) {
        cout << endl << "cwb_merge.C : Error - --nthread argument '" 
             << option << "' is not a number" << endl << endl;
        exit(1);
      } else {
        cwb_nthreads  = option.Atoi();
      }
    }
    // get input dir (default is output_dir defined in $CWB_UPARAMETERS_FILE)
    option = TB.getParameter(cwb_merge_opts,"--idir");
    if(option!="") strcpy(output_dir,option.Data());
    // get user tag : is added to the merged files .U_utag
    option = TB.getParameter(cwb_merge_opts,"--utag");
    if(option!="") cwb_utag = option;
  }

  //vector<TString> fileList = TB.getFileListFromDir(output_dir,".root");
  //cout << "getDirFileLists - nfiles : " << fileList.size() << endl;
  char tag[256];sprintf(tag,"wave_%s",data_label);
  vector<TString> fileList = TB.getFileListFromDir(merge_dir, ".root", tag,"",true);

  TString mtag;
  if(cwb_utag!="") mtag = TString::Format(".%s.U_%s.root",cwb_merge_label.Data(),cwb_utag.Data());
  else             mtag = TString::Format(".%s.root",cwb_merge_label.Data());
  if(cwb_merge_label.CompareTo("M0")!=0) {
    // check if label already exist
    for(int i=0;i<fileList.size();i++) {
      if(fileList[i].Contains(mtag)) {
        cout << "Ecwb_merge.C : Error - merge label " << cwb_merge_label << " already exist" << endl << endl; 
        gSystem->Exit(1); 
      }
    }
  }
  // compute new version for label
  int iversion=0;
  int ipos = (cwb_utag=="") ? 2 : 3;
  for(int i=0;i<fileList.size();i++) {
    cout << i << " " << fileList[i].Data() << endl;
    if(cwb_utag!="") if(!fileList[i].Contains(TString::Format(".U_%s",cwb_utag.Data()))) continue;
    if(fileList[i].Contains(".root")) {
      TObjArray* token = TString(fileList[i]).Tokenize(TString("."));
      TString srescueID = ((TObjString*)token->At(token->GetEntries()-ipos))->GetString();
      //srescueID.ReplaceAll(".root","");
      srescueID.ReplaceAll("M","");
      if(srescueID.IsDigit()) {
        //cout << i << " " << fileList[i].Data() << endl;
        int rescueID = srescueID.Atoi();
        if(iversion<rescueID) iversion=rescueID;
      }
    }
  }
  iversion++;

  if(cwb_merge_number==0 || cwb_merge_number>=iversion) {	// iversion replaced with the input user merge version
     if(cwb_merge_number!=0) iversion=cwb_merge_number;
  } else {
     cout << endl << "cwb_merge.C : Error - the input merge version (M" << cwb_merge_number 
          << ") must be greater of the most recent merge version (M" << iversion-1 << ")" << endl << endl;
     gSystem->Exit(1);
  } 

  // merge mdc log files (if produced by "On the Fly MDC")
#ifdef _USE_ROOT6
  char cmd[128]; sprintf(cmd,"iversion = %d",iversion);
  EXPORT(int,iversion,cmd)
#endif
  gROOT->Macro(gSystem->ExpandPathName("$CWB_MACROS/cwb_merge_log.C"));

  char new_data_label[256]; 
  if(cwb_utag=="") sprintf(new_data_label,"%s.M%d",data_label,iversion);
  else             sprintf(new_data_label,"%s.M%d.U_%s",data_label,iversion,cwb_utag.Data());
  cout << "Start merging to file : " << new_data_label << endl << endl;
  vector<TString> mergeList; 
  if(cwb_nthreads<=1) {
    mergeList=TB.mergeCWBTrees(output_dir,simulation,merge_dir,new_data_label,cwb_brms,cwb_bvar,cwb_bpsm);
  } else {
    // multi threads
    mergeList=TB.mergeCWBTrees(cwb_nthreads,output_dir,simulation,merge_dir,new_data_label,cwb_brms,cwb_bvar,cwb_bpsm);
  }
  char merge_list_file[1024]; 
  if(cwb_utag=="") sprintf(merge_list_file,"%s/merge_%s.M%d.lst",merge_dir,data_label,iversion);
  else             sprintf(merge_list_file,"%s/merge_%s.M%d.U_%s.lst",merge_dir,data_label,iversion,cwb_utag.Data());
  TB.dumpFileList(mergeList, merge_list_file);

  if(nfactor<=0) nfactor=1;	// fix nfactor when nfactor is not defined

  // get the number of job submit by condor
  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);
  char condor_dag_file[1024];
  sprintf(condor_dag_file,"%s/%s%s.dag",full_condor_dir,data_label,"");
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(condor_dag_file,&id,&size,&flags,&mt);
  vector<int> jobList;
  if (estat==0) {
    jobList=TB.getCondorJobList(full_condor_dir, data_label);
    int ncondor_jobs = jobList.size();
    cout << endl << "Merged files : " << mergeList.size() << "/" << nfactor*ncondor_jobs << endl << endl;
  }

  // for standard merge labels the merge/cwb_user_parameters.C is created
  // this file is extract from the merge/wave_'label'.root
  // it is the merge between $CWB_PARAMETERS_FILE & CWB_UPARAMETERS_FILE
  // it is used in the post production stage

  if(cwb_utag!="") gSystem->Exit(0);

  char wave_name[1024]; sprintf(wave_name,"%s/wave_%s.root",merge_dir,new_data_label);
  TFile *ifile = TFile::Open(wave_name);
  if(ifile==NULL) {cout << "Failed to open " << wave_name << endl;exit(1);}

  CWB::History* ihistory = (CWB::History*)ifile->Get("history");
  if(ihistory==NULL) {
    cout << "Error : history is not present!!!" << endl;exit(1);
  }

  TString config;
  TString git_version;
  int nStages = cwb::GetStageSize();
  if(ihistory) {
    TList* stageList = ihistory->GetStageNames();    	// get stage list
    for(int i=0;i<stageList->GetSize();i++) {   	// get lib versions
      TObjString* stageObjString = (TObjString*)stageList->At(i);
      TString stageName = stageObjString->GetString();
      // check if stage belong to processing stages (skip pp stages)
      bool isProcessingStage=false;
      for(int j=0;j<=nStages;j++) {
        if(cwb::GetStageString((CWB_STAGE)j).Contains(stageName)) isProcessingStage=true;
      }
      if(!isProcessingStage) continue;
      char* stage = const_cast<char*>(stageName.Data());
      TString slog = ihistory->GetHistory(stage,const_cast<char*>("PARAMETERS"));
      if(slog!="") {
        config = slog;
        git_version = ihistory->GetHistory(stage,const_cast<char*>("GITVERSION"));
      }
    }
  }

  char configFile[1024]; sprintf(configFile,"%s/cwb_user_parameters.C",merge_dir);
  ofstream out;
  out.open(configFile,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << configFile << endl;exit(1);}
  out << config.Data();
  out.close();

  gSystem->Exit(0);
}
