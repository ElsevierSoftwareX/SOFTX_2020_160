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


// get the cwb parameters.C used in the last stage : used by the cwb_dump command

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));

  TString cwb_dump_hist_file_name;
  if(gSystem->Getenv("CWB_DUMP_HIST_FILE_NAME")==NULL) {
    cout << "Error : environment CWB_DUMP_HIST_FILE_NAME is not defined!!!" << endl;exit(1);
  } else {
    cwb_dump_hist_file_name=TString(gSystem->Getenv("CWB_DUMP_HIST_FILE_NAME"));
  }
  if(cwb_dump_hist_file_name.Contains(".root")==0) {
    cout << "Error : " << cwb_dump_hist_file_name.Data() << " is not a root file!!!" << endl;exit(1);
  }

  TString cwb_dump_hist_mode="view";
  if(gSystem->Getenv("CWB_DUMP_HIST_MODE")!=NULL) {
    cwb_dump_hist_mode=TString(gSystem->Getenv("CWB_DUMP_HIST_MODE"));
  }

  TFile *ifile = TFile::Open(cwb_dump_hist_file_name);
  if(ifile==NULL) {cout << "Failed to open " << cwb_dump_hist_file_name.Data() << endl;exit(-1);}

  CWB::History* ihistory = (CWB::History*)ifile->Get("history");
  if(ihistory==NULL) {
    cout << "Error : history is not present!!!" << endl;exit(1);
  }

  TString config;
  TString config_md5;
  int nStages = cwb::GetStageSize();
  if(ihistory) {
    TList* stageList = ihistory->GetStageNames();    // get stage list
    for(int i=0;i<stageList->GetSize();i++) {   // get lib versions
      TObjString* stageObjString = (TObjString*)stageList->At(i);
      TString stageName = stageObjString->GetString();
      // check if stage belong to processing stages (skip pp stages)
      bool isProcessingStage=false;
      for(int j=0;j<=nStages;j++) {
        if(cwb::GetStageString((CWB_STAGE)j).Contains(stageName)) isProcessingStage=true;
      }
      if(!isProcessingStage) continue;
      char* stage = const_cast<char*>(stageName.Data());
      TString log = ihistory->GetHistory(stage,const_cast<char*>("PARAMETERS"));
      if(log!="") {
        config = log;
        config_md5 = ihistory->GetHistory(stage,const_cast<char*>("PARAMETERS_MD5"));
      }
    }
  }

  if((cwb_dump_hist_mode=="view")||(cwb_dump_hist_mode=="dump")) {

    TString label = (cwb_dump_hist_mode=="dump") ? "config" : "view";

    char configFile[512];
    TObjArray* token = TString(cwb_dump_hist_file_name).Tokenize(TString("/"));
    
    TString odir = (cwb_dump_hist_mode=="view") ? "/tmp" : dump_dir;
    sprintf(configFile,"%s/%s_%s",odir.Data(),label.Data(),
            TString(((TObjString*)token->At(token->GetEntries()-1))->GetString()).ReplaceAll(".root",".C").Data());

    ofstream out;
    out.open(configFile,ios::out);
    if (!out.good()) {cout << "Error Opening File : " << configFile << endl;exit(1);}
    out << config.Data();
    out.close();

    if(cwb_dump_hist_mode=="view") {
      char cmd[1024];
      sprintf(cmd,"vim %s",configFile);
      cout << cmd << endl;
      gSystem->Exec(cmd);
      sprintf(cmd,"rm %s",configFile);
      cout << cmd << endl;
      gSystem->Exec(cmd);
    } else {
      cout << "Write : " << configFile << endl;
    }

  } else {
    cout << endl << "Config MD5 : " << config_md5.Data() << endl << endl;
  }

  exit(0);
}
