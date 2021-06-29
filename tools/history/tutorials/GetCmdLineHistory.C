//
// Get Command Line infos from the OUTPUT CWB ROOT FILES contained in the 'idir' directory
// Author : Gabriele Vedovato
//
// Example : root 'GetCmdLineHistory.C("output")'
// Output  : /home/waveburst/soft/root/root-v5-32-04.patched/bin/root.exe -splash -n -l -b ...
//

void GetCmdLineHistory(TString idir) {

  TString data_label="";

  vector<TString> fileList = CWB::Toolbox::getFileListFromDir(idir,".root","",data_label,true);
  for(int n=0;n<fileList.size();n++) {

    // Get STOP JOB info from history
    TFile *ifile = TFile::Open(fileList[n]);
    if(ifile==NULL) {cout << "Failed to open " << fileList[n].Data() << endl;exit(-1);}
    CWB::History* ihistory = (CWB::History*)ifile->Get("history");
    if(ihistory==NULL) { cout << "Error : history is not present!!!" << endl;exit(1); }

    TString cmd_line;
    int nStages = cwb::GetStageSize();
    if(ihistory) {
      TList* stageList = ihistory->GetStageNames();    // get stage list
      for(int i=0;i<stageList->GetSize();i++) {        // loop over the stage list
        TObjString* stageObjString = (TObjString*)stageList->At(i);
        TString stageName = stageObjString->GetString();
        char* stage = const_cast<char*>(stageName.Data());
        TString info = ihistory->GetHistory(stage,const_cast<char*>("CMDLINE"));
        if(info!="") cmd_line=info;
      }
    }
    cout << cmd_line << endl;
    ifile->Close();
  }

  exit(0);
}

