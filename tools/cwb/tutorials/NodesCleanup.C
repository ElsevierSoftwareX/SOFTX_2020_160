// This macro is used to remove temporary cwb files from the remote nodes
// Works only in ATLAS cluster
// WARNING : no jobs must be running !!!
// Author : Gabriele Vedovato
//
// Macro scan all the availables nodes and remove :
// wave_*,  job_*,  CWB_Plugin_Config_*   files
// 

//#define GET_SIZE 		// get total file size (MB)
#define FILE_CLEANUP		// remove temporary files
#define DUMMY			// file are not removed
//#define VERBOSE		// print remove command

int GetSize(TString dir);
int FilesCleanup(TString dir, int& fsize);

void NodesCleanup() {

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  gROOT->Macro(gSystem->ExpandPathName("$CWB_PARAMETERS_FILE"));

  // ---------------------------------------------------
  // get node list
  // ---------------------------------------------------
  vector<TString> nodeList;

  cout << "NodesCleanup.C : get node list ..." << endl;
  char cmd[1024]; strcpy(cmd,"condor_status -format \"%s\n\" Machine");
  FILE* fp = gSystem->OpenPipe(cmd, "r");
  if(fp != NULL) {
    char cmdout[1024];
    while(fgets(cmdout, sizeof(cmdout)-1, fp)!=NULL) {
      cmdout[strlen(cmdout)-1]=0; // strip new line char
      if(TString(cmdout).BeginsWith("a")) {
        bool check=false;
        // check if node is already in the list
        for(int j=0;j<(int)nodeList.size();j++) {
          if(TString(nodeList[j])==cmdout) {check=true;break;}
        }
        if(!check) nodeList.push_back(cmdout);
      }
    }
  }
  gSystem->ClosePipe(fp);

  UserGroup_t* uinfo = gSystem->GetUserInfo();
  TString uname = uinfo->fUser;

  int fsize=0;
  int tot_fsize=0;
  int tot_fsize2=0;
  int nfiles=0;
  int tot_nfiles=0;
  int nodes = nodeList.size();
  for(int i=0;i<nodes;i++) {
    //cout << i << " " << nodeList[i] << endl;
    char node_user_path[1024];
    sprintf(node_user_path,"/atlas/node/%s/user/%s",nodeList[i].Data(),uname.Data());
    cout << "NODE : " << node_user_path << " (" << i << "/" << nodes << ")" << endl;
    Long_t id,size,flags,mt;
    int estat = gSystem->GetPathInfo(node_user_path,&id,&size,&flags,&mt);
    if (estat==0) {	// dir exist
#ifdef GET_SIZE
      fsize = GetSize(node_user_path);		// get total file size (MB)
      tot_fsize+=fsize;
#endif
#ifdef FILE_CLEANUP
      nfiles = FilesCleanup(node_user_path,fsize);	// remove temporary files
      tot_fsize2+=fsize;
      tot_nfiles+=nfiles;
#endif
    }
    cout << endl;
  }

#ifdef GET_SIZE
  cout << "total file size : " << tot_fsize << " MB" << endl;
#endif
#ifdef FILE_CLEANUP
  cout << "tot deleted files : " << tot_nfiles << " - tot size " << tot_fsize2 << " MB" << endl;
#endif

  gSystem->Exit(0);
}

int FilesCleanup(TString dir, int& fsize) {

  char cmd[1024];
  int  size_tot=0;
  int  nfiles=0;  

  // remove temporary files
  vector<TString> fileList = CWB::Toolbox::getFileListFromDir(dir);
  for(int j=0;j<fileList.size();j++) {
    bool toBeRemoved=false;
    if(fileList[j].Contains("/wave_")) toBeRemoved=true;
    if(fileList[j].Contains("/job_"))  toBeRemoved=true;
    if(fileList[j].Contains("/CWB_Plugin_")) toBeRemoved=true;
    if(toBeRemoved) {
      Long_t id,size,flags,mt;
      int estat = gSystem->GetPathInfo(fileList[j].Data(),&id,&size,&flags,&mt);
      if (estat==0) {	// file exist
        //cout << size << endl;
        size_tot+=size;
        sprintf(cmd,"rm %s",fileList[j].Data());
#ifdef VERBOSE
        cout << cmd << endl;
#endif
#ifndef DUMMY
        gSystem->Exec(cmd);
#endif
        nfiles++;
      }
    }
  }

  fsize=(int)size_tot/(1024.*1024.);

  cout << "deleted files : " << nfiles << " - size " << fsize << " MB" << endl;

  return nfiles;
}

int GetSize(TString dir) {

  int size=0;

  char cmd[1024]; 
  sprintf(cmd,"du -sm %s",dir.Data());
  //cout << cmd << endl;
  FILE* fp = gSystem->OpenPipe(cmd, "r");
  if(fp != NULL) {
    char cmdout[512];
    if(fgets(cmdout, sizeof(cmdout)-1, fp)!=NULL) {
      //cout << cmdout << endl;
      TString SIZE = TString(cmdout);
      // size is separated from path by a tab (0x9)
      // 9711    /atlas/node/a0134.atlas.local/user/vedovato
      SIZE.Resize(SIZE.First(0x9));
      if(SIZE.IsDigit()) size=SIZE.Atoi();
      //cout << size << endl;
    }
  }
  gSystem->ClosePipe(fp);

  cout << "total file size " << size << " MB" << endl;

  return size;
}
