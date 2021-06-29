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


// check job finished (look history saved in the output root file) : used by the cwb_condor command

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

/*
  vector<slag> slagList;
  vector<slag> rslagList;

  vector<TString> ifos(nIFO);
  for(int n=0;n<nIFO;n++) ifos[n]=ifo[n];

  char* slagFile = new char[256];
  sprintf(slagFile,"%s/%s.slag",dump_dir,data_label);

  vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, CWB_CAT1);
  int slagSegs=TB.getSlagJobList(cat1List, segLen).size();

  slagList=TB.getSlagList(nIFO, slagSize, slagSegs, slagOff, slagMin, slagMax, slagSite);
  rslagList=TB.getSlagList( slagList, ifos, segLen, segMLS, segEdge, nDQF, DQF, CWB_CAT1);
  rslagList=TB.getSlagList(rslagList, ifos, segLen, segTHR, segEdge, nDQF, DQF, CWB_CAT2);

  bool* bslag = new bool[slagList.size()]; 
  for(int i=0;i<slagList.size();i++) bslag[i]=false;
  int* jobStatus = new int[slagList.size()]; 
  for(int i=0;i<slagList.size();i++) jobStatus[i]=0;
  for(int i=0;i<rslagList.size();i++) {
//    printf("%14d %14d", rslagList[i].jobId, rslagList[i].slagId[0]);
//    for (int n=0; n<nIFO; n++) printf("%14d",rslagList[i].segId[n]);
//    printf("\n");
    bslag[rslagList[i].jobId]=true;
    jobStatus[rslagList[i].jobId]=1;
  }
*/

  if(nfactor<=0) nfactor=1;     // fix nfactor when nfactor is not defined
  // get the number of job submit by condor
  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);
  char condor_dag_file[1024];
  sprintf(condor_dag_file,"%s/%s%s.dag",full_condor_dir,data_label,"");
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(condor_dag_file,&id,&size,&flags,&mt);
  vector<int> jobList;
  int* jobStatus = NULL; 
  int ncondor_jobs = 0;
  if (estat==0) {
    jobList=TB.getCondorJobList(full_condor_dir, data_label);
    ncondor_jobs = jobList.size();
    jobStatus = new int[ncondor_jobs+1]; 
    for(int i=0;i<=ncondor_jobs;i++) jobStatus[i]=1;
  } else {
    cout << "cwb_condor_check: condor dag file not exist, exit" << endl;
    exit(1);
  }

  cout << "Starting reading output directory ..." << endl;
  vector<TString> fileList = TB.getFileListFromDir(output_dir,".root","","wave_");
  for(int n=0;n<fileList.size();n++) {
    //cout << n << " " << fileList[n].Data()<< endl;

    if (n%1000==0) cout << "cwb_condor check - " << n << "/" << fileList.size() << " files" << endl;

    TFile ifile(fileList[n]);
    if(!ifile.IsOpen()) {cout << "Failed to open " << fileList[n].Data() << endl;exit(-1);}

    CWB::History* ihistory = (CWB::History*)ifile.Get("history");
    if(ihistory==NULL) ihistory = (CWB::History*)ifile.Get("CWB::History"); // for back compatibility
    if(ihistory==NULL) { cout << "Error : history is not present!!!" << endl;exit(1); }

    int log_size = ihistory->GetLogSize((char*)"FULL");
    TString log = ihistory->GetLog((char*)"FULL",log_size-1);

    int jobId = TB.getJobId(fileList[n]);    // Get JOB ID
    //cout << jobId << " " << fileList[n].Data() << endl;
    if(jobId>ncondor_jobs) continue;
    if(ifile.IsZombie()) {
      jobStatus[jobId]=2;
    } else {
      // Check if "STOP JOB" is in the log history
      if(log=="STOP JOB") jobStatus[jobId]=3;
    }
  }

  int njobCondor=0;
  int njobZombie=0;
  int njobProcessed=0;
  for(int i=1;i<=ncondor_jobs;i++) {
    if(jobStatus[i]==1) njobCondor++;
    if(jobStatus[i]==2) njobZombie++;
    if(jobStatus[i]==3) njobProcessed++;
  }
//  cout << "njobZombie    : " << njoZombie << endl;
  cout << endl;
  cout << "Number of Jobs                 : " << nfactor*ncondor_jobs << endl;
  cout << "Number of Processed Jobs       : " << njobProcessed << endl;
  cout << "Number of Jobs to be Processed : " << njobCondor << endl;
  cout << endl;

  delete [] jobStatus; 

  exit(0);
}
