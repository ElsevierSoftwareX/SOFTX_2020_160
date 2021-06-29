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


// used by the cwb_condor_rescue

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);

  vector<int> jobList=TB.getCondorJobList(full_condor_dir, data_label);

  int max_jobs = jobList.size();

  int jobStart[max_jobs+1];
  int jobStop[max_jobs+1];

  bool jobIdStatus[max_jobs+1];
  for (int i=0;i<max_jobs+1;i++) jobIdStatus[i]=false;

  char tag[256];sprintf(tag,"%s.dag.rescue.",data_label);
  vector<TString> fileList = TB.getFileListFromDir(condor_dir, tag, "", "_wave");
  int iversion=0;
  for(int i=0;i<fileList.size();i++) {
    //cout << i << " " << fileList[i].Data() << endl;
    TObjArray* token = TString(fileList[i]).Tokenize(TString("."));
    TObjString* srescueID = (TObjString*)token->At(token->GetEntries()-1);
    if(srescueID->GetString().IsDigit()) {
      cout << i << " " << fileList[i].Data() << endl;
      int rescueID = srescueID->GetString().Atoi();
      if(iversion<rescueID) iversion=rescueID;
    }
  }
  iversion++;

  char ofile[1024];
  sprintf(ofile,"%s/%s.dag.rescue.%d",condor_dir,data_label,iversion);

  // Check if file exist
  Long_t id,size,flags,mt;
  int estat = gSystem->GetPathInfo(ofile,&id,&size,&flags,&mt);
  if (estat==0) {
    char answer[256];
    strcpy(answer,"");
    do {
      cout << "File \"" << ofile << "\" already exist" << endl;
      cout << "Do you want to overwrite the file ? (y/n) ";
      cin >> answer;
      cout << endl << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"n")==0) {
      exit(0);
    }
  }

  int nrescue=0;
  char ifile_name[1024];
  for (int i=1;i<=max_jobs;i++) {
    sprintf(ifile_name,"%s/%s/%d_%s.err",work_dir,log_dir,jobList[i],data_label);
    //cout << ifile_name << endl;
    Long_t id,size,flags,mt;
    int estat = gSystem->GetPathInfo(ifile_name,&id,&size,&flags,&mt);
    if (estat==0) {
      if (size>0) cout << ifile_name << endl;
      //cout << size << endl;
      if (size>0) {jobIdStatus[i]=true;nrescue++;}
    }
  }
  if(nrescue==0) {
    cout << "Unfinished Error Jobs : " << 0 << "/" << max_jobs << endl;
    exit(0);
  }
  cout << endl;
  cout << "New Rescue File " << endl;
  cout << ofile << endl;

  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);

  ofstream out;
  out.open(ofile,ios::out);
  int cnt = 0;
  for (int i=0;i<=max_jobs;i++) {
    if (jobIdStatus[i]) {
      cnt++;
      char ostring[256];
      int jobID=jobList[i];
      sprintf(ostring,"JOB A%i %s/%s.sub",jobID,full_condor_dir,data_label);
      out << ostring << endl;
      sprintf(ostring,"VARS A%i PID=\"%i\"",jobID,jobID);
      out << ostring << endl;
      sprintf(ostring,"RETRY A%i 3000",jobID);
      out << ostring << endl;
    }
  }
  out.close();

  cout << "Unfinished Error Jobs : " << cnt << "/" << max_jobs << endl;
  cout << endl;
  cout << "To submit condor rescued jobs type :" << endl;
  cout << "cd " << condor_dir << endl;
  sprintf(ofile,"%s.dag.rescue.%d",data_label,iversion);
  cout << "condor_submit_dag " << ofile << endl;
  cout << endl;

  exit(0);
}
