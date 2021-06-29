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


// compare list of jobs included in the dag file and check job finished (from history) produce condor dag file : used by the cwb_condor command

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  if(TString(condor_tag)=="") {
    cout << endl;
    cout << "cwb_condor_recovery.C : Error - the accounting_group is not defined !!!" << endl;
    cout << "The accounting_group must be defined in the user_parameters.C file" << endl;
    cout << "See the following link:" << endl;
    cout <<" https://ldas-gridmon.ligo.caltech.edu/accounting/condor_groups/determine_condor_account_group.html" << endl;
    cout << "Examples : " << endl;
    cout << "strcpy(condor_tag,\"ligo.dev.o2.burst.allsky.cwboffline\");" << endl;
    cout << "strcpy(condor_tag,\"ligo.prod.o2.burst.allsky.cwboffline\");" << endl;
    cout << "If you don't need it set : strcpy(condor_tag,\"disabled\");" << endl << endl;
    exit(1);
  }
  if(TString(condor_tag)=="disabled") strcpy(condor_tag,"");

  TString cwb_uparameters_file;
  if(gSystem->Getenv("CWB_UPARAMETERS_FILE")==NULL) {
    cout << "Error : environment CWB_UPARAMETERS_FILE is not defined!!!" << endl;exit(1);
  } else {
    cwb_uparameters_file=TString(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  }

  // get cwb stage name
  TString cwb_stage_name="CWB_STAGE_FULL";
  if(gSystem->Getenv("CWB_STAGE_NAME")!=NULL) {
    cwb_stage_name=TString(gSystem->Getenv("CWB_STAGE_NAME"));
  }
  if(cwb_stage_name=="CWB_STAGE_FULL") cwb_stage_name="CWB_STAGE_LIKELIHOOD";
  // convert stage name to value
  TString cwb_resume_label=output_dir;	// used with resume 
  CWB_STAGE cwb_stage;
  if(cwb_stage_name=="CWB_STAGE_FULL")         {cwb_stage=CWB_STAGE_LIKELIHOOD;	 cwb_resume_label+="/supercluster_";}
  if(cwb_stage_name=="CWB_STAGE_INIT")         {cwb_stage=CWB_STAGE_INIT;	 cwb_resume_label+="";}
  if(cwb_stage_name=="CWB_STAGE_STRAIN")       {cwb_stage=CWB_STAGE_STRAIN;	 cwb_resume_label+="/init_";}
  if(cwb_stage_name=="CWB_STAGE_CSTRAIN")      {cwb_stage=CWB_STAGE_CSTRAIN;	 cwb_resume_label+="/strain_";}
  if(cwb_stage_name=="CWB_STAGE_COHERENCE")    {cwb_stage=CWB_STAGE_COHERENCE;	 cwb_resume_label+="/cstrain_";}
  if(cwb_stage_name=="CWB_STAGE_SUPERCLUSTER") {cwb_stage=CWB_STAGE_SUPERCLUSTER;cwb_resume_label+="/coherence_";}
  if(cwb_stage_name=="CWB_STAGE_LIKELIHOOD")   {cwb_stage=CWB_STAGE_LIKELIHOOD;	 cwb_resume_label+="/supercluster_";}
  if(gSystem->Getenv("CWB_STAGE_NAME")!=NULL) {
    cwb_stage_name=TString(gSystem->Getenv("CWB_STAGE_NAME"));
  }

  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);

  // check if condor dag file is present, otherwise it is created 
  // the dag file could be not present in the second stage analysis
  bool exists = TB.isFileExisting(TString::Format("%s/%s.dag",full_condor_dir,data_label));
  if(!exists) {
    TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
    TString exec_cmd = TString::Format("%s/cwb_condor.csh create",cwb_scripts.Data());
    int ret=gSystem->Exec(exec_cmd);
    if(ret) {cout << "Error  while executing cwb_condor create !!!" << endl;exit(1);}
  }

  // read condor job list
  vector<int> jobList=TB.getCondorJobList(full_condor_dir, data_label);

  int max_jobs = 0;
  for(int i=0;i<jobList.size();i++) if(jobList[i]>max_jobs) max_jobs=jobList[i];

  int* jobStart = new int[max_jobs];
  int* jobStop = new int[max_jobs];

  CWB_STAGE* jobStage = new CWB_STAGE[max_jobs];
  for (int i=0;i<max_jobs;i++) jobStage[i]=(CWB_STAGE)-1;	             // excluded jobs
  for (int i=0;i<jobList.size();i++) jobStage[jobList[i]-1]=(CWB_STAGE)0;    // jobs in the dag file

  char tag[256];sprintf(tag,"%s.dag.recovery.",data_label);
  vector<TString> fileList = TB.getFileListFromDir(condor_dir, "", tag);
  int iversion=0;
  for(int i=0;i<fileList.size();i++) {
    //cout << i << " " << fileList[i].Data() << endl;
    TObjArray* token = TString(fileList[i]).Tokenize(TString("."));
    TObjString* srecoveryID = (TObjString*)token->At(token->GetEntries()-1);
    if(srecoveryID->GetString().IsDigit()) {
      cout << i << " " << fileList[i].Data() << endl;
      int recoveryID = srecoveryID->GetString().Atoi();
      if(iversion<recoveryID) iversion=recoveryID;
    }
  }
  iversion++;

  char dagfile[1024];
  sprintf(dagfile,"%s/%s.dag.recovery.%d",condor_dir,data_label,iversion);

  // Check if dag file already exist
  Long_t id,size,flags,mt;
  int estat = gSystem->GetPathInfo(dagfile,&id,&size,&flags,&mt);
  if (estat==0) {
    char answer[256];
    strcpy(answer,"");
    do {
      cout << "File \"" << dagfile << "\" already exist" << endl;
      cout << "Do you want to overwrite the file ? (y/n) ";
      cin >> answer;
      cout << endl << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"n")==0) {
      exit(0);
    }
  }

  // get cwb_stage_resume 
  // if true then recovery done only if previous cwb_stage is present in the output dir
  TString cwb_stage_resume="FALSE";
  if(gSystem->Getenv("CWB_STAGE_RESUME")!=NULL) {
    cwb_stage_resume=TString(gSystem->Getenv("CWB_STAGE_RESUME"));
  }

  // get cwb stage input dir (input files produced by the previous stage)
  TString cwb_stage_input=output_dir;
  if(gSystem->Getenv("CWB_STAGE_INPUT")!=NULL) {
    cwb_stage_input=TString(gSystem->Getenv("CWB_STAGE_INPUT"));
  }
  if(cwb_stage_input=="") cwb_stage_input=output_dir;
  TB.checkFile(cwb_stage_input);

  // factor label : extract the last factor value
  char sfactor[32]="";
  if(simulation) {
    if(simulation==3) {
      if(factor<0)  sprintf(sfactor,"_n%g",fabs(factors[nfactor-1]));
      if(factor==0) sprintf(sfactor,"_z%g",factors[nfactor-1]);
      if(factor>0)  sprintf(sfactor,"_p%g",factors[nfactor-1]);
    } else if(simulation==4) {
      int ioffset = int(factors[0])<=0 ? 1 : int(factors[0]);
      ioffset+=nfactor-1;
      sprintf(sfactor,"_%i",ioffset);
    } else          sprintf(sfactor,"_%g",factors[nfactor-1]);
  }
  char job_label[512];sprintf(job_label,"%s%s",data_label,sfactor);

  cout << "Starting reading output directory ..." << endl;
  vector<TString> jobFiles(max_jobs);
  for(int i=0;i<max_jobs;i++) jobFiles[i]=cwb_uparameters_file; 
  fileList = TB.getFileListFromDir(cwb_stage_input,".root","",data_label,true);
  for(int n=0;n<fileList.size();n++) {

    int jobId = TB.getJobId(fileList[n]);    // Get JOB ID
    jobId-=1;
     
    if(fileList[n].BeginsWith(cwb_stage_input+"/init_")) 
      if(CWB_STAGE_INIT>jobStage[jobId]) 
        {jobStage[jobId]=CWB_STAGE_INIT;jobFiles[jobId]=fileList[n];continue;}  
    if(fileList[n].BeginsWith(cwb_stage_input+"/strain_")) 
      if(CWB_STAGE_STRAIN>jobStage[jobId])       
        {jobStage[jobId]=CWB_STAGE_STRAIN;jobFiles[jobId]=fileList[n];continue;}  
    if(fileList[n].BeginsWith(cwb_stage_input+"/cstrain_")) 
      if(CWB_STAGE_CSTRAIN>jobStage[jobId])      
        {jobStage[jobId]=CWB_STAGE_CSTRAIN;jobFiles[jobId]=fileList[n];continue;}  
    if(fileList[n].BeginsWith(cwb_stage_input+"/coherence_")) 
      if(CWB_STAGE_COHERENCE>jobStage[jobId])    
        {jobStage[jobId]=CWB_STAGE_COHERENCE;jobFiles[jobId]=fileList[n];continue;}  
    if(fileList[n].BeginsWith(cwb_stage_input+"/supercluster_")) 
      if(CWB_STAGE_SUPERCLUSTER>jobStage[jobId]) 
        {jobStage[jobId]=CWB_STAGE_SUPERCLUSTER;jobFiles[jobId]=fileList[n];continue;}  
    if(fileList[n].BeginsWith(cwb_stage_input+"/wave_")&&fileList[n].Contains(job_label)) 
      if(CWB_STAGE_LIKELIHOOD>jobStage[jobId])   
        {jobStage[jobId]=CWB_STAGE_LIKELIHOOD;jobFiles[jobId]=fileList[n];continue;}  

/*
    // Get STOP JOB info from history
    TFile *ifile = TFile::Open(fileList[n]);
    if(ifile==NULL) {cout << "Failed to open " << fileList[n].Data() << endl;exit(-1);}
    CWB::History* ihistory = (CWB::History*)ifile->Get("history");
    if(ihistory==NULL) { cout << "Error : history is not present!!!" << endl;exit(1); }
    int log_size = ihistory->GetLogSize("FULL");
    TString log = ihistory->GetLog("FULL",log_size-1);
    ifile->Close();
    if(log!="STOP JOB") nrecovery++;
*/
  }
  //for (int i=0;i<max_jobs;i++) cout << i << " " << jobStage[i] << " " << jobFiles[i].Data() << " " << cwb_stage << endl;

  int nrecovery=0;
  for (int i=0;i<max_jobs;i++)
    if ((jobStage[i]>=cwb_stage)||(jobStage[i]>=CWB_STAGE_LIKELIHOOD)) nrecovery++; 
  nrecovery=jobList.size()-nrecovery;
  if(nrecovery==0) {
    cout << "No Jobs to be recovered" << endl;
    gSystem->Exit(0);
  }
  cout << endl;
  cout << "New Recovey File " << endl;
  cout << dagfile << endl;

  // condor log dirs
  char full_condor_out_dir[1024];
  char full_condor_err_dir[1024];
  sprintf(full_condor_out_dir,"%s/%s",work_dir,log_dir);
  sprintf(full_condor_err_dir,"%s/%s",work_dir,log_dir);

  // create dag condor file
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);

  ofstream out;
  out.open(dagfile,ios::out);
  int cnt = 0;
  for (int i=0;i<max_jobs;i++) {
    if (i%1000==0) cout << i << "/" << max_jobs << endl;
    if ((jobStage[i]!=(CWB_STAGE)-1)&&(jobStage[i]<cwb_stage)&&(jobStage[i]<CWB_STAGE_LIKELIHOOD)) {
      if(cwb_stage_resume=="TRUE") if(!jobFiles[i].BeginsWith(cwb_resume_label)) continue;
      cnt++;
      char ostring[256];
      int jobID=i+1;
      sprintf(ostring,"JOB A%i %s/%s.sub.recovery.%d",jobID,full_condor_dir,data_label,iversion);
      out << ostring << endl;
      sprintf(ostring,"VARS A%i PID=\"%i\" CWB_UFILE=\"%s\" CWB_STAGE=\"%s\"",
                      jobID,jobID,jobFiles[i].Data(),cwb_stage_name.Data());
      out << ostring << endl;
      sprintf(ostring,"RETRY A%i 3000",jobID);
      out << ostring << endl;
      // remove broken symbolic links of condor log files (avoid init condor failure)
      TString path;
      char symlink[1024];
      Long_t id,size,flags,mt;
      sprintf(symlink,"%s/%d_%s_%s.out",full_condor_out_dir,jobID,data_label,cwb_stage_name.Data());
      path = CWB::Toolbox::getFileName(symlink);
      if(path!="") {
        int estat = gSystem->GetPathInfo(path.Data(),&id,&size,&flags,&mt); 
        if(estat!=0) {	// condor log out symbolic link is broken
          char cmd[1024]; sprintf(cmd,"rm -f %s",symlink);
          gSystem->Exec(cmd);
        }
      }
      sprintf(symlink,"%s/%d_%s_%s.err",full_condor_err_dir,jobID,data_label,cwb_stage_name.Data());
      path = CWB::Toolbox::getFileName(symlink);
      if(path!="") {
        int estat = gSystem->GetPathInfo(symlink,&id,&size,&flags,&mt); 
        if(estat!=0) {	// condor log err symbolic link is broken
          char cmd[1024]; sprintf(cmd,"rm -f %s",symlink);
          gSystem->Exec(cmd);
        }
      }
    }
  }
  out.close();

  if(gSystem->Getenv("_USE_LSF")!=NULL) {

    // make lsf label
    char lsf_label[1024];
    if(cwb_stage_name=="CWB_STAGE_FULL") {
      sprintf(lsf_label,"%s",data_label);
    } else {
      sprintf(lsf_label,"%s_%s",data_label,cwb_stage_name.Data());
    }

    // create tgz of the working dir
    TString exec_cmd = TString::Format("tar -czf  %s/%s.tgz  %s %s %s %s/*.sh --exclude='*/.svn'",
                                       condor_dir, lsf_label, input_dir, config_dir, macro_dir, condor_dir);
    gSystem->Exec(exec_cmd);
    cout << endl << "Created tgz file : " << condor_dir<<"/"<<lsf_label<<".tgz" << endl;

    // create LSF file from the DAG file
    TString lsfFile=TB.DAG2LSF(dagfile, data_label, nodedir, data_dir, condor_dir, log_dir, output_dir, work_dir);
    if(lsfFile!="") {
      cout << endl << "Unfinished Jobs  : " << cnt << "/" << jobList.size() << endl;
      cout << endl << "Created LSF file : " << lsfFile << endl << endl;
      cout << "To submit LSF recovered jobs, type :" << endl;
      cout << "cwb_lsf submit " << lsfFile << endl << endl;
    } else {
      cout << endl << "No jobs to be submitted !!!" << endl << endl;
    }
    gSystem->Exit(0);
  }

  if(gSystem->Getenv("_USE_PEGASUS")!=NULL) {
    // create in tgz file
    ofstream out;
    char infile[1024];
    sprintf(infile,"%s/%s.in.recovery.%d",condor_dir,data_label,iversion);
    out.open(infile,ios::out);
    out << "../" << config_dir << "/" << endl;
    out << "../" << input_dir  << "/" << endl;
    out << "../" << macro_dir  << "/" << endl;
/*
    for (int i=0;i<max_jobs;i++) {
      if ((jobStage[i]!=-1)&&(jobStage[i]<cwb_stage)&&(jobStage[i]<CWB_STAGE_LIKELIHOOD)) {
        if(jobStage[i]!=CWB_STAGE_FULL) out << "../" << jobFiles[i] << endl;
      }
    }
*/
    out.close();
    // execute cwb_pegasus_create.sh
    sprintf(dagfile,"%s.dag.recovery.%d",data_label,iversion);
    TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
    TString exec_cmd = TString::Format("cd %s;%s/cwb_pegasus_create.sh %s",
                                       condor_dir,cwb_scripts.Data(),dagfile);
    int ret=gSystem->Exec(exec_cmd);
    if(ret) {cout << "Error  while executing cwb_pegasus_create !!!" << endl;exit(1);}
 
    cout << endl;
    cout << "Unfinished Jobs : " << cnt << "/" << jobList.size() << endl;
    cout << endl;
    sprintf(dagfile,"%s/%s.dag.recovery.%d",condor_dir,data_label,iversion);
    cout << "To submit pegasus recovered jobs, type :" << endl;
    cout << "cwb_pegasus submit " << dagfile << endl;
  } else {
    // create sub condor file
    char extention[1024];
    sprintf(extention,"recovery.%d",iversion);
    TB.createSubFile(data_label, full_condor_dir, full_condor_out_dir, 
                     full_condor_err_dir, condor_log, extention, condor_tag);
    cout << endl;
    cout << "Unfinished Jobs : " << cnt << "/" << jobList.size() << endl;
    cout << endl;
    sprintf(dagfile,"%s/%s.dag.recovery.%d",condor_dir,data_label,iversion);
    cout << "To submit condor recovered jobs, type :" << endl;
    cout << "cwb_condor submit " << dagfile << endl;
  }
  cout << endl;

  if(gSystem->Getenv("_USE_LSF")!=NULL) {

    TString cwb_stage_label="supercluster_";                   
    if( cwb_stage_input=="FULL")         cwb_stage_label="wave_";
    if( cwb_stage_input=="INIT")         cwb_stage_label="init_";
    if( cwb_stage_input=="STRAIN")       cwb_stage_label="strain_";
    if( cwb_stage_input=="CSTRAIN")      cwb_stage_label="cstrain_";
    if( cwb_stage_input=="COHERENCE")    cwb_stage_label="coherence_";
    if( cwb_stage_input=="SUPERCLUSTER") cwb_stage_label="supercluster_";
    if( cwb_stage_input=="LIKELIHOOD")   cwb_stage_label="wave_";

    int jobID=1;   // jobID must be defined -> TO BE FIXED !!!!
    TString exec_cmd = TString::Format("export file_n_st=""$(ls %s*_job%i.root)""",cwb_stage_label.Data(),jobID);
    gSystem->Exec(exec_cmd);
    gSystem->Exec("echo $file_n_st");
    if(gSystem->Getenv("file_n_st")!=NULL) {
      char *file_tmp = gSystem->ExpandPathName("output/$file_n_st");
      //cout<<file_tmp<<endl;
      exec_cmd = TString::Format("tar -czf  %s/%s.tgz  %s %s %s %s %s %s --exclude='*/.svn' --exclude='%s/*' --exclude='%s/*'",
                                 condor_dir, data_label, input_dir, config_dir, output_dir, 
                                 log_dir, macro_dir, output_dir, log_dir, file_tmp);
    } else {
      // create tgz of the working dir
      exec_cmd = TString::Format("tar -czf  %s/%s.tgz  %s %s %s %s %s --exclude='*/.svn' --exclude='%s/*' --exclude='%s/*'",
                                 condor_dir, data_label, input_dir, config_dir, output_dir, 
                                 log_dir, macro_dir, output_dir, log_dir);
      gSystem->Exec(exec_cmd);
      cout << endl << "Created tgz file : " << condor_dir<<"/"<<data_label<<".tgz" << endl;
    }
  }

  delete [] jobStage;
  delete [] jobStart;
  delete [] jobStop;

  gSystem->Exit(0);
}
