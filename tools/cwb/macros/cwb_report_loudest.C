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


// produce condor dag file of loudest event list : used by the cwb_condor command

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  if(simulation) {
    cout << endl << "cwb_report_loudest.C : Error - "
         << "this command can not be used in simulation mode !!!" << endl << endl;
    exit(1);
  }

  if(TString(condor_tag)=="") {
    cout << endl;
    cout << "cwb_condor_loudest.C : Error - the accounting_group is not defined !!!" << endl;
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

  TString cwb_stage_name="CWB_STAGE_LIKELIHOOD";
  CWB_STAGE cwb_stage=CWB_STAGE_LIKELIHOOD;

  // extract options
  float     cwb_loudest_rho    = 0;  
  int       cwb_loudest_nevt   = pp_max_nloudest_list;  
  int       cwb_loudest_idmin  = 0;  
  int       cwb_loudest_idmax  = 0;  
  TString   cwb_loudest_ufile  = cwb_uparameters_file;
  TString   cwb_loudest_ced    = "FALSE";
  TString   cwb_loudest_veto   = "FALSE";
  TString   cwb_loudest_odir   = "loudest";
  TString   cwb_loudest_stage  = "AUTO";
  TString cwb_report_options = TString(gSystem->Getenv("CWB_REPORT_OPTIONS"));
  if(cwb_report_options!="") {
    TString option="";
    // get the stage
    option = TB.getParameter(cwb_report_options,"--stage");
    if(option!="") {
       option.ToUpper();
       if(option=="FULL") cwb_loudest_stage = option;
    }
    // get the loudest output dir
    option = TB.getParameter(cwb_report_options,"--odir");
    if(option!="") cwb_loudest_odir  = option;
    // get the min id number of the loudest events to be processed
    option = TB.getParameter(cwb_report_options,"--idmin");
    if(option!="") {
       if(option.IsDigit()) cwb_loudest_idmin = option.Atoi();
       else { cout << endl << "cwb_report_loudest.C : Error - "
                   << "the option --idmin value is not a number !!!" 
                   << endl << endl; exit(1);}
    }
    // get the max id number of the loudest events to be processed
    // if not defined then idmax=idmin
    option = TB.getParameter(cwb_report_options,"--idmax");
    if(option!="") {
       if(option.IsDigit()) cwb_loudest_idmax = option.Atoi();
       else { cout << endl << "cwb_report_loudest.C : Error - "
                   << "the option --idmax value is not a number !!!" 
                   << endl << endl; exit(1);}
    } else cwb_loudest_idmax=cwb_loudest_idmin;
    if(cwb_loudest_idmax<cwb_loudest_idmin) {
       cout << endl << "cwb_report_loudest.C : Error - "
            << "the option --idmax value is < of --idmin value !!!" << endl << endl; 
       exit(1);
    }
    // get number of loudest events to be processed
    option = TB.getParameter(cwb_report_options,"--nevt");
    if(option!="") {
       if(option.IsDigit()) cwb_loudest_nevt = option.Atoi();
       else { cout << endl << "cwb_report_loudest.C : Error - "
                   << "the option --nevt value is not a number !!!" 
                   << endl << endl; exit(1);}
    }
    if(cwb_loudest_nevt==0) cwb_loudest_nevt=100000; 
    // get the user parameters file option 
    option = TB.getParameter(cwb_report_options,"--ufile");
    if(option!="") cwb_loudest_ufile = option;
    // get the veto option 
    option = TB.getParameter(cwb_report_options,"--veto");
    if(option!="") cwb_loudest_veto = option;
    cwb_loudest_veto.ToUpper();
    // get the ced option 
    option = TB.getParameter(cwb_report_options,"--ced");
    if(option!="") cwb_loudest_ced = option;
    cwb_loudest_ced.ToUpper();
    // get the rho threshold option 
    option = TB.getParameter(cwb_report_options,"--rho");
    if(option!="") {
       if(option.IsFloat()) cwb_loudest_rho = option.Atof();
       else { cout << endl << "cwb_report_loudest.C : Error - "
                   << "the option --rho value is not a number !!!" 
                   << endl << endl; exit(1);}
    }
  }

  // creates dir for loudest events
  TString pp_loudest_dir = TString(pp_dir)+TString("/")+cwb_loudest_odir ;
  CWB::Toolbox::mkDir(pp_loudest_dir,!pp_batch);

  // read condor job list
  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);
  vector<int> jobList=TB.getCondorJobList(full_condor_dir, data_label);
  int max_jobs = 0;
  for(int i=0;i<jobList.size();i++) if(jobList[i]>max_jobs) max_jobs=jobList[i];

  // read loudest event list
  TB.checkFile(netdir);
  char events_sorted[1024];
  sprintf(events_sorted,"%s/events_sorted.txt",netdir);
  CWB::Toolbox::checkFile(events_sorted);

  char pm[8];
  char c3[8];
  float icc, icc2, icc3, irho, iacor, ilag, islag, ilik, ipen, icHH, ivHH, ivED;
  int ifreq, ilow, ihigh;
  float idur;
  int isize, irate, irun;
  float phi, theta, psi;

  double* itime = new double[NIFO_MAX];
  double* iSNR  = new double[NIFO_MAX];
  double* ihrss = new double[NIFO_MAX];

  vector<int>* lagLoudest = new vector<int>[max_jobs]; // contains the lag list of the loudest events for each job

  int iline=0;	// line counter
  int ievt=0;	// event select counter
  vector<TString> JTAG;
  ifstream f_ev(events_sorted);
  while(ievt<cwb_loudest_nevt) {
    f_ev>>pm>>c3>>irho>>icc>>icc2>>icc3>>iacor>>ilag>>islag>>ilik;
    f_ev>>ipen>>icHH>>ifreq>>ilow>>ihigh>>idur>>isize>>irate>>irun;
    for(int i=0;i<nIFO;i++) f_ev>>itime[i];                                                                                 
    for(int i=0;i<nIFO;i++) f_ev>>iSNR[i];                                                                                  
    for(int i=0;i<nIFO;i++) f_ev>>ihrss[i];                                                                                 
    f_ev>>phi>>theta>>psi;                                                                                                  
    if (!f_ev.good()) break;                                                                                                
    iline++;
    if (!TString(pm).Contains("+")) continue;                                                                               
    if ((cwb_loudest_veto=="TRUE")&&(strcmp(pm,"+")||strcmp(c3,"-"))) continue;
    if (cwb_loudest_idmin&&(iline<cwb_loudest_idmin)) continue;
    if (cwb_loudest_idmax&&(iline>cwb_loudest_idmax)) continue;
    if (irho<cwb_loudest_rho) continue;
    ievt++;                                                                                                                 
    //cout << irun << " " << 0 << " " << "true" << " " << 0 << " " << ilag << " " << ced_dir << endl;
    //out << irun << " " << 0 << " " << "true" << " " << 0 << " " << ilag << " " << ced_dir << endl;
    char jtag[1024];sprintf(jtag,"%i_%i_%i",irun,(int)islag,(int)ilag);
    bool bjtag=false;for(int j=0;j<JTAG.size();j++) if(JTAG[j]==jtag) bjtag=true;
    if(!bjtag) JTAG.push_back(jtag); else continue;     // skip job if already created
    if(irun>max_jobs) {
      cout << endl << "cwb_report_loudest.C : Error - "
           << "job number is greater than max allowed jobs, check your yours dag and user_parameter.C files" << endl << endl;
      exit(1);
    }
    lagLoudest[irun-1].push_back(ilag);
  }
  f_ev.close();

  delete [] itime;
  delete [] iSNR;
  delete [] ihrss;

  if(ievt==0) {
    cout << endl << "cwb_report_loudest.C : Warning - "
         << "no events has been selected !!!" << endl << endl;
    exit(1);
  }

  // check stages of output root files
  int* jobStart = new int[max_jobs];
  int* jobStop  = new int[max_jobs];

  CWB_STAGE* jobStage = new CWB_STAGE[max_jobs];
  for (int i=0;i<max_jobs;i++) jobStage[i]=(CWB_STAGE)-1;	                        // excluded jobs
  for (int i=0;i<jobList.size();i++) jobStage[jobList[i]-1]=(CWB_STAGE)0;               // jobs in the dag file

  char tag[1024];sprintf(tag,"%s.dag.loudest.",data_label);
  vector<TString> fileList = TB.getFileListFromDir(condor_dir, "", tag);
  int iversion=0;
  for(int i=0;i<fileList.size();i++) {
    //cout << i << " " << fileList[i].Data() << endl;
    TObjArray* token = TString(fileList[i]).Tokenize(TString("."));
    TObjString* sloudestID = (TObjString*)token->At(token->GetEntries()-1);
    if(sloudestID->GetString().IsDigit()) {
      cout << i << " " << fileList[i].Data() << endl;
      int loudestID = sloudestID->GetString().Atoi();
      if(iversion<loudestID) iversion=loudestID;
    }
  }
  iversion++;

  char ofile[1024];
  sprintf(ofile,"%s/%s.dag.loudest.%d",condor_dir,data_label,iversion);

  // Check if file exist
  Long_t id,size,flags,mt;
  int estat = gSystem->GetPathInfo(ofile,&id,&size,&flags,&mt);
  if (estat==0) {
    char answer[1024];
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

  // get cwb stage input dir (input files produced by the previous stage)
  TString cwb_stage_input=output_dir;
  if(gSystem->Getenv("CWB_STAGE_INPUT")!=NULL) {
    cwb_stage_input=TString(gSystem->Getenv("CWB_STAGE_INPUT"));
  }
  if(cwb_stage_input=="") cwb_stage_input=output_dir;
  TB.checkFile(cwb_stage_input);

  char job_label[1024];sprintf(job_label,"%s",data_label);

  cout << "Starting reading output directory ..." << endl;
  vector<TString> jobFiles(max_jobs);
  for(int i=0;i<max_jobs;i++) jobFiles[i]=cwb_uparameters_file; 
  vector<TString> fileList2 = TB.getFileListFromDir(cwb_stage_input,".root","",data_label,true);
  if(cwb_loudest_stage=="AUTO") for(int n=0;n<fileList2.size();n++) {

    int jobId = TB.getJobId(fileList2[n]);    // Get JOB ID
    jobId-=1;
     
    if(fileList2[n].BeginsWith(cwb_stage_input+"/init_")) 
      if(CWB_STAGE_INIT>jobStage[jobId]) 
        {jobStage[jobId]=CWB_STAGE_INIT;jobFiles[jobId]=fileList2[n];continue;}  
    if(fileList2[n].BeginsWith(cwb_stage_input+"/strain_")) 
      if(CWB_STAGE_STRAIN>jobStage[jobId])       
        {jobStage[jobId]=CWB_STAGE_STRAIN;jobFiles[jobId]=fileList2[n];continue;}  
    if(fileList2[n].BeginsWith(cwb_stage_input+"/cstrain_")) 
      if(CWB_STAGE_CSTRAIN>jobStage[jobId])      
        {jobStage[jobId]=CWB_STAGE_CSTRAIN;jobFiles[jobId]=fileList2[n];continue;}  
    if(fileList2[n].BeginsWith(cwb_stage_input+"/coherence_")) 
      if(CWB_STAGE_COHERENCE>jobStage[jobId])    
        {jobStage[jobId]=CWB_STAGE_COHERENCE;jobFiles[jobId]=fileList2[n];continue;}  
    if(fileList2[n].BeginsWith(cwb_stage_input+"/supercluster_")) 
      if(CWB_STAGE_SUPERCLUSTER>jobStage[jobId]) 
        {jobStage[jobId]=CWB_STAGE_SUPERCLUSTER;jobFiles[jobId]=fileList2[n];continue;}  
  }
  //for (int i=0;i<max_jobs;i++) cout << i << " " << jobStage[i] << " " << jobFiles[i].Data() << " " << cwb_stage << endl;

  // if jobFile is a configuration file and user has provided an auxiliary config file
  // then jobFile is initialized with the auxiliary config file
  for(int i=0;i<max_jobs;i++) {
    if(jobFiles[i]==cwb_uparameters_file) { 
      if(cwb_loudest_ufile!=cwb_uparameters_file) { 
        jobFiles[i]=cwb_loudest_ufile;
      }
    }
  }

  int nloudest=0;
  for (int i=0;i<max_jobs;i++) nloudest+=lagLoudest[i].size();
  if(nloudest==0) {
    cout << "No loudest events to be processed !!!" << endl;
    gSystem->Exit(0);
  }
  cout << endl;
  cout << "New Loudest File " << endl;
  cout << ofile << endl;

  // condor log dirs
  char full_condor_out_dir[1024];
  char full_condor_err_dir[1024];
  sprintf(full_condor_out_dir,"%s/%s",work_dir,log_dir);
  sprintf(full_condor_err_dir,"%s/%s",work_dir,log_dir);

  // create dag condor file
  char full_condor_dir2[1024];
  sprintf(full_condor_dir2,"%s/%s",work_dir,condor_dir);

  ofstream out;
  out.open(ofile,ios::out);
  int jID = 0;
  for (int i=0;i<max_jobs;i++) {
    for (int j=0;j<lagLoudest[i].size();j++) {
      jID++;
      char ostring[1024];
      int jobID=i+1;
      char jtag[1024];sprintf(jtag,"%i_%i",jobID,lagLoudest[i][j]);
      sprintf(ostring,"JOB A%i_%s %s/%s.sub.loudest.%d",jID,jtag,full_condor_dir2,data_label,iversion);
      out << ostring << endl;
      sprintf(ostring,"VARS A%i_%s PID=\"%i\" CWB_UFILE=\"%s\" CWB_STAGE=\"%s\" ",
                      jID,jtag,jobID,jobFiles[i].Data(),cwb_stage_name.Data());
      out << ostring;
      sprintf(ostring,"CWB_MDC_FACTOR=\"0\" CWB_JOB_LAG=\"%i\" CWB_CED_DIR=\"%s\" ",
                      lagLoudest[i][j],pp_loudest_dir.Data());
      out << ostring;
      sprintf(ostring,"CWB_INET_OPTIONS=\"%s\" CWB_UPARAMETERS_FILE=\"%s\" CWB_BATCH=\"true\"",
                      cwb_loudest_ced.Data(),cwb_loudest_ufile.Data());
      out << ostring << endl;
      sprintf(ostring,"RETRY A%i_%s 3000",jID,jtag);
      out << ostring << endl;
      // remove broken symbolic links of condor log files (avoid init condor failure)
/*    COMMENTED BECAUSE TAKES TIME
      TString path;
      char symlink[1024];
      Long_t id,size,flags,mt;
      sprintf(symlink,"%s/%d_%s_%s.out",full_condor_out_dir,jobID,data_label,cwb_stage_name.Data());
      path = CWB::Toolbox::getFileName(symlink);
      if(path!="") {
        int estat = gSystem->GetPathInfo(path.Data(),&id,&size,&flags,&mt); 
        if(estat!=0) { // condor log out symbolic link is broken
          char cmd[1024]; sprintf(cmd,"rm -f %s",symlink);
          gSystem->Exec(cmd);
        }
      }
      sprintf(symlink,"%s/%d_%s_%s.err",full_condor_err_dir,jobID,data_label,cwb_stage_name.Data());
      path = CWB::Toolbox::getFileName(symlink);
      if(path!="") {
        int estat = gSystem->GetPathInfo(symlink,&id,&size,&flags,&mt); 
        if(estat!=0) { // condor log err symbolic link is broken
          char cmd[1024]; sprintf(cmd,"rm -f %s",symlink);
          gSystem->Exec(cmd);
        }
      }
*/
    }
  }
  out.close();

  // check if loudest.sh script is present in the condor directory
  // if not the condor_dir/loudest.sh is linked to $CWB_SCRIPTS/cwb_loudest.sh 
  TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
  estat = gSystem->GetPathInfo(TString(condor_dir)+"/loudest.sh",&id,&size,&flags,&mt);
  if (estat!=0) {	// file not found
    char cmd[1024];
    sprintf(cmd,"ln -sf %s/cwb_loudest.sh %s/loudest.sh",cwb_scripts.Data(),condor_dir);
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  // create sub condor file
  char ofile_condor_sub[1024];
  sprintf(ofile_condor_sub,"%s/%s.sub.loudest.%d",full_condor_dir2,data_label,iversion);

  FILE *fP=NULL;
  if((fP = fopen(ofile_condor_sub, "w")) == NULL) {
    cout << "cwb_report_loudest.C : Error - cannot open file " << ofile_condor_sub << endl;
    exit(1);                                                                                  
  }                                                                                           
  cout << ofile_condor_sub << endl;                                                           

  fprintf(fP,"universe = vanilla\n");
  fprintf(fP,"getenv = true\n");     
  fprintf(fP,"priority = $(PRI)\n"); 
  fprintf(fP,"on_exit_hold = ( ExitCode != 0 )\n");
  fprintf(fP,"request_memory = 2000\n");           
  fprintf(fP,"executable = loudest.sh\n");             
  fprintf(fP,"job_machine_attrs = Machine\n");     
  fprintf(fP,"job_machine_attrs_history_length = 5\n");
  fprintf(fP,"requirements = target.machine =!= MachineAttrMachine1 && target.machine =!= MachineAttrMachine2 && target.machine =!= MachineAttrMachine3 && target.machine =!= MachineAttrMachine4 && target.machine =!= MachineAttrMachine5\n");          
  fprintf(fP,"environment = CWB_JOBID=$(PID);CWB_UFILE=$(CWB_UFILE);CWB_STAGE=$(CWB_STAGE);CWB_MDC_FACTOR=$(CWB_MDC_FACTOR);CWB_JOB_LAG=$(CWB_JOB_LAG);CWB_CED_DIR=$(CWB_CED_DIR);CWB_INET_OPTIONS=$(CWB_INET_OPTIONS);CWB_UPARAMETERS_FILE=$(CWB_UPARAMETERS_FILE);CWB_BATCH=$(CWB_BATCH)\n");
  if(TString(condor_tag)!="") fprintf(fP,"accounting_group = %s\n",condor_tag);
  fprintf(fP,"output = %s/$(PID)_$(CWB_JOB_LAG)_%s_%s.out\n",full_condor_out_dir,data_label,cwb_loudest_odir.Data());
  fprintf(fP,"error = %s/$(PID)_$(CWB_JOB_LAG)_%s_%s.err\n",full_condor_err_dir,data_label,cwb_loudest_odir.Data());
  fprintf(fP,"log = %s/%s_%s.log\n",condor_log,data_label,cwb_loudest_odir.Data());
  fprintf(fP,"notification = never\n");                                                                                      
  fprintf(fP,"rank=memory\n");                                                                                               
  fprintf(fP,"queue\n");                                                                                                     
  fclose(fP);


  cout << "Number of Jobs : " << jID << "/" << jobList.size() << endl;
  cout << endl;
  cout << "To submit condor recovered jobs type :" << endl;
  sprintf(ofile,"%s/%s.dag.loudest.%d",condor_dir,data_label,iversion);
  cout << "cwb_condor submit " << ofile << endl;
  cout << endl;

  delete [] lagLoudest;
  delete [] jobStart;
  delete [] jobStop;
  delete [] jobStage;
 
  gSystem->Exit(0);
}
