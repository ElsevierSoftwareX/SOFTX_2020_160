/*
# Copyright (C) 2019 Gabriele Vedovato, Claudia Lazzaro
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


// submit condor jobs : used by the cwb_condor command

{
  CWB::Toolbox TB;
  char cmd[1024];

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  // get condor dag
  TString cwb_condor_dag="";
  if(gSystem->Getenv("CWB_CONDOR_DAG")!=NULL) {
    // user define condor tag
    cwb_condor_dag=TString(gSystem->Getenv("CWB_CONDOR_DAG"));
    //strip '/' from path
    TObjArray* token = cwb_condor_dag.Tokenize(TString('/'));
    TObjString* stoken =(TObjString*)token->At(token->GetEntries()-1);
    cwb_condor_dag = stoken->GetString();
  }

  // get pegasus site
  TString cwb_pegasus_usite="";
  if(gSystem->Getenv("CWB_PEGASUS_USITE")!=NULL) {
    // user define pegasus site
    TString cwb_pegasus_usite=TString(gSystem->Getenv("CWB_PEGASUS_USITE"));
  }

  // get workflow file (used by pegasus)
  TString cwb_condor_workflow=data_label;
  cwb_condor_workflow.ReplaceAll(".","_");  // pegasus do not accept '.'
  char workflowfile[1024];
  sprintf(workflowfile,"%s/%s",condor_dir,cwb_condor_workflow.Data());

  // get job batch command
  TString batch_cmd = "condor_submit_dag";
  TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
  if(gSystem->Getenv("_USE_PEGASUS")!=NULL) batch_cmd = cwb_scripts+"/cwb_pegasus_submit.sh";
  if(gSystem->Getenv("_USE_LSF")!=NULL) batch_cmd = cwb_scripts+"/cwb_lsf_submit.sh";

  if(gSystem->Getenv("_USE_PEGASUS")!=NULL) {

    // Check if pagasus work file already exist
    Long_t id,size,flags,mt;
    int estat = gSystem->GetPathInfo(workflowfile,&id,&size,&flags,&mt);
    if (estat==0) {
      char answer[256];
      strcpy(answer,"");
      do {
        cout << "File \"" << workflowfile << "\" already exist" << endl;
        cout << "Do you want to submit again ? (y/n) ";
        cin >> answer;
        cout << endl << endl;
      } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
      if (strcmp(answer,"n")==0) {
        exit(0);
      } else {
        sprintf(cmd,"rm %s",workflowfile);
        int ret=gSystem->Exec(cmd);
        if(ret) {cout << "Error  while executing " << cmd << " !!!" << endl;exit(1);}
      }
    }
  }

  if(cwb_condor_dag=="") {
    sprintf(cmd,"cd %s/%s;%s %s.dag %s",
            work_dir,condor_dir,batch_cmd.Data(),data_label,cwb_pegasus_usite.Data());
  } else {                  
    sprintf(cmd,"cd %s/%s;%s %s %s",
            work_dir,condor_dir,batch_cmd.Data(),cwb_condor_dag.Data(),cwb_pegasus_usite.Data());
  }
  cout << cmd << endl;
  int ret=gSystem->Exec(cmd);
  if(ret) {cout << "Error  while executing " << batch_cmd << " !!!" << endl;exit(1);}

  if(gSystem->Getenv("_USE_LSF")!=NULL) {
    cout << endl << endl;
    cout << "Your LSF jobs has been submitted" << endl;
    cout << "To monitor the jobs do                       : cwb_lsf status" << endl;
    cout << "To monitor the queue                         : cwb_lsf queue" << endl;
    cout << "To kill the all jobs do                      : cwb_lsf kill"   << endl;
    cout << "To resubmit paused jobs do                   : cwb_lsf resume" << endl;
    cout << "To suspend all jobs do                       : cwb_lsf stop"   << endl << endl;
  }

  if(gSystem->Getenv("_USE_PEGASUS")!=NULL) {
    sprintf(cmd, "cd %s;ls -1 -trd workflows/*/pegasus/cwb/* | tail -n 1 | awk 'BEGIN { OFS = \"\"; ORS = \"\" } ; {print $1}; {print \" %s\"}' | xargs ln -sf",condor_dir,data_label);
    cout << cmd << endl;
    int ret=gSystem->Exec(cmd);
    if(ret) {cout << "Error  while executing " << cmd << " !!!" << endl;exit(1);}
    sprintf(cmd, "ln -sfn %s/%s pegasus",condor_dir,data_label);
    cout << cmd << endl;
    ret=gSystem->Exec(cmd);
    cout << endl << endl;
    cout << "Your workflow has been started and is running in the base directory: " << endl;
    cout <<  workflowfile << endl << endl;
    cout << "To monitor the workflow do                   : cwb_pstatus"   << endl;
    cout << "To remove the workflow do                    : cwb_premove"   << endl;
    cout << "To resubmit an aborted or failed workflow do : cwb_prun"      << endl;
    cout << "To analyze the workflow do                   : cwb_panalyzer" << endl << endl;
  }

  gSystem->Exit(0);
}
