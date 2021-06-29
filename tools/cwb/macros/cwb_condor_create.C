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


// creates dag/sub files under the condor dir : used by the cwb_condor command

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  if(TString(condor_tag)=="") {
    cout << endl;
    cout << "cwb_condor_create.C : Error - the accounting_group is not defined !!!" << endl;
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

  // get cwb stage name
  TString cwb_stage_name="CWB_STAGE_FULL";
  if(gSystem->Getenv("CWB_STAGE_NAME")!=NULL) {
    cwb_stage_name=TString(gSystem->Getenv("CWB_STAGE_NAME"));
  }

  char full_condor_dir[1024];
  char full_condor_out_dir[1024];
  char full_condor_err_dir[1024];

  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);
  sprintf(full_condor_out_dir,"%s/%s",work_dir,log_dir);
  sprintf(full_condor_err_dir,"%s/%s",work_dir,log_dir);


  char dagfile[1024];
  sprintf(dagfile,"%s/%s.dag",condor_dir,data_label);

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


  if(nIFO==1) { // Single Detector Mode
    CWB::config config;
    config.Import();
    config.SetSingleDetectorMode();
    config.Export();
  }

  vector<slag> rslagList;
  vector<TString> ifos(nIFO);
  for(int n=0;n<nIFO;n++) {
    if(strlen(ifo[n])!=0) ifos[n]=ifo[n];           // built in detector
    else                  ifos[n]=detParms[n].name; // user define detector
  }

  vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, CWB_CAT1);

  // extract detector's segments range from dq cat 1
  if(slagSize==0) {
    vector<waveSegment> cat2List=TB.readSegList(nDQF, DQF, CWB_CAT2);
    // segments are standard style !!!
    vector<waveSegment> jobList=TB.getJobList(cat1List, cat2List, segLen, segMLS, segTHR, segEdge);
    cout << "Standard Segments : " << jobList.size() << endl;
    // create dag condor file
    vector<TString> jobFiles;
    TB.createDagFile(jobList, full_condor_dir, data_label, jobFiles, cwb_stage_name);
  } else {
    // segments are slag style
    int slagSegs=TB.getSlagJobList(cat1List, segLen).size();
    // full slag list
    vector<slag> slagList=TB.getSlagList(nIFO, slagSize, slagSegs, slagOff, slagMin, slagMax, slagSite, slagFile);  
    cout << "slagList size : " << slagList.size() << endl;
    // reduced slag list (removed slag with max segLength < segMSL after cat1 & segLength < segTHR after cat2)
    cout << endl << "Start segments selection from dq cat1 list ..." << endl << endl; 
    rslagList=TB.getSlagList( slagList, ifos, segLen, segMLS, segEdge, nDQF, DQF, CWB_CAT1);
    cout << "Number of selected jobs after cat1 : " << rslagList.size() << endl;
    cout << endl << "Start segments selection from dq cat2 list ..." << endl << endl; 
    rslagList=TB.getSlagList(rslagList, ifos, segLen, segTHR, segEdge, nDQF, DQF, CWB_CAT2);
    cout << "Number of selected jobs after cat2 : " << rslagList.size() << endl;
    // create dag condor file
    vector<TString> jobFiles;
    TB.createDagFile(rslagList, full_condor_dir, data_label, jobFiles, cwb_stage_name);
  }

  if(gSystem->Getenv("_USE_OSG")!=NULL) {
    // create tgz of the working dir
    TString exec_cmd = TString::Format("tar -czf  %s/%s.tgz  %s %s %s %s %s --exclude='*/.svn' --exclude='%s/*' --exclude='%s/*'",
                       condor_dir, data_label, input_dir, config_dir, output_dir, log_dir, macro_dir, output_dir, log_dir);
    gSystem->Exec(exec_cmd);
    cout << endl << "Created tgz file : " << condor_dir<<"/"<<data_label<<".tgz" << endl;
  }

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
    char dagFile[1024];
    sprintf(dagFile,"%s/%s.dag",condor_dir,data_label);
    TString lsfFile=TB.DAG2LSF(dagFile, data_label, nodedir, data_dir, condor_dir, log_dir, output_dir, work_dir); 
    if(lsfFile!="") {
      cout << endl << "Created LSF file : " << lsfFile << endl << endl;
      cout << "To submit lsf jobs, type : cwb_lsf submit" << endl;
    } else {
      cout << endl << "No jobs to be submitted !!!" << endl << endl;
    }
    gSystem->Exit(0);
  }

  if(gSystem->Getenv("_USE_PEGASUS")!=NULL) {
    // create in tgz file
    ofstream out;
    char infile[1024];
    sprintf(infile,"%s/%s.in",full_condor_dir,data_label);
    out.open(infile,ios::out);
    out << "../" << config_dir << "/" << endl;
    out << "../" << input_dir  << "/" << endl;
    out << "../" << macro_dir  << "/" << endl;
    out.close();
    // execute cwb_pegasus_create.sh 
    sprintf(dagfile,"%s.dag",data_label);
    TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
    TString exec_cmd = TString::Format("cd %s;%s/cwb_pegasus_create.sh %s",
                                       condor_dir,cwb_scripts.Data(),dagfile);
    int ret=gSystem->Exec(exec_cmd);
    if(ret) {cout << "Error  while executing cwb_pegasus_create !!!" << endl;exit(1);}

    cout << endl;
    cout << "To submit pegasus jobs, type : cwb_pegasus submit" << endl;
    cout << endl;
  } else {
    // create sub condor file
    TB.createSubFile(data_label, full_condor_dir, full_condor_out_dir, 
                     full_condor_err_dir, condor_log, "", condor_tag);

    cout << endl;
    cout << "To submit condor jobs, type : cwb_condor submit" << endl;
    cout << endl;
  }

  gSystem->Exit(0);
}
