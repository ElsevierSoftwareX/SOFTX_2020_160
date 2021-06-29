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


// cwb_xnet : interface to the cwb class for 1G/2G pipelines
// 
// fName : for 1G/(2G jstage=1) analysis it is the user_parameters.C file name
//         for 2G jstage=2 analysis :
//            if fName is the (1G jstage=1) output cluster file name 
//               then runID is retrieved from file name  
// jstage : for 1G jstage=CWB_STAGE_FULL
//          for 2G jstage=CWB_STAGE_FULL / CWB_STAGE_INIT / CWB_STAGE_STRAIN / CWB_STAGE_CSTRAIN
//                        CWB_STAGE_COHERENCE / CWB_STAGE_SUPERCLUSTER / CWB_STAGE_LIKELIHOOD            
// uName : used only for (2G jstage=2) -> is an auxiliary user_parameters.C file
// batch : true/false -> batch/interactive mode
// ced   : true/false -> if true the easy ced script cwb_eced.C is executed (used only for easy ced)
//                       

#include <unistd.h>

void SystemExec(char* command, int maxtry=3);
void GetJobFromGPS();

void 
cwb_xnet(TString fName, CWB_STAGE jstage=CWB_STAGE_FULL, TString uName="", bool batch=false, bool eced=false) {

  GetJobFromGPS(); // if srunID="@" then jobID is automatically computed using the environmental CWB_GPS_EVENT

  int  runID=0;
  char cmd[1024];

  char tmpOut[1024];	// temporary out log file name
  char outLog[1024];	// final     out log file name
  char tmpErr[1024];	// temporary err log file name
  char errLog[1024];	// final     err log file name
  if(batch) {	
    // job stdout is redirect to local file when job is executed in batch mode (condor)
    // this is done to avoid network traffic between the execution node and the working user directory
    // when the job is finished the log file is moved to the log working directory

    // import nodedir from CINT (used in batch mode in multistages analysis)
    char nodedir_cint[1024];
    TGlobal* GLOBAL = (TGlobal*)gROOT->GetGlobal("nodedir",true); 
    if(GLOBAL!=NULL)  
      memcpy((void*)&nodedir_cint,(void*)GLOBAL->GetAddress(),1024*sizeof(char));
    // makes the nodedir_cint readable from the group
    sprintf(cmd,"chmod g+r %s",nodedir_cint);
    cout << cmd << endl;
    SystemExec(cmd);

    unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files

    TString srunID = TString(gSystem->Getenv("CWB_JOBID"));	// get jobID
    TString slabel = gSystem->BaseName(work_dir);		// get label
    TString sstage = cwb::GetStageString(jstage);          	// get stage name

    // create tmpOut & outLog
    sprintf(tmpOut,"%s/%s_%s_%s_%d.out", nodedir_cint,srunID.Data(),slabel.Data(),sstage.Data(),Pid);
    sprintf(outLog,"%s/%s/%s_%s_%s.out", work_dir,log_dir,srunID.Data(),slabel.Data(),sstage.Data());
    cout << tmpOut << endl; cout << outLog << endl;

    // create tmpErr & errLog
    sprintf(tmpErr,"%s/%s_%s_%s_%d.err", nodedir_cint,srunID.Data(),slabel.Data(),sstage.Data(),Pid);
    sprintf(errLog,"%s/%s/%s_%s_%s.err", work_dir,log_dir,srunID.Data(),slabel.Data(),sstage.Data());
    cout << tmpErr << endl; cout << errLog << endl;

    // NOTE: The creation of symbolic links has been removed because in general the working directory in the remote node is not visible from head node 
    //       for condor batch systems the output and error log files can be see at run time using the remote access to the the node using condor_ssh_to_job
    //       only at the end of jobs the output and error log files are copied to the user log directory in the head node 

/*  
    if(gSystem->Getenv("_USE_PEGASUS")==NULL) {
      // create in log dir a symbolic link to the temporary log files
      // the remove is necessary because the command 'ln' fails if the link points to a unavailable node
      // the command touch it is necessary because condor fails if xxxLog not exist
      sprintf(cmd,"rm -f %s;ln -sf %s %s;touch %s",outLog,tmpOut,outLog,outLog);
      cout << cmd << endl;
      SystemExec(cmd);
      sprintf(cmd,"rm -f %s;ln -sf %s %s;touch %s",errLog,tmpErr,errLog,errLog);
      cout << cmd << endl;
      SystemExec(cmd);
    }
*/

    freopen(tmpOut,"w",stdout);  	// redirect stdout to file
    freopen(tmpErr,"w",stderr);  	// redirect stderr to file

    // Print CWB logo
    PrintLogoCWB(GetLALVersion());
  }

  if(TString(analysis)=="1G") {
    cwb1G CWB(fName);
    if(!batch) {
      CWB::config* cfg = CWB.GetConfig();
      if(eced) cfg->Import(gSystem->ExpandPathName("$CWB_MACROS/cwb_eced.C"));  // easy ced
      cfg->Import(gSystem->ExpandPathName("$CWB_MACROS/cwb_inet.C"));
      //runID=1;
    }
    //cfg->Print();
    CWB.run(runID);
  } else
  if(TString(analysis)=="2G") {
    CWB::Toolbox TB;

    // import nodedir from CINT (used in batch mode in multistages analysis)
    char nodedir_cint[1024];
    TGlobal* GLOBAL = (TGlobal*)gROOT->GetGlobal("nodedir",true); 
    if(GLOBAL!=NULL)  
      memcpy((void*)&nodedir_cint,(void*)GLOBAL->GetAddress(),1024*sizeof(char));

    // create cwb2G object
    cwb2G CWB(fName,uName,jstage);

    // set configuration
    CWB::config* cfg = CWB.GetConfig();
    if(!batch) {
      if(eced) gROOT->Macro(gSystem->ExpandPathName("$CWB_MACROS/cwb_eced.C"));  // easy ced
      cfg->Import(gSystem->ExpandPathName("$CWB_MACROS/cwb_inet.C"));
      //runID=1;
    } else {
      // nodedir must be updated because in multistage the nodedir change
      sprintf(cfg->nodedir,"%s",nodedir_cint);	
    }

    // set save options for job file
    CWB.SetupStage(jstage);
    //cfg->Print();

    // run analysis
    CWB.run(runID);
  } else {
    cout << "cwb_xnet - Error : analysis must be 1G or 2G" << endl;
    gSystem->Exit(1);
  }

  if(batch) {	
    if(ftell(stdout)!=-1) fclose(stdout);  // close stdout
    if(ftell(stderr)!=-1) fclose(stderr);  // close stderr

    // move condor temporary log files to the working log directory 
    // if size>1000000 it is not copied (to avoid network traffic)
    int estat;
    Long_t id,size,flags,mt;

    estat = gSystem->GetPathInfo(tmpOut,&id,&size,&flags,&mt);
    if (estat==0 && size<1000000) {
      sprintf(cmd,"mv %s %s",tmpOut,outLog);
      cout << cmd << endl;
      SystemExec(cmd);
    }

    estat = gSystem->GetPathInfo(tmpErr,&id,&size,&flags,&mt);
    if (estat==0 && size<1000000) {
      sprintf(cmd,"mv %s %s",tmpErr,errLog);
      cout << cmd << endl;
      SystemExec(cmd);
    }
  }

  gSystem->Exit(0);
}

void
SystemExec(char* command, int maxtry) {
//
// Execute System Command
//
// command : system command
// maxtry  : number of retry before to exit with error (default=3)
//

  int ntry=0;
  int ecommand=1;
  while(ecommand&&(ntry<maxtry)) {
    ecommand=gSystem->Exec(command);
    if(ecommand) gSystem->Sleep(int(gRandom->Uniform(1000,3000)));  // wait [1,3] sec
    ntry++;
  }
  if(ecommand) {cout << "cwb_xnet.C - Error -> " << command << endl;gSystem->Exit(1);}
  return;
}

void 
GetJobFromGPS() {

  TString srunID = TString(gSystem->Getenv("CWB_JOBID"));	// get jobID
  TString slabel = gSystem->BaseName(work_dir);			// get label

  if(srunID!="@") return;

  int gps_event=0;
  TString cwb_gps_event=TString(gSystem->Getenv("CWB_GPS_EVENT"));
  if(cwb_gps_event.CompareTo("")!=0) {
    if(!cwb_gps_event.IsFloat()) {cout<< "Error : CWB_GPS_EVENT is not a number" << endl;exit(1);}
    if(cwb_gps_event.Atoi()>0) gps_event=cwb_gps_event.Atoi();
  }

  TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
  TString exec_cmd = TString::Format("%s/cwb_dump.csh sjob",cwb_scripts.Data());
  int ret=gSystem->Exec(exec_cmd);
  if(ret) {cout << "Error  while executing cwb_dum sjob !!!" << endl;exit(1);}

  TString sjob_file = TString("report/dump/")+slabel+".sjob";

  CWB::Toolbox::checkFile(sjob_file);

  // Open chunk list
  ifstream in;
  in.open(sjob_file.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << sjob_file << endl;exit(1);}

  int job=0;
  double start,stop;
  while(true) {
    in >> start >> stop;
    if(!in.good()) break;
    job++;
    if(gps_event>=start && gps_event<=stop) {
      TString sjob = TString::Format("%d",job);
      gSystem->Setenv("CWB_JOBID",sjob.Data());
    }
  }

  return;
}

