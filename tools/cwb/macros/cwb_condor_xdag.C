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


// This macro is used to generate the dag file when CWB_Plugin_Periodic_Frames.C is used
// Selects jobs with range overlapping 2 periods -> do not work with CWB_Plugin_Periodic_Frames.C
// Read the standard *.dag file produced with cwb_condor create and produced a new dag file *.xdag
// Creates also *.idag which permits to run the selected jobs in interactive mode

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  if(slagSize != 1) {
    cout << "cwb_condor_xdag.C - Error : This procedure do not works with segments with variable lenght (slagSize must be =1)!!!" << endl;exit(1);
  }

  int cwb_condor_fstart=-1;  
  if(gSystem->Getenv("CWB_CONDOR_FSTART")==NULL) {
    cout << "cwb_condor_xdag.C - Error : environment CWB_CONDOR_FSTART is not defined!!!" << endl;exit(1);
  }
  if(TString(gSystem->Getenv("CWB_CONDOR_FSTART")).IsDigit()) {
    cwb_condor_fstart=TString(gSystem->Getenv("CWB_CONDOR_FSTART")).Atoi();
  }
  if(cwb_condor_fstart<=1) {
    cout << "cwb_condor_xdag.C - Error : FRAME GPS START must be defined" << endl;
    gSystem->Exit(1);
  }
  cout << cwb_condor_fstart << endl;

  int cwb_condor_flen=-1;  
  if(gSystem->Getenv("CWB_CONDOR_FLEN")==NULL) {
    cout << "Error : environment CWB_CONDOR_FLEN is not defined!!!" << endl;exit(1);
  }
  if(TString(gSystem->Getenv("CWB_CONDOR_FLEN")).IsDigit()) {
    cwb_condor_flen=TString(gSystem->Getenv("CWB_CONDOR_FLEN")).Atoi();
  }
  if(cwb_condor_flen<=1) {
    cout << "cwb_condor_xdag.C - Error : FRAME LENGHT must be defined" << endl;
    gSystem->Exit(1);
  }
  cout << cwb_condor_flen << endl;

  // get list start, stop for each job
  double xtime;
  char jobListFile[256];
  sprintf(jobListFile,"%s/%s.sjob",dump_dir,data_label);
  char dq1ListFile[256];
  sprintf(dq1ListFile,"%s/%s.cat1",dump_dir,data_label);

  cout<<endl<<"-------------------------------------------------------------------------------------"<< endl<<endl;

  CWB_CAT dqcat = CWB_CAT1;
  vector<waveSegment> dq1List=TB.readSegList(nDQF, DQF, dqcat);
  TB.dumpSegList(dq1List,dq1ListFile, false);
  cout << "Dump file : " << dq1ListFile << endl;
  xtime=TB.getTimeSegList(dq1List);
  cout << "cat1 livetime : " << int(xtime) << " sec "
       << xtime/3600. << " h " << xtime/86400. << " day" << endl;

  cout<<endl<<"-------------------------------------------------------------------------------------"<< endl<<endl;

  // find jobs with range overlapping 2 periods -> do not work with CWB_Plugin_Strain_Frames.C
  vector<waveSegment> jobSeg=TB.getSlagJobList(dq1List, segLen);   // get  slag job list
  cout << "Dump file : " << jobListFile << endl;
  cout << "nJob = " << jobSeg.size() << endl;
  cout << endl;
  vector<bool> jobSelected(jobSeg.size());
  for(int i=0;i<jobSeg.size();i++) {
     double frame_offset = fmod(cwb_condor_fstart, cwb_condor_flen);
     double seg_offset = fmod(jobSeg[i].start-frame_offset, cwb_condor_flen);
     double seg_len    = jobSeg[i].stop-jobSeg[i].start;
     //if(seg_offset<segEdge || (cwb_condor_flen-seg_offset)<seg_len+2*segEdge) 
     //  cout << i+1 << " " << (int)seg_offset << " " << cwb_condor_flen-(int)seg_offset << endl;
     if(seg_offset<segEdge || (cwb_condor_flen-seg_offset)<seg_len+2*segEdge) jobSelected[i]=false; else jobSelected[i]=true;
  }

  // get condor dag
  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);
  TString cwb_condor_dag = TString::Format("%s/%s.dag",condor_dir,data_label);
  TB.checkFile(cwb_condor_dag,false);
  if(gSystem->Getenv("CWB_CONDOR_DAG")!=NULL) {
    // user define condor tag
    TString cwb_condor_dag=TString(gSystem->Getenv("CWB_CONDOR_DAG"));
  }
  if(!cwb_condor_dag.Contains(".dag")) {
    cout << "cwb_condor_xdag.C - Error : dag file name must contains '.dag'" << endl;
    gSystem->Exit(1);
  }

  vector<int> jobList=TB.getCondorJobList(full_condor_dir, data_label);

  TString cwb_condor_xdag = cwb_condor_dag;
  cwb_condor_xdag.ReplaceAll(".dag",".xdag");
  cout << cwb_condor_xdag << endl;
  //TB.checkFile(cwb_condor_xdag,false);

  // creates interactive job list
  TString cwb_condor_idag = cwb_condor_dag;
  cwb_condor_idag.ReplaceAll(".dag",".idag");
  cout << cwb_condor_idag << endl;
  //TB.checkFile(cwb_condor_idag,false);

  ofstream xout;
  xout.open(cwb_condor_xdag.Data(),ios::out);
  if (!xout.good()) {cout << "cwb_condor_xdag.C - Error Opening File : " << cwb_condor_xdag << endl;exit(1);}

  ofstream iout;
  iout.open(cwb_condor_idag.Data(),ios::out);
  if (!iout.good()) {cout << "cwb_condor_xdag.C - Error Opening File : " << cwb_condor_idag << endl;exit(1);}

  ifstream in;
  in.open(cwb_condor_dag.Data(),ios::in);
  if (!in.good()) {cout << "cwb_condor_dag.C - Error Opening File : " << cwb_condor_dag << endl;exit(1);}

  // copy dag to xdag
  char istr[1024];
  int njob=jobList.size()-1;
  int nxjob=0;		// index of selected jobs
  int nlines=1;
  while(true) {
    in.getline(istr,1024);
    if (!in.good()) break;
    if(njob<=0) break;
    int jobId = jobList[njob];
    if(jobSelected[jobId-1]) xout << istr << endl;
    if(jobSelected[jobId-1] && nlines%3==2) {
      iout << "cwb_inet " << jobId << endl;
      nxjob++;
    }
    if(nlines%3==0) njob--;
    nlines++;
  }
  in.close();
  xout.close();
  iout.close();
  cout << endl << "Created new xdag file : " << cwb_condor_xdag << endl << endl;
  cout << endl << "Created new idag file : " << cwb_condor_idag << endl << endl;

  cout << "Number of  job files : " << jobList.size() << endl;
  cout << "Number of xjob files : " << nxjob << endl << endl;

  exit(0);
}
