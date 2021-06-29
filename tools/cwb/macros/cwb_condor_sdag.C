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


// read the standard *.dag file produced with cwb_condor create and produced a new dag file *.sdag
// The *.sdag permits to submit in each node N jobs in sequential mode
// This solution speedup the submission of jobs when the job execution time is too short  

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  int cwb_condor_njobs=-1;  // slag>0
  if(gSystem->Getenv("CWB_CONDOR_NJOBS")==NULL) {
    cout << "Error : environment CWB_CONDOR_NJOBS is not defined!!!" << endl;exit(1);
  }
  if(TString(gSystem->Getenv("CWB_CONDOR_NJOBS")).IsDigit()) {
    cwb_condor_njobs=TString(gSystem->Getenv("CWB_CONDOR_NJOBS")).Atoi();
  }
  if(cwb_condor_njobs<=1) {
    cout << "cwb_condor_sdag.C - Error : nJOBS must be >1" << endl;
    gSystem->Exit(1);
  }

  cout << cwb_condor_njobs << endl;

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
    cout << "cwb_condor_sdag.C - Error : dag file name must contains '.dag'" << endl;
    gSystem->Exit(1);
  }

  vector<int> jobList=TB.getCondorJobList(full_condor_dir, data_label);

  TString cwb_condor_sdag = cwb_condor_dag;
  cwb_condor_sdag.ReplaceAll(".dag",".sdag");
  cout << cwb_condor_sdag << endl;

  TB.checkFile(cwb_condor_sdag,true);

  ofstream out;
  out.open(cwb_condor_sdag.Data(),ios::out);
  if (!out.good()) {cout << "cwb_condor_sdag.C - Error Opening File : " << cwb_condor_sdag << endl;exit(1);}

  ifstream in;
  in.open(cwb_condor_dag.Data(),ios::in);
  if (!in.good()) {cout << "cwb_condor_dag.C - Error Opening File : " << cwb_condor_dag << endl;exit(1);}

  // copy dag to sdag
  char istr[1024];
  while(true) {
    in.getline(istr,1024);
    if (!in.good()) break;
    out << istr << endl;
  }
  in.close();

  // add super jobs instructions to sdag file
  out << endl;
  int nsjob=0;
  for(int i=0;i<jobList.size();i++) {
    for(int j=1;j<cwb_condor_njobs && i<jobList.size()-1;j++) {
      out << "PARENT A" << jobList[i] << " CHILD A" << jobList[i+1] << endl;
      //cout << i << " " << j << " jobId : " << jobList[i] << endl;
      i++;
    }
    out << endl;
    nsjob++;
  }

  out.close();

  cout << endl << "Created new sdag file : " << cwb_condor_sdag << endl << endl;

  cout << "Number of  job files : " << jobList.size() << endl;
  cout << "Number of sjob files : " << nsjob << endl << endl;

  exit(0);
}
