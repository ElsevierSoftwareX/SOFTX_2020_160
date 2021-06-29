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


{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  int cwb_jobid=0;
  if(gSystem->Getenv("CWB_JOBID")==NULL) {
    cout << "Error : environment CWB_JOBID is not defined!!!" << endl;exit(1);
  } else {
    if(TString(gSystem->Getenv("CWB_JOBID")).IsDigit()) {
      cwb_jobid=TString(gSystem->Getenv("CWB_JOBID")).Atoi();
    } else {
      cout << "Error : environment CWB_JOBID is not defined!!!" << endl;exit(1);
    }
  }
  cout << "cwb_jobid : " << cwb_jobid << endl;

  TString cwb_ifo="";
  if(gSystem->Getenv("CWB_IFO")==NULL) {
    cout << "Error : environment CWB_IFO is not defined!!!" << endl;exit(1);
  } else {
    cwb_ifo=TString(gSystem->Getenv("CWB_IFO"));
  }
  cout << "cwb_ifo : " << cwb_ifo.Data() << endl;

  TString home_frdisplay="";
  if(gSystem->Getenv("HOME_FRDISPLAY")==NULL) {
    cout << "Error : environment HOME_FRDISPLAY is not defined!!!" << endl;exit(1);
  } else {
    home_frdisplay=TString(gSystem->Getenv("HOME_FRDISPLAY"));
  }
  cout << "home_frdisplay : " << home_frdisplay.Data() << endl;

  // check if ifo is declared in user_parameters.C
  bool icheck=false;
  for(int i=0;i<nIFO;i++) {
    if(cwb_ifo.CompareTo(ifo[i])==0) icheck=true;
  }
  if(!icheck) {
    cout << "Error - ifo : " << cwb_ifo.Data() << " is not declared in " 
         << gSystem->Getenv("CWB_UPARAMETERS_FILE") << endl;
    cout << endl << "List of allowed ifos " << endl << endl; 
    for(int i=0;i<nIFO;i++) cout << ifo[i] << endl;
    cout << endl;
    exit(1);
  }

  int ifoID=0;
  for(int i=0;i<nIFO;i++) if(cwb_ifo.CompareTo(ifo[i])==0) ifoID=i;
  cout << "ifoID : " << ifoID << endl;
 
  vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, CWB_CAT1);
  vector<waveSegment> jobList;
  if(slagSize==0) {
    // get  standard job list
    jobList=TB.getJobList(cat1List, segLen, segMLS, segEdge);        
  } else {
    // get  slag job list
    jobList=TB.getSlagJobList(cat1List, segLen);   
  }

  int job_start = jobList[cwb_jobid].start;
  int job_stop  = jobList[cwb_jobid].stop;
  
  cout << "job_start : " << job_start << " job_stop : " << job_stop << endl;
 
  // prepare frame list format for frdisplay
  gRandom->SetSeed(0);
  int baudline_rnID = gRandom->Uniform(0,10000000);   // random name ID

  UserGroup_t* uinfo = gSystem->GetUserInfo();

  char cmd[256];
  sprintf(cmd,"mkdir -p /dev/shm/%s",uinfo->fUser.Data());
  gSystem->Exec(cmd);

  char baudline_FFL[512];
  sprintf(baudline_FFL,"/dev/shm/%s/cwb_frdisplay_%d.ffl",uinfo->fUser.Data(),baudline_rnID);
  cout << "baudline_FFL : " << baudline_FFL << endl;

  ofstream out;
  out.open(baudline_FFL,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << baudline_FFL << endl;exit(1);}
  
  ifstream in;
  in.open(frFiles[ifoID],ios::in);
  if (!in.good()) {cout << "Error Opening File : " << frFiles[ifoID] << endl;exit(1);}

  TString pfile_path="";
  char istring[1024];
  while (1) {
    in >> istring;
    if (!in.good()) break;
    TString file_path = istring;
    file_path.ReplaceAll("file://localhost",""); 

    TString file_path_tmp=file_path;
    file_path_tmp.ReplaceAll(".gwf","");
    TObjArray* token = TString(file_path_tmp).Tokenize(TString("-"));
    int frfile_start = ((TObjString*)token->At(token->GetEntries()-2))->GetString().Atoi();
    int frfile_len   = ((TObjString*)token->At(token->GetEntries()-1))->GetString().Atoi();
    int frfile_stop = frfile_start+frfile_len;
    if(frfile_stop<job_start) continue;
    if(frfile_start>job_stop) continue;
    out << file_path.Data() << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
    cout << file_path.Data() << endl;
  }
  in.close();
  out.close();

  cout << "out file : " << baudline_FFL << endl;

  sprintf(cmd,"%s/FrDisplay -d 5 -proc -t %s -i %s -k \"-Bu -Hp -o 6 -a 50\"",home_frdisplay.Data(),channelNamesRaw[ifoID],baudline_FFL);
  cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"ps | grep baudline | awk '{print $1}' | xargs kill -9");
  cout << cmd << endl;
  gSystem->Exec(cmd);

  exit(0);
}
