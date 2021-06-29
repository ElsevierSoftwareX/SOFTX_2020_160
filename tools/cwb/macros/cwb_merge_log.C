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


// merge log injection files : this script is called from cwb_merge.C 
// iversion is defined in cwb_merge.C
// WARNING!!!  -> For unamed macro all variables must be declared before any return statement !!! 

{

  #include <vector>

  bool first=true;

  char _cmd[1024];
  char _log[1024];
  char logFile[1024];
  char lstFile[1024];

  vector<TString> logList;
  vector<TString> _fileList;

  ifstream in;
  ofstream _out;

  TString sortFile;
  TString rootFile;
  TString logHeaderFile;

  detector* pD[NIFO_MAX];                                     //! pointers to detectors
  for(int i=0; i<nIFO; i++) {
    if(strlen(ifo[i])>0) pD[i] = new detector(ifo[i]);        // built in detector
    else                 pD[i] = new detector(detParms[i]);   // user define detector
  }
  CWB::mdc MDC(nIFO,pD);

  sprintf(logFile,"%s/log_%s.M%d.txt",merge_dir,data_label,iversion);

  // read log file list
  _fileList = CWB::Toolbox::getFileListFromDir(output_dir, ".txt","log_","",true);
  if(_fileList.size()==0) { // no log file are present in the directory
    cout << "No Log files are present on the directory : " << output_dir << endl;
    return;	
  }
  for(int j=0;j<_fileList.size();j++) {
    if(_fileList[j].Contains(data_label)) {
      logList.push_back(_fileList[j].Data());
      //cout << j << " " << _fileList[j].Data() << endl;
    }
  }
  if(logList.size()==0) { // no log file are present in the directory
    cout << "No Log files are present on the directory : " << output_dir << endl;
    return;	
  }

  // write lstFile to nodedir
  sprintf(lstFile,"%s/%s-Log.lst",tmp_dir,data_label);
  _out.open(lstFile,ios::out);
  if (!_out.good()) {cout << "cwb_merge_log.C - Error Opening File : " << lstFile << endl;exit(1);}
  for(int i=0;i<logList.size();i++) _out << logList[i].Data() << endl;
  _out.close();

  // Sort lstFile
  sortFile = lstFile;
  sortFile.ReplaceAll(".lst",".sort");
  sprintf(_cmd,"sort %s > %s",lstFile,sortFile.Data());
  cout << _cmd << endl;
  gSystem->Exec(_cmd);

  // write merged log
  in.open(sortFile.Data(),ios::in);
  if (!in.good()) {cout << "cwb_merge_log.C - Error Opening Sorted File : " << sortFile.Data() << endl;exit(1);}
  while(true) {
    in.getline(_log,1024);
    if (!in.good()) break;
    if(first) {
      sprintf(_cmd,"cat %s >  %s",_log,logFile);
      first=false;
    } else {
      sprintf(_cmd,"cat %s >> %s",_log,logFile);
    }
    //cout << _cmd << endl;
    gSystem->Exec(_cmd);
  }
  in.close();

  // remove lstFile
  sprintf(_cmd,"rm %s %s",lstFile,sortFile.Data());
  //cout << _cmd << endl;
  gSystem->Exec(_cmd);

  cout << endl << "Merged LogFile : " << logFile << endl << endl;

  // write log root
  rootFile = logFile;
  rootFile.ReplaceAll(".txt",".root");
  //cout << TString(rootFile).Data() << endl;
  MDC.SetSkyDistribution(MDC_LOGFILE,logFile,0);
  MDC.DumpLog(rootFile);

  // write log header
  logHeaderFile = logFile;
  logHeaderFile.ReplaceAll("log_","hlog_");
  MDC.DumpLogHeader(logHeaderFile,data_label);

  // merge log and log header
  sprintf(_cmd,"cat %s >> %s",logFile,logHeaderFile.Data());
  //cout << _cmd << endl;
  gSystem->Exec(_cmd);
  sprintf(_cmd,"rm %s",logFile);
  //cout << _cmd << endl;
  gSystem->Exec(_cmd);
  sprintf(_cmd,"mv %s %s",logHeaderFile.Data(),logFile);
  //cout << _cmd << endl;
  gSystem->Exec(_cmd);
}
