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


// list the jobs used for slags (fixed length) : used by the cwb_dump command

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

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

  vector<waveSegment> jobList=TB.getSlagJobList(dq1List, segLen);   // get  slag job list
  TB.dumpSegList(jobList, jobListFile);           // dump slag job list to file  (ONLY FOR DEBUG !!!)
  cout << "Dump file : " << jobListFile << endl;
  //cout << "nJob = " << jobList.size() << endl;
  xtime=TB.getTimeSegList(jobList);
  //cout << "job livetime   : " << int(xtime) << " sec " 
  //     << xtime/3600. << " h " << xtime/86400. << " day" << endl;
  cout << endl;

  exit(0);
}
