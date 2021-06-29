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


// list the standard jobs (variable length) : used by the cwb_dump command

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  // output job list file
  char ojobListFile[256];
  sprintf(ojobListFile,"%s/%s.job",dump_dir,data_label);

  cout<<endl<<"-------------------------------------------------------------------------------------"<< endl<<endl;

  vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, CWB_CAT1);
  double cat1_time=TB.getTimeSegList(cat1List);
  cout << "total cat1 livetime : " << int(cat1_time) << " sec " 
       << cat1_time/3600. << " h " << cat1_time/86400. << " day" << endl;
  cout << endl;

  vector<waveSegment> cat2List=TB.readSegList(nDQF, DQF, CWB_CAT2);
  double cat2_time=TB.getTimeSegList(cat2List);

  vector<waveSegment> jobList=TB.getJobList(cat1List, cat2List, segLen, segMLS, segTHR, segEdge);
  double job_time=TB.getTimeSegList(jobList);

  cout << endl;
  cout << "cat1      livetime (zero lag) of the standard job list : " << int(job_time) << " sec " 
       << job_time/3600. << " h " << job_time/86400. << " day" << endl;

  double job_time_cat2 = TB.getLiveTime(jobList,cat2List);
  cout << "cat1+cat2 livetime (zero lag) of the standard job list : " << int(job_time_cat2) << " sec " 
       << job_time_cat2/3600. << " h " << job_time_cat2/86400. << " day" << endl;

  cout<<endl<<"-------------------------------------------------------------------------------------"<< endl<<endl;

  cout<<"Final number of standard jobs : " << jobList.size() <<endl<<endl;
  cout<<"Dump job list : includes jobs discarted by the condition livetime<segTHR"<<endl;
  TB.dumpJobList(cat1List, ojobListFile, segLen, segMLS, segEdge);
  cout << endl;

  exit(0);
}
