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


// dump data quality duty times : used by the cwb_dump command

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  char dq1ListFile[256];  
  cout << endl;
  sprintf(dq1ListFile,"%s/%s.dq1",dump_dir,data_label);
  vector<waveSegment> dq1List=TB.readSegList(nDQF, DQF, CWB_CAT1);
  TB.dumpSegList(dq1List,dq1ListFile, false);
  double ctime_dq1=TB.getTimeSegList(dq1List);
  cout << "dq cat 1 ctime : " << int(ctime_dq1) << " sec " << ctime_dq1/3600. << " h " << ctime_dq1/86400. << " day" << endl;

  char dq2ListFile[256];  
  cout << endl;
  sprintf(dq2ListFile,"%s/%s.dq2",dump_dir,data_label);
  vector<waveSegment> dq2List=TB.readSegList(nDQF, DQF, CWB_CAT2);
  TB.dumpSegList(dq2List,dq2ListFile, false);
  double ctime_dq2=TB.getTimeSegList(dq2List);
  cout << "dq cat 2 ctime : " << int(ctime_dq2) << " sec " << ctime_dq2/3600. << " h " << ctime_dq2/86400. << " day" << endl;
  cout << endl;

  exit(0);
}
