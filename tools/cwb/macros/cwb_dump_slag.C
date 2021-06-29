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


// list superlag used in production : used by the cwb_dump command

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  // check input user configuration
  CWB::config cfg;
  cfg.Import();
  cfg.Check();

  vector<slag> slagList;
  vector<slag> rslagList;

  vector<TString> ifos(nIFO);
  for(int n=0;n<nIFO;n++) ifos[n]=ifo[n];

  char* lagslagFile = new char[256];
  sprintf(lagslagFile,"%s/%s.lagslag",dump_dir,data_label);

  vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, CWB_CAT1);
  int slagSegs=TB.getSlagJobList(cat1List, segLen).size();

  slagList=TB.getSlagList(nIFO, slagSize, slagSegs, slagOff, slagMin, slagMax, slagSite, slagFile);

  rslagList=TB.getSlagList(slagList, ifos, segLen, segMLS, segEdge, nDQF, DQF, CWB_CAT1);
  cout << endl << "Number of selected jobs after cat1 : " << rslagList.size() << endl << endl;
  rslagList=TB.getSlagList(rslagList, ifos, segLen, segTHR, segEdge, nDQF, DQF, CWB_CAT2);
  cout << endl << "Number of selected jobs after cat2 : " << rslagList.size() << endl << endl;

  TB.dumpSlagList(rslagList, lagslagFile);

  // dump only slag
  char* _slagFile = new char[256];
  sprintf(_slagFile,"%s/%s.slag",dump_dir,data_label);
  TB.dumpSlagList(slagList, _slagFile, true);

  cout << endl;
  printf("%14s ","SLAG");
  for(int n=0; n<nIFO; n++) printf("%8sifo[%d]","",n);
  printf("\n");
  for(int i=0;i<slagList.size();i++) {
    if(slagList[i].slagId[1]==1) {
      printf("%14d", slagList[i].slagId[0]);
      for (int n=0; n<nIFO; n++) printf("%14d",slagList[i].slagId[n+2]);
      printf("\n");
    }
  }
  cout << endl;
  cout << "Write lag+slag list : " << lagslagFile << endl;
  cout << "Write     slag list : " << _slagFile << endl;
  cout << endl;

  exit(0);
}
