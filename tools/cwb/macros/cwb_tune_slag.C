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
  #include <vector>

  CWB::Toolbox TB;
  vector<slag> slagList;

  slagFile = new char[256];
  sprintf(slagFile,"report/%s.slg",data_label);

  vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, CWB_CAT1);
  int slagSegs=TB.getSlagJobList(cat1List, segLen).size();

  slagList=TB.getSlagList(nIFO, slagSize, slagSegs, slagOff, slagMin, slagMax, slagSite);  
  slagList=TB.dumpSlagList(slagList, slagFile);  
 
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
  cout << "Total number of slag jobs : " << slagList.size() << endl;
  cout << endl;

  exit(0);
}
