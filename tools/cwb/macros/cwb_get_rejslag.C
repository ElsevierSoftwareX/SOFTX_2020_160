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

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));


  // create super lag list
  CWB::Toolbox slagTB;		       
  CWB_CAT dqcat=CWB_CAT1;
  vector<TString> ifos(nIFO);
  for(int n=0;n<nIFO;n++) ifos[n]=ifo[n];

  vector<waveSegment> cat1List=slagTB.readSegList(nDQF, DQF, CWB_CAT1);
  int slagSegs=slagTB.getSlagJobList(cat1List, segLen).size();

  vector<slag> slagList=slagTB.getSlagList(nIFO, slagSize, slagSegs, slagOff, slagMin, slagMax, slagSite);
  vector<slag> rslagList=slagTB.getSlagList(slagList, ifos, segLen, segMLS, segEdge, nDQF, DQF, dqcat);
  
  cout << "Rejected slags = " << slagList.size()-rslagList.size() << "/" << slagList.size() << endl;
  exit(0);
}
