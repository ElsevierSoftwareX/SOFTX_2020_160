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


// This macro is used to generate periodic DQ when CWB_Plugin_Periodic_Frames.C is used

void MakePeriodicDQ(TString dq_idir, double period, int nperiod=10, TString dq_odir="") {
//
// dq_idir: input Data Quality directory
// period:  period to generate periodc DQ
// nperiod: number of periods to be generated   
// dq_odir: output Data Quality directory

  #include <vector>

  CWB::Toolbox TB;

  vector<TString> dqfList = TB.getFileListFromDir(dq_idir, ".txt","","_cat");

  dqfile* DQF = new dqfile[dqfList.size()];
  for(int i=0;i<dqfList.size();i++) {
    cout << i << " MakePeriodicDQ - Input File: " << dqfList[i] << endl;
    vector<waveSegment> dqList = TB.readSegments(dqfList[i]);
//    cout<<endl<<"-------------------------------------------------------------------------------------"<< endl;
//    for(int j=0;j<dqList.size();j++) cout << j << " " << (int)dqList[j].start << " " << (int)dqList[j].stop << endl;

    vector<waveSegment> xdqList = dqList;
    for(int j=0;j<dqList.size();j++) {
      for(int k=-nperiod; k<=nperiod; k++) {
        waveSegment seg = dqList[j];
        seg.start += k*period;
        seg.stop  += k*period;
        xdqList.push_back(seg);
      }
    }
    xdqList = TB.unionSegments(xdqList);
//    cout<<endl<<"-------------------------------------------------------------------------------------"<< endl;
//    for(int j=0;j<xdqList.size();j++) cout << j << " " << (int)xdqList[j].start << " " << (int)xdqList[j].stop << endl;

    if(dq_odir!="") {
      TString fName = TString::Format("%s/%s",dq_odir.Data(),gSystem->BaseName(dqfList[i]));
      TB.dumpSegList(xdqList, fName, false);
    }
  }

  exit(0);
}
