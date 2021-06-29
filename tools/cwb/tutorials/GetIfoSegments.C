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

//
// Read livetime tree from input root file  and produce segments.txt -> used for DetChar
//
// Example: root -l -b 'GetIfoSegments.C("merge/live_O3_K19_C01_LV_IMBHB_BKG_xrun1.M1.root",1,"K19_segments.txt")'
//

void GetIfoSegments(TString liveFileName, int ifoID, TString ofile) {
//
// liveFileName : input livetime root file
// ifoID        : detectior ID
// ofile	: output file name
//

// Read livetime ROOT file

  TFile *_file0 = TFile::Open(liveFileName);
  TTree* tree = (TTree *) gROOT->FindObject("liveTime");
  if(tree==NULL) {cout << "ScanLIVE : liveTime tree not found !!!" << endl;exit(1);}
  int nentries = tree->GetEntries();
  cout << endl << "in nentries = " << nentries << endl << endl;

  if(ifoID<0 || ifoID>=NIFO_MAX) {cout << "Error - input ifoID not correct, max available ID is: " << NIFO_MAX-1 << endl; exit(1);}

// Read start/stop entries from livetime ROOT file
// and fill the vector segmnts list iseg

  int index=0;
  waveSegment seg;
  vector<waveSegment> iseg;

  double* start = new double[NIFO_MAX];
  tree->SetBranchAddress("start",start);
  double* stop = new double[NIFO_MAX];
  tree->SetBranchAddress("stop",stop);

  for(int i=0;i<nentries;i++) {                     // loop over the detected events
    tree->GetEntry(i);                          // load entry #i
    seg.index=index++; seg.start=start[ifoID]; seg.stop=stop[ifoID]; iseg.push_back(seg);
  }

  // --------------------------------------------------------------
  // merge segments
  // --------------------------------------------------------------
  vector<waveSegment> oseg = CWB::Toolbox::unionSegments(iseg);
  cout << endl << "out size = " << oseg.size() << endl << endl;

//  cout << "List of output segments" << endl;
//  for(int n=0;n<oseg.size();n++) cout << n << " " << oseg[n].start << " " << oseg[n].stop << endl;

  // --------------------------------------------------------------
  // Write output file segment list
  // --------------------------------------------------------------
  ofstream out;
  out.open(ofile.Data(),ios::out);
  out.precision(16);
  double liveTime=0.;
  for(int n=0;n<oseg.size();n++) {
    out << oseg[n].start << " " << oseg[n].stop << endl;
    liveTime+=oseg[n].stop-oseg[n].start;
  }
  out.close();

  cout << endl << "out liveTime = " << liveTime << " (sec) " << liveTime/(24*3600.) << " (days) " << endl << endl;

  cout << "output file list : " << ofile << endl;

  exit(0);
}
