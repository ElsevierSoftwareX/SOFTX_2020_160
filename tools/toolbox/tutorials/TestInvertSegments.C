//
// Test Invert Segments
// Author : Gabriele Vedovato


void TestInvertSegments() {


  vector<waveSegment> ilist; 
  vector<waveSegment> olist; 

  waveSegment seg;

  // fill ilist1
  seg.index=0; seg.start=0; seg.stop =20; ilist.push_back(seg);
  seg.index=1; seg.start=40; seg.stop =60; ilist.push_back(seg);
  cout << "List of input segments" << endl;
  for(int n=0;n<ilist.size();n++) cout << n << " " << ilist[n].start << " " << ilist[n].stop << endl;

  olist = CWB::Toolbox::invertSegments(ilist);

  cout << "List of inverted segments" << endl;
  for(int n=0;n<olist.size();n++) cout << n << " " << olist[n].start << " " << olist[n].stop << endl;

  exit(0);
}

