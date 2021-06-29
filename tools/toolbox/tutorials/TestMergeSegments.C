//
// Test Merge Segments
// Author : Gabriele Vedovato


void TestMergeSegments() {


  vector<waveSegment> ilist1; 
  vector<waveSegment> ilist2; 
  vector<waveSegment> olist; 

  waveSegment seg;

  // fill ilist1
  seg.index=0; seg.start=0; seg.stop =20; ilist1.push_back(seg);
  seg.index=1; seg.start=40; seg.stop =60; ilist1.push_back(seg);
  cout << "List of input segments 1" << endl;
  for(int n=0;n<ilist1.size();n++) cout << n << " " << ilist1[n].start << " " << ilist1[n].stop << endl;

  // fill ilist2
  seg.index=0; seg.start=10; seg.stop =30; ilist2.push_back(seg);
  seg.index=1; seg.start=30; seg.stop =50; ilist2.push_back(seg);
  seg.index=2; seg.start=70; seg.stop =80; ilist2.push_back(seg);
  cout << "List of input segments 2" << endl;
  for(int n=0;n<ilist2.size();n++) cout << n << " " << ilist2[n].start << " " << ilist2[n].stop << endl;

  olist = CWB::Toolbox::mergeSegLists(ilist1, ilist2);

  cout << "List of merged segments" << endl;
  for(int n=0;n<olist.size();n++) cout << n << " " << olist[n].start << " " << olist[n].stop << endl;

  exit(0);
}

