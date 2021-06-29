//
// Test Union Segments
// Author : Gabriele Vedovato


void TestUnionSegments() {


  vector<waveSegment> iseg; 
  vector<waveSegment> oseg; 

  waveSegment seg;
  seg.index=0; seg.start=40; seg.stop =60; iseg.push_back(seg);
  seg.index=1; seg.start=30; seg.stop =40; iseg.push_back(seg);
  seg.index=2; seg.start=0; seg.stop =20; iseg.push_back(seg);

  cout << "List of input segments" << endl;
  for(int n=0;n<iseg.size();n++) cout << n << " " << iseg[n].start << " " << iseg[n].stop << endl;

  oseg = CWB::Toolbox::unionSegments(iseg);

  cout << "List of output segments" << endl;
  for(int n=0;n<oseg.size();n++) cout << n << " " << oseg[n].start << " " << oseg[n].stop << endl;

  exit(0);
}

