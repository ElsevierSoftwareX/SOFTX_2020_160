//
// Merge Segments
// this macro read a list of start stop segments from a file
// and merge the overlapping segments
// Example : 
// root -b -l 'MergeSegments.C("ifile.txt","ofile.txt")'
//
// Author : Gabriele Vedovato


void MergeSegments(TString ifile, TString ofile) {

  // --------------------------------------------------------------
  // Open file segment list
  // --------------------------------------------------------------
  ifstream in;
  in.open(ifile.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << ifile << endl;gSystem->Exit(1);}
  cout << "input file list : " << ifile << endl;

  char str[1024];
  int fpos=0;
  int index=0;
  double start;
  double stop;
  waveSegment seg;
  vector<waveSegment> iseg; 
  while (1) {
    fpos=in.tellg();
    in.getline(str,1024);
    if(str[0] == '#') continue;
    in.seekg(fpos, ios::beg);

    in >> start >> stop;
    if (!in.good()) break;
    fpos=in.tellg();
    in.seekg(fpos+1, ios::beg);

    seg.index=index++; seg.start=start; seg.stop=stop; iseg.push_back(seg);
  }

  in.close();

//  cout << "List of input segments" << endl;
//  for(int n=0;n<iseg.size();n++) cout << n << " " << iseg[n].start << " " << iseg[n].stop << endl;

  // --------------------------------------------------------------
  // merge segments
  // --------------------------------------------------------------
  vector<waveSegment> oseg = CWB::Toolbox::unionSegments(iseg);

//  cout << "List of output segments" << endl;
//  for(int n=0;n<oseg.size();n++) cout << n << " " << oseg[n].start << " " << oseg[n].stop << endl;

  // --------------------------------------------------------------
  // Write output file segment list
  // --------------------------------------------------------------
  ofstream out;
  out.open(ofile.Data(),ios::out);
  out.precision(16);
  for(int n=0;n<oseg.size();n++) {
    out << oseg[n].start << " " << oseg[n].stop << endl;
  }
  out.close();

  cout << "output file list : " << ofile << endl;

  exit(0);
}

