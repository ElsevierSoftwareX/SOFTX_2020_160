{

  #define OFILE  "Segments/SN2016B_s15a3o15_segments.lst"

  #define GPS_START	1135238643
  #define GPS_LENGTH	600
  #define NSEGMENTS	100

  ofstream out;
  out.open(OFILE,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << OFILE << endl;exit(1);}

  for(int i=1;i<=NSEGMENTS;i++) {
    int start = GPS_START+i*GPS_LENGTH;
    int stop  = start+GPS_LENGTH;
    cout << i << "\t" << start << "\t" << stop << "\t" << GPS_LENGTH << endl;
    out << i << "\t" << start << "\t" << stop << "\t" << GPS_LENGTH << endl;
  }

  out.close();

  exit(0);
}
