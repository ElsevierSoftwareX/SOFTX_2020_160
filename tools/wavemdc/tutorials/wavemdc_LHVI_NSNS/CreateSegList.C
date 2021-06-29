#define GPS_START 931158000
#define GPS_STOP  (931158000+100000)
#define GPS_DELTA 900

#define OFILE_NAME "Segments/NSNS_SegmentsList.txt"

{

  TString ofile_name = OFILE_NAME;

  ofstream out;
  out.open(ofile_name,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << ofile_name << endl;exit(1);}

  int index = 1;
  int start = GPS_START;
  int stop  = start+GPS_DELTA;
  while (1) {
    out << index++ << " " << start << " " << stop << " " << stop-start << endl;
    start+=GPS_DELTA;
    stop+=GPS_DELTA;
    if(start>GPS_STOP) break;
  }
  out.close();
  exit(0);
}
