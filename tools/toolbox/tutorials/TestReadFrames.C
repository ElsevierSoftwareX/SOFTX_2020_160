//
// Test Read Frame Files
// Author : Gabriele Vedovato


//#define WRITE_OFILE
#define SORT_FILE_NAME "root/GHLTV-HBRST14_S6D_R1-Sorted.root"
#define JOBS_FILE_NAME "S6D_R9_segments/L1H1V1_S6D_R9_jobs.txt"

{

  int nDQF=12;
  dqfile DQF[12]={
                  {0 ,"S6D_R9_segments/S6D_OFFLINE_L1SCIENCE.txt"        , 0, 0., false, false},
                  {0 ,"S6D_R9_segments/S6D_OFFLINE_L1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {0 ,"S6D_R9_segments/S6D_OFFLINE_L1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {0 ,"S6D_R9_segments/S6D_OFFLINE_L1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                  {1 ,"S6D_R9_segments/S6D_OFFLINE_H1SCIENCE.txt"        , 0, 0., false, false},
                  {1 ,"S6D_R9_segments/S6D_OFFLINE_H1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {1 ,"S6D_R9_segments/S6D_OFFLINE_H1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {1 ,"S6D_R9_segments/S6D_OFFLINE_H1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                  {2 ,"S6D_R9_segments/S6D_OFFLINE_V1SCIENCE.txt"        , 0, 0., false, false},
                  {2 ,"S6D_R9_segments/S6D_OFFLINE_V1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {2 ,"S6D_R9_segments/S6D_OFFLINE_V1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {2 ,"S6D_R9_segments/S6D_OFFLINE_V1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                 };

  cwbtb tb;

  int cat = 1;
  std::vector<waveSegment> olist;
  double ctime=tb.readSegList(nDQF, DQF, cat, olist);
  cout << "ctime : " << int(ctime) << " sec " << ctime/3600. << " h " << ctime/86400. << " day" << endl;

#ifdef WRITE_OFILE
  tb.dumpSegList(olist,"olist.txt", false);
#endif

  TString sort_file_name = SORT_FILE_NAME;

  char ifile_name[256];
  sprintf(ifile_name,"%s",JOBS_FILE_NAME);
  cout << "Opening File : " << ifile_name << endl;

  ifstream in;
  in.open(ifile_name,ios::in);
  if (!in.good()) {cout << "Error Opening File : " << ifile_name << endl;exit(1);}

  frfile frf;
  int start,stop;
  int jobId=0;
  while (1) {
    in >> start >> stop;
    if (!in.good()) break;
    jobId++;
    if(jobId==3991) {
      cout << "jobId : " << jobId << " " <<  "start " << start << " stop " << stop << endl;
      frf=tb.getFrList(sort_file_name, start, stop, 8);
      tb.dumpFrList(frf, "lists/HBRST14_S6D.lst.3991");
    }
  }
  cout << "OUTPUT" << endl;
  cout << "start  : " << frf.start << endl;
  cout << "stop   : " << frf.stop << endl;
  cout << "length : " << frf.length << endl;
  for(int n=0;n<frf.file.size();n++) {
    cout << "file  : " << frf.file[n].Data() << endl;
  }

  olist.clear();

  exit(0);
}
