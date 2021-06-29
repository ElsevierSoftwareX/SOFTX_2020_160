//
// Test Merge Segment Lists
// Author : Gabriele Vedovato


#define WRITE_OFILE

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

  int dqcat = 1;
  std::vector<waveSegment> olist;
  double ctime=tb.readSegList(nDQF, DQF, dqcat, olist);
  cout << "ctime : " << int(ctime) << " sec " << ctime/3600. << " h " << ctime/86400. << " day" << endl;

#ifdef WRITE_OFILE
  tb.dumpSegList(olist,"olist.txt", false);
#endif

  //olist.clear();

  exit(0);
}
