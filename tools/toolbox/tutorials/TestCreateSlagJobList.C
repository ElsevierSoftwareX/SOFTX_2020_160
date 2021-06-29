//
// Test Create SuperLags Job List
// Author : Gabriele Vedovato

#define SORT_FILE_NAME "root/GHLTV-HBRST14_S6D_R1-Sorted.root"
#define CHANNEL_NAME "L1:GW-H"

{

  int nDQF=12;
  dqfile DQF[12]={
                  {0 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_L1SCIENCE.txt"        , 0, 0., false, false},
                  {0 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_L1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {0 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_L1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {0 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_L1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                  {1 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_H1SCIENCE.txt"        , 0, 0., false, false},
                  {1 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_H1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {1 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_H1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {1 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_H1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                  {2 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_V1SCIENCE.txt"        , 0, 0., false, false},
                  {2 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_V1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {2 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_V1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {2 ,"SEGMENTS/S6D_R10_segments/S6D_OFFLINE_V1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                 };

  cwbtb tb;

  // create merged dq file list
  int dqcat = 1;
  std::vector<waveSegment> olist;
  olist=tb.readSegList(nDQF, DQF, dqcat);
  double ctime=tb.getTimeSegList(olist);
  cout << "ctime : " << int(ctime) << " sec " << ctime/3600. << " h " << ctime/86400. << " day" << endl;
  cout << "Total SEG Time          : " << tb.getTimeSegList(olist) << endl;

  // create slag job list
  std::vector<waveSegment> jlist;
  jlist=tb.getSlagJobList(olist, 600);

  cout.precision(14);
  waveSegment SEG;
  waveSegment MAXSEG;
  std::vector<waveSegment> slist;
  std::vector<waveSegment> mlist;
  for(int i=0;i<jlist.size();i++) {
    //cout << "jlist : " << i << " " << jlist[i].start << " " << jlist[i].stop << " " << jlist[i].stop-jlist[i].start << endl;
    mlist.clear();
    SEG.start=jlist[i].start-8;
    SEG.stop=jlist[i].stop+8;
    mlist.push_back(SEG);
    mlist=tb.mergeSegLists(olist, mlist);
    MAXSEG=tb.getMaxSeg(mlist);
    for(int j=0;j<mlist.size();j++) {
    //  cout << "mlist : " << j << " " << mlist[j].start << " " << mlist[j].stop << " " << mlist[j].stop-mlist[j].start << endl;
    }
    MAXSEG.start+=8;
    MAXSEG.stop-=8;
    //cout << "maxseg: " << i << "/" << jlist.size() << " " << MAXSEG.start << " " << MAXSEG.stop << " " << MAXSEG.stop-MAXSEG.start << endl;
    //cout << endl;
    if((MAXSEG.stop-MAXSEG.start)>=300) slist.push_back(MAXSEG); 
  }
  cout << "Total SLAG JOB SEG Time : " << tb.getTimeSegList(slist) << endl;

  cout << "Total LOST JOB SEG Time : " << tb.getTimeSegList(olist)-tb.getTimeSegList(slist) << endl;

  olist.clear();
  jlist.clear();

  exit(0);
}
