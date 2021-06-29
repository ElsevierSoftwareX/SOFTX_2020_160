//
//Test Create FrList
// Author : Gabriele Vedovato

#define LIST_FILE_NAME "SEGMENTS/lists/GHLTV-HBRST14_S6D_R1.frames"
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

  // set frame file list
  int nfiles=tb.frl2FrTree(LIST_FILE_NAME);
  cout << "nfiles : " << nfiles << endl;

  // create merged dq file list
  int dqcat = 1;
  std::vector<waveSegment> olist;
  olist=tb.readSegList(nDQF, DQF, dqcat);
  double ctime=tb.getTimeSegList(olist);
  cout << "ctime : " << int(ctime) << " sec " << ctime/3600. << " h " << ctime/86400. << " day" << endl;

  // create job file list
  std::vector<waveSegment> jlist;
  jlist=tb.getJobList(olist, 300, 600, 8);

  frfile frf;
  int start,stop;
  int jobId=0;
  char ofile[256];
  wavearray<double> w;
  cout << "nJobs : " << jlist.size() << endl;
  //for(int j=0;j<jlist.size();j++) {
  //for(int j=4030;j<jlist.size();j++) {
  for(int j=0;j<2;j++) {
    // read frame file list
    start=jlist[j].start;
    stop=jlist[j].stop;
    jobId=j+1;
    cout << "jobId : " << jobId << " " <<  "start " << start << " stop " << stop << endl;
    frf=tb.getFrList(start, stop, 8);
    // read frames
    tb.readFrames(frf, CHANNEL_NAME, w);
    cout << "Data : " << w.size() << " " << w.rate() << endl;
  }

  olist.clear();

  exit(0);
}
