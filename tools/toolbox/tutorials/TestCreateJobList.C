//
// Test Create Job List
// Author : Gabriele Vedovato

//#define WRITE_OFILE
#define SORT_FILE_NAME "root/GHLTV-HBRST14_S6D_R1-Sorted.root"
#define JOBS_FILE_NAME "SEGMENTS/S6D_R10_segments/L1H1V1_S6D_R10_jobs.txt"

#define OJOBS_FILE_NAME "ojob.lst"

#define JOB_DAG_LABEL "jobDagFile"

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

  int dqcat = 1;
  std::vector<waveSegment> olist;
  olist=tb.readSegList(nDQF, DQF, dqcat);
  double ctime=tb.getTimeSegList(olist);
  cout << "ctime : " << int(ctime) << " sec " << ctime/3600. << " h " << ctime/86400. << " day" << endl;

  tb.dumpSegList(olist,"olist.txt", false);
  tb.dumpJobList(olist, OJOBS_FILE_NAME, 300, 600, 8);
  cout << "Total SEG Time      : " << tb.getTimeSegList(olist) << endl;

  std::vector<waveSegment> jlist;
  jlist=tb.getJobList(olist, 300, 600, 8);
  cout << "Total JOB SEG Time  : " << tb.getTimeSegList(jlist) << endl;

  cout << "Total SEG Lost Time : " << tb.getTimeSegList(olist)-tb.getTimeSegList(jlist) << endl;

  tb.createDagFile(jlist, JOB_DAG_LABEL);

  exit(0);
}
