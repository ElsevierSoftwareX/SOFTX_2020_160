//
// Show how to create a super lags list
// Author : Gabriele Vedovato


#define NIFO 2
#define NUMBER_OF_SEGMENTS 142
#define NUMBER_OF_SLAG 40  // (0 are the the default lag)
#define SEGMENT_LENGTH 600
#define SLAG_OFFSET 0
#define SLAG_MAX 0
//#define SLAG_LIST_FILE_NAME "slagList_O1N1_Nseg142_Nslag40.txt"
#define SLAG_LIST_FILE_NAME NULL

#define SLAG_DAG_LABEL "slagFileList"

{

  cwbtb tb;
  vector<slag> slist;

  slist=tb.getSlagList(NIFO, NUMBER_OF_SLAG, NUMBER_OF_SEGMENTS, 
                       SEGMENT_LENGTH, SLAG_OFFSET, SLAG_MAX, SLAG_LIST_FILE_NAME);  

  cout << "slag size : " << slist.size() << endl;

  tb.createSlagDagFile(slist, SLAG_DAG_LABEL);

  exit(0);
}
