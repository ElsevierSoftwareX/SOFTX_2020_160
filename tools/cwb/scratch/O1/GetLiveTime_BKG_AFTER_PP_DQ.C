#define L1_FILE_CAT0 "/home/vedovato/O1//DQvetos/ER8b_12Sep20Oct_C0101/L1Cat0.txt"
#define H1_FILE_CAT0 "/home/vedovato/O1//DQvetos/ER8b_12Sep20Oct_C0101/H1Cat0.txt"

#define L1_FILE_CAT1 "/home/vedovato/O1//DQvetos/ER8b_12Sep20Oct_C0101/L1Cat1.txt"
#define H1_FILE_CAT1 "/home/vedovato/O1//DQvetos/ER8b_12Sep20Oct_C0101/H1Cat1.txt"

#define L1_FILE_CAT2 "/home/vedovato/O1/DQvetos/ER8b_12Sep20Oct_C0101/L1Cat2.txt"
#define H1_FILE_CAT2 "/home/vedovato/O1/DQvetos/ER8b_12Sep20Oct_C0101/H1Cat2.txt"

//#define L1_FILE_CAT2 "L1Cat2.txt"
//#define H1_FILE_CAT2 "H1Cat2.txt"

#define L1_FILE_CAT3 "/home/vedovato/O1//DQvetos/ER8b_12Sep20Oct_C0101/HVETO_L1_SEP12OCT20_MERGED.txt"
#define H1_FILE_CAT3 "/home/vedovato/O1//DQvetos/ER8b_12Sep20Oct_C0101/HVETO_H1_SEP12OCT20_MERGED.txt"

//#define L1_FILE_CAT3 "HVETO_L1_SEP12OCT20_MERGED.txt"
//#define H1_FILE_CAT3 "HVETO_H1_SEP12OCT20_MERGED.txt"

#define L1_FILE_CAT4 "/home/vedovato/O1//DQvetos/ER8b_12Sep20Oct_C0101/L1Cat4.txt"
#define H1_FILE_CAT4 "/home/vedovato/O1//DQvetos/ER8b_12Sep20Oct_C0101/H1Cat4.txt"

//#define L1_FILE_CAT4 "L1Cat4.txt"
//#define H1_FILE_CAT4 "H1Cat4.txt"

// livetime background run0a
//#define LIVE_TIME_FILE "/home/vedovato/O1/ER8b_12Sep20Oct_C0101/ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run0a/merge/live_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run0a.M1.root"
// livetime background runa1to10a
#define LIVE_TIME_FILE "/home/vedovato/O1/ER8b_12Sep20Oct_C0101/ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a/merge/live_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run1ato10a.M1.root"

// ------------------------------------------------------------------------------------------------------
// This macro compute a conservative percentage of lost live time after the CAT2/4 Hveto Post-Production DQ
// We compute only the percentage for each single detector : L1,H1 
// We do not consider the combined percentage obtaine from the intersection of L1,H1 but only the sum of L1,H1
// 1) extract from the livetime tree all segment-jobs(600s) used in the analysis (only 1 lag)
// 2) extract all segment-dq after DQ CAT0/1/2/4 Hveto  
// 3) for each segment-job and each detector we compute the intersection of segment-job & segment-dq = segment-job-dq
// 4) for each detector the livetime after the CAT2/4 Hveto Post-Production DQ is the sum of all segment-job-dq livetimes
// 5) the conservative livetime is sum of L1+H1 after the step 4) 
//
// NOTE : works only foe O1 : L1,H1
// NOTE : this procedure is correct only when the following user_parameters parameters are used :
// segLen = 600;
// segMLS = 600;
// All segments have SegLen=600 sec !!! 
// ------------------------------------------------------------------------------------------------------

#define SEGLEN  600

vector<waveSegment> readSegment(TString ifile);

void GetLiveTime_BKG_AFTER_PP_DQ() {

  TFile *live = TFile::Open(LIVE_TIME_FILE);
  if(live==NULL) {
    cout << "CreateSymbolicLinksC01 - Error : File " << LIVE_TIME_FILE << " not exist !!!" << endl;
    gSystem->Exit(1);
  }

  // get list of jobId which belong to C01
  TTree* tree = (TTree *) gROOT->FindObject("liveTime");
  if(tree==NULL) {
    cout << "CheckCAT2 - Error : liveTime tree not found !!!" << endl;
    gSystem->Exit(1);
  }

  long int isize = tree->GetEntries();
  isize=2000000000;
  cout << "isize : " << isize << endl;
  tree->SetEstimate(isize);
  //tree->Draw("start[0]:start[1]","","goff",isize);
  tree->Draw("start[0]:start[1]","lag[2]==0","goff",isize);	// select only one lag for each slag
  long int nseg = (Int_t)tree->GetSelectedRows();
  cout << "nseg : " << nseg << endl;
  double* start_l1 = tree->GetV1();
  double* start_h1 = tree->GetV2();

  double mintime_h1=1e20;
  for(long int i=0;i<nseg;i++) if(start_h1[i]<mintime_h1) mintime_h1=start_h1[i];
  double maxtime_h1=0;
  for(long int i=0;i<nseg;i++) if(start_h1[i]>maxtime_h1) maxtime_h1=start_h1[i];
  int size=(maxtime_h1-mintime_h1)/SEGLEN;
  cout << "h1 size : " << size << endl;
  wavearray<int> seg_h1(size+1);
  seg_h1=0;
  for(long int i=0;i<nseg;i++) {
    double start=start_h1[i]-mintime_h1;
    int jobid=start/SEGLEN;
    seg_h1[jobid]++;
  }


  double mintime_l1=1e20;
  for(long int i=0;i<nseg;i++) if(start_l1[i]<mintime_l1) mintime_l1=start_l1[i];
  double maxtime_l1=0;
  for(long int i=0;i<nseg;i++) if(start_l1[i]>maxtime_l1) maxtime_l1=start_l1[i];
  size=(maxtime_l1-mintime_l1)/SEGLEN;
  cout << "l1 size : " << size << endl;
  wavearray<int> seg_l1(size+1);
  seg_l1=0;
  for(long int i=0;i<nseg;i++) {
    double start=start_l1[i]-mintime_l1;
    int jobid=start/SEGLEN;
    seg_l1[jobid]++;
  }

  // ---------------------------------------------------------------------
  // CAT0
  // ---------------------------------------------------------------------
  vector<waveSegment> segCat0_L1 = CWB::Toolbox::readSegments(L1_FILE_CAT0);
  vector<waveSegment> segCat0_H1 = CWB::Toolbox::readSegments(H1_FILE_CAT0);

  double h1_time_cat0=CWB::Toolbox::getTimeSegList(segCat0_H1);
  cout << "h1_time_cat0 : " << (int)h1_time_cat0 << endl;
  double l1_time_cat0=CWB::Toolbox::getTimeSegList(segCat0_L1);
  cout << "l1_time_cat0 : " << (int)l1_time_cat0 << endl;

  vector<waveSegment> L1H1_cat0 = CWB::Toolbox::mergeSegLists(segCat0_L1,segCat0_H1);
  double l1h1_time_cat0=CWB::Toolbox::getTimeSegList(L1H1_cat0);
  cout << "l1h1_time_cat0 : " << (int)l1h1_time_cat0 << " " << l1h1_time_cat0/(24.*3600.) << endl;

  // ---------------------------------------------------------------------
  // CAT1
  // ---------------------------------------------------------------------
  vector<waveSegment> segCat1_L1 = CWB::Toolbox::readSegments(L1_FILE_CAT1);
  vector<waveSegment> segCat1_H1 = CWB::Toolbox::readSegments(H1_FILE_CAT1);

  double h1_time_cat1=CWB::Toolbox::getTimeSegList(segCat1_H1);
  cout << "h1_time_cat1 : " << (int)h1_time_cat1 << endl;
  double l1_time_cat1=CWB::Toolbox::getTimeSegList(segCat1_L1);
  cout << "l1_time_cat1 : " << (int)l1_time_cat1 << endl;

  vector<waveSegment> isegCat1_L1 = CWB::Toolbox::invertSegments(segCat1_L1);
  vector<waveSegment> isegCat1_H1 = CWB::Toolbox::invertSegments(segCat1_H1);

  vector<waveSegment> H1_cat0_and_cat1 = CWB::Toolbox::mergeSegLists(segCat0_H1,isegCat1_H1);
  double h1_time_cat0_cat1=CWB::Toolbox::getTimeSegList(H1_cat0_and_cat1);
  cout << "h1_time_cat0_cat1 : " << (int)h1_time_cat0_cat1 << endl;

  vector<waveSegment> L1_cat0_and_cat1 = CWB::Toolbox::mergeSegLists(segCat0_L1,isegCat1_L1);
  double l1_time_cat0_cat1=CWB::Toolbox::getTimeSegList(L1_cat0_and_cat1);
  cout << "l1_time_cat0_cat1 : " << (int)l1_time_cat0_cat1 << endl;

  vector<waveSegment> L1H1_cat0_and_cat1 = CWB::Toolbox::mergeSegLists(H1_cat0_and_cat1,L1_cat0_and_cat1);
  double l1h1_time_cat0_cat1=CWB::Toolbox::getTimeSegList(L1H1_cat0_and_cat1);
  cout << "l1h1_time_cat0_cat1 : " << (int)l1h1_time_cat0_cat1 << " " << l1h1_time_cat0_cat1/(24.*3600.) << endl;

  // ---------------------------------------------------------------------
  // JOBS -> CAT0+CAT1
  // ---------------------------------------------------------------------

  vector<waveSegment> L1_jobs = L1_cat0_and_cat1;
  vector<waveSegment> H1_jobs = H1_cat0_and_cat1;

  l1_live_time=CWB::Toolbox::getTimeSegList(L1_jobs);
  cout << "l1_live_time " << (int)l1_live_time << " " << l1_live_time/(24.*3600.) << " days" << endl;
  h1_live_time=CWB::Toolbox::getTimeSegList(H1_jobs);
  cout << "h1_live_time " << (int)h1_live_time << " " << l1_live_time/(24.*3600.) << " days" << endl;

  vector<waveSegment> L1H1_jobs = CWB::Toolbox::mergeSegLists(H1_jobs,L1_jobs); 
  l1h1_live_time=CWB::Toolbox::getTimeSegList(L1H1_jobs);
  cout << "l1h1_live_time : " << (int)l1h1_live_time << " " << l1h1_live_time/(24.*3600.) << " days" << endl;

  // ---------------------------------------------------------------------
  // CAT2
  // ---------------------------------------------------------------------
  vector<waveSegment> segCat2_L1 = CWB::Toolbox::readSegments(L1_FILE_CAT2);
  vector<waveSegment> segCat2_H1 = CWB::Toolbox::readSegments(H1_FILE_CAT2);

  double h1_time_cat2=CWB::Toolbox::getTimeSegList(segCat2_H1);
  cout << "h1_time_cat2 : " << (int)h1_time_cat2 << endl;
  double l1_time_cat2=CWB::Toolbox::getTimeSegList(segCat2_L1);
  cout << "l1_time_cat2 : " << (int)l1_time_cat2 << endl;

  vector<waveSegment> isegCat2_L1 = CWB::Toolbox::invertSegments(segCat2_L1);
  vector<waveSegment> isegCat2_H1 = CWB::Toolbox::invertSegments(segCat2_H1);

  vector<waveSegment> H1_jobs_and_cat2 = CWB::Toolbox::mergeSegLists(H1_jobs,isegCat2_H1); 
  double h1_time_job_cat2=CWB::Toolbox::getTimeSegList(H1_jobs_and_cat2);
  cout << "h1_time_job_cat2 : " << (int)h1_time_job_cat2 
       << " Vetoed (%) : " << 100*(h1_live_time-h1_time_job_cat2)/h1_live_time << endl;

  vector<waveSegment> L1_jobs_and_cat2 = CWB::Toolbox::mergeSegLists(L1_jobs,isegCat2_L1); 
  double l1_time_job_cat2=CWB::Toolbox::getTimeSegList(L1_jobs_and_cat2);
  cout << "l1_time_job_cat2 : " << (int)l1_time_job_cat2 
       << " Vetoed (%) : " << 100*(l1_live_time-l1_time_job_cat2)/l1_live_time << endl;

  vector<waveSegment> L1H1_jobs_and_cat2 = CWB::Toolbox::mergeSegLists(H1_jobs_and_cat2,L1_jobs_and_cat2); 
  double l1h1_time_job_cat2=CWB::Toolbox::getTimeSegList(L1H1_jobs_and_cat2);
  cout << "l1h1_time_job_cat2 : " << (int)l1h1_time_job_cat2 <<  " " << l1h1_time_job_cat2/(24.*3600.) << " days" << endl;

  cout << "CAT2 vetoed time : " << (int)(l1h1_live_time-l1h1_time_job_cat2) << "/" << (int)l1h1_live_time 
       << " Vetoed (%) : " << 100*(l1h1_live_time-l1h1_time_job_cat2)/l1h1_live_time << endl;

  // ---------------------------------------------------------------------
  // CAT3
  // ---------------------------------------------------------------------
  vector<waveSegment> segCat3_L1 = CWB::Toolbox::readSegments(L1_FILE_CAT3);
  vector<waveSegment> segCat3_H1 = CWB::Toolbox::readSegments(H1_FILE_CAT3);

  double h1_time_cat3=CWB::Toolbox::getTimeSegList(segCat3_H1);
  cout << "h1_time_cat3 : " << (int)h1_time_cat3 << endl;
  double l1_time_cat3=CWB::Toolbox::getTimeSegList(segCat3_L1);
  cout << "l1_time_cat3 : " << (int)l1_time_cat3 << endl;

  vector<waveSegment> isegCat3_L1 = CWB::Toolbox::invertSegments(segCat3_L1);
  vector<waveSegment> isegCat3_H1 = CWB::Toolbox::invertSegments(segCat3_H1);

  vector<waveSegment> H1_jobs_and_cat3 = CWB::Toolbox::mergeSegLists(H1_jobs,isegCat3_H1); 
  double h1_time_job_cat3=CWB::Toolbox::getTimeSegList(H1_jobs_and_cat3);
  cout << "h1_time_job_cat3 : " << (int)h1_time_job_cat3 
       << " Vetoed (%) : " << 100*(h1_live_time-h1_time_job_cat3)/h1_live_time << endl;

  vector<waveSegment> L1_jobs_and_cat3 = CWB::Toolbox::mergeSegLists(L1_jobs,isegCat3_L1); 
  double l1_time_job_cat3=CWB::Toolbox::getTimeSegList(L1_jobs_and_cat3);
  cout << "l1_time_job_cat3 : " << (int)l1_time_job_cat3 
       << " Vetoed (%) : " << 100*(l1_live_time-l1_time_job_cat3)/l1_live_time << endl;

  vector<waveSegment> L1H1_jobs_and_cat3 = CWB::Toolbox::mergeSegLists(H1_jobs_and_cat3,L1_jobs_and_cat3); 
  double l1h1_time_job_cat3=CWB::Toolbox::getTimeSegList(L1H1_jobs_and_cat3);
  cout << "l1h1_time_job_cat3 : " << (int)l1h1_time_job_cat3 <<  " " << l1h1_time_job_cat3/(24.*3600.) << " days" << endl;

  cout << "CAT3 vetoed time : " << (int)(l1h1_live_time-l1h1_time_job_cat3) << "/" << (int)l1h1_live_time 
       << " Vetoed (%) : " << 100*(l1h1_live_time-l1h1_time_job_cat3)/l1h1_live_time << endl;

  // ---------------------------------------------------------------------
  // CAT4
  // ---------------------------------------------------------------------
  vector<waveSegment> segCat4_L1 = CWB::Toolbox::readSegments(L1_FILE_CAT4);
  vector<waveSegment> segCat4_H1 = CWB::Toolbox::readSegments(H1_FILE_CAT4);

  double h1_time_cat4=CWB::Toolbox::getTimeSegList(segCat4_H1);
  cout << "h1_time_cat4 : " << (int)h1_time_cat4 << endl;
  double l1_time_cat4=CWB::Toolbox::getTimeSegList(segCat4_L1);
  cout << "l1_time_cat4 : " << (int)l1_time_cat4 << endl;

  vector<waveSegment> isegCat4_L1 = CWB::Toolbox::invertSegments(segCat4_L1);
  vector<waveSegment> isegCat4_H1 = CWB::Toolbox::invertSegments(segCat4_H1);

  vector<waveSegment> H1_jobs_and_cat4 = CWB::Toolbox::mergeSegLists(H1_jobs,isegCat4_H1); 
  double h1_time_job_cat4=CWB::Toolbox::getTimeSegList(H1_jobs_and_cat4);
  cout << "h1_time_job_cat4 : " << (int)h1_time_job_cat4
       << " Vetoed (%) : " << 100*(h1_live_time-h1_time_job_cat4)/h1_live_time << endl;

  vector<waveSegment> L1_jobs_and_cat4 = CWB::Toolbox::mergeSegLists(L1_jobs,isegCat4_L1); 
  double l1_time_job_cat4=CWB::Toolbox::getTimeSegList(L1_jobs_and_cat4);
  cout << "l1_time_job_cat4 : " << (int)l1_time_job_cat4 
       << " Vetoed (%) : " << 100*(l1_live_time-l1_time_job_cat4)/l1_live_time << endl;

  vector<waveSegment> L1H1_jobs_and_cat4 = CWB::Toolbox::mergeSegLists(H1_jobs_and_cat4,L1_jobs_and_cat4); 
  double l1h1_time_job_cat4=CWB::Toolbox::getTimeSegList(L1H1_jobs_and_cat4);
  cout << "l1h1_time_job_cat4 : " << (int)l1h1_time_job_cat4 <<  " " << l1h1_time_job_cat4/(24.*3600.) << " days" << endl;

  cout << "CAT4 vetoed time : " << (int)(l1h1_live_time-l1h1_time_job_cat4) << "/" << (int)l1h1_live_time 
       << " Vetoed (%) : " << 100*(l1h1_live_time-l1h1_time_job_cat4)/l1h1_live_time << endl;

  // ---------------------------------------------------------------------
  // L1 CAT2+CAT3+CAT4
  // ---------------------------------------------------------------------
  vector<waveSegment> L1_jobs_and_cat23  = CWB::Toolbox::mergeSegLists(L1_jobs_and_cat2,L1_jobs_and_cat3); 
  vector<waveSegment> L1_jobs_and_cat234 = CWB::Toolbox::mergeSegLists(L1_jobs_and_cat23,L1_jobs_and_cat4); 
  double l1_time_job_cat234=CWB::Toolbox::getTimeSegList(L1_jobs_and_cat234);
  cout << "l1_time_job_cat234 : " << (int)l1_time_job_cat234 <<  " " << l1_time_job_cat234/(24.*3600.) << " days" << endl;

  cout << "L1 CAT234 vetoed time : " << (int)(l1_live_time-l1_time_job_cat234) << "/" << (int)l1_live_time 
       << " Vetoed (%) : " << 100*(l1_live_time-l1_time_job_cat234)/l1_live_time << endl;


  // ---------------------------------------------------------------------
  // H1 CAT2+CAT3+CAT4
  // ---------------------------------------------------------------------
  vector<waveSegment> H1_jobs_and_cat23  = CWB::Toolbox::mergeSegLists(H1_jobs_and_cat2,H1_jobs_and_cat3); 
  vector<waveSegment> H1_jobs_and_cat234 = CWB::Toolbox::mergeSegLists(H1_jobs_and_cat23,H1_jobs_and_cat4); 
  double h1_time_job_cat234=CWB::Toolbox::getTimeSegList(H1_jobs_and_cat234);
  cout << "h1_time_job_cat234 : " << (int)h1_time_job_cat234 <<  " " << h1_time_job_cat234/(24.*3600.) << " days" << endl;

  cout << "H1 CAT234 vetoed time : " << (int)(h1_live_time-h1_time_job_cat234) << "/" << (int)h1_live_time 
       << " Vetoed (%) : " << 100*(h1_live_time-h1_time_job_cat234)/h1_live_time << endl;


  // ---------------------------------------------------------------------
  // L1H1 CAT2+CAT3+CAT4
  // ---------------------------------------------------------------------
  vector<waveSegment> L1H1_jobs_and_cat23  = CWB::Toolbox::mergeSegLists(L1H1_jobs_and_cat2,L1H1_jobs_and_cat3); 
  vector<waveSegment> L1H1_jobs_and_cat234 = CWB::Toolbox::mergeSegLists(L1H1_jobs_and_cat23,L1H1_jobs_and_cat4); 
  double l1h1_time_job_cat234=CWB::Toolbox::getTimeSegList(L1H1_jobs_and_cat234);
  cout << "l1h1_time_job_cat234 : " << (int)l1h1_time_job_cat234 <<  " " << l1h1_time_job_cat234/(24.*3600.) << " days" << endl;

  cout << "L1H1 CAT234 vetoed time : " << (int)(l1h1_live_time-l1h1_time_job_cat234) << "/" << (int)l1h1_live_time 
       << " Vetoed (%) : " << 100*(l1h1_live_time-l1h1_time_job_cat234)/l1h1_live_time << endl;

  // ---------------------------------------------------------------------
  // L1 CAT2+CAT3+CAT4+SEG
  // ---------------------------------------------------------------------

  vector<waveSegment> seg(1);
  seg.index=1; 

  double l1_tot_time=0;
  double l1_tot_time_after_dq=0;
  for(long int i=0;i<seg_l1.size();i++) {
    if(seg_l1[i]==0) continue;
    l1_tot_time+=seg_l1[i]*SEGLEN; 
    seg[0].start=mintime_l1+i*SEGLEN; seg[0].stop=seg[0].start+SEGLEN;
    vector<waveSegment> segs = CWB::Toolbox::mergeSegLists(L1_jobs_and_cat234,seg); 
    double segs_time=CWB::Toolbox::getTimeSegList(segs);
    //cout << i << " segs_time : " << segs_time << endl;
    l1_tot_time_after_dq+=seg_l1[i]*segs_time; 
  }
  cout << "L1 CAT234 vetoed time (ALL) : " << (l1_tot_time-l1_tot_time_after_dq) << "/" << l1_tot_time 
       << " Vetoed (%) : " << 100*(l1_tot_time-l1_tot_time_after_dq)/l1_tot_time << endl;


  // ---------------------------------------------------------------------
  // H1 CAT2+CAT3+CAT4+SEG
  // ---------------------------------------------------------------------

  double h1_tot_time=0;
  double h1_tot_time_after_dq=0;
  for(long int i=0;i<seg_h1.size();i++) {
    if(seg_h1[i]==0) continue;
    h1_tot_time+=seg_h1[i]*SEGLEN; 
    seg[0].start=mintime_h1+i*SEGLEN; seg[0].stop=seg[0].start+SEGLEN;
    vector<waveSegment> segs = CWB::Toolbox::mergeSegLists(H1_jobs_and_cat234,seg); 
    double segs_time=CWB::Toolbox::getTimeSegList(segs);
    //cout << i << " segs_time : " << segs_time << endl;
    h1_tot_time_after_dq+=seg_h1[i]*segs_time; 
  }
  cout << "H1 CAT234 vetoed time (ALL) : " << (h1_tot_time-h1_tot_time_after_dq) << "/" << h1_tot_time 
       << " Vetoed (%) : " << 100*(h1_tot_time-h1_tot_time_after_dq)/h1_tot_time << endl;

  exit(0);
}

vector<waveSegment> readSegment(TString ifile) {

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

  return iseg;
}

