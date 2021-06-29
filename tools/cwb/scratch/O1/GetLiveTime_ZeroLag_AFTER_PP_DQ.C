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

#define LIVE_TIME_FILE "/home/vedovato/O1/ER8b_12Sep20Oct_C0101/ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run0a/merge/live_ER8b_12Sep20Oct_C0101_BKG_LF_rMRA_run0a.M1.root"

#define ZERO_LAG_LIVETIME_JOBS	// get live time using only the job segments (600sec)

vector<waveSegment> readSegment(TString ifile);

void GetLiveTime_ZeroLag_AFTER_PP_DQ() {

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
  tree->Draw("start[0]:stop[0]:start[1]:stop[1]","lag[2]==0 && slag[2]==0","goff");
  int nseg = (Int_t)tree->GetSelectedRows();
  cout << "nseg : " << nseg << endl;
  double* start_l1 = tree->GetV1();
  double* stop_l1  = tree->GetV2();
  double* start_h1 = tree->GetV3();
  double* stop_h1  = tree->GetV4();

  vector<waveSegment> L1_jobs;
  int index=0;
  waveSegment seg;
  for(int i=0;i<nseg;i++) {
    seg.index=index++; seg.start=start_l1[i]; seg.stop=stop_l1[i]; L1_jobs.push_back(seg);
  }
  L1_jobs = CWB::Toolbox::sortSegments(L1_jobs);
  vector<waveSegment> H1_jobs;
  index=0;
  for(int i=0;i<nseg;i++) {
    seg.index=index++; seg.start=start_h1[i]; seg.stop=stop_h1[i]; H1_jobs.push_back(seg);
  }
  H1_jobs = CWB::Toolbox::sortSegments(H1_jobs);

  cout << "L1_jobs : " << L1_jobs.size() << endl;
  cout << "H1_jobs : " << H1_jobs.size() << endl;

  double l1_live_time=CWB::Toolbox::getTimeSegList(L1_jobs);
  cout << "l1_live_time " << (int)l1_live_time << endl;
  double h1_live_time=CWB::Toolbox::getTimeSegList(H1_jobs);
  cout << "h1_live_time " << (int)h1_live_time << endl;

  double l1h1_live_time = l1_live_time; 

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

#ifndef ZERO_LAG_LIVETIME_JOBS
  L1_jobs = L1_cat0_and_cat1;
  H1_jobs = H1_cat0_and_cat1;
#endif

  l1_live_time=CWB::Toolbox::getTimeSegList(L1_jobs);
  cout << "l1_live_time " << (int)l1_live_time << " " << l1_live_time/(24.*3600.) << " days" << endl;
  h1_live_time=CWB::Toolbox::getTimeSegList(H1_jobs);
  cout << "h1_live_time " << (int)h1_live_time << " " << l1_live_time/(24.*3600.) << " days" << endl;

  vector<waveSegment> L1H1_jobs = CWB::Toolbox::mergeSegLists(H1_jobs,L1_jobs); 
  l1h1_live_time=CWB::Toolbox::getTimeSegList(L1H1_jobs);
  cout << "----> BEFORE PP DQ : l1h1_live_time : " << (int)l1h1_live_time << " " << l1h1_live_time/(24.*3600.) << " days" << endl;

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
  cout << "----> AFTER PP DQ : l1h1_time_job_cat234 : " << (int)l1h1_time_job_cat234 <<  " " << l1h1_time_job_cat234/(24.*3600.) << " days" << endl;

  cout << "L1H1 CAT234 vetoed time : " << (int)(l1h1_live_time-l1h1_time_job_cat234) << "/" << (int)l1h1_live_time 
       << " Vetoed (%) : " << 100*(l1h1_live_time-l1h1_time_job_cat234)/l1h1_live_time << endl;

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

