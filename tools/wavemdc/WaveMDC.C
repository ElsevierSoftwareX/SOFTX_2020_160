{
  #include <vector>

  int seg_start_id = TString(gSystem->Getenv("WMDC_SEG_START")).Atoi();
  int seg_stop_id  = TString(gSystem->Getenv("WMDC_SEG_STOP")).Atoi();

  cout << "seg_start_id : " << seg_start_id << " seg_stop_id : " << seg_stop_id << endl;

  int nsegs = seg_stop_id-seg_start_id+1;

  int segs_start[nsegs]; 
  int segs_length[nsegs]; 

  ifstream in;
  in.open(segmentList,ios::in);
  if (!in.good()) {cout << "Error Opening Segments File : " << segmentList << endl;exit(1);}

  int seg_id;
  int seg_start;
  int seg_stop;
  int seg_length;
  int cnt=0;
  while(true) {
    in >> seg_id >> seg_start >> seg_stop >> seg_length;
    if (!in.good()) break;
    //cout << " " << seg_id << " " << seg_start << " " << seg_stop << " " << seg_length << endl;
    if((seg_id>=seg_start_id)&&(seg_id<=seg_stop_id)) {
      segs_start[cnt]  = seg_start;
      segs_length[cnt] = seg_stop-seg_start;
      cnt++;
    }
  }

  in.close();


  for(int n=0;n<nsegs;n++) {
    cout << n << " " << segs_start[n] << " " << segs_length[n] << endl;
  }

  bool log = true;

  for(int n=0;n<nsegs;n++) {
    cout << n << " " << segs_start[n] << " " << segs_length[n] << endl;
    size_t gps      = segs_start[n];
    size_t length   = segs_length[n];
    if(chName.size()) MDC.WriteFrameFile(frDir, frLabel, gps, length, log, chName);
    else              MDC.WriteFrameFile(frDir, frLabel, gps, length, log);
  }

  exit(0);
}
