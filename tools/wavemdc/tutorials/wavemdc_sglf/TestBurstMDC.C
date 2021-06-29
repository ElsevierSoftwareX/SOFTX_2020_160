#define SEGMENTS  "Segments/S6a_LS-segs.txt"

{

  #include <vector>

  #define NWAVEFORM 15

  double F[NWAVEFORM];
  double Q[NWAVEFORM];

  F[0]=70;      Q[0]=8.9;
  F[1]=100;     Q[1]=8.9;
  F[2]=153;     Q[2]=8.9;
  F[3]=235;     Q[3]=8.9;
  F[4]=361;     Q[4]=8.9;
  F[5]=554;     Q[5]=8.9;
  F[6]=849;     Q[6]=8.9;
  F[7]=945;     Q[7]=9;
  F[8]=1053;    Q[8]=9;
  F[9]=1172;    Q[9]=9;
  F[10]=1304;   Q[10]=9;
  F[11]=1451;   Q[11]=9;
  F[12]=1615;   Q[12]=9;
  F[13]=1797;   Q[13]=9;
  F[14]=2000;   Q[14]=9;


  int seg_start_id = TString(gSystem->Getenv("MDC_SEG_START")).Atoi();
  int seg_stop_id  = TString(gSystem->Getenv("MDC_SEG_STOP")).Atoi();

  cout << "seg_start_id : " << seg_start_id << " seg_stop_id : " << seg_stop_id << endl;

  int nsegs = seg_stop_id-seg_start_id+1;

  int segs_start[nsegs]; 
  int segs_length[nsegs]; 

  ifstream in;
  in.open(SEGMENTS,ios::in);
  if (!in.good()) {cout << "Error Opening Segments File : " << SEGMENTS << endl;exit(1);}

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

  int nIFO=3;
  TString ifo[nIFO]={"L1","H1","V1"};
  CWB::mdc MDC(nIFO,ifo); 

  char wf_name[256];
  waveform wf;
  for (int n=0;n<NWAVEFORM;n++) {
    sprintf(wf_name,"SG%dQ%1.1f",F[n],Q[n]);
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    wf.name.ReplaceAll("d0","");
    wf.hp = MDC.GetSGQ(F[n],Q[n]);
    wf.hx = wf.hp;
    wf.hx = 0;
    MDC.AddWaveform(wf);
  }
  MDC.PrintWaveformsList();

  TArrayD parms(3);parms[0]=100000;parms[1]=1;parms[2]=100;
  MDC.SetSkyDistribution(MDC_RANDOM,"",&parms);

  //TArrayD parms(1);parms[0]=100000;  // number of entries
  //MDC.SetSkyDistribution(MDC_MNGD,"",&parms);

  //TArrayD parms(1);parms[0]=100;  // distance max = 100 Mpc
  //MDC.SetSkyDistribution(MDC_GWGC,"../data/GWGCCatalog_Rev1d2.txt",&parms);

  //TArrayD parms(5);parms[0]=1111111;parms[1]=30;parms[2]=60;parms[3]=90;parms[4]=1;
  //MDC.SetSkyDistribution(MDC_EARTH_FIX,"",&parms);

  MDC.SetInjRate(0.01);
  MDC.SetInjTimeJitter(10.0);

  TString frDir   = "frames";
  TString frLabel = "TestBurstMDC";
  bool log        = true;

  for(int n=0;n<nsegs;n++) {
    cout << n << " " << segs_start[n] << " " << segs_length[n] << endl;
    size_t gps      = segs_start[n];
    size_t length   = segs_length[n];
    MDC.WriteFrameFile(frDir, frLabel, gps, length, log);
  }

  exit(0);
}
