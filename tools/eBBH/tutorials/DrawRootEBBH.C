#define DRAW_TIME
//#define DRAW_FFT
//#define SAVE_PLOT

void DrawRootEBBH(TString fName, int entry) {  

  TFile* efile = new TFile(fName);
  if(efile==NULL) {
    cout << "CWB::mdc::AddWaveform - Error opening root file : " << fName.Data() << endl;
    exit(1);
  }

  int id;
  double m1,m2,rp0,e0;
  wavearray<double>* hp = new wavearray<double>;
  wavearray<double>* hx = new wavearray<double>;

  TTree* etree = (TTree *) efile->Get("ebbh");
  if(etree==NULL) {
    cout << "CWB::mdc::AddWaveform - file : " << fName.Data()
         << " not contains tree ebbh" << endl;
    exit(1);
  }
  etree->SetBranchAddress("id",&id);
  etree->SetBranchAddress("m1",&m1);
  etree->SetBranchAddress("m2",&m2);
  etree->SetBranchAddress("rp0",&rp0);
  etree->SetBranchAddress("e0",&e0);
  etree->SetBranchAddress("hp",&hp);
  etree->SetBranchAddress("hx",&hx);

  int esize = etree->GetEntries();

  if(entry>=esize) {
    cout << "entry not present in the tree" << endl;
    return;
  }

  etree->GetEntry(entry);
  cout << id << " " << m1 << " " << m2 << " " << rp0 << " " << e0 << endl;

  gwavearray<double>* gw = new gwavearray<double>(hp);

#ifdef DRAW_TIME
  gw->Draw();
  gw->Draw(hx,GWAT_TIME,"SAME",kRed);
#endif

#ifdef DRAW_FFT
  gw->Draw(GWAT_FFT);
  gw->Draw(hx,GWAT_FFT,"SAME",kRed);
#endif

#ifdef SAVE_PLOT
  watplot* plot = gw->GetWATPLOT();

  char gtitle[256];
  sprintf(gtitle,"eBBH : hp(black) - hx(red)");
#ifdef DRAW_TIME
  plot->gtitle(gtitle,"time(sec)","amplitude");
  TString gfile="eBHH_time_plot.png";
#endif
#ifdef DRAW_FFT
  plot->gtitle(gtitle,"freq(hz)","amplitude");
  TString gfile="eBHH_freq_plot.png";
#endif

  // save plot to file
  (*plot) >> gfile;

  exit(0);
#endif

  delete hp;
  delete hx;
  delete efile;

}
