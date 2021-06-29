//
// Read/Display Objects from strain job root file
// Author : Gabriele Vedovato

{

  #define IFILE  "strain_931158216_184_ADV_SIM_NSNS_L1H1V1_2G_run3_job1.root"

  #define DRAW_IFO
  #define DRAW_FFT
  //#define SAVE_PLOT

  //#define DATA_TYPE "mdc/L1"		// read mdc L1 data
  #define DATA_TYPE "strain/L1"		// read strain L1 data

  //#define PRINT_CFG
  //#define PRINT_HISTORY
  //#define PRINT_NETWORK
 
  TFile* ifile = new TFile(IFILE);
  if(ifile==NULL) {cout << "Error : file " << IFILE << " not found" <<  endl;exit(1);}
  ifile->ls();

  // read config object
  CWB::config* cfg = (CWB::config*)ifile->Get("config");
#ifdef PRINT_CFG
  if(cfg) cfg->Print();
#endif

  // read history object
  CWB::History* history = (CWB::History*)ifile->Get("history");
#ifdef PRINT_HISTORY
  if(history) history->Print();
#endif

  // read network object
  network* net = (network*)ifile->Get("network");
#ifdef PRINT_NETWORK
  if(net) net->print();
#endif

#ifdef DRAW_IFO
  wavearray<double>* x = (wavearray<double>*)ifile->Get(DATA_TYPE);
  if(x==NULL) {cout << "Error : wavearray not found" <<  endl;exit(1);}
  cout.precision(14);
  cout << "start : " << x->start() << " (sec)" << endl;
  cout << "len   : " << x->size()/x->rate() << " (sec)" << endl;
  x->start(0);

  gwavearray<double> gx(x);

#ifdef DRAW_FFT
  gx.Draw(GWAT_FFT);
#else
  gx.Draw(GWAT_TIME);
#endif
  watplot* plot = gx.GetWATPLOT();

#ifdef SAVE_PLOT
  plot->gtitle(IFILE,"time(sec)","amplitude");

  // save plot to file
  TString gfile="gwavearray_plot.png";
  (*plot) >> gfile;

  exit(0);
#endif
#endif

}
