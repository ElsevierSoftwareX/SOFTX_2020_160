{

  #define IFILE_NAME "sasi3DL68randgwraynew.dat"
  #define TITLE "Kotake"
  #define WF_FILE_NAME "kotake_STRAIN.txt"
  #define WF_PNG_FILE_NAME "kotake_STRAIN.png"
  #define FFT_PNG_FILE_NAME "kotake_FFT.png"

  // time is in ms

  //#define WRITE_WF
  #define DRAW_FFT

  #define OUTPUT_SAMPLE_RATE 16384.
  //#define OUTPUT_SAMPLE_RATE 1000.

  CWB::Toolbox::checkFile(IFILE_NAME);

  TTree tt("ott","ott");
  //tt.ReadFile(IFILE_NAME,"t/D:hp/D");
  tt.ReadFile(IFILE_NAME,"t/D:hp1/D:hp2/D:hp3/D:hp4/D:hp5/D:hp6/D:hp7/D:hp/D:hp9/D:hp10/D");
//  tt.StartViewer();
  //tv__tree->Draw("hp:t/1000.","","");
//  tv__tree->Draw("hp:t","","");
//return;
  //tv__tree->Draw("hp:t/1000.","","goff");
  tt.Draw("hp:t.","","goff");
  double* hp = tt.GetV1(); 
  double* t  = tt.GetV2(); 
  int size = tt.GetSelectedRows();
  cout << "size : " << size << endl;

  // check if input sample rate is constant
  double mdt=0.;
  for(int i=1;i<size;i++) mdt+=t[i]-t[i-1];
  mdt/=(size-1);
  double isample_rate=1./mdt;
  cout << "Estimate input sample rate : " << isample_rate << endl;
  bool check=true;
  for(int i=1;i<size;i++) {
    if(fabs((t[i]-t[i-1])-mdt)/mdt>0.01) {
      cout << "Warning : " << t[i]-t[i-1] << " " << mdt << " " << 1.-(t[i]-t[i-1])/mdt << endl;
      check=false;
    }
  }
  if(!check) exit(1);

  TCanvas *c1 = new TCanvas("c","C",0,0,800,500);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetGridx();
  c1->SetGridy();
#ifdef DRAW_FFT
  c1->SetLogy();
#endif
  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.65);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleX(0.30);
  gStyle->SetTitleFont(72);
  gStyle->SetTitleFillColor(kWhite);

  TGraph gr(size,t,hp);
  //gr.Draw("ALP");

  double sample_rate = OUTPUT_SAMPLE_RATE;
  int nsize = sample_rate;
  double* nhp = new double[nsize]; 
  double* nt = new double[nsize]; 
  for (int i=0;i<nsize;i++) {
    int j = i;
    nt[j] = i*(1./sample_rate);
    if(nt[j]<t[size-1]) 
      nhp[j] = gr.Eval(nt[j]);
    else 
      nhp[j] = 0.;
  }

  int mkf_argc = 11;  
  const char *mkf_argv[11] = {"mkf", "-Bu","-Hp", "-o", "4", "-a", "200", "-s", "16384","-l",NULL};  
  CWB::Filter filter(mkf_argc, mkf_argv);
//  filter.usage();exit(0);

  // cut frequencies > isample_rate
  wavearray<double> w(nsize); // time series for injections
//  for (int i=0;i<w.size();i++) cout << i << " " << nhp[i] << endl;
  for (int i=0;i<w.size();i++) w.data[i]=filter.Arma(nhp[i]);
  //for (int i=0;i<w.size();i++) w.data[i]=filter.Arma(i==0 ? 1:0..);
//  for (int i=0;i<w.size();i++) cout << i << " " << w.data[i] << endl;
//exit(0);
//  for (int i=0;i<w.size();i++) w.data[i]=nhp[i];
/*
  w.FFT(1);
  double df=(double)sample_rate/(double)nsize;
  for (int i=0;i<nsize;i++) {
    double freq = i*df;
    if(freq>isample_rate) w.data[i]=0;
  }
  w.FFT(-1);
*/
  for (int i=0;i<w.size();i++) nhp[i]=w.data[i];

  int wf_size = sample_rate;  // 1sec

  double dt=1./sample_rate;
  double hrss = 0;
  for (int n=0;n<wf_size;n++) hrss+=dt*nhp[n]*nhp[n];
  hrss=sqrt(hrss);
  cout << "hrss : " << hrss << endl;

  // normalize to hrss=1
  for (int n=0;n<wf_size;n++) nhp[n]/=hrss;
  hrss = 0;
  for (int n=0;n<wf_size;n++) hrss+=dt*nhp[n]*nhp[n];
  hrss=sqrt(hrss);
  cout << "norm hrss : " << hrss << endl;


#ifdef DRAW_FFT
  double* xa = new double[wf_size];
  double* f = new double[wf_size];
  double df=(double)sample_rate/(double)wf_size;

  wavearray<double> x(wf_size); // time series for injections
  for (int i=0;i<x.size()/2;i++) f[i]=(double)i*df;
  for (int i=0;i<x.size();i++) x.data[i]=nhp[i];
  x.FFT(1);
  for (int i=0;i<x.size()/2;i++) xa[i]=sqrt(x.data[2*i]*x.data[2*i]+x.data[2*i+1]*x.data[2*i+1]);
  for (int i=0;i<x.size()/2;i++) xa[i]*=sqrt(1/df);  // one side psd

  TGraph* grf;
  grf = new TGraph(x.size()/2, f, xa);
  grf->SetLineWidth(1);
  grf->SetMarkerColor(kRed);
  grf->SetLineColor(kRed);
  grf->SetMarkerStyle(6);
  //grf->Draw("ALP");
#else
  TGraph* ngr = new TGraph(wf_size,nt,nhp);
  ngr->SetMarkerColor(kRed);
  ngr->SetLineColor(kRed);
  //ngr->Draw("ALP");
#endif


  TMultiGraph *multigraph = new TMultiGraph();
  multigraph->SetName("");
  multigraph->SetTitle(TITLE);

#ifdef DRAW_FFT
    multigraph->Add(grf);
#else
    multigraph->Add(ngr);
#endif

  multigraph->Paint("AP");
  multigraph->GetHistogram()->SetMarkerColor(kBlack);
  multigraph->GetHistogram()->SetMarkerStyle(6);
  multigraph->GetHistogram()->SetMarkerSize(1.0);
  multigraph->GetHistogram()->GetYaxis()->SetTitle("");
#ifdef DRAW_FFT
    multigraph->GetHistogram()->GetXaxis()->SetTitle("Frequency [Hz]");
#else
    multigraph->GetHistogram()->GetXaxis()->SetTitle("Time [sec]");
#endif
  multigraph->GetHistogram()->GetXaxis()->CenterTitle(true);
  multigraph->GetHistogram()->GetXaxis()->SetLabelFont(42);
  multigraph->GetHistogram()->GetXaxis()->SetTitleOffset(1.1);
  multigraph->GetHistogram()->GetXaxis()->SetTitleFont(42);
  multigraph->GetHistogram()->GetYaxis()->CenterTitle(true);
  multigraph->GetHistogram()->GetYaxis()->SetLabelFont(42);
  multigraph->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
  multigraph->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
  multigraph->GetHistogram()->GetYaxis()->SetTitleFont(42);
  multigraph->Draw("alp");

#ifdef DRAW_FFT
  c1->Print(FFT_PNG_FILE_NAME);
#else
  c1->Print(WF_PNG_FILE_NAME);
#endif

#ifdef WRITE_WF
  ofstream out;
  out.open(WF_FILE_NAME,ios::out);
  for (int n=0;n<wf_size;n++) out << nhp[n] << endl;
  out.close();
#endif

}
