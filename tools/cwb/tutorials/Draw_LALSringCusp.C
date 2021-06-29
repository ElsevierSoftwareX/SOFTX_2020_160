//
// Show comparison between LAL Cosmic String Cusps waveformss vs High frequency cut-off
//
// Author : Gabriele Vedovato
//

#define nMDC	4

#define PLOT_TIME
//#define PLOT_FFT
//#define PLOT_TF	0

//#define ZOOM

#define AMPLITUDE	1.

#define SAVE_PLOT

CWB::mdc* MDC;

void Draw_LALSringCusp() {

  MDC = new CWB::mdc();

  mdcid mdcID[nMDC];

  vector<mdcpar> par;
  par.resize(2);
  par[0].name="frequency"; 
  par[1].name="amplitude"; par[1].value=AMPLITUDE;

  for(int n=0;n<nMDC;n++) {
    if(n==3) par[0].value=50;
    if(n==2) par[0].value=150;
    if(n==1) par[0].value=500;
    if(n==0) par[0].value=1500;
    mdcID[n] = MDC->AddWaveform(MDC_SC_LAL, par);
    cout << n << " " << mdcID[n].name << endl;
  }

  MDC->Print();

  waveform wf[nMDC];
  for(int n=0;n<nMDC;n++) {

    //wf[n] = MDC->GetWaveform(mdcID[n].name,0);
    wf[n] = MDC->GetWaveform(n,0);

    if(wf[n].status==false) {
      cout << "Error : Waveform " << mdcID[n].name << " not exist in the MDC pool !!!" << endl << endl;
      gSystem->Exit(1);
    }

    double dt=1./wf[n].hp.rate();
    double hrss=0.;
    for(int i=0;i<wf[n].hp.size();i++) hrss+=pow(wf[n].hp[i],2);
    hrss=sqrt(hrss*dt);
    wf[n].hx=0.;

    cout << wf[n].name << "\tfreq : " << wf[n].par[0].value << "\tsize : " << wf[n].hp.size() 
         << "\trate : " << wf[n].hp.rate() << "\tstart : " << (int)wf[n].hp.start() << "\thrss : " << hrss << endl;

    wf[n].hp.start(0);	// set start to 0 (needed by draw Method)
    wf[n].hx.start(0);
  }

  int color[nMDC] = {kGreen, kBlue, kRed, kBlack};

  watplot* plot = NULL;

  TLegend* leg;
  double hleg = 0.8-nMDC*0.05;

  TString gfile;

  char title[256];sprintf(title,"Cosmic String Cusps vs High frequency cut-off (amplitude=%g)",AMPLITUDE);

#ifdef PLOT_TIME
  for(int n=0;n<nMDC;n++) {
    if(n==0) {
      plot=MDC->Draw(wf[n].hp,MDC_TIME,"ALP ZOOM",color[n]);
      plot->graph[0]->SetTitle(title);
#ifdef ZOOM
      plot->graph[0]->GetHistogram()->GetXaxis()->SetRangeUser(4.49,4.51);
      plot->graph[0]->GetHistogram()->GetYaxis()->SetRangeUser(1.5,4.5);
      gfile = "CosmicStringCusps_TimeZoom.png";
#else
      plot->graph[0]->GetHistogram()->GetXaxis()->SetRangeUser(4,5);
      plot->graph[0]->GetHistogram()->GetYaxis()->SetRangeUser(-1,5);
      gfile = "CosmicStringCusps_Time.png";
#endif
    } else MDC->Draw(wf[n].hp,MDC_TIME,"same",color[n]);
    plot->graph[n]->SetLineWidth(2);
  }
#endif

#ifdef PLOT_FFT
  for(int n=0;n<nMDC;n++) {
    if(n==0) {
      plot = MDC->Draw(wf[n].hp,MDC_FFT,"ALP ZOOM LOGX",color[n]);	// draw hp in frequency domain
      plot->graph[0]->GetHistogram()->GetXaxis()->SetRangeUser(10,2048);
      plot->graph[0]->GetHistogram()->GetYaxis()->SetRangeUser(1e-5,0.1);
      plot->graph[0]->SetTitle(title);
    } else MDC->Draw(wf[n].hp,MDC_FFT,"same LOGX",color[n]);
    plot->graph[n]->SetLineWidth(2);
    gfile = "CosmicStringCusps_Frequency.png";
  }
#endif

  leg = new TLegend(0.6120401,hleg,0.9615385,0.8721805,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextAlign(22);
  leg->SetTextFont(12);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.04);
  leg->SetLineColor(kBlack);
  leg->SetFillColor(kWhite);
  for(int n=nMDC-1;n>=0;n--) {
    char legLabel[256];
    sprintf(legLabel,"%s",wf[n].name.Data());
    leg->AddEntry(plot->graph[n],legLabel,"lp");
  }
  leg->Draw();

#ifdef SAVE_PLOT
  if(plot) plot->canvas->SaveAs(gfile.Data());
#endif

#ifdef PLOT_TF
  MDC->Draw(wf[PLOT_TF].hp,MDC_TF);	// draw hp in frequency domain
  CWB::STFT* stft = MDC->GetSTFT(); 
  stft->GetHistogram()->GetXaxis()->SetRangeUser(4.4,4.6);
  stft->GetHistogram()->GetYaxis()->SetRangeUser(0,200);
  stft->GetCanvas()->SetLogy(false);
#endif

}
