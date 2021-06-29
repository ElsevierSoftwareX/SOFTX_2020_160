//
// Show comparison between SGE/WNB CWB vs LAL waveforms
//
// MDC_SGE_LAL use CG instead of SG !!!
// MDC_WNB_LAL has a gaussian envelope also in frequency
//
// Author : Gabriele Vedovato
//

//#define TYPE_MDC	MDC_CGE		// CWB
//#define TYPE_MDC	MDC_SGE		// CWB
#define TYPE_MDC	MDC_SGE_LAL	// LAL

//#define TYPE_MDC	MDC_WNB		// CWB
//#define TYPE_MDC	MDC_WNB_LAL	// LAL

#define PLOT_TIME
//#define PLOT_FFT

//#define SAVE_PLOT

CWB::mdc* MDC;

void Draw_LAL_vs_CWB() {

  MDC = new CWB::mdc();

  vector<mdcpar> par;

  if(TYPE_MDC==MDC_SGE || TYPE_MDC==MDC_CGE || TYPE_MDC==MDC_SGE_LAL) {
    par.resize(3);
    par[0].name="frequency"; par[0].value=831.3633225811645;
    par[1].name="Q";         par[1].value=6.239809618564323;
    par[2].name="decimals";  par[2].value=2;
  }

  if(TYPE_MDC==MDC_WNB || TYPE_MDC==MDC_WNB_LAL) {
    par.resize(6);
    par[0].name="frequency"; par[0].value=539.1519966628402;
    par[1].name="bandwidth"; par[1].value=196.8777744751424;
    par[2].name="duration";  par[2].value=0.07535462749423459;
    par[3].name="lseed";     par[3].value=1549043712;
    par[4].name="hseed";     par[4].value=523230030;
    par[5].name="decimals";  par[5].value=2;
  }


  MDC_TYPE mdc_type = TYPE_MDC;

  mdcid wf_mdcid = MDC->AddWaveform(mdc_type, par);
  TString wf_name = wf_mdcid.name;

  MDC->Print();

  // Get the first waveform hp,hx components 
  // NOTE : The current list contains only 1 waveform
  waveform wf = MDC->GetWaveform(wf_name,0);

  if(wf.status==false) {
    cout << "Error : Waveform " << wf_name << " not exist in the MDC pool !!!" << endl << endl;
    gSystem->Exit(1);
  }

  cout << "size : " << wf.hp.size() << " rate : " << wf.hp.rate() 
       << " start : " << (int)wf.hp.start() << endl;

  wf.hp.start(0);	// set start to 0 (needed by draw Method)
  wf.hx.start(0);

  watplot* plot = NULL;

#ifdef PLOT_TIME
  plot=MDC->Draw(wf.hp,MDC_TIME,"ALP ZOOM");
  MDC->Draw(wf.hx,MDC_TIME,"same",kRed);
#endif

#ifdef PLOT_FFT
  plot = MDC->Draw(wf.hp,MDC_FFT,"ALP ZOOM NOLOGX");	// draw hp in frequency domain
  MDC->Draw(wf.hx,MDC_FFT,"same NOLOGX",kRed);
  plot->graph[0]->GetHistogram()->GetXaxis()->SetRangeUser(0,1100);
#endif

  if(plot) plot->graph[0]->SetTitle(wf_name);

#ifdef SAVE_PLOT
  if(plot) plot->canvas->SaveAs("Plot_LAL_vs_CWB.png");
#endif

}
