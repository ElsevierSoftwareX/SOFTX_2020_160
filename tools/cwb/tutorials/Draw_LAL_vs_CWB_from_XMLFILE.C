//
// Show comparison between SGE/WNB CWB vs LAL waveforms 
// The waveform parameters are loaded from LAL XML file
//
// MDC_SGE_LAL use CG instead of SG !!!
// MDC_WNB_LAL has a gaussian envelope also in frequency
//
// Author : Gabriele Vedovato
//

/*
// SGE
#define XMLFILE_NAME	"SG_LAL.xml"
// defines the gps range used to fill the MDC pool
#define GPS_START	 968386896
#define GPS_STOP	 968387096
*/


// WNB
#define XMLFILE_NAME	"WNB_LAL.xml"
// defines the gps range used to fill the MDC pool
#define GPS_START	 966769371
#define GPS_STOP	 966769571


//#define PLOT_TIME
#define PLOT_FFT

// defines the ID used to extract the waveform from the MDC pool
// see ID from the MDC->Print() list
#define MDC_ID	1

// defines the MDC engine used to generate the waveforms
// the parameters are read from XMLFILE_NAME 
#define MDC_ENGINE	0	// 0/1 -> LAL/CWB

//#define SAVE_PLOT

CWB::mdc* MDC;

void Draw_LAL_vs_CWB_from_XMLFILE() {

  MDC = new CWB::mdc();

  vector<mdcpar> par;

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  par.resize(4);
  par[0].name="start";   par[0].value=GPS_START;
  par[1].name="stop";    par[1].value=GPS_STOP;
  par[2].name="engine";  par[2].value=MDC_ENGINE; 
  par[3].name="auto";    par[3].value=1; 

  TString xmlfile = XMLFILE_NAME;

  // fill the MDC waveform pool 
  MDC->SetSkyDistribution(MDC_XMLFILE,xmlfile,par,0);
  // print the MDC waveform pool 
  MDC->Print();
  // Get the MDC_ID waveform hp,hx components from the MDC pool
  waveform wf = MDC->GetWaveform(MDC_ID,0);

  if(wf.status==false) {
    cout << "Error : ID is not presents in the MDC pool !!!" << endl << endl;
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

  if(plot) plot->graph[0]->SetTitle(wf.name);

#ifdef SAVE_PLOT
  if(plot) plot->canvas->SaveAs("Plot_LAL_vs_CWB_from_XMLFILE.png");
#endif

}
