CWB::mdc* MDC; 
watplot* WTS;
watplot* plot;

void DrawWDM() {
  //
  // Instantiation of object CWB::mdc 
  // Get Waveform , WDM Transform - Plot WDM Scalogram using watplot
  // Author : Gabriele Vedovato

  //#define WNB_MDC
  //#define RD_MDC
  #define SGQ_MDC

  #define WDM_TF
  #define WDM_XCHECK

  MDC = new CWB::mdc; 

#ifdef WNB_MDC
  wavearray<double> x = MDC->GetWNB(1000., 1000., 0.01);
#endif
#ifdef RD_MDC
  wavearray<double> x = MDC->GetRD(1000., 0.2, 10.);
#endif
#ifdef SGQ_MDC
  wavearray<double> x = MDC->GetSGQ(1053., 9);
#endif

  //MDC->Draw(x);		// Draw Waveform in time domain
  //MDC->Draw(x,MDC_FFT);	// Draw Waveform in frequency domain
  //MDC->Draw(x,MDC_TF);		// Draw Waveform in time/frequency domain

  // compute x energy
  double E=0;
  double dt=1./x.rate();
  for(int i=0;i<x.size();i++) E+=x[i]*x[i]*dt;
  cout << "E = " << E << endl;

  // perform WDM transform level=512
  double R = x.rate();
  int level = 512;

  double dF = R/(2*level);	// frequency resolution
  double dT = level/R;		// time resolution

  WDM<double> wdtf(level, level, 4, 10);
  WSeries<double> w;
  w.Forward(x, wdtf);

  cout << "WDM Decomposition Level : " << w.getLevel() 
       << " dF : " << dF << " Hz " << " dT : " << 1000*dT << " ms " <<endl;

  // Compute Energy 
  WDM<double>* wdm = (WDM<double>*) w.pWavelet;
  int M = w.getLevel();
  int L = wdm->maxLayer();
  cout << "WDM Decomposition Level : " << M << " Layers " << L << endl; 

  int mF = int(w.size()/wdm->nSTS);
  int nTC = w.size()/(M+1)/mF;                   // # of Time Coefficients
  // Extract map00 & map90 : size = (M+1)*nTC
  double* map00 = wdm->pWWS;
  double* map90 = map00 + (mF-1)*(M+1)*nTC;
  
  wavearray<double> xx;
  wavearray<double> XX;
  E=0;
  // this example show the corrispondence between map00,map90 and xx,XX
  for(int j=0;j<M+1;j++) {	// loop over the layers
    w.getLayer(xx,j);		// extract phase 00 @ layer j
    w.getLayer(XX,-j);		// extract phase 90 @ layer j
    // for phase 00 -> (xx[i] = map00[i*(M+1)+j])
    // for phase 90 -> (XX[i] = map90[i*(M+1)+j])
    for(int i=0;i<xx.size();i++) {
//      cout << j << " P00 " << xx[i] << " " << map00[i*(M+1)+j] 
//                << " P90 " << XX[i] << " " << map90[i*(M+1)+j] << endl;
    }
    for(int i=0;i<xx.size();i++) E+=xx[i]*xx[i];
    //for(int i=0;i<XX.size();i++) E+=XX[i]*XX[i];
  }
  cout << "E = " << E << endl;

#ifdef WDM_TF
  // Plot WDM Scalogram
  WTS = new watplot(const_cast<char*>("wtswrc"));
  //scalogram maps
  double start = w.start();
  double stop  = w.start()+w.size()/w.rate();
  double flow  = 64;
  double fhigh = 2048;
  WTS->plot(&w, 2, start, stop,const_cast<char*>("COLZ"));
  WTS->hist2D->GetYaxis()->SetRangeUser(flow, fhigh);
  // dump spectrum
  char fname[1024];
  sprintf(fname,"wdm_scalogram_%d_lev%d.root", int(w.start()),w.getLevel());
  cout << endl << "Dump WDM Scalogram : " << fname << endl << endl;
  WTS->canvas->Print(fname);
#endif

#ifdef WDM_XCHECK
  // Compare Signal before/after WDM Forward/Inverse Tranform
  w.Inverse();
  wavearray<double> y = w;

  plot = new watplot(const_cast<char*>("plot"),200,20,800,500);
  char gtitle[256];
  sprintf(gtitle,"WDM : Original(black) - After WDM Forward/Inverse Tranform(Red)");
  plot->gtitle(gtitle,"time(sec)","amplitude");
  plot->goptions(const_cast<char*>("alp"), 1, 0., 0.);
  // draw signal
  //y-=x;	// difference
  x >> *plot;
  y >> *plot;

  // save plot to file
  TString gfile="WDM_time_xcheck.png";
  *plot >> gfile;
#endif

}
