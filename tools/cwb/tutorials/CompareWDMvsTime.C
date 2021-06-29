// Read waveform from file
// resample data
// WDM transfom
// compare energy in time domain vs energy in time-freq domain
// plot WDM
// plot waveform in time/fft domain 
// save plots to disk
// compare Time/WDM energies for each time bin

#define WF_HP	"Wolfy_hplus.txt" 	// Wolfy waveform
#define NZEROES	50000

#define USE_WDM_00_90			// use 00 + 90 components
					// if commented then only 00 component is used 

//#define SAVE_PLOTS			// save plots to disk

gwavearray<double>* gx;
gwavearray<double>* geT;
gwavearray<double>* greT;
gwavearray<double>* geWDM;
watplot* WTS;

void CompareWDMvsTime() {

  TString hp = WF_HP;
 
  // ---------------------------------------------------
  // read waveform and fill x wavearray
  // ---------------------------------------------------
  wavearray<double> x;
  wavearray<double> time;

  double X,T;
  ifstream in;
  in.open(hp.Data());
  for(int i=0;i<NZEROES;i++) x.append(0);	// padding with zeroes
  while(1) {
    in >> X >> T;
    if(!in.good()) break;
    time.append(X);
    x.append(T);
  }
  for(int i=0;i<NZEROES;i++) x.append(0);	// padding with zeroes

  double rate = (time.size()-1)/(time[time.size()-1]-time[0]);
  x.rate(rate);
  cout << "rate time series  " << x.rate() << endl;
  cout << "size time series  " << x.size() << endl;
  cout << "start time series " << x.start()<< endl;
  cout << "stop time series  " << x.stop() << endl;

  x.resample(8192);				// resample data to 8192 Hz

  // ---------------------------------------------------
  // compute energy in time domain
  // ---------------------------------------------------
  double E_TIME=0;
  double dt=1./x.rate();
  geT = new gwavearray<double>(x);			// energy array
  for(int i=0;i<x.size();i++) {geT->data[i]=x[i]*x[i]*dt;E_TIME+=geT->data[i];}
  cout << "E_TIME = " << E_TIME << endl; 

  // ---------------------------------------------------
  // perform WDM transform level=512
  // ---------------------------------------------------
  double R  = x.rate();
  int level = 32;

  double dF = R/(2*level);	// frequency resolution
  double dT = level/R;		// time resolution

  cout<< "R="<< R << " dF= " << dF << " dT= " << dT << endl;

  WDM<double> wdtf(level, level, 4, 10);
  WSeries<double> w; 
  w.Forward(x, wdtf); 		// WDM Transform

  cout << "size WDM  " << w.size() << endl;
  cout << "start WDM " << w.start() << endl;
  cout << "stop WDM  " << w.stop() << endl;

  // ---------------------------------------------------
  // compute energy in time-freq domain
  // ---------------------------------------------------
  int slices = w.sizeZero();
  float df2  = w.resolution();              
  float dt2  = 1./(2*df2);        
  int layers = w.maxLayer()+1;
  double E_WDM=0;
  geWDM = new gwavearray<double>(slices); 	// energy array
  geWDM->rate(w.rate()/level);geWDM->start(w.start());	
  for(int i=1;i<=slices;i++) {
    geWDM->data[i-1] = 0;
    for(int j=1;j<=layers;j++) { // loop over the frequencies with the same time
#ifdef USE_WDM_00_90
      // en = (E00+E90)/2
      double en = pow(w.getSample(i-1,j-1),2)+pow(w.getSample(i-1,-(j-1+0.01)),2);
      en/=2*R;
#else 
      // en = E00
      double en = pow(w.getSample(i-1,j-1),2)/R;
#endif
      geWDM->data[i-1] += en;
      E_WDM += en;
    }
  }
  cout << "slices "<< slices << " layers"<< layers << endl;
  cout << "E_WDM = " << E_WDM << endl; 
  cout << "E_WDM/E_TIME = " << E_WDM/E_TIME << endl; 

  // ---------------------------------------------------
  // Plot Waveform in time/fft domain
  // ---------------------------------------------------
  gx = new gwavearray<double>(&x);
  gx->Draw(GWAT_TIME);
  //gx->Draw(GWAT_TIME,"FULL");
  //gx->Draw(GWAT_FFT);

#ifdef SAVE_PLOTS
  // save plot to file
  watplot* plot = gx->GetWATPLOT();	// get pointer to watplot object
  char gtitle[256]; sprintf(gtitle,"Wolfy waveform");
  plot->gtitle(gtitle,"time(sec)","amplitude");
  //plot->gtitle(gtitle,"frequency(hz)","amplitude");
  TString gfile="Wolfy_Time.png";
  (*plot) >> gfile;
#endif

  // ---------------------------------------------------
  // Plot Waveform in time-freq domain 
  // Note : Energy in the plot is not normalize : (Norm is 1/R)
  // ---------------------------------------------------
  WTS = new watplot(const_cast<char*>("wtswrc"));
  double start,stop;
  gx->GetTimeRange(start,stop);
  double flow  = 64;
  double fhigh = 4096;
#ifdef USE_WDM_00_90
  WTS->plot(&w, 2, start, stop,const_cast<char*>("COLZ"));
#else
  WTS->plot(&w, 4, start, stop,const_cast<char*>("COLZ"));
#endif
  WTS->hist2D->GetYaxis()->SetRangeUser(flow, fhigh);
#ifdef SAVE_PLOTS
  WTS->canvas->SaveAs("Wolfy_WDM.png"); // dump plot
#endif


  // ---------------------------------------------------
  // Compare energy eT vs eWDM
  // ---------------------------------------------------
  greT = new gwavearray<double>(geWDM->size());	// energy array
  greT->start(geWDM->start());		 
  greT->rate(geWDM->rate());
  // energy of time domain signal is integrated over the WDM bin size (=level)
  for(int i=0;i<greT->size()-1;i++) {
    greT->data[i]=0;
    for(int j=0;j<level;j++) greT->data[i]+=geT->data[i*level+j];
  }

//  geT->Draw(GWAT_TIME);
//  geWDM->Draw(GWAT_TIME);		
  greT->Draw(GWAT_TIME);			// black is greT
  // the left edge of the bin is shifted in time by +half_sample because  
  // the WDM pixel time is the central time of the pixel (the same is for the frequency)		 
  geWDM->start(geWDM->start()-1./(2.*geWDM->rate()));	
  greT->Draw(geWDM,GWAT_TIME,"SAME",kRed);	// red is geWDM
#ifdef SAVE_PLOTS
  // save plot to file
  plot = greT->GetWATPLOT();			// get pointer to watplot object
  gfile="Wolfy_Time_Freq_Comp.png";
  (*plot) >> gfile;
#endif
}
