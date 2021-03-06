
// get number of WDM layers (first & last bins have df/2)
inline int GetLayers(WSeries<double> *WS) {return 1<<WS->pWavelet->m_Level;}
// get number of samples in wavelet layer
inline int GetSlices(WSeries<double> *WS) {return WS->size()/(1<<WS->pWavelet->m_Level);}

// get frequency resolution (Hz)
inline float GetFreqResolution(WSeries<double> *WS) {return (float)WS->rate()/2/(float)(1<<WS->pWavelet->m_Level);}
// get time resolution (sec)
inline float GetTimeResolution(WSeries<double> *WS) {return (float)(1<<WS->pWavelet->m_Level)/(float)WS->rate();}

void PlotWSeries(WSeries<double> *WS, TString fname) {

  if(WS->pWavelet->m_WaveType==WDMT) {
    cout << "PlotWSeries Error : PlotWSeries don't works for WDMT wavelet types" << endl;
    gSystem->Exit(1);
  }

  TCanvas* canvas;

  canvas= new TCanvas("WSeries", "WDMT", 200, 20, 800, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);
  canvas->SetRightMargin(0.10);
  canvas->SetLeftMargin(0.10);
  canvas->SetBottomMargin(0.13);
  canvas->SetBorderMode(0);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);


  TH2F* hist2D;

  int layers = GetLayers(WS);            // numbers of frequency bins (first & last bins have df/2)
  int slices = GetSlices(WS);            // number of time bins
  float dt = GetTimeResolution(WS);      // time bin resolution (sec)
  float df = GetFreqResolution(WS);      // frequency bin resolution (Hz)
  int rate = int(1./dt);

  cout << "layers : " << layers << "\t slices : " << slices << "\t rate : " << rate
       << "\t dt : " << dt << "\t df : " << df << endl;

  int nT = slices;
  int nF = layers;
  cout << "Time Range : " << nT*dt << " Freq Range " << nF*df;

  hist2D=new TH2F("vSS","", nT, 0, nT*dt, nF, 0., nF*df);

  wavearray<double> wl;
  for(int j=1;j<=nF;j++) {
    WS->getLayer(wl,j-1);
    for(int i=1;i<=nT;i++) {
      float EE = wl.data[i-1];
      hist2D->SetBinContent(i,j,EE);
    }
  }

  hist2D->SetXTitle("time, sec");
  hist2D->SetYTitle("frequency, Hz");

  hist2D->SetStats(kFALSE);
  hist2D->SetTitleFont(12);
  hist2D->SetFillColor(kWhite);

  char title[256];
  sprintf(title,"Scalogram (%s)","energy");
  hist2D->SetTitle(title);

  hist2D->GetXaxis()->SetNdivisions(506);
  hist2D->GetXaxis()->SetLabelFont(42);
  hist2D->GetXaxis()->SetLabelOffset(0.014);
  hist2D->GetXaxis()->SetTitleOffset(1.4);
  hist2D->GetYaxis()->SetTitleOffset(1.2);
  hist2D->GetYaxis()->SetNdivisions(506);
  hist2D->GetYaxis()->SetLabelFont(42);
  hist2D->GetYaxis()->SetLabelOffset(0.01);
  hist2D->GetZaxis()->SetLabelFont(42);
  hist2D->GetZaxis()->SetNoExponent(false);
  hist2D->GetZaxis()->SetNdivisions(506);

  hist2D->GetXaxis()->SetTitleFont(42);
  hist2D->GetXaxis()->SetTitle("Time (sec)");
  hist2D->GetXaxis()->CenterTitle(true);
  hist2D->GetYaxis()->SetTitleFont(42);
  hist2D->GetYaxis()->SetTitle("Frequency (Hz)");
  hist2D->GetYaxis()->CenterTitle(true);

  hist2D->GetZaxis()->SetTitleOffset(0.6);
  hist2D->GetZaxis()->SetTitleFont(42);
  hist2D->GetZaxis()->CenterTitle(true);

  hist2D->GetXaxis()->SetLabelSize(0.03);
  hist2D->GetYaxis()->SetLabelSize(0.03);
  hist2D->GetZaxis()->SetLabelSize(0.03);

  canvas->cd();

  hist2D->Draw("COLZ");

  // dump spectrum
  cout << endl << "Dump SS Scalogram : " << fname.Data() << endl << endl;
  canvas->Print(fname);

  delete hist2D;
  delete canvas;
}

