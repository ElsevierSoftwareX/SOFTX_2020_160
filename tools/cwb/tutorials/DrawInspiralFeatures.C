{
  //
  // Show how to use Toolbox to get signal envelop / instantaneous frequency, WignerVille Transform, ...
  // Author : Gabriele Vedovato

  #define DRAW_WIGNER_VILLE
  //#define DRAW_IFREQUENCY
  //#define DRAW_ENVELOPE
  //#define DRAW_ENVELOPE_EXTREME
  //#define DRAW_FREQUENCY_EXTREME


  #include <vector>

  CWB::mdc MDC; 

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------
  TString inspOptions="";
  inspOptions = "--time-step 60.0 --time-interval 3 ";
  inspOptions+= "--dir . ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--gps-start-time 931072130 --gps-end-time 933491330 ";
  inspOptions+= "--d-distr volume --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 10.000000 ";
  inspOptions+= "--min-mass1 25.000000 --max-mass1 225.000000 ";
  inspOptions+= "--min-mass2 25.000000 --max-mass2 225.000000 ";
  inspOptions+= "--min-mtotal 50.000000 --max-mtotal 250.000000 ";
  inspOptions+= "--min-mratio 0.25 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 1000000.0 --max-distance 1500000.0 ";
  inspOptions+= "--approximant EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789 ";
  inspOptions+= "--output inspirals.xml ";		// set output xml file

  // set and write xml file
  MDC.SetInspiral("EOBNRv2",inspOptions);

  // set inspiral using xml file (speedup GetInspiral method)
  TString inspOptions="--xml inspirals.xml";
  MDC.SetInspiral("EOBNRv2",inspOptions);

  // Get the first waveform hp,hx components starting from gps = 931072130
  wavearray<double> hp = MDC.GetInspiral("hp",931072130,931072230);
  wavearray<double> hx = MDC.GetInspiral("hx",931072130,931072230);
  cout << "size : " << hp.size() << " rate : " << hp.rate() << " start : " << (int)hp.start() << endl;
  hp.start(0);		// set start to 0 (needed by draw Method)
  hx.start(0);
  MDC.Draw(hp,MDC_TIME);
  //MDC.Draw(hx,MDC_TIME,"same",kRed);

  //MDC.Draw(hp,MDC_FFT);	// draw hp in frequency domain
  //MDC.Draw(hp,MDC_TF);        // draw hp in time-frequency domain;

#ifdef DRAW_ENVELOPE
  // Get hp envelope
  wavearray<double> ehp = CWB::Toolbox::getHilbertEnvelope(hp);
  MDC.Draw(ehp,MDC_TIME,"same",kRed);
#endif

#ifdef DRAW_IFREQUENCY
  // Get hp instantaneous frequency
  wavearray<double> fhp = CWB::Toolbox::getHilbertIFrequency(hp);
  MDC.Draw(fhp,MDC_TIME);
#endif

#ifdef DRAW_WIGNER_VILLE
  // draw Wigner-Ville transform
  wavearray<double> xhp(hp.size()/64);
  xhp.rate(hp.rate()/64); xhp.start(hp.start());
  for(int i=0;i<xhp.size();i++) xhp[i]=hp[i*64];    // decimated hp
  wavearray<double> wvhp = CWB::Toolbox::getWignerVilleTransform(xhp);
  int N = sqrt(wvhp.size());
  double dt = 1./wvhp.rate();
  double df = wvhp.rate()/N;

  TH2D* h2 = new TH2D("h2","h2",N,0,N*dt,N,0,N*df);
  h2->SetStats(kFALSE);
  h2->SetTitle("Wigner Ville Transform");
  h2->GetXaxis()->SetTitle("Time (sec)");
  h2->GetYaxis()->SetTitle("freq (Hz)");

  for(int n=0; n<N; n++) {
   for(int m=0; m<N; m++) {
      if(wvhp[n*N+m]>0) h2->Fill(m*dt,n*df,wvhp[n*N+m]);
    }
  }

  h2->Draw("colz");
#endif

  // Draw envelope/frequency extreme of signal
  wavearray<double> t(hp.size());
  wavearray<double> f(hp.size());
  wavearray<double> a(hp.size());
  int size=0;
  double dt = 1./hp.rate();
  double amplitude,omega,phase;
  for (int i=1;i<hp.size()-1;i++) {
    CWB::Toolbox::getSineFittingParams(hp[i-1],hp[i],hp[i+1],hp.rate(),amplitude,omega,phase);
    if((amplitude!=0)&&(omega!=0)) {
      t[size]=dt*(i-1)-(phase/omega);
      a[size]=fabs(amplitude);
      f[size]=omega/TMath::TwoPi();
      size++;
    }
  }

#ifdef DRAW_ENVELOPE_EXTREME
  TGraph* agr = new TGraph(size,t.data,a.data);
  agr->SetLineColor(kBlue);
  agr->SetMarkerColor(kBlue);
  agr->SetMarkerStyle(20);
  agr->SetMarkerSize(1);
  fgr->SetTitle("Envelope");
  agr->GetHistogram()->GetXaxis()->SetTitle("Time (sec)");
  agr->GetHistogram()->GetYaxis()->SetTitle("amplitude (Hz)");
  //agr->Draw("ALP");
  agr->Draw("SAME");
#endif

#ifdef DRAW_FREQUENCY_EXTREME
  TGraph* fgr = new TGraph(size,t.data,f.data);
  fgr->SetLineColor(kRed);
  fgr->SetMarkerColor(kRed);
  fgr->SetMarkerStyle(20);
  fgr->SetMarkerSize(1);
  fgr->SetTitle("Frequency");
  fgr->GetHistogram()->GetYaxis()->SetRangeUser(0,256);
  fgr->GetHistogram()->GetXaxis()->SetTitle("Time (sec)");
  fgr->GetHistogram()->GetYaxis()->SetTitle("freq (Hz)");
  fgr->Draw("ALP");
  //fgr->Draw("SAME");
#endif

}
