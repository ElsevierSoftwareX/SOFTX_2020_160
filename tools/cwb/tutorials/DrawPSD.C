{
  #define ADV_LIGO_PSD    "$HOME_WAT/tools/cwb/plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt"
  #define BANDWIDTH    	8192.	// Hz
  #define DELTA_FREQ	0.125	// Hz

  //
  // generate Read PSD from file & display
  // Author : Gabriele Vedovato

  #include <vector>

  TString psdName  = gSystem->ExpandPathName(ADV_LIGO_PSD);

  CWB::Toolbox TB;

  double fWidth = BANDWIDTH;	// bandwidth
  double dFreq  = DELTA_FREQ;	// frequency resolution
  wavearray<double> psd = TB.GetDetectorPSD(psdName,fWidth,dFreq);
  double df=psd.rate()/psd.size();
  cout << "PSD : " << " rate : " << psd.rate() << " size : " << psd.size() << " df " << df << endl;

  // define frequency array
  wavearray<double> f=psd;
  for(int i=0;i<f.size();i++) f[i]=i*df;


  TObjArray* token = TString(psdName).Tokenize(TString('/'));
  TObjString* sname = (TObjString*)token->At(token->GetEntries()-1);
  TString Title = sname->GetString();

  // display PSD
  TCanvas canvas;
  canvas.SetLogx();
  canvas.SetLogy();
  canvas.SetGridx();
  canvas.SetGridy();
  TGraph gr(psd.size(),f.data,psd.data);
  gr.SetTitle(Title);
  gr.SetMarkerColor(kRed);
  gr.SetLineColor(kRed);
  gr.GetHistogram()->GetXaxis()->SetRangeUser(8,8192);
  gr.GetHistogram()->GetYaxis()->SetRangeUser(1e-24,1e-21);
  gr.GetHistogram()->GetXaxis()->SetTitle("Hz");
  gr.GetHistogram()->GetYaxis()->SetTitle("PSD");
  gr.GetHistogram()->GetYaxis()->SetTitleOffset(1.2);

  gr.Draw("ALP");
}

