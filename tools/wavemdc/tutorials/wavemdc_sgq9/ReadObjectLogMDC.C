{

  TFile *jfile = TFile::Open("frames/L1H1V1-TestWaveMDC-Log.root");

  jfile->ls();

  // read mdc object
  CWB::mdc MDC = *(CWB::mdc*)jfile->Get("TestWaveMDC");

  MDC.PrintWaveformsList();

  MDC.Draw("SG945Q9", 0, "hp", MDC_DRAW_TIME);
}
