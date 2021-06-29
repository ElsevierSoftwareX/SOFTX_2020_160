{

  TFile *jfile = TFile::Open("frames/L1H1V1-RDE_TST1-Log.root");

  jfile->ls();

  // read mdc object
  CWB::mdc MDC = *(CWB::mdc*)jfile->Get("RDE_TST1");

  MDC.PrintWaveformsList();

  //MDC.Draw("RDE_1590_0d2", 0, "hp", MDC_DRAW_TF);
  MDC.Draw("RDE_2090_0d2", 0, "hp", MDC_DRAW_TIME);
  MDC.Draw("RDE_2090_0d2", 0, "hx", MDC_DRAW_TIME,"SAME",2);
}
