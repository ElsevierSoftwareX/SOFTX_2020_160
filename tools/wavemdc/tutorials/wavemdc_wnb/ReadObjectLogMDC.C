{

  TFile *jfile = TFile::Open("frames/L1H1V1-TestWNB-Log.root");

  jfile->ls();

  // read mdc object
  CWB::mdc MDC = *(CWB::mdc*)jfile->Get("TestWNB");

  MDC.Print();

  MDC.Draw("WNB1000_100_0d100", 0, "hp", MDC_DRAW_TF);
  //MDC.Draw("WNB1000_1000_0d001", 0, "hp", MDC_DRAW_TIME);
  //MDC.Draw("WNB1000_1000_0d001", 0, "hx", MDC_DRAW_TIME,"SAME",2);
}
