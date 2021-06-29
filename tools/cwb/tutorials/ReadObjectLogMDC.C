{
  //
  // Draw Waveform from mdc object saved to root file
  // Author : Gabriele Vedovato

  TFile *jfile = TFile::Open("frames/L1H1V1-RDE_TST1-Log.root");
  if(jfile==NULL) {
    cout << "Input ROOT File  not exist !!!" << endl;
    gSystem->Exit(1);
  }

  jfile->ls();

  // read mdc object
  CWB::mdc MDC = *(CWB::mdc*)jfile->Get("RDE_TST1");

  MDC.Print();

  //MDC.Draw("RDE_1590_0d2", 0, "hp", MDC_TF);
  MDC.Draw("RDE_2090_0d2", 0, "hp", MDC_TIME);
  MDC.Draw("RDE_2090_0d2", 0, "hx", MDC_TIME,"SAME",2);
}
