{
  //
  // Read & Draw Waveforms 
  // Author : Gabriele Vedovato

  CWB::mdc MDC; 

  wavearray<double> x;
  //MDC.ReadWaveform(x,"Waveforms/s11.2.h.dat");
  //MDC.ReadWaveform(x,"Waveforms/data.GRW.md.h10kpc.15_3.2.txt");
  MDC.ReadWaveform(x,"Waveforms/SG554Q8d9.txt",16384.);

  cout << x.size() << " " << x.start() << " " << x.rate() << endl;

  double hrss=0;
  for(int i=0;i<x.size();i++) hrss+=x[i]*x[i];
  hrss=sqrt(hrss/x.rate());
  for(int i=0;i<x.size();i++) x[i]/=hrss;


  //MDC.Dump("TEST_WF_DUMP.txt",x);

  //MDC.Draw(x);
  //MDC.Draw(x,MDC_FFT);
  MDC.Draw(x,MDC_TF);

  //exit(0);
}
