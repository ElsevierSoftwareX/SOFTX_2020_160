{
  cout << "-----> plugins/CWB_Plugin_TestClassCBC_Config.C" << endl;

  CWB::mdc MDC(net);

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 60.0 --time-interval 3 ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--gps-start-time 931072130 --gps-end-time 933491330 ";
  inspOptions+= "--d-distr volume --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 10.000000 ";
  inspOptions+= "--min-mass1 25.000000 --max-mass1 225.000000 ";
  inspOptions+= "--min-mass2 25.000000 --max-mass2 225.000000 ";
  inspOptions+= "--min-mtotal 50.000000 --max-mtotal 250.000000 ";
  inspOptions+= "--min-mratio 0.25 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 1000000.0 --max-distance 1500000.0 ";
  inspOptions+= "--waveform EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  MDC.SetInspiral("EOBNRv2",inspOptions);

}
