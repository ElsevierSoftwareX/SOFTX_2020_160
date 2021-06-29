{
  cout << "-----> plugins/CWB_Plugin_TestClassCBC_Config.C" << endl;

  CWB::mdc MDC(net);

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="";
  //inspOptions = "--time-step 100.0 --time-interval 0 ";
  inspOptions = "--time-step 100.0 --time-interval 0 ";
  inspOptions+= "--gps-start-time 931158300 --gps-end-time 931158700 ";
  inspOptions+= "--dir "+TString(cfg->nodedir)+" ";
  //inspOptions+= "--l-distr random ";
  inspOptions+= "--l-distr fixed --longitude 45 --latitude 45 ";
  //inspOptions+= "--polarization 0 ";
  inspOptions+= "--d-distr uniform --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 32.000000 ";
  inspOptions+= "--min-mass1 1.4 --max-mass1 1.4 ";
  inspOptions+= "--min-mass2 1.4 --max-mass2 1.4 ";
  inspOptions+= "--min-mtotal 2.8 --max-mtotal 2.8 ";
  inspOptions+= "--min-mratio 1.0 --max-mratio 1.0 ";
  inspOptions+= "--min-distance 200000.0 --max-distance 200000.0 ";
  inspOptions+= "--waveform EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  MDC.SetInspiral("NSNS",inspOptions);
}
