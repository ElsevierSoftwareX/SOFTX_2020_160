{
  //!NOISE_MDC_SIMULATION
  // Config Plugin to generate EOBNRv2 MDC OTF (used by Overlap example)

  cout << "-----> CWB_Plugin_EOBNRv2_Config.C" << endl;

  CWB::mdc MDC(net);

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 30.0 --time-interval 3 ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--gps-start-time 931158205 --gps-end-time 931158800 ";
  inspOptions+= "--d-distr volume --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 24.000000 ";
  inspOptions+= "--min-mass1 15.0 --max-mass1 25.0 ";
  inspOptions+= "--min-mass2 15.0 --max-mass2 25.0 ";
  inspOptions+= "--min-mtotal 30. --max-mtotal 50. ";
  inspOptions+= "--min-mratio 0.20 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 30000.0 --max-distance 45000.0 ";
  inspOptions+= "--waveform EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  MDC.SetInspiral("EOBNRv2",inspOptions);

}
