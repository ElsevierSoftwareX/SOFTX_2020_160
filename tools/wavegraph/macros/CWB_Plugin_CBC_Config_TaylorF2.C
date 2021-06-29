{
  cout << "-----> CWB_Plugin_NN_Config.C" << endl;

  CWB::mdc MDC(net);

  // -------------------------------->nodedir)+" ");
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 70.0 --time-interval 3 ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--dir "+TString(cfg->nodedir)+" ";
  inspOptions+= "--gps-start-time 931081124 --gps-end-time 935147777 ";
  inspOptions+= "--d-distr volume --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 40.000000 ";
  // inspOptions+= "--min-mass1 2.5 --max-mass1 10.0 ";  // incompatible with --m-distr totalMassRatio with lalsuite v6.13.2
  // inspOptions+= "--min-mass2 2.5 --max-mass2 10.0 ";  // incompatible with --m-distr totalMassRatio with lalsuite v6.13.2
  inspOptions+= "--min-mtotal 3.75 --max-mtotal 20. ";
  inspOptions+= "--min-mratio 0.5 --max-mratio 1.0 ";
  inspOptions+= "--min-distance 30000.0 --max-distance 30000.0 ";
  //inspOptions+= "--waveform TaylorT1newtonian --disable-spin ";
  inspOptions+= "--waveform TaylorT2threePointFivePN --disable-spin ";
  //inspOptions+= "--waveform EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  MDC.SetInspiral("TaylorF2",inspOptions);
  //MDC.SetInspiral("EOBNRv2",inspOptions);

}
