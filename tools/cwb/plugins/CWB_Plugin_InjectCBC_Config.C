{
  //!NOISE_MDC_SIMULATION
  // Config Plugin to injected 'on the fly' CBC MDC

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  // ---------------------------------
  // set LAL inspiral parameters
  // ---------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 33.333333 --time-interval 10 ";
  inspOptions+= "--dir "+TString(cfg->nodedir)+" ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--gps-start-time 931158092 --gps-end-time 931999916 ";
  inspOptions+= "--d-distr uniform --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 10.000000 ";
  inspOptions+= "--min-mass1 25.000000 --max-mass1 25.000000 ";
  inspOptions+= "--min-mass2 25.000000 --max-mass2 25.000000 ";
  inspOptions+= "--min-mtotal 50.000000 --max-mtotal 50.000000 ";
  inspOptions+= "--min-mratio 1.000000 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 5000.000000 --max-distance 5000.000000 ";
  inspOptions+= "--waveform EOBNRv2HMpseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  // Injections can be read from xml file
  // TString inspOptions="--xml injections.xml";

  MDC->SetInspiral("inspNameTEST",inspOptions);

}
