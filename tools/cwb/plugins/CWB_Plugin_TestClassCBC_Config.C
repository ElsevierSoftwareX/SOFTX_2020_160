{
  //!NOISE_MDC_SIMULATION
  // Config Plugin to generate simulated gaussian noise and injected 'on the fly' CBC MDC

  cout << "-----> plugins/CWB_Plugin_TestClassCBC_Config.C" << endl;

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 33.333333 --time-interval 10 ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--gps-start-time 931158092 --gps-end-time 931999916 ";
  inspOptions+= "--dir "+TString(cfg->nodedir)+" ";
//  inspOptions+= "--l-distr fixed --longitude 30 --latitude 76 ";
  inspOptions+= "--d-distr uniform --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 10.000000 ";
  inspOptions+= "--min-mass1 25.000000 --max-mass1 25.000000 ";
  inspOptions+= "--min-mass2 25.000000 --max-mass2 25.000000 ";
  inspOptions+= "--min-mtotal 50.000000 --max-mtotal 50.000000 ";
  inspOptions+= "--min-mratio 1.000000 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 5000.000000 --max-distance 5000.000000 ";
  inspOptions+= "--waveform EOBNRv2HMpseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  MDC->SetInspiral("inspNameTEST",inspOptions);

}
