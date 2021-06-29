{                                                                                                    
  //!NOISE_MDC_SIMULATION
  // Config Plugin to generate injected 'on the fly' EOBNRv2 from LAL

  cout << "Execute CWB_Plugin_MDC_OTF_Config_EOBNRv2pseudoFourPN.C ..." << endl;                                 

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  // ---------------------------------------------------------------------------
  // set inspiral parms
  // waveforms are produced by the LAL library
  // to dump all the available options do:
  // $LALINSPINJ_EXEC --help
  // for any details refer to the LAL documentation.
  // there are some special options added only for the mdc class
  // --approximant       : is used as alternative to --waveform to force 
  //                       the use of the new XLALSimInspiralChooseWaveformFromSimInspiral
  // --output "file.xml" : write mdc injection's parameters to xml file. default=/tmp   
  // WARNING : write a space characters at the end of each inspOptions line !!!
  // ---------------------------------------------------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 60.0 --time-interval 3 ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--gps-start-time 931081124 --gps-end-time 932377124 ";
  inspOptions+= "--dir "+TString(cfg->nodedir)+" ";
  inspOptions+= "--d-distr volume --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 5.000000 ";
  inspOptions+= "--min-mass1 170.000000 --max-mass1 680.000000 ";
  inspOptions+= "--min-mass2 170.000000 --max-mass2 680.000000 ";
  inspOptions+= "--min-mtotal 650.000000 --max-mtotal 850.000000 ";
  inspOptions+= "--min-mratio 0.25 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 1000000.0 --max-distance 1500000.0 ";
//  inspOptions+= "--waveform EOBNRv2pseudoFourPN ";	// use the old LALFindChirpInjectSignals
  inspOptions+= "--approximant EOBNRv2HMpseudoFourPN ";	// use the new XLALSimInspiralChooseWaveformFromSimInspiral
  inspOptions+= "--disable-spin ";
  inspOptions+= "--taper-injection start --seed 125431123 ";

  // if user wants provide its own xml injection list {injections.xml}
  // all previous declaration must be replaced with :
  /*
  TString inspOptions="";
  inspOptions+= "--xml injection.xml "
  inspOptions+= "--dir "+TString(cfg->nodedir)+" ";
  inspOptions += "--approximant EOBNRv2HMpseudoFourPN ";
  */

  // the first parameter a the name of MDC defined by the user
  MDC->SetInspiral("EOBNRv2HMpseudoFourPN",inspOptions);
}

