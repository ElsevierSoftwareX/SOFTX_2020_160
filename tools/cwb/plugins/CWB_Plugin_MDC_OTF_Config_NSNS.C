// THIS IS THE NEW VERSION OF CONFIG PLUGIN !!!
// THE NEW VERSION IS REQUIRED BY ROOT6
// IN ROOT5 IS STILL POSSIBLE TO USE OLD & NEW VERSION

#include "CWB_Plugin.h"

void CWB_PluginConfig() {

  //!NOISE_MDC_SIMULATION
  // Config Plugin to generate injected 'on the fly' NSNS from LAL

  cout << "Execute CWB_Plugin_MDC_OTF_Config_NSNS.C ..." << endl;                                 

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int* gIFACTOR;
  CWB_PLUGIN_IMPORT(int*,gIFACTOR);

  int seed = 100*((*net)->nRun)+*gIFACTOR; 
  char SEED[64]; sprintf(SEED,"%d",seed);

  cout << "gIFACTOR : " << *gIFACTOR << " net->nRun : " << (*net)->nRun << " SEED " << SEED << endl; 

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  CWB::config** cfg;
  CWB_PLUGIN_IMPORT(CWB::config**,cfg);


  // ---------------------------------------------------------------------------
  // set inspiral parms
  // waveforms are produced by the LAL library
  // to dump all the available options do:
  // $LALINSPINJ_EXEC --help
  // for any details refer to the LAL documentation.
  // there are some special options added only for the mdc class
  // --approximant       : is used as alternative to --waveform to force 
  //                       the use of the new XLALSimInspiralChooseWaveformFromSimInspiral
  // --output "file.xml" : write mdc injection's parameters to xml file   
  // --dir "tmp dir"     : directory used to store the temporary xml file, default=/tmp         
  // WARNING : write a space characters at the end of each inspOptions line !!!
  // ---------------------------------------------------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 600.0 --time-interval 0 ";
  inspOptions+= "--gps-start-time 931158300 --gps-end-time 931158700 ";
  inspOptions+= "--dir "+TString((*cfg)->nodedir)+" ";
  inspOptions+= "--l-distr fixed --longitude 45 --latitude 45 ";
  inspOptions+= "--d-distr uniform --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 32.000000 ";
//OLD_LAL  inspOptions+= "--min-mass1 1.4 --max-mass1 1.4 ";
//OLD_LAL  inspOptions+= "--min-mass2 1.4 --max-mass2 1.4 ";
  inspOptions+= "--min-mtotal 2.8 --max-mtotal 2.8 ";
  inspOptions+= "--min-mratio 1.0 --max-mratio 1.0 ";
  inspOptions+= "--min-distance 200000.0 --max-distance 200000.0 ";
  inspOptions+= "--waveform EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  // the first parameter a the name of MDC defined by the user
  MDC->SetInspiral("NSNS",inspOptions);

}

