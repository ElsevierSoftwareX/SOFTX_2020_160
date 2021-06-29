// THIS IS THE NEW VERSION OF CONFIG PLUGIN !!!
// THE NEW VERSION IS REQUIRED BY ROOT6
// IN ROOT5 IS STILL POSSIBLE TO USE OLD & NEW VERSION

#include "CWB_Plugin.h"

void CWB_PluginConfig() {

  // this config plugin must be used with the CWB_Plugin_MDC_OTF.C plugin 
  // to generate a uniform in volume distribution of BBH waveforms

  cout << "Execute : CWB_Plugin_BBH_low_UnVol_Config_ER8.C" << endl;

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int* gIFACTOR;
  CWB_PLUGIN_IMPORT(int*,gIFACTOR);

  int* xstart;
  int* xstop;
  CWB_PLUGIN_IMPORT(int*,xstart);
  CWB_PLUGIN_IMPORT(int*,xstop);

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
  // --output "file.xml" : write mdc injection's parameters to xml file. default=/tmp   
  // WARNING : write a space characters at the end of each inspOptions line !!!
  // ---------------------------------------------------------------------------

  TString inspOptions="";
  inspOptions = "--gps-start-time 1124005012 --gps-end-time 1125234224 ";
  inspOptions+= "--taper-injection start --seed "+TString(SEED)+" ";
  inspOptions+= "--dir "+TString((*cfg)->nodedir)+" ";
  inspOptions+= "--f-lower 24.000000 ";
  inspOptions+= "--time-step 30 ";

  inspOptions+= "--waveform IMRPhenomBthreePointFivePN ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--i-distr uniform ";
  inspOptions+= "--min-mass1 15. --max-mass1 25. ";
  inspOptions+= "--min-mass2 15. --max-mass2 25. ";
  inspOptions+= "--min-mtotal 30. --max-mtotal 50 ";
  inspOptions+= "--m-distr componentMass ";
  inspOptions+= "--enable-spin --aligned ";
  inspOptions+= "--min-spin1 0. --max-spin1 0.9 ";
  inspOptions+= "--min-spin2 0. --max-spin2 0.9 ";
  inspOptions+= "--amp-order 0 ";
  inspOptions+= "--d-distr volume ";
  //inspOptions+= "--min-z 1.e-4 --max-z 0.33 ";
  inspOptions+= "--min-distance 100 --max-distance 1100000 ";

  // the first parameter a the name of MDC defined by the user
  MDC->SetInspiral("IMRPhenomBthreePointFivePN",inspOptions);

}
