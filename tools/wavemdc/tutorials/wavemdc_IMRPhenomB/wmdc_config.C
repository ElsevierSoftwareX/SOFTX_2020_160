{
  #include <vector>

  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "wavemdc_IMRPhenomB";

  // define frame segment list
  TString segmentList = "Segments/S6a_LS-segs.txt";

  // select segments
  int jobmin  = 210;    // start segment
  int jobmax  = 312;    // end segment
  int jobstep =  10;    // frames per job 

  // -------------------------------------------------------- 
  // define network
  // -------------------------------------------------------- 
  int nIFO=3;
  TString ifo[nIFO]={"L1","H1","V1"};
  CWB::mdc MDC(nIFO,ifo); 

  // --------------------------------------------------------
  // define channel names
  // --------------------------------------------------------
  vector<TString> chName;

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 60 --time-interval 3 ";
  inspOptions+= "--d-distr volume --l-distr random --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 8.000000 ";
  inspOptions+= "--min-mass1 25.0 --max-mass1 400.0 ";
  inspOptions+= "--min-mass2 25.0 --max-mass2 400.0 ";
  inspOptions+= "--min-mtotal 50.0 --max-mtotal 450.0 ";
  inspOptions+= "--min-mratio 0.25 --max-mratio 1.0 ";
  inspOptions+= "--min-distance 100000.0 --max-distance 150000.0 ";
  inspOptions+= "--waveform IMRPhenomBpseudoFourPN ";
  inspOptions+= "--enable-spin --aligned ";
  inspOptions+= "--min-spin1 0.000000 --max-spin1 0.8 ";
  inspOptions+= "--min-spin2 0.000000 --max-spin2 0.8 ";
  inspOptions+= "--taper-injection start --seed 123451234";

  MDC.SetInspiral("IMRPhenomBpseudoFourPN",inspOptions);

}
