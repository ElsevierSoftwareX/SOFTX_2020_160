{
  #include <vector>

  // --------------------------------------------------------
  // define job parameters
  // --------------------------------------------------------
  TString frDir       = "frames";
  TString frLabel     = "wavemdc_EOBNRv2HMpseudoFourPN";

  // define frame segment list
  TString segmentList = "Segments/Chris-segs.txt";

  // select segments
  int jobmin  =   1;    // start segment
  int jobmax  =   1;    // end segment
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
  vector<TString> chName(nIFO);
  chName[0]="L1:FAKE-STRAIN_INJ";
  chName[1]="H1:FAKE-STRAIN_INJ";
  chName[2]="V1:FAKE-STRAIN_INJ";

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="--xml input/GHLTV-UNIFORM_v1-931072130-841.xml ";
  inspOptions+= "--approximant EOBNRv2HMpseudoFourPN ";
  //inspOptions+= "--dump Test.xml ";
  MDC.SetInspiral("EOBNRv2HMpseudoFourPN",inspOptions);

}
