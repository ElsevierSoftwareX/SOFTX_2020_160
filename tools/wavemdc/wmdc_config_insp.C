{
  #include <vector>

  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "TestInsp";

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

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="";
//  inspOptions = "--gps-start-time 968650000 --gps-end-time 968660000 --time-step 60 --time-interval 5 ";
  inspOptions = "--time-step 40 --time-interval 5 ";
  inspOptions+= "--l-distr random ";
//  inspOptions+= "--l-distr fixed --longitude 30 --latitude 76 ";
  inspOptions+= "--d-distr uniform --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 10.000000 ";
  inspOptions+= "--min-mass1 25.000000 --min-mass2 25.000000 --max-mass1 25.000000 --max-mass2 25.000000 ";
  inspOptions+= "--min-mtotal 50.000000 --max-mtotal 50.000000 --min-mratio 1.000000 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 5000.000000 --max-distance 5000.000000 ";
//  inspOptions+= "--output GHLTV-UNIFORM_v1-968650000-10000.xml ";
  inspOptions+= "--waveform EOBNRv2HMpseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  MDC.SetInspiral(frLabel,inspOptions);

}
