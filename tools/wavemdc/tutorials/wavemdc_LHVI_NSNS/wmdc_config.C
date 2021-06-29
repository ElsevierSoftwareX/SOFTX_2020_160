{
  #include "wat.hh"
  #include <vector>

  // -------------------------------------------------------- 
  // define job parameters
  // -------------------------------------------------------- 
  TString frDir       = "frames";
  TString frLabel     = "LHVI_NSNS";

  // define frame segment list
  TString segmentList = "Segments/NSNS_SegmentsList.txt";

  // select segments
  int jobmin  = 1;    // start segment
  int jobmax  = 100;    // end segment
  int jobstep =  2;    // frames per job 

  // -------------------------------------------------------- 
  // define network
  // -------------------------------------------------------- 
  int nIFO=4;
  TString ifo[nIFO]={"L1","H1","V1",""};

  detectorParams detParms[nIFO] = {
                          {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716), 0, (   -197.716 )},
                          {"H1", 46.4551, -119.408,  0.0, 0, ( +90-125.998), 0, (   -125.998 )},
                          {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},
                          {"I1", 14.4,      76.4,    0.0, 0, ( +90        ), 0, (      0.0   )}  // INDIGO
                         };

  detector* pD[NIFO_MAX];          // pointers to detectors

  for(int i=0; i<nIFO; i++) {
    if(ifo[i]!="") pD[i] = new detector(ifo[i].Data()); // built in detector
    else           pD[i] = new detector(detParms[i]);   // user define detector
  }

  CWB::mdc MDC(nIFO,pD); 

  // --------------------------------------------------------
  // define channel names
  // --------------------------------------------------------
  vector<TString> chName;

  // ---------------------------------
  // set inspiral parms
  // ---------------------------------

  TString inspOptions="";
  inspOptions = "--time-step 300.0 --time-interval 0 ";
  inspOptions+= "--gps-start-time 931158295 --gps-end-time 931258000 ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--d-distr uniform --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 16.000000 ";
  inspOptions+= "--min-mass1 1.4 --max-mass1 1.4 ";
  inspOptions+= "--min-mass2 1.4 --max-mass2 1.4 ";
  inspOptions+= "--min-mtotal 2.8 --max-mtotal 2.8 ";
  inspOptions+= "--min-mratio 1.0 --max-mratio 1.0 ";
  inspOptions+= "--min-distance 100000.0 --max-distance 100000.0 ";
  inspOptions+= "--waveform EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789";

  MDC.SetInspiral("NSNS",inspOptions);
}
