{
  //
  // Draw Inspiral Waveform with the mdc class (for user defined detectors)
  // Author : Gabriele Vedovato

  #include <vector>

  // --------------------------------------------------------
  // define network
  // --------------------------------------------------------

  int nIFO = 4;
  char ifo[NIFO_MAX][4];
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"");    // if empty then select user detectors
  strcpy(ifo[2],"");    // if empty then select user detectors
  strcpy(ifo[3],"V1");

  detectorParams detParms[4] = {
                          {"L1", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, (    -197.716)},

                          {"Y2", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, ( +45-125.999)},  // H2 LCI
                          {"Y3", 46.4551, -119.408, 0.0, 0, ( +45-125.999), 0, (    -125.999)},  // H3 LCI

                          {"V1", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},
                         };


  network NET;
  detector* pD[NIFO_MAX];
  for(int i=0; i<nIFO; i++) {
    if(strlen(ifo[i])>0) pD[i] = new detector(ifo[i]);        // built in detector
    else                 pD[i] = new detector(detParms[i]);   // user define detector
  }
  for(int i=0; i<nIFO; i++) NET.add(pD[i]);

  bool ifoBuiltin=true;
  for(int n=0; n<nIFO; n++) {
    if(pD[n]->isBuiltin()) {
     cout << pD[n]->Name << " -> Builtin " << endl;
    } else { 
     cout << pD[n]->Name << " -> User Defined " << endl;
    }
  }


  CWB::mdc MDC(&NET);
 
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

  MDC.SetInspiral("TEST",inspOptions);

  MDC.SetInjLength(10);
  int id=0;
  MDC.Draw("Y2",968650200,968650300,id,MDC_TIME);
  //MDC.Draw("Y2",968650000,968650100,id,MDC_TIME);
  //MDC.Draw("L1",968650000,968650100,id,MDC_TIME);
  //MDC.Draw("L1",968650000,968650100,id,MDC_FFT);
  //MDC.Draw("L1",968650000,968650100,id,MDC_TF);
}
