//#define READ_SETUP_FROM_XML_FILE

{
  //
  // Show how to use mdc class to get & draw LAL inpiral waveforms
  // Author : Gabriele Vedovato

  #include <vector>

  CWB::mdc MDC; 

#ifndef READ_SETUP_FROM_XML_FILE
  // ---------------------------------
  // set inspiral parms
  // ---------------------------------
  TString inspOptions="";
  inspOptions = "--time-step 60.0 --time-interval 3 ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--gps-start-time 931072130 --gps-end-time 933491330 ";
  inspOptions+= "--d-distr volume --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 10.000000 ";
//  inspOptions+= "--min-mass1 25.000000 --max-mass1 225.000000 ";
//  inspOptions+= "--min-mass2 25.000000 --max-mass2 225.000000 ";
  inspOptions+= "--min-mtotal 50.000000 --max-mtotal 250.000000 ";
  inspOptions+= "--min-mratio 0.25 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 1000000.0 --max-distance 1500000.0 ";
  inspOptions+= "--approximant EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789 ";
  inspOptions+= "--dir ./ ";
  inspOptions+= "--output inspirals.xml ";		// set output xml file

  // set and write xml file
  MDC.SetInspiral("EOBNRv2",inspOptions);

#else 
  // set inspiral using xml file (speedup GetInspiral method)
  TString inspOptions="--xml inspirals.xml ";
  inspOptions+= "--dir ./ ";
  MDC.SetInspiral("EOBNRv2",inspOptions);
#endif

  // Get the first waveform hp,hx components starting from gps = 931072130
  wavearray<double> hp = MDC.GetInspiral("hp",931072130,931072230);
  wavearray<double> hx = MDC.GetInspiral("hx",931072130,931072230);
  cout << "size : " << hp.size() << " rate : " << hp.rate() << " start : " << (int)hp.start() << endl;
  hp.start(0);		// set start to 0 (needed by draw Method)
  hx.start(0);
  MDC.Draw(hp,MDC_TIME);
  MDC.Draw(hx,MDC_TIME,"same",kRed);

  //MDC.Draw(hp,MDC_FFT);	// draw hp in frequency domain
  //MDC.Draw(hp,MDC_TF);        // draw hp in time-frequency domain;

}
