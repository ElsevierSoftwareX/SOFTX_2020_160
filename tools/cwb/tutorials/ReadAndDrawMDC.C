// Draw Waveform from text file
// Author : Gabriele Vedovato
//

CWB::mdc* MDC; 

void ReadAndDrawMDC() {

  MDC = new CWB::mdc; 

  TString fName = "../waveforms/ott-burrow/s25WW.h.dat";
  wavearray<double> x;
  MDC->SetInjLength(1.5);
  MDC->ReadWaveform(x, fName);
  //MDC->Draw(x);
  //MDC->Draw(x,MDC_FFT);
  MDC->Draw(x,MDC_TF);

  //exit(0);
}
