//
// Test Create FrList
// Author : Gabriele Vedovato

{
  #define FRLIST_NAME  "input/H1_LDAS_C02_L2.frl"
  #define CHANNEL_NAME "H1:LDAS-STRAIN"

  double start = 942449664;
  double stop = start+600;
  wavearray<double> xx;
  wavearray<double> yy;

  double sRate = 8192./2;  // if > 0 && > (original srate) data are resampled to sRate

  CWB::frame fr(FRLIST_NAME);

  int nfiles=fr.getNfiles();
  cout << "nfiles : " << nfiles << endl;
  frfile FRF = fr.getFrList(start, stop, 0);
  fr.readFrames(FRF,const_cast<char*>(CHANNEL_NAME),xx);
  fr.setSRIndex(12);
  fr.readFrames(FRF,const_cast<char*>(CHANNEL_NAME),yy);

  cout << "srate " << yy.rate() << endl;

  watplot plot(const_cast<char*>("plot"),200,20,800,500);

  char gtitle[256];  
  sprintf(gtitle,"Original PSD (Black) Reasampled PSD (Red)");
  plot.gtitle(gtitle,"frequency (Hz)","strain/#sqrt{Hz}");
  plot.goptions(const_cast<char*>("alp logy"), 1, 0, 0, true, 0,xx.rate()/2., true, 32);

  xx >> plot; yy >> plot; 
/*
  char ofName[256];  
  sprintf(ofName,"ResamplePSD.png");
  cout << "write results to " << ofName << endl;
  TString gfile=ofName;
  plot >> gfile;

  exit(0);
*/
}
