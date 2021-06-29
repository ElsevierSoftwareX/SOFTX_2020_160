//
// Test Create FrList
// Author : Gabriele Vedovato

{
  #define FRLIST_NAME  "input/H1_LDAS_C02_L2.frl"
  #define CHANNEL_NAME "H1:LDAS-STRAIN"

  double start = 942449664;
  double stop = start+600;
  wavearray<double> xx;

  CWB::frame fr(FRLIST_NAME);

  int nfiles=fr.getNfiles();
  cout << "nfiles : " << nfiles << endl;
  frfile FRF = fr.getFrList(start, stop, 0);
  fr.readFrames(FRF,const_cast<char*>(CHANNEL_NAME),xx);

  exit(0);
}
