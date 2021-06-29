//
// Test Create FrList
// Author : Gabriele Vedovato

#define FRLIST_NAME  "input/H1_LDAS_C02_L2.frl"
#define CHANNEL_NAME "H1:LDAS-STRAIN"

using namespace CWB;

void TestFrame3() {

  wavearray<double> x;
  x.start(942449664);
  x.stop(x.start()+600);

  frame fr(FRLIST_NAME,CHANNEL_NAME,"READ");
  int nfiles=fr.getNfiles();
  cout << "nfiles : " << nfiles << endl;
  fr >> x;
  fr.close();

  fr.open("TEST.gwf",CHANNEL_NAME,"WRITE");
  fr << x;
  fr.close();

  exit(0);
}
