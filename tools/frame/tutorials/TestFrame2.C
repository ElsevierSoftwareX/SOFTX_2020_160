//
// Test Create FrList
// Author : Gabriele Vedovato

#define FRLIST_NAME  "input/H1_LDAS_C02_L2.frl"
#define CHANNEL_NAME "H1:LDAS-STRAIN"

using namespace CWB;

void TestFrame2() {

  wavearray<double> x;
  x.start(942449664);
  x.stop(x.start()+600);

  frame ifr(FRLIST_NAME,CHANNEL_NAME,"READ");
  int nfiles=ifr.getNfiles();
  cout << "nfiles : " << nfiles << endl;

  ifr >> x;

  cout << x.start() << " " << x.rate() << " " << x.size() << endl;

  frame ofr("TEST.gwf","CHNAME","WRITE");
//  ofr.setFrName("FRNAME");
  ofr << x;
  ofr.close();

  exit(0);
}
