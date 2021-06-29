//
// Test Create FrList
// Author : Gabriele Vedovato

#define FRLIST_NAME  "input/H1_LDAS_C02_L2.frl"
#define CHANNEL_NAME "H1:LDAS-STRAIN"

using namespace CWB;

watplot* plot;

void TestFrame4() {

  wavearray<double> x;
  x.start(942449664);
  x.stop(x.start()+600);

  plot = new watplot(const_cast<char*>("plot"),200,20,800,500);
  //plot->goptions(const_cast<char*>("alp logx logy"), 1, 0., 0., true, 0., 0.);
  //plot->goptions(const_cast<char*>("alp logx logy"), 1, x.start(), x.start()+100., true, 10., 0.);
  plot->goptions(const_cast<char*>("alp logx logy"), 1, x.start(), x.start()+100., true, 10., 0., true, 8.);
  //plot->goptions(const_cast<char*>("alp logy"), 1, x.start(), x.stop(), true, 170., 190., true, 64.);
  plot->gtitle(CHANNEL_NAME);

  int mkf_argc = 11;
  const char *mkf_argv[11] = {"mkf", "-Bu","-Hp", "-o", "6", "-a", "200", "-s", "16384","-l",NULL};
  Filter filter(mkf_argc, mkf_argv);

  frame ifr(FRLIST_NAME,CHANNEL_NAME);
  int nfiles=ifr.getNfiles();
  cout << "nfiles : " << nfiles << endl;

  char frName[512];
  sprintf(frName,"%s-%s-%d-%d.gwf","H1","Filtered",int(x.start()),int(x.stop()-x.start()));
  cout << frName << endl;

  frame ofr(frName,CHANNEL_NAME,"WRITE");
  ofr.setFrName("H1:GW-H");

  //ifr >> x >> filter >> x >> *plot;
  //ifr >> x >> filter >> x >> ofr;
  //ifr >> x >> *plot >> x >> filter >> x >> ofr;
  //ifr >> x >> filter >> x >> *plot >> x >> ofr;
  //ifr >> x >> *plot >> x >> filter >> x >> *plot >> x >> ofr;

  char plName[512];
  sprintf(plName,"%s-%s-%d-%d.png","H1","Filtered",int(x.start()),int(x.stop()-x.start()));
  cout << plName << endl;
  //plot->print(plName);
  *plot >> plName;

  ifr.close();
  ofr.close();

  //exit(0);
}
