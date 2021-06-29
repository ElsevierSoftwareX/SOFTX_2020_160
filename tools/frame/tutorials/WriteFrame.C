//
// Write Frame File
// Author : Gabriele Vedovato


{
  char frFile[256];
  char frName[256];
  char chName[256];

  wavearray<double> x;
  x.rate(16386);
  x.start(123456789);
  x.resize(10*x.rate());
  x=1.;

  int gps = int(x.start());
  int length = int(x.size()/x.rate());

  sprintf(frFile,"frfile-%d-%d.gwf",gps,length);
  sprintf(frName,"FRNAME");
  sprintf(chName,"CHNAME");

  CWB::frame fr(frFile,chName,"WRITE");
  fr.writeFrame(x, frName, chName);

  exit(0);
}
