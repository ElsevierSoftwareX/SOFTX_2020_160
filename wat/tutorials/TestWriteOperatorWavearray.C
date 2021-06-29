{
  wavearray<double> x(100);
  for(int i=0;i<x.size();i++) x[i]=i;
  // save wavearray in binary format
  x >> const_cast<char*>("wavearray.dat");
  // save wavearray in ascii format
  x >> const_cast<char*>("wavearray.txt");
  // save object wavearray in root format
  x >> const_cast<char*>("wavearray.root");
  exit(0);
}
