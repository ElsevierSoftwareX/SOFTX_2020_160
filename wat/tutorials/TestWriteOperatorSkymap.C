{
  // cwb sky segmentation
  skymap sm(0.4);

  int L = sm.size();
  cout << "skymap size : " << L << endl; 

  for(int l=0;l<L;l++) sm.set(l,l); 
  // save skymap in binary format
  sm >> const_cast<char*>("skymap.dat");
  // save skymap object to root file
  sm >> const_cast<char*>("skymap.root");

  // healpix sky segmentation
  skymap esm(7);

  int eL = esm.size();
  cout << "healpix skymap size : " << eL << endl; 

  for(int l=0;l<L;l++) esm.set(l,l); 
  // save skymap in fits format
  esm >> const_cast<char*>("skymap.fits");

  exit(0);
}
