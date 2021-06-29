{
  //
  // Example which shown how to rescale the strain of eBBH according to the distance of the source
  // Author : Gabriele Vedovato


  // the distance of source is G*M/c^2  meters
  double G  = watconstants::GravitationalConstant();
  double M  = watconstants::SolarMass();
  double c  = watconstants::SpeedOfLightInVacuo();
  double pc = watconstants::Parsec();

  // source parameters
  double m1  = 19.8523;
  double m2  = 23.8574;
  double rp0 = 13.955599;
  double e0  = 0.3749609;
  double distance = 4509.8798;

  wavearray<double> hp, hx;
  // get hp,hx 
  getEBBH(m1, m2, rp0, e0, hp, hx);

  // distance of the source computed with getEBBH
  double distance_source_Kpc = (m1+m2)*G*M/(c*c)/pc/1.e3;  	// = 2.09169e-15 Kpc
  cout << "distance source Kpc : " << distance_source_Kpc << " (Kpc) " << endl; 

  // sample time
  double dt=1./hp.rate();

  // compute hp hrss @ distance=distance_source_Kpc
  double hrss_hp=0;
  for(int i=0;i<hp.size();i++) hrss_hp+=hp[i]*hp[i];
  hrss_hp*=dt; hrss_hp=sqrt(hrss_hp);

  // compute hx hrss @ distance=distance_source_Kpc
  double hrss_hx=0;
  for(int i=0;i<hx.size();i++) hrss_hx+=hx[i]*hx[i];
  hrss_hx*=dt; hrss_hx=sqrt(hrss_hx);

  // rescale hp,hx to 10 Kpc
  hrss_hp*=distance_source_Kpc/10.;
  hrss_hx*=distance_source_Kpc/10.;

  cout << "hrss_hp @ 10Kpc : " << hrss_hp << endl;		// = 1.52606e-17
  cout << "hrss_hx @ 10Kpc : " << hrss_hx << endl;		// = 1.52422e-17

  // strain=sqrt(sqrt(hrss_hp*hrss_hp+hrss_hx*hrss_hx) @ 10Kpc
  double strain_10kpc=sqrt(hrss_hp*hrss_hp+hrss_hx*hrss_hx);
  cout << "strain @ 10Kpc : " << strain_10kpc << endl;		// = 2.15687e-17

  // rescale strain @ 4509.8798 Kpc
  double strain = 10*strain_10kpc/distance;
  cout << "strain @ 4509.8798 Kpc : " << strain << endl;	// = 4.78254e-20

  exit(0);

}
