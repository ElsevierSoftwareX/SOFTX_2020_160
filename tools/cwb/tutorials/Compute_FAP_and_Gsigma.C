// Compute Gaussian Sigma & FAP from OBSERVATIONAL_TIME, BACKGROUND_TIME, TRIALS_FACTOR

{
  #define OBSERVATIONAL_TIME	16.		// days
  #define BACKGROUND_TIME	67800.		// years
  #define TRIALS_FACTOR		3		// number of bins

  double N = OBSERVATIONAL_TIME/365.*(1./BACKGROUND_TIME);
  double FAP = 1-exp(-N*TRIALS_FACTOR);
  double Gsigma = sqrt(2)*TMath::ErfcInverse(2.*FAP);

  double FAP_from_Gsigma = TMath::Erfc(Gsigma*1./sqrt(2))/2;	// xcheck

  cout << endl;
  cout << "-----------------------------------------------" << endl;
  cout << "OBSERVATIONAL_TIME : " << OBSERVATIONAL_TIME << " days" << endl;
  cout << "BACKGROUND_TIME    : " << BACKGROUND_TIME << " years" << endl;
  cout << "TRIALS_FACTOR      : " << TRIALS_FACTOR << endl;
  cout << "-----------------------------------------------" << endl;
  cout << "FAP                : " << FAP << endl;
  cout << "Gaussian sigma     : " << Gsigma << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;

  exit(0);
}
