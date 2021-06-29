#define MSMBH 3.5e6
#define MMIN  10
#define MMAX  25

#define FREQ_CUTOFF 20  // Hz

#define NEVENTS 1000

//#define BUILTIN_WRITE

{
  double mSMBH = MSMBH;
  double mmin  = MMIN;
  double mmax  = MMAX;
  double beta  = 2;

  GNGen eGen(mSMBH, mmin, mmax, beta);
  eGen.setFreqCutoff(FREQ_CUTOFF);

  char fName[256];
  sprintf(fName,"macro/eBBH_%g_%g_%g_%g.lst",mSMBH,mmin,mmax,beta);
  TString sName = fName;
  sName.ReplaceAll("+","");		// remove '+' character
  sprintf(fName,"%s",sName.Data());

#ifdef BUILTIN_WRITE
  cout << "Starting eBBH parameters generation ..." << endl;
  eGen.generateEvents(NEVENTS,fName);
#else
  double m1, m2, rp, e;
  FILE* f = stdout;
  f=fopen(fName, "w");
  int nevents=0;
  while(nevents<NEVENTS) {

    eGen.generateEvent(m1, m2, rp, e);

    if(e<0.5) {
      //printf("%d %lf %lf %lf %.8lf\n", nevents, m1, m2, rp, e);
      fprintf(f, "%d %lf %lf %lf %.8lf\n", nevents++, m1, m2, rp, e);
    }
  }
  fclose(f);
#endif

  gSystem->Exit(0);
}
