
#define GPS_START_TIME	1116917000
#define GPS_END_TIME	1116918000

#define TMP_BUF_LENGTH  100	// temporary buffer length in sec used for MDC generation

#define XML_FILE_NAME	"burst.xml"
#define LOG_FILE_NAME	"burst.txt"	// create BurstMDC log file

void MakeBurstXML() {
  //
  // Show how to use mdc class to create a LAL Burst XML 
  // Author : Gabriele Vedovato

  // MDC class 
  #define N_IFO 3
  TString ifo[N_IFO]={"L1","H1","V1"};
  CWB::mdc MDC(N_IFO,ifo);

  // add SGE waveforms 
  // uniform distribution in frequency [40:1000] Hz
  // uniform distribution in Q [3:30] Hz
  vector<mdcpar> par(3);

  double value;
  int    NWAVE 	 = 100;
  float  minfreq = 40.;
  float  maxfreq = 1500.;
  float  minq    = 3;
  float  maxq    = 30.;
  int    seed    = 12345;

  gRandom->SetSeed(seed);
  for(int n=0;n<NWAVE;n++) {
    value = gRandom->Uniform(minfreq,maxfreq);
    par[0].name="frequency"; par[0].value=value;
    value = gRandom->Uniform(minq,maxq);
    par[1].name="Q"; par[1].value=value;
    par[2].name="decimals"; par[2].value=3;
    MDC.AddWaveform(MDC_SGE, par);
  }
  MDC.Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------

  MDC.SetInjHrss(2.5e-21);    // injection hrss (the factors array in user_parameter.C file is used to rescale it)
  MDC.SetInjRate(0.1);        // injection rate (injections per sec)
  MDC.SetInjJitter(1.0);      // time jitter in sec

  // ---------------------------------------------------------------------------------------
  // define sky distribution
  // this example show how to set a random sky distribution
  // for a full description of the parameters see the CWB::mdc::SetSkyDistribution in mdc.cc
  // ---------------------------------------------------------------------------------------

  par.resize(4);
  par[0].name="entries"; par[0].value=100000;  // number of entries randomly generated
  par[1].name="rho_min"; par[1].value=0;       // the source hrss is not rescaled
  par[2].name="rho_max"; par[2].value=0;       // the source hrss is not rescaled
  par[3].name="iota";    par[3].value=-1;      // iota : random distribution, iota < 0: random distribution
  MDC.SetSkyDistribution(MDC_RANDOM,par,seed); // random sky distribution (theta,phi,psi)

  par.resize(10);
  par[0].name="gps-start-time"; par[0].value=GPS_START_TIME; 
  par[1].name="gps-end-time";   par[1].value=GPS_END_TIME; 
  par[2].name="f-distr";  	par[2].svalue="uniform"; 
  par[3].name="min-frequency";  par[3].value=minfreq; 
  par[4].name="max-frequency";  par[4].value=maxfreq; 
  par[5].name="q-distr";  	par[5].svalue="uniform"; 
  par[6].name="min-q";  	par[6].value=minq; 
  par[7].name="max-q";  	par[7].value=maxq; 
  par[8].name="seed";  		par[8].value=seed; 
  par[9].name="population";  	par[9].svalue="all_sky_sinegaussian"; 
  MDC.CreateBurstXML(XML_FILE_NAME,par);
#ifdef LOG_FILE_NAME	
  MDC.DumpLogHeader(LOG_FILE_NAME,"",1);
#endif
  wavearray<double> x(MDC.GetSampleRate()*TMP_BUF_LENGTH);
  x.rate(MDC.GetSampleRate());
  bool verbose = true;
  double start=GPS_START_TIME;
  int loops = TMath::CeilNint((GPS_END_TIME-GPS_START_TIME)/TMP_BUF_LENGTH);
  for(int i=0;i<loops;i++) {
    if(i%10==0) cout << "loop : " << i << endl;
    start += x.size()/x.rate(); 
    x.start(start);
    MDC.GetBurst(x,"H1");
    MDC.FillBurstXML(verbose);
#ifdef LOG_FILE_NAME	
    MDC.DumpLog(LOG_FILE_NAME,"",true);
#endif
  }
  MDC.CloseBurstXML();

  gSystem->Exit(0);
}
