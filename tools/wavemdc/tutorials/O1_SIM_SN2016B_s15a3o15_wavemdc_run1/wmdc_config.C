{
  // --------------------------------------------------------
  // define job parameters
  // --------------------------------------------------------
  TString frDir       = "frames";
  TString frLabel     = "SN2016B_s15a3o15";

  // define frame segment list
  TString segmentList = "Segments/SN2016B_s15a3o15_segments.lst";

  // select segments
  int jobmin  = 1;     // start segment
  int jobmax  = 100;   // end segment
  int jobstep = 10;    // frames per job

  // --------------------------------------------------------
  // define network
  // --------------------------------------------------------
  int nIFO=2;
  TString ifo[nIFO]={"L1","H1"};
  CWB::mdc MDC(nIFO,ifo);

  // ------------------------------------------------------------------------------
  // Ott2013(?): rotbar
  // https://wiki.ligo.org/viewauth/Bursts/O1Waveforms
  // ------------------------------------------------------------------------------ 

  #define SN1  "s15a2o09"
  //#define NTH   18     
  //#define NPH   27     

  // --------------------------------------------------------
  // define channel names
  // --------------------------------------------------------
  vector<TString> chName(nIFO);
  for(int n=0;n<nIFO;n++) chName[n] = TString::Format("%s:%s",ifo[n].Data(),SN1);

  // --------------------------------------------------------
  // define waveforms
  // --------------------------------------------------------

  int seed = 0;

  // read SN waveform file list
  TString hp_name = "Waveforms/processed_signal_s15a3o15_ls-plus.txt";
  TString hc_name = "Waveforms/processed_signal_s15a3o15_ls-cross.txt";
  vector<mdcpar> sn_par(2);                                                                        
  sn_par[0].name="name";  sn_par[0].svalue=SN1;                                                    
  sn_par[1].name="hrss";  sn_par[1].value=-1;     // do not normalize hrss                         
  MDC->AddWaveform(SN1,hp_name,hc_name,16384,sn_par);                                              

  //MDC->Print(1);

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  // if inj_hrss=0 the hrss @10Kpc is the one which is defined by hp,hc waveforms
  MDC->SetInjHrss(0);
  MDC->SetInjRate(0.0333333);
  MDC->SetInjJitter(10.0);

  // --------------------------------------------------------
  // define MNGD sky distribution
  // --------------------------------------------------------
  #define SN2016B_DEC    1.719
  #define SN2016B_RA     178.7677

  vector<mdcpar> par(5);
  par[0].name="entries"; par[0].value=1000;
  //par[1].name="theta";   par[1].value=SN2016B_DEC;
  //par[2].name="phi";     par[2].value=SN2016B_RA;
  par[1].name="psi";   par[1].value=0;
  par[2].name="iota";     par[2].value=0;
  par[3].name="rho_min";     par[3].value=0.01;
  par[4].name="rho_max";     par[4].value=10.0;
  //par[3].name="rho";     par[3].value=10.0;//kpc
  MDC->SetSkyDistribution(MDC_RANDOM,"0.13127*3*x*x - 2.4894*2*x + 20.56",par,seed);
  //MDC->SetSkyDistribution(MDC_CELESTIAL_FIX,"0.13127*3*x*x - 2.4894*2*x + 20.56",par,seed);
  //MDC->SetSkyDistribution(MDC_RANDOM,"x+100./x",par,seed);
  //MDC->SetSkyDistribution(MDC_RANDOM,par,seed);
}

