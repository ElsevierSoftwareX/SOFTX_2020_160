#define GWGCCatalog "$CWB_GWAT/data/GWGCCatalog_Rev1d8.txt"

{

  int seed = 10;

  int nIFO=3;
  TString ifo[nIFO]={"L1","H1","V1"};
  CWB::mdc MDC(nIFO,ifo);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par; 

  par.resize(5);
  par[0].name="frequency"; par[0].value=250.;
  par[1].name="bandwidth"; par[1].value=100.;
  par[2].name="duration";  par[2].value=0.1;
  for(int n=0;n<30;n++) {
    par[3].name="pseed";     par[3].value=seed+100000+n;
    par[4].name="xseed";     par[4].value=seed+100001+n;
    MDC.AddWaveform(MDC_WNB, par);
  }

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=3.;
  MDC.AddWaveform(MDC_SG, par);

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=8.9;
  MDC.AddWaveform(MDC_SG, par);

  par.resize(2);
  par[0].name="frequency"; par[0].value=235.;
  par[1].name="Q"; par[1].value=9.;
  MDC.AddWaveform(MDC_SGC, par);

  MDC.Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC.SetInjHrss(2.5e-21);
  MDC.SetInjRate(0.0333333);
  MDC.SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  par.resize(1);
  par[0].name="distance_thr"; par[0].value=50000;
  MDC.SetSkyDistribution(MDC_GWGC,gSystem->ExpandPathName(GWGCCatalog),par);

  MDC.DrawSkyDistribution("skyplot","","geographic",2,true);
  //MDC.DrawSkyDistribution("skyplot","","celestial",2,true);
  //MDC.DrawSkyDistribution("skyplot","","celestial",2,false);
}
