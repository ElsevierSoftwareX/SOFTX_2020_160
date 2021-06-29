{
  //
  // Create waveforms and write to frame file
  // Author : Gabriele Vedovato

  #define N_IFO 3
  #define NWAVEFORM 15

  double F[NWAVEFORM];
  double Q[NWAVEFORM];

  F[0]=70;      Q[0]=8.9;
  F[1]=100;     Q[1]=8.9;
  F[2]=153;     Q[2]=8.9;
  F[3]=235;     Q[3]=8.9;
  F[4]=361;     Q[4]=8.9;
  F[5]=554;     Q[5]=8.9;
  F[6]=849;     Q[6]=8.9;
  F[7]=945;     Q[7]=9;
  F[8]=1053;    Q[8]=9;
  F[9]=1172;    Q[9]=9;
  F[10]=1304;   Q[10]=9;
  F[11]=1451;   Q[11]=9;
  F[12]=1615;   Q[12]=9;
  F[13]=1797;   Q[13]=9;
  F[14]=2000;   Q[14]=9;



  TString ifo[N_IFO]={"L1","H1","V1"};
  CWB::mdc MDC(N_IFO,ifo); 

  char wf_name[256];
  waveform wf;
  for (int n=0;n<NWAVEFORM;n++) {
    sprintf(wf_name,"SG%fQ%1.1f",F[n],Q[n]);
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    wf.name.ReplaceAll("d0","");
    wf.hp = MDC.GetSGQ(F[n],Q[n]);
    wf.hx = wf.hp;
    wf.hx = 0;
    MDC.AddWaveform(wf);
  }
  MDC.Print();


  vector<mdcpar> par(3);
  par[0].name="entries"; par[0].value=100000;
  par[1].name="rho_min"; par[1].value=5;
  par[2].name="rho_max"; par[2].value=10;
  MDC.SetSkyDistribution(MDC_RANDOM,par);

  //par.resize(1);
  //par[0].name="entries"; par[0].value=100000;
  //MDC.SetSkyDistribution(MDC_MNGD,par);

  //vector<mdcpar> par(1);
  //par[0].name="distance_thr"; par[0].value=100;
  //MDC.SetSkyDistribution(MDC_GWGC,"../data/GWGCCatalog_Rev1d7.txt",par);

/*
  par.resize(3);
  par[0].name="theta"; par[0].value=30;
  par[1].name="phi"; par[0].value=60;
  par[2].name="rho"; par[0].value=1;
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par);
*/

  MDC.SetInjRate(0.1);
  MDC.SetInjJitter(1.0);

  TString frDir   = "frames";
  TString frLabel = "TEST";
  size_t gps      = 123456789;
  size_t length   = 1000;
  bool xlog       = true;

  TString ofName = MDC.WriteFrameFile(frDir, frLabel, gps, length, xlog);
/*
  for(int n=0;n<3;n++) {
    MDC.WriteFrameFile(frDir, frLabel, gps, length, xlog);
    gps+=length;
  }
*/
  //cout << ofName.Data() << endl;

  exit(0);
}
