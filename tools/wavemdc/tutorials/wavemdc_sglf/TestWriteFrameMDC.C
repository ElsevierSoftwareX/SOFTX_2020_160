{

  #include <vector>

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



  int nIFO=3;
  TString ifo[nIFO]={"L1","H1","V1"};
  CWB::mdc MDC(nIFO,ifo); 

  char wf_name[256];
  waveform wf;
  for (int n=0;n<NWAVEFORM;n++) {
    sprintf(wf_name,"SG%dQ%1.1f",F[n],Q[n]);
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    wf.name.ReplaceAll("d0","");
    wf.hp = MDC.GetSGQ(F[n],Q[n]);
    wf.hx = wf.hp;
    wf.hx = 0;
    MDC.AddWaveform(wf);
  }
  MDC.PrintWaveformsList();

  TArrayD parms(3);parms[0]=100000;parms[1]=1;parms[2]=100;
  MDC.SetSkyDistribution(MDC_RANDOM,"",&parms);

  //TArrayD parms(1);parms[0]=100000;  // number of entries
  //MDC.SetSkyDistribution(MDC_MNGD,"",&parms);

  //TArrayD parms(1);parms[0]=100;  // distance max = 100 Mpc
  //MDC.SetSkyDistribution(MDC_GWGC,"../data/GWGCCatalog_Rev1d2.txt",&parms);

  //TArrayD parms(5);parms[0]=1111111;parms[1]=30;parms[2]=60;parms[3]=90;parms[4]=1;
  //MDC.SetSkyDistribution(MDC_EARTH_FIX,"",&parms);

  MDC.SetInjRate(0.1);
  MDC.SetInjTimeJitter(1.0);

  TString frDir   = "frames";
  TString frLabel = "TEST";
  size_t gps      = 123456789;
  size_t length   = 1000;
  bool log        = true;

  //TString ofName = MDC.WriteFrameFile(frDir, frLabel, gps, length, log);
  for(int n=0;n<3;n++) {
    MDC.WriteFrameFile(frDir, frLabel, gps, length, log);
    gps+=length;
  }

  //cout << ofName.Data() << endl;

  exit(0);
}
