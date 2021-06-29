{
  //
  // Instantiation of object CWB::mdc and setup of waveforms, sky distribution
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
    wf.type = MDC_USER;
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

  //vector<mdcpar> par(1);
  //par[0].name="entries"; par[0].value=100000;
  //MDC.SetSkyDistribution(MDC_MNGD,par);

  //vector<mdcpar> par(1);
  //par[0].name="distance_thr"; par[0].value=100;
  //MDC.SetSkyDistribution(MDC_GWGC,"../data/GWGCCatalog_Rev1d7.txt",par);

/*
  vector<mdcpar> par(3);
  par[0].name="theta"; par[0].value=30;
  par[1].name="phi"; par[1].value=60;
  par[2].name="rho"; par[2].value=1;
  MDC.SetSkyDistribution(MDC_EARTH_FIX,par);
*/

  double theta; double phi; double psi; double rho;
  double xiota; double hrss; int ID; int id;
  for(int n=0;n<10;n++) {
    MDC.GetSourceCoordinates(theta, phi, psi, rho,xiota, hrss, ID, id);
    cout << "Source : " << theta << " " << phi << " " << psi << " " << rho
         << xiota << " " << hrss << " " << ID << " " << id << endl;
  }

  //MDC.DrawSkyDistribution("skyplot","","geographic",2,true);

  //MDC.AddWaveform("SG554Q8d9","Waveforms/SG554Q8d9.txt");

  //waveform wf = MDC.GetWaveform("SG554Q8d9",0);
  //wf.name="SG554Q8d9_CLONE";
  //MDC.AddWaveform(wf);

  //wavearray<double> wnb = MDC.GetWNB(1000., 1000., 0.01);
  //MDC.Draw(wnb);
  //MDC.Draw(wnb,MDC_FFT);
  //MDC.Draw(wnb,MDC_TF);

  //wavearray<double> rd = MDC.GetRD(1000., 0.2, 10.);
  //MDC.Draw(rd);
  //MDC.Draw(rd,MDC_TF);

  //wavearray<double> sgq = MDC.GetSGQ(554., 8.9);
  //MDC.Draw(sgq);
  //MDC.Draw(sgq,MDC_TF);

  //cout << MDC.wfList[0].hp.size() << endl;

  MDC.wfList[14].hx=MDC.wfList[14].hp;
  MDC.TimeShift(MDC.wfList[14].hx,0.001);
  //MDC.PhaseShift(MDC.wfList[0].hx,90);

  //MDC.Print();

  //MDC.Draw(14,0,"hp",MDC_TIME);
  //MDC.Draw(14,0,"hx",MDC_TIME,"same",kRed);

  // resampling waveform
  int resampling_factor = 8;
  wavearray<double> hp = MDC.wfList[14].hp;
  wavearray<double> hx = MDC.wfList[14].hx;
  int size = hp.size();
  int rate = hp.rate();
  hp.FFT(1);
  hp.resize(resampling_factor*size);
  hp.rate(resampling_factor*rate);
  for(int i=size;i<resampling_factor*size;i++) hp[i]=0;
  hp.FFT(-1);
  hx.FFT(1);
  hx.resize(resampling_factor*size);
  hx.rate(resampling_factor*rate);
  for(int i=size;i<resampling_factor*size;i++) hx[i]=0;
  hx.FFT(-1);
  MDC.Draw(hp,MDC_TIME,"zoom");
  MDC.Draw(hx,MDC_TIME,"same",kRed);

  //MDC.Draw(0,0,"hp",MDC_FFT);
  //MDC.Draw(0,0,"hx",MDC_FFT,"same",kRed);

  //MDC.Draw(0,0,"hp",MDC_TF);

/*

  MDC.SetInjRate(0.1);
  MDC.SetInjTimeJitter(1.0);

  double srate=16384,;
  wavearray<double> x(25*srate);
  x.rate(srate);
  x.start(10000000.);
  MDC.Get(x,"L1");
  //for(int i=0;i<x.size();i++) x[i]=i; 
  //MDC.Draw(x);
  //MDC.Draw(x,MDC_FFT);
  //MDC.Draw(x,MDC_TF);
*/

  //exit(0);
}
