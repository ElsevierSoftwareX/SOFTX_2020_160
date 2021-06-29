//!NOISE_MDC_SIMULATION
// Config Plugin to injected 'on the fly' MDC & to set SkyMaskCC 'on the fly'

{
  #define HEN_LIST "config/HEN_list_2009_2010.txt"
  #define ONSOURCE_TIME   1000   // sec

  cout << "-----> CWB_Plugin_HEN_SIM_Config.C" << endl;

  network** net;
  CWB_PLUGIN_IMPORT(network**,net);

  int* xstart;
  int* xstop;
  CWB_PLUGIN_IMPORT(int*,xstart);
  CWB_PLUGIN_IMPORT(int*,xstop);


  int seed = (*net)->nRun;  // WARNING : seed must be the same for each detector within the same job

  double GPS=0;
  float RA;
  float DEC;
  float RADIUS;

  ifstream in;
  in.open(HEN_LIST,ios::in);
  if (!in.good()) {
    cout << "CWB_Plugin_HEN_SIM_Config.C - Error Opening File : " << HEN_LIST << endl;
    gSystem->Exit(1);
  }

  char str[1024];
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] == '#') continue;
    TObjArray* token = TString(str).Tokenize(TString(" "));
    if(token->GetEntries()!=16) {
      cout << "CWB_Plugin_HEN_SIM_Config.C - bad line format : " << str << endl;
      gSystem->Exit(1);
    }
    TObjString* otoken;
    TString stoken;

    otoken = (TObjString*)token->At(2);
    stoken = otoken->GetString();
    GPS = stoken.Atof();

    otoken = (TObjString*)token->At(3);
    stoken = otoken->GetString();
    RA = stoken.Atof();

    otoken = (TObjString*)token->At(4);
    stoken = otoken->GetString();
    DEC = stoken.Atof();

    otoken = (TObjString*)token->At(5);
    stoken = otoken->GetString();
    RADIUS = stoken.Atof();

    int GPS1 = GPS-ONSOURCE_TIME/2;
    if(GPS1>xstart && GPS1<=xstop) break;
    
    int GPS2 = GPS+ONSOURCE_TIME/2;
    if(GPS2>xstart && GPS2<=xstop) break;
    
    delete token;
  }
  in.close();

  if(GPS==0) {
    cout << "CWB_Plugin_HEN_SIM_Config.C - no trigger in the range : " << xtart << " " << xstop << endl;
    gSystem->Exit(1);
  } else {
    cout << "CWB_Plugin_HEN_SIM_Config.C - trigger found at gps : " 
         << GPS << " ra : " << RA << " dec : " << DEC << " radius : " << RADIUS << endl;
  }

  // --------------------------------------------------------
  // define MDC injections
  // --------------------------------------------------------

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  char wf_name[256];
  waveform wf;
  vector<mdcpar> par;

  par.resize(2);

  par[0].name="frequency"; par[0].value=100.;
  par[1].name="Q"; par[1].value=8.9;
  MDC->AddWaveform(MDC_SGC, par);

  par[0].name="frequency"; par[0].value=153.;
  par[1].name="Q"; par[1].value=8.9;
  MDC->AddWaveform(MDC_SGC, par);

  par[0].name="frequency"; par[0].value=1053.;
  par[1].name="Q"; par[1].value=9;
  MDC->AddWaveform(MDC_SGC, par);

  par[0].name="frequency"; par[0].value=1304.;
  par[1].name="Q"; par[1].value=9;
  MDC->AddWaveform(MDC_SGC, par);

  MDC->Print();

  // --------------------------------------------------------
  // define injection parameters
  // --------------------------------------------------------
  MDC->SetInjHrss(2.5e-21);
  MDC->SetInjRate(0.02);
  MDC->SetInjJitter(10.0);

  // --------------------------------------------------------
  // define sky distribution
  // --------------------------------------------------------
  par.resize(3);
  par[0].name="theta";   par[0].value=DEC;
  par[1].name="phi";     par[1].value=RA;
  par[2].name="entries"; par[2].value=1000;
  MDC->SetSkyDistribution(MDC_CELESTIAL_FIX,par,seed);


  // --------------------------------------------------------
  // define SetSkyMaskCC
  // --------------------------------------------------------
  char options[256];
  sprintf(options,"--dec %f --ra %f --radius %f",DEC,RA,RADIUS);
  cwb* CWB = new cwb;
  CWB->SetSkyMaskCC(net,cfg,options);
  delete CWB;
}
