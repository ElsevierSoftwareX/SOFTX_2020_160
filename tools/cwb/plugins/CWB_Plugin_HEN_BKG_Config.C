//!NOISE_MDC_SIMULATION
// Config Plugin to set SkyMaskCC 'on the fly'

{
  #define HEN_LIST "config/HEN_list_2009_2010.txt"
  #define ONSOURCE_TIME   1000   // sec

  cout << "-----> CWB_Plugin_HEN_BKG_Config.C" << endl;

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
    cout << "CWB_Plugin_HEN_BKG_Config.C - Error Opening File : " << HEN_LIST << endl;
    gSystem->Exit(1);
  }

  char str[1024];
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] == '#') continue;
    TObjArray* token = TString(str).Tokenize(TString(" "));
    if(token->GetEntries()!=16) {
      cout << "CWB_Plugin_HEN_BKG_Config.C - bad line format : " << str << endl;
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
    cout << "CWB_Plugin_HEN_BKG_Config.C - no trigger in the range : " << xtart << " " << xstop << endl;
    gSystem->Exit(1);
  } else {
    cout << "CWB_Plugin_HEN_BKG_Config.C - trigger found at gps : " 
         << GPS << " ra : " << RA << " dec : " << DEC << " radius : " << RADIUS << endl;
  }

  // --------------------------------------------------------
  // define SetSkyMaskCC
  // --------------------------------------------------------
  sprintf(cfg->skyMaskCCFile,"--theta %f --phi %f --radius %f",DEC,RA,RADIUS);
  cwb* CWB = new cwb;
  CWB->SetSkyMask(net,cfg,cfg->skyMaskCCFile,'c');
  delete CWB;
}
