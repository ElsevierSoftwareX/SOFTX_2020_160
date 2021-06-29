{

  printf("Loading WaveMDC macros...\n");

  TString home_wavemdc;
  if(gSystem->Getenv("HOME_WAVEMDC")==NULL) {
    cout << "Error : environment HOME_WAVEMDC is not defined!!!" << endl;exit(1);
  } else {
    home_wavemdc=TString(gSystem->Getenv("HOME_WAVEMDC"));
  }

  
  gROOT->LoadMacro(home_wavemdc+"/CreateSubWaveMDC.C");
  gROOT->LoadMacro(home_wavemdc+"/CreateDagWaveMDC.C");

  // condor log
  char condor_log[512];
  UserGroup_t* uinfo = gSystem->GetUserInfo();
  sprintf(condor_log,"/local/user/%s",uinfo->fUser.Data());

  cout << "CreateSubWaveMDC ..." << endl;
  TString condor_dir   = gSystem->WorkingDirectory();
  TString condor_out   = condor_dir+"/log";
  TString condor_err   = condor_dir+"/log";
  CreateSubWaveMDC(frLabel,condor_dir,condor_out,condor_err,condor_log);

  cout << "CreateDagWaveMDC ..." << endl;
  CreateDagWaveMDC(condor_dir, frLabel, jobmin, jobmax, jobstep);

  // create wavemdc stuff 
  char cmd[256];
  sprintf(cmd,"ln -sf %s/WaveMDC.sh",home_wavemdc.Data());
  gSystem->Exec(cmd);
  gSystem->Exec("mkdir -p log frames");

  exit(0);
}
