{
  //
  // Write & Read config object to/from root file
  // Author : Gabriele Vedovato

  CWB::config config;
  config.Import("$CWB_PARAMETERS_FILE");

  char*    xlagFile = new char[1024]; strcpy(xlagFile,"LAG");
  size_t*  xlagSite = new size_t[2];  xlagSite[0]=0;xlagSite[1]=1;
  char*   xslagFile = new char[1024]; strcpy(xslagFile,"SLAG");
  size_t* xslagSite = new size_t[2];  xslagSite[0]=2;xslagSite[1]=3;

  config.lagFile  = xlagFile;
  config.lagSite  = xlagSite;
  config.slagFile = xslagFile;
  config.slagSite = xslagSite;

  config.Import();
  if(config.lagFile) cout << "config.lagFile : " << config.lagFile << endl;
  if(config.lagSite) cout << "config.lagSite : " << config.lagSite[0] << " " << config.lagSite[1] << endl;
  if(config.slagFile) cout << "config.slagFile : " << config.slagFile << endl;
  if(config.slagSite) cout << "config.slagSite : " << config.slagSite[0] << " " << config.slagSite[1] << endl;

  TFile *froot = new TFile("test.root", "RECREATE");
  config.Write("CC");
  froot->Close();

  TFile *f = new TFile("test.root");

  f->ls();

  CWB::config *iconfig = (CWB::config*)f->Get("CC");
  if(iconfig->lagFile) cout << "iconfig.lagFile : " << iconfig->lagFile << endl;
  if(iconfig->lagSite) cout << "iconfig.lagSite : " << iconfig->lagSite[0] << " " << iconfig->lagSite[1] << endl;
  if(iconfig->slagFile) cout << "iconfig.slagFile : " << iconfig->slagFile << endl;
  if(iconfig->slagSite) cout << "iconfig.slagSite : " << iconfig->slagSite[0] << " " << iconfig->slagSite[1] << endl;

  f->Close();

  exit(0);
}
