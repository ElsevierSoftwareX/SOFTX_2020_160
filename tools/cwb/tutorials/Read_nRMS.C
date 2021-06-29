// This macro read the nRMS from the cstrain root files produced with 
// the plugin CWB_Plugin_nRMS.C 
// Author : Gabriele Vedovato

//#define SAVE_PLOT	// save plots

void DumpToFile(WSeries<double>& nRMS, double gps, CWB::config* cfg, network* net, TString ifo, TString fName);
void ReadFile(TString fName); 

void Read_nRMS() {

  CWB::Toolbox TB;

  cout << "Starting reading output directory ..." << endl;
  vector<TString> fileList = TB.getFileListFromDir("output",".root","cstrain_","",true);
  int nfile = fileList.size();
  for(int n=0;n<nfile;n++) {
    ReadFile(fileList[n]); 
  }

  exit(0);
}

void ReadFile(TString fName) {

  network NET;
  CWB::config  CFG;
  detector* pD[NIFO_MAX];                       //! pointers to detectors

  // open root file
  TFile* ifile = new TFile(fName);
  if(ifile==NULL) {cout << "Error opening root file : " << fName.Data() << endl;exit(1);}
  //ifile->ls();

  // read network object from job file
  network* pnet = (network*)ifile->Get("network");
  if(pnet!=NULL) {
    // copy network object into local NET object
    NET = *pnet;
    delete pnet;
  } else {
    cout << "cwb::InitNetwork - Error : net is not contained in root file " << fName.Data() << endl;
    exit(1);
  }

  // read config object
  if(ifile->Get("config")!=NULL) {
    // read config object
    CWB::config* pcfg = (CWB::config*)ifile->Get("config");
    CFG = *pcfg;
    delete pcfg;
  } 

  int nIFO = NET.ifoListSize();

  // save detector pointers into local structures pD
  for(int n=0; n<nIFO; n++) pD[n] = NET.getifo(n);

  char cdrms_name[32];sprintf(cdrms_name,"rms-f%d",0);

  // read strain rms
  for(int n=0; n<nIFO; n++) {
    //cout << n << " " << pD[n]->Name << endl;

    ifile->cd();
    WSeries<double>* pws = (WSeries<double>*)ifile->Get(TString("rms/")+cdrms_name+"/"+pD[n]->Name);
    if(pws==NULL)
      {cout << "Error : strain rms not present, job terminated!!!" << endl;exit(1);}
    pD[n]->nRMS=*pws;;
    delete pws;
  }

  ifile->Close();

  // dump to file
  for(int n=0; n<nIFO; n++) {
    DumpToFile(pD[n]->nRMS, -1., &CFG, &NET, pD[n]->Name, fName);
  }

}

void DumpToFile(WSeries<double>& nRMS, double gps, CWB::config* cfg, network* net, TString ifo, TString fName) {

  // gps<0 : the rms is the the average over the all rms samples in the buffer
  // gps=0 : the selected rms is in the middle of the buffer
  // gps>0 : only the rms at that time is selected 

  double GPS = gps;

  bool oneside = true;

  int levels = nRMS.getLevel()+1;               // number of levels
  int slices = nRMS.size()/levels;              // number of nRMS samples
  double length = slices*cfg->whiteStride;      // nRMS len in sec

  // nRMS do not come from a TF transform
  // nRMS params must be fixed before to pass to watplot
  double rate = nRMS.rate();
  double fNinquist = rate/2.;
  nRMS.rate(nRMS.size()/length);
  WDM<double>* wdm = (WDM<double>*) nRMS.pWavelet;
  wdm->nSTS=nRMS.size();

  // Plot nRMS Scalogram
  double start = nRMS.start();
  double stop  = nRMS.start()+nRMS.size()/nRMS.rate();
  double flow  = cfg->fLow;
  double fhigh = cfg->fHigh;

  int M = nRMS.getLevel();
  const double scale = 1./nRMS.wrate();

  int igps;
  if(GPS<=0) {
    igps = (nRMS.stop()-nRMS.start())/scale/2;    // time index half buffer
    gps  = nRMS.start()+igps*scale;
  } else {
    igps = (gps-nRMS.start())/scale + 1;    	  // time index at gps time
  }

  int N = nRMS.size()/(M+1);			  // number of time indexes
/*
  // dump nRMS samples at the same frequency
  for(int j=0; j<N; j++) {
    cout << j*scale << " " << nRMS.data[j*(M+1)+100] << endl;
  }
*/

  // dump nRMS to file
  fName.ReplaceAll("cstrain_","nrms_");
  fName.ReplaceAll(".root",".txt");
  cout << endl << "Dump nRMS " << " : " << fName << endl;

  ofstream out;
  out.open(fName.Data(),ios::out);
  if(!out.good()) {cout << "Error Opening File : " << fName << endl;exit(1);}

  double inRate = cfg->fResample>0 ? cfg->fResample : cfg->inRate;
  double rescale = oneside ? sqrt(2.)/sqrt(inRate) : 1./sqrt(inRate);

  double df=fNinquist/M;
  for(int j=1; j<M; ++j) {
    double f = j*df+df/2.;
    if(f<flow || f>fhigh) continue;
    double rms = 0;
    if(GPS<0) {
      for(int i=0; i<N; i++) rms += nRMS.data[i*(M+1)+j]*rescale;
      rms/=N;
    } else {
      rms = nRMS.data[igps*(M+1)+j]*rescale;
    }
    out << j*df+df/2. << " " << rms << endl;
  }

  out.close();

#ifdef SAVE_PLOT
  char cmd[1024];
  // execute cwb_draw_sensitivity
  sprintf(cmd,"export CWB_SENSITIVITY_FILE_NAME=\"%s\";",fName.Data());
  sprintf(cmd,"%s export CWB_SENSITIVITY_SAVE_PLOT=%d;",cmd,1);
  sprintf(cmd,"%s export CWB_SENSITIVITY_RANGE_FIX=%d;",cmd,1);
  sprintf(cmd,"%s root -n -l -b ${CWB_MACROS}/cwb_draw_sensitivity.C",cmd);
  cout << cmd << endl; gSystem->Exec(cmd);
#endif

}
