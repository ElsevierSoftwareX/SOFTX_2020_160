#include <vector>

#define NET_FILE_NAME "data/net_931158208_ADV_SIM_NSNS_L1H1V1_2G_2PHS_64_2048_1_564770044.root"
#define MDC_FILE_NAME "NSNS/MDC_L1_931158200_208_ADV_SIM_NSNS_L1H1V1_2G_2PHS_job1.root"

#define NIFO 4
#define RATE 4096.
#define LAG 0
#define CID 2

#define WDM_SPECTROGRAM
#define WDM_TF
#define WDM_T
#define WDM_F

#define READ_MDC

void ReadConfig(TString fName, network* NET, CWB::config* cfg);
void PrintConfig(network* NET, CWB::config* cfg);
void DrawSignal(int ifo, vector<int> level);
void DrawMonster(int ifo, vector<int> leveli, int iRATE=0, bool like=true);
void ReadMDC(TString fName, wavearray<double>& x);

wavearray<float>  fp[NIFO];                    // + antenna pattern
wavearray<float>  fx[NIFO];                    // x antenna pattern
wavearray<float>  am[NIFO];                    // 00 phase response
wavearray<float>  AM[NIFO];                    // 90 phase response
wavearray<float>   u[NIFO];                    // unity vector
wavearray<float>   v[NIFO];                    // unity vector
wavearray<double> nv[NIFO];                    // noise variance
wavearray<int>      layers;                    // layers
wavearray<int>        time;                    // time index
wavearray<int>   frequency;                    // frequency index
wavearray<int>  likelihood;                    // likelihood 
wavearray<int>        null;                    // null 

void MonsterEventDisplay() {

  CWB::config cfg;
  network NET;

  ReadConfig(NET_FILE_NAME, &NET, &cfg);
  PrintConfig(&NET, &cfg);

  int nIFO = NET.ifoListSize();

  float _fp[4*NIFO],_fx[4*NIFO],_am[4*NIFO],_AM[4*NIFO],_u[4*NIFO],_v[4*NIFO];

  double a,b,h,H,psi,gg,E,ee,EE,fp2,fx2;

  wavearray<double> cid;                        // buffers for cluster ID
  wavearray<double> cTo;                        // buffers for cluster time
  netpixel* pix;
  std::vector<int>* vint;
  std::vector<int> level;

//+++++++++++++++++++++++++++++++++++++++
// loop over clusters
//+++++++++++++++++++++++++++++++++++++++

  double cF = NET.MRA ? 1. : 2.;     			 // NDM normalization Coefficient

  netcluster* pwc = &NET.wc_List[LAG];
  pwc->print();

  cid = pwc->get((char*)"ID",  0,'S',0);                 // get cluster ID
  cTo = pwc->get((char*)"time",0,'L',0);                 // get cluster time

  int K = cid.size();
  for(int k=0; k<K; k++) {                               // loop over clusters

    int id = size_t(cid.data[k]+0.1);
    if(id!=CID) continue;

    if(pwc->sCuts[id-1] != 0) continue;                  // skip rejected/processed clusters

    vint = &(pwc->cList[id-1]);                          // pixel list

    int V = vint->size();
    if(!V) continue;

    layers.resize(V);  
    time.resize(V);  
    frequency.resize(V);  
    likelihood.resize(V);  
    null.resize(V);  
    for(int k=0;k<NIFO;k++) {
      am[k].resize(V); AM[k].resize(V);
      fp[k].resize(V); fx[k].resize(V);
       u[k].resize(V);  v[k].resize(V);
      nv[k].resize(V);  
    }

    int L=0;
    for(int j=0; j<V; j++) {
      netpixel* pix = pwc->getPixel(id,j);
      if(!pix->core) continue;
      //cout << j << " likelihood : " << pix->likelihood << endl;

      for(int n=0; n<nIFO; n++) {
        b    = pix->getdata('N',n);        // noise rms
        gg  += 1./b/b;                     // noise normalization
      }
      gg = sqrt(gg);

      for(int n=0; n<nIFO; n++) {                 
        _am[n] = pix->getdata('S',n);             		// snr amplitude
        _AM[n] = pix->getdata('P',n);       			// snr 90 degrees amplitude
        _u[n]  = pix->getdata('W',n)/pix->getdata('N',n);
        _v[n]  = pix->getdata('U',n)/pix->getdata('N',n);
         b     = pix->getdata('N',n)*gg;    			// noise rms
        nv[n][L] = b; 
        //cout << j << " " << n << " " << _am[n] << " " << _AM[n] << " " 
        //     << _u[n]*pix->getdata('N',n) << " " << _v[n]*pix->getdata('N',n) << endl;
         pd    = NET.getifo(n);
        _fp[n] = pd->mFp.get(pix->theta,pix->phi);
        _fx[n] = pd->mFx.get(pix->theta,pix->phi);
        _fp[n]/=b;
        _fx[n]/=b;
      }

      // DPF and DSP transformations
      //NET.dspx(_fp, _fx, _am, _AM, _u, _v);

      // save data
      for(int n=0; n<nIFO; n++) {              
        fp[n][L] = _fp[n]; fx[n][L] = _fx[n];
        am[n][L] = _am[n]; AM[n][L] = _AM[n];
        // u[n][L] = _u[n];   v[n][L] = _v[n];

         u[n][L] = _u[n]*pix->getdata('N',n);   
         v[n][L] = _v[n]*pix->getdata('N',n);
      }

      ee=EE=0;
      for(int n=0; n<nIFO; n++) {
        ee += am[n][L]*am[n][L]; 
        EE += AM[n][L]*AM[n][L];
      }

      layers[L]     = pix->layers;
      time[L]       = pix->time;
      frequency[L]  = pix->frequency;
      likelihood[L] = pix->likelihood;
      null[L]       = (ee+EE)/cF-pix->likelihood;

      //cout << j << " " << pix->layers << endl;
      bool check=false;
      for(int k=0;k<level.size();k++) if(level[k]==pix->layers) {check=true;break;}
      if(!check) level.push_back(pix->layers);

      L++;
    }

    layers.resize(L);  
    time.resize(L);  
    frequency.resize(L);  
    likelihood.resize(L);  
    null.resize(L);  
    for(int k=0;k<NIFO;k++) {
      am[k].resize(L); AM[k].resize(L);
      fp[k].resize(L); fx[k].resize(L);
       u[k].resize(L);  v[k].resize(L);
      nv[k].resize(L);  
    }
  }

  for(int k=0;k<level.size();k++) level[k]-=1;
  for(int k=0;k<level.size();k++) cout << k << " " << level[k] << endl;

  int ifo=0;
  DrawSignal(ifo,level);
  DrawMonster(ifo,level);

  //exit(0);
}

void ReadConfig(TString fName, network* NET, CWB::config* cfg) {

  TFile* rfile = new TFile(fName);
  if(rfile==NULL) {cout << "Error opening root file : " << fName.Data() << endl;exit(1);}
  // read network object
  if(rfile->Get("network")!=NULL) {
    *NET = *(network*)rfile->Get("network");
  } else {
    cout << "Error : net is not contained in root file " << fName.Data() << endl; exit(1);
  }
  // read config object
  if(rfile->Get("config")!=NULL) {
    cfg = *(CWB::config*)rfile->Get("config");
  } else {
    cout << "Error : config is not contained in root file " << fName.Data() << endl; exit(1);
  }
  rfile->Close();

  int nIFO = NET->ifoListSize();
  detector* pD[NIFO_MAX]; //! pointers to detectors
  for(int n=0; n<nIFO; n++) pD[n] = NET->getifo(n);

  // restore skymaps
  if(cfg->healpix) NET->setSkyMaps(int(cfg->healpix));
  else             NET->setSkyMaps(cfg->angle,cfg->Theta1,cfg->Theta2,cfg->Phi1,cfg->Phi2);
  NET->setAntenna();
  NET->setDelay(pD[0]->Name);

  // there is a issue in network class (no TClass for char*)
  NET->ifoName.clear();
  for(int n=0; n<nIFO; n++) NET->ifoName.push_back(pD[n]->Name);
}

void PrintConfig(network* NET, CWB::config* cfg) {

  // Print Configuration
  cout << "\n Analysis : " << cfg->analysis << endl;

  if(cfg->search=='E' || cfg->search=='E') cout<<"\n un-modeled search (Energy): "<<cfg->search<<endl;
  if(cfg->search=='b' || cfg->search=='B') cout<<"\n un-modeled search (single stream): "<<cfg->search<<endl;
  if(cfg->search=='r' || cfg->search=='R') cout<<"\n un-modeled search (dual stream): "<<cfg->search<<endl;
  if(cfg->search=='i' || cfg->search=='I') cout<<"\n elliptical polarisation: "<<cfg->search<<endl;
  if(cfg->search=='g' || cfg->search=='G') cout<<"\n circular polarisation: "<<cfg->search<<endl;
  if(cfg->search=='s' || cfg->search=='S') cout<<"\n linear polarisation: "<<cfg->search<<endl;

  double mTau=NET->getDelay(const_cast<char*>("MAX"));  // maximum time delay
  double dTau=NET->getDelay(const_cast<char*>(""));     // time delay difference
  cout<<"maximum time delay between detectors: "<<mTau<<endl;
  cout<<"       maximum time delay difference: "<<dTau<<endl;
  cout<<"                  skymap search mode: "<<cfg->mode<<endl;
  if(cfg->healpix) {
    cout<<"                       HEALPix order: "<<cfg->healpix<<endl;
  } else {
    cout<<"           skymap angular resolution: "<<cfg->angle<<endl;
    cout<<"          skymap size in polar angle: "<<cfg->Theta1<<", "<<cfg->Theta2<<endl;
    cout<<"      skymap size in azimuthal angle: "<<cfg->Phi1<<", "<<cfg->Phi2<<endl;
  }
  cout<<"                    netRHO and netCC: "<<cfg->netRHO<<", "<<cfg->netCC<<endl;
  cout<<"              regulator delta, local: "<<cfg->delta<<" "<<NET->local<<endl;
  cout<<endl;
  if(cfg->mask>0.)
    cout<<"              earth sky mask applied; "<<cfg->skyMaskFile<<endl;
  if(strlen(cfg->skyMaskCCFile)>0)
    cout<<"          celestial sky mask applied; "<<cfg->skyMaskCCFile<<endl;
  cout<<endl;

  detector* pD[NIFO_MAX];                       //! pointers to detectors
  int nIFO = NET->ifoListSize();
  for(int n=0; n<nIFO; n++) pD[n] = NET->getifo(n);

  cout<<" network of ";
  for(int i=0; i<cfg->nIFO; i++) cout<<pD[i]->Name<<" ";
  cout<<" detectors\n\n";

  cout << endl << "Print list of injected MDC" << endl << endl;
  cout.precision(14);
  for(int k=0;k<(int)NET->mdcList.size();k++) cout << k << " mdcList " << NET->mdcList[k] << endl;
  for(int k=0;k<(int)NET->mdcTime.size();k++) cout << k << " mdcTime " << NET->mdcTime[k] << endl;
  for(int k=0;k<(int)NET->mdcType.size();k++) cout << k << " mdcType " << NET->mdcType[k] << endl;

  // fill injection object

  injection INJ(nIFO);
  size_t mdcID, ID, M;
  double T, injTime, TAU, pcc;
  wavearray<double> time_net;
  int ifo=0;

  double gps = NET->getifo(ifo)->getTFmap()->start();
  netcluster* p = NET->getwc(LAG);       // pointer to netcluster
  if(!p->size()) continue;

  time_net.append(p->get((char*)"time",ifo,'L',0));

  for(int n=0;n<time_net.size();n++) {

    double time = time_net.data[n] + gps;

    injTime = 1.e12;
    injID   = -1;
    M = NET->mdc__IDSize();
    for(int m=0; m<M; m++) {
      mdcID = NET->getmdc__ID(m);
      T = fabs(time - NET->getmdcTime(mdcID));
      if(T<injTime && INJ.fill_in(NET,mdcID)) {
        injTime = T;
        injID = mdcID;
      }
      printf("Injection ID : %d  %12.3f  %12.3f\n",mdcID,NET->getmdcTime(mdcID),T);
    }
  }
}

void DrawSignal(int ifo, vector<int> level) {

  int nRES = level.size();

  WSeries<double> w00;
  WSeries<double> w90;
  WDM<double>* wdm00;
  WDM<double>* wdm90;
  WDM<double>* wdtf;

  wavearray<double> x,z;
  z.rate(RATE); z.resize(200*RATE); z.start(0); z=0;
  x = z;

  for(int n=0;n<nRES;n++) {

    x=0;

    double dF = RATE/(2*level[n]);         // frequency resolution
    double dT = level[n]/RATE;             // time resolution

    cout << "WDM Decomposition Level : " << level[n]
         << " dF : " << dF << " Hz " << " dT : " << 1000*dT << " ms " <<endl;

    wdtf = new WDM<double>(level[n], level[n], 4, 10);
    w00.Forward(x, *wdtf);
    w90.Forward(x, *wdtf);

    double flength = wdtf->m_H/RATE;
    cout << "WDM Filter length : " << wdtf->m_H << " samples " << flength << " sec " << endl; 

    // Set pixels
    wdm00 = (WDM<double>*) w00.pWavelet;
    wdm90 = (WDM<double>*) w90.pWavelet;
    int M = w00.getLevel();
    int L = wdm00->maxLayer();

    int mF = int(w00.size()/wdm00->nSTS);
    int nTC = w00.size()/(M+1)/mF;                   // # of Time Coefficients
    // Extract map00 & map90 : size = (M+1)*nTC
    double* map00 = wdm00->pWWS;
    double* map90 = map00 + (mF-1)*(M+1)*nTC;

    double* MAP00 = wdm90->pWWS;
    double* MAP90 = MAP00 + (mF-1)*(M+1)*nTC;

    int cnt=0;
    for(int m=0;m<u[ifo].size();m++) {
      if(layers[m]!=(level[n]+1)) continue;
      map00[time[m]]=u[ifo][m];
      map90[time[m]]=v[ifo][m];
      MAP00[time[m]]=u[ifo][m];
      MAP90[time[m]]=v[ifo][m];
      cnt++;
    }
    cout << "level : " << level[n] << " entries " << cnt << endl;

    w00.Inverse();
    z+=(wavearray<double>)w00;
    w90.Inverse(-2);
    z+=(wavearray<double>)w90;

    delete wdtf;
  }

#ifdef WDM_SPECTROGRAM

  int nfact=4;
  int nfft=nfact*512;
  int noverlap=nfft-10;
  double fparm=nfact*6;
  int ystart = 85*z.rate();
  int ystop  = 105*z.rate();
  ystart-=nfft;
  ystop+=nfft;
  int ysize=ystop-ystart;
  wavearray<double> y;y.resize(ysize);y.rate(z.rate());y.start(ystart/z.rate());

  for(int i=0;i<(int)y.size();i++) y.data[i]=z.data[i+ystart];

  CWB::STFT stft(y,nfft,noverlap,"energy","gauss",fparm);

  TCanvas* canvas;
  double tstart = nfft/z.rate()+ystart/z.rate();
  double tstop = (ysize-nfft)/z.rate()+ystart/z.rate();
  stft.SetTitle("Monster Event Spectrogram (Energy)");
  stft.Draw(tstart,tstop,0,1024,0,0,1);
  TString fName = "monster_event_spectrogram.png";
  canvas = stft.GetCanvas();
  stft.Print(fName); 

#endif

#ifdef WDM_TF
  int Mlevel=0;
  for(int n=0;n<nRES;n++) if(Mlevel<level[n]) Mlevel=level[n];
  wdtf = new WDM<double>(Mlevel, Mlevel, 4, 10);
  w00.Forward(z, *wdtf);

  // Plot WDM Scalogram
  watplot WTS(const_cast<char*>("wts"));
  //scalogram maps
  double start = w00.start();
  double stop  = w00.start()+w00.size()/w00.rate();
  double flow  = 64;
  double fhigh = 2048;
  WTS.plot(&w00, 2, start, stop,const_cast<char*>("COLZ"));
  WTS.hist2D->GetXaxis()->SetRangeUser(85, 105);
  WTS.hist2D->GetYaxis()->SetRangeUser(32, 1000);
  // dump spectrum
  char fname[1024];
  sprintf(fname,"monster_event_scalogram.root");
  cout << endl << "Dump WDM Scalogram : " << fname << endl << endl;
  WTS.canvas->Print(fname);
  sprintf(fname,"monster_event_scalogram.png");
  WTS.canvas->Print(fname);
#endif

  char gtitle[256];
  TString gfile;

#ifdef READ_MDC
  ReadMDC(MDC_FILE_NAME, x);
#endif

#ifdef WDM_T
  watplot tplot(const_cast<char*>("tplot"),200,20,800,500);
  sprintf(gtitle,"WDM : Reconstructed Signal");
  tplot.gtitle(gtitle,"time(sec)","amplitude");
  tplot.goptions("alp", 1, 85., 105.);
  // draw signal
#ifdef READ_MDC
  x >> tplot;
#endif
  z >> tplot;
  // save plot to file
  gfile="monster_event_time.png";
  tplot >> gfile;
  gfile="monster_event_time.root";
  tplot >> gfile;
#endif

#ifdef WDM_F
  watplot fplot(const_cast<char*>("fplot"),200,20,800,500);
  fplot.gtitle(gtitle,"frequency (Hz)","strain/#sqrt{Hz}");
  fplot.goptions("alp logy", 1, 85., 105., true, 32,1024);
#ifdef READ_MDC
  x >> fplot;
#endif
  z >> fplot;
  // save plot to file
  gfile="monster_event_freq.png";
  fplot >> gfile;
  gfile="monster_event_freq.root";
  fplot >> gfile;
#endif

}

void DrawMonster(int ifo, vector<int> level, int iRATE, bool like) {

  int nRES = level.size();

  int Mlevel=0;
  for(int n=0;n<nRES;n++) if(Mlevel<level[n]) Mlevel=level[n];
  int mlevel=1000;
  for(int n=0;n<nRES;n++) if(mlevel>level[n]) mlevel=level[n];
  int mrate=RATE/Mlevel;
  int Mrate=RATE/mlevel;

  cout << "mrate  : " << mrate <<  " Mrate  : " << Mrate << endl;
  cout << "mlevel : " << mlevel << " Mlevel : " << Mlevel << endl;

  TCanvas* canvas;
  canvas = new TCanvas("Monster Event Display", "LVC experiment", 300,40, 800, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);
  canvas->SetRightMargin(0.10);
  canvas->SetLeftMargin(0.10);
  canvas->SetBottomMargin(0.13);
  canvas->SetBorderMode(0);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);

  TH2F* hist2D = NULL;
  hist2D=new TH2F("WTF", "WTF", 105*Mrate, 0, 105, Mlevel+1, 0, RATE/2);
  hist2D->SetXTitle("time, sec");
  hist2D->SetYTitle("frequency, Hz");

  hist2D->SetStats(kFALSE);
  hist2D->SetTitleFont(12);
  hist2D->SetFillColor(kWhite);

  hist2D->GetXaxis()->SetNdivisions(506);
  hist2D->GetXaxis()->SetLabelFont(42);
  hist2D->GetXaxis()->SetLabelOffset(0.014);
  hist2D->GetXaxis()->SetTitleOffset(1.4);
  hist2D->GetYaxis()->SetTitleOffset(1.2);
  hist2D->GetYaxis()->SetNdivisions(506);
  hist2D->GetYaxis()->SetLabelFont(42);
  hist2D->GetYaxis()->SetLabelOffset(0.01);
  hist2D->GetZaxis()->SetLabelFont(42);
  hist2D->GetZaxis()->SetNoExponent(false);
  hist2D->GetZaxis()->SetNdivisions(506);

  hist2D->GetXaxis()->SetTitleFont(42);
  hist2D->GetXaxis()->SetTitle("Time (sec)");
  hist2D->GetXaxis()->CenterTitle(true);
  hist2D->GetYaxis()->SetTitleFont(42);
  hist2D->GetYaxis()->SetTitle("Frequency (Hz)");
  hist2D->GetYaxis()->CenterTitle(true);

  hist2D->GetZaxis()->SetTitleOffset(0.6);
  hist2D->GetZaxis()->SetTitleFont(42);
  //hist2D->GetZaxis()->SetTitle(ztitle);
  hist2D->GetZaxis()->CenterTitle(true);

  hist2D->GetXaxis()->SetLabelSize(0.03);
  hist2D->GetYaxis()->SetLabelSize(0.03);
  hist2D->GetZaxis()->SetLabelSize(0.03);

  //hist2D->GetXaxis()->SetRangeUser(mtime-0.1,Mtime+0.1);
  //hist2D->GetYaxis()->SetRangeUser(0,Mfrequency+32);

  hist2D->GetXaxis()->SetRangeUser(85,105);
  hist2D->GetYaxis()->SetRangeUser(0,1000);

  hist2D->SetBinContent(1,1,1);
  hist2D->SetBinContent(100*Mrate,1,1);

  double Null=0;
  double Likelihood=0;
  for(int n=0;n<u[ifo].size();n++) {
    int irate = RATE/(layers[n]-1); 
    int M=Mrate/irate;
    int K=Mlevel/(layers[n]-1);
    double itime=((double)time[n]/(double)irate)/layers[n];
    int i=itime*Mrate;
    int j=frequency[n]*K;
    //cout << n << " " << i << " " << j << " " << time[n] << " " << frequency[n] << " " << likelihood[n] << endl;
    if(irate!=iRATE && iRATE!=0) continue;
    Null+=null[n];
    Likelihood+=likelihood[n];
    int L=0;int R=1;while (R < irate) {R*=2;L++;}
    for(int m=0;m<M;m++) {
      for(int k=0;k<K;k++) {
        if(like) hist2D->SetBinContent(i+1+m,2*j+1+k,likelihood[n]);
        else     hist2D->SetBinContent(i+1+m,2*j+1+k,null[n]);
        //hist2D->SetBinContent(i+1+m,2*2*j+1+k,L);
      }
    }
  }

  char title[256];
  sprintf(title,"Monster Event Display - Rate : [%d:%d] - Likelihood %3.0f - Null %3.0f",mrate,Mrate,Likelihood,Null);
  hist2D->SetTitle(title);

  cout << "Event Likelihood : " << Likelihood << endl;
  cout << "Event Null       : " << Null << endl;

  hist2D->Draw("COLZ");

  TString gfileName("monster_event_likelihood.root");
  cout << gfileName.Data() << endl;
  canvas->Print(gfileName);
  gfileName.ReplaceAll(".root",".gif");
  canvas->Print(gfileName);
  TString pfileName=gfileName;
  pfileName.ReplaceAll(".gif",".png");
  char cmd[1024];
  sprintf(cmd,"convert %s %s",gfileName.Data(),pfileName.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"rm %s",gfileName.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);

}

void ReadMDC(TString fName, wavearray<double>& x) {

  // READ MDC DATA
  // data are loaded from root file
  TFile* jfile = new TFile(fName,"UPDATE");
  if(jfile==NULL||!jfile->IsOpen())
    {cout << "ReadMDC - Error opening root file : " << fName.Data() << endl;exit(1);}
  wavearray<double>* mdc = (wavearray<double>*)jfile->Get("L1");
  if(mdc==NULL)
    {cout << "ReadMDC - Error : mdc not present in file : " << fName.Data() << endl;exit(1);}

  x=*mdc;
  delete mdc;

  cout << "MDC " << x.size() << " " << x.rate() << " "
                 << x.size()/x.rate() << " " << int(x.start()) << endl;

  x.start(0);

  double rateFactor=2;  	// sample rate = 16384
//  double rateFactor=1;                // sample rate = 4096
  double factor=60;
  double mdcFactor=0.17207748345771;
  x*=factor*mdcFactor*rateFactor;;
  cout << "MDC factor " << factor*mdcFactor*rateFactor << endl;

  return;
}

