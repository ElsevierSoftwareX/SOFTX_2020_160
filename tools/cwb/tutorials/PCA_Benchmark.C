#include <vector>

// --------------------------------------------------------
// PRINCIPAL COMPONENT ANALYSIS
// --------------------------------------------------------

#define PCA_NRES	7	// number of resolutions
#define PCA_IRES	4	// initial resolution

#define PCA_MAX		60	// number of max PC extractions
#define PCA_PREC	0.01 	// precision = Eresidual/Esignal

#define ACORE		1.7	// set selection threshold : Ethr = 2*ACORE*ACORE

//#define APPLY_TSHIFT		// enable/disable (uncomment/comment) find best wavelet time shift

//#define OUTPUT_FILE	"pca_benchmark.txt"	
				// enable/disable (uncomment/comment) write precision vs nPC to ascii file

// --------------------------------------------------------
// SIGNAL
// --------------------------------------------------------

#define WAVE_TYPE	"INSPIRAL"	// inspiral
//#define WAVE_TYPE	"SGE"		// sine gaussian
//#define WAVE_TYPE	"WNB"		// white noise burst
//#define WAVE_TYPE	"GA"		// gaussian

#define WAVE_LENGTH     6               // Is effective only for WAVE_TYPE = SGE,WNB,GA (def=1sec)

#define WAVE_POL	"hp"	// select hp,hx component

#define WAVE_SNR	30	// signal SNR

// --------------------------------------------------------
// WDM
// --------------------------------------------------------

#define WDM_NU		6	// WDM mayer order 
				// NOTE : wp_xtalk coefficients are defined only for WDM_NU=2,4,6
#define WDM_PREC	10	// WDM precision
#define WDM_TDSIZE	12	// Time Delay filter size
//#define WDM_SCRATCH	14.0	// WDM TF segment scratch (sec)
#define WDM_SCRATCH	27.0	// WDM TF segment scratch (sec)
//#define WDM_SCRATCH	54.0	// WDM TF segment scratch (sec)

//#define LOAD_XTALK_CATALOG	// enable/disable (uncomment/comment) xcatalog
//#define XTALK_CATALOG "wdmXTalk/OverlapCatalog16-1024.bin"			// WDM_NU = 4
#define XTALK_CATALOG "wdmXTalk/OverlapCatalog-ilLev4-hLev10-iNu6-P10.bin"	// WDM_NU = 6

// --------------------------------------------------------
// WAVELET PACKET
// --------------------------------------------------------

#define NWP		4	// number of WP types used in the PCA

#define WP_SINGLE 		// enable/disable (uncomment/comment) single wavelet
//#define WP_DIAG 		// enable/disable (uncomment/comment) diagonal WP
//#define WP_HORIZ		// enable/disable (uncomment/comment) horizontal WP
//#define WP_VERT 		// enable/disable (uncomment/comment) vertical WP
//#define WP_BDIAG		// enable/disable (uncomment/comment) back-diagonal WP

				// enable/disable (uncomment/comment) wavelet packet analysis
    				// NOTE : for more infos see InitXTALK function
				// ----------------------------------------------------------
                                // Phase 00        Phase 90
                                // - - c*w         - - C*W
                                // - b*v -         - B*V -
                                // a*u - -         A*U - -
				//
				// WAVELET PACKET BASE
        			// parity = odd  -> z=U+v+W  , Z=u-V+w
        			// parity = even -> z=U-v+W  , Z=u+V+w
				// (v,U) = (v,W) = -(V,u) = -(V,w) = WP_XTALK
                                // (W,u) = (U,w) = 0
        			// (z,z) = (Z,Z) = 3+4*WP_XTALK
				// (z,Z) = 0
				//
				// DATA
        			// x=a*u+b*v+c*w 
				// A00 = (x,z)/(z,z)	// A00 is the amplitude of x respect to z
				// A90 = (x,Z)/(Z,Z)	// A90 is the amplitude of x respect to Z
        			// parity = odd 
				// A00 = (x,z) = (x,U+v+W) = A+b+C
				// A90 = (x,Z) = (x,u-V+w) = a-B+c
        			// parity = even 
				// A00 = (x,z) = (x,U-v+W) = A-b+C
				// A90 = (x,Z) = (x,u+V+w) = a+B+c
				// ----------------------------------------------------------

// --------------------------------------------------------
// PLOTS
// --------------------------------------------------------

//#define PLOT_MONSTER		// enable/disable (uncomment/comment) show likelihood TFmap
#define PLOT_REC_VS_INJ     	// enable/disable (uncomment/comment) show time rec vs inj waveform 

//#define PLOT_TFMAP 		// enable/disable (uncomment/comment) show TF map @ PLOT_TFRES resolution
#define PLOT_TFRES	6	// TF frequency resolution showed
#define PLOT_MINFREQ	0	// min freq shoed in the TF plots
#define PLOT_MAXFREQ	256	// max freq shoed in the TF plots
#define PLOT_TFTYPE	2
				// 0 : sqrt((E00+E90)/2)
				// 1 : (E00+E90)/2
				// 2 : sqrt((E00+E90)/2)
				// 3 : amplitude:00
				// 4 : energy:00
				// 5 : |amplitude:00|
				// 6 : amplitude:90
				// 7 : energy:90
				// 8 : |amplitude:90|

// --------------------------------------------------------
// DECLARATIONS
// --------------------------------------------------------

WDM<double>* wdm[PCA_NRES];   	// define a WDM transform containers
WSeries<double> wsig[PCA_NRES]; // waveform TF map containers
wavearray<double> sig;		// signal
gwavearray<double>* grec;	// used to plot rec
gwavearray<double>* gsig;	// used to plot sig
watplot *WTS[3];		// plot objects
double wp_xtalk[4];
bool wp_mask[5];

#ifdef LOAD_XTALK_CATALOG
monster wdmMRA;                 // wdm multi-resolution analysis
void CheckXTalkCatalog(int l_low, int l_high);
double GetPrincipalFactor(int nmax, int mmax, int rmax);
#endif

int GetPrincipalComponent(double &amax, int &nmax, int &mmax, int &rmax); 

void InitXTALK();

void PCA_Benchmark() {

  //
  // Principal Component Analysis : Test Macro
  // Extraction of Pincipal Components from multiresolution TF maps
  // Author : Gabriele Vedovato

  if(TString("hp")!=WAVE_POL && TString("hx")!=WAVE_POL) {
    cout << "Error - Bad WAVE_POL value : must be hp or hx" << endl;
    exit(1);
  }
#ifdef PLOT_TFMAP
  if(PLOT_TFRES<0 || PLOT_TFRES>=PCA_NRES) {
    cout << "Error - Bad PLOT_TFRES value : must be >=0 && < " << PCA_NRES << endl;
    exit(1);
  }
#endif
#ifdef OUTPUT_FILE
  ofstream out;
  out.open(OUTPUT_FILE,ios::out); 
#endif

  for(int l=0;l<=NWP;l++) wp_mask[l] = false;
#ifdef WP_SINGLE
  wp_mask[0] = true;
#endif
#ifdef WP_DIAG
  wp_mask[1] = true;
#endif
#ifdef WP_HORIZ
  wp_mask[2] = true;
#endif
#ifdef WP_VERT
  wp_mask[3] = true;
#endif
#ifdef WP_BDIAG
  wp_mask[4] = true;
#endif

  InitXTALK();	// init wp_xtalk coefficients

  CWB::mdc MDC; 
  MDC.SetInjLength(WAVE_LENGTH);

  if(WAVE_TYPE == "INSPIRAL") {
    // ---------------------------------
    // set inspiral parameters
    // ---------------------------------
    TString inspOptions="";
    inspOptions = "--time-step 60.0 --time-interval 0 ";
    inspOptions+= "--l-distr random ";
    inspOptions+= "--gps-start-time 931072160 --gps-end-time 931072235 ";
    inspOptions+= "--d-distr volume --m-distr totalMassRatio --i-distr uniform ";

//  inspOptions+= "--f-lower 80.000000 ";	// 0.5 sec
//  inspOptions+= "--f-lower 50.000000 ";	// 1.5 sec
//  inspOptions+= "--f-lower 40.000000 ";	// 3   sec
    inspOptions+= "--f-lower 30.000000 ";	// 6.5 sec
//  inspOptions+= "--f-lower 20.000000 ";	// 19  sec
//  inspOptions+= "--f-lower 10.000000 ";	// 122 sec
    inspOptions+= "--min-mass1 5.0 --max-mass1 5.0 ";
    inspOptions+= "--min-mass2 5.0 --max-mass2 5.0 ";
    inspOptions+= "--min-mtotal 10.0 --max-mtotal 10.0 ";

/*
    inspOptions+= "--f-lower 40.000000 ";
    inspOptions+= "--min-mass1 10.0 --max-mass1 10.0 ";
    inspOptions+= "--min-mass2 10.0 --max-mass2 10.0 ";
    inspOptions+= "--min-mtotal 20.0 --max-mtotal 20.0 ";
*/
/*
    inspOptions+= "--f-lower 24.000000 ";
    inspOptions+= "--min-mass1 42.1 --max-mass1 42.1 ";
    inspOptions+= "--min-mass2 29.6 --max-mass2 29.6 ";
    inspOptions+= "--min-mtotal 71.7 --max-mtotal 71.7 ";
*/
/*
    inspOptions+= "--f-lower 10.000000 ";
    inspOptions+= "--min-mass1 75.000000 --max-mass1 75.000000 ";
    inspOptions+= "--min-mass2 75.000000 --max-mass2 75.000000 ";
    inspOptions+= "--min-mtotal 150.000000 --max-mtotal 150.000000 ";
*/

    inspOptions+= "--m-distr componentMass ";
    inspOptions+= "--min-distance 1000000.0 --max-distance 1500000.0 ";
    inspOptions+= "--approximant EOBNRv2pseudoFourPN --disable-spin ";
    inspOptions+= "--taper-injection start --seed 123456789 ";
    inspOptions+= "--dir ./ ";
    inspOptions+= "--output inspirals.xml ";		// set output xml file

    MDC.SetInspiral("EOBNRv2",inspOptions);

    // Get the first waveform hp,hx components starting from gps = 931072160
    sig = MDC.GetInspiral(WAVE_POL,931072160,931072235);
  }

  if(WAVE_TYPE == "SGE") {
    vector<mdcpar> par(2);
    par[0].name="frequency"; par[0].value=100.;
    par[1].name="Q"; par[1].value=8.9;
    MDC.AddWaveform(MDC_SGE, par);
    MDC.Print();
    waveform wf = MDC.GetWaveform(0,0);
    if(TString(WAVE_POL)=="hp") sig = wf.hp;
    if(TString(WAVE_POL)=="hx") sig = wf.hx;
  }

  if(WAVE_TYPE == "WNB") {
    vector<mdcpar> par(6);
    par[0].name="frequency"; par[0].value=250;
    par[1].name="bandwidth"; par[1].value=100;
    par[2].name="duration";  par[2].value=0.1;
    par[3].name="pseed";     par[3].value=111;
    par[4].name="xseed";     par[4].value=222;
    par[5].name="mode";      par[5].value=1;  // simmetric
    MDC.AddWaveform(MDC_WNB, par);
    MDC.Print();
    waveform wf = MDC.GetWaveform(0,0);
    if(TString(WAVE_POL)=="hp") sig = wf.hp;
    if(TString(WAVE_POL)=="hx") sig = wf.hx;
  }

  if(WAVE_TYPE == "GA") {
    vector<mdcpar> par(1);
    par[0].name="duration"; par[0].value=0.004;
    MDC.AddWaveform(MDC_GA, par);
    MDC.Print();
    waveform wf = MDC.GetWaveform(0,0);
    if(TString(WAVE_POL)=="hp") sig = wf.hp;
    if(TString(WAVE_POL)=="hx") {
      cout << "Error - WAVE=GA -> hx=0" << endl;
      exit(1);  
    }
  }

  cout << endl << "sig size : " << sig.size() << " sig rate : " 
       << sig.rate() << " sig start : " << (int)sig.start() 
       << " sig length : " << sig.size()/sig.rate() << " sec" << endl;
  sig.start(0);		// set start to 0 (needed by draw method)

  // add scratch array
  int iscratch = WDM_SCRATCH*sig.rate();
  wavearray<double> X(sig.size()+2*iscratch);
  X.rate(sig.rate());
  X=0.; for(int i=0;i<sig.size();i++) X[i+iscratch]=sig[i]; 
  // apply time shift
  //int jshift = 0.05*sig.rate();
  //for(int i=jshift;i<sig.size();i++) X[i+iscratch]=sig[i-jshift]; 
  sig=X;

  // resample data 16384 Hz -> 2048 Hz
  sig.FFTW(1);
  sig.resize(sig.size()/8);
  sig.FFTW(-1);
  sig.rate(sig.rate()/8.);

  // normalize sig energy to WAVE_SNR^2  
  double Esig=0.; for(int i=0;i<sig.size();i++) Esig+=sig[i]*sig[i];
  for(int i=0;i<sig.size();i++) sig[i]/=sqrt(Esig)/sqrt(WAVE_SNR*WAVE_SNR);
  Esig=WAVE_SNR*WAVE_SNR;

  // print pixel resolutions
  cout << endl;
  for(int level=PCA_NRES+PCA_IRES-1; level>=PCA_IRES; level--) {
    int rateANA = int(sig.rate()); 
    int layers = level>0 ? 1<<level : 0;
    int rate  = rateANA>>level;
    cout << "level : " << level << "\t rate(hz) : " << rate << "\t layers : " << layers
         << "\t df(hz) : " << rateANA/2./double(1<<level)
         << "\t dt(ms) : " << 1000./rate << endl;
  }
  cout << endl;

#ifdef APPLY_TSHIFT
  cout << endl << "find best wavelet time shift ENABLED" << endl << endl;
#else
  cout << endl << "find best wavelet time shift DISABLED" << endl << endl;
#endif
  bool wpEnabled=false;
  for(int l=0;l<=NWP;l++) if(wp_mask[l]) wpEnabled=true;
  if(wpEnabled) {
    cout << endl;
    if(wp_mask[0]) cout << "wavelet packet       single  ENABLED" << endl;
    if(wp_mask[1]) cout << "wavelet packet      diagonal ENABLED" << endl;
    if(wp_mask[2]) cout << "wavelet packet    horizontal ENABLED" << endl;
    if(wp_mask[3]) cout << "wavelet packet      vertical ENABLED" << endl;
    if(wp_mask[4]) cout << "wavelet packet back-diagonal ENABLED" << endl;
    cout << endl;
  } else {
    cout << endl << "Error : all WP options are disabled !!!" << endl << endl;
    exit(1);
  }

  // compute WDM with PCA_NRES TF resolutions
  cout << "Init WDM ..." << endl;
  for(int i=PCA_NRES-1;i>=0;i--) {
    int j=i+PCA_IRES;
    wdm[i] = new WDM<double>(1<<j, 1<<j, WDM_NU, WDM_PREC);   // define a WDM transform 

    // check if filter lenght is less than cwb scratch length
    double wdmFLen = double(wdm[i]->m_H)/sig.rate();    // sec
    if(wdmFLen > WDM_SCRATCH) {
       cout << endl;
       cout << "Error - filter length must be <= WDM_SCRATCH !!!" << " RES = " << j << endl;
       cout << "filter length : " << wdmFLen << " sec" << endl;
       cout << "cwb   scratch : " << WDM_SCRATCH << " sec" << endl;
       exit(1);
    }
  }

#ifdef LOAD_XTALK_CATALOG

  char filter_dir[1024];
  if(gSystem->Getenv("HOME_WAT_FILTERS")==NULL) {
    cout << "Error : environment HOME_WAT_FILTERS is not defined!!!" << endl;exit(1);
  } else {
    strcpy(filter_dir,TString(gSystem->Getenv("HOME_WAT_FILTERS")).Data());
  }

  char MRAcatalog[1024];
  sprintf(MRAcatalog,"%s/%s",filter_dir,XTALK_CATALOG);
  cout << "cwb2G::Init - Loading catalog of WDM cross-talk coefficients ... " << endl;
  cout << MRAcatalog << endl;
  CWB::Toolbox::checkFile(MRAcatalog);

  wdmMRA.read(MRAcatalog);

  // check if analysis layers are contained in the MRAcatalog
  // level : is the decomposition level
  // layes : are the number of layers along the frequency axis rateANA/(rateANA>>level)

  int l_low  = PCA_IRES;
  int l_high = PCA_IRES+PCA_NRES-1;

  int check_layers=0;
  for(int level=l_high; level>=l_low; level--) {
    int layers = level>0 ? 1<<level : 0;
    for(int j=0;j<wdmMRA.nRes;j++) if(layers==wdmMRA.layers[j]) check_layers++;
  }

  if(check_layers!=PCA_NRES) {
    cout << "Error : analysis layers do not match the MRA catalog" << endl;
    cout << endl << "analysis layers : " << endl;
    for(int level=l_high; level>=l_low; level--) {
      int layers = level>0 ? 1<<level : 0;
      cout << "level : " << level << " layers : " << layers << endl;
    }
    cout << endl << "MRA catalog layers : " << endl;
    for(int i=0;i<wdmMRA.nRes;i++)
       cout << "layers : " << wdmMRA.layers[i] << endl;
    exit(1);
  }

  //CheckXTalkCatalog(l_low, l_high);
#endif

  // compute TF of sig signal
  cout << "Apply WDM Transform to signal ..." << endl;
  for(int i=0;i<PCA_NRES;i++) {
    wsig[i].Forward(sig, *wdm[i]);             	// apply the WDM to the time series
  }

  // display TF sig
#ifdef PLOT_TFMAP		
  double dt0 = 1./wsig[PLOT_TFRES].rate();
  double tbeg0 = WDM_SCRATCH-2*dt0;
  double tend0 = dt0*wsig[PLOT_TFRES].size()/2-WDM_SCRATCH;
  WTS[0] = new watplot(const_cast<char*>("wts0"));
  WTS[0]->plot(&wsig[PLOT_TFRES], PLOT_TFTYPE, tbeg0, tend0, const_cast<char*>("COLZ"));
  WTS[0]->hist2D->GetYaxis()->SetRangeUser(PLOT_MINFREQ, PLOT_MAXFREQ);
  double hist2D_MAX = WTS[0]->hist2D->GetMaximum();
#endif

  // define time array of the PC reconstructed waveform 
  wavearray<double> rec(sig.size());    
  rec.rate(sig.rate());
  rec=0;
  double Erec=0.;

  // define netcluster for monster display
  netcluster* pwc = new netcluster;
  netpixel* pix = new netpixel;
  pixdata pdata;
  float crate = sig.rate();
  float ctime;
  float cfreq;
  std::vector<int> list;                            // cluster list
  std::vector<int> vtof(1);                         // time configuration array
  std::vector<int> vtmp;                            // sky index array
  std::vector<float> varea;                         // sky error regions array
  std::vector<float> vpmap;                         // sky pixel map array
  clusterdata cd;                                   // dummy cluster data

  // start PC extraction
  cout << "Start PC extraction ..." << endl;
  double Ethr = 2*ACORE*ACORE;		// selection pixel threshold
  double precision=1;
  int npix=0;
  int nn,mm;
  int nWP[NWP+1]; for(int l=0;l<=NWP;l++) nWP[l]=0;
  int nPC=0;
  for(int k=0;k<PCA_MAX;k++) { 

    //if(fabs(Esig-Erec)/Esig<PCA_PREC) break;
    if(precision<PCA_PREC) break;

    // find pixel with max energy -> principal component
    double amax;	// max pixel amplitude
    int    nmax;	// time index @ amax
    int    mmax;	// freq index @ amax
    int    rmax;	// resolution @ amax
    int isWP=GetPrincipalComponent(amax, nmax, mmax, rmax); 
    if(amax*amax<Ethr) break;	// break loop if the PC energy is below the threshold
    double factor = 1.;
#ifdef LOAD_XTALK_CATALOG
    factor = GetPrincipalFactor(nmax, mmax, rmax);
#endif
    nPC++; nWP[isWP]++;

    cout << endl << "--------------------------------------------------------" << endl;
    cout <<         "Principal Component # " << k << endl;
    cout <<         "--------------------------------------------------------" << endl;

    if(isWP) cout << "---> WAVELET PACKET : " << isWP << endl;
    else     cout << "---> wavelet pixel" << endl;
    cout << " amax : " << amax << " nmax : " << nmax 
         << " mmax : " << mmax << " rmax : " << rmax << endl;

    double a00[3];
    double a90[3];
    a00[0] = wsig[rmax].getSample(nmax,mmax);
    a90[0] = wsig[rmax].getSample(nmax,-mmax);

    if(isWP) {	// wavelet packet : select 2 more pixels

      if(isWP==1) {nn=nmax+1; mm=mmax+1;}
      if(isWP==2) {nn=nmax-1; mm=mmax;  }
      if(isWP==3) {nn=nmax;   mm=mmax+1;}
      if(isWP==4) {nn=nmax-1; mm=mmax+1;}

      a00[1] = wsig[rmax].getSample(nn,mm);
      a90[1] = wsig[rmax].getSample(nn,-mm);

      if(isWP==1) {nn=nmax-1; mm=mmax-1;}
      if(isWP==2) {nn=nmax+1; mm=mmax;  }
      if(isWP==3) {nn=nmax;   mm=mmax-1;}
      if(isWP==4) {nn=nmax+1; mm=mmax-1;}

      a00[2] = wsig[rmax].getSample(nn,mm);
      a90[2] = wsig[rmax].getSample(nn,-mm);
    }

    double ishift=0;
#ifdef APPLY_TSHIFT
    // get time delayed amplitutes A00,A90
    WDM<double>* pWDM = (WDM<double>*)wsig[rmax].pWavelet;
    pWDM->setTDFilter(WDM_TDSIZE);           // computes TD filters
    int nsmp = (1<<rmax+PCA_IRES)/2;
    double Emax=a00[0]*a00[0]+a90[0]*a90[0];
    if(isWP) {	// wavelet packet
      Emax+=(a00[1]*a00[1]+a90[1]*a90[1]);
      Emax+=(a00[2]*a00[2]+a90[2]*a90[2]);
    }
    double A00[3];
    double A90[3];
    A00[0]=a00[0]; A90[0]=a90[0];
    if(isWP) {
      A00[1]=a00[1]; A90[1]=a90[1];
      A00[2]=a00[2]; A90[2]=a90[2];
    }
    double d00[3];
    double d90[3];
    double EMAX=0;
    for(int i=-nsmp;i<=nsmp;i++) {
      d00[0] = (double)pWDM->getPixelAmplitude(mmax, nmax,  i, false);  // delay a00 by i samples
      d90[0] = (double)pWDM->getPixelAmplitude(mmax, nmax,  i, true);   // delay a90 by i samples
      EMAX = d00[0]*d00[0]+d90[0]*d90[0];
      if(EMAX>Emax) {Emax=EMAX;ishift=i;A00[0]=d00[0];A90[0]=d90[0];}

      if(isWP) {	// wavelet packet : select 2 more pixels

        if(isWP==1) {nn=nmax+1; mm=mmax+1;}
        if(isWP==2) {nn=nmax-1; mm=mmax;  }
        if(isWP==3) {nn=nmax;   mm=mmax+1;}
        if(isWP==4) {nn=nmax-1; mm=mmax+1;}

        d00[1] = (double)pWDM->getPixelAmplitude(mm, nn,  i, false);  // delay a00 by i samples
        d90[1] = (double)pWDM->getPixelAmplitude(mm, nn,  i, true);   // delay a90 by i samples

        if(isWP==1) {nn=nmax-1; mm=mmax-1;}
        if(isWP==2) {nn=nmax+1; mm=mmax;  }
        if(isWP==3) {nn=nmax;   mm=mmax-1;}
        if(isWP==4) {nn=nmax+1; mm=mmax-1;}

        d00[2] = (double)pWDM->getPixelAmplitude(mm, nn,  i, false);  // delay a00 by i samples
        d90[2] = (double)pWDM->getPixelAmplitude(mm, nn,  i, true);   // delay a90 by i samples

	// parity
        int sign;
        if(isWP==1) sign = mmax%2 ?  1 :-1;
        if(isWP==2) sign = mmax%2 ? -1 : 1;
        if(isWP==3) sign = mmax%2 ? -1 : 1;
        if(isWP==4) sign = mmax%2 ? -1 : 1;

        // get D00,D90 with respect to z,Z 
        // x = a90[1]*U +/- a00[0]*v + a90[2]*W
        // X = a00[1]*u -/+ a90[0]*V + a00[2]*w
        double D00 = 0;
        double D90 = 0;
        if(sign==1) {
          D00 =  a00[0]+a90[1]+a90[2]; 	// D00 = (x,z) 
          D90 = -a90[0]+a00[1]+a00[2]; 	// D00 = (x,z)
        } else {
          D00 = -a00[0]+a90[1]+a90[2]; 	// D00 = (x,z)
          D90 =  a90[0]+a00[1]+a00[2]; 	// D00 = (x,z)
        }

        // D00 = (x,z)/(z,z)
        // D90 = (X,Z)/(Z,Z)
        // (z,z) = (Z,Z) = 3+4*WP_XTALK
        D00/=(3+4*wp_xtalk[isWP-1]);
        D90/=(3+4*wp_xtalk[isWP-1]);

        // z = z/sqrt(3+4*WP_XTALK)
        // Z = Z/sqrt(3+4*WP_XTALK)
        // (z,z) = (Z,Z) = 1
        D00*=sqrt(3+4*wp_xtalk[isWP-1]);
        D90*=sqrt(3+4*wp_xtalk[isWP-1]);

        EMAX = (D00*D00+D90*D90);

        if(EMAX>Emax) {
          Emax=EMAX;ishift=i;
          A00[0]=d00[0];A90[0]=d90[0];
          A00[1]=d00[1];A90[1]=d90[1];
          A00[2]=d00[2];A90[2]=d90[2];
        }
      }
    }
    cout << 0      << " a00[0] : " << a00[0] << " a90[0] : " << a90[0] << endl;
    cout << ishift << " A00[0] : " << A00[0] << " A90[0] : " << A90[0] << endl;
    cout << "tshift : " << 1000*ishift/wsig[rmax].rate() << " msec" << endl;
    cout << "tmax   : " << 1000*nsmp/wsig[rmax].rate() << " msec" << endl;
    a00[0]=A00[0];a90[0]=A90[0];
    if(isWP) {
      a00[1]=A00[1];a90[1]=A90[1];
      a00[2]=A00[2];a90[2]=A90[2];
    }
#endif

    // get max pixels time representation

    int layers = wsig[rmax].maxLayer()+1;   // numbers of frequency bins (first & last bins have df/2)
    int slices = wsig[rmax].sizeZero();     // number of time bins
    //cout << "layers " << layers << " slices " << slices << endl;

    int itime = nmax*layers + mmax;

    wavearray<double> w00[3];   //time series container
    wavearray<double> w90[3];   //time series container
    int j00[3];
    int j90[3];

    j00[0] = wdm[rmax]->getBaseWave(itime,w00[0],false)-ishift;
    j90[0] = wdm[rmax]->getBaseWave(itime,w90[0],true)-ishift;

    wavearray<double> x(sig.size());   // PC time series container
    x.rate(sig.rate());
    x=0;

    double D00 = 0;
    double D90 = 0;
    if(isWP) {	// wavelet packet : select 2 more pixels

      if(isWP==1) {nn=nmax+1; mm=mmax+1;}
      if(isWP==2) {nn=nmax-1; mm=mmax;  }
      if(isWP==3) {nn=nmax;   mm=mmax+1;}
      if(isWP==4) {nn=nmax-1; mm=mmax+1;}

      itime = nn*layers + mm;
      j00[1] = wdm[rmax]->getBaseWave(itime,w00[1],false)-ishift;
      j90[1] = wdm[rmax]->getBaseWave(itime,w90[1],true)-ishift;

      if(isWP==1) {nn=nmax-1; mm=mmax-1;}
      if(isWP==2) {nn=nmax+1; mm=mmax;  }
      if(isWP==3) {nn=nmax;   mm=mmax-1;}
      if(isWP==4) {nn=nmax+1; mm=mmax-1;}

      itime = nn*layers + mm;
      j00[2] = wdm[rmax]->getBaseWave(itime,w00[2],false)-ishift;
      j90[2] = wdm[rmax]->getBaseWave(itime,w90[2],true)-ishift;

      // parity	
      // sign= 1 -> z=U+v+W  , Z=u-V+w
      // sign=-1 -> z=U-v+W  , Z=u+V+w
      // v=w00[0], u=w00[1], w=w00[2]
      // V=w90[0], U=w90[1], W=w90[2]
      // (z,z) = (Z,Z) = 3+4*WP_XTALK
      // (z,Z) = 0
      int sign;
      if(isWP==1) sign = mmax%2 ?  1 :-1;
      if(isWP==2) sign = mmax%2 ? -1 : 1;
      if(isWP==3) sign = mmax%2 ? -1 : 1;
      if(isWP==4) sign = mmax%2 ? -1 : 1;
      cout << " parity : " << sign << endl;

      // get D00,D90 with respect to z,Z 
      // x = a90[1]*U +/- a00[0]*v + a90[2]*W
      // X = a00[1]*u -/+ a90[0]*V + a00[2]*w
      // (x,X)=0
      if(sign==1) {
        D00 =  a00[0]+a90[1]+a90[2]; 	// D00 = (x,z) 
        D90 = -a90[0]+a00[1]+a00[2]; 	// D00 = (x,z)
      } else {
        D00 = -a00[0]+a90[1]+a90[2]; 	// D00 = (x,z)
        D90 =  a90[0]+a00[1]+a00[2]; 	// D00 = (x,z)
      }
      D00/=(3+4*wp_xtalk[isWP-1]);	// D00 = (x,z)/(z,z)
      D90/=(3+4*wp_xtalk[isWP-1]);	// D90 = (x,Z)/(Z,Z)

      // D00*z
      // sign= 1 -> z=U+v+W 
      // sign=-1 -> z=U-v+W 
      for(int i=0; i<w90[1].size(); i++){	// U=w90[1]
        if(j90[1]+i<0 || j90[1]+i>=x.size()) continue;
        x.data[j90[1]+i] +=      D00*w90[1][i];
      }
      for(int i=0; i<w00[0].size(); i++){	// v=w00[0]
        if(j00[0]+i<0 || j00[0]+i>=x.size()) continue;
        x.data[j00[0]+i] += sign*D00*w00[0][i];
      }
      for(int i=0; i<w90[2].size(); i++){	// v=w00[0]
        if(j90[2]+i<0 || j90[2]+i>=x.size()) continue;
        x.data[j90[2]+i] +=      D00*w90[2][i];
      }

      // D90*Z
      // sign= 1 -> Z=u-V+w
      // sign=-1 -> Z=u+V+w
      for(int i=0; i<w00[1].size(); i++){	// u=w00[1]
        if(j00[1]+i<0 || j00[1]+i>=x.size()) continue;
        x.data[j00[1]+i] +=      D90*w00[1][i];
      }
      for(int i=0; i<w90[0].size(); i++){	// V=w90[0]
        if(j90[0]+i<0 || j90[0]+i>=x.size()) continue;
        x.data[j90[0]+i] -= sign*D90*w90[0][i];
      }
      for(int i=0; i<w00[2].size(); i++){	// w=w00[2]
        if(j00[2]+i<0 || j00[2]+i>=x.size()) continue;
        x.data[j00[2]+i] +=      D90*w00[2][i];
      }

    } else {	// single pixel

      for(int i=0; i<w00[0].size(); i++){
        if(j00[0]+i<0 || j00[0]+i>=x.size()) continue;
        x.data[j00[0]+i] += factor*a00[0]*w00[0][i];
      }
      for(int i=0; i<w90[0].size(); i++){
        if(j90[0]+i<0 || j90[0]+i>=x.size()) continue;
        x.data[j90[0]+i] += factor*a90[0]*w90[0][i];
      }
      
      D00 = a00[0];
      D90 = a90[0];
    }

    // fill pixel infos for monster display
    npix++;
    pix->core = 1;
    pix->clusterID = 0;                       // setup new ID
    pix->rate = x.rate()/(1<<(rmax+PCA_IRES));
    pix->layers = layers;
    pix->time = itime;
    pix->frequency = mmax;
    pix->data.clear();
    pdata.asnr = D00;
    pdata.a_90 = D90;
    pix->data.push_back(pdata);
    pwc->pList.push_back(*pix);               // update pList and counter

    // PC waveform
    for(int i=0; i<x.size(); i++) rec.data[i] += x[i];
    Erec=0.; for(int i=0; i<rec.size(); i++) Erec+=rec[i]*rec[i];

    // subtract max pixels from TF maps
    WSeries<double> wx;   		// TF map container
    for(int i=0;i<PCA_NRES;i++) {
      wx.Forward(x, *wdm[i]);            // apply the WDM to the max pixel time series
      wsig[i]-=wx;			// subtract max pixels from TF maps
    }

    // compute precision = Eresidual/Esignal
    wavearray<double> y = rec;
    y-=sig; y*=y;
    wavearray<double> y2 = sig;
    y2*=y2;
    precision = y.mean()/y2.mean();
    cout << " Erec " << Erec << " Esig " << Esig 
         << " fabs(Esig-Erec)/Esig : " << fabs(Esig-Erec)/Esig  
         << " precision " << precision << endl;
#ifdef OUTPUT_FILE
    out << Erec << " " << Esig << " " << fabs(Esig-Erec)/Esig << " " << precision << endl;
#endif
  }
#ifdef OUTPUT_FILE
  out.close();
#endif

  cout << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << "number of               SW/nPC used in the PCA : " << nWP[0] << "/" << nPC << endl; 
  cout << "number of      diagonal_WP/nPC used in the PCA : " << nWP[1] << "/" << nPC << endl; 
  cout << "number of    horizontal_WP/nPC used in the PCA : " << nWP[2] << "/" << nPC << endl; 
  cout << "number of      vertical_WP/nPC used in the PCA : " << nWP[3] << "/" << nPC << endl; 
  cout << "number of back-diagonal_WP/nPC used in the PCA : " << nWP[4] << "/" << nPC << endl; 
#ifdef OUTPUT_FILE
  cout << "output file : " << OUTPUT_FILE << endl;
#endif
  cout << "-----------------------------------------------------" << endl;
  cout << endl;

  // display tf (residuals)
#ifdef PLOT_TFMAP		
  double dt1 = 1./wsig[PLOT_TFRES].rate();
  double tbeg1 = WDM_SCRATCH-2*dt1;
  double tend1 = dt1*wsig[PLOT_TFRES].size()/2-WDM_SCRATCH;
  WTS[1] = new watplot(const_cast<char*>("wts1"));
  WTS[1]->plot(&wsig[PLOT_TFRES], PLOT_TFTYPE, tbeg1, tend1, const_cast<char*>("COLZ"));
  WTS[1]->hist2D->GetYaxis()->SetRangeUser(PLOT_MINFREQ, PLOT_MAXFREQ);
  WTS[1]->hist2D->GetZaxis()->SetRangeUser(0, hist2D_MAX);
#endif

#ifdef PLOT_MONSTER
  // monster display
  for(int i=0;i<npix;i++) list.push_back(i);
  // fill pixel infos for monster display
  pwc->rate = sig.rate();
  pwc->cList.push_back(list);
  pwc->cRate.push_back(crate);
  pwc->cTime.push_back(ctime);
  pwc->cFreq.push_back(cfreq);
  pwc->sArea.push_back(varea);                    // recreate sky error regions array
  pwc->p_Map.push_back(vpmap);                    // recreate sky pixel map array
  pwc->nTofF.push_back(vtof);                     // recreate time configuration array
  pwc->p_Ind.push_back(vtmp);                     // recreate sky index array
  pwc->cData.push_back(cd);             	  // recreate dummy cluster data

  WTS[2] = new watplot(const_cast<char*>("wts2"));
  WTS[2]->plot(pwc, 1, 1, 'L', 0, const_cast<char*>("COLZ"));
#endif

#ifdef PLOT_REC_VS_INJ
  // display rec vs sig waveform in time domain
  grec = new gwavearray<double>;
  *grec = rec;
  gsig = new gwavearray<double>;
  *gsig = sig;
  gsig->SetTimeRange(WDM_SCRATCH, gsig->size()/gsig->rate()-WDM_SCRATCH);
  gsig->Draw(GWAT_TIME,"ALPCUSTOM");
  gsig->Draw(grec,GWAT_TIME,"SAME",kRed);
#endif
}

int GetPrincipalComponent(double &amax, int &nmax, int &mmax, int &rmax) {
 
  rmax=0;
  nmax=0;
  mmax=0;
  amax=0;

  int Rmax=0;
  double Nmax=0;
  double Mmax=0;
  double Amax=0;

  double Ethr = 2*ACORE*ACORE;		// selection pixel threshold
  int isWP=0;
  int nn,mm;
  double a00[3];
  double a90[3];
  double A00 = 0;
  double A90 = 0;
 
  for(int i=0;i<PCA_NRES;i++) {

    int layers = wsig[i].maxLayer()+1;   // numbers of frequency bins (first & last bins have df/2)
    int slices = wsig[i].sizeZero();     // number of time bins
    //cout << "layers " << layers << " slices " << slices << endl;

    int iscratch = WDM_SCRATCH*(wsig[i].rate()/(1<<(i+PCA_IRES)));
    for(int n=iscratch;n<slices-iscratch;n++) {	// exclude scratch data
      for(int m=2;m<layers-2;m++) {
        a00[0] = wsig[i].getSample(n,m);
        a90[0] = wsig[i].getSample(n,-m);

        double E = (a00[0]*a00[0]+a90[0]*a90[0]);

        double A = sqrt(E);
        if(wp_mask[0]) {
          if(A>amax) {amax=A; nmax=n; mmax=m; rmax=i; isWP=0;}
          if(A>Amax) {Amax=A; Nmax=n; Mmax=m; Rmax=i; A00=a00[0]; A90=a90[0];}
        }

        for(int l=0;l<NWP;l++) {

          if(!wp_mask[l+1]) continue;
       
          if(l==0) {nn=n+1; mm=m+1;}	// diagonal WP
          if(l==1) {nn=n-1; mm=m;  }	// horizontal WP
          if(l==2) {nn=n;   mm=m+1;}	// vertical WP
          if(l==3) {nn=n-1; mm=m+1;}	// back-diagonal WP

          a00[1] = wsig[i].getSample(nn,mm);
          a90[1] = wsig[i].getSample(nn,-mm);
       
          if(l==0) {nn=n-1; mm=m-1;}	// diagonal WP
          if(l==1) {nn=n+1; mm=m;  }	// horizontal WP
          if(l==2) {nn=n;   mm=m-1;}	// vertical WP
          if(l==3) {nn=n+1; mm=m-1;}	// back-diagonal WP

          a00[2] = wsig[i].getSample(nn,mm);
          a90[2] = wsig[i].getSample(nn,-mm);

          // parity
          int sign;
          if(l==0) sign = m%2 ?  1 :-1;
          if(l==1) sign = m%2 ? -1 : 1;
          if(l==2) sign = m%2 ? -1 : 1;
          if(l==3) sign = m%2 ? -1 : 1;

          // sign= 1 -> z=U+v+W  , Z=u-V+w
          // sign=-1 -> z=U-v+W  , Z=u+V+w
          // v=w00[0], u=w00[1], w=w00[2]
          // V=w90[0], U=w90[1], W=w90[2]
          // (v,U) = (v,W) = -(V,u) = -(V,w) = WP_XTALK
          // (z,Z) = 0

          // get D00,D90 with respect to z,Z 
          // x = a90[1]*U +/- a00[0]*v + a90[2]*W
          // X = a00[1]*u -/+ a90[0]*V + a00[2]*w
          double D00 = 0;
          double D90 = 0;
          if(sign==1) {
            D00 =  a00[0]+a90[1]+a90[2]; 	// D00 = (x,z) 
            D90 = -a90[0]+a00[1]+a00[2]; 	// D00 = (x,z)
          } else {
            D00 = -a00[0]+a90[1]+a90[2]; 	// D00 = (x,z)
            D90 =  a90[0]+a00[1]+a00[2]; 	// D00 = (x,z)
          }
          // (z,z) = (Z,Z) = 3+4*WP_XTALK
          D00/=(3+4*wp_xtalk[l]);
          D90/=(3+4*wp_xtalk[l]);

          // get D00,D90 with respect to the normalized z,Z 
          // z = z/sqrt(3+4*WP_XTALK)
          // Z = Z/sqrt(3+4*WP_XTALK)
          // (z,z) = (Z,Z) = 1
          D00*=sqrt(3+4*wp_xtalk[l]);
          D90*=sqrt(3+4*wp_xtalk[l]);

          double E = (D00*D00+D90*D90);

          A = sqrt(E);

          if(A>amax) {amax=A; nmax=n; mmax=m; rmax=i; isWP=l+1;}
        }
        // set isWP=0 if the difference between 1 pixels and WP is below the noise threshold
        // 2*Ethr is the noise of the 2 WP extra pixels
//        if(fabs(amax*amax-Amax*Amax)<2*Ethr) {amax=Amax; nmax=Nmax; mmax=Mmax; rmax=Rmax; isWP=0;}
      }
    }
  }
  return isWP;
}

#ifdef LOAD_XTALK_CATALOG
void CheckXTalkCatalog(int l_low, int l_high) {

  for(int level=l_high; level>=l_low; level--) {
    int layers = level>0 ? (1<<level)+1 : 0;
    cout << "layers : " << layers << endl;

    // v,V
    int _n = 100;
    int _m = 9;		// parity = m%2 ?  even : odd;	m=9 -> parity=even
    // w,W
    int _nn=_n+1; 
    int _mm=_m+1;

/*
    for(int ii=-1;ii<=1;ii++) {
      for(int jj=-1;jj<=1;jj++) {
        // CC[0]  CC[2]
        // CC[1]  CC[3] 
        _nn=_n+ii; 
        _mm=_m+jj;
        struct xtalk xt = wdmMRA.getXTalk(layers, _n*layers+_m, layers, _nn*layers+_mm);
        cout << ii << " " << jj << " (v,w) : " << xt.CC[0] << " (v,W) : " << xt.CC[1] 
                                << " (V,w) : " << xt.CC[2] << " (V,W) : " << xt.CC[3] << endl;
      }
    }
*/

    struct xtalk xt = wdmMRA.getXTalk(layers, _n*layers+_m, layers, _nn*layers+_mm);
    cout << "xtalk catalog    : " << " (v,w) : " << xt.CC[0] << " (v,W) : " << xt.CC[1] 
                                  << " (V,w) : " << xt.CC[2] << " (V,W) : " << xt.CC[3] << endl;

    int j00[2];
    int j90[2];
    wavearray<double> w00[2];   //time series container
    wavearray<double> w90[2];   //time series container

    j00[0] = wdm[level-l_low]->getBaseWave( _n*layers+_m, w00[0],false);	// v
    j90[0] = wdm[level-l_low]->getBaseWave( _n*layers+_m, w90[0],true);	// V
    j00[1] = wdm[level-l_low]->getBaseWave(_nn*layers+_mm,w00[1],false);	// w
    j90[1] = wdm[level-l_low]->getBaseWave(_nn*layers+_mm,w90[1],true);	// W

    wavearray<double> W00[2];   //time series container
    wavearray<double> W90[2];   //time series container
    for(int k=0;k<2;k++) {
      W00[k]=sig; W00[k]=0;      
      for(int i=0;i<w00[k].size();i++) {
        if(i+j00[k]<0 || i+j00[k]>W00[k].size()) continue;
        W00[k][i+j00[k]] = w00[k][i];
      }
      W90[k]=sig; W90[k]=0;      
      for(int i=0;i<w90[k].size();i++) {
        if(i+j90[k]<0 || i+j90[k]>W90[k].size()) continue;
        W90[k][i+j90[k]] = w90[k][i];
      }
    }
 
    double vw=0.;
    double vW=0.;
    double Vw=0.;
    double VW=0.;
    for(int i=0;i<W90[0].size();i++) {
      vw+=W00[0][i]*W00[1][i];		// (v,w)
      vW+=W00[0][i]*W90[1][i];		// (v,W)
      Vw+=W90[0][i]*W00[1][i];		// (V,w)
      VW+=W90[0][i]*W90[1][i];		// (V,W)
    }
    // NOTE : for more infos see InitXTALK function
    cout << "xtalk base waves : " << " (v,w) : " << vw << " (v,W) : " << vW
                                  << " (V,w) : " << Vw << " (V,W) : " << VW << endl << endl;
    cout << "wp_xtalk         : " << " (v,w) : " << 0            << " (v,W) : " << wp_xtalk[0]
                                  << " (V,w) : " << -wp_xtalk[0] << " (V,W) : " << 0 << endl << endl;
  }
}
#endif

double GetPrincipalFactor(int nmax, int mmax, int rmax) {
 
  double Ethr = 2*ACORE*ACORE;		// selection pixel threshold
  int nn,mm;
  vector<double> vb;
  vector<double> vB;
  vector<int> vN;
  vector<int> vM;
  vector<int> vL;

  struct xtalk xt;
 
  for(int i=0;i<PCA_NRES;i++) {

    int layers = wsig[i].maxLayer()+1;   // numbers of frequency bins (first & last bins have df/2)
    int slices = wsig[i].sizeZero();     // number of time bins
    //cout << "layers " << layers << " slices " << slices << endl;

    int iscratch = WDM_SCRATCH*(wsig[i].rate()/(1<<(i+PCA_IRES)));
    for(int n=iscratch;n<slices-iscratch;n++) {	// exclude scratch data
      for(int m=2;m<layers-2;m++) {
        double a00 = wsig[i].getSample(n,m);
        double a90 = wsig[i].getSample(n,-m);

        double E = (a00*a00+a90*a90);

        if(E<Ethr) continue;

        //cout << "r:n:m " << i << ":" << n << ":" << m << " E : " << E << endl;
        //cout << "lmax:nmax:mmax " << lmax << ":" << nmax << ":" << mmax << endl; 
        
        vb.push_back(a00); 
        vB.push_back(a90); 
        vN.push_back(n); 
        vM.push_back(m); 
        vL.push_back(layers); 
      }
    }
  }

  // PC is X = a*Y00 + A*Y90
  // x - q*X	 
  // sum_i((x,Yi)-q*(X,Yi))^2 = sum_i(bi-q*(X,Yi))^2 + sum_i(Bi-q*(X,Yi))^2
  // q = (sum_i(bi*(X,Yi)) + sum_i(Bi*(X,Yi))) / (sum_i((X,Yi)^2) + sum_i((X,Yi)^2))

  int lmax = wsig[rmax].maxLayer()+1;

  double a = wsig[rmax].getSample(nmax,mmax);
  double A = wsig[rmax].getSample(nmax,-mmax);

  double num=0;
  double den=0;
  for(int k=0;k<vb.size();k++) {
    double b = vb[k];
    double B = vB[k];
    int n = vN[k];
    int m = vM[k];
    int l = vL[k];
    xt = wdmMRA.getXTalk(lmax, nmax*lmax+mmax, l, n*l+m);
    if(xt.CC[0]>2) continue;
    //cout << "r:n:m " << l << ":" << n << ":" << m << " b : " << b << " B : " << B << endl;
    //cout << "lmax:nmax:mmax " << lmax << ":" << nmax << ":" << mmax << endl; 
    //cout << k << " (u,w) : " << xt.CC[0] << " (u,W) : " << xt.CC[1] 
    //          << " (U,w) : " << xt.CC[2] << " (U,W) : " << xt.CC[3] << endl;

    num+=b*(xt.CC[0]*a+xt.CC[2]*A)+B*(xt.CC[1]*a+xt.CC[3]*A);
    den+=pow(xt.CC[0]*a+xt.CC[2]*A,2)+pow(xt.CC[1]*a+xt.CC[3]*A,2);
  }

  double q = num/den;

  cout << "q = " << q << endl;

  return q;
}

void InitXTALK() {

// -------------------------------------------------------------------------
// XTALK DEFINITIONS
// NOTE : All xtalk values are computed with 'v,V' parity = even (see below)
// -------------------------------------------------------------------------
//
// WAVELET PACKET BASE
// parity = m%2 ?  even : odd; (n,m) are the 'v,V' time,frequency index
// parity = odd  -> z=U+v+W  , Z=u-V+w
// parity = even -> z=U-v+W  , Z=u+V+w
// (z,z) = (Z,Z) = 3+4*WP_XTALK
// (z,Z) = 0
//
// diagonal WP
// - - w         - - W 
// - v -         - V - 
// u - -         U - - 
//                                                                                                              
// horizontal WP      
// - - -         - - -
// u v w         U V W 
// - - -         - - -
//                     
// vertical WP         
// - w -         - W - 
// - v -         - V -
// - u -         - U -  
//                                                                                                              
// back-diagonal WP    
// w - -         W - -
// - v -         - V - 
// - - u         - - U 
//                     
// wp_xtalk = (v,W) = (W,v) = -(V,w) = -(w,V)  
// wp_xtalk = (v,U) = (U,v) = -(V,u) = -(u,V) 
//                                           

#if   WDM_NU == 2

  wp_xtalk[0] =  0.207345;      // diagonal WP          (v,W)=(n:m, n+1:m+1) (v,U)=(n:m, n-1:m-1)
  wp_xtalk[1] =  0.561945;      // horizontal WP        (v,U)=(n:m, n-1:m  ) (v,W)=(n:m, n+1:m  )
  wp_xtalk[2] =  0.243038;      // vertical WP          (v,W)=(n:m, n  :m+1) (v,U)=(n:m, n  :m-1)
  wp_xtalk[3] = -0.207345;      // back-diagonal WP     (v,W)=(n:m, n-1:m+1) (v,U)=(n:m, n+1:m-1)

#elif WDM_NU == 4

  wp_xtalk[0] =  0.162725;      // diagonal WP          (v,W)=(n:m, n+1:m+1) (v,U)=(n:m, n-1:m-1)
  wp_xtalk[1] =  0.597494;      // horizontal WP        (v,U)=(n:m, n-1:m  ) (v,W)=(n:m, n+1:m  )
  wp_xtalk[2] =  0.178856;      // vertical WP          (v,W)=(n:m, n  :m+1) (v,U)=(n:m, n  :m-1)
  wp_xtalk[3] = -0.162725;      // back-diagonal WP     (v,W)=(n:m, n-1:m+1) (v,U)=(n:m, n+1:m-1)

#elif WDM_NU == 6

  wp_xtalk[0] =  0.138369;      // diagonal WP          (v,W)=(n:m, n+1:m+1) (v,U)=(n:m, n-1:m-1)
  wp_xtalk[1] =  0.610111;      // horizontal WP        (v,U)=(n:m, n-1:m  ) (v,W)=(n:m, n+1:m  )
  wp_xtalk[2] =  0.148011;      // vertical WP          (v,W)=(n:m, n  :m+1) (v,U)=(n:m, n  :m-1)
  wp_xtalk[3] = -0.138369;      // back-diagonal WP     (v,W)=(n:m, n-1:m+1) (v,U)=(n:m, n+1:m-1)

#else

  bool wpCheck=false;
  for(int l=1;l<=NWP;l++) if(wp_mask[l]) wpCheck=true;
  if(wpCheck) {
    if(WDM_NU!=2) {
      cout << endl << "Error : wp_xtalk coefficients are defined only for WDM_NU=2,4,6 !!!" << endl << endl;
      exit(1);
    }
  }

#endif

}

