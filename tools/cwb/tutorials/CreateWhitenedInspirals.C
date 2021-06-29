//#define GENERATE_XML

#define ADV_LIGO_PSD "detectors_PSD/LIGO_zero_det_HP.txt"
#define SRATE     16384

#define FREQ_MIN  32
#define FREQ_MAX  2048

#define OFILE	"inspirals.root"

#define nWAVEFORM 1000

#define WHITENING

wavearray<double> ReadDetectorPSD(TString fName);
double GetBoundaries(wavearray<double> x, double P, double& bT, double& eT);

void 
CreateWhitenedInspirals() {
  //
  // generate LAL whitened inspiral waveforms
  // Author : Gabriele Vedovato

  wavearray<double> psd = ReadDetectorPSD(ADV_LIGO_PSD);
  cout << "PSD : " << " rate : " << psd.rate() << " size : " << psd.size() << endl;

  CWB::mdc MDC; 

#ifdef GENERATE_XML
  // ---------------------------------
  // set inspiral parms
  // ---------------------------------
  TString inspOptions="";
  inspOptions = "--time-step 60.0 --time-interval 3 ";
  inspOptions+= "--l-distr random ";
  inspOptions+= "--gps-start-time 931072130 --gps-end-time 933491330 ";
  inspOptions+= "--d-distr volume --m-distr totalMassRatio --i-distr uniform ";
  inspOptions+= "--f-lower 10.000000 ";
  inspOptions+= "--min-mass1 25.000000 --max-mass1 225.000000 ";
  inspOptions+= "--min-mass2 25.000000 --max-mass2 225.000000 ";
  inspOptions+= "--min-mtotal 50.000000 --max-mtotal 250.000000 ";
  inspOptions+= "--min-mratio 0.25 --max-mratio 1.000000 ";
  inspOptions+= "--min-distance 1000000.0 --max-distance 1500000.0 ";
  inspOptions+= "--approximant EOBNRv2pseudoFourPN --disable-spin ";
  inspOptions+= "--taper-injection start --seed 123456789 ";
  inspOptions+= "--output inspirals.xml ";		// set output xml file

  // set and write xml file
  MDC.SetInspiral("EOBNRv2",inspOptions);
#endif

  // set inspiral using xml file (speedup GetInspiral method)
  TString inspOptions="--xml inspirals.xml";
  MDC.SetInspiral("EOBNRv2",inspOptions);

  int gps_start_time = 931072131;
  int gps_end_time   = gps_start_time+100;

  wavearray<double> w(8*SRATE);
  w.start(0);
  w.rate(SRATE);
  wavearray<double> hp;

  TFile ofile(OFILE,"RECREATE");

  // Get the first waveform hp,hx components starting from gps = 931072130
  for(int n=0;n<nWAVEFORM;n++) {
    cout << n << " gps_start_time " << gps_start_time << " gps_end_time " << gps_end_time << endl;
    hp = MDC.GetInspiral("hp",gps_start_time,gps_end_time);
    cout << "size : " << hp.size() << " rate : " << hp.rate() << " start : " << (int)hp.start() << endl;

    if(hp.size()>w.size() || (hp.size()==0)) {
      gps_start_time += 50;
      gps_end_time    = gps_start_time+100;
      continue;
    }

    for(int i=0;i<hp.size();i++)        w[i]=hp[i];
    for(int i=hp.size();i<w.size();i++) w[i]=0;

#ifdef WHITENING
    // waveform whitening 
    w.FFTW(1);
    double df = w.rate()/w.size();
    for(int i=0;i<w.size();i++) {
      double f=i*df;
      if(f<FREQ_MIN || f>FREQ_MAX) {
        w[i]=0.;
      } else {
        w[i]/=psd[int(i/2)];
      }
    }
    w.FFTW(-1);
    w*=1./w.rms();	// normalization
#endif

    double bT,eT;
    GetBoundaries(w, 0.99, bT, eT);
    cout << w.start() << " " << bT << " " << eT << " " << eT-bT << endl;

    double dt = 1./w.rate();
    int size = (eT-bT)/dt;
    int os = bT/dt;
    wavearray<double> W(size);
    W.start(0);
    W.rate(w.rate());
    for(int i=0;i<size;i++) W[i]=w[i+os];

    // save whitened waveform to root file
    W.Write("Inspiral");

//    MDC.Draw(W,MDC_TIME);break;
//    MDC.Draw(W,MDC_FFT);	// draw hp in frequency domain
//    MDC.Draw(W,MDC_TF);        // draw hp in time-frequency domain;

    //gps_start_time = hp.start()+ hp.size()/hp.rate();
    gps_start_time += 50;
    gps_end_time    = gps_start_time+100;
  }

  ofile.Close();
}

wavearray<double>
ReadDetectorPSD(TString fName) {

  #define FREQ_RES       0.125	// Hz

  CWB::Toolbox TB;

  // read PSD
  ifstream in;
  in.open(fName.Data(),ios::in);
  if (!in.good()) {cout << "ReadDetectorPSD - Error Opening File : " << fName.Data() << endl;exit(1);}
  cout << "ReadDetectorPSD - Read File : " << fName.Data() << endl;

  int size=0;
  char str[1024];
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') size++;
  }
  //cout << "size " << size << endl;
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);

  wavearray<double> ifr(size);
  wavearray<double> ish(size);

  int cnt=0;
  while (1) {
    in >> ifr.data[cnt] >> ish.data[cnt];
    if (!in.good()) break;
    if(ish.data[cnt]<=0)
      {cout << "ReadDetectorPSD - input sensitivity file : " << fName.Data()
            << " contains zero at frequency : " << ifr.data[cnt] << " Hz " << endl;exit(1);}
    cnt++;
  }
  in.close();

  // convert frequency sample
  size=int((SRATE/2)/FREQ_RES);
  wavearray<double> ofr(size);
  wavearray<double> osh(size);

  for(int i=0;i<(int)ofr.size();i++) ofr[i]=i*FREQ_RES;
  TB.convertSampleRate(ifr,ish,ofr,osh);

  osh.rate(SRATE/2);
  osh*=1./sqrt(osh.size()/(SRATE/2));  // normalization

  return osh;
}

double
GetBoundaries(wavearray<double> x, double P, double& bT, double& eT) {

  if(P<0) P=0;
  if(P>1) P=1;

  int N = x.size();

  double E = 0;                                                 // signal energy
  double avr = 0;                                               // average
  for(int i=0;i<N;i++) {avr+=i*x[i]*x[i]; E+=x[i]*x[i];}
  int M=int(avr/E);                                             // central index

  // search range which contains percentage P of the total energy E
  int jB=0;
  int jE=N-1;
  double a,b;
  double sum = ((M>=0)&&(M<N)) ? x[M]*x[M] : 0.;
  for(int j=1; j<N; j++) {
    a = ((M-j>=0)&&(M-j<N)) ? x[M-j] : 0.;
    b = ((M+j>=0)&&(M+j<N)) ? x[M+j] : 0.;
    if(a) jB=M-j;
    if(b) jE=M+j;
    sum += a*a+b*b;
    if(sum/E > P) break;
  }

  bT = x.start()+jB/x.rate();
  eT = x.start()+jE/x.rate();

  return eT-bT;
}


