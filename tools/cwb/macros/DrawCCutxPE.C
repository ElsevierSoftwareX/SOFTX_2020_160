/*
# Copyright (C) 2019 Gabriele Vedovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

// PURPOSE
// This is a macro used to visualize the effect of the "ccut" time-frequency cut implemented in the macro cwb_pereport.C
// The "ccut" time-frequency cut is used for inspiral higher modes investigations
// NOTE: most of the functions are extracted from the macro cwb_pereport.C

#define PE_MAX_EVENTS	100
#define N_IFO		3

std::vector<wavearray<double> > wREC[PE_MAX_EVENTS];    // reconstructed signal in the trials
std::vector<wavearray<double> > wINJ[PE_MAX_EVENTS];    // iwhitened injected posteriors in the trials

double wTCOA[PE_MAX_EVENTS];                    // tcoa   injected posteriors in the trials
float  wTHETA[PE_MAX_EVENTS];                   // theta  injected posteriors in the trials
float  wPHI[PE_MAX_EVENTS];                     // phi    injected posteriors in the trials
float  wM1[PE_MAX_EVENTS];                      // m1     injected posteriors in the trials
float  wM2[PE_MAX_EVENTS];                      // m1     injected posteriors in the trials
float  wS1[PE_MAX_EVENTS];                      // spin1        OffSource injected posteriors
float  wS2[PE_MAX_EVENTS];                      // spin2        OffSource injected posteriors
float  wRHO[PE_MAX_EVENTS];                     // rho    injected posteriors in the trials

int                gEVENTS;                     // number of read events
static TString     gIFO[3] = {"L1","H1","V1"};  // network ifos. WARNING: this setup is hard coded !!!  

double GetInjTcoa(double geocentric_tcoa, network* NET, TString ifo, double theta, double phi);

TString GetParameterFromInjLog(TString log, TString param);

wavearray<double> GetCCut(wavearray<double>* ts, double tcoa, double m1, double m2, double s1z, double s2z);
void GetCCutParms(TString options);
void DrawCCut(wavearray<double>* ts, double tcoa, double m1, double m2, double s1z, double s2z); 
void DrawCCut(WSeries<double> W, double tcoa, double m1, double m2, double s1z, double s2z);

struct uoptions {

  TString   ccut;               // Is used to select TF area using chirp TF cut. If="" then it is not applied.
                                // The selected TF region is defined by the parameters declared in the ccut string.
                                // format: wdm_fres:mc_lower_factor:mc_upper_factor:left_time:right_time
                                // default: 8:1.25:1.75:0.5:0.0
 
  // the following parameters are setup parsing the user ccut parameter. Are initialized with the default values
  int       ccut_wdm_fres;
  double    ccut_bchirp;
  double    ccut_uchirp;
  double    ccut_ltime;
  double    ccut_rtime;
};

uoptions     gOPT;         // global User Options
WDM<double>* gWDM;         // wdm used for chirp TF cut 


// chirp tf cut default parameters

#define CCUT_WDM_FRES           16      // the WDM freq resolution used to decompose the signal in TF domain
#define CCUT_BCHIRP             1.40    // mchirp factor used to select the upper region boundary Mc=Mc.ccut_bchirp
#define CCUT_UCHIRP             1.60    // mchirp factor used to select the upper region boundary Mc=Mc.ccut_uchirp
#define CCUT_LTIME              0.5     // left time from tcoa (sec) used to select the the time region: [tcoa-left_time:tcoa-right_time]
#define CCUT_RTIME              0.0     // right time from tcoa (sec) used to select the the time region: [tcoa-left_time:tcoa-right_time]


CWB::mdc *MDC;
TLine* ftcoa;		// tcoa time        function
TLine* flchirp;		// left   chirp cut function
TLine* frchirp;		// right  chirp cut function
TF1* fbchirp;		// bottom chirp cut function
TF1* fuchirp;		// upper  chirp cut function

// options

//#define SHOW_FFT
//#define SHOW_TIME
//#define SHOW_ENVELOPE

#define APPLY_CCUT
//#define APPLY_BANDPASS
//#define APPLY_SYNC

//#define SHOW_WDM_INJ
//#define SHOW_WDM_REC
#define SHOW_WDM_DIF

#define DEBUGGING	

void DrawCCutxPE(TString ifname, int wf_id=0, int ifo_id=0, TString options="16:1.40:1.60:0.5:0.0") {

  gOPT.ccut=options;

  GetCCutParms(gOPT.ccut);

  // read waveform from output root file
  TFile* froot = new TFile(ifname,"READ");
  if(froot==NULL) {
    cout << "Error : Failed to open file : " <<  ifname << endl;
    gSystem->Exit(1);
  }
  TTree* itree = (TTree*)froot->Get("waveburst");
  if(itree==NULL) {
    cout << "Error : Failed to open tree waveburst from file : " <<  ifname << endl;
    //continue;
    gSystem->Exit(1);
  }

  gnetwork* NET = new gnetwork(N_IFO,gIFO);

  wavearray<double>* sample_wREC[N_IFO];
  for(int i=0;i<N_IFO;i++) sample_wREC[i] = new wavearray<double>;

  wavearray<double>* sample_wINJ[N_IFO];
  for(int i=0;i<N_IFO;i++) sample_wINJ[i] = new wavearray<double>;

  double* time = new double[2*N_IFO];
  float likelihood;
  double tcoa=-1;               // geocentric tcoa
  float rho[2];
  float theta[4];
  float phi[4];
  float ifar=0;
  string* log = new string;     // injection log

  cout.precision(14);
  for(int k=0;k<itree->GetEntries();k++) {

    for(int i=0;i<N_IFO;i++) itree->SetBranchAddress(TString::Format("wREC_%d",i).Data(),&sample_wREC[i]);
    for(int i=0;i<N_IFO;i++) itree->SetBranchAddress(TString::Format("wINJ_%d",i).Data(),&sample_wINJ[i]);
    itree->SetBranchAddress("time",time);
    itree->SetBranchAddress("theta",theta);
    itree->SetBranchAddress("phi",phi);
    itree->SetBranchAddress("rho",rho);
    itree->SetBranchAddress("likelihood",&likelihood);
    itree->SetBranchAddress("wf_tcoa",&tcoa);
    itree->SetBranchAddress("log",&log);

    itree->GetEntry(k);                                               // read wREC,skyprob objects

    for(int i=0;i<N_IFO;i++) {
      wREC[gEVENTS].push_back(*sample_wREC[i]);             // store reconstructed posterior sample waveform
      wINJ[gEVENTS].push_back(*sample_wINJ[i]);             // store whitened posterior sample waveform
    }
    wTCOA[gEVENTS]=tcoa;                                    // store injected geocentric tcoa
    wTHETA[gEVENTS]=theta[1];                               // store OnSource sample theta
    wPHI[gEVENTS]=phi[1];                                   // store OnSource sample phi
    wRHO[gEVENTS]=rho[1];                                   // store OnSource sample rho

    ifar/=(24.*3600.*365.);                                 // sec -> year
    //cout << log->c_str() << endl;
    TString m1 = GetParameterFromInjLog(log->c_str(), "mass1");
    TString m2 = GetParameterFromInjLog(log->c_str(), "mass2");
    wM1[gEVENTS]=m1.Atof();                                 // store injected m1
    wM2[gEVENTS]=m2.Atof();                                 // store injected m2
    TString s1 = GetParameterFromInjLog(log->c_str(), "spin1");
    TString s2 = GetParameterFromInjLog(log->c_str(), "spin2");
    wS1[gEVENTS]=s1.Atof();                                 // store injected s1
    wS2[gEVENTS]=s2.Atof();                                 // store injected s2

    double mc = pow(m1.Atof()*m2.Atof(),3./5.)/pow(m1.Atof()+m2.Atof(),1./5.);	    // chirp mass
    cout << gEVENTS << " m1 " << m1 << " m2 " << m2 << " mc " << mc << " ifar " << ifar << " (years) " << endl;

    gEVENTS++;
  }
  froot->Close();

  if(wf_id>=gEVENTS) {
    cout << "Error : max wf_id = " <<  gEVENTS-1 << endl;
    gSystem->Exit(1);
  }

  cout << endl << "RHO  " << wRHO[wf_id] << endl << endl;

  // from wdm_fres and input data rate we compute the WDM levels and init WDM used for TF transform
  int ccut_wdm_levels = int(wREC[wf_id][ifo_id].rate()/gOPT.ccut_wdm_fres/2);
  cout << endl << "RATE " << wREC[wf_id][ifo_id].rate() << " WDM levels " << ccut_wdm_levels << endl << endl;
  gWDM = new WDM<double>(ccut_wdm_levels, ccut_wdm_levels, 6, 10);        // init a WDM transform 

  vector<double> vtcoa(N_IFO);
  for(int n=0;n<N_IFO;n++) {
    vtcoa[n] = GetInjTcoa(wTCOA[wf_id], NET, gIFO[n], wTHETA[wf_id], wPHI[wf_id]);
    cout << n << "\t" << int(wINJ[wf_id][n].start()) << "\t" << vtcoa[n]-wINJ[wf_id][n].start() << endl;
  }

  MDC = new CWB::mdc(N_IFO,gIFO);

#ifdef APPLY_SYNC
  double sync_phase=0, sync_time=0, sync_xcor=0;
  //sync_xcor = CWB::mdc::TimeSync(wREC[wf_id][ifo_id], wINJ[wf_id][ifo_id], sync_time);
  sync_xcor = CWB::mdc::TimePhaseSync(wREC[wf_id][ifo_id], wINJ[wf_id][ifo_id], sync_time, sync_phase);
  cout << "sync_time " << sync_time << " sync_phase " << sync_phase << " sync_xcor " << sync_xcor << endl;
#endif

/*
  wINJ[wf_id][ifo_id].start(0);
  wREC[wf_id][ifo_id].start(0);
  MDC->DrawTime(wINJ[wf_id][ifo_id],"ALP ZOOM");
  MDC->DrawTime(wREC[wf_id][ifo_id],"SAME",kRed);
  return;
*/
//  CWB::mdc::Align(wREC[wf_id][ifo_id], wINJ[wf_id][ifo_id]);
  wavearray<double> wDIF = CWB::mdc::GetDiff(&wREC[wf_id][ifo_id], &wINJ[wf_id][ifo_id]);
/*
  MDC->DrawFFT(wINJ[wf_id][ifo_id],"ALP ZOOM");
  MDC->DrawFFT(wREC[wf_id][ifo_id],"SAME",kRed);
  MDC->DrawFFT(wDIF,"SAME",kGreen+2);return;
*/

//  MDC->DrawTF(wREC[wf_id][ifo_id]); return;

  std::vector<wavearray<double> > vREC(N_IFO);
  std::vector<wavearray<double> > vINJ(N_IFO);
  for(int n=0;n<N_IFO;n++) {
    vINJ[n] = wINJ[wf_id][ifo_id];
    vREC[n] = wREC[wf_id][ifo_id];
  }
  vector<double> tstart,tstop; 
  tstart=vtcoa;
  tstop=vtcoa;
  for(int n=0;n<N_IFO;n++) tstart[n]-=0.25;
  double FF = CWB::mdc::GetMatchFactor("ff",vREC,vINJ,tstart,tstop);
  cout << "FF -> " << FF << endl;
  double RE = CWB::mdc::GetMatchFactor("re",vREC,vINJ,tstart,tstop);
  cout << "RE -> " << RE << endl;

#ifdef APPLY_BANDPASS
  vINJ[ifo_id] = CWB::mdc::GetBandpass(vINJ[ifo_id],80,400);
  vREC[ifo_id] = CWB::mdc::GetBandpass(vREC[ifo_id],80,400);
#endif

#ifdef APPLY_CCUT

#ifdef SHOW_WDM_INJ
  cout << endl;
  wINJ[wf_id][ifo_id] = GetCCut(&wINJ[wf_id][ifo_id], vtcoa[ifo_id]-wINJ[wf_id][ifo_id].start(), wM1[wf_id], wM2[wf_id], wS1[wf_id], wS2[wf_id]); 
  cout << "-------> APPLY_CCUT wINJ Energy (time domain phase 00) -------> " << pow(wINJ[wf_id][ifo_id].rms(),2)*wINJ[wf_id][ifo_id].size() << endl;
  cout << endl;
#endif

#ifdef SHOW_WDM_REC
  wREC[wf_id][ifo_id] = GetCCut(&wREC[wf_id][ifo_id], vtcoa[ifo_id]-wREC[wf_id][ifo_id].start(), wM1[wf_id], wM2[wf_id], wS1[wf_id], wS2[wf_id]); 
  cout << "-------> APPLY_CCUT wREC Energy (time domain phase 00) -------> " << pow(wREC[wf_id][ifo_id].rms(),2)*wREC[wf_id][ifo_id].size() << endl;
  cout << endl;
#endif

#ifdef SHOW_WDM_DIF
  wDIF = GetCCut(&wDIF, vtcoa[ifo_id]-wREC[wf_id][ifo_id].start(), wM1[wf_id], wM2[wf_id], wS1[wf_id], wS2[wf_id]); 
  cout << "-------> APPLY_CCUT wDIF Energy (time domain phase 00) -------> " << pow(wDIF.rms(),2)*wDIF.size() << endl;
  cout << endl;
#endif

#ifdef DEBUGGING	
  cout << "END " << endl;
  return;
#endif

#endif

#ifdef SHOW_FFT
  MDC->DrawFFT(wINJ[wf_id][ifo_id],"ALP ZOOM");
  MDC->DrawFFT(wREC[wf_id][ifo_id],"SAME",kRed);
#endif

#ifdef SHOW_TIME
  wINJ[wf_id][ifo_id].start(0);
  wREC[wf_id][ifo_id].start(0);
  MDC->DrawTime(wINJ[wf_id][ifo_id],"ALP ZOOM");
return;
  MDC->DrawTime(wREC[wf_id][ifo_id],"SAME",kRed);
#endif

#ifdef SHOW_ENVELOPE
  for(int n=0;n<N_IFO;n++) {
    vINJ[n] = CWB::mdc::GetEnvelope(&(wINJ[wf_id][ifo_id]));
    vREC[n] = CWB::mdc::GetEnvelope(&(wREC[wf_id][ifo_id]));
  }
  vINJ[ifo_id].start(0);
  vREC[ifo_id].start(0);

  MDC->DrawTime(vINJ[ifo_id],"ALP ZOOM");
  MDC->DrawTime(vREC[ifo_id],"SAME",kRed);
#endif

#ifdef SHOW_WDM_INJ
  DrawCCut(&wINJ[wf_id][ifo_id], vtcoa[ifo_id]-wINJ[wf_id][ifo_id].start(), wM1[wf_id], wM2[wf_id], wS1[wf_id], wS2[wf_id]); 
#endif
#ifdef SHOW_WDM_REC
  DrawCCut(&wREC[wf_id][ifo_id], vtcoa[ifo_id]-wINJ[wf_id][ifo_id].start(), wM1[wf_id], wM2[wf_id], wS1[wf_id], wS2[wf_id]); 
#endif
#ifdef SHOW_WDM_DIF
  DrawCCut(&wDIF, vtcoa[ifo_id]-wINJ[wf_id][ifo_id].start(), wM1[wf_id], wM2[wf_id], wS1[wf_id], wS2[wf_id]); 
#endif

}

double
GetInjTcoa(double geocentric_tcoa, network* NET, TString ifo, double theta, double phi) {
// compute coalescence time

  // Add Time Delay respect to geocenter
  CWB::mdc MDC(NET);
  double tdelay = MDC.GetDelay(ifo,"",phi,theta);
  double ifo_tcoa = geocentric_tcoa+tdelay;

  return ifo_tcoa;
}

TString GetParameterFromInjLog(TString log, TString param) {
// get params from log string                                                     

  TObjArray* token = log.Tokenize(TString(" "));
  for(int i=0;i<token->GetEntries();i++) {
    TObjString* otoken = (TObjString*)token->At(i);
    TString stoken = otoken->GetString();
    // exctract value
    if(stoken==param) {
      if(i<token->GetEntries()-1) {
        if(stoken=="spin1" || stoken=="spin2") {
          otoken = (TObjString*)token->At(i+3);         // get spinz
          return otoken->GetString();
        } else {
          otoken = (TObjString*)token->At(i+1);
          return otoken->GetString();
        }
      }
    }
  }
  if(token) delete token;

  return "";
}

void DrawCCut(wavearray<double>* ts, double tcoa, double m1, double m2, double s1, double s2) {  // apply transformation on white noise and plot it

  // produce the TF map:
  WSeries<double> W;           			// TF map container
  W.Forward(*ts, *gWDM);         		// apply the WDM to the time series 

  DrawCCut(W, tcoa, m1, m2, s1, s2);
}

void GetCCutParms(TString options) {

  // initialize with the default values
  gOPT.ccut_wdm_fres   = CCUT_WDM_FRES;
  gOPT.ccut_bchirp     = CCUT_BCHIRP;
  gOPT.ccut_uchirp     = CCUT_UCHIRP;
  gOPT.ccut_ltime      = CCUT_LTIME;
  gOPT.ccut_rtime      = CCUT_RTIME;

  // read the user parameter ccut
  TObjArray* token = TString(options).Tokenize(TString(':'));
  if(token->GetEntries()!=5) {
    cout << endl;
    cout << "cwb_pereport - Error : ccut format - must be: wdm_fres:bchirp:uchirp:ltime:rtime" <<  endl;
    cout << "               Default setup is: 16:1.35:1.65:0.5:0.0" << endl;
    cout << "               To setup a default value use *" << endl;
    cout << "               Example: *:1.25:1.75:*:*" << endl;
    cout << endl;
    exit(1);
  } 
  for(int j=0;j<token->GetEntries();j++) {
    TObjString* tok = (TObjString*)token->At(j);
    TString stok = tok->GetString();

    if(stok=="*") continue;

    if(j==0) {
      TString wdm_fres=stok;
      if(wdm_fres.IsDigit()) gOPT.ccut_wdm_fres=wdm_fres.Atoi();
      else {cout<<"cwb_pereport - Error : wdm_fres is not integer"<<endl;exit(1);}                   
    }
    if(j==1) {
      TString mc_lower_factor=stok;
      if(mc_lower_factor.IsFloat()) gOPT.ccut_bchirp=mc_lower_factor.Atof();
      else {cout<<"cwb_pereport - Error : mc_lower_factor is not float"<<endl;exit(1);}                   
    }
    if(j==2) {
      TString mc_upper_factor=stok;
      if(mc_upper_factor.IsFloat()) gOPT.ccut_uchirp=mc_upper_factor.Atof();
      else {cout<<"cwb_pereport - Error : mc_upper_factor is not float"<<endl;exit(1);}                   
    }
    if(j==3) {
      TString left_time=stok;
      if(left_time.IsFloat()) gOPT.ccut_ltime=left_time.Atof();
      else {cout<<"cwb_pereport - Error : left_time is not float"<<endl;exit(1);}                   
    }
    if(j==4) {
      TString right_time=stok;
      if(right_time.IsFloat()) gOPT.ccut_rtime=right_time.Atof();
      else {cout<<"cwb_pereport - Error : right_time is not float"<<endl;exit(1);}                   
    }
  }

  int x=gOPT.ccut_wdm_fres;        // must be power of 2
  if(!((x != 0) && ((x & (~x + 1)) == x)) || gOPT.ccut_wdm_fres<=0) {
    cout<<"cwb_pereport : upTDF  parameter non valid : must be power of 2 : "<<gOPT.ccut_wdm_fres<<endl;
    exit(1);
  }
  if(gOPT.ccut_bchirp < 0) {
    cout<<"cwb_pereport - Error :.ccut_bchirp must be >= 0"<<endl;
    exit(1);
  }
  if(gOPT.ccut_bchirp > gOPT.ccut_uchirp) {
    cout<<"cwb_pereport - Error :.ccut_bchirp must be < ccut_uchirp"<<endl;
    exit(1);
  }
  if(gOPT.ccut_rtime < 0) {
    cout<<"cwb_pereport - Error :.ccut_rtime >=0"<<endl;
    exit(1);
  }
  if(gOPT.ccut_ltime < gOPT.ccut_rtime) {
    cout<<"cwb_pereport - Error :.ccut_ltime must be >.ccut_rtime"<<endl;
    exit(1);
  }
}

wavearray<double> GetCCut(wavearray<double>* ts, double tcoa, double m1, double m2, double s1z, double s2z) {
// apply chirp TF cuts and compute total energy inside the chirp cut region

  const auto& fchirp = static_cast<double(*)(double,double,double,double,double)>(CWB::mdc::SimIMRSEOBNRv4ROMFrequencyOfTime);
  const auto& tchirp = static_cast<double(*)(double,double,double,double,double)>(CWB::mdc::SimIMRSEOBNRv4ROMTimeOfFrequency);

  // produce the TF map:
  WSeries<double> W;                            // TF map container
  W.Forward(*ts, *gWDM);                        // apply the WDM to the time series 

  double fmin=16;
  double fmax= ts->rate()/2.>320 ? 320. : ts->rate()/2.;

  //#define nPx	100
  #define nPx	30
  double T[nPx],Fb[nPx],Fu[nPx];
  double dF = (fmax-fmin)/nPx;
  for(int i=0;i<nPx;i++) {
    double F = fmin+i*dF;
    Fb[i] = gOPT.ccut_bchirp*F;			// alpha bottom
    Fu[i] = gOPT.ccut_uchirp*F;			// alpha up
    T[i]  = tcoa-tchirp(F,m1,m2,s1z,s2z);
  }
  TSpline3 tbchirp("tbchirp",Fb,T,nPx);		// bottom chirp line 	(time vs freq)
  TSpline3 tuchirp("tuchirp",Fu,T,nPx);		// up chirp line	(time vs freq)
  TSpline3 fbchirp("fbchirp",T,Fb,nPx);		// bottom chirp line	(freq vs time)
  TSpline3 fuchirp("fuchirp",T,Fu,nPx);		// up chirp line	(freq vs time)

  double df = W.resolution();                   // frequency pixel resolution (hz)
  double dt = 1./(2*df);                        // time pixel resolution (sec)

  int layers = W.maxLayer()+1;                  // numbers of frequency bins (first & last bins have df/2)
  int slices = W.sizeZero();                    // number of time bins

  double tl = tcoa-gOPT.ccut_ltime;		// left time cut
  double tr = tcoa-gOPT.ccut_rtime;		// right time cut

  // rescale the phase 00/90 pixel amplitudes

  for(int i=1;i<slices;i++) {                   // loop over time bins

    double time  = i*dt;			// pixel central time

    double tm = time-dt/2;			// pixel minimum time 
    double tM = time+dt/2;			// pixel maximum time

    if(tm<tl-dt || tM>tr+dt) {			// remove pixels outside the time interval [tl-dt:tr+dt]
      for(int j=1;j<layers;j++) {               // loop over frequency bins
        W.putSample(0,i,j+0.01);    		// zero a00 pixel amplitude
        W.putSample(0,i,-(j+0.01)); 		// zero a90 pixel amplitude
      }
      continue;
    }

    double fb = fbchirp.Eval(time);		// fbchirp frequency at pixel central time
    double fu = fuchirp.Eval(time);		// fuchirp frequency at pixel central time

    if(TMath::IsNaN(fb)) fb = ts->rate()/2.;	// check if fb is NaN
    if(TMath::IsNaN(fu)) fu = ts->rate()/2.;	// check if fu is NaN

    double f1 = fbchirp.Eval(tm);		// fbchirp frequency at pixel minimum time
    double f2 = fbchirp.Eval(tM);		// fbchirp frequency at pixel maximum time
    double f3 = fuchirp.Eval(tm);		// fuchirp frequency at pixel minimum time
    double f4 = fuchirp.Eval(tM);		// fuchirp frequency at pixel maximum time

    for(int j=1;j<layers;j++) {                 // loop over frequency bins

      double freq = j*df;			// pixel central frequency

      double fm = freq-df/2;			// pixel minimum frequency
      double fM = freq+df/2;			// pixel maximum frequency

      if(fm<fmin || fM>fmax) {
        W.putSample(0,i,j+0.01);    		// zero a00 pixel amplitude
        W.putSample(0,i,-(j+0.01)); 		// zero a90 pixel amplitude
        continue;
      }

      double t1 = tbchirp.Eval(fm);		// fbchirp time at pixel minimum frequency
      double t2 = tbchirp.Eval(fM);		// fbchirp time at pixel maximum frequency
      double t3 = tuchirp.Eval(fm);		// fuchirp time at pixel minimum frequency
      double t4 = tuchirp.Eval(fM);		// fuchirp time at pixel maximum frequency

      // check if pixel is inside the time cut range tl,tr
      bool ctime = (tM>tl && tm<tr) ? true : false; 

      // check if pixel is inside the frequency cut range fb,fu
      bool cfreq = (fM>fb && fm<fu) ? true : false; 

      double x1=t1,x2=t2,x3=t3,x4=t4,xm=tm,xM=tM;	// assign x* = t*
      double y1=f1,y2=f2,y3=f3,y4=f4,ym=fm,yM=fM;	// assign y* = f*

      // for pixels crossed by the cut times tl,tr assign new time pixel boundaries xm,xM and new f*
      if(tl>tm && tl<tM) {xm=tl; f1 = fbchirp.Eval(xm); f3 = fuchirp.Eval(xm);}
      if(tr>tm && tr<tM) {xM=tr; f2 = fbchirp.Eval(xM); f4 = fuchirp.Eval(xM);} 

      y1=f1,y2=f2,y3=f3,y4=f4;				// update y*

      // check if pixel is crossed by fbchirp,fuchirp lines
      // a pixel is crossed by the chirp lines if at least one of t*,f* points are inside the time/freq pixel ranges

      bool cline=false;							// false/true = not-crossed/crossed

      if(t1>xm && t1<xM) cline=true; 
      if(t2>xm && t2<xM) cline=true; 
      if(t3>xm && t3<xM) cline=true; 
      if(t4>xm && t4<xM) cline=true; 

      if(f1>ym && f1<yM) cline=true; 
      if(f2>ym && f2<yM) cline=true; 
      if(f3>ym && f3<yM) cline=true; 
      if(f4>ym && f4<yM) cline=true; 

      double Ap=(xM-xm)*(yM-ym);					// area pixel (is <=0.5 = total pixel area) 

      if(ctime&&cline) {	// pixels crossed by fbchirp,fuchirp lines and inside the time cut range tl,tr

        // for pixels with frequency outside the frequency pixel range assign ym,yM
        if(y1<ym) y1=ym; if(y1>yM) y1=yM;
        if(y2<ym) y2=ym; if(y2>yM) y2=yM;
        if(y3<ym) y3=ym; if(y3>yM) y3=yM;
        if(y4<ym) y4=ym; if(y4>yM) y4=yM;

	// compute area Ab for pixels crossed by fbchirp line 

        double Ab=0;							// area pixel above the fbchirp line

        // for pixels with times t1,t2 outside the time pixel range assign xm,xM and compute the corresponding new frequencies y1,y2
        if(t1<xm) {x1=xm; y1 = fbchirp.Eval(x1);}
        if(t1>xM) {x1=xM; y1 = fbchirp.Eval(x1);}
        if(t2<xm) {x2=xm; y2 = fbchirp.Eval(x2);}
        if(t2>xM) {x2=xM; y2 = fbchirp.Eval(x2);}

        // for pixels with frequency outside the frequency pixel range assign ym,yM
        if(y1<ym) y1=ym; if(y1>yM) y1=yM;
        if(y2<ym) y2=ym; if(y2>yM) y2=yM;

        // compute area pixel above the fbchirp line
        if(y1 >ym && y2==yM) Ab = (x2-xm)*(yM-y1)/2;
        if(y1 >ym && y2 <yM) Ab = (xM-xm)*(yM-y2)+(xM-xm)*(y2-y1)/2; 
        if(y1==ym && y2==yM) Ab = (x1-xm)*(yM-ym)+(x2-x1)*(yM-ym)/2;
        if(y1==ym && y2 <yM) Ab = Ap-(xM-x1)*(y2-ym)/2;

	// compute area Au for pixels crossed by fuchirp line 

        double Au=0;							// area pixel below the fuchirp line

        // for pixels with times t3,t4 outside the time pixel range assign xm,xM and compute the corresponding new frequencies y3,y4
        if(t3<xm) {x3=xm; y3 = fuchirp.Eval(x3);}
        if(t3>xM) {x3=xM; y3 = fuchirp.Eval(x3);}
        if(t4<xm) {x4=xm; y4 = fuchirp.Eval(x4);}
        if(t4>xM) {x4=xM; y4 = fuchirp.Eval(x4);}

        // for pixels with frequency outside the frequency pixel range assign ym,yM
        if(y3<ym) y3=ym; if(y3>yM) y3=yM;
        if(y4<ym) y4=ym; if(y4>yM) y4=yM;

        // compute area pixel below the fuchirp line
        if(y3 >ym && y4==yM) Au = Ap-(x4-xm)*(yM-y3)/2;
        if(y3 >ym && y4 <yM) Au = (xM-xm)*(y3-ym)+(xM-xm)*(y4-y3)/2; 
        if(y3==ym && y4==yM) Au = (xM-x4)*(yM-ym)+(x4-x3)*(yM-ym)/2;
        if(y3==ym && y4 <yM) Au = (xM-x3)*(y4-ym)/2;

        if(Ab<0) Ab=0; if(Ab>Ap) Ab=Ap;                                 // fix area inconsistency due to precision errors
        if(Au<0) Au=0; if(Au>Au) Au=Au;                                 // fix area inconsistency due to precision errors
        double A = Ab+Au; A-=Ap;                                        // combine Ab and Au when pixel is crossed by fbchirp fuchirp lines

        if(A<0)   {cout << "GetCCut Warning: pixel area < 0  " << endl; A=0;}
        if(A>Ap)  {cout << "GetCCut Warning: pixel area > Ap " << endl; A=Ap;}
        double R=sqrt(A/0.5);						// compute the rescaling amplitude factor (0.5 is the total pixel area)
#ifdef DEBUGGING
        TString tag = R==0 ? "*" : "";
        TString sR = TString::Format("A - GetCCut Rescaling Factor: Ap %0.3f %0.3f %0.3f %0.3f \ttime %0.3f \tfreq %0.3f \tR %0.9f \t%s",Ap,Ab,Au,A,time,freq,R,tag.Data());
        cout << sR << endl;
#endif
        double a00 = W.getSample(i,j+0.01); W.putSample(a00*R,i,j+0.01);	// rescaling a00 pixel amplitude
        double a90 = W.getSample(i,-(j+0.01));W.putSample(a90*R,i,-(j+0.01));	// rescaling a90 pixel amplitude
      }

      if((ctime && cfreq) && !cline) { 	// pixels inside/outside the chirp cut area not belonging to the cline pixels 

        double R = (x4<=tl || x1>=tr) ? 0 : sqrt(Ap/0.5);	// compute the rescaling amplitude factor (0.5 is the total pixel area)
								// remove pixels outside the time interval [tl,tr]

#ifdef DEBUGGING
        TString tag = "+";
        TString sR = TString::Format("B - GetCCut Rescaling Factor: Ap %0.3f \ttime %0.3f \tfreq %0.3f \tR %0.3f \t%s",Ap,time,freq,R,tag.Data());
        cout << sR << endl;
#endif
        double a00 = W.getSample(i,j+0.01); W.putSample(a00*R,i,j+0.01);	// rescaling a00 pixel amplitude
        double a90 = W.getSample(i,-(j+0.01));W.putSample(a90*R,i,-(j+0.01));	// rescaling a90 pixel amplitude
      }

      //if(ctime && cfreq) { 
      //if(ctime && cline) { 
      //if((ctime && cfreq) && !cline) { 
      if((ctime && cline)||(ctime && cfreq)) { 
      } else {
        W.putSample(0,i,j+0.01);
        W.putSample(0,i,-(j+0.01));
      }
    }
  }

  double Et=0;                                  // Energy in time-frequency chirp cut region
  for(int i=1;i<slices;i++) {                   // loop over time bins
    for(int j=1;j<layers;j++) {                 // loop over frequency bins
      double a00 = W.getSample(i,j+0.01);       // read a00 pixel amplitude
      double a90 = W.getSample(i,-(j+0.01));    // read a90 pixel amplitude
      Et+=(a00*a00+a90*a90)/2;
    }
  }

  cout << "GetCCut Energy chirp cut (TF domain): " << Et << endl;

#ifdef DEBUGGING
  DrawCCut(W, tcoa, m1, m2, s1z, s2z); 
#endif

  // return to time domain
  WSeries<double> xW = W;
  W.Inverse();
  xW.Inverse(-2);
  W += xW;
  W *= 0.5;

  return W;
}

void DrawCCut(WSeries<double> W, double tcoa, double m1, double m2, double s1z, double s2z) {  // apply transformation on white noise and plot it

  const auto& fchirp = static_cast<double(*)(double,double,double,double,double)>(CWB::mdc::SimIMRSEOBNRv4ROMFrequencyOfTime);
  const auto& tchirp = static_cast<double(*)(double,double,double,double,double)>(CWB::mdc::SimIMRSEOBNRv4ROMTimeOfFrequency);

  #define nPX2	100

  double G  = watconstants::GravitationalConstant();
  double SM = watconstants::SolarMass();
  double C  = watconstants::SpeedOfLightInVacuo();
  double D  = G*SM/C/C/C;

 double tmax=tcoa;
 double tmin=tmax-0.5;

  double fmin=16;
  double fmax= W.rate()/2.>320 ? 320. : W.rate()/2.;

  #define nPx2	100
  //#define nPx2	30
  double T[nPx2],Fb[nPx2],Fu[nPx2];
  double dF = (fmax-fmin)/nPx2;
  for(int i=0;i<nPx2;i++) {
    double F = fmin+i*dF;
    Fb[i] = gOPT.ccut_bchirp*F;			// alpha bottom
    Fu[i] = gOPT.ccut_uchirp*F;			// alpha up
    T[i]  = tcoa-tchirp(F,m1,m2,s1z,s2z);
  }
  TSpline3 tbchirp("tbchirp",Fb,T,nPx2);		// bottom chirp line 	(time vs freq)
  TSpline3 tuchirp("tuchirp",Fu,T,nPx2);		// up chirp line	(time vs freq)
  TSpline3 fbchirp("fbchirp",T,Fb,nPx2);		// bottom chirp line	(freq vs time)
  TSpline3 fuchirp("fuchirp",T,Fu,nPx2);		// up chirp line	(freq vs time)

  // plot the TF map as energy average of 0 and 90 degree phases (quadratures)
  TString wplot_id=TString::Format("%d",int(gRandom->Uniform(0,10000000)));
  watplot* wplot = new watplot(const_cast<char*>(wplot_id.Data()));
  wplot->canvas->SetLogy();
  wplot->plot(W);
  wplot->hist2D->GetXaxis()->SetRangeUser(6.4,7.2);
  wplot->hist2D->GetYaxis()->SetRangeUser(0,600);
  wplot->hist2D->GetYaxis()->SetRangeUser(30,600);

  double xb[nPX2],yb[nPX2];
  double dxb=(tmax-tmin)/nPX2;
  for(int i=0;i<nPX2;i++) {
    xb[i]=tmin+i*dxb;
    //yb[i]=fbchirp->GetX(xb[i]);
    yb[i]=fbchirp.Eval(xb[i]);
  }
  TGraph* grb = new TGraph(nPX2,xb,yb);
  grb->SetLineColor(kWhite);
  grb->SetLineWidth(3);
  grb->SetLineStyle(2);
  grb->Draw("same");

  double xu[nPX2],yu[nPX2];
  double dxu=(tmax-tmin)/nPX2;
  for(int i=0;i<nPX2;i++) {
    xu[i]=tmin+i*dxu;
    //yu[i]=fuchirp->GetX(xu[i]);
    yu[i]=fuchirp.Eval(xu[i]);
  }
  TGraph* gru = new TGraph(nPX2,xu,yu);
  gru->SetLineColor(kWhite);
  gru->SetLineWidth(3);
  gru->SetLineStyle(2);
  gru->Draw("same");

  double p0;

  // left time chirp cut line
  p0 = tcoa-gOPT.ccut_ltime;		// left time cut
  flchirp = new TLine(p0,0,p0,600);
  flchirp->SetLineColor(kBlack);
  flchirp->SetLineWidth(3);
  flchirp->SetLineStyle(2);
  flchirp->Draw("same");

  // right time chirp cut line
  p0 = tcoa-gOPT.ccut_rtime;		// left time cut
  frchirp = new TLine(p0,0,p0,600);
  frchirp->SetLineColor(kBlack);
  frchirp->SetLineWidth(3);
  frchirp->SetLineStyle(2);
  frchirp->Draw("same");

  // right time chirp cut line
  p0 = tcoa;				// tcoa time
  ftcoa = new TLine(p0,0,p0,600);
  ftcoa->SetLineColor(kRed);
  ftcoa->SetLineWidth(2);
  ftcoa->SetLineStyle(1);
  ftcoa->Draw("same");
}

