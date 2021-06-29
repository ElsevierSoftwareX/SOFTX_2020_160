#define XIFO 5

#pragma GCC system_header

#include "cwb.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TMath.h"
#include "mdc.hh"
#include <vector>

//#define DUMP_DATA_SPECTRUM
//#define DUMP_STRAIN_AND_MDC_SPECTRUM
//#define CHECK_SNR

void GetGlitchList(TString ifo, std::vector<double>&  gpsList, std::vector<double>&  snrList);
void SetSNR(wavearray<double>* x, std::vector<double>  gpsList, std::vector<double>  snrList, CWB::config* cfg);
void CheckSNR(wavearray<double>* x, std::vector<double>  gpsList, std::vector<double>  snrList, CWB::config* cfg);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  cout << endl;
  cout << "-----> macro/CWB_Plugin_OnlineGlitches.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

#ifdef DUMP_DATA_SPECTRUM
  if(type==CWB_PLUGIN_DATA) {
    CWB::Toolbox TB;
    char file[1024];
    sprintf(file,"%s/sensitivity_%s_%d_%s_job%lu.txt",
            cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 2, cfg->segEdge);

    int nIFO=net->ifoListSize();
    if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) gSystem->Exit(0);  // last ifo
  }
#endif

#ifdef DUMP_STRAIN_AND_MDC_SPECTRUM
  if(type==CWB_PLUGIN_STRAIN_AND_MDC) {  
    CWB::Toolbox TB;
    int level=-1;
    if(TString(cfg->analysis)=="1G") {
      level=x->getLevel();
      x->Inverse(-1);
    }
    // dump spectrum
    char file[1024];
    sprintf(file,"%s/sensitivity_strain_mdc_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 8, cfg->segEdge);
    int nIFO=net->ifoListSize();
    if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) gSystem->Exit(0);  // last ifo
    if(TString(cfg->analysis)=="1G") x->Forward(level);
  }
#endif

  if(type==CWB_PLUGIN_STRAIN_AND_MDC) {  

    int level=-1;
    if(TString(cfg->analysis)=="1G") {
      level=x->getLevel();
      x->Inverse(-1);
    }

    CWB::frame fr;
    int nfrFiles;
    frfile FRF;

    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    fr.open(cfg->frFiles[id]);
    fr.setVerbose();
    fr.setRetryTime(cfg->frRetryTime);
    fr.setSRIndex(14);   	// force readFrames to resample to 16384 Hz
    nfrFiles=fr.getNfiles();
    cout << ifo << " -> nfrFiles : " << nfrFiles << endl;
    double start = x->start();
    double stop  = start + x->size()/x->rate();
    FRF = fr.getFrList(start, stop, 0);

    wavearray<double> y;
    char channelNamesGlitch[256];sprintf(channelNamesGlitch,"%s:FAKE-GLITCHES",ifo.Data());   
    fr.readFrames(FRF,channelNamesGlitch,y);

    // resample glitch channel
    WSeries<double> wM;              // WSeries
    Meyer<double> B(1024);           // set wavelet for resampling
    wM.Forward(y,B,cfg->levelR);
    wM.getLayer(y,0);

    if(y.rate()!=x->rate())
      {cout << "CWB_Plugin_OnlineGlitches.C - glitch rate " << y.rate()
            << " do not match the mdc rate : " << x->rate() << endl;gSystem->Exit(1);}

    // set snr glitches (glitches in channelNamesGlitch have snr=1)
    std::vector<double>  gpsList;
    std::vector<double>  snrList;
    GetGlitchList(ifo, gpsList, snrList);
    SetSNR(&y, gpsList, snrList, cfg);

    *x+=y; 			// add glitches to mdc data
    if(TString(cfg->analysis)=="1G") x->Forward(level); // restore wavelet transform
  }

#ifdef CHECK_SNR
  if(type==CWB_PLUGIN_WHITE) {
    int level=-1;
    if(TString(cfg->analysis)=="1G") {
      level=x->getLevel();
      x->Inverse(-1);
    }

    // check snr glitches in whitened data
    std::vector<double>  gpsList;
    std::vector<double>  snrList;
    GetGlitchList(ifo, gpsList, snrList);
    CheckSNR(x, gpsList, snrList, cfg);

    if(TString(cfg->analysis)=="1G") x->Forward(level); // restore wavelet transform
  }
#endif

  return;
}

void
GetGlitchList(TString ifo, std::vector<double>&  gpsList, std::vector<double>&  snrList) {

  TString fName;

  if(ifo=="L1") fName="input/L1_glitches.lst";
  if(ifo=="H1") fName="input/H1_glitches.lst";
  if(ifo=="V1") fName="input/V1_glitches.lst";

  ifstream in;
  in.open(fName.Data(),ios::in);
  if (!in.good()) {cout << "CWB_Plugin_OnlineGlitches.C - Error Opening File : " << fName.Data() << endl;gSystem->Exit(1);}

  char str[1024];
  char name[128];
  double gps,theta,phi,psi,rho,iota,snr;  // geographic coordinates
  int ID,id;
  int fpos=0;
  while (1) {
    fpos=in.tellg();
    in.getline(str,1024);
    if(str[0] == '#') continue;
    if (!in.good()) break;

    std::stringstream linestream(str);
    if(!(linestream >> gps >> name >> theta >> phi >> psi >> rho >> iota >> snr >> ID >> id)) {
       cout << "CWB_Plugin_OnlineGlitches.C - Wrong Format for File : " << fName.Data() << endl;
       cout << "input line : " << endl;
       cout << str << endl;
       cout << "must be : " << endl;
       cout << "gps " << "name " <<  "theta " <<  "phi " <<  "psi " << "rho " << "iota " << "snr " << "..." << endl;
       gSystem->Exit(1);
    }
    gpsList.push_back(gps);
    snrList.push_back(snr);
    //cout << gps << " " <<  name << " " <<  theta << " " <<  phi << " " <<  psi << " " <<  rho << endl;
  }
  in.close();

  return;
}

void 
SetSNR(wavearray<double>* x, std::vector<double>  gpsList, std::vector<double>  snrList, CWB::config* cfg) {

  int K = gpsList.size();
  double rate = x->rate();
  int N = x->size();      
  //double dT=fabs(cfg->gap/2);            
  double dT=0.5;            

  // isolate injections
  for(int k=0; k<K; k++) {

    double T = gpsList[k] - x->start();

    int nstrt = int((T - dT)*rate);
    int nstop = int((T + dT)*rate);
    if(nstrt<=0) nstrt = 0;        
    if(nstop>=int(N)) nstop = N;   
    if(nstop<=0) continue;                     // outside of the segment
    if(nstrt>=int(N)) continue;                // outside of the segment

    double snr = snrList[k];              
    for(int i=nstrt;i<nstop;i++) x->data[i]*=snr;	// set snr

    // compute hrss
    double hrss=0;
    for(int i=nstrt;i<nstop;i++) hrss+=pow(x->data[i],2);	
    hrss*=1./rate;
    hrss=sqrt(hrss);
    cout.precision(14);
    cout << "glitch : " << k << " gps " << gpsList[k] << " hrss : " << hrss << endl;
  }                                                        
}                                                          


void 
CheckSNR(wavearray<double>* x, std::vector<double>  gpsList, std::vector<double>  snrList, CWB::config* cfg) {

  int K = gpsList.size();
  double rate = x->rate();
  int N = x->size();      
  //double dT=fabs(cfg->gap/2);            
  double dT=0.5;            

  // isolate injections
  for(int k=0; k<K; k++) {

    double T = gpsList[k] - x->start();

    int nstrt = int((T - dT)*rate);
    int nstop = int((T + dT)*rate);
    if(nstrt<=0) nstrt = 0;        
    if(nstop>=int(N)) nstop = N;   
    if(nstop<=0) continue;                     // outside of the segment
    if(nstrt>=int(N)) continue;                // outside of the segment

    // compute snr
    double snr=0;
    for(int i=nstrt;i<nstop;i++) snr+=pow(x->data[i],2);	
    snr*=1./rate;
    snr=sqrt(snr);
    cout.precision(14);
    cout << "Check glitch snr : " << k << " gps " << gpsList[k] << " snr : " << snr << " inj snr : " << snrList[k] << endl;
  }                                                        
}                                                          

