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


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  cout << endl;
  cout << "-----> macro/CWB_Plugin_MergeStrainMDC.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_STRAIN) {  
    CWB::frame   fr;
    frfile       FRF;
    int          nIFO = cfg->nIFO;

    int id=-1;
    for(int n=0;n<nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    // set frame file list
    fr.open(cfg->frFiles[id]);
    fr.setVerbose();
    fr.setRetryTime(cfg->frRetryTime);
    fr.setSRIndex(14);   // force readFrames to resample to 16384 Hz
    int nfrFiles=fr.getNfiles();
    cout << "MDC " << " -> nfrFiles : " << nfrFiles << endl;
    FRF = fr.getFrList(int(x->start()), int(x->stop()), 0.);

    wavearray<double> mdc;
    fr.readFrames(FRF,cfg->channelNamesMDC[id],mdc);

    mdc*=cfg->factors[0];
    *x+=mdc;
  }

  return;
}

