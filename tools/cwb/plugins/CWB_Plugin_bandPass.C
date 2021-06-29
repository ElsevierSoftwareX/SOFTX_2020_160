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


#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "Toolbox.hh"

#define DUMP_WHITE_PSD

struct fBand {
  double flow;
  double fhigh;
};

#define nBAND 6 

fBand fband[nBAND]={
                    {77,83},
                    {116,125},
                    {175,185},
                    {196,204},
                    {491,497},
                    {573,580}
                   };

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!DATA_CONDITIONING
// Plugin which implement the bandFilter data conditioning on whitened data

  cout << endl;
  cout << "-----> CWB_Plugin_bandPass.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_WHITE) {
    CWB::Toolbox TB;
    int level=-1;
    if(TString(cfg->analysis)=="1G") {
      level=x->getLevel();
      x->Inverse(-1);
    }

    wavearray<double>* px = (wavearray<double>*) x;

    int layers = px->rate()/2;
    WDM<double> wdm(layers,layers,6,10);  			// df = (px->rate()/2)/layers
    detector D;
    D.getTFmap()->Forward(*px,wdm);
    for(int i=0;i<nBAND;i++) D.bandCut(fband[i].flow,fband[i].fhigh);	// band cut  filter
    D.getTFmap()->Inverse();              			// return to time domain
    *px = *D.getTFmap();  					// get band filtered data

#ifdef DUMP_WHITE_PSD
    // dump spectrum
    char file[1024];
    sprintf(file,"%s/sensitivity_white_%s_%d_%s_job%lu.txt",
            cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x);
    int nIFO=net->ifoListSize();
    if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) gSystem->Exit(0);  // last ifo
#endif

    if(TString(cfg->analysis)=="1G") x->Forward(level);
  }

  return;
}
