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

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!DATA_CONDITIONING
// Plugin to produce frames of whitened data
// To display data with FrDisplay do (Ex) :
// FrDisplayPROC -t L1:LDAS-STRAIN \
// -i report/dump/white_H1_S6A_BKG_L1H1V1_HF_2G_run2_job535-932501639-600.gwf 
//

  cout << endl;
  cout << "-----> CWB_Plugin_MakeWhiteFrame.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {

    // get ifo id
    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    // transform to time domain 
    int level=-1;
    if(TString(cfg->analysis)=="1G") {
      level=x->getLevel();
      x->Inverse(-1);
    }

    // save cleaned data into frame
    wavearray<double> X = *x;
    // remove scratch
    int os = cfg->segEdge*x->rate();
    X.start(x->start()+cfg->segEdge);
    X.stop(x->stop()-cfg->segEdge);
    X.resize(x->size()-2*os);
    for(int i=0;i<X.size()-2*os;i++) X[i]=X[i+os];

    // create frame file
    char frFile[1024];
    char frName[1024];
    char chName[1024]; 
    sprintf(frFile,"%s/white_%s_%s_job%lu-%lu-%lu.gwf",
            cfg->dump_dir,ifo.Data(),cfg->data_label,net->nRun,
            int(X.start()),(int)(X.stop()-X.start()));
    sprintf(frName,"WFRAME");
    sprintf(chName,cfg->channelNamesRaw[id]);
    // open frame file
    CWB::frame fr(frFile,chName,"WRITE");
    // write frame to file
    fr.writeFrame(X, frName, chName);
    cout << "CWB_Plugin_MakeWhiteFrame.C : write " << frFile << endl;

    fr.close();

    int nIFO=net->ifoListSize();
    if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) gSystem->Exit(0);  // last ifo
    if(TString(cfg->analysis)=="1G") x->Forward(level);
  }

  return;
}
