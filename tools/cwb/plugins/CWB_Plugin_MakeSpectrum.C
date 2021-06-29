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
//!DISPLAY_SPECTRA
// Plugin to produce spectrum of input/whitened data

  cout << endl;
  cout << "-----> CWB_Plugin_MakeSpectrum.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;
/*
  if(type==CWB_PLUGIN_DATA) {
    CWB::Toolbox TB;
    char file[1024];
    sprintf(file,"%s/sensitivity_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    //sprintf(file,"%s/sensitivity_%s_%d_%s_job%lu.txt",".",ifo.Data(),int(x->start()),"prova",net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 8, cfg->segEdge);

    int nIFO=net->ifoListSize();
    if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) gSystem->Exit(0);  // last ifo
  }
*/
  if(type==CWB_PLUGIN_WHITE) {
    CWB::Toolbox TB;
    int level=-1;
    if(TString(cfg->analysis)=="1G") {
      level=x->getLevel();
      x->Inverse(-1);
    }
    // dump spectrum
    char file[1024];
    sprintf(file,"%s/sensitivity_white_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 8, cfg->segEdge);
    int nIFO=net->ifoListSize();
    if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) gSystem->Exit(0);  // last ifo
    if(TString(cfg->analysis)=="1G") x->Forward(level);
  }

  return;
}
