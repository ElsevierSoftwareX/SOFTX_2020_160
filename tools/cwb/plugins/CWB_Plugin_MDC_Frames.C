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


//!NOISE_MDC_SIMULATION
// Plugin used to setup frame files and log MDC 'on the fly' using simulation=4

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
#include "TComplex.h"
#include "TMath.h"
#include "cwb2G.hh"
#include "mdc.hh"
#include <vector>


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->mdcPlugin=true;         // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_MDC) {  

    cout << "Execute CWB_Plugin_MDC_Frames.C : Inject On The Fly MDC ..." << endl;

    int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {ifoID=n;break;}

    int nIFO     = net->ifoListSize();          // number of detectors
    int gIFACTOR = -1; IMPORT(int,gIFACTOR)
    int xstart   = (int)x->start();
    int xstop    = (int)x->stop();
    int segEdge  = int(cfg->segEdge);

    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)

    // Export parms to the config plugin
    static bool cfg_init = false;
    if(!cfg_init) {	// declaration is made only one time to avoid cling error
      CWB_PLUGIN_EXPORT(cfg)
      cfg_init=true;
    }
    CWB_PLUGIN_EXPORT(gIFACTOR)

    // Export MDC frame files & log from plugin config
    std::vector<TString> mdcPar;
    CWB_PLUGIN_EXPORT(mdcPar)
    // Export mdcType
    std::vector<std::string> mdcType;
    CWB_PLUGIN_EXPORT(mdcType)
    // Export channelName
    std::vector<std::string> channelName;
    CWB_PLUGIN_EXPORT(channelName)

    // read plugin config
    cfg->configPlugin.Exec();

    TString injectionList = mdcPar[0];
    TString frFiles = mdcPar[1];
    cout << "injectionList : " << injectionList << endl;
    cout << "frFiles       : " << frFiles << endl;

    // open frame file
    CWB::frame fr;
    fr.open(frFiles.Data());
    fr.setVerbose();
    fr.setRetryTime(cfg->frRetryTime);
    fr.setSRIndex(TMath::Nint(TMath::Log2(cfg->inRate))); // force data to be resampled to inRate
    int nfrFiles = fr.getNfiles();
    cout << "MDC " << " -> nfrFiles : " << nfrFiles << endl;
    if(cfg->mdc_shift.startMDC<0) {  // the mdc shift is automatically selected
      // read mdc range from mdc frl files
      waveSegment mdc_range  = fr.getFrRange();
      cout << "mdc_range : " << mdc_range.start << " " << mdc_range.stop << endl;
      cfg->mdc_shift.startMDC=mdc_range.start;
      cfg->mdc_shift.stopMDC=mdc_range.stop;
    }
    double mdcShift = CWB::Toolbox::getMDCShift(cfg->mdc_shift, xstart);
    cout << "mdcShift : " << mdcShift << endl;

    // read injection list
    net->mdcList.clear();
    net->mdcType.clear();
    net->mdcTime.clear();
    vector<TString> ifos(nIFO);
    for(int n=0;n<nIFO;n++) ifos[n] = net->getifo(n)->Name;
    int i=net->readMDClog(const_cast<char*>(injectionList.Data()),double(long(gCWB2G->Tb))-mdcShift);
    for(int i=0;i<(int)net->mdcTime.size();i++) net->mdcTime[i]+=mdcShift;
    printf("GPS: %16.6f saved,  injections: %d\n",double(long(gCWB2G->Tb)),i);
    net->mdcType=mdcType;

    // read frame file 
    frfile FRF = fr.getFrList(xstart-mdcShift+segEdge, xstop-mdcShift-segEdge, segEdge);
    fr.readFrames(FRF,const_cast<char*>(channelName[ifoID].c_str()),*x);
    x->start(x->start()+mdcShift);
    if(x->rate()!=cfg->inRate)
      {cout << "CWB_Plugin_MDC_Frames.C - input rate from frame " << x->rate()
            << " do not match the one defined in config : " << cfg->inRate << endl;gSystem->Exit(1);}

    // print MDC injections list 
    cout.precision(14);
    if(ifo.CompareTo(net->ifoName[0])==0) {
      //for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << net->mdcList[k] << endl;
      //for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << net->mdcTime[k] << endl;
      for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << net->mdcType[k] << endl;
    }

  }

  return;
}
