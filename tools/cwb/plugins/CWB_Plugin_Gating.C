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
#include "cwb2G.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "Toolbox.hh"

//!DATA_CONDITIONING

// This plugin implements the time gating in the pixel's selection stage.
// Gating is a veto of the pixels belonging to a time interval.
// This plugin is used to exclude from the analysis the pixels 
// in a time interval where there is a huge glitch.
// Huge glitches are harmful because they affect the correct estimation 
// of the selection threshold and could mask the nearby events at lower SNR. 
// Warning : events with high SNR can be rejected by this procedure (see SETHR)

#define SETHR	1000000	// Is the threshold which define the time slices to be cutted
                        // Warning : this value depends on the frequency interval [fHigh:fLow] 

#define TWIN    0.5     // Is the window (sec) used to integrate the energies in time
                        // TWIN must be a multiple of the greatest time resolution used in the analysis 

#define TEDG	1.5	// Is the time window (sec) to be cutted when the time integrated energy is > SETHR

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// this plugin is called in cwb2G.cc to the whitened data at the beginning of the COHERENT stage

//  cout << endl;
//  cout << "-----> CWB_Plugin_Gating.C" << endl;
//  cout << "ifo " << ifo.Data() << endl;
//  cout << "type " << type << endl;
//  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {
    cout << endl;
    cout << "-----> CWB_Plugin_Gating.C" << endl;
    cout << endl;
  }

  if(type==CWB_PLUGIN_ICOHERENCE) {

    // get data range
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
    double Tb = gCWB2G->GetSegBegin()-cfg->segEdge;
    double Te = gCWB2G->GetSegEnd()+cfg->segEdge;


    //cout << "CWB_Plugin_Gating.C - " << "Tb : " << int(Tb) << " Te : " << int(Te) << endl; 

    int nIFO = (gCWB2G->IsSingleDetector()) ? 1 : net->ifoListSize();   // number of detectors
    int    size = net->getifo(0)->getHoT()->size();    			// number of time bins
    double rate = net->getifo(0)->getHoT()->rate();    			// rate
    double dt = 1./rate;              					// time bin resolution (sec)

    //cout << "-----> CWB_Plugin_Gating.C - size : " << size << "\t rate : " << rate << "\t dt : " << dt << " length(sec) " << size*dt << endl;                            

    wavearray<double>* hot[NIFO_MAX];	
    for(int n=0; n<nIFO; n++) hot[n] = net->getifo(n)->getHoT();	// whitened data

    // A new array 'SE' is obtained from the arrays 'hot whitened data'
    // It contains the time integrated energies over the time window TWIN
    // A time sample is marked as to be vetoed when the energy SE>SETHR
    // When a time sample 'i' is vetoed also the nearby M time size are vetoed [i-M:i+M]
    // where M is the number of samples contained in TEDG sec 

    int X = int(cfg->segEdge/dt+0.5); 	// samples in segEdge : this the is scratch time
    int M = int(TEDG/dt)+1; 		// samples in TEDG sec
    int N = int(TWIN/dt);		// number of time samples in TWIN
    if(N<1) N=1;
    vector<waveSegment> glist;		// gating segment list
    wavearray<double> SE(size);
    wavearray<int> sVeto(size);		// array which contains samples to be vetoed

    int gsize=0;
    for(int n=0; n<nIFO; n++) {		// loop over detectors
      SE=0;
//      for(int i=N-1;i<size;i++) for(int j=0;j<N;j++) SE[i]+=pow(hot[n]->data[i-j],2);

      for(int j=0;j<N;j++) SE[N-1]+=pow(hot[n]->data[N-1-j],2);
      for(int i=N;i<size;i++) {
        SE[i]  = SE[i-1];
        SE[i] -= pow(hot[n]->data[i-N],2);
        SE[i] += pow(hot[n]->data[i],2);
      }

      sVeto=0;
      for(int i=0;i<size;i++) {
        if(SE[i]>SETHR) { 		
          int is = i-M>X ? i-M : X; 
          int es = i+M<(size-X) ? i+M : size-X; 
          for(int k=is;k<es;k++) sVeto[k]=1;
        }
      }
      // array derivative -> gating_start=1, gating_stop=-1
      sVeto[0]=0;sVeto[size-1]=0;
      for(int i=0;i<size-1;i++) sVeto[i]=sVeto[i+1]-sVeto[i];

      // build the gating segment list
      waveSegment gseg;					// gating segment
      gseg.index = 0;
      for(int i=0;i<size;i++) {
        if(sVeto[i]== 1) gseg.start = Tb+int(i*dt);	// round to nearest lower integer
        if(sVeto[i]==-1) {
          gseg.stop  = Tb+int(i*dt+0.5);		// round to nearest higher integer
          gseg.index++;
          glist.push_back(gseg);
        }
      } 

      if(glist.size()>0) {
        cout << endl;
        cout.precision(10);
        for(int j=gsize;j<glist.size();j++) {
          cout << j << " -----> CWB_Plugin_Gating.C - " << net->getifo(n)->Name << " gating time : start="
               << glist[j].start << " stop=" << glist[j].stop << " len=" << glist[j].stop-glist[j].start << endl;
        }
        cout << endl;
        gsize=glist.size();
      }
    }                                                                                               

    double gating_time = 0;						// total gating time
    if(glist.size()) { 
      glist = CWB::Toolbox::unionSegments(glist);  			// Join & sort a waveSegment list
      gating_time = CWB::Toolbox::getTimeSegList(glist);		// get total gating time
      glist = CWB::Toolbox::invertSegments(glist);  			// Invert waveSegment list
      net->segList = CWB::Toolbox::mergeSegLists(glist,net->segList);  	// merge gating veto list to the net->segList
    }

    cout<< "-----> CWB_Plugin_Gating.C - gating time : " << gating_time << endl << endl; 

    // add infos to history
    char info[256];
    sprintf(info,"-GT:%g",gating_time);
    gCWB2G->PrintAnalysisInfo(CWB_STAGE_COHERENCE,"cwb2G::Coherence",info,false);
  }

  return;
}
