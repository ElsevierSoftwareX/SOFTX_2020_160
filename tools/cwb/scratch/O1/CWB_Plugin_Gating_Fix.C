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

// this plugin is called in cwb2G.cc after the production of the TF maps with max over the sky energy (TFmax)
// The procedure is applied for each detector and for each resolution level using the TFmax.

  cout << endl;
  cout << "-----> CWB_Plugin_Gating_Fix.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_ECOHERENCE) {

    // WSeries contains energy (it is stored in the 0 phase amplitudes)

    // import resolution level
    size_t gILEVEL=-1; IMPORT(size_t,gILEVEL)
    char slevel[256];sprintf(slevel,",%d,",gILEVEL); 

    int nIFO = net->ifoListSize();
    detector* pD[NIFO_MAX];             // pointers to detectors
    for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);
    WSeries<double>* WS[NIFO_MAX];
    for(int n=0; n<nIFO; n++) WS[n] = pD[n]->getTFmap();

    int layers = WS[0]->maxLayer()+1;  	// numbers of frequency bins (first & last bins have df/2)
    int slices = WS[0]->sizeZero();    	// number of time bins

    double df = WS[0]->resolution();    // frequency bin resolution (hz)	//FIX
    double dt = 1./(2*df);              // time bin resolution (sec)		//FIX

    int rate = int(1./dt);                                                                           

    cout << "layers : " << layers << "\t slices : " << slices << "\t rate : " << rate
         << "\t dt : " << dt << "\t df : " << df << endl;                            

    double rTWIN = fabs(TWIN/dt-TMath::Nint(TWIN/dt));
    if(TWIN<dt || rTWIN>1e-7) {	// FIX!!!	fix precision !!!
      cout << "-----> CWB_Plugin_Gating_Fix.C : Error-> " << " TWIN=" << TWIN 
           << " is not a multiple of dt=" << dt << endl;
      gSystem->Exit(1);
    }      

    double rTEDG = fabs(TEDG/dt-TMath::Nint(TEDG/dt));
    if(TEDG<dt || rTEDG>1e-7) {	// FIX!!!       fix precision !!!
      cout << "-----> CWB_Plugin_Gating_Fix.C : Error-> " << " TEDG=" << TEDG 
           << " is not a multiple of dt=" << dt << endl;
      gSystem->Exit(1);
    }      

    // For each time slice (time index) we compute the sum of the pixel 
    // energy over all frequency layers (freq index)
    // The sum is stored in the array 'se[NIFO_MAX]' 
    // se = projection of TF map energy on the time axis
    wavearray<double> se[NIFO_MAX];	
    for(int n=0; n<nIFO; n++) {
      se[n].resize(slices); se[n]=0;
      for(int i=0;i<slices;i++) {
        for(int j=0;j<layers;j++) se[n][i] += WS[n]->getSample(i,j);
      }                                                                                                     
    }              

    // A new array 'SE[NIFO_MAX]' is obtained from the arrays 'se[NIFO_MAX]'
    // It contains the time integrated energies over the time window TWIN
    // NOTE : this procedure makes the energies independents from the resolution level
    int N = int(TWIN/dt);		// number of time samples in TWIN
    if(N<1) N=1;
    wavearray<double> SE[NIFO_MAX];
    for(int n=0; n<nIFO; n++) {
      SE[n].resize(slices); SE[n]=0;
      for(int i=N-1;i<slices;i++) for(int j=0;j<N;j++) SE[n][i]+=se[n][i-j];
    }                                                                                               

    // find the time slices to be excluded :
    // 1) from the computation of the pixel's selection threshold
    // 2) from the pixel's selection
    // A time slice is marked as to be cutted when the energy SE>SETHR
    // When a time slice is cutted also the nearby M time slices are cutted [slice-M:slice+M]
    // where M is the number of samples contained in TEDG sec 
    int X = int(cfg->segEdge/dt+0.5); 	// samples in segEdge : this the is scratch time
    int M = int(TEDG/dt)+1; 		// samples in TEDG sec
    wavearray<int> sCut[NIFO_MAX];	// array which contains the slices to be cut
    for(int n=0; n<nIFO; n++) {
      sCut[n].resize(slices); sCut[n]=0;
      for(int i=0;i<slices;i++) {
        if(SE[n][i]>SETHR) { 		
          int is = i-M>X ? i-M : X; 
          int es = i+M<(slices-X) ? i+M : slices-X; 
          for(int ii=is;ii<es;ii++) sCut[n][ii]=1;
        }
        //cout << i << " " << SE[n][i] << endl;  
      }
    }

    // All pixels in the layers which time slice is marked as cutted are filled with a negative energy = -1e20
    // -> these pixels do not partecipate to the computation of the pixel's selection thereshold
    //    and moreover these pixels are excluded from the selection
    double gating_time=0.;	// total time vetoed by gating
    for(int n=0; n<nIFO; n++) {
      for(int i=0;i<slices;i++) {	// set to -1e20 the energy
        if(sCut[n][i]) {
          gating_time+=dt;
          for(int j=0;j<layers;j++) WS[n]->putSample(-1e20,i,j);
        }
      }
    } 

    // add infos to history
    char info[256];
    sprintf(info,"-RES:%d-GT:%g",cfg->l_high-gILEVEL,gating_time);
    cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G)
    gCWB2G->PrintAnalysisInfo(CWB_STAGE_COHERENCE,"cwb2G::Coherence",info,false);
  }

  return;
}
