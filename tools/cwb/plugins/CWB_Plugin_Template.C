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
#include "TComplex.h"
#include "TMath.h"
#include "TSystem.h"
#include "mdc.hh"
#include "WDM.hh"
#include "regression.hh"
#include "frame.hh"
#include "watplot.hh"
#include "Biorthogonal.hh"
#include <fstream>
#include <vector>


using namespace CWB;

//______________________________________________________________________________
void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// Plugin Template
//
//
// jfile - pointer to job file
//         this file is used by pipeline to store temporary informations like noise,mdc data
//
//   cfg - pointer to the user configuration setup (see config.hh file)
//         it can be used to read or change setting 
//
//   net - pointer to network object (see network.hh file)
//         this is the class which perform tha main anaalysis
//         this object contains the detector definitions  
//
//     x - this arrays contains the ifo noise or mdc data (see wseries.hh/wavearray.hh files)  
//
//   ifo - detector name (Ex: L1,H1,V1)
//
//  type - is the stage of the cwb analysis  
//	   this is the list of stage types 
// 
//         CWB_PLUGIN_CONFIG             = 0
//         CWB_PLUGIN_NETWORK            = 1
//         CWB_PLUGIN_INIT_JOB           = 9
//         CWB_PLUGIN_DATA               = 2
//         CWB_PLUGIN_MDC                = 3
//         CWB_PLUGIN_RMDC               = 11
//         CWB_PLUGIN_OREADDATA          = 14
//         CWB_PLUGIN_IDATA_CONDITIONING = 4
//         CWB_PLUGIN_DATA_CONDITIONING  = 15
//         CWB_PLUGIN_ODATA_CONDITIONING = 5
//         CWB_PLUGIN_ICOHERENCE         = 6
//         CWB_PLUGIN_XCOHERENCE         = 16
//         CWB_PLUGIN_OCOHERENCE         = 13
//         CWB_PLUGIN_ISUPERCLUSTER      = 7
//         CWB_PLUGIN_OSUPERCLUSTER      = 17
//         CWB_PLUGIN_ILIKELIHOOD        = 12
//         CWB_PLUGIN_OLIKELIHOOD        = 8
//         CWB_PLUGIN_CLOSE_JOB          = 10
//


  cout << endl;
  cout << "-----> CWB_Plugin_Template.C : " << ifo.Data() << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {
    //
    // at this stage the config contains the configuration 
    // user can read or change the configuration parameters 
    //
    // jfile : NULL
    // cfg   : configured
    // net   : not configured
    // x     : NULL
    // ifo   : ""
    // 
    // Examples
    //
    // cfg->Print();	    // print configuration
    // cfg->bpp = 0.0001;   // change probability for black pixel selection
  }

  if(type==CWB_PLUGIN_NETWORK) {
    //
    // at this stage the network object has been created
    // Ex: user can read/change the detectors objects 
    //
    // jfile : NULL
    // cfg   : configured
    // net   : configured
    // x     : NULL
    // ifo   : ""
    // 
    // Examples
    //
    // CWB_Plugin_SGW.C : Set Scalar Gravitational Wave Model
  }

  if(type==CWB_PLUGIN_INIT_JOB) {
    //
    // at this stage the data quality list has been read
    // the job segment is computed according the DQ and run ID
    // this stage is performed for each detector, ifo parameter
    //
    // jfile : NULL
    // cfg   : configured
    // net   : configured
    // x     : is the pointer to a dummy wavearray 
    //         used only to get x.rate(),x.start(),x.size()
    // ifo   : ifo label (L1,H1,V1,...)
    // 
  }

  if(type==CWB_PLUGIN_DATA) {
    //
    // at this stage the detector strain has been read and stored in x array
    // user can modify or inject its own noise (Ex: on the fly noise)
    // user can read/change noise data (Ex: Calibration)
    //
    // jfile : opened in update mode
    //         When noise is read it is stored to job file :
    //         TDirectoryFile*               strain  
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //
    // cfg   : configured
    // net   : configured
    // x     : is the pointer to noise data wavearray of detector ifo
    // ifo   : ifo label (L1,H1,V1,...)
    // 
    // Examples
    //
    // CWB_Plugin_SimNoise.C : generate simulated gaussian noise
    //
    // for(int i=0;i<x->size();i++) x->data[i]=0;  // access to ifo noise data array
  }

  if(type==CWB_PLUGIN_MDC) {
    //
    // at this stage the mdc has been read and stored in x array
    // user can modify or inject its own mdc (Ex: on the fly mdc)
    //
    // jfile : opened in update mode
    //         When mdc is read it is stored to job file :
    //         TDirectoryFile*               mdc    
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //
    // cfg   : configured
    // net   : configured
    // x     : is the pointer to MDC data wavearray of detector ifo
    // ifo   : ifo label (L1,H1,V1,...)
    // 
    // Examples
    //
    // CWB_Plugin_MDC_OTF.C : injected MDC 'on the fly'
    //
    // for(int i=0;i<x->size();i++) x->data[i]=0;  // access to mdc ifo data array
  }

  if(type==CWB_PLUGIN_RMDC) {
    //
    // at this stage the mdc has been read and stored in x array
    // user can modify or inject its own mdc (Ex: on the fly mdc)
    //
    // jfile : opened in update mode
    //         job file contains :
    //         TDirectoryFile*               strain 
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //         TDirectoryFile*               mdc   
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //
    // cfg   : configured
    // net   : configured
    // x     : is the pointer to MDC data wavearray of detector ifo
    // ifo   : ifo label (L1,H1,V1,...)
    // 
  }

  if(type==CWB_PLUGIN_OREADDATA) {
    //
    // at this stage the mdc has been read 
    // this plugin is called only if simulation=2
    // when simulation=2 (fixed input SNR) the MDC are rescaled to snr network = 1
    // user can modify or inject its own mdc (Ex: on the fly mdc)
    // 
    // jfile : opened in update mode
    //         job file contains :
    //         TDirectoryFile*               strain 
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //         TDirectoryFile*               mdc   
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //
    // cfg   : configured
    // net   : configured
    //         - MDC rescaled to snr network = 1 : stored in net->getifo(ifoID)->HoT
    // x     : NULL
    // ifo   : ""
    // 
    // Examples
    //
    // CWB_Plugin_SNR.C : compute and save the mdc SNR
  }

  if(type==CWB_PLUGIN_IDATA_CONDITIONING) {
    //
    // called at the begining of DATA_CONDITIONING stage
    // at this stage the rescaled mdc is added to noise data and stored in x array
    //
    // jfile : opened in update mode
    //         job file contains :
    //         TDirectoryFile*               strain
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //         TDirectoryFile*               mdc  
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //         
    //         if analysis is done in 2 stages job file contains also :
    //         TDirectoryFile*               waveform                            
    //          TDirectoryFile*              waveform-f0                      
    //           KEY: detector       L1;1                                                
    //           KEY: detector       H1;1                                                
    //           KEY: detector       V1;1                                                
    //
    // cfg   : configured
    // net   : configured
    // x     : is the pointer to MDC+noise data wavearray of detector ifo
    // ifo   : ifo label (L1,H1,V1,...)
    // 
    // Examples
    //
    // CWB_Plugin_MakeScalogram.C : produce scalograms of input data
  }

  if(type==CWB_PLUGIN_DATA_CONDITIONING) {
    //
    // at this stage the rescaled mdc is added to noise data and stored in x array
    // this stage is executed only if cfg.dcPlugin=true
    // if cfg.dcPlugin=true : user can apply data conditioning (Ex: regression analysis)
    //
    // jfile : opened in update mode
    //         job file contains :
    //         TDirectoryFile*               strain
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //         TDirectoryFile*               mdc  
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //
    // cfg   : configured
    // net   : configured
    // x     : is the pointer to MDC+noise data wavearray of detector ifo
    // ifo   : ifo label (L1,H1,V1,...)
    // 
    // Examples
    //
    // CWB_Plugin_DataConditioning.C : shows how to implement a custom Data Conditioning 
  }

  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {
    //
    // called at the end of DATA_CONDITIONING stage
    // at this stage the noise+mdc data is whitened and stored in x array
    // Ex : user can monitor if whitened is correct
    //
    // jfile : opened in update mode
    //         job file contains :
    //         TDirectoryFile*               strain 
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //         TDirectoryFile*               mdc  
    //          KEY: wavearray<double>       L1;1
    //          KEY: wavearray<double>       H1;1
    //          KEY: wavearray<double>       V1;1
    //         
    //         if analysis is done in 2 stages job file contains also :
    //         TDirectoryFile*               waveform                            
    //          TDirectoryFile*              waveform-f0                      
    //           KEY: detector       L1;1                                                
    //           KEY: detector       H1;1                                                
    //           KEY: detector       V1;1                                                
    //           .. for each factor
    //
    // cfg   : configured
    // net   : configured
    // x     : is the pointer to MDC+noise whitened data wavearray of detector ifo
    // ifo   : ifo label (L1,H1,V1,...)
    // 
    // Examples
    //
    // CWB_Plugin_WDM_freqCuts.C : Implements the 2G frequency cuts after the whitening
    // CWB_Plugin_MakeWhiteFrame.C : produces frames of whitened data  
  }

  if(type==CWB_PLUGIN_ICOHERENCE) {
    //
    // called at the begining of COHERENCE stage
    // at this stage the noise+mdc data is whitened and stored in x array
    // Ex : user can monitor if whitened is correct
    //
    // jfile : opened in update mode
    //         
    //         if analysis is done in 2 stages job file contains also :
    //         TDirectoryFile*               waveform                            
    //
    // cfg   : configured
    // net   : configured
    // x     : NULL
    // ifo   : factor index
    // 
  }

  if(type==CWB_PLUGIN_XCOHERENCE) {
    //
    // called at the inside the resolution level loop in the COHERENCE stage
    //
    // to get the factor the variable gIFACTOR must be imported
    // int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    // to get the resolution level the variable gILEVEL must be imported
    // size_t gILEVEL=-1; IMPORT(size_t,gILEVEL)
    // 
    // at this stage the max energy TF has been computed
    // max energy TF are stored in net->getifo(ifoID)->getTFmap()
    // Ex : user can monitor the max energy TF
    //
    // jfile : opened in update mode
    //         TDirectoryFile*               coherence     
    //          OBJ: TTree   clusters-cycle:0        
    //         
    //         if analysis is done in 2 stages job file contains also :
    //         TDirectoryFile*               waveform                            
    //
    // cfg   : configured
    // net   : configured
    //         - TF max energy is stored net->getifo(ifoID)->getTFmap()
    // x     : NULL
    // ifo   : factor index
    // 
    // Examples
    //
    // CWB_Plugin_fCuts.C : Implements the 2G frequency cuts in pixel selection stage
  }

  if(type==CWB_PLUGIN_OCOHERENCE) {
    //
    // called at the end of COHERENCE stage
    // at this stage the coherence analysis has been done
    // selected pixels are stored in jfile
    // user can read/modify pixels
    //
    // jfile : opened in update mode
    //         TDirectoryFile*               coherence     
    //          OBJ: TTree   clusters-cycle:0        
    //         
    //         if analysis is done in 2 stages job file contains also :
    //         TDirectoryFile*               waveform                            
    //
    // cfg   : configured
    // net   : configured 
    //         - clusters are stored in the net 
    // x     : NULL
    // ifo   : factor index
    // 
  }

  if(type==CWB_PLUGIN_ISUPERCLUSTER) {
    //
    // at this stage the supercluster analysis started
    // 
  }

  if(type==CWB_PLUGIN_OSUPERCLUSTER) {
    //
    // at this stage the supercluster analysis has been done and stored in net object
    // user can read clusters
    // 
    // jfile : opened in update mode
    //         TDirectoryFile*               coherence      
    //          KEY: TTree                   clusters-cycle:0;1      
    //         TDirectoryFile*               supercluster  
    //          OBJ: TTree                   clusters-cycle:0        
    //         
    //         if analysis is done in 2 stages job file contains also :
    //         TDirectoryFile*               waveform                            
    //
    // cfg   : configured
    // net   : configured
    //         - clusters are stored in the net 
    // x     : NULL
    // ifo   : ""
    // 
  }

  if(type==CWB_PLUGIN_ILIKELIHOOD) {
    //
    // this stage is the bebinning of likelihood analysis 
    // sparse maps are loaded
    // reconstructed events are stored net object
    // user can extract reconstructed events (waveform, params) and compute parameters estimation
    // 
    // jfile : opened in update mode
    //         TDirectoryFile*               sparse  
    //          KEY: SSeries<double> L1-level:0:8;1
    //          KEY: SSeries<double> H1-level:0:8;1
    //          KEY: SSeries<double> V1-level:0:8;1
    //          .. for all resolution levels 
    //         KEY: TDirectoryFile   supercluster;1 
    //          KEY: TTree   clusters-cycle:0;1    
    //
    //         if analysis is done in 2 stages job file contains also :
    //         KEY: TDirectoryFile   waveform;1                               
    //         KEY: CWB::config      config;1                                            
    //         KEY: network          network;1                                                   
    //         KEY: CWB::History     history;1                                           
    //         KEY: cwb2G            cwb;1                                                       
    //
    // cfg   : configured
    // net   : configured
    //         - sparse map is stored net->getifo(ifoID)->vSS[ifoID]
    //         - clusters are stored in the net 
    // x     : NULL
    // ifo   : factor-id or lag id
    // 
    // Examples
    //
    // CWB_Plugin_cwb_inet.C : plot sparse TF maps
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {
    //
    // at this stage the likelihood analysis has been done
    // reconstructed events are stored net object
    // user can extract reconstructed events (waveform, params) and compute parameters estimation
    // user can store custom results to jfile 
    // 
    // jfile : opened in update mode
    //         TDirectoryFile*               sparse  
    //          KEY: SSeries<double> L1-level:0:8;1
    //          KEY: SSeries<double> H1-level:0:8;1
    //          KEY: SSeries<double> V1-level:0:8;1
    //          .. all resolution levels 
    //         KEY: TDirectoryFile   supercluster;1 
    //          KEY: TTree   clusters-cycle:0;1    0
    //
    //         if analysis is done in 2 stages job file contains also :
    //         KEY: TDirectoryFile   waveform;1                                
    //         KEY: CWB::config      config;1                                            
    //         KEY: network          network;1                                                   
    //         KEY: CWB::History     history;1                                           
    //         KEY: cwb2G            cwb;1                                                       
    //
    // cfg   : configured
    // net   : configured
    //         - reconstructed clusters are contained in net object net->getwc(lag)
    // x     : NULL
    // ifo   : ""
    // 
    // Examples
    //
    // CWB_Plugin_NN.C : includes netcluster event stricture into the output wave root file 
    // CWB_Plugin_SkyProb.C : save SkyProb into the output wave root file
    // CWB_Plugin_pixeLHood.C : dump/plot the likelihood/null pixels of the detected/reconstructed event
  }

  if(type==CWB_PLUGIN_CLOSE_JOB) {
    //
    // at this stage job is finished and 
    //
    // jfile : opened in update mode
    //         if analysis is done in 2 stages job file contains also :
    //         This is the trigger file used as input for the final stage analysis
    //         KEY: TDirectoryFile   waveform;1                                 
    //         KEY: TDirectoryFile   supercluster;1                         
    //         KEY: TDirectoryFile   sparse;1                                    
    //         KEY: CWB::config      config;1                                            
    //         KEY: network          network;1                                                   
    //         KEY: CWB::History     history;1                                           
    //         KEY: cwb2G            cwb;1                                                       
    //
    // cfg   : configured
    // net   : configured
    // x     : is the pointer to a dummy wavearray 
    //         used only to get x.rate(),x.start(),x.size()
    // ifo   : ifo label (L1,H1,V1,...)
    // 
  }

  return;
}
