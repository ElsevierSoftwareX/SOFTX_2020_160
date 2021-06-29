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
#include "mdc.hh"
#include "watplot.hh"
#include "gwavearray.hh"
#include <vector>

//#define SAVE_MRA_PLOT			// enable event MRA plots 
//#define SAVE_WHT_PLOT			// enable event WHITE plots 
//#define SAVE_STR_PLOT			// enable event STRAIN plots 
//#define SAVE_TF_PLOT			// enable event WHITE TF plots 

#define MDC_OTF				// enable MDC On The Fly
#define NOISE_OTF			// enable NOISE On The Fly

//#define SAVE_EVT_CLUSTER		// save cluster to the event output tree		
//#define SAVE_EVT_SNR			// save iSNR,oSNR,ioSNR 

void ComputeResidualEnergyOptTF(wavearray<double>* wfINJ, wavearray<double>* wfREC, 
                                network* NET, netevent* EVT, int lag, int ID, int ifoID);
void ComputeResidualEnergy(wavearray<double>* wfINJ, wavearray<double>* wfREC, 
                           double& enINJ, double& enREC, double& xcorINJ_REC);
double ComputeEnergy(WSeries<double> *WS);
double ComputeResidualEnergy(WSeries<double> *WS1, WSeries<double> *WS2);
void ApplyFrequencyCuts(wavearray<double>* wf, CWB::config* cfg);
void PlotWaveform(TString ifo, wavearray<double>* wfINJ, wavearray<double>* wfREC,
                  CWB::config* cfg, bool fft=false, bool strain=false);
void PlotMRA(netcluster* pwc, int nIFO, int ID);


void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!NOISE_MDC_SIMULATION
// Extract whitened/strain injected & reconstructed waveforms, compute residual energy 
// Save residual energy to the output wave root file 
// This plugin can be use for Waveform Reconstruction Challenge (WRC) studies

  cout << endl;
  cout << "-----> CWB_Plugin_WRC.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
#ifdef NOISE_OTF
    cfg->dataPlugin=true; 				// disable read data from frames
#endif
#ifdef MDC_OTF
    cfg->mdcPlugin=true;  				// disable read mdc from frames
#endif
    cfg->outPlugin=true;  				// disable built-in output root file
  }


#ifdef NOISE_OTF
  if(type==CWB_PLUGIN_DATA) {  

    CWB::Toolbox TB;

    int seed;
    if(ifo.CompareTo("L1")==0) seed=1000;
    if(ifo.CompareTo("H1")==0) seed=2000;
    if(ifo.CompareTo("V1")==0) seed=3000;
    if(ifo.CompareTo("J1")==0) seed=4000;
    if(ifo.CompareTo("A2")==0) seed=5000;
    if(ifo.CompareTo("Y2")==0) seed=6000;
    if(ifo.CompareTo("Y3")==0) seed=7000;

    TString fName;
    if(ifo.CompareTo("L1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("V1")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";
    if(ifo.CompareTo("J1")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";            
    if(ifo.CompareTo("A2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";            
    if(ifo.CompareTo("Y2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";            
    if(ifo.CompareTo("Y3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";            

    int size=x->size();
    double start=x->start();
    TB.getSimNoise(*x, fName, seed, NET->nRun);
    x->resize(size);                           
    x->start(start);                           
  }                                            
#endif

#ifdef MDC_OTF
  if(type==CWB_PLUGIN_MDC) {  

    char cmd[128];
    sprintf(cmd,"network* net = (network*)%p;",NET);
    gROOT->ProcessLine(cmd);

    CWB::mdc MDC(NET);

    // ---------------------------------
    // read plugin config 
    // ---------------------------------

    cfg->configPlugin.Exec();

    // ---------------------------------
    // set list of mdc waveforms
    // ---------------------------------
   
    IMPORT(CWB::mdc,MDC) 
    MDC.Print();

    // ---------------------------------
    // get mdc data
    // ---------------------------------

    MDC.Get(*x,ifo);

    // ---------------------------------
    // set mdc list in the network class 
    // ---------------------------------

    if(ifo.CompareTo(NET->ifoName[0])==0) {
      NET->mdcList.clear();
      NET->mdcType.clear();
      NET->mdcTime.clear();
      NET->mdcList=MDC.mdcList;
      NET->mdcType=MDC.mdcType;
      NET->mdcTime=MDC.mdcTime;
    }

    cout.precision(14);
    for(int k=0;k<(int)NET->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)NET->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)NET->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;
  }
#endif

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_WRC.C -> CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_WRC.C -> " << " gIFACTOR : " << gIFACTOR << endl;

    // import slagShift
    float* gSLAGSHIFT=NULL; IMPORT(float*,gSLAGSHIFT)

    int nIFO = NET->ifoListSize();			// number of detectors
    int K = NET->nLag;  				// number of time lag          
    netevent* EVT;
    wavearray<double> id;
    double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];                 
    int rate = 0;					// select all resolutions

    // search output root file in the system list
    TFile* froot = NULL;                         
    TList *files = (TList*)gROOT->GetListOfFiles();
    TString outDump="";
    if (files) {                                   
      TIter next(files);                           
      TSystemFile *file;                           
      TString fname;                               
      bool check=false;                            
      while ((file=(TSystemFile*)next())) {        
         fname = file->GetName();                  
         // set output root file as the current file
         if(fname.Contains("wave_")) {
           froot=(TFile*)file;froot->cd();
           outDump=fname;
           outDump.ReplaceAll(".root.tmp",".txt");
           cout << "output file name : " << fname << endl;
         }
      }                                                               
      if(!froot) {                                                    
        cout << "CWB_Plugin_WRC.C : Error - output root file not found" << endl;
        gSystem->Exit(1);                                                                             
      }                                                                                      
    } else {                                                                                 
      cout << "CWB_Plugin_WRC.C : Error - output root file not found" << endl;  
      gSystem->Exit(1);                                                                               
    }                                                                                        

    netcluster* pwc = new netcluster();  
    char ciSNR[16]; sprintf(ciSNR,"_iSNR[%1d]/F",nIFO);
    char coSNR[16]; sprintf(coSNR,"_oSNR[%1d]/F",nIFO);
    char cioSNR[16]; sprintf(cioSNR,"_ioSNR[%1d]/F",nIFO);
    float* iSNR  = new float[nIFO];
    float* oSNR  = new float[nIFO];
    float* ioSNR = new float[nIFO];

    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree!=NULL) {
      EVT = new netevent(net_tree,nIFO);
#ifdef SAVE_EVT_CLUSTER
      net_tree->SetBranchAddress("netcluster",&pwc); 
#endif
#ifdef SAVE_EVT_SNR
      net_tree->SetBranchAddress("_iSNR",iSNR);
      net_tree->SetBranchAddress("_oSNR",oSNR);
      net_tree->SetBranchAddress("_ioSNR",ioSNR);
#endif
    } else {
      EVT = new netevent(nIFO);
      net_tree = EVT->setTree();
#ifdef SAVE_EVT_CLUSTER
      // add new netcluster leaf
      net_tree->Branch("netcluster","netcluster",&pwc,32000,0); 
#endif
#ifdef SAVE_EVT_SNR
      net_tree->Branch("_iSNR",iSNR,ciSNR);
      net_tree->Branch("_oSNR",oSNR,coSNR);
      net_tree->Branch("_ioSNR",ioSNR,cioSNR);
#endif
    }
    EVT->setSLags(gSLAGSHIFT);                          // set slags into netevent

    for(int k=0; k<K; k++) {  				// loop over the lags

      id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);

      for(int j=0; j<(int)id.size(); j++) {  		// loop over cluster index

        int ID = size_t(id.data[j]+0.5);

        if(NET->getwc(k)->sCuts[ID-1]!=-1) continue;    // skip rejected/processed clusters

        double ofactor=0;
        if(cfg->simulation==4)      ofactor=-factor;
        else if(cfg->simulation==3) ofactor=-gIFACTOR;
        else                        ofactor=factor;

        if(k==0) {					// only for zero lag

          EVT->output2G(NULL,NET,ID,k,ofactor);		// get reconstructed parameters

          double recTime = EVT->time[0];		// rec event time det=0
          double injTh   = EVT->theta[1];		// inj theta
          double injPh   = EVT->phi[1];			// inj phi
          double recTh   = EVT->theta[0];		// rec theta
          double recPh   = EVT->phi[0];			// rec phi

          injection INJ(nIFO);
          // get inj ID
          int M = NET->mdc__IDSize();
          double injTime = 1.e12;
          int injID   = -1;
          for(int m=0; m<M; m++) {
             int mdcID = NET->getmdc__ID(m);		
             double T = fabs(recTime - NET->getmdcTime(mdcID));
             if(T<injTime && INJ.fill_in(NET,mdcID)) {
                injTime = T;
                injID = mdcID;
             }
          }

          if(INJ.fill_in(NET,injID)) {                 // get inj parameters

             wavearray<double>** pwfINJ = new wavearray<double>*[nIFO];
             wavearray<double>** pwfREC = new wavearray<double>*[nIFO];
             detector* pd = NET->getifo(0);
             int idSize = pd->RWFID.size();
             int wfIndex=-1;
             for(int mm=0; mm<idSize; mm++) if (pd->RWFID[mm]==ID) wfIndex=mm;

             // extract whitened injected & reconstructed waveforms
             for(int n=0; n<nIFO; n++) {

                pd = NET->getifo(n);

                pwfINJ[n] = INJ.pwf[n];
                if (pwfINJ[n]==NULL) {
                   cout << "CWB_Plugin_WRC.C : Error : Injected waveform not saved !!! : detector "
                        << NET->ifoName[n] << endl;
                   continue;
                }
                if (wfIndex<0) {
                   cout << "CWB_Plugin_WRC.C : Error : Reconstructed waveform not saved !!! : ID -> "
                        << ID << " : detector " << NET->ifoName[n] << endl;
                   continue;
                }

                if (wfIndex>=0) pwfREC[n] = pd->RWFP[wfIndex];
                double R = pd->TFmap.rate();
                double rFactor = 1.;
                rFactor *= factor;
                wavearray<double>* wfINJ = pwfINJ[n];	// array of injected waveforms
                *wfINJ*=rFactor;			// injected wf is multiplied by the custom factor
                wavearray<double>* wfREC = pwfREC[n];	// array of reconstructed waveforms

                double enINJ, enREC, xcorINJ_REC;
                // compute residual energy in time domain
                ComputeResidualEnergy(wfINJ, wfREC, enINJ, enREC, xcorINJ_REC);
                // compute residual energy in TF domain at optimal resolution
                ComputeResidualEnergyOptTF(wfINJ, wfREC, NET, EVT, k, ID, n);
#ifdef SAVE_WHT_PLOT
                PlotWaveform(NET->ifoName[n], wfINJ, wfREC, cfg, false, false);
                PlotWaveform(NET->ifoName[n], wfINJ, wfREC, cfg, true, false);
#endif

                iSNR[n]    = enINJ;			// inj whitened SNR^2 [fLow,fHigh]
                oSNR[n]    = enREC;			// rec whitened SNR^2 [fLow,fHigh]
                ioSNR[n]   = xcorINJ_REC;		// inj*rec whitened energy [fLow,fHigh]

                *wfINJ*=1./rFactor;			// restore injected amplitude
             }
             delete [] pwfINJ;
             delete [] pwfREC;
          }

          // extract strain injected & reconstructed waveforms
          if(INJ.fill_in(NET,-(injID+1))) {                 // set injections

             wavearray<double>** pwfINJ = new wavearray<double>*[nIFO];
             wavearray<double>** pwfREC = new wavearray<double>*[nIFO];
             detector* pd = NET->getifo(0);
             int idSize = pd->RWFID.size();
             int wfIndex=-1;
             for(int mm=0; mm<idSize; mm++) if (pd->RWFID[mm]==-ID) wfIndex=mm;

             for(int n=0; n<nIFO; n++) {

                pd = NET->getifo(n);

                pwfINJ[n] = INJ.pwf[n];
                if (pwfINJ[n]==NULL) {
                   cout << "CWB_Plugin_WRC.C : Error : Injected waveform not saved !!! : detector "
                        << NET->ifoName[n] << endl;
                   continue;
                }
                if (wfIndex<0) {
                   cout << "CWB_Plugin_WRC.C : Error : Reconstructed waveform not saved !!! : ID -> "
                        << ID << " : detector " << NET->ifoName[n] << endl;
                   continue;
                }

                if (wfIndex>=0) pwfREC[n] = pd->RWFP[wfIndex];
                double R = pd->TFmap.rate();
                double rFactor = 1.;
                rFactor *= factor;
                wavearray<double>* wfINJ = pwfINJ[n];	// array of injected waveforms
                *wfINJ*=rFactor;			// injected wf is multiplied by the custom factor
                wavearray<double>* wfREC = pwfREC[n];	// array of reconstructed waveforms

                ApplyFrequencyCuts(wfINJ, cfg);
#ifdef SAVE_STR_PLOT
                PlotWaveform(NET->ifoName[n], wfINJ, wfREC, cfg, false, true);
                PlotWaveform(NET->ifoName[n], wfINJ, wfREC, cfg, true, true);
#endif

                double enINJ, enREC, xcorINJ_REC;
                // compute residual energy in time domain
                ComputeResidualEnergy(wfINJ, wfREC, enINJ, enREC, xcorINJ_REC);
                // compute residual energy in TF domain at optimal resolution
                ComputeResidualEnergyOptTF(wfINJ, wfREC, NET, EVT, k, ID, n);

                iSNR[n]    = enINJ;			// inj whitened SNR^2 [fLow,fHigh]
                oSNR[n]    = enREC;			// rec whitened SNR^2 [fLow,fHigh]
                ioSNR[n]   = xcorINJ_REC;		// inj*rec whitened energy [fLow,fHigh]

                *wfINJ*=1./rFactor;			// restore injected amplitude
             }
             delete [] pwfINJ;
             delete [] pwfREC;
          }
        }

        std::vector<int> sCuts = NET->getwc(k)->sCuts; 	// save sCuts
        // set sCuts=1 for the events which must be not copied with cps to pwc
        for(int i=0; i<(int)sCuts.size(); i++) if(i!=ID-1) NET->getwc(k)->sCuts[i]=1;  

        // after cpf, pwc contains only one event
        // ID can not be used to get the event, to get event use ID=1 (ex: for watplot)
        pwc->cpf(*(NET->getwc(k)));			// note: likelihood2G delete tdAmp 
        NET->getwc(k)->sCuts = sCuts; 			// restore sCuts

        if(cfg->dump) EVT->dopen(outDump.Data(),const_cast<char*>("a"),false);
        EVT->output2G(net_tree,NET,ID,k,ofactor);	// get reconstructed parameters
        if(cfg->dump) EVT->dclose();
        if(!cfg->cedDump) NET->getwc(k)->sCuts[ID-1]=1; // mark as processed                              

#ifdef SAVE_MRA_PLOT					// monster event display
        PlotMRA(pwc, nIFO, ID);
#endif
      }
    }

    jfile->cd();

    if(pwc)   delete pwc;
    if(iSNR)  delete [] iSNR;
    if(oSNR)  delete [] oSNR;
    if(ioSNR) delete [] ioSNR;
    if(EVT)   delete EVT;
  }

  return;
}

void ComputeResidualEnergy(wavearray<double>* wfINJ, wavearray<double>* wfREC, 
                           double& enINJ, double& enREC, double& xcorINJ_REC) {

  double bINJ = wfINJ->start();
  double eINJ = wfINJ->stop();
  double bREC = wfREC->start();
  double eREC = wfREC->stop();
  //cout.precision(14);
  //cout << "bINJ : " << bINJ << " eINJ : " << eINJ << endl;
  //cout << "bREC : " << bREC << " eREC : " << eREC << endl;

  double R=wfINJ->rate();

  int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
  int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);
  //cout << "oINJ : " << oINJ << " oREC : " << oREC << endl;

  double startXCOR = bINJ>bREC ? bINJ : bREC;
  double endXCOR   = eINJ<eREC ? eINJ : eREC;
  int sizeXCOR     = int((endXCOR-startXCOR)*R+0.5);
  //cout << "startXCOR : " << startXCOR << " endXCOR : " 
  //     << endXCOR << " sizeXCOR :" << sizeXCOR << endl;

  enINJ=0;
  enREC=0;
  xcorINJ_REC=0;
  if (sizeXCOR<=0) return;

  // the enINJ, enREC, xcorINJ_REC are computed in the INJ range
  for (int i=0;i<wfINJ->size();i++) enINJ+=wfINJ->data[i]*wfINJ->data[i];
  for (int i=0;i<sizeXCOR;i++) enREC+=wfREC->data[i+oREC]*wfREC->data[i+oREC];
  for (int i=0;i<sizeXCOR;i++) xcorINJ_REC+=wfINJ->data[i+oINJ]*wfREC->data[i+oREC];

  double erINJ_REC = enINJ+enREC-2*xcorINJ_REC;
  cout << "enINJ : " << enINJ << " enREC : " << enREC 
       << " xcorINJ_REC : " << xcorINJ_REC << " enINJREC : " << enINJ+enREC-2*xcorINJ_REC << endl;
  //cout << "erINJ_REC/enINJ : " << erINJ_REC/enINJ << endl;

  return;
}

void ApplyFrequencyCuts(wavearray<double>* wf, CWB::config* cfg) {

  double f_low  = cfg->fLow;
  double f_high = cfg->fHigh;
  cout << "f_low : " << f_low << " f_high : " << f_high << endl;
  wf->FFTW(1);
  double Fs=((double)wf->rate()/(double)wf->size())/2.;
  for (int j=0;j<wf->size();j+=2) {
    double f=j*Fs;
    if((f<f_low)||(f>f_high)) {wf->data[j]=0.;wf->data[j+1]=0.;}
  }
  wf->FFTW(-1);
  wf->resetFFTW();

  return;
}

void PlotMRA(netcluster* pwc, int nIFO, int ID) {

  watplot WTS(const_cast<char*>("wts"));
  WTS.canvas->cd();
  char fname[1024];
  sprintf(fname, "l_tfmap_scalogram_%d.png",ID);
  cout << "write " << fname << endl;
  WTS.plot(pwc, ID, nIFO, 'L', 0, const_cast<char*>("COLZ")); 
  WTS.canvas->Print(fname);
  WTS.clear();
}

void PlotWaveform(TString ifo, wavearray<double>* wfINJ, wavearray<double>* wfREC, 
                  CWB::config* cfg, bool fft, bool strain) {

  watplot PTS(const_cast<char*>("ptswrc"),200,20,800,500);

  //cout << "Print " << fname << endl;
  double tmin = wfINJ->start()<wfREC->start() ? wfINJ->start() : wfREC->start();
  wfINJ->start(wfINJ->start()-tmin);
  wfREC->start(wfREC->start()-tmin);
  if(fft) {
    PTS.plot(wfINJ, const_cast<char*>("ALP"), 1, 0, 0, true, cfg->fLow, cfg->fHigh);
  } else {
    PTS.plot(wfINJ, const_cast<char*>("ALP"), 1, 0, 0);
  }
  PTS.graph[0]->SetLineWidth(1);
  if(fft) {
    PTS.plot(wfREC, const_cast<char*>("SAME"), 2, 0, 0, true, cfg->fLow, cfg->fHigh);
  } else {
    PTS.plot(wfREC, const_cast<char*>("SAME"), 2, 0, 0);
  }
  PTS.graph[1]->SetLineWidth(2);
  wfINJ->start(wfINJ->start()+tmin);
  wfREC->start(wfREC->start()+tmin);

  char label[64]="";
  if(fft) sprintf(label,"%s","fft");
  else    sprintf(label,"%s","time");
  if(strain) sprintf(label,"%s_%s",label,"strain");
  else       sprintf(label,"%s_%s",label,"white"); 

  char fname[1024];
  sprintf(fname, "%s_wf_%s_inj_rec_gps_%d.root",ifo.Data(),label,int(tmin));
  PTS.canvas->Print(fname); 
  cout << "write : " << fname << endl;
  //PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
}

void ComputeResidualEnergyOptTF(wavearray<double>* wfINJ, wavearray<double>* wfREC, 
                           network* NET, netevent* EVT, int lag, int ID, int ifoID) {

  double enINJ=0;
  double enREC=0; 
  double enINJREC=0;

  // find TF optimal resolution level
  netcluster* pwc = NET->getwc(lag); 
  double optRate  = (pwc->cRate[ID-1])[0];
  double optLayer = pwc->rate/optRate;
  double optLevel = int(log10(optLayer)/log10(2));

  double bINJ = wfINJ->start();
  double eINJ = wfINJ->stop();
  double bREC = wfREC->start();
  double eREC = wfREC->stop();
  //cout.precision(14);
  //cout << "bINJ : " << bINJ << " eINJ : " << eINJ << endl;
  //cout << "bREC : " << bREC << " eREC : " << eREC << endl;

  double R=wfINJ->rate();

  int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
  int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);
  //cout << "oINJ : " << oINJ << " oREC : " << oREC << endl;

  double startXCOR = bINJ>bREC ? bINJ : bREC;
  double endXCOR   = eINJ<eREC ? eINJ : eREC;
  int sizeXCOR     = int((endXCOR-startXCOR)*R+0.5);
  //cout << "startXCOR : " << startXCOR << " endXCOR : " 
  //     << endXCOR << " sizeXCOR :" << sizeXCOR << endl;

  if (sizeXCOR<=0) return;

  watplot WTS(const_cast<char*>("wtswrc"));
  WTS.canvas->cd();

  detector* pd = NET->getifo(0);

  // length of temporary buffer for tf plots
  int xcor_length = sizeXCOR/wfINJ->rate();
  int wts_size = xcor_length<8 ? 16 : 16*int(1+xcor_length/8.);
  wts_size*=wfINJ->rate();

  char fname[1024];
  Meyer<double> S(1024,2);     // set wavelet for production
  double wts_start,wts_stop;

  // ------------------------------------------------------------------------
  // INJECTED WAVEFORM
  // ------------------------------------------------------------------------

  wavearray<double> xINJ(wts_size);
  xINJ.start(wfINJ->start()-EVT->gps[0]+double(oINJ+wts_size/2)/wfINJ->rate());
  xINJ.rate(wfINJ->rate());
  xINJ=0.;
  for(int m=0;m<sizeXCOR;m++)
    if(m<(int)xINJ.size()/2) xINJ.data[m+wts_size/2]=wfINJ->data[m+oINJ];
  WSeries<double> wINJ(xINJ,S);
  xINJ.resize(0);

  if(NET->getwdm(optLayer+1)) wINJ.Forward(wINJ,*NET->getwdm(optLayer+1));

#ifdef SAVE_TF_PLOT
  //scalogram maps
  sprintf(fname, "%s_wf_white_inj_tf.png", NET->ifoName[ifoID]);
  cout << "Print " << fname << endl;
  //PlotWSeries(&wINJ, fname);
  wts_start = wINJ.start()+(double)(wts_size/2)/wINJ.rate();
  wts_stop  = sizeXCOR<wts_size/2 ? wts_start+sizeXCOR/wINJ.rate() : wts_start+(wts_size/2)/wINJ.rate();
  WTS.plot(&wINJ, 2, wts_start, wts_stop,const_cast<char*>("COLZ"));
  WTS.hist2D->GetYaxis()->SetRangeUser(EVT->low[ifoID], EVT->high[ifoID]);
  WTS.canvas->Print(fname); 
  WTS.clear();
#endif

  // ------------------------------------------------------------------------
  // RECONSTRUCTED WAVEFORM
  // ------------------------------------------------------------------------

  wavearray<double> xREC(wts_size);
  xREC.start(wfREC->start()-EVT->gps[0]+double(oREC+wts_size/2)/wfREC->rate());
  xREC.rate(wfREC->rate());
  xREC=0.;
  for (int m=0;m<sizeXCOR;m++)
    if(m<(int)xREC.size()/2) xREC.data[m+wts_size/2]=wfREC->data[m+oREC];
  WSeries<double> wREC(xREC,S);
  xREC.resize(0);

  if(NET->getwdm(optLayer+1)) wREC.Forward(wREC,*NET->getwdm(optLayer+1));

#ifdef SAVE_TF_PLOT
  //scalogram maps
  sprintf(fname, "%s_wf_white_rec_tf.png", NET->ifoName[ifoID]);
  cout << "Print " << fname << endl;
  wts_start = wREC.start()+(double)(wts_size/2)/wREC.rate();
  wts_stop  = sizeXCOR<wts_size/2 ? wts_start+sizeXCOR/wREC.rate() : wts_start+(wts_size/2)/wREC.rate();
  WTS.plot(&wREC, 2, wts_start, wts_stop,const_cast<char*>("COLZ"));
  WTS.hist2D->GetYaxis()->SetRangeUser(EVT->low[ifoID], EVT->high[ifoID]);
  WTS.canvas->Print(fname); 
  WTS.clear();
#endif

  // ------------------------------------------------------------------------
  // INJECTED-RECONSTRUCTED WAVEFORM
  // ------------------------------------------------------------------------

  wavearray<double> xDIF(wts_size);
  xDIF.start(wfREC->start()-EVT->gps[0]+double(oREC+wts_size/2)/wfREC->rate());
  xDIF.rate(wfREC->rate());
  xDIF=0.;
  for (int m=0;m<sizeXCOR;m++)
    if(m<(int)xDIF.size()/2) xDIF.data[m+wts_size/2]=wfREC->data[m+oREC]-wfINJ->data[m+oINJ];
  WSeries<double> wDIF(xDIF,S);
  xDIF.resize(0);

#ifdef SAVE_TF_PLOT
  //scalogram maps
  sprintf(fname, "%s_wf_white_dif_tf.png", NET->ifoName[ifoID]);
  cout << "Print " << fname << endl;
  if(NET->getwdm(optLayer+1)) wDIF.Forward(wDIF,*NET->getwdm(optLayer+1));
  wts_start = wDIF.start()+(double)(wts_size/2)/wDIF.rate();
  wts_stop  = sizeXCOR<wts_size/2 ? wts_start+sizeXCOR/wDIF.rate() : wts_start+(wts_size/2)/wDIF.rate();
  WTS.plot(&wDIF, 2, wts_start, wts_stop,const_cast<char*>("COLZ"));
  WTS.hist2D->GetYaxis()->SetRangeUser(EVT->low[ifoID], EVT->high[ifoID]);
  WTS.canvas->Print(fname); 
  WTS.clear();
#endif

  // compute the energy 
  enINJ = ComputeEnergy(&wINJ);
  enREC = ComputeEnergy(&wREC);
  enINJREC = ComputeResidualEnergy(&wINJ,&wREC);

  cout << "TF : " << "enINJ : " << enINJ << " enREC : " << enREC 
       << " enINJREC : " << enINJREC << endl;

  return;
}

double ComputeEnergy(WSeries<double> *WS) {

  int layers = WS->maxLayer()+1;  // numbers of frequency bins (first & last bins have df/2)
  int slices = WS->sizeZero();    // number of time bins

  float df = WS->resolution();    // frequency bin resolution (hz)
  float dt = 1./(2*df);           // time bin resolution (sec)

  int rate = int(1./dt);

//  cout << "layers : " << layers << "\t slices : " << slices << "\t rate : " << rate
//       << "\t dt : " << dt << "\t df : " << df << endl;

  double EE=0;
  double pixel_frequency;
  double pixel_time;
  for(int i=1;i<=slices;i++) {
    pixel_time = i*dt;
    for(int j=1;j<=layers;j++) {
      if(j==1) pixel_frequency = df/2;
      if(j==layers) pixel_frequency = (layers-1)*df;
      if((j!=1)&&(j!=layers)) pixel_frequency = df/2 + (j-1)*df;
      double ER = pow(WS->getSample(i-1,j-1),2); 
      if(j==1) EE += ER/2.;
      if(j==layers) EE += ER/2.;
      if((j!=1)&&(j!=layers)) EE += ER;
      //cout << pixel_time << " : " << pixel_frequency << endl;  
    }
  }

  return EE;
}

double ComputeResidualEnergy(WSeries<double> *WS1, WSeries<double> *WS2) {

  int layers = WS1->maxLayer()+1;  // numbers of frequency bins (first & last bins have df/2)
  int slices = WS1->sizeZero();    // number of time bins

  if((layers != WS2->maxLayer()+1)||(slices != WS2->sizeZero())) {
    cout << "ComputeResidualEnergy - Error : WS1 & WS2 are not compatible !!!" << endl;
    exit(1);
  }

  float df = WS1->resolution();    // frequency bin resolution (hz)
  float dt = 1./(2*df);           // time bin resolution (sec)

  int rate = int(1./dt);

//  cout << "layers : " << layers << "\t slices : " << slices << "\t rate : " << rate
//       << "\t dt : " << dt << "\t df : " << df << endl;

  double EE=0;
  double pixel_frequency;
  double pixel_time;
  for(int i=1;i<=slices;i++) {
    pixel_time = i*dt;
    for(int j=1;j<=layers;j++) {
      if(j==1) pixel_frequency = df/2;
      if(j==layers) pixel_frequency = (layers-1)*df;
      if((j!=1)&&(j!=layers)) pixel_frequency = df/2 + (j-1)*df;
      double ER = pow(WS1->getSample(i-1,j-1)-WS2->getSample(i-1,j-1),2); 
      if(j==1) EE += ER/2.;
      if(j==layers) EE += ER/2.;
      if((j!=1)&&(j!=layers)) EE += ER;
      //cout << pixel_time << " : " << pixel_frequency << endl;  
    }
  }

  return EE;
}
