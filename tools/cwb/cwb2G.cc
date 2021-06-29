/*
# Copyright (C) 2019 Gabriele Vedovato, Sergey Klimenko
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


#include "cwb2G.hh"
#include "TKey.h"
#include "regression.hh"
#include "watplot.hh"
#include "sseries.hh"
#include "gwavearray.hh"
#include <iomanip>

// minimun skymap resolution used for subNetCut 
#define MIN_SKYRES_HEALPIX	4	
#define MIN_SKYRES_ANGLE	3	

// regression parameters
#define REGRESSION_FILTER_LENGTH	8
#define REGRESSION_MATRIX_FRACTION	0.95
#define REGRESSION_SOLVE_EIGEN_THR	0.
#define REGRESSION_SOLVE_EIGEN_NUM	10
#define REGRESSION_SOLVE_REGULATOR	'h'
#define REGRESSION_APPLY_THR		0.8

// select function parameter
#define SELECT_SUBRHO		5.0
#define SELECT_SUBNET		0.1

// WDM default parameters
#define WDM_BETAORDER		6	// beta function order for Meyer
#define WDM_PRECISION		10	// wavelet precision

#define EXIT(ERR) gSystem->Exit(ERR)   // better exit handling for ROOT stuff

ClassImp(cwb2G)

cwb2G::~cwb2G() {
//
// Destructor
//

  for(int i=0; i<NRES_MAX; i++) {
    if(pwdm[i]!=NULL) delete pwdm[i];
  }
}

void cwb2G::Init() {
//
// - Initialize the WDM<double> tranforms for each resolution level used in the analysis
// - Load the xTalk Catalog (network::setMRAcatalog) used for the Multi Resolution Analysis
//

  // export this (cwb2G) object (can be used by plugins)
  char cmd[256];
  sprintf(cmd,"gCWB2G = (void*)%p;",this); EXPORT(void*,gCWB2G,cmd);

  // check & set analysis stage
  if(TString(cfg.analysis)!="2G")
    {cout << "cwb2G::Init - Error : analysis must be 2G" << endl;EXIT(1);}
  if(istage!=CWB_STAGE_FULL && iname=="")
    {cout << "cwb2G::Init - Error : if initial stage!=CWB_STAGE_FULL "
          << "then input file name must be a job file " << endl;EXIT(1);}

  nRES = cfg.l_high-cfg.l_low+1;     // number of frequency resolution levels

  for(int i=0; i<NIFO_MAX; i++) hot[i]=NULL;
  for(int i=0; i<NRES_MAX; i++) pwdm[i]=NULL;

  // loading MRA catalog
  char MRAcatalog[1024];
  sprintf(MRAcatalog,"%s/%s",cfg.filter_dir,cfg.wdmXTalk);
  cout << "cwb2G::Init - Loading catalog of WDM cross-talk coefficients ... " << endl;
  cout << MRAcatalog << endl;
  TB.checkFile(MRAcatalog);
  NET.setMRAcatalog(MRAcatalog);

  int BetaOrder = WDM_BETAORDER;	// beta function order for Meyer
  int precision = WDM_PRECISION;	// wavelet precision
  if(NET.wdmMRA.tag!=0) {		// new catalog format : read BetaOrder,precision from catalog
    BetaOrder=NET.wdmMRA.BetaOrder; 
    precision=NET.wdmMRA.precision;
  }
  // print catalog infos 
  char info[256];
  sprintf(info,"-Tag:%d-BetaOrder:%d-precision:%d",NET.wdmMRA.tag,BetaOrder,precision);
  PrintAnalysisInfo(CWB_STAGE_INIT,"cwb2G::Init",info,true,false);

  //cout << "cwb2G::Init - Create WDM ... " << endl;
  for(int level=cfg.l_high; level>=cfg.l_low; level--) {
    int layers = level>0 ? 1<<level : 0;
    pwdm[cfg.l_high-level] = new WDM<double>(layers,layers,BetaOrder,precision);
    // check if filter lenght is less than cwb scratch length
    double wdmFLen = double(pwdm[cfg.l_high-level]->m_H)/rateANA;    // sec
    if(wdmFLen > cfg.segEdge+0.001) {
       cout << endl;
       cout << "cwb2G::Init : Error - filter length must be <= segEdge !!!" << endl;
       cout << "filter length : " << wdmFLen << " sec" << endl;
       cout << "cwb   scratch : " << cfg.segEdge << " sec" << endl;
       EXIT(1);
    } else {
       cout << "Filter length = " << wdmFLen << " (sec)" << endl;
    }
    // check if the length for time delay amplitudes is less than cwb scratch length
    // the factor 1.5 is used to avoid to use pixels on the border which could be distorted
    double rate  = rateANA>>level;
    if(cfg.segEdge<int(1.5*(cfg.TDSize/rate)+0.5)) {
       cout << endl;
       cout << "cwb2G::Init : Error - segEdge must be > " 
            << "1.5x the length for time delay amplitudes!!!" << endl;
       cout << "TD length : " << cfg.TDSize/rate << " sec" << endl;
       cout << "segEdge   : " << cfg.segEdge << " sec" << endl;
       cout << "Select segEdge > " << int(1.5*(cfg.TDSize/rate)+0.5) << endl << endl;
       EXIT(1);
    } 
    // add WDM to network
    NET.add(pwdm[cfg.l_high-level]); // network vector must be filled starting from max resolution level
  }

  // check if analysis layers are contained in the MRAcatalog
  // level : is the decomposition level
  // layes : are the number of layers along the frequency axis rateANA/(rateANA>>level)
  int check_layers=0;
  for(int level=cfg.l_high; level>=cfg.l_low; level--) {
    int layers = level>0 ? 1<<level : 0;
    for(int j=0;j<NET.wdmMRA.nRes;j++) if(layers==NET.wdmMRA.layers[j]) check_layers++;
  }

  if(check_layers!=nRES) {
    cout << "cwb2G::Init - Error : analysis layers do not match the MRA catalog" << endl;
    cout << endl << "analysis layers : " << endl;
    for(int level=cfg.l_high; level>=cfg.l_low; level--) {
      int layers = level>0 ? 1<<level : 0;
      cout << "level : " << level << " layers : " << layers << endl;
    }
    cout << endl << "MRA catalog layers : " << endl;
    for(int i=0;i<NET.wdmMRA.nRes;i++) 
       cout << "layers : " << NET.wdmMRA.layers[i] << endl;
    EXIT(1);
  } 
  else {
    cout << endl;
    for(int level=cfg.l_high; level>=cfg.l_low; level--) {
      int layers = level>0 ? 1<<level : 0;
      int rate  = rateANA>>level;
      cout << "level : " << level << "\t rate(hz) : " << rate << "\t layers : " << layers
           << "\t df(hz) : " << rateANA/2./double(1<<level)
           << "\t dt(ms) : " << 1000./rate << endl;
    }
    cout << endl;
  }

  // Check if lagStep compatible with WDM parity 
  // This condition is necessary to avoid mixing between odd 
  // and even pixels when circular buffer is used for lag shift
  // The MRAcatalog distinguish odd and even pixels
  int rate_min  = rateANA>>cfg.l_high;
  double dt_max = 1./rate_min;
  if(fmod(rate_min,1.)) {
    cout << "cwb2G::Init - Error : rate min=" << rate_min << "(Hz) is not integer" << endl << endl;
    EXIT(1);  
  }
  if(int(cfg.lagStep*rate_min+0.001)&1) {
    cout << "cwb2G::Init - Error : lagStep=" << cfg.lagStep << "(sec)"
         << " is not a multple of 2*max_time_resolution=" << 2*dt_max << "(sec)" << endl << endl;
    EXIT(1);  
  }
  if(int(cfg.segEdge*rate_min+0.001)&1) {
    cout << "cwb2G::Init - Error : segEdge=" << cfg.segEdge << "(sec)"
         << " is not a multple of 2*max_time_resolution=" << 2*dt_max << "(sec)" << endl << endl;
    EXIT(1);  
  }
  if(int(cfg.segMLS*rate_min+0.001)&1) {
    cout << "cwb2G::Init - Error : segMLS=" << cfg.segMLS << "(sec)"
         << " is not a multple of 2*max_time_resolution=" << 2*dt_max << "(sec)" << endl << endl;
    EXIT(1);  
  }

  // time-delay filter rate
  if(cfg.fResample>0) {                                 // RESAMPLING
    TDRate = (cfg.fResample>>cfg.levelR)*cfg.upTDF;	
  } else {
    TDRate = (cfg.inRate>>cfg.levelR)*cfg.upTDF;	
  }

  return;
}

double cwb2G::ReadData(double mdcShift, int ifactor) {
//
// Read Noise & MDC data from frame file or "On The Fly" from plugin
//
// Loop over detectors
// - Read noise from frames or "On The Fly" from plugin
//   - Resampling data
// - Read injections from frames or "On The Fly" from plugin (config::simulation>0)
//   - Resampling data
// - if(simulation==2) MDC are rescaled (detector::setsnr) with a fixed 
//                     network SNR according to the config::factors
// - Store noise & MDC to job file
//

  // if 2G analysis istage>=CWB_STAGE_STRAIN then skip read data
  if(istage>=CWB_STAGE_STRAIN) return 0.; 

  PrintStageInfo(CWB_STAGE_STRAIN,"cwb2G::ReadData");

  // data are loaded from root file
  if((istage!=CWB_STAGE_INIT)&&(iname!="")) return cwb::ReadData(iname);

  Meyer<double> B(1024);           // set wavelet for resampling
  wavearray<double> x,z,w;         // temporary time series
  std::vector<wavearray<double> > y; // temporary time series for snr mode
  y.resize(nIFO);
  WSeries<double> wM;              // mdc WSeries
  wavearray<double>* px;

  int layers_high = 1<<cfg.l_high;
  WDM<double> WDMwhite(layers_high,layers_high,6,10);         // set whitening WDM

  // check if wdm filter lenght is less than cwb scratch
  double wdmFLen = double(WDMwhite.m_H)/rateANA; // sec
  if(wdmFLen > cfg.segEdge+0.001) {
     cout << endl;
     cout << "cwb2G::ReadData : Error - filter scratch must be <= cwb scratch!!!" << endl;
     cout << "filter length : " << wdmFLen << " sec" << endl;
     cout << "cwb   scratch : " << cfg.segEdge << " sec" << endl;
     EXIT(1);
  } else {
     cout << "WDM filter length for regression = " << wdmFLen << " (sec)" << endl;
  }

  int    xsize=0.;
  double xrate=0.;
  double xstart=0.;

  // data are loaded from frame files
  jfile = new TFile(jname,"UPDATE");
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb2G::ReadData - Error opening root file : " << jname << endl;EXIT(1);}
  TDirectory* cdstrain = (TDirectory*)jfile->Get("strain");
  if(cdstrain==NULL) cdstrain = jfile->mkdir("strain");

  for(int i=0; i<nIFO; i++) {

    if(cfg.simulation==4) {       	  		  // sim4 -> read ifo strain from job file
      px = (wavearray<double>*)jfile->Get(TString("strain/")+ifo[i]);
    }
    if((cfg.simulation==4)&&(px!=NULL)) {       	  // sim4 -> use strain from job file
      x = *px; delete px;
    } else {
      if(cfg.dataPlugin) {  	                          // data are provided by the user plugin
        x.rate(cfg.inRate); x.start(FRF[i].start); x.resize(int(x.rate()*(FRF[i].stop-FRF[i].start)));
      } else {  			                  // data are read from frame
        fr[i].readFrames(FRF[i],cfg.channelNamesRaw[i],x);
      }
      if(bplugin) CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)&x,ifo[i],CWB_PLUGIN_STRAIN);
      if(TMath::IsNaN(x.mean())) 
        {cout << "cwb2G::ReadData - Error : found NaN in strain data !!!" <<  endl;EXIT(1);}

      if(x.rate()!=cfg.inRate)
        {cout << "cwb2G::ReadData - input rate from frame " << x.rate()
              << " do not match the one defined in config : " << cfg.inRate << endl;EXIT(1);}
  
      x.start(x.start()+cfg.dataShift[i]);                // dataShift
      x.start(x.start()-cfg.segLen*(segID[i]-segID[0]));  // SLAG
      if(singleDetector) TB.resampleToPowerOfTwo(x);
      if(cfg.dcCal[i]>0.) x*=cfg.dcCal[i];                // DC correction
      if(cfg.fResample>0) x.Resample(cfg.fResample);   	  // RESAMPLING
      x.Resample(x.rate()/(1<<cfg.levelR));		  // resampling
      x*=sqrt(1<<cfg.levelR);				  // rescaling

      if(cfg.simulation==2) w = x;                        // snr mode - save strain

      // save ifo data to temporary job file
      cdstrain->cd();gwavearray<double> gx(x);gx.Write(ifo[i],TObject::kOverwrite);
    }

    if(i==0) {xrate=x.rate();xstart=x.start();xsize=x.size();}

    fprintf(stdout,"start=%f duration=%f rate=%f\n", x.start(),x.size()/x.rate(),x.rate());
    if(i>0 && xstart != x.start()) {
      cout << "cwb2G::ReadData - Error : ifo noise data not synchronized" << endl;
      cout << ifo[i] << " " << x.start() << " != " << ifo[0] << " " << xstart << endl;
      EXIT(1);
    }
    if(i>0 && xrate != x.rate()) {
      cout << "cwb2G::ReadData - Error : ifo noise data have different rates" << endl;
      cout << ifo[i] << " " << x.rate() << " != " << ifo[0] << " " << xrate << endl;
      EXIT(1);
    }

    if(cfg.simulation) {
      TDirectory* cdmdc = (TDirectory*)jfile->Get("mdc");
      if(cdmdc==NULL) cdmdc = jfile->mkdir("mdc");

      if(cfg.mdcPlugin) {	// mdc are provided by the user plugin
        x.rate(cfg.inRate); x.start(FRF[i+nIFO].start); 
        x.resize(int(x.rate()*(FRF[i+nIFO].stop-FRF[i+nIFO].start)));
      } else {			// mdc are read from frame
        fr[nIFO].readFrames(FRF[i+nIFO],cfg.channelNamesMDC[i],x);
        if(x.rate()!=cfg.inRate)
          {cout << "cwb2G::ReadData - input rate from frame " << x.rate()
                << " do not match the one defined in config : " << cfg.inRate << endl;EXIT(1);}
      }
      if(bplugin) CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)&x,ifo[i],CWB_PLUGIN_MDC);
      if(TMath::IsNaN(x.mean())) 
        {cout << "cwb2G::ReadData - Error : found NaN in MDC data !!!" <<  endl;EXIT(1);}

      x.start(x.start()+mdcShift);
      if(cfg.fResample>0) x.Resample(cfg.fResample);   	// RESAMPLING
      y[i] = x;                   			// save inj for snr mode
      x.Resample(x.rate()/(1<<cfg.levelR));		// resampling
      x*=sqrt(1<<cfg.levelR);			        // rescaling

      if(cfg.simulation==2) { 				// snr mode

        // calculate noise rms 
        pTF[i] = pD[i]->getTFmap();
        pTF[i]->Forward(w,WDMwhite);            
        pTF[i]->setlow(cfg.fLow);
        pTF[i]->sethigh(cfg.fHigh);
        pD[i]->white(cfg.whiteWindow,0,cfg.segEdge,cfg.whiteStride);  	

        // compute snr
        wM.Forward(x,WDMwhite);            
        wM.setlow(cfg.fLow);
        wM.sethigh(cfg.fHigh);
        // zero f<fLow to avoid whitening issues when psd noise is not well defined for f<fLow
        int layers  = wM.maxLayer();             	
        for(int j=0;j<layers;j++) if(wM.frequency(j)<cfg.fLow) {
          double layer = j+0.01;				// -epsilon select 0 layer for 90 phase
          wM.getLayer(z, layer);z=0;wM.putLayer(z, layer);	//  0 phase
          wM.getLayer(z,-layer);z=0;wM.putLayer(z,-layer);	// 90 phase
        }
 
        // compute mdc snr 
        pD[i]->setsim(wM,NET.getmdcTime(),cfg.iwindow/2.,cfg.segEdge,false);   
      } 
      else {                                // save mdc data to temporary job file
        cdmdc->cd();gwavearray<double> gx(x);gx.Write(ifo[i],TObject::kOverwrite);
      }

      fprintf(stdout,"start=%f duration=%f rate=%f\n", x.start(),x.size()/x.rate(),x.rate());
      if(xstart != x.start()) {
        cout << "cwb2G::ReadData - Error : mdc/noise data with different start time" << endl; 
        printf("start time : noise = %10.6f - mdc = %10.6f\n",xstart,x.start());
        EXIT(1);
      }
      if(xrate != x.rate()) {
        cout << "cwb2G::ReadData - Error : mdc/noise data with different rate" << endl; 
        printf("rate : noise = %10.6f - mdc = %10.6f\n",xrate,x.rate());
        EXIT(1);
      }
      if(xsize != x.size()) {
        cout << "cwb2G::ReadData - Error : mdc/noise data with different buffer size" << endl; 
        printf("buffer size : noise = %d - mdc = %lu\n",xsize,x.size());
        EXIT(1);
      }
    }
    if(singleDetector) break;
  }

  // for snr mode the factors parameters set to be the mdc snr
  if(cfg.simulation==2) {
    TDirectory* cdmdc = (TDirectory*)jfile->Get("mdc");
    if(cdmdc==NULL) cdmdc = jfile->mkdir("mdc");

    // compute rescale factor -> snr=1
    std::vector<double>  mdcFactor;
    for (int k=0;k<(int)pD[0]->ISNR.size();k++) {
      double snr=0;
      if(singleDetector) {
        for(int i=0; i<nIFO; i++) snr+=pD[0]->ISNR.data[k];
      } else {
        for(int i=0; i<nIFO; i++) snr+=pD[i]->ISNR.data[k];
      }
      snr=sqrt(snr);
      if(snr>0) mdcFactor.push_back(1./snr); else mdcFactor.push_back(0.);
    }

    size_t K = mdcFactor.size();
    for(int k=0; k<(int)K; k++) {
      if(mdcFactor[k]) cout << k << " mdcFactor : " << mdcFactor[k] << endl;
    }

    // rescale mdc snr to 1
    for(int i=0; i<nIFO; i++) {
      pD[i]->setsnr(y[i],NET.getmdcTime(),&mdcFactor,cfg.iwindow/2.,cfg.segEdge);
      if(bplugin) CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)&y[i],ifo[i],CWB_PLUGIN_RMDC);

      wM.Forward(y[i],B,cfg.levelR);
      wM.getLayer(y[i],0);
      cdmdc->cd();gwavearray<double> gyi(y[i]);gyi.Write(ifo[i],TObject::kOverwrite);

      // rescaled saved injected waveforms
      for(int k=0; k<(int)pD[i]->IWFP.size(); k++) {
        wavearray<double>* pwf = pD[i]->IWFP[k];
        int iwfid = pD[i]->IWFID[k];
        for(int j=0;j<(int)pwf->size();j++) pwf->data[j]*=mdcFactor[iwfid];
      }
      for(int k=0; k<(int)K; k++) {
        pD[i]->HRSS[k]*=mdcFactor[k];
        pD[i]->ISNR[k]*=mdcFactor[k];
      }

      if(singleDetector) break;
    }
    // rescale amplitudes stored in the mdcList
    for(int k=0; k<(int)K; k++) {
      int ilog[5] = {1,3,12,13,14};
      for(int l=0;l<5;l++) {
        double mfactor = l<2 ? mdcFactor[k] : mdcFactor[k]*mdcFactor[k];
        TString slog = TB.GetMDCLog(NET.mdcList[k], ilog[l]);
        NET.mdcList[k]=TB.SetMDCLog(NET.mdcList[k], ilog[l], mfactor*slog.Atof());
      }
    }
  }

  // mdc data are copied to detectors HoT to be accessible by the user plugin
  for(int i=0; i<nIFO; i++) pD[i]->HoT = y[i];	
  if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_OREADDATA);
  // clean mdc from detectors HoT
  for(int i=0; i<nIFO; i++) {pD[i]->HoT.resize(0);y[i].resize(0);}	

  jfile->Close();

  x.resize(0);
  z.resize(0);
  w.resize(0);
  wM.resize(0);

  return x.rate();
}

void
cwb2G::DataConditioning(int ifactor) {
//
// Apply regression to remove lines & whiten data
//
// Loop over detectors
// - read ifo strain from job file
// - read MDC data from temporary job file (config::simulation>0)
// - if(config::simulation==1) MDC are rescaled according to the config::factors
// - Add MDC to noise
// - Apply regression to remove lines 
// - Use detector::white to estimate noise (detector::nRMS) 
// - Use the estimated noise to whiten data (WSeries<double>::white)
// - Store injected waveforms (SaveWaveforms)
// - Store whitened data (detector::HoT) to job file (jfile)
// - Store estimated noise to job file (detector::nRMS) 
//

  // data are loaded from root file
  if((istage>CWB_STAGE_STRAIN)&&(istage<CWB_STAGE_SUPERCLUSTER)) 
    return DataConditioning(iname, ifactor);
  // if 2G analysis istage==CWB_STAGE_SUPERCLUSTER then skip data conditioning
  if(istage==CWB_STAGE_SUPERCLUSTER) return; 

  PrintStageInfo(CWB_STAGE_CSTRAIN,"cwb2G::DataConditioning");

  char info[256];
  WSeries<double> wM;              		
  wavearray<double> xM;
  wavearray<double> x;
  wavearray<double>* px;
  TDirectory* cdrms=NULL;
  TDirectory* cdcstrain=NULL;
  double factor=cfg.factors[ifactor];

  int layers_high = 1<<cfg.l_high;
  WDM<double> WDMwhite(layers_high,layers_high,6,10);         // set whitening WDM
  int layers = rateANA/8;   
  WDM<double> WDMlpr(layers,layers,6,10);           // set LPE filter WDM

  // check if whitening WDM filter lenght is less than cwb scratch
  double wdmFLen = double(WDMwhite.m_H)/rateANA;    // sec
  if(wdmFLen > cfg.segEdge+0.001) {
     cout << endl;
     cout << "cwb2G::DataConditioning : Error - filter scratch must be <= cwb scratch!!!" << endl;
     cout << "filter length : " << wdmFLen << " sec" << endl;
     cout << "cwb   scratch : " << cfg.segEdge << " sec" << endl;
     EXIT(1);
  } else {
     cout << "WDM filter max length = " << wdmFLen << " (sec)" << endl;
  }

  // open input job file
  TFile* ifile = ((istage==CWB_STAGE_STRAIN)&&(iname!="")) ? new TFile(iname) : NULL;
  // open temporary job file
  jfile = new TFile(jname, "UPDATE");
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb2G::DataConditioning - Error : file " << jname << " not found" <<  endl;EXIT(1);}

  // create cstrain,rms root dirs
  if(jobfOptions&CWB_JOBF_SAVE_CSTRAIN) {
    TDirectory* dcstrain = (TDirectory*)jfile->Get("cstrain");
    if(dcstrain==NULL) dcstrain=jfile->mkdir("cstrain"); 
    char cdcstrain_name[32];sprintf(cdcstrain_name,"cstrain-f%d",ifactor);
    cdcstrain=dcstrain->mkdir(cdcstrain_name);

    TDirectory* drms = (TDirectory*)jfile->Get("rms");
    if(drms==NULL) drms=jfile->mkdir("rms"); 
    char cdrms_name[32];sprintf(cdrms_name,"rms-f%d",ifactor);
    cdrms=drms->mkdir(cdrms_name);
  }

  for(int i=0; i<nIFO; i++) {
    // read ifo strain from temporary job file
    if(ifile!=NULL) px = (wavearray<double>*)ifile->Get(TString("strain/")+ifo[i]);
    else            px = (wavearray<double>*)jfile->Get(TString("strain/")+ifo[i]);
    if(TMath::IsNaN(px->mean())) 
      {cout << "cwb2G::DataConditioning - Error : found NaN in strain data !!!" <<  endl;EXIT(1);}
    hot[i] = pD[i]->getHoT(); 
    *hot[i] = *px; delete px;

    if(cfg.simulation) {
      // read mdc data from temporary job file
      if(ifile!=NULL) px = (wavearray<double>*)ifile->Get(TString("mdc/")+ifo[i]);
      else            px = (wavearray<double>*)jfile->Get(TString("mdc/")+ifo[i]);
      if(TMath::IsNaN(px->mean())) 
        {cout << "cwb2G::DataConditioning - Error : found NaN in MDC data !!!" <<  endl;EXIT(1);}

      if(cfg.simulation==3) {                   // time shift   : factor is the shift time
        int nshift = int(factor*px->rate());   	// number of shifted samples
        wavearray<double> pX = *px;(*px)=0;	// temporary array
        int jstart = nshift<0 ? -nshift : 0;
        int jstop  = nshift<0 ? px->size() : px->size()-nshift;
        for(int j=jstart;j<jstop;j++) px->data[j+nshift] = pX.data[j];

        double tshift=nshift/px->rate();	// time shift (sec)
        // take into account of the previous applied time shift  	
        tshift = ifactor==0 ? tshift : tshift-int(cfg.factors[ifactor-1]*px->rate())/px->rate();	

        // tshift saved injected waveforms 
        for(int k=0; k<(int)pD[i]->IWFP.size(); k++) {
          wavearray<double>* pwf = pD[i]->IWFP[k];
          pwf->start(pwf->start()+tshift);	
        }
        // tshift saved central times
        for(int k=0; k<(int)pD[i]->TIME.size(); k++) pD[i]->TIME[k]+=tshift;

        // shift times stored in the NET.mdcList & NET.mdcTime
        if(i==0) {
          vector<TString> ifos(nIFO);
          for(int n=0;n<nIFO;n++) ifos[n]=ifo[n];
          TB.shiftBurstMDCLog(NET.mdcList, ifos, tshift);       			
          for(int k=0;k<(int)NET.mdcTime.size();k++) NET.mdcTime[k]+=tshift;
        } 
        xM = *px;				// copy MDC to temporary wavearray
      } else if(cfg.simulation==4) {
        xM = *px;                              // copy MDC to temporary wavearray
      } else {
        xM = *px;				// copy MDC to temporary wavearray
        (*px)*=factor;
      }
      hot[i]->add(*px);
      delete px;
    }

    if(bplugin) 
      CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)hot[i],ifo[i],CWB_PLUGIN_IDATA_CONDITIONING);

    //  2G Data Conditioning & Whitening

    if(!cfg.dcPlugin) {    	// built in data conditioning 
      pTF[i] = pD[i]->getTFmap();

      // regression
      pTF[i]->Forward(*hot[i],WDMlpr);
      regression rr(*pTF[i],const_cast<char*>("target"),1.,cfg.fHigh);
      rr.add(*hot[i],const_cast<char*>("target"));
      rr.setFilter(REGRESSION_FILTER_LENGTH);
      rr.setMatrix(NET.Edge,REGRESSION_MATRIX_FRACTION);
      rr.solve(REGRESSION_SOLVE_EIGEN_THR,REGRESSION_SOLVE_EIGEN_NUM,REGRESSION_SOLVE_REGULATOR);
      rr.apply(REGRESSION_APPLY_THR);
      *hot[i] = rr.getClean();

      // whitening
      pTF[i]->Forward(*hot[i],WDMwhite);            
      pTF[i]->setlow(cfg.fLow);
      pTF[i]->sethigh(cfg.fHigh);
      pD[i]->white(cfg.whiteWindow,0,cfg.segEdge,cfg.whiteStride);  	// calculate noise rms 
      pD[i]->nRMS.bandpass(16.,0.,1); 	                                // high pass filtering at 16Hz
      pTF[i]->white(pD[i]->nRMS,1);               			// whiten  0 phase WSeries
      pTF[i]->white(pD[i]->nRMS,-1);              			// whiten 90 phase WSeries
      if(cfg.simulation) {						// estimated whitened MDC parms
        wM.Forward(xM,WDMwhite);            
        wM.setlow(cfg.fLow);
        wM.sethigh(cfg.fHigh);
        // zero f<fLow to avoid whitening issues when psd noise is not well defined for f<fLow
        int layers  = wM.maxLayer();      
        for(int j=0;j<layers;j++) if(wM.frequency(j)<cfg.fLow) {
          double layer = j+0.01;				// -epsilon select 0 layer for 90 phase
          wM.getLayer(x, layer);x=0;wM.putLayer(x, layer);	//  0 phase
          wM.getLayer(x,-layer);x=0;wM.putLayer(x,-layer);	// 90 phase
        }
	// compute mdc snr and save whiten waveforms
        pD[i]->setsim(wM,NET.getmdcTime(),cfg.iwindow/2.,cfg.segEdge,true);   
      }
      WSeries<double> wtmp = *pTF[i];
      pTF[i]->Inverse();
      wtmp.Inverse(-2);
      *hot[i] = *pTF[i];
      *hot[i] += wtmp;
      *hot[i] *= 0.5;
      // add infos to history
      sprintf(info,"-IFO:%d-RMS:%g",i,hot[i]->rms());
      PrintAnalysisInfo(CWB_STAGE_CSTRAIN,"cwb2G::DataConditioning",info,false);
    } else { 			// data conditioning is provided by the user plugin
      char cmd[128];
      // export to CINT variables 
      sprintf(cmd,"gRATEANA = %lu;",rateANA); EXPORT(size_t,gRATEANA,cmd);
      sprintf(cmd,"gMDC = %p;",&xM); EXPORT(void*,gMDC,cmd);
      if(bplugin) 
        CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)hot[i],ifo[i],CWB_PLUGIN_DATA_CONDITIONING);
    }

    if(bplugin) 
       CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)hot[i],ifo[i],CWB_PLUGIN_ODATA_CONDITIONING);

    if((jstage!=CWB_STAGE_FULL)&&(jstage<CWB_STAGE_LIKELIHOOD)) { 
       SaveWaveforms(jfile, pD[i], ifactor); 	// write inj waveforms to output job file
    }

    if(jobfOptions&CWB_JOBF_SAVE_CSTRAIN) {
       cdcstrain->cd();gwavearray<double> ghi(*hot[i]);ghi.Write(ifo[i]); 
       cdrms->cd();pD[i]->nRMS.Write(ifo[i]);
    }

    cout<<"After "<<ifo[i]<<" data conditioning"<<endl;
    gSystem->Exec("/bin/date"); GetProcInfo();

    if(singleDetector) {
      *pD[1]=*pD[0];
      // copy detector data not implemented in the copy operator
      pD[1]->HRSS  = pD[0]->HRSS;
      pD[1]->ISNR  = pD[0]->ISNR;
      pD[1]->FREQ  = pD[0]->FREQ;
      pD[1]->BAND  = pD[0]->BAND;
      pD[1]->TIME  = pD[0]->TIME;
      pD[1]->TDUR  = pD[0]->TDUR;
      pD[1]->IWFID = pD[0]->IWFID;
      pD[1]->IWFP  = pD[0]->IWFP;
      pD[1]->RWFID = pD[0]->RWFID;
      pD[1]->RWFP  = pD[0]->RWFP;
      break;
    }
  }
  if(ifile!=NULL) ifile->Close();
  jfile->Close();

  xM.resize(0); wM.resize(0); x.resize(0);

  // strains and mdc data are removed if not set in the jobfOptions (only for the last factor)
  int ioffset = (cfg.simulation==4) ? int(cfg.factors[0]) : 0; // ifactor offset
  if(ifactor==ioffset+cfg.nfactor-1) {  // the last factor
    vector<TString> delObjList;
    if(!(jobfOptions&CWB_JOBF_SAVE_STRAIN)) delObjList.push_back("strain");
    if(cfg.simulation && !(jobfOptions&CWB_JOBF_SAVE_MDC)) delObjList.push_back("mdc");
    FileGarbageCollector(jname,"",delObjList);
  }

  return;
}

void
cwb2G::DataConditioning(TString fName, int ifactor) {
//
// Read Conditioned Data from Job File (from the previous processed stage)
//
// Loop over detectors
// - Read conditioned data
//

  PrintStageInfo(CWB_STAGE_CSTRAIN,"cwb2G::DataConditioning from file");

  char cdrms_name[32];sprintf(cdrms_name,"rms-f%d",ifactor);
  char cdcstrain_name[32];sprintf(cdcstrain_name,"cstrain-f%d",ifactor);

  // data are loaded from input root file
  TFile* ifile = new TFile(fName);
  if(ifile==NULL)
    {cout << "cwb2G::DataConditioning - Error opening root file : " << fName << endl;EXIT(1);}
  // open temporary job file
  jfile = new TFile(jname,"UPDATE");
  if(jfile==NULL||!jfile->IsOpen())
    {cout << "cwb::DataConditioning - Error opening root file : " << jname << endl;EXIT(1);}

  // open job name
  TDirectory* jcdrms = NULL;
  TDirectory* jcdcstrain = NULL;
  if(jobfOptions&CWB_JOBF_SAVE_CSTRAIN) {
    TDirectory* dcstrain = (TDirectory*)jfile->Get("cstrain");
    if(dcstrain==NULL) dcstrain=jfile->mkdir("cstrain"); 
    jcdcstrain=dcstrain->mkdir(cdcstrain_name);

    TDirectory* drms = (TDirectory*)jfile->Get("rms");
    if(drms==NULL) drms=jfile->mkdir("rms"); 
    jcdrms=drms->mkdir(cdrms_name);
  }

  WSeries<double>* pws=NULL;
  for(int i=0; i<nIFO; i++) {
    pTF[i] = pD[i]->getTFmap();
    pTF[i]->setlow(cfg.fLow);
    pTF[i]->sethigh(cfg.fHigh);
    pTF[i]->w_mode=1;

    if((istage>=CWB_STAGE_CSTRAIN)&&(jstage<=CWB_STAGE_LIKELIHOOD)) { 
      LoadWaveforms(ifile, pD[i], ifactor);	// read inj waveforms from input job file
      SaveWaveforms(jfile, pD[i], ifactor);	// write inj waveforms to output job file
    }

    // restore ctrain
    ifile->cd();
    pws = (WSeries<double>*)ifile->Get(TString("cstrain/")+cdcstrain_name+"/"+pD[i]->Name);
    if(pws==NULL)
      {cout << "cwb2G::DataConditioning - Error : cstrain not present, job terminated!!!" << endl;EXIT(1);}
    if(jobfOptions&CWB_JOBF_SAVE_CSTRAIN) 
      {jfile->cd();jcdcstrain->cd();pws->Write(ifo[i]);}
    hot[i] = pD[i]->getHoT(); 
    *hot[i]=*pws;
    *pTF[i]=*hot[i];
    delete pws;

    // restore strain rms
    ifile->cd();
    pws = (WSeries<double>*)ifile->Get(TString("rms/")+cdrms_name+"/"+pD[i]->Name);
    if(pws==NULL)
      {cout << "cwb2G::DataConditioning - Error : strain rms not present, job terminated!!!" << endl;EXIT(1);}
    if(jobfOptions&CWB_JOBF_SAVE_CSTRAIN) 
      {jfile->cd();jcdcstrain->cd();jcdrms->cd();pws->Write(ifo[i]);}
    pD[i]->nRMS=*pws;; 
    delete pws;

    if(singleDetector) {
      *pD[1]=*pD[0];
      // copy detector data not implemented in the copy operator
      pD[1]->HRSS  = pD[0]->HRSS;
      pD[1]->ISNR  = pD[0]->ISNR;
      pD[1]->FREQ  = pD[0]->FREQ;
      pD[1]->BAND  = pD[0]->BAND;
      pD[1]->TIME  = pD[0]->TIME;
      pD[1]->TDUR  = pD[0]->TDUR;
      pD[1]->IWFID = pD[0]->IWFID;
      pD[1]->IWFP  = pD[0]->IWFP;
      pD[1]->RWFID = pD[0]->RWFID;
      pD[1]->RWFP  = pD[0]->RWFP;
      break;
    }
  }

  if(jfile!=NULL) jfile->Close();
  ifile->Close();
  return;
}

void
cwb2G::Coherence(int ifactor) {
//
// Select the significant pixels
//
// Loop over resolution levels (nRES)
// - Loop over detectors (cwb::nIFO)
//   - Compute the maximum energy of TF pixels (WSeries<double>::maxEnergy)
// - Set pixel energy selection threshold (network::THRESHOLD)
// - Loop over time lags (network::nLag)
//   - Select the significant pixels (network::getNetworkPixels)
//   - Single resolution clustering (network::cluster)
//   - Store selected pixels to job file (netcluster::write)
//

  // if 2G analysis istage>=CWB_STAGE_COHERENCE then skip coherence
  if(istage>=CWB_STAGE_COHERENCE) return; 

  if(bplugin) {
    char sfactor[8];sprintf(sfactor,"%d",ifactor);
    CWB_Plugin(jfile,&cfg,&NET,NULL,sfactor,CWB_PLUGIN_BCOHERENCE);
  }

  PrintStageInfo(CWB_STAGE_COHERENCE,"cwb2G::Coherence");

  char info[256];
  double Eo;
  netcluster* pwc;
  netcluster wc;

  // used to build sparse map
  SSeries<double> SS;
  std::vector<SSeries<double> > vSS;
  for(int n=0; n<nIFO; n++) vSS.push_back(SS);
  WSeries<double> WS[NIFO_MAX];

  // open temporary job file
  jfile = new TFile(jname,"UPDATE");
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb2G::Coherence - Error : file " << jname << " not found" <<  endl;EXIT(1);}

  // create sparse directories in the job root file to store sparse TF maps
  TDirectory* sdir = (TDirectory*)jfile->Get("csparse");
  if(sdir==NULL) sdir = jfile->mkdir("csparse");

  if(bplugin) {
    char sfactor[8];sprintf(sfactor,"%d",ifactor);
    CWB_Plugin(jfile,&cfg,&NET,NULL,sfactor,CWB_PLUGIN_ICOHERENCE);
  }

  if(!cfg.cohPlugin) {    				// built in coherence stage 
    TStopwatch watchCoh;					// coherence benchmark
    int upN = rateANA/1024; if(upN<1) upN=1;                // calculate upsample factor

    for(int i=0; i<nRES; i++) {  				// loop over TF resolutions
      watchCoh.Start();
      // print level infos
      int level=cfg.l_high-i;
      int layers = level>0 ? 1<<level : 0;
      int rate  = rateANA>>level;
      cout << "level : " << level << "\t rate(hz) : " << rate 
           << "\t layers : " << layers << "\t df(hz) : " << rateANA/2./double(1<<level) 
           << "\t dt(ms) : " << 1000./rate << endl;
  
      // produce TF maps with max over the sky energy
      double alp=0;
      for(int n=0; n<nIFO; n++) {
        alp+=NET.getifo(n)->getTFmap()->maxEnergy(*hot[n],*pwdm[i],mTau,upN,NET.pattern);
        // restore the frequency boundaries changed by the maxEnergy call
        NET.getifo(n)->getTFmap()->setlow(cfg.fLow);
        NET.getifo(n)->getTFmap()->sethigh(cfg.fHigh);
        if(singleDetector) {
          *(NET.getifo(1)->getTFmap()) = *(NET.getifo(0)->getTFmap());
          break;
        }
      }
  
      if(bplugin) {					// call user plugin
        char cmd[128];
        // export resolution level to CINT  
        sprintf(cmd,"gILEVEL = %lu;",level); EXPORT(size_t,gILEVEL,cmd);
        CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_ECOHERENCE);
      }

      // threshold on pixel energy
      alp /= nIFO;
      if(NET.pattern!=0) { 
        Eo = NET.THRESHOLD(cfg.bpp,alp);          
      } else {
        Eo = NET.THRESHOLD(cfg.bpp); 
      }
      cout.precision(5);
      cout<<"thresholds in units of noise variance: Eo="<<Eo<<" Emax="<<Eo*2<<endl;
      // add infos to history
      sprintf(info,"-RES:%d-THR:%g",i,Eo);
      PrintAnalysisInfo(CWB_STAGE_COHERENCE,"cwb2G::Coherence",info,false);
  
      double TL = NET.setVeto(cfg.iwindow);
      cout<<"live time in zero lag: "<<TL<<endl;        // set veto array
      if(TL <= 0.) {froot->Close();EXIT(1);}  	        // exit if live time is zero
  
      if(bplugin) {					// call user plugin
        char cmd[128];
        // export resolution level to CINT  
        sprintf(cmd,"gILEVEL = %lu;",level); EXPORT(size_t,gILEVEL,cmd);
        CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_XCOHERENCE);
      }


      // init sparse table (used in supercluster stage : set the TD filter size)
      pwdm[i]->setTDFilter(cfg.TDSize, 1); 
      for(int n=0; n<nIFO; n++) {
         WS[n].Forward(*hot[n],*pwdm[i]);
         vSS[n].SetMap(&WS[n]);
         vSS[n].SetHalo(mTau);
         if(singleDetector) {
           WS[1]=WS[0];
           vSS[1].SetMap(&WS[1]);
           vSS[1].SetHalo(mTau);
           break;
         }
      }

      // select pixels
      if(cfg.simulation) {cout<<"ifactor|clusters|pixels ";cout.flush();}
      else               {cout<<"lag|clusters|pixels ";    cout.flush();}
      int csize_tot=0;int psize_tot=0; 
      for(int j=0; j<(int)NET.nLag; j++) {

         NET.getNetworkPixels(j,Eo);
         pwc = NET.getwc(j);
         if(NET.pattern!=0) {
            NET.cluster(2,3);
            wc.cpf(*(pwc),false);
            wc.select(const_cast<char*>("subrho"),SELECT_SUBRHO);
            wc.select(const_cast<char*>("subnet"),SELECT_SUBNET);
            pwc->cpf(wc,false);
         } else NET.cluster(1,1);
         // store cluster into temporary job file
         int cycle = cfg.simulation ? ifactor : Long_t(pwc->shift);
         pwc->write(jfile,"coherence","clusters",0,cycle);
         pwc->write(jfile,"coherence","clusters",-1,cycle,-(rateANA>>(cfg.l_high-i)));
         cout<<cycle<<"|"<<pwc->csize()<<"|"<<pwc->size()<<" ";cout.flush();
         csize_tot+=pwc->csize(); psize_tot+=pwc->size(); 

         // add core pixels to sparse table
         for(int n=0; n<nIFO; n++) vSS[n].AddCore(n,pwc);

         pwc->clear();
      }
      // write clusters to the job file & free trees memory
      jfile->Write();
      TList* fList = gDirectory->GetList();
      TObject *obj;
      TIter nextobj(fList);
      while ((obj = (TObject *) nextobj())) {
        if(TString(obj->GetName()).Contains("clusters")) delete obj;
      }
      // update sparse table (halo added to core pixels)
      for(int n=0; n<nIFO; n++) {vSS[n].UpdateSparseTable();}
      // write sparse table to job file
      sdir->cd();
      for(int n=0; n<nIFO; n++) {
        vSS[n].Clean(); 
        char wsname[32];
        if(cfg.simulation) sprintf(wsname,"%s-level:%d:%d",ifo[n],ifactor,cfg.l_high-i);
        else               sprintf(wsname,"%s-level:%d",ifo[n],cfg.l_high-i);
        vSS[n].Write(wsname,TObject::kWriteDelete);
      }
      // print benchmark infos
      watchCoh.Stop(); cout<<endl;
      PrintElapsedTime(watchCoh.RealTime(),"Coherence Elapsed Time for this level : ");
      cout<<endl; watchCoh.Reset();
      // add infos to history
      sprintf(info,"-RES:%d-CSIZE:%d-PSIZE:%d",i,(int)csize_tot/NET.nLag,(int)psize_tot/NET.nLag);
      PrintAnalysisInfo(CWB_STAGE_COHERENCE,"cwb2G::Coherence",info,false);
    }
  }

  if(bplugin) {
    char sfactor[8];sprintf(sfactor,"%d",ifactor);
    CWB_Plugin(jfile,&cfg,&NET,NULL,sfactor,CWB_PLUGIN_OCOHERENCE);
  }

  jfile->Close();

  return;
}

void
cwb2G::SuperCluster(int ifactor) {
//
// Multi resolution clustering & Rejection of the sub-threshold clusters
//
// Loop over time lags 
// - Read clusters from job file (netcluster::read)
// - Multi resolution clustering (netcluster::supercluster)
// - Compute for each pixel the time delay amplitudes (netcluster::loadTDampSSE)
// - Rejection of the sub-threshold clusters (network::subNetCut)
// - Defragment clusters (netcluster::defragment)
// - Store superclusters to job file (netcluster::write)
// Build & Write to job file the sparse TF maps (WriteSparseTFmap)
//

  // if 2G analysis istage>=CWB_STAGE_SUPERCLUSTER then skip super cluster analysis
  if(istage>=CWB_STAGE_SUPERCLUSTER) return; 

  PrintStageInfo(CWB_STAGE_SUPERCLUSTER,"cwb2G::SuperCluster");

  char info[256];

  // open input job file
  TFile* ifile = ((istage==CWB_STAGE_COHERENCE)&&(iname!="")) ? new TFile(iname) : NULL;
  // open temporary job file
  jfile = new TFile(jname,"UPDATE");
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb2G::SuperCluster - Error opening : " << jname <<  endl;EXIT(1);}

  // create prod/sim/lag directories in the output root file to store dignostic histograms
  if(froot==NULL) {cout << "cwb2G::SuperCluster - Error opening : " << froot->GetPath() << endl;EXIT(1);}
  TDirectory* wdir = (TDirectory*)froot->Get("histogram");
  if(wdir==NULL) wdir = froot->mkdir("histogram");

  // decrease skymap resolution to improve subNetCut performances
  double skyres=0;
  if(cfg.healpix) skyres = cfg.healpix>MIN_SKYRES_HEALPIX ? MIN_SKYRES_HEALPIX : 0;
  else            skyres = cfg.angle<MIN_SKYRES_ANGLE ? MIN_SKYRES_ANGLE : 0;
  if(skyres) {
    if(cfg.healpix) NET.setSkyMaps(int(skyres));
    else            NET.setSkyMaps(skyres,cfg.Theta1,cfg.Theta2,cfg.Phi1,cfg.Phi2);
    NET.setAntenna();
    NET.setDelay(cfg.refIFO);
    // the down resampling of the skymask works only for the built-in skymask
    if(strlen(cfg.skyMaskFile)>0)   SetSkyMask(&NET,&cfg,cfg.skyMaskFile,'e',skyres);
    if(strlen(cfg.skyMaskCCFile)>0) SetSkyMask(&NET,&cfg,cfg.skyMaskCCFile,'c',skyres);
  }

  for(int i=0; i<nIFO; i++) pTF[i] = pD[i]->getTFmap();
  // set low-rate TD filters 
  for(int k=0;k<nRES;k++) pwdm[k]->setTDFilter(cfg.TDSize, 1); 
  // read sparse map from job file
  cout << "Loading sparse TF map ... " << endl;
  for(int n=0; n<nIFO; n++) {
    pD[n]->sclear();   // clear vector with sparse maps
    for(int i=0; i<nRES; i++) {
      char swname[32];
      if(cfg.simulation) sprintf(swname,"csparse/%s-level:%d:%d",ifo[n],ifactor,i+cfg.l_low);
      else               sprintf(swname,"csparse/%s-level:%d",ifo[n],i+cfg.l_low);
      SSeries<double>* psw;
      if(ifile!=NULL) psw = (SSeries<double>*)ifile->Get(swname);	
      else            psw = (SSeries<double>*)jfile->Get(swname);	
      if(psw==NULL) {
        cout << "cwb2G::SuperCluster : sparse map " << swname
             << " not exist in job file" << endl;EXIT(1);
      }
      SSeries<double> SS = *psw;
      pD[n]->vSS.push_back(SS);
      delete psw;
    }
    cout<<endl;
  }

  // export to CINT ifile (used by plugins)
  char cmd[128];
  if(ifile!=NULL) sprintf(cmd,"gIFILE = (void*)%p;",ifile);
  else            sprintf(cmd,"gIFILE = NULL;");
  EXPORT(void*,gIFILE,cmd);
  // in supercluster plugin
  if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_ISUPERCLUSTER);

  if(!cfg.scPlugin) {    		 // built in supercluster function 

    int nevt = 0;
    int nnn = 0;
    int mmm = 0;
    size_t count = 0;
    netcluster  wc;
    netcluster* pwc;

    for(int j=0; j<(int)lags; j++) {

      int cycle = cfg.simulation ? ifactor : Long_t(NET.wc_List[j].shift);

      // read cluster metadata
      if(ifile!=NULL) wc.read(ifile,"coherence","clusters",0,cycle);	
      else            wc.read(jfile,"coherence","clusters",0,cycle);	
      // read clusters from temporary job file, loop over TF resolutions
      if(ifile!=NULL) {
        for(int i=nRES-1; i>=0; i--)     // reverse loop is faster loading cluster (?)
          wc.read(ifile,"coherence","clusters",-2,cycle,-(rateANA>>(i+cfg.l_low))); 
      } else {           
        for(int i=nRES-1; i>=0; i--)     // reverse loop is faster loading cluster (?)
          wc.read(jfile,"coherence","clusters",-2,cycle,-(rateANA>>(i+cfg.l_low))); 
      }
      cout<<"-----------------------------------------------------"<<endl; 
      char cycle_name[32];
      if(cfg.simulation) sprintf(cycle_name," factor[%d]=%g",ifactor,cfg.factors[ifactor]);
      else               sprintf(cycle_name," lag=%d",cycle); 
      cout<<"-> Processing "   <<cycle_name<<" ..."<<endl;
      cout<<"   --------------------------------------------------"<<endl; 
      cout<<"   coher  clusters|pixels      : "
          <<setfill(' ')<<setw(6)<<wc.csize()<<"|"<<wc.size()<<endl;

      // supercluster analysis
      if(cfg.l_high==cfg.l_low) wc.pair=false;		// if only one resolution is used pair is false 
      if(NET.pattern!=0) wc.pair=false;                 // if other than pattern=0 - allow one resolution cluster
      wc.supercluster('L',NET.e2or,cfg.TFgap,false);  	// likehood2G
      cout<<"   super  clusters|pixels      : "
          <<setfill(' ')<<setw(6)<<wc.esize(0)<<"|"<<wc.psize(0)<<endl;

      // defragmentation for pattern != 0
      if(NET.pattern!=0) {
         wc.defragment(cfg.Tgap,cfg.Fgap);                                
         cout<<"   defrag clusters|pixels      : "
             <<setfill(' ')<<setw(6)<<wc.esize(0)<<"|"<<wc.psize(0)<<"\n";
      }

      // copy selected clusters to network
      pwc = NET.getwc(j);
      pwc->cpf(wc, false);

      // apply subNetCut() only for pattern=0 || cfg.subnet>0 || cfg.subcut>0
      if(NET.pattern==0 || cfg.subnet>0 || cfg.subcut>0) {
         NET.setDelayIndex(hot[0]->rate());
         pwc->setcore(false);   
         int psel = 0;
         while(1) {
           count = pwc->loadTDampSSE(NET, 'a', cfg.BATCH, cfg.LOUD);
           psel += NET.subNetCut((int)j,cfg.subnet,cfg.subcut,NULL);
           int ptot = pwc->psize(1)+pwc->psize(-1);
           double pfrac = ptot>0 ? double(psel)/double(ptot) : 0.;
           //cout<<"selected pixels: "<<psel<<", fraction: "<<pfrac<< endl;
           if(count<10000) break;
         }
         cout<<"   subnet clusters|pixels      : "
             <<setfill(' ')<<setw(6)<<NET.events()<<"|"<<pwc->psize(-1)<<"\n";
      }
      if(NET.pattern==0) {
         // defragmentation
         pwc->defragment(cfg.Tgap,cfg.Fgap);    
         cout<<"   defrag clusters|pixels      : "
             <<setfill(' ')<<setw(6)<<NET.events()<<"|"<<pwc->psize(-1)<<"\n";
      }

      nevt += NET.events();
      nnn  += pwc->psize(-1);
      mmm  += pwc->psize(1)+pwc->psize(-1);
  
      // store cluster into temporary job file [NEWSS]
      pwc->write(jfile,"supercluster","clusters",0,cycle);
      pwc->write(jfile,"supercluster","clusters",-1,cycle);
      //cout<<cycle<<"|"<<pwc->csize()<<"|"<<pwc->size()<<" ";cout.flush();
  
      pwc->clear();
      cout<<endl;cout.flush();
    }

    // print final statistic
    cout<<endl<<"Supercluster done"<<endl;  
    if(mmm) cout<<"total  clusters|pixels|frac : "
                <<setfill(' ')<<setw(6)<<nevt<<"|"<<nnn<<"|"<<nnn/double(mmm)<<"\n";
    else    cout<<"total  clusters             : "<<nevt<<"\n";
    cout<<endl;cout.flush();

    // add infos to history
    if(mmm) sprintf(info,"-NEVT:%d-PSIZE:%d-FRAC:%g",nevt,nnn,nnn/double(mmm));
    else    sprintf(info,"-NEVT:%d",nevt);
    PrintAnalysisInfo(CWB_STAGE_SUPERCLUSTER,"cwb2G::SuperCluster",info,false);
  }

  // out supercluster plugin
  if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_OSUPERCLUSTER);

  // restore skymap resolution 
  if(skyres) {
    if(cfg.healpix) NET.setSkyMaps(int(cfg.healpix));
    else            NET.setSkyMaps(cfg.angle,cfg.Theta1,cfg.Theta2,cfg.Phi1,cfg.Phi2);
    NET.setAntenna();
    NET.setDelay(cfg.refIFO);
    if(strlen(cfg.skyMaskFile)>0)   SetSkyMask(&NET,&cfg,cfg.skyMaskFile,'e');
    if(strlen(cfg.skyMaskCCFile)>0) SetSkyMask(&NET,&cfg,cfg.skyMaskCCFile,'c');
  }

  if(ifile!=NULL) {
    for(int i=0; i<nIFO; i++) {
      if((istage>=CWB_STAGE_CSTRAIN)&&(jstage<=CWB_STAGE_LIKELIHOOD)) { 
        LoadWaveforms(ifile, pD[i], ifactor);		// read inj waveforms from input job file
        SaveWaveforms(jfile, pD[i], ifactor);    	// write inj waveforms to output job file
      }
    }
  }

  jfile->Write();
  jfile->Close();

  // clear vector with sparse maps
  for(int n=0; n<nIFO; n++) pD[n]->sclear();   
  // write sparse table to job file
  jfile = new TFile(jname,"UPDATE","",9);	// compression level 9
  WriteSparseTFmap(jfile, ifactor, "sparse", "supercluster");
  jfile->Close();

  // job root file garbage collector
  vector<TString> delObjList;
  // coherence clusters is removed if not set in the jobfOptions
  if(!(jobfOptions&CWB_JOBF_SAVE_COHERENCE)) {
    delObjList.push_back("coherence");
  }
  // cstrains and rms data are removed if not set in the jobfOptions 
  if(!(jobfOptions&CWB_JOBF_SAVE_CSTRAIN)) {
    char cdrms_name[32];sprintf(cdrms_name,"rms/rms-f%d",ifactor);
    char cdcstrain_name[32];sprintf(cdcstrain_name,"cstrain/cstrain-f%d",ifactor);
    delObjList.push_back(cdrms_name);
    delObjList.push_back(cdcstrain_name);
  }
  // sparse maps are removed if not set in the jobfOptions (only for the last factor)
  if(ifactor==cfg.nfactor-1) {  // the last factor
    if(!(jobfOptions&CWB_JOBF_SAVE_CSPARSE)) delObjList.push_back("csparse");
  }
  FileGarbageCollector(jname,"",delObjList);

  return;
}

bool
cwb2G::Likelihood(int ifactor, char* ced_dir, netevent* netburst, TTree* net_tree, char* outDump) {
//
// Event reconstruction & parameters estimation
//
// Read sparse map from job file
// Loop over time lags
// - Read cluster list from job file (netcluster::read)
// - Loop over cluster list
//   - Read pixels (netcluster::read)
//   - Compute for each pixel the time delay amplitudes (netcluster::loadTDampSSE)
//   - Event reconstruction+parameter estimation (network::likelihood2G)
//   - Store event parameters to job file (netevent::output2G)
//   - If(config::cedDump>0) Generate Coherent Event Display (CWB::ced) 
//

  // if 2G analysis istage>=CWB_STAGE_LIKELIHOOD then skip likelihood analysis
  if(istage>=CWB_STAGE_LIKELIHOOD) return false; 

  PrintStageInfo(CWB_STAGE_LIKELIHOOD,"cwb2G::Likelihood");

  char info[256];
  netcluster* pwc;
  int ceddir = 0;  // flag if ced directory exists
  double factor = cfg.factors[ifactor];

  NET.setDelayIndex(TDRate);
  if(singleDetector) NET.setIndexMode(cfg.mode);  // when nIFO=1 exclude duplicate delay configurations

  if(bplugin) {
    // export to CINT variables
    EXPORT(float*,gSLAGSHIFT,TString::Format("gSLAGSHIFT = (float*)%p;",slagShift).Data());
  }

  // open temporary job file
  TString xname = ((istage==CWB_STAGE_SUPERCLUSTER)&&(iname!="")) ? iname : jname;
  jfile = (xname==jname) ? new TFile(xname,"UPDATE") : new TFile(xname);
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb2G::Likelihood - Error : file " << xname.Data() << " not found" <<  endl;EXIT(1);}

  if((istage>=CWB_STAGE_CSTRAIN)&&(jstage<=CWB_STAGE_LIKELIHOOD)) { 
    // restore in detector object the inj waveforms from input job file
    for(int i=0; i<nIFO; i++) LoadWaveforms(jfile, pD[i], ifactor); 	
  }
  if((jobfOptions&CWB_JOBF_SAVE_WFINJ)||(jobfOptions&CWB_JOBF_SAVE_WFREC)) { 
    // save to job file whitened inj/rec waveforms
    TFile* ifile = jfile;
    if(xname!=jname) {
      ifile = new TFile(jname,"UPDATE");
      if(ifile==NULL||!ifile->IsOpen()) {
        cout << "cwb2G::Likelihood - Error : file " << jname << " not found" <<  endl; EXIT(1); }
    }
    int wfSAVE = 0;
    if(jobfOptions&CWB_JOBF_SAVE_WFINJ) wfSAVE+=1;
    if(jobfOptions&CWB_JOBF_SAVE_WFREC) wfSAVE+=2;
    for(int i=0; i<nIFO; i++) SaveWaveforms(ifile, pD[i], ifactor, wfSAVE);   
    ifile->Write();
    if(xname!=jname) ifile->Close();
  }

  // set low-rate TD filters
  for(int k=0; k<nRES; k++) pwdm[k]->setTDFilter(cfg.TDSize, cfg.upTDF);
  NET.setDelayIndex(TDRate);

  // read sparse map from job file
  cout << "Loading sparse TF map ... " << endl;
  for(int n=0; n<nIFO; n++) {
    pD[n]->sclear();   // clear vector with sparse maps 
    for(int i=0; i<nRES; i++) {
      char swname[32];
      if(cfg.simulation) sprintf(swname,"sparse/%s-level:%d:%d",ifo[n],ifactor,i+cfg.l_low);
      else               sprintf(swname,"sparse/%s-level:%d",ifo[n],i+cfg.l_low);
      SSeries<double>* psw = (SSeries<double>*)jfile->Get(swname); 
      if(psw==NULL) {
        cout << "cwb2G::Likelihood : sparse map " << swname 
             << " not exist in job file" << endl;EXIT(1);
      }
      SSeries<double> SS = *psw;
      pD[n]->vSS.push_back(SS);
      delete psw; 
    }
    cout<<endl;
  }

  if(cfg.dump && (!cfg.outPlugin)) // write header to ascii dump file 
    {netburst->dopen(outDump,const_cast<char*>("a"),true); netburst->dclose();}

  int LRETRY  = 0;                 // number of likelihood2G retry for each event 
  int NEVENTS = 0;		   // total recontructed events for all lags
  for(int j=0; j<(int)lags; j++) {

    // init likelihood plugin
    EXPORT(int,gLRETRY,"gLRETRY = 0;")   // needed to avoid error msg when gLRETRY is not defined in the plugin 
    if(bplugin) {
      EXPORT(int,gRET,"gRET = 0;")       // needed to avoid error msg when gRET    is not defined in the plugin 
      char strcycle[8];sprintf(strcycle,"%d",cfg.simulation==0 ? j : ifactor);
      CWB_Plugin(jfile,&cfg,&NET,NULL,strcycle,CWB_PLUGIN_ILIKELIHOOD);
      // get gRET from plugin
      int gRET=0; IMPORT(int,gRET) if(gRET==-1) continue;
    }
    int gLRETRY=0; IMPORT(int,gLRETRY) LRETRY=gLRETRY;                                                          

    int cycle = cfg.simulation ? ifactor : Long_t(NET.wc_List[j].shift);

    pwc = NET.getwc(j);
    pwc->cData.clear();

    // read cluster list & metadata netcluster object  
    vector<int> clist = pwc->read(jfile,"supercluster","clusters",0,cycle);
    //pwc->print();

    // print header
    char cycle_name[32];
    if(cfg.simulation) sprintf(cycle_name," factor[%d]=%g",ifactor,cfg.factors[ifactor]);
    else               sprintf(cycle_name," lag=%d",cycle); 
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-> Processing " << clist.size() << " clusters in" << cycle_name<<endl;
    cout<<"   ----------------------------------------------------"<<endl;cout.flush();

    //int nmax = 1000;              // load tdAmp associated to the first nmax loudest pixels
    int nmax = -1;                  // load all tdAmp
    int npixels = 0;		    // total loaded pixels per lag
    int nevents = 0;		    // total recontructed events per lag
    int nselected_core_pixels = 0;
    int nrejected_weak_pixels = 0;  // remove weak glitches
    int nrejected_loud_pixels = 0;  // remove loud glitches
    for(int k=0;k<(int)clist.size();k++) {	// loop over the cluster list 

      // read pixels & tdAmp into netcluster pwc
      pwc->read(jfile,"supercluster","clusters",nmax,cycle,0,clist[k]);

      // likelihood plugin                                                                                      
      if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_XLIKELIHOOD);                                   

      wavearray<double> cid = pwc->get((char*)"ID",  0,'S',0); // get cluster ID
      int id = size_t(cid.data[cid.size()-1]+0.1);
      pwc->setcore(false,id);
      pwc->loadTDampSSE(NET, 'a', cfg.BATCH, cfg.BATCH);  // attach TD amp to pixels

      int lag = j;

      NET.MRA = true;
      int ID = cfg.cedDump ? -id : 0;
      int selected_core_pixels = 0;
      if(NET.pattern>0) { 
        selected_core_pixels = NET.likelihoodWP(cfg.search, lag, ID, NULL);
      } else { 
        selected_core_pixels = NET.likelihood2G(cfg.search, lag, ID, NULL);
      }
      if(!cfg.outPlugin) { 	// if true then output to root file is provided by the user plugin
        double ofactor=0;
        if(cfg.simulation==4)      ofactor=-factor;
        else if(cfg.simulation==3) ofactor=-ifactor;
        else                       ofactor=factor;
        if(cfg.dump) netburst->dopen(outDump,const_cast<char*>("a"),false);
        netburst->output2G(net_tree,&NET,id,lag,ofactor);
        if(cfg.dump) netburst->dclose();
      } 
      int rejected_weak_pixels = 0;
      int rejected_loud_pixels = 0;

      bool detected = (bool)(NET.getwc(j)->sCuts[k] == -1);

      // print reconstructed event  
      cout<<"   cluster-id|pixels: "<<setfill(' ')<<setw(5)<<clist[k]<<"|"<<pwc->size()-npixels;  
      if(detected) cout << "\t -> SELECTED !!!" << endl;
      else 	   cout << "\t <- rejected    " << endl;
      cout.flush();

      // likelihood clusters are stored into temporary job file if set in the jobfOptions
      if(((k==0)||detected)&&(jobfOptions&CWB_JOBF_SAVE_LIKELIHOOD)) {
        TFile* ifile = jfile;
        if(xname!=jname) {
          ifile = new TFile(jname,"UPDATE");
          if(ifile==NULL||!ifile->IsOpen()) {
            cout << "cwb2G::Likelihood - Error : file " << jname << " not found" <<  endl; EXIT(1); }
        }
        pwc->write(ifile,"likelihood","clusters",0,cycle);
        pwc->write(ifile,"likelihood","clusters",-1,cycle,0,k+1);
        if(detected) cout<<"saved"<<endl;cout.flush();
        ifile->Write();
        if(xname!=jname) ifile->Close();
      }

      if(detected) nevents++;
      if(gLRETRY==0) npixels=pwc->size();
      if(!detected && gLRETRY==LRETRY) {
        gLRETRY=0;
        EXPORT(int,gLRETRY,TString::Format("gLRETRY = %d;",gLRETRY).Data())
      }

      // likelihood plugin
      if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_OLIKELIHOOD);

      // ced dump
      if(detected&&(cfg.cedDump)) {
        TFile* ifile = jfile;
        if(xname!=jname) {
          ifile = new TFile(jname,"UPDATE");
          if(ifile==NULL||!ifile->IsOpen()) {
            cout << "cwb2G::Likelihood - Error : file " << jname << " not found" <<  endl; EXIT(1); }
        }
        ced = NULL;
        if(jobfOptions&CWB_JOBF_SAVE_CED) {
          // save ced to temporary job file
          cout<<"dump ced to job file ..." <<endl;
          TDirectory* cdced = NULL;
          cdced = (TDirectory*)ifile->Get("ced");
          if(cdced == NULL) cdced = ifile->mkdir("ced");
          ced = new CWB::ced(&NET,netburst,cdced);
        } else {
          cout<<"dump ced to disk ..." <<endl;
          ced = new CWB::ced(&NET,netburst,ced_dir);
        }
        switch(istage) { 
        case CWB_STAGE_FULL: 
        case CWB_STAGE_INIT: 
          // use TF map & inj signals
          ced->SetOptions(cfg.simulation,cfg.cedRHO,cfg.inRate);
          for(int n=0; n<nIFO; n++) {pTF[n]->setlow(cfg.fLow); pTF[n]->sethigh(cfg.fHigh);}
          break;
        default:                              
          // use sparse map & inj signals
          ced->SetOptions(cfg.simulation,cfg.cedRHO,cfg.inRate,true);
        }
        if(singleDetector) ced->SetChannelName(cfg.channelNamesRaw[0]);
        bool fullCED = singleDetector ? false : true;
        double ofactor=0;
        // CED plugin
        if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_CED);
        if(cfg.simulation==4)      ofactor=-factor;
        else if(cfg.simulation==3) ofactor=-ifactor;
        else                       ofactor=factor;
        if(ced->Write(ofactor,fullCED)) ceddir = 1;
        delete ced;
        if(xname!=jname) ifile->Close();
      }

      //pwc->print();
      pwc->clean(k+1); // clean time delay amplitudes

      IMPORT(int,gLRETRY)
      if(gLRETRY-- > 0) {      // clear last cluster in pwc
        clist = pwc->read(jfile,"supercluster","clusters",0,cycle);
        if(pwc->cList.size()) {
          int npix = pwc->cList[id-1].size();
          for(int p=0;p<npix;p++) pwc->pList.pop_back();
          pwc->cList.pop_back(); pwc->cData.pop_back(); pwc->sCuts.pop_back();
          pwc->cRate.pop_back(); pwc->cTime.pop_back(); pwc->cFreq.pop_back();
          pwc->sArea.pop_back(); pwc->p_Map.pop_back(); pwc->p_Ind.pop_back();
        }
        k--;
      } else gLRETRY=LRETRY;   // every new event gLRETRY is reinitialize to LRETRY
      EXPORT(int,gLRETRY,TString::Format("gLRETRY = %d;",gLRETRY).Data())
    }  // end loop over the cluster list
    /*
    // print final lag statistic
    if(nevents==0) 
      cout << "--------------------------------------------------------------" << endl;
    cout<<endl; 
    cout<<"   total pixels         : "<<pwc->size()<<endl;
    cout<<"   selected_core_pixels : "<<nselected_core_pixels<<endl;
    cout<<"   rejected_weak_pixels : "<<nrejected_weak_pixels<<endl;
    cout<<"   rejected_loud_pixels : "<<nrejected_loud_pixels <<endl;
    cout<<"-> reconstruct events   : "<<nevents<<endl;
    cout<<endl; 
    */
    cout<<endl<<endl;cout.flush();
    NEVENTS+=nevents;
  } // end loop over lags  

  // add infos to history
  sprintf(info,"-NEVT:%d",NEVENTS);
  PrintAnalysisInfo(CWB_STAGE_LIKELIHOOD,"cwb2G::Likelihood",info,false);

  jfile->Close();

  // sparse maps & supercluster clusters are removed if not set in the jobfOptions (only for the last factor)
  if(ifactor==cfg.nfactor-1) {  // the last factor
    vector<TString> delObjList;
    // supercluster clusters are removed if not set in the jobfOptions
    if(!(jobfOptions&CWB_JOBF_SAVE_SUPERCLUSTER)) delObjList.push_back("supercluster");
    if(!(jobfOptions&CWB_JOBF_SAVE_SPARSE)) delObjList.push_back("sparse");
    // If CED is saved to the job file than FileGarbageCollector is not applied
    // FileGarbageCollector do not preserve the style of the plots
    if(!((jobfOptions&CWB_JOBF_SAVE_CED) && cfg.cedDump)) FileGarbageCollector(jname,"",delObjList);
  }

  return ceddir;
}

void
cwb2G::WriteSparseTFmap(TFile* jfile, int ifactor, TString tdir, TString tname) {
//
// Fill sparse map with cluster pixels and write to job file 
//
// jfile   : pointer to job file
// ifactor : factor ID
// tdir    : directory name which contains the sparse map
// tname   : tree name : Ex "coherence" or "supercluster"
//

  // create sparse directories in the job root file to store sparse TF maps
  TDirectory* sdir = (TDirectory*)jfile->Get(tdir);
  if(sdir==NULL) sdir = jfile->mkdir(tdir);

  int ncore_tot=0;
  int ncluster_tot=0;
  int ccluster_tot=0;
  netcluster  wc;
  SSeries<double> ss;
  std::vector<SSeries<double> > vSS;
  for(int n=0; n<nIFO; n++) vSS.push_back(ss);

  for(int i=0; i<nRES; i++) {
    // init sparse table
    pwdm[i]->setTDFilter(cfg.TDSize, 1); 
    for(int n=0; n<nIFO; n++) {
       pTF[n] = pD[n]->getTFmap();
       pTF[n]->Forward(*hot[n],*pwdm[i]);
       vSS[n].SetMap(pTF[n]);
       vSS[n].SetHalo(mTau);
       if(singleDetector) {
         pTF[1] = pD[1]->getTFmap();
         *pTF[1] = *pTF[0];
         vSS[1].SetMap(pTF[1]);
         vSS[1].SetHalo(mTau);
         break;
       }
    }

    for(int j=0; j<(int)lags; j++) {

      int cycle = cfg.simulation ? ifactor : Long_t(NET.wc_List[j].shift);

      // read cluster metadata
      wc.read(jfile,tname,"clusters",0,cycle);

      // read pixels & tdAmp into netcluster pwc                             
      if(tname=="coherence") {
        wc.read(jfile,tname,"clusters",-2,cycle,-vSS[0].wrate()); 
      } else {
        wc.read(jfile,tname,"clusters",-1,cycle,vSS[0].wrate()); 
      }
      //cout<<"loaded clusters|pixels: "<<wc.csize()<<"|"<<wc.size()<<endl;  

      // add core pixels to sparse table
      for(int n=0; n<nIFO; n++) vSS[n].AddCore(n,&wc);
      wc.clear();                                                           
    }                                                                     

    // update sparse table (halo added to core pixels)
    for(int n=0; n<nIFO; n++) {vSS[n].UpdateSparseTable();} 

    // compute sparse statistic
    int ncore=0;	// core pixels
    int ncluster=0;	// core+halo pixels
    int ccluster=0;  // core+halo pixels associated to each core pixels          
    for(int n=0; n<nIFO; n++) {ncore+=vSS[n].GetSparseSize();ncluster+=vSS[n].GetSparseSize(0);}
    ccluster = 2*(vSS[0].GetHaloSlice()+vSS[0].GetHaloSlice(true))+1;	
    ccluster*= 2*vSS[0].GetHaloLayer()+1;
    ncore_tot+=ncore; ncluster_tot+=ncluster; ccluster_tot+=2*ccluster*ncore;
    int rate = rateANA/(rateANA>>(cfg.l_high-i));
    /*
    // print sparse statistic (per resolution)                                                        
    cout << "----------------------- SPARSE STATISTIC / DETECTOR / RES -----------------------" << endl;
    cout << "rate|sSlice|sLayer|eSlice|ccluster : " << rate << "|" << vSS[0].GetHaloSlice() << "|" 
         << vSS[0].GetHaloLayer() << "|" <<vSS[0].GetHaloSlice(1) << "|" << ccluster << endl;
    cout << "npix_core|npix_cluster|3*npix_cluster|2*ccluster*npix_core|ratio|sfilesize(byte) : " << endl;
    cout << ncore << " | " << ncluster  << " | " << 3*ncluster                            
         << " | " << 2*ccluster*ncore << " | " << ((ccluster*ncore)?double(3*ncluster)/(2*ccluster*ncore):0)
         << " | " << 3*ncluster*4 << endl;                                                              
    cout<<endl;
    */
    // write sparse table to job file
    sdir->cd();
    for(int n=0; n<nIFO; n++) {
      if(!(cfg.cedDump && jstage==CWB_STAGE_FULL)) vSS[n].Clean(); // TF map not cleaned if used by ced
      char wsname[32];
      if(cfg.simulation) sprintf(wsname,"%s-level:%d:%d",ifo[n],ifactor,cfg.l_high-i);
      else               sprintf(wsname,"%s-level:%d",ifo[n],cfg.l_high-i);
      vSS[n].Write(wsname,TObject::kWriteDelete);
    }
  }            

  // print final sparse statistic
  cout << "----------------- FINAL SPARSE STATISTIC / DETECTOR -------------------" << endl;
  cout << "npix_core_tot|npix_cluster_tot|3*npix_cluster_tot|ccluster_tot|ratio : " << endl;
  cout << ncore_tot << " | " << ncluster_tot  << " | " << 3*ncluster_tot                          
       << " | " << ccluster_tot << " | " << (ccluster_tot?double(3*ncluster_tot)/(ccluster_tot):0) << endl;        

  return;                      
}                              

void
cwb2G::FillSparseTFmap(TFile* jfile, int ifactor, TString tname) {
//
// Fill sparse map with cluster pixels  
//
// jfile   : pointer to job file
// ifactor : factor ID
// tname   : tree name : Ex "coherence" or "supercluster"
//

  netcluster  wc;

  SSeries<double> SS; 
  for(int n=0; n<nIFO; n++) pD[n]->sclear();   // clear vector with sparse maps

  cout << "cwb2G::FillSparseTFmap ..." << endl;
  for(int i=0; i<nRES; i++) {

    // init sparse table
    for(int n=0; n<nIFO; n++) {
       pTF[n] = pD[n]->getTFmap();
       pTF[n]->Forward(*hot[n],*pwdm[i]);
       pD[n]->vSS.push_back(SS);
       pD[n]->vSS[i].SetMap(pTF[n]);
       pD[n]->vSS[i].SetHalo(mTau);
       if(singleDetector) {
         pTF[1] = pD[1]->getTFmap();
         *pTF[1] = *pTF[0];
         pD[1]->vSS.push_back(SS); 
         pD[1]->vSS[i].SetMap(pTF[1]);
         pD[1]->vSS[i].SetHalo(mTau);
         break;
       }
    }

    for(int j=0; j<(int)lags; j++) {

      int cycle = cfg.simulation ? ifactor : Long_t(NET.wc_List[j].shift);

      // read cluster metadata
      wc.read(jfile,"coherence","clusters",0,cycle);

      // read pixels & tdAmp into netcluster pwc                             
      if(tname=="coherence") {
        wc.read(jfile,tname,"clusters",-2,cycle,-pD[0]->vSS[i].wrate()); 
      } else {
        wc.read(jfile,tname,"clusters",-1,cycle,pD[0]->vSS[i].wrate()); 
      }
      //cout<<"loaded clusters|pixels: "<<wc.csize()<<"|"<<wc.size()<<endl;  

      // add core pixels to sparse table
      for(int n=0; n<nIFO; n++) pD[n]->vSS[i].AddCore(n,&wc);

      wc.clear();                                                           
    }                                                                     
    // update sparse table (halo added to core pixels)
    for(int n=0; n<nIFO; n++) {pD[n]->vSS[i].UpdateSparseTable();} 
  }            
  return;                      
}                              

void   
cwb2G::SaveWaveforms(TFile* jfile, detector* pD, int ifactor, int wfSAVE) {
//
// write injected waveforms to job file (only for simulation>0)
//
// jfile   : pointer to job file
// pD      : pointer to detector object
// ifactor : factor ID
// wfSAVE  : 1/2/3
//           1 : save injected waveforms
//           2 : save reconstructed waveforms
//           3 : save injected & reconstructed waveforms
// 

  if(cfg.simulation==0) return;

  if(jfile==NULL) 
    {cout << "cwb2G::SaveWaveforms - NULL input root file : " << jname << endl;EXIT(1);}
  if(wfSAVE==0) 
    {cout << "cwb2G::SaveWaveforms - save mode must be 1/2/3 " << endl;EXIT(1);}

  // create wf directory
  TDirectory* wf = (TDirectory*)jfile->Get("waveform");
  if(wf==NULL) wf=jfile->mkdir("waveform"); 
  char cdwf_name[32];sprintf(cdwf_name,"waveform-f%d",ifactor);
  TDirectory* cdwf = (TDirectory*)wf->Get(cdwf_name);
  if(cdwf==NULL) cdwf=wf->mkdir(cdwf_name);
  cdwf->cd();

  pD->wfsave(wfSAVE);  // save waveforms
  pD->Write(pD->Name,TObject::kSingleKey);
  pD->wfsave(0);       // restore default value

  return;
}

void   
cwb2G::LoadWaveforms(TFile* ifile, detector* pD, int ifactor, int wfSAVE) {
//
// restore in detector object the wf waveforms from input job file (only for simulation>0)
//
// ifile   : pointer to job file
// pD      : pointer to detector object
// ifactor : factor ID
// wfSAVE  : 1/2/3
//           1 : restore injected waveforms
//           2 : restore reconstructed waveforms
//           3 : restore injected & reconstructed waveforms
// 

  if(cfg.simulation==0) return;

  if(ifile==NULL) 
    {cout << "cwb2G::LoadWaveforms - NULL input root file : " << iname << endl;EXIT(1);}
  if(wfSAVE==0) 
    {cout << "cwb2G::LoadWaveforms - save mode must be 1/2/3 " << endl;EXIT(1);}

  // read wf waveforms
  ifile->cd();
  char cdwf_name[32];sprintf(cdwf_name,"waveform/waveform-f%d",ifactor);
  detector* pd = (detector*)ifile->Get(TString(cdwf_name)+TString("/")+pD->Name);
  if(pd==NULL) return;

  pd->wfsave(wfSAVE);  // load waveforms
  *pd >> *pD; 

  delete pd;
}

