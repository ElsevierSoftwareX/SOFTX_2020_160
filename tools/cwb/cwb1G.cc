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


#include "cwb1G.hh"

#define EXIT(ERR) gSystem->Exit(ERR)   // better exit handling for ROOT stuff

ClassImp(cwb1G)

cwb1G::~cwb1G() {
//
// Destructor
//
}

void
cwb1G::Init() {
//
// Initialize & Check variables
//

  // check & set analysis stage
  if(TString(cfg.analysis)!="1G")
    {cout << "cwb1G::Init - Error : analysis must be 1G" << endl;EXIT(1);}
  if(jstage!=CWB_STAGE_FULL)
    {cout << "cwb1G::Init - Error : if analysis=1G then stage must be CWB_STAGE_FULL" << endl;EXIT(1);}

  // Check if lagStep is a multiple of the max time resolution
  // This condition is necessary to ensure a shift of an integer
  // number of pixels when circular buffer is used for lag shift
  int rate_min  = rateANA>>cfg.l_high;
  double dt_max = 1./rate_min;
  if((cfg.lagStep*rate_min-TMath::Nint(cfg.lagStep*rate_min))>1e-12) {
    cout << "cwb1G::Init - Error : lagStep is not a multple of max-time-resolution" << endl;
    cout << "lagStep(sec) : " << cfg.lagStep << "\t max dt(sec) : " << dt_max << endl << endl;
    EXIT(1);
  }

  return;
}

double
cwb1G::ReadData(double mdcShift, int ifactor) {
//
// Read Noise & MDC data from frame file or from the "On The Fly" generator
//

  // data are loaded from root file
  if(iname!="") return cwb::ReadData(iname);

  PrintStageInfo(CWB_STAGE_STRAIN,"cwb1G::ReadData");

  // data are loaded from root file
  if(iname!="") return cwb::ReadData(iname);

  Meyer<double> B(1024);           // set wavelet for resampling
  Meyer<double> S(1024,2);         // set wavelet for production
 
  wavearray<double> x,z;           // temporary time series
  std::vector<wavearray<double> > y; // temporary time series for snr mode
  y.resize(nIFO);
  WSeries<double> wM;              // mdc WSeries

  int    xsize=0.;
  double xrate=0.;
  double xstart=0.;

  // data are loaded from frame files
  jfile = new TFile(jname,"UPDATE");
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb1G::ReadData - Error opening root file : " << jname << endl;EXIT(1);}
  TDirectory* cdstrain = (TDirectory*)jfile->Get("strain");
  if(cdstrain==NULL) cdstrain = jfile->mkdir("strain");

  for(int i=0; i<nIFO; i++) {
    if(cfg.dataPlugin) {
      x.rate(cfg.inRate); x.start(FRF[i].start); x.resize(int(x.rate()*(FRF[i].stop-FRF[i].start)));
    } else {
      fr[i].readFrames(FRF[i],cfg.channelNamesRaw[i],x);  
      if(x.rate()!=cfg.inRate) 
        {cout << "cwb1G::ReadData - input rate from frame " << x.rate() 
              << " do not match the one defined in config : " << cfg.inRate << endl;EXIT(1);}
    }
    if(bplugin) CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)&x,ifo[i],CWB_PLUGIN_STRAIN);
    x.start(x.start()+cfg.dataShift[i]);                // dataShift
    x.start(x.start()-cfg.segLen*(segID[i]-segID[0]));  // SLAG
    if(singleDetector) TB.resampleToPowerOfTwo(x);
    if(cfg.dcCal[i]>0.) x*=cfg.dcCal[i];                // DC correction
    if(cfg.fResample>0) {                               // RESAMPLING
      x.FFTW(1); 
      x.resize(cfg.fResample/x.rate()*x.size()); 
      x.FFTW(-1); 
      x.rate(cfg.fResample); 
    }
    pTF[i] = pD[i]->getTFmap();
    pTF[i]->Forward(x,B,cfg.levelR);
    pTF[i]->getLayer(x,0);
    pTF[i]->Forward(x,S,cfg.levelD);

    // save ifo data to temporary job file
    cdstrain->cd(); pTF[i]->Write(ifo[i]);

    if(i==0) {xrate=x.rate();xstart=x.start();xsize=x.size();}

    fprintf(stdout,"start=%f duration=%f rate=%f\n", x.start(),x.size()/x.rate(),x.rate());
    if(i>0 && xstart != x.start()) {
      cout << "cwb1G::ReadData - Error : ifo noise data not synchronized" << endl;
      cout << ifo[i] << " " << x.start() << " != " << ifo[0] << " " << xstart << endl;
      EXIT(1);
    }
    if(i>0 && xrate != x.rate()) {
      cout << "cwb1G::ReadData - Error : ifo noise data have different rates" << endl;
      cout << ifo[i] << " " << x.rate() << " != " << ifo[0] << " " << xrate << endl;
      EXIT(1);
    }
 
    if(cfg.simulation) {
      TDirectory* cdmdc = (TDirectory*)jfile->Get("mdc");
      if(cdmdc==NULL) cdmdc = jfile->mkdir("mdc");

      if(cfg.mdcPlugin) {
        x.rate(cfg.inRate); x.start(FRF[i+nIFO].start); x.resize(int(x.rate()*(FRF[i+nIFO].stop-FRF[i+nIFO].start)));
      } else {
        fr[nIFO].readFrames(FRF[i+nIFO],cfg.channelNamesMDC[i],x);  
        if(x.rate()!=cfg.inRate) 
          {cout << "cwb1G::ReadData - input rate from frame " << x.rate() 
                << " do not match the one defined in config : " << cfg.inRate << endl;EXIT(1);}
      }

      if(bplugin) CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)&x,ifo[i],CWB_PLUGIN_MDC);
      x.start(x.start()+mdcShift);  
      if(cfg.fResample>0) {                              // RESAMPLING
        x.FFTW(1); 
        x.resize(cfg.fResample/x.rate()*x.size()); 
        x.FFTW(-1); 
        x.rate(cfg.fResample); 
      }
      if(cfg.simulation==2) y[i] = x;  	 	// snr mode
      wM.Forward(x,B,cfg.levelR); 
      wM.getLayer(x,0); 
      wM.Forward(x,S,cfg.levelD);
      if(cfg.simulation==2) {			// snr mode
        pTF[i]->lprFilter(2,0,cfg.Tlpr,4.);
        pTF[i]->setlow(cfg.fLow);
        pD[i]->white(cfg.whiteWindow,1,cfg.segEdge,cfg.whiteStride);
        // set to 0 f<fLow to avoid whitening issues when psd noise is not well defined for f<fLow 
        int layers  = wM.maxLayer()+1;
        for(int j=0;j<layers;j++) if(wM.frequency(j)<cfg.fLow) {wM.getLayer(z,j);z=0;wM.putLayer(z,j);}
        // compute snr
        pD[i]->setsim(wM,NET.getmdcTime(),cfg.iwindow/2.,cfg.segEdge,false);   
      } else {
        cdmdc->cd();wM.Write(ifo[i]);
      }

      fprintf(stdout,"start=%f duration=%f rate=%f\n", x.start(),x.size()/x.rate(),x.rate());
      if(xstart != x.start()) {
        cout << "cwb1G::ReadData - Error : mdc/noise data with different start time" << endl;
        printf("start time : noise = %10.6f - mdc = %10.6f\n",xstart,x.start());
        EXIT(1);
      }
      if(xrate != x.rate()) {
        cout << "cwb1G::ReadData - Error : mdc/noise data with different rate" << endl;
        printf("rate : noise = %10.6f - mdc = %10.6f\n",xrate,x.rate());
        EXIT(1);
      }
      if(xsize != x.size()) {
        cout << "cwb1G::ReadData - Error : mdc/noise data with different buffer size" << endl;
        printf("buffer size : noise = %d - mdc = %lu\n",xsize,x.size());
        EXIT(1);
      }
    }
    if(singleDetector) break;
  }

  // if simulation==2 the factors parameters set the mdc snr
  if(cfg.simulation==2) {
    TDirectory* cdmdc = (TDirectory*)jfile->Get("mdc");
    if(cdmdc==NULL) cdmdc = jfile->mkdir("mdc");

    // compute rescale factor -> snr network=1 
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
    for(int k=0; k<K; k++) {
      if(mdcFactor[k]) cout << k << " mdcFactor : " << mdcFactor[k] << endl;
    }

    // rescale mdc snr network to 1 
    for(int i=0; i<nIFO; i++) {
      pD[i]->setsnr(y[i],NET.getmdcTime(),&mdcFactor,cfg.iwindow/2.,cfg.segEdge);
      if(bplugin) CWB_Plugin(jfile,&cfg,&NET,(WSeries<double>*)&y[i],ifo[i],CWB_PLUGIN_RMDC);

      wM.Forward(y[i],B,cfg.levelR); 
      wM.getLayer(y[i],0); 
      wM.Forward(y[i],S,cfg.levelD);
      cdmdc->cd();wM.Write(ifo[i]);

      y[i].resize(0);
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

  if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_OREADDATA);

  jfile->Close();

  x.resize(0);
  z.resize(0);

  return x.rate();
}

void
cwb1G::DataConditioning(int ifactor) {
//
// Apply line predictor filter to remove lines & whiten data
//

  PrintStageInfo(CWB_STAGE_CSTRAIN,"cwb1G::DataConditioning");

  WSeries<double> wM;              // mdc WSeries
  WSeries<double>* pWS;
  wavearray<double> x;
  TDirectory* cdcstrain=NULL;
  double factor=cfg.factors[ifactor];

  // data are loaded from root file
  //if(!cfg.simulation && iname!="") return cwb::DataConditioning(iname);

  jfile = new TFile(jname, "UPDATE"); 
  if(jfile!=NULL && (cfg.jobfOptions&CWB_JOBF_SAVE_CSTRAIN)) cdcstrain=jfile->mkdir("cstrain");
   
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb1G::DataConditioning - Error : file " << jname << " not found" <<  endl;EXIT(1);}

  for(int i=0; i<nIFO; i++) {
    pWS = (WSeries<double>*)jfile->Get(TString("strain/")+ifo[i]);
    *pTF[i] = *pWS;
    delete pWS;
    
    if(cfg.simulation) {
      pWS = (WSeries<double>*)jfile->Get(TString("mdc/")+ifo[i]); 
      if(cfg.simulation==3) { 			// time shift	: factor is the shift time
        int nshift = int(factor*pWS->rate()); 	// number of shifted samples
        int level = pWS->getLevel();
        pWS->Inverse();
        wM = *pWS; wM=0;
        int jstart = nshift<0 ? -nshift : 0; 
        int jstop  = nshift<0 ? pWS->size() : pWS->size()-nshift; 
        for(int j=jstart;j<jstop;j++) wM.data[j+nshift] = pWS->data[j];
        wM.Forward(level);			// return to the original decomposition level
        pTF[i]->add(wM);

        double tshift=nshift/pWS->rate();        // time shift (sec)
        // take into account of the previous applied time shift
        tshift = ifactor==0 ? tshift : tshift-int(cfg.factors[ifactor-1]*pWS->rate())/pWS->rate();

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
      } else {
        wM = *pWS;
        (*pWS)*=factor;
        pTF[i]->add(*pWS);
      }
      delete pWS;
    }

    if(bplugin) CWB_Plugin(jfile,&cfg,&NET,pTF[i],ifo[i],CWB_PLUGIN_IDATA_CONDITIONING);

    if(!cfg.dcPlugin) {     	// built in data conditioning
      pTF[i]->lprFilter(2,0,cfg.Tlpr,4.);
      pTF[i]->setlow(cfg.fLow);
      pD[i]->white(cfg.whiteWindow,1,cfg.segEdge,cfg.whiteStride);
      if(cfg.simulation) {
        // set to 0 f<fLow to avoid whitening issues when psd noise is not well defined for f<fLow 
        int layers  = wM.maxLayer()+1;
        for(int j=0;j<layers;j++) if(wM.frequency(j)<cfg.fLow) {wM.getLayer(x,j);x=0;wM.putLayer(x,j);}
        // compute mdc params & save whiten mdc 
        pD[i]->setsim(wM,NET.getmdcTime(),cfg.iwindow/2.,cfg.segEdge,true);   
      }
      pTF[i]->Inverse(cfg.levelD-cfg.levelF);
      pTF[i]->lprFilter(2,0,cfg.Tlpr,4.);
      pTF[i]->Forward(cfg.levelD-cfg.levelF);
      pTF[i]->sethigh(cfg.fHigh);
      v[i] = pTF[i]->variability();
      pD[i]->bandPass1G();      // band pass filtering
    } else {                    // data conditioning is provided by the user plugin
      char cmd[128];
      // export to CINT variables
      sprintf(cmd,"gMDC = %p;",&wM); EXPORT(void*,gMDC,cmd);
      if(bplugin) CWB_Plugin(jfile,&cfg,&NET,pTF[i],ifo[i],CWB_PLUGIN_DATA_CONDITIONING);
    }

    if(bplugin) CWB_Plugin(jfile,&cfg,&NET,pTF[i],ifo[i],CWB_PLUGIN_ODATA_CONDITIONING);

    if(cfg.jobfOptions&CWB_JOBF_SAVE_CSTRAIN) {cdcstrain->cd();pTF[i]->Write(ifo[i]);}
          
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
  jfile->Close();

  x.resize(0);

  // strains and mdc data are removed if not set in the jobfOptions (only for the last factor)
  if(ifactor==cfg.nfactor-1) {  // the last factor
    vector<TString> delObjList;
    if(!(cfg.jobfOptions&CWB_JOBF_SAVE_STRAIN)) delObjList.push_back("strain");
    if(cfg.simulation && !(cfg.jobfOptions&CWB_JOBF_SAVE_MDC)) delObjList.push_back("mdc");
    FileGarbageCollector(jname,"",delObjList);
  }

  return;
}

void
cwb1G::Coherence(int ifactor) {
//
// Set pixel energy threshold
// Select the significant pixels
// Single level clustering
//

  PrintStageInfo(CWB_STAGE_COHERENCE,"cwb1G::Coherence");

  int n,m;
  char tdf00[1024];
  double Ao;
  netcluster wc;

  double TL = NET.setVeto(cfg.iwindow);
  cout<<"live time in zero lag: "<<TL<<endl<<endl;  // set veto array
  if(TL <= 0.) {froot->Close();EXIT(1);}  	    // exit if live time is zero

  jfile = new TFile(jname,"UPDATE"); 
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb1G::Coherence - Error : file " << jname << " not found" <<  endl;EXIT(1);}

  if(bplugin) {
    char sfactor[8];sprintf(sfactor,"%d",ifactor);
    CWB_Plugin(jfile,&cfg,&NET,NULL,sfactor,CWB_PLUGIN_ICOHERENCE);
  }

  for(int i=cfg.levelD; i>=cfg.l_low; i--) {  // loop over TF resolutions
  	  
    if(i<=cfg.l_high) {
  	    	    	    
      sprintf(tdf00,"%s/data64_wat-4.8.2/Meyer1024wat482_00_L%1d.dat",cfg.filter_dir,i);
      NET.setDelayFilters(tdf00);
      if(i==cfg.l_high) {
        NET.setDelayIndex();
        NET.setIndexMode(1);
      }
  
      Ao = NET.threshold(cfg.bpp,dTau);
      NET.set2or(cfg.x2or*Ao*Ao);
      cout<<"pixel threshold in units of noise rms: "<<Ao<<endl;
      cout<<"2 OR  threshold in units of noise var: "<<cfg.x2or*Ao*Ao<<endl;
  	
      cout<<"total    pixels: "<<NET.coherence(Ao)<<"  ";
  	    
      n = size_t(2.*cfg.Tgap*pD[0]->getTFmap()->resolution(0)+0.1);
      m = size_t(cfg.Fgap/pD[0]->getTFmap()->resolution(0)+0.0001);
  	    
      cout<<"clusters: "<<NET.cluster(n,m)<<"  ";
      cout<<"selected pixels: "<<NET.likelihood('E',cfg.Acore)<<"\n";
  
      for(int j=0; j<(int)NET.nLag; j++) {
  	wc = *(NET.getwc(j));
        // write cluster data
        int cycle = cfg.simulation ? ifactor : Long_t(wc.shift);
        wc.write(jfile,"coherence","clusters",0,cycle);
        wc.write(jfile,"coherence","clusters",-1,cycle);
        cout<<wc.csize()<<"|"<<wc.size()<<" ";cout.flush(); 
        wc.clear();
      }
      cout<<endl;
    }
	  
    if(i>cfg.l_low) NET.Inverse(1); 
    gSystem->Exec("/bin/date"); GetProcInfo();
  
  }

  if(bplugin) {
    char sfactor[8];sprintf(sfactor,"%d",ifactor);
    CWB_Plugin(jfile,&cfg,&NET,NULL,sfactor,CWB_PLUGIN_OCOHERENCE);
  }

  jfile->Write();
  jfile->Close();

  return;
}

void 
cwb1G::SuperCluster(int ifactor) {
//
// Multi level clustering 
//

  PrintStageInfo(CWB_STAGE_SUPERCLUSTER,"cwb1G::SuperCluster");

  netcluster wc;

  jfile = new TFile(jname); 
  if(jfile==NULL||!jfile->IsOpen()) 
    {cout << "cwb1G::SuperCluster - Error : file " << jname << " not found" <<  endl;EXIT(1);}

  if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_ISUPERCLUSTER);

  for(int j=0; j<(int)lags; j++) {
    // read clusters from temporary job file
    int cycle = cfg.simulation ? ifactor : Long_t(NET.wc_List[j].shift);
    // read metadata netcluster object
    wc.read(jfile,"coherence","clusters",0,cycle);
    // read cluster objects
    for(int i=cfg.l_low; i<=cfg.l_high; i++) {
      wc.read(jfile,"coherence","clusters",-1,cycle,rateANA>>i);
    }
    if(cfg.l_high==cfg.l_low) wc.pair=false;  // if only one resolution is used pair is false
    int m = wc.supercluster('L',NET.e2or,true);
    netcluster* pwc = NET.getwc(j); pwc->cpf(wc,true);
    cout<<m<<"|"<<pwc->size()<<" ";
    wc.clear();
  }

  if(bplugin) CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_OSUPERCLUSTER);

  jfile->Close();

  // coherence clusters are removed if not set in the jobfOptions (only for the last factor)
  if(ifactor==cfg.nfactor-1) {  // the last factor
    vector<TString> delObjList;
    // coherence clusters are removed if not set in the jobfOptions
    if(!(cfg.jobfOptions&CWB_JOBF_SAVE_COHERENCE)) delObjList.push_back("coherence");
    FileGarbageCollector(jname,"",delObjList);
  }

  return;
}

bool
cwb1G::Likelihood(int ifactor, char* ced_dir, netevent* netburst, TTree* net_tree, char* outDump) {
//
// event reconstruction
// event parameters estimation
//

  PrintStageInfo(CWB_STAGE_LIKELIHOOD,"cwb1G::Likelihood");

  char tdf00[1024], tdf90[1024];

  int ceddir = 0;  // flag if ced directory exists

  for(int i=cfg.l_low; i<=cfg.l_high; i++) {
    sprintf(tdf00,"%s/data64_wat-4.8.2/Meyer1024wat482_00%s_L%1d.dat",cfg.filter_dir,cfg.filter,i);
    sprintf(tdf90,"%s/data64_wat-4.8.2/Meyer1024wat482_90%s_L%1d.dat",cfg.filter_dir,cfg.filter,i);
    NET.setDelayFilters(tdf00,tdf90);

    if(i==cfg.l_low) {
      NET.setDelayIndex();
      NET.setIndexMode(cfg.mode);
    }

    cout<<"selected core pixels: "<<NET.likelihood(cfg.search,cfg.Acore)<<" for level "<<i<<"\n";
    cout<<"rejected weak pixels: "<<NET.netcut(cfg.netRHO,'r',0,1)<<"\n";  // remove weak glitches
    cout<<"rejected loud pixels: "<<NET.netcut(cfg.netCC,'c',0,1)<<"\n";   // remove loud glitches
    cout<<"events in the buffer: "<<NET.events()<<"\n";

    if(cfg.cedDump) {
      CWB::ced *ced = NULL;
      if(cfg.jobfOptions&CWB_JOBF_SAVE_CED) {
        // save ced to temporary job file 
        cout<<"dump ced into "<<jname<<"\n";
        jfile = new TFile(jname,"UPDATE");
        if(jfile==NULL||!jfile->IsOpen()) 
          {cout << "cwb1G::Likelihood - Error : file " << jname << " not found" <<  endl;EXIT(1);}
        TDirectory* cdced = NULL; 
        cdced = (TDirectory*)jfile->Get("ced");
        if(cdced == NULL) cdced = jfile->mkdir("ced"); 
        ced = new CWB::ced(&NET,netburst,cdced);
      } else {
        cout<<"dump ced into "<<ced_dir<<"\n";
        ced = new CWB::ced(&NET,netburst,ced_dir);
      }
      ced->SetOptions(cfg.simulation,cfg.cedRHO,cfg.inRate);
      if(singleDetector) ced->SetChannelName(cfg.channelNamesRaw[0]);
      bool fullCED = singleDetector ? false : true; 
      if(ced->Write(cfg.factors[ifactor],fullCED)) ceddir = 1;
      if(cfg.jobfOptions&CWB_JOBF_SAVE_CED) jfile->Close();
      delete ced;
    }

    if(bplugin) {
      jfile = new TFile(jname);
      if(jfile==NULL||!jfile->IsOpen()) 
        {cout << "cwb1G::Likelihood - Error : file " << jname << " not found" <<  endl;EXIT(1);}
      CWB_Plugin(jfile,&cfg,&NET,NULL,"",CWB_PLUGIN_OLIKELIHOOD);
      jfile->Close();
    }

    if(i<cfg.l_high) NET.Forward(1);
  }

  return ceddir;
}
