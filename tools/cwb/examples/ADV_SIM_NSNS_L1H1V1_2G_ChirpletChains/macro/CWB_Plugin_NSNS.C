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
#include <vector>

#include "macro/PlotWSeries.C"

//#define DUMP_WHITE_WDM_TF_IFO
#define USER_COHERENCE

void PlotData(WSeries<double>* WS, double factor, int mode, TString tag, TString odir, double segEdge=0);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// Plugin to generate simulated gaussian noise and injected 'on the fly' NSNS

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_NSNS.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
//    cfg->dataPlugin=true; // disable read data from frames
    cfg->mdcPlugin=true;  // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_DATA) {  
/*
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
    TB.getSimNoise(*x, fName, seed, net->nRun);
    x->resize(size);
    x->start(start);
*/
/*
    // dump spectrum
    char file[1024];
    sprintf(file,"%s/sensitivity_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 8, cfg->segEdge);
    if(TString(ifo).CompareTo("V1")==0) exit(0);
*/
  }

  if(type==CWB_PLUGIN_MDC || type==CWB_PLUGIN_INIT_JOB) {  

    char cmd[128];
    sprintf(cmd,"network* net = (network*)%p;",net);
    gROOT->ProcessLine(cmd);
    sprintf(cmd,"CWB::config* cfg = (CWB::config*)%p;",cfg);
    gROOT->ProcessLine(cmd);

    CWB::mdc MDC(net);

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

    if(type==CWB_PLUGIN_INIT_JOB) {  
      if(ifo.CompareTo(net->ifoName[0])==0) {
        net->mdcList.clear();
        net->mdcType.clear();
        net->mdcTime.clear();
        net->mdcList=MDC.mdcList;
        net->mdcType=MDC.mdcType;
        net->mdcTime=MDC.mdcTime;
        double Tb = x->start();
        double dT = x->size()/x->rate();
        char log_label[512];
        char tmpFile[1024];
        sprintf(log_label,"%d_%d_%s_job%d",int(Tb),int(dT),cfg->data_label,net->nRun);
        sprintf(tmpFile,"%s/%s-LogTMP.txt",cfg->log_dir,log_label);
        cout << "Write MDC Log : " << tmpFile << endl;
        MDC.DumpLog(tmpFile);  
      }
    }

    cout.precision(14);
    for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;
  }

  if(type==CWB_PLUGIN_CLOSE_JOB) {  
    if(ifo.CompareTo(net->ifoName[0])==0) {
      double Tb = x->start();
      double dT = x->size()/x->rate();
      char log_label[512];
      char tmpFile[1024];
      char outFile[1024];
      sprintf(log_label,"%d_%d_%s_job%d",int(Tb),int(dT),cfg->data_label,net->nRun);
      sprintf(tmpFile,"%s/%s-LogTMP.txt",cfg->log_dir,log_label);
      sprintf(outFile,"%s/%s-Log.txt",cfg->log_dir,log_label);
      char command[1024];
      sprintf(command,"/bin/mv %s %s", tmpFile, outFile);
      gSystem->Exec(command);
    }
  }

  if(type==CWB_PLUGIN_WHITE) {

#ifdef DUMP_WHITE_WDM_TF_IFO

    if(TString(cfg->analysis)!="2G") {
      cout << "plugins/CWB_Plugin_SparseData.C : Error - analysis type must be 2G !!!" << endl;
      gSystem->Exit(1);                                                                        
    }

    wavearray<double> y = *x;                   // copy original data to y (contains the final handled data)
    WSeries<double> w;                          // temporary array for data manipulation                    
    WDM<double>* pwdm = NULL;                                                                               

    for(int level=cfg->l_high; level>=cfg->l_low; level--) {    // loop over levels

      int layers = level>0 ? 1<<level : 0;      // get layers
      cout << "level : " << level << " layers : " << layers << endl;

      pwdm = net->getwdm(layers+1);             // get pointer to wdm transform
      if(pwdm==NULL) {
        cout << "CWB_Plugin_HandleWhiteDataWDM.C : Error - WDM not defined !!!" << endl;
        gSystem->Exit(1);
      }

      w.Forward(y,*pwdm);                       // apply wdm transformation

      WSeries<double>* WS = (WSeries<double>*)(&w);
      char tag[32];sprintf(tag,"%s_TF",ifo.Data());
      PlotData(WS, 0, 2, tag, cfg->dump_dir, cfg->segEdge);

      w.Inverse();                              // inverse transform
      y = w;                                    // copy manipulated data to y
    }
    gSystem->Exit(0);

#endif
  }

  if(type==CWB_PLUGIN_ICOHERENCE) {

#ifdef USER_COHERENCE

    cout << "type==CWB_PLUGIN_ICOHERENCE" << endl;
    int rateANA=cfg->inRate>>cfg->levelR;
    int nIFO = net->ifoListSize();
    detector* pD[NIFO_MAX];                              // pointers to detectors
    for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);
    double mTau=net->getDelay(const_cast<char*>("MAX")); // maximum time delay

    char wname[1024];			
    sprintf(wname,"%s/wdm_res_ifo.root", cfg->dump_dir);
    TFile *wroot = new TFile(wname, "RECREATE");	 // create output root file for wdm TF maps

    int nRES = net->wdmListSize();			 // number of resolution levels
    for(int i=0; i<nRES; i++) {                          // loop over TF resolutions
      // print level infos
      int level=cfg->l_high-i;
      int layers = level>0 ? 1<<level : 0;
      int rate  = rateANA>>level;
      cout << "level : " << level << "\t rate(hz) : " << rate
           << "\t layers : " << layers << "\t df(hz) : " << rateANA/2./double(1<<level)
           << "\t dt(ms) : " << 1000./rate << endl;

      for(int n=0; n<nIFO; n++) {                        // produce TF maps with max over the sky energy
        WDM<double>* pwdm = net->wdmList[i];
        wavearray<double>* hot = pD[n]->getHoT();
        pD[n]->getTFmap()->maxEnergy(*hot,*pwdm,mTau,4);

        char wlabel[32];sprintf(wlabel,"%s:%d",net->ifoName[n],level);
        pD[n]->getTFmap()->Write(wlabel);

        char tag[32];sprintf(tag,"%s_MAXEE_TF",net->ifoName[n]);
//        PlotData(pD[n]->getTFmap(), cfg->factors[ifo.Atoi()], 0, tag, cfg->dump_dir, cfg->segEdge);
      }

      double Eo = net->THRESHOLD(cfg->bpp);              // threshold on pixel energy
      cout<<"thresholds in units of noise variance: Eo="<<Eo<<" Emax="<<Eo*2<<endl;

      // sum ifo energies
      WSeries<double> WSE = *pD[0]->getTFmap();
      for(int n=1; n<nIFO; n++) {                        // produce TF maps with sum of energy over the ifo 
        WSE.add(*pD[n]->getTFmap());
      }
//      PlotData(&WSE, cfg->factors[ifo.Atoi()], 0, "NET_MAXEE_TF", cfg->dump_dir, cfg->segEdge);

      double TL = net->setVeto(cfg->iwindow);
      cout<<"live time in zero lag: "<<TL<<endl<<endl;   // set veto array
      if(TL <= 0.) {gSystem->Exit(1);}                   // exit if live time is zero

      for(int lag=0; lag<(int)net->nLag; lag++) {

        net->getNetworkPixels(lag,Eo,Eo*2);

        netcluster* pwc = net->getwc(lag);
        //net->cluster(1,1);
        //cout<<lag<<"|"<<pwc->csize()<<"|"<<pwc->size()<<" ";cout.flush();

        // plot selected pixels
        WSeries<double> WSE = *pD[0]->getTFmap();
        WSE=0;
        int V = net->wc_List[lag].pList.size();
        cout << "SELECTED PIXELS : " << V << endl;
        for(int j=0;j<V;j++) {			// pixel loop
          netpixel* pix = &(net->wc_List[lag].pList[j]);
          WSE.data[pix->time]=pix->likelihood;
/*
          for(int n=0;n<nIFO;n++) {
            size_t index = size_t(pix->getdata('I',n)+0.1);
            WSE.data[index]+=pix->getdata('S',n);		
          }
*/
        }
//        PlotData(&WSE, cfg->factors[ifo.Atoi()], 0, "NET_SEL_MAXEE_TF", cfg->dump_dir, cfg->segEdge);
        char fname[1024];
        sprintf(fname,"%s/%s_wdm_scalogram_%d_lev%d_%g.root", 
                cfg->dump_dir, "NET_SEL_MAXEE_TF", int(WSE.start()),WSE.getLevel(),cfg->factors[ifo.Atoi()]);
        cout << endl << "Dump WDM Scalogram : " << fname << endl << endl;
//        PlotWSeries(&WSE, fname);

        pwc->clear();
      }
      cout<<endl;
      gSystem->Exec("/bin/date"); 
    }
    wroot->Close();
    
    gSystem->Exit(0);

#endif
  }

  return;
}


void PlotData(WSeries<double>* WS, double factor, int mode, TString tag, TString odir, double segEdge) {

  // if(mode==0) plot "amplitude";
  // if(mode==1) plot "energy";
  // if(mode==2) plot "|amplitude|";

  // Plot WDM Scalogram
  watplot WTS(const_cast<char*>("WTS"));

  //scalogram maps
  double start = WS->start()+segEdge;
  double stop  = WS->start()+WS->size()/WS->rate()-segEdge;
  double flow  = 0;                                
  double fhigh = 2048;                             
  WTS.plot(*WS, mode, start, stop,const_cast<char*>("COLZ"));
  WTS.hist2D->GetYaxis()->SetRangeUser(flow, fhigh);      

  // fill WTS.hist2D with nearest pixels
  int xsize=WTS.hist2D->GetNbinsX();
  int ysize=WTS.hist2D->GetNbinsY();
  cout << "xsize : " << xsize << endl;
  cout << "ysize : " << ysize << endl;

  // dump spectrum
  char fname[1024];
  sprintf(fname,"%s/%s_wdm_scalogram_%d_lev%d_%g.root", odir.Data(), tag.Data(), int(start),WS->getLevel(),factor);
  cout << endl << "Dump WDM Scalogram : " << fname << endl << endl;
  WTS.canvas->cd();
  WTS.canvas->Print(fname);

  return;
}


