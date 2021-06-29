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

//#define USER_COHERENCE
//#define DUMP_IFO_EE_TF
//#define DUMP_NET_SEL_MAXEE_TF
#define DUMP_IFO_SEL_MAXEE_TF

std::vector<source> srcList;

void PlotData(WSeries<double>* WS, double factor, int mode, TString tag, TString odir, double segEdge=0);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// Plugin to generate simulated gaussian noise and injected 'on the fly' CBC MDC

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_BRST2.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->dataPlugin=false; // enable read data from frames
    cfg->mdcPlugin=true;   // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_DATA) {  

    CWB::Toolbox TB;

    int seed;
    if(ifo.CompareTo("L1")==0) seed=1001;
    if(ifo.CompareTo("L2")==0) seed=1002;
    if(ifo.CompareTo("L3")==0) seed=1003;

    if(ifo.CompareTo("H1")==0) seed=2001;
    if(ifo.CompareTo("H2")==0) seed=2002;
    if(ifo.CompareTo("H3")==0) seed=2003;

    if(ifo.CompareTo("V1")==0) seed=3001;
    if(ifo.CompareTo("V2")==0) seed=3002;
    if(ifo.CompareTo("V3")==0) seed=3003;

    if(ifo.CompareTo("I1")==0) seed=4001;
    if(ifo.CompareTo("I2")==0) seed=4002;
    if(ifo.CompareTo("I3")==0) seed=4003;

    if(ifo.CompareTo("J1")==0) seed=5001;
    if(ifo.CompareTo("J2")==0) seed=5002;
    if(ifo.CompareTo("J3")==0) seed=5003;

    TString fName;
    if(ifo.CompareTo("L1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("L2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("L3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";

    if(ifo.CompareTo("H1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";

    if(ifo.CompareTo("V1")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";
    if(ifo.CompareTo("V2")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";
    if(ifo.CompareTo("V3")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";

    if(ifo.CompareTo("I1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("I2")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("I3")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";

    if(ifo.CompareTo("J1")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";
    if(ifo.CompareTo("J2")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";
    if(ifo.CompareTo("J3")==0) fName="plugins/strains/LCGT_sensitivity_8khz_one_side.txt";

    int size=x->size();
    double start=x->start();
    TB.getSimNoise(*x, fName, seed, net->nRun);
    x->resize(size);
    x->start(start);

/*
    // dump spectrum
    char file[1024];
    sprintf(file,"%s/sensitivity_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 8, cfg->segEdge);
    if(TString(ifo).CompareTo("V1")==0) exit(0);
*/
  }

  if(type==CWB_PLUGIN_MDC) {  

    char cmd[128];
    sprintf(cmd,"network* net = (network*)%p;",net);
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

    if(ifo.CompareTo(net->ifoName[0])==0) {
      net->mdcList.clear();
      net->mdcType.clear();
      net->mdcTime.clear();
      net->mdcList=MDC.mdcList;
      net->mdcType=MDC.mdcType;
      net->mdcTime=MDC.mdcTime;
    }

    // save input mdc list
    if(ifo==cfg->ifo[0]) srcList = MDC.srcList;

    cout.precision(14);
    for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;
  }

  if(type==CWB_PLUGIN_WHITE) {  
/*
    CWB::Toolbox TB;
    int level=x->getLevel(); 
    x->Inverse(-1);
    // dump spectrum
    char file[1024];
    sprintf(file,"%s/sensitivity_white_%s_%d_%s_job%lu.txt",cfg->dump_dir,ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump Sensitivity : " << file << endl << endl;
    TB.makeSpectrum(file, *x, 8, cfg->segEdge);
    if(TString(ifo).CompareTo("V1")==0) exit(0);
    x->Forward(level);
*/
  }

  if(type==CWB_PLUGIN_STRAIN_AND_MDC) {  

    if(ifo!=cfg->ifo[0]) return;

    CWB::Toolbox TB;

    // get rescaled strain from mdcList & update srcList
    // to be used when simulation=2
    int K = net->mdcList.size();
    for(int k=0; k<(int)K; k++) {
      TString slog = TB.GetMDCLog(net->mdcList[k], 1);
      double strain = slog.Atof();
      srcList[k].hrss = strain*cfg->factors[0];
    }

    // init MDC object
    char cmd[128];
    sprintf(cmd,"network* net = (network*)%p;",net);
    gROOT->ProcessLine(cmd);
    CWB::mdc MDC(net);
    cfg->configPlugin.Exec();
    IMPORT(CWB::mdc,MDC) 
    MDC.srcList = srcList;

    // write MDC list file
    char logFile[512];
    int runID = net->nRun;
    int Tb=x->start()+cfg->segEdge;
    int dT=x->size()/x->rate()-2*cfg->segEdge;
    sprintf(logFile,"%s/log_%d_%d_%s_job%d.lst",cfg->output_dir,int(Tb),int(dT),cfg->data_label,runID);
    cout << "Dump : " << logFile << endl;
    MDC.DumpLog(logFile);
  }


  if(type==CWB_PLUGIN_ICOHERENCE) {

#ifdef USER_COHERENCE

    cout << "type==CWB_PLUGIN_ICOHERENCE" << endl;
    int nIFO = net->ifoListSize();                
    detector* pD[NIFO_MAX];                              // pointers to detectors
    for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);                              
    double dTau=net->getDelay(const_cast<char*>(""));    // time delay difference

    char tdf00[512];
    for(int i=cfg->levelD; i>=cfg->l_low; i--) {  	 // loop over TF resolutions

#ifdef DUMP_IFO_EE_TF
      for(int n=0; n<nIFO; n++) {                        // produce TF maps with max over the sky energy
        char tag[32];sprintf(tag,"%s_EE_TF",net->ifoName[n]);
        PlotData(pD[n]->getTFmap(), cfg->factors[ifo.Atoi()], 1, tag, cfg->dump_dir, cfg->segEdge);
      }                                                                                              
#endif
                                                                         
      if(i<=cfg->l_high) {                                                  
                                                                         
        sprintf(tdf00,"%s/data64_wat-4.8.2/Meyer1024wat482_00_L%1d.dat",cfg->filter_dir,i);
        net->setDelayFilters(tdf00);                                                       
        if(i==cfg->l_high) {                                                               
          net->setDelayIndex();                                                            
          net->setIndexMode(1);                                                            
        }
  
        double Ao = net->threshold(cfg->bpp,dTau);
        net->set2or(cfg->x2or*Ao*Ao);      
        cout<<"pixel threshold in units of noise rms: "<<Ao<<endl;
        cout<<"2 OR  threshold in units of noise var: "<<cfg->x2or*Ao*Ao<<endl;
                                                                            
        cout<<"total    pixels: "<<net->coherence(Ao)<<"  ";                   
                                                                            
//        int nn = size_t(2.*cfg->Tgap*pD[0]->getTFmap()->resolution(0)+0.1);         
//        int mm = size_t(cfg->Fgap/pD[0]->getTFmap()->resolution(0)+0.0001);         
//        cout<<"clusters: "<<net->cluster(nn,mm)<<"  ";                           
//        cout<<"selected pixels: "<<net->likelihood('E',cfg->Acore)<<"\n";       
  
        for(int lag=0; lag<(int)net->nLag; lag++) {

          netcluster* pwc = net->getwc(lag);

          // plot selected pixels
          WSeries<double> WSE = *pD[0]->getTFmap();
          WSE=0;                                   
          int V = net->wc_List[lag].pList.size();  
          cout << "SELECTED PIXELS : " << V << endl;
          for(int j=0;j<V;j++) {                  // pixel loop
            netpixel* pix = &(net->wc_List[lag].pList[j]);     
            WSE.data[pix->time]=pix->likelihood;               
            //cout << j << " " << pix->time << " " << pix->likelihood << endl;
/*                                                           
            for(int n=0;n<nIFO;n++) {                          
              size_t index = size_t(pix->getdata('I',n)+0.1);  
              WSE.data[index]+=pix->getdata('S',n);            
            }                                                  
*/                                                           
          }                                                    

#ifdef DUMP_NET_SEL_EE_TF
          char fname[1024];                                                                              
          sprintf(fname,"%s/%s_wdm_scalogram_%d_lev%d_%g.root",                                          
                  cfg->dump_dir, "NET_SEL_MAXEE_TF", int(WSE.start()),WSE.getLevel(),cfg->factors[ifo.Atoi()]);
          cout << endl << "Dump WDM Scalogram : " << fname << endl << endl;                                    
          PlotWSeries(&WSE, fname);                                                                            

//          PlotData(&WSE, cfg->factors[ifo.Atoi()], 2, "NET_SEL_MAXEE_TF", cfg->dump_dir, cfg->segEdge);
#endif

#ifdef DUMP_IFO_SEL_MAXEE_TF
          WSE=0;                                   
          for(int n=0;n<nIFO;n++) {                          
            for(int j=0;j<V;j++) {                  // pixel loop
              netpixel* pix = &(net->wc_List[lag].pList[j]);     
              size_t index = size_t(pix->getdata('I',n)+0.1);  
              WSE.data[index]+=pix->getdata('S',n);            
            }                                                  
            char tag[32];sprintf(tag,"%s_SEL_MAXEE_TF",net->ifoName[n]);
            PlotData(&WSE, cfg->factors[ifo.Atoi()], 2, tag, cfg->dump_dir, cfg->segEdge);
          }                                                    
#endif

          pwc->clear();
        }
        cout<<endl;
      }

      if(i>cfg->l_low) net->Inverse(1);
      gSystem->Exec("/bin/date"); 

    }

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


