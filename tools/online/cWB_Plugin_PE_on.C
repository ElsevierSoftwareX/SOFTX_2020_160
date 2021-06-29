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

#define SAVE_TF_RES_ENERGY              //save TF residual energy (reconstructed -injected)


#define SAVE_EVT_CLUSTER		// save cluster to the event output tree	

void ComputeResidualEnergyOptTF( wavearray<double>* wfREC, network* NET, netevent* EVT, int lag, int ID, int ifoID, double& entf_REC,  double& TFtime_REC, double& sig_TFtime_REC, double& TFfreq_REC, double& sig_TFfreq_REC, double& TFduration, double& TFband);

double ComputeEnergy(WSeries<double> *WS);

void ComputeTFtimefreq(WSeries<double> *WS,  double& TFtime, double& sig_TFtime,  double& TFfreq, double& sig_TFfreq,double& TFduration, double& TFband);


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
    cfg->outPlugin=true;  				// disable built-in output root file

    //    gCUTS.resize(cfg->lagSize);				// resize gCUTS
  }

 if(type==CWB_PLUGIN_ILIKELIHOOD) {

    NET->wfsave=true;                                   // enable waveform save

  }


 cfg->mdcPlugin=true;  
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
 for(int k=0;k<(int)NET->mdcList.size();k++) cout << k << " mdcList " << MDC\
					       .mdcList[k] << endl;
 for(int k=0;k<(int)NET->mdcTime.size();k++) cout << k << " mdcTime " << MDC\
					       .mdcTime[k] << endl;
 for(int k=0;k<(int)NET->mdcType.size();k++) cout << k << " mdcType " << MDC\
					       .mdcType[k] << endl;
}




  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(TString(cfg->analysis)!="2G") {
      cout << "CWB_Plugin_WRC.C -> CWB_PLUGIN_OLIKELIHOOD implemented only for 2G" << endl;
      gSystem->Exit(1);
    }


    //    bool sim = false;
    //if(cfg->simulation>0)  sim=true; 


    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    cout << "-----> CWB_Plugin_WRC.C -> " << " gIFACTOR : " << gIFACTOR << endl;

    //    if(gIFACTOR!=sIFACTOR) {	// reset gCUTS when gIFACTOR change
    //  for(int i=0;i<gCUTS.size();i++) gCUTS[i].clear();
    //  sIFACTOR=gIFACTOR;
    // }

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
    char cEnTf_REC[24]; sprintf(cEnTf_REC,"_EnTf_REC[%1d]/D",nIFO);
    double* EnTf_REC  = new double[nIFO];
    char cEnTfWh_REC[24]; sprintf(cEnTfWh_REC,"_EnTfWh_REC[%1d]/D",nIFO);
    double* EnTfWh_REC  = new double[nIFO];


    // leaf                        
    char cTFtauREC_tot[24]; sprintf(cTFtauREC_tot,"_TFtauREC_tot[%1d]/D",nIFO);
    char cTFfreqREC_tot[24]; sprintf(cTFfreqREC_tot,"_TFfreqREC_tot[%1d]/D",nIFO);
    char cTFtauWhREC_tot[24]; sprintf(cTFtauWhREC_tot,"_TFtauWhREC_tot[%1d]/D",nIFO);
    char cTFfreqWhREC_tot[24]; sprintf(cTFfreqWhREC_tot,"_TFfreqWhREC_tot[%1d]/D",nIFO);

    char csig_TFtauREC_tot[48]; sprintf(csig_TFtauREC_tot,"_sig_TFtauREC_tot[%1d]/D",nIFO);
    char csig_TFfreqREC_tot[48]; sprintf(csig_TFfreqREC_tot,"_sig_TFfreqREC_tot[%1d]/D",nIFO);
    char csig_TFtauWhREC_tot[48]; sprintf(csig_TFtauWhREC_tot,"_sig_TFtauWhREC_tot[%1d]/D",nIFO);
    char csig_TFfreqWhREC_tot[48]; sprintf(csig_TFfreqWhREC_tot,"_sig_TFfreqWhREC_tot[%1d]/D",nIFO);

    

    char cTFbandREC_tot[24]; sprintf(cTFbandREC_tot,"_TFbandREC_tot[%1d]/D",nIFO);
    char cTFdurationREC_tot[24]; sprintf(cTFdurationREC_tot,"_TFdurationREC_tot[%1d]/D",nIFO);
    char cTFbandWhREC_tot[24]; sprintf(cTFbandWhREC_tot,"_TFbandWhREC_tot[%1d]/D",nIFO);
    char cTFdurationWhREC_tot[24]; sprintf(cTFdurationWhREC_tot,"_TFdurationWhREC_tot[%1d]/D",nIFO);


    double* TFtauREC_tot =new double[nIFO];
    double* TFfreqREC_tot =new double[nIFO];
    double* TFtauWhREC_tot =new double[nIFO];
    double* TFfreqWhREC_tot =new double[nIFO];

    double* sig_TFtauREC_tot =new double[nIFO];
    double* sig_TFfreqREC_tot =new double[nIFO];
    double* sig_TFtauWhREC_tot =new double[nIFO];
    double* sig_TFfreqWhREC_tot =new double[nIFO];



    double* TFbandREC_tot =new double[nIFO];
    double* TFdurationREC_tot =new double[nIFO];
    double* TFbandWhREC_tot =new double[nIFO];
    double* TFdurationWhREC_tot =new double[nIFO];


    TTree* net_tree = (TTree *) froot->Get("waveburst");
    if(net_tree!=NULL) {
      EVT = new netevent(net_tree,nIFO);
#ifdef SAVE_EVT_CLUSTER
      net_tree->SetBranchAddress("netcluster",&pwc); 
#endif

#ifdef SAVE_TF_RES_ENERGY
      net_tree->SetBranchAddress("_EnTf_REC", EnTf_REC);     
      net_tree->SetBranchAddress("_EnTfWh_REC", EnTfWh_REC);     

      net_tree->SetBranchAddress("_TFtauREC_tot",TFtauREC_tot);
      net_tree->SetBranchAddress("_TFfreqREC_tot",TFfreqREC_tot);
      net_tree->SetBranchAddress("_TFtauWhREC_tot",TFtauWhREC_tot);
      net_tree->SetBranchAddress("_TFfreqWhREC_tot",TFfreqWhREC_tot);

      net_tree->SetBranchAddress("_sig_TFtauREC_tot",sig_TFtauREC_tot);
      net_tree->SetBranchAddress("_sig_TFfreqREC_tot",sig_TFfreqREC_tot);
      net_tree->SetBranchAddress("_sig_TFtauWhREC_tot",sig_TFtauWhREC_tot);
      net_tree->SetBranchAddress("_sig_TFfreqWhREC_tot",sig_TFfreqWhREC_tot);

      net_tree->SetBranchAddress("_TFbandREC_tot",TFbandREC_tot);
      net_tree->SetBranchAddress("_TFdurationREC_tot",TFdurationREC_tot);
      net_tree->SetBranchAddress("_TFbandWhREC_tot",TFbandWhREC_tot);
      net_tree->SetBranchAddress("_TFdurationWhREC_tot",TFdurationWhREC_tot);
 #endif

     } else {
      EVT = new netevent(nIFO);
       net_tree = EVT->setTree();
 #ifdef SAVE_EVT_CLUSTER
       // add new netcluster leaf
       net_tree->Branch("netcluster","netcluster",&pwc,32000,0); 
 #endif

 #ifdef SAVE_TF_RES_ENERGY
       net_tree->Branch("_EnTf_REC", EnTf_REC,cEnTf_REC);
       net_tree->Branch("_EnTfWh_REC", EnTfWh_REC,cEnTfWh_REC);

       net_tree->Branch("_TFtauREC_tot",TFtauREC_tot, cTFtauREC_tot);
       net_tree->Branch("_TFfreqREC_tot",TFfreqREC_tot,cTFfreqREC_tot);
       net_tree->Branch("_TFtauWhREC_tot",TFtauWhREC_tot,cTFtauWhREC_tot);
       net_tree->Branch("_TFfreqWhREC_tot",TFfreqWhREC_tot,cTFfreqWhREC_tot);

       net_tree->Branch("_sig_TFtauREC_tot",sig_TFtauREC_tot, csig_TFtauREC_tot);
       net_tree->Branch("_sig_TFfreqREC_tot",sig_TFfreqREC_tot,csig_TFfreqREC_tot);
       net_tree->Branch("_sig_TFtauWhREC_tot",sig_TFtauWhREC_tot,csig_TFtauWhREC_tot);
       net_tree->Branch("_sig_TFfreqWhREC_tot",sig_TFfreqWhREC_tot,csig_TFfreqWhREC_tot);

       net_tree->Branch("_TFbandREC_tot",TFbandREC_tot, cTFbandREC_tot);
       net_tree->Branch("_TFdurationREC_tot",TFdurationREC_tot,cTFdurationREC_tot);
       net_tree->Branch("_TFbandWhREC_tot",TFbandWhREC_tot,cTFbandWhREC_tot);
       net_tree->Branch("_TFdurationWhREC_tot",TFdurationWhREC_tot,cTFdurationWhREC_tot);


#endif


    }
  
  
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
          double recTh   = EVT->theta[0];		// rec theta
          double recPh   = EVT->phi[0];			// rec phi



	   wavearray<double>** pwfREC = new wavearray<double>*[nIFO];
	    detector* pd = NET->getifo(0);
	   int idSize = pd->RWFID.size();
	   int wfIndex=-1;
	   for(int mm=0; mm<idSize; mm++) if (pd->RWFID[mm]==ID) wfIndex=mm;
	   // extract whitened injected & reconstructed waveforms
	   for(int n=0; n<nIFO; n++) {
	     if (wfIndex<0) {
	  	 cout << "CWB_Plugin_WRC.C : Error : Reconstructed waveform not saved !!! : ID -> "
	  	      << ID << " : detector " << NET->ifoName[n] << endl;
	  	 continue;
	     }
	       
	        pd = NET->getifo(n);
	        if (wfIndex>=0) pwfREC[n] = pd->RWFP[wfIndex];
	        double R = pd->TFmap.rate();
	      wavearray<double>* wfREC = pwfREC[n];	// array of reconstructed waveforms
	       
	      double entf_REC, tftime_REC,sig_tftime_REC,tffreq_REC,sig_tffreq_REC, tfband_rec, tfduration_rec;
	       
	       //compute residual energy in time domain
	      ComputeResidualEnergyOptTF(  wfREC, NET,EVT, k, ID, n,  entf_REC, tftime_REC,  sig_tftime_REC,  tffreq_REC,  sig_tffreq_REC,tfduration_rec, tfband_rec);
	       EnTfWh_REC[n] = entf_REC;
	       TFtauWhREC_tot[n] = wfREC->start()+tftime_REC;
	       TFfreqWhREC_tot[n] = tffreq_REC;
               sig_TFtauWhREC_tot[n] = sig_tftime_REC;
               sig_TFfreqWhREC_tot[n] = sig_tffreq_REC;
               TFbandWhREC_tot[n] = tfband_rec;
               TFdurationWhREC_tot[n] = tfduration_rec;
	       cout<<"WHITE "<<tfduration_rec<<"  "<<tfband_rec<<endl;
	     }
	   delete [] pwfREC;
	       
	   
	   wavearray<double>** pwfREC_s = new wavearray<double>*[nIFO];
	   //	   detector* pd = NET->getifo(0);
	   // int idSize = pd->RWFID.size();
	   //int wfIndex=-1;
	   for(int mm=0; mm<idSize; mm++) if (pd->RWFID[mm]==-ID) wfIndex=mm;
	   
             for(int n=0; n<nIFO; n++) {
               pd = NET->getifo(n);
	       if (wfIndex<0) {
		 cout << "CWB_Plugin_WRC.C : Error : Reconstructed waveform not \
saved !!! : ID -> "
		      << ID << " : detector " << NET->ifoName[n] << endl;
		 continue;
	       }
	       
	       if (wfIndex>=0) pwfREC_s[n] = pd->RWFP[wfIndex];
	       double R = pd->TFmap.rate();
	      wavearray<double>* wfREC_s = pwfREC_s[n];	// array of reconstructed waveforms
	       
	      double entf_REC, tftime_REC,sig_tftime_REC,tffreq_REC,sig_tffreq_REC, tfband_rec, tfduration_rec;
	      cout<<"STRAIN "<<endl;
	       
	       //compute residual energy in time domain
	      ComputeResidualEnergyOptTF(  wfREC_s, NET,EVT, k, ID, n,  entf_REC, tftime_REC,  sig_tftime_REC,  tffreq_REC,  sig_tffreq_REC,tfduration_rec, tfband_rec);

	      
               EnTf_REC[n] = entf_REC;
               TFtauREC_tot[n] = wfREC_s->start()+tftime_REC;
               TFfreqREC_tot[n] = tffreq_REC;
               sig_TFtauREC_tot[n] = sig_tftime_REC;
               sig_TFfreqREC_tot[n] = sig_tffreq_REC;
               TFbandREC_tot[n] = tfband_rec;
               TFdurationREC_tot[n] = tfduration_rec;
	     }
	    
	 
	     delete [] pwfREC_s;
   
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
      }

    }
    jfile->cd();

    if(pwc)   delete pwc;
     if(EVT)   delete EVT;


        if( TFtauREC_tot) delete []  TFtauREC_tot;
     if( TFfreqREC_tot) delete []  TFfreqREC_tot;
    if( TFtauWhREC_tot) delete []  TFtauWhREC_tot;
    if( TFfreqWhREC_tot) delete []  TFfreqWhREC_tot;

  }
    
       
 return;

  
 

}



void ComputeResidualEnergyOptTF( wavearray<double>* wfREC, network* NET, netevent* EVT, int   lag, int ID, int ifoID, double& entf_REC,  double& TFtime_REC, double& sig_TFtime_REC, double& TFfreq_REC, double& sig_TFfreq_REC, double& TFduration, double& TFband)
 {

  entf_REC=0; 
 // find TF optimal resolution level
  netcluster* pwc = NET->getwc(lag); 
  double optRate  = (pwc->cRate[ID-1])[0];
  double optLayer = pwc->rate/optRate;
  double optLevel = int(log10(optLayer)/log10(2));

  double bREC = wfREC->start();
  double eREC = wfREC->stop();

  double R=wfREC->rate();
  detector* pd = NET->getifo(0);

  //  int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
  //int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);
  //cout << "oINJ : " << oINJ << " oREC : " << oREC << endl;

  // double startXCOR = bINJ>bREC ? bINJ : bREC;
  // double endXCOR   = eINJ<eREC ? eINJ : eREC;
   int sizeXCOR     = int((eREC-bREC)*R+0.5);

  // length of temporary buffer for tf plots
  int xcor_length = sizeXCOR/wfREC->rate();
  int wts_size =xcor_length<8 ? 16 :  16*int(1+xcor_length/8.);
   wts_size*=wfREC->rate();
  //wts_size=wfREC->size();
  cout<<"SIXX "<<sizeXCOR<<"  "<<wts_size<<endl;
 
  char fname[1024];
  Meyer<double> S(1024,2);     // set wavelet for production
  double wts_start,wts_stop;


  //  cout<<wfREC->start()-EVT->gps[0]<<" bbb  "<<bREC<< "   "<<double(bREC+wts_size/2)<<"  "<<wfREC->rate()<<endl;

  // ------------------------------------------------------------------------
    // RECONSTRUCTED WAVEFORM
  // ------------------------------------------------------------------------
  
  wavearray<double> xREC(wts_size);
  xREC.start(wfREC->start()-EVT->gps[0]);
  // xREC.start(wfREC->start()-EVT->gps[0]+double(bREC+wts_size/2)/wfREC->rate());
  xREC.rate(wfREC->rate());
  xREC=0.;
  for (int m=0;m<sizeXCOR;m++)
   if(m<(int)xREC.size()/2) xREC.data[m+wts_size/2]=wfREC->data[m];
  WSeries<double> wREC(xREC,S);
 xREC.resize(0);
  
  if(NET->getwdm(optLayer+1)) wREC.Forward(wREC,*NET->getwdm(optLayer+1));
  cout<<"LLLLL  "<< wREC.start()+(double)(wts_size/2)/wREC.rate() <<"  "<<wts_start+sizeXCOR/wREC.rate()<<endl;
  
  entf_REC = ComputeEnergy(&wREC);
  cout << "TF : " << "enINJ : "<< " enREC : " << entf_REC 
       << endl;
  double tm, fq,stm,sfq, bd, dr;
  ComputeTFtimefreq(&wREC,tm,stm,fq,sfq,dr,bd);
  TFtime_REC= tm;
  sig_TFtime_REC =stm;
  TFfreq_REC =fq;
  sig_TFfreq_REC=sfq;

  TFduration=dr;
  TFband =bd;

  cout<<"TEST freqtime"<<TFfreq_REC<<"  "<< sig_TFfreq_REC<<"  TFtime " <<TFtime_REC<<" band   "<<TFband<< endl;

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



void ComputeTFtimefreq(WSeries<double> *WS,  double& TFtime, double& sig_TFtime, double& TFfreq, double& sig_TFfreq, double& TFduration, double& TFband) {

   int layers = WS->maxLayer()+1;  // numbers of frequency bins (first & last bins  have df/2)                                                                         
   int slices = WS->sizeZero();    // number of time bins                           

   float df = WS->resolution();    // frequency bin resolution (hz)                 
   float dt = 1./(2*df);           // time bin resolution (sec)                     

   int rate = int(1./dt);

   //  cout << "layers : " << layers << "\t slices : " << slices << "\t rate : " << r\
ate                                                                                
  //       << "\t dt : " << dt << "\t df : " << df << endl;                          
  double timeEn =0; double freqEn=0;
 double EE=0;  double EE4=0;
 double pixel_frequency;
 double pixel_time;

 for(int i=1;i<=slices;i++) {
   pixel_time = i*dt;
   for(int j=1;j<=layers;j++) {
     double f = 0;
     double t = 0;
     double en=0;
     if(j==1) {pixel_frequency = df/4; }
     if(j==layers) {pixel_frequency = (layers-1)*df -df/4;}
     if((j!=1)&&(j!=layers)) pixel_frequency = df/2 + (j-1)*df  +df/2;

     en = pow(WS->getSample(i-1,j-1),2);
     f = en*pixel_frequency;
     t = en*pixel_time;
     if(j==1) { EE += en/2.; EE4 += (pow(en,2))/2.; timeEn += t/2; freqEn +=t/2;}
     if(j==layers){ EE += en/2.;EE4 += (pow(en,2))/2.; timeEn +=t/2; freqEn +=t/2 ;}
     if((j!=1)&&(j!=layers)){ EE += en; EE4 += (pow(en,2));timeEn +=t; freqEn +=f;}
     //      if(en !=0) cout<<"AAIUTO "<< pixel_time << " : " << pixel_frequency << "  "<<f<<" "<<t<<" "<<en<<" Sig "<<EE4<<"  "<<endl;                             
   }
 }

 // cout<<"??? "<<timeEn<<"  "<<EE<<"  "<<freqEn<<"  "<<EE4<<endl;
 TFtime = timeEn/EE;
 sig_TFtime = dt*sqrt(1./12.)*sqrt(EE4)/EE;
 TFfreq = freqEn/EE;
 sig_TFfreq= df*sqrt(1./12.)*sqrt(EE4)/EE;

double bandEn =0; double durationEn=0;

 for(int i=1;i<=slices;i++) {
   pixel_time = i*dt;
   for(int j=1;j<=layers;j++) {
     double band = 0;
     double duration = 0;
     double en=0;
     if(j==1) {pixel_frequency = df/4; }
     if(j==layers) {pixel_frequency = (layers-1)*df -df/4;}
     if((j!=1)&&(j!=layers)) pixel_frequency = df/2 + (j-1)*df  +df/2;

     en = pow(WS->getSample(i-1,j-1),2);
     band  = en*(pixel_frequency-TFfreq)*(pixel_frequency-TFfreq);
     duration = en*(pixel_time-TFtime)*(pixel_time-TFtime);
     if(j==1) { EE += en/2.;durationEn +=duration/2;  bandEn += band/2; }
     if(j==layers){ EE += en/2.; durationEn +=duration/2; bandEn +=band/2    ;}
     if((j!=1)&&(j!=layers)){ EE += en; durationEn +=duration; bandEn +=band;}
     //      if(en !=0) cout<<"AAIUTO "<< pixel_time << " : " << pixel_frequency << "  "<<f<<" "<<t<<" "<<en<<" Sig "<<EE4<<"  "<<endl;                             
   }
 }

 TFband = sqrt(bandEn/EE);
 TFduration= sqrt(durationEn/EE);
 cout<<"  BBB DDF "<<TFband <<"  "<<TFduration<<endl;




 return;
  }

