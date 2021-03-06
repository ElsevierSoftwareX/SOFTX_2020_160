#define XIFO 5

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

wavearray<double>* strain[NIFO_MAX];
void WriteFrame(wavearray<double>* strain, wavearray<double>* mdc, TString ifo, CWB::config* cfg, network* net, int frame);

int  FRLEN = 4;   // output frame len (sec) - must be a submultiple of segLen

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

  cout << endl;
  cout << "-----> macro/CWB_Plugin_BRST.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cfg->dataPlugin=true; 	// disable read data from frames
    cfg->mdcPlugin=true;   	// disable read mdc from frames
    for(int n=0;n<cfg->nIFO;n++) strain[n] = new wavearray<double>;
  }

  if(type==CWB_PLUGIN_DATA) {  
    CWB::Toolbox TB;

    int seed;
    if(ifo.CompareTo("L1")==0) seed=1000;
    if(ifo.CompareTo("H1")==0) seed=2000;
    if(ifo.CompareTo("V1")==0) seed=3000;

    TString fName;
    if(ifo.CompareTo("L1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("H1")==0) fName="plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt";
    if(ifo.CompareTo("V1")==0) fName="plugins/strains/advVIRGO_sensitivity_12May09_8khz_one_side.txt";

    int size=x->size();
    double start=x->start();
    TB.getSimNoise(*x, fName, seed, net->nRun);
    x->resize(size);
    x->start(start);

    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    // save strain to tamporary wavearray
    *strain[id] = *x;
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

    IMPORT(int,FRLEN)		// read FRLEN from config plugin

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

    cout.precision(14);
    for(int k=0;k<(int)net->mdcList.size();k++) cout << k << " mdcList " << MDC.mdcList[k] << endl;
    for(int k=0;k<(int)net->mdcTime.size();k++) cout << k << " mdcTime " << MDC.mdcTime[k] << endl;
    for(int k=0;k<(int)net->mdcType.size();k++) cout << k << " mdcType " << MDC.mdcType[k] << endl;

    int nIFO=net->ifoListSize();
    char ifoLabel[64]="";
    for(int i=0;i<nIFO;i++) sprintf(ifoLabel,"%s%s",ifoLabel,net->ifoName[i]);

    // write log when last ifo
    if(ifo==cfg->ifo[cfg->nIFO-1]) {
      char logFile[256];
      sprintf(logFile,"%s/logs/Log-%s-%s-job%d.txt",cfg->output_dir,ifoLabel,cfg->injectionList,net->nRun);
      cout << logFile << endl;
      MDC.DumpLog(logFile);
    }
  }

  if(type==CWB_PLUGIN_RMDC) {  

    wavearray<double>* mdc = x;
    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    double segStart = mdc->start();
    double segEdge  = cfg->segEdge;
    double segLen   = mdc->size()/mdc->rate()-2*segEdge;
    cout.precision(14);
    if(fmod(segStart+segEdge,FRLEN)!=0) 
      {cout << "segment start " << segStart << " is not mod " << FRLEN << " sec !!!" << endl; gSystem->Exit(-1);}
    if(fmod(segLen,FRLEN)!=0) 
      {cout << "segment len " << segLen << " is not mod " << FRLEN << " sec !!!" << endl; gSystem->Exit(-1);}

    // write original segment into nFrame frames
    int nFrame = int(segLen/FRLEN);
    int nEdge  = int(segEdge*mdc->rate());
    int nSize  = mdc->size()-2*nEdge;
    int size4s = FRLEN*mdc->rate();
 
    wavearray<double> strain4s(size4s);
    strain4s.rate(mdc->rate());
    wavearray<double> mdc4s;
    for(int i=0;i<nFrame;i++) {
      //cout << i << " " << mdc->start()+i*FRLEN << endl;
      strain4s.start(mdc->start()+segEdge+i*FRLEN);  
      mdc4s=strain4s;
      for(int j=0;j<size4s;j++) {
         strain4s[j] = strain[id]->data[nEdge+i*size4s+j];
         mdc4s[j]    = mdc->data[nEdge+i*size4s+j];
      }
      int frame =  (net->nRun-1)*nFrame+i;
      WriteFrame(&strain4s, &mdc4s, ifo, cfg, net, frame);
    }

    // exit when last ifo
    // end job - clean temporary files
    if(ifo==cfg->ifo[cfg->nIFO-1]) {
      cout << "remove temporary file ..." << endl;
      TString jname = jfile->GetPath();
      jname.ReplaceAll(":/","");
      cout << jname.Data() << endl;
      gSystem->Exec(TString("rm "+jname).Data());
      cout << "end job" << endl;
      gSystem->Exit(0);
    }
  }

  return;
}

void 
WriteFrame(wavearray<double>* strain, wavearray<double>* mdc, TString ifo, CWB::config* cfg, network* net, int frame) {

  //cout << "WriteFrame " << ifo.Data() << endl;
  double deg2rad = TMath::Pi()/180.;

  int id=-1;
  for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
  if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

  // if frFiles[i+nIFO]==frFiles[i] then strain=starin+factor*mdc
  if(TString(cfg->frFiles[id+cfg->nIFO])=="ADD_TO_STRAIN") {
    wavearray<double> rmdc = *mdc; 
    rmdc *= cfg->factors[0];     // rescaled mdc
    *strain+=rmdc;
  }

  int size=strain->size();
  double start=strain->start();
  double lenght=strain->size()/strain->rate();

  //------------------ create output directory ----------
  char oDir[256];
  sprintf(oDir,"%s/%s",cfg->output_dir,ifo.Data());
  char command[256];
  sprintf(command,"mkdir -p %s",oDir);
  gSystem->Exec(command);

  //------------------ Create a new frame file ----------
  char ofName[256];
  //sprintf(ofName,"%s/%c-%s_llhoft-%d-%d.gwf",oDir,ifo[0],ifo.Data(),int(start),int(lenght));
  sprintf(ofName,"%s/%s-%d-%d.gwf",oDir,cfg->frFiles[id],int(start),int(lenght));
  cout << ofName << endl;
  FrFile *ofp = FrFileONew(ofName,0);

  //------------------ Create a new frame with snr=1 ----
  char frName[256];
  sprintf(frName,"%s-CWB",ifo.Data());
  FrameH* simFrame = FrameNew(frName);
  simFrame->frame = frame;
  simFrame->run = -1;
  simFrame->dt = lenght;
  simFrame->GTimeS = start;
  simFrame->GTimeN = 0;

  //------------------ Set Detector ---------------------
  detector* det = net->getifo(id); 
  detectorParams detParms = det->getDetectorParams();

  FrDetector* detectProc   = simFrame->detectProc;
  FrStrCpy(&detectProc->name, detParms.name);
  detectProc->longitude    = deg2rad*detParms.longitude;
  detectProc->latitude     = deg2rad*detParms.latitude;
  detectProc->elevation    = detParms.elevation;
  detectProc->armXazimuth  = detParms.AzX;
  detectProc->armYazimuth  = detParms.AzY;
  detectProc->armXaltitude = detParms.AltX;
  detectProc->armYaltitude = detParms.AltY;
  detectProc->armXmidpoint = 0;
  detectProc->armYmidpoint = 0;

  //------------------ Set STRAIN channel --------------
  FrProcData* procStrain = FrProcDataNew(simFrame,cfg->channelNamesRaw[id],strain->rate(),strain->size(),-64);
  if(procStrain == NULL) {cout << "Cannot create procStrain" << endl;gSystem->Exit(-1);}
  procStrain->timeOffset = 0;
  procStrain->tRange = simFrame->dt;
  procStrain->type = 1;   // Time Serie
  char strain_unitX[16]="Time";
  FrVectSetUnitX(procStrain->data,strain_unitX);
  char strain_unitY[16]="h";
  FrVectSetUnitY(procStrain->data,strain_unitY);
  for (int i=0;i<procStrain->data->nData;i++) procStrain->data->dataD[i] = strain->data[i];

  //------------------ Set MDC channel -----------------
  FrProcData* procMDC = FrProcDataNew(simFrame,cfg->channelNamesMDC[id],mdc->rate(),mdc->size(),-64);
  if(procMDC == NULL) {cout << "Cannot create procMDC" << endl;gSystem->Exit(-1);}
  procMDC->timeOffset = 0;
  procMDC->tRange = simFrame->dt;
  procMDC->type = 1;   // Time Serie
  char MDC_unitX[16]="Time";
  FrVectSetUnitX(procMDC->data,MDC_unitX);
  char MDC_unitY[16]="Counts";
  FrVectSetUnitY(procMDC->data,MDC_unitY);
  for (int i=0;i<procMDC->data->nData;i++) procMDC->data->dataD[i] = mdc->data[i];

  //------------------ Set DQ channel ------------------
  char dqName[256];
  sprintf(dqName,"%s:LLD-DQ_VECTOR",ifo.Data());
  double dqRate = 1.; 	// 1Hz
  double dqSize = int(simFrame->dt);
  FrAdcData* adcDQ = FrAdcDataNew(simFrame,dqName,1.,dqSize,32);
  if(adcDQ == NULL) {cout << "Cannot create adcDQ" << endl;gSystem->Exit(-1);}
  adcDQ->timeOffset = 0;
  char DQ_unitX[16]="Time";
  FrVectSetUnitX(adcDQ->data,DQ_unitX);
  char DQ_unitY[16]="Counts";
  FrVectSetUnitY(adcDQ->data,DQ_unitY);
  for (int i=0;i<adcDQ->data->nData;i++) adcDQ->data->dataI[i] = 7;	// ifo in science mode

  //------------------ Write Frame --------------------
  int err=FrameWrite(simFrame,ofp);
  if (err) {cout << "Error writing frame" << endl;gSystem->Exit(1);}
  FrameFree(simFrame);

  //------------------ Close File Frame ---------------
  if (ofp!=NULL) FrFileOEnd(ofp);

  return;
}
