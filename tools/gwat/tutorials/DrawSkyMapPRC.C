//
// Draw FOM plots for Position Reconstruction 
// Author : Gabriele Vedovato
// Note : run with ACLiC -> root -l 'DrawSkyMapPRC.C+("IN_FILE_NAME","median50",1)'
// Note : to plot the skymap root files produced by this macro use : 
//        DrawLarsHistogramPRC.C 
//

#define XIFO 4

#pragma GCC system_header

#include "gnetwork.hh"
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
#include "gwat.hh"
#include "mdc.hh"
#include <vector>


//------------------------------------------------------
// General Settings
//------------------------------------------------------

#define WRITE_PLOT
#define DATA_LABEL ""

#define WRITE_SKYMAP

//------------------------------------------------------
// Plot Settings
//------------------------------------------------------

#define DISPLAY_WORLD_MAP

#define DRAW_SITES

//------------------------------------------------------
// Thresholds Settings
//------------------------------------------------------

#define TREE_NAME "waveburst"

#define RHO_THR 2.0
#define CC_THR  0.5
//#define ED_THR  0.35   // ONLY FOR SEARCH I
#define SNR_THRESHOLD  100000000.

//------------------------------------------------------
// Smoothing Settings
//------------------------------------------------------

#define MIN_TRIGGERS_PER_PIXEL 2

//#define MAX_ENTRIES_PER_PIXEL 2000
//#define MAX_ENTRIES_PER_PIXEL 500
//#define MAX_ENTRIES_PER_PIXEL 400
#define MAX_ENTRIES_PER_PIXEL 300
//#define MAX_ENTRIES_PER_PIXEL 200

//------------------------------------------------------
// SKYMAP Settings
//------------------------------------------------------

//char COORDINATES[256]="cWB";
char COORDINATES[256]="Geographic";

//char PROJECTION[256]="";
char PROJECTION[256]="hammer";
//char PROJECTION[256]="sinusoidal";

// default  (for Lars Histogram)
#define ANGULAR_RESOLUTION (1./3.)
#define NPIX_SMOOTH 13

//#define ANGULAR_RESOLUTION 1.
//#define NPIX_SMOOTH 4

//#define ANGULAR_RESOLUTION 3. 
//#define NPIX_SMOOTH 1

//#define ANGULAR_RESOLUTION 9.
//#define NPIX_SMOOTH 0

//------------------------------------------------------


void DrawSkyMapPRC(TString ifName, TString pType="MEDIAN50", int type=0, TString Title="", TString Tag="") {

#ifdef __CINT__
  printf("Please run this script in compiled mode by running \".x DrawSkyMapPRC.C+\"\n");
  return;
#endif

  pType.ToUpper();
  if((pType!="MEDIAN50")&&(pType!="MEDIAN90")&&(pType!="WRC50")&&(pType!="EFFICIENCY")) {
    cout << "Plot Type not available !!!" << endl;
    cout << "Select Plot Type : MEDIAN50/MEDIAN90/WRC50/EFFICIENCY" << endl;
    exit(1);
  }

  // open input tree 
  TFile *file0 = TFile::Open(ifName);
  cout << ifName.Data() << endl;
  if (file0==NULL) {cout << "Error Opening file" << endl;exit(1);}
  TTree* itree = (TTree *) gROOT->FindObject(TREE_NAME);
  if (itree==NULL) {cout << "Error Opening tree" << endl;exit(1);}

  // get detector list
  TList* list = itree->GetUserInfo();
  int nIFO=list->GetSize();
  if (nIFO==0) {cout << "Error : no ifo present in the tree" << endl;exit(1);}
  detector* pDetector[nIFO];
  for (int n=0;n<list->GetSize();n++) {
    pDetector[n] = (detector*)list->At(n);
    detectorParams dParams = pDetector[n]->getDetectorParams();
    pDetector[n]->print();
  }

  // define selection cuts
  char cut[1024];
  sprintf(cut,"run>=0");
#ifdef CC_THR
  sprintf(cut,"%s && netcc[0]>%f",cut,CC_THR);
#endif
#ifdef RHO_THR
  sprintf(cut,"%s && rho[1]>%f",cut,RHO_THR);
#endif
#ifdef ED_THR
  sprintf(cut,"%s && neted[0]/ecor<%f",cut,ED_THR);
#endif
  sprintf(cut,"%s && abs(time[0]-time[%d])<0.1",cut,nIFO);  // WAVE file

  if(type>0) sprintf(cut,"%s && type[1]==%d",cut,type);

  if(SNR_THRESHOLD>0.) {
    sprintf(cut,"%s && sqrt(likelihood)<%f",cut,SNR_THRESHOLD);  // WAVE file
  } else {
    sprintf(cut,"%s && sqrt(likelihood)>%f",cut,-SNR_THRESHOLD);  // WAVE file
  }

  if (pType=="WRC50") {
    char selection[256]="theta[1]:phi[1]";

    char selection_wrc_num[256]="iSNR[0]+oSNR[0]-2*ioSNR[0]";
    for(int i=1;i<nIFO;i++) sprintf(selection_wrc_num,"%s+iSNR[%d]+oSNR[%d]-2*ioSNR[%d]",selection_wrc_num,i,i,i);

    char selection_wrc_den[256]="oSNR[0]";
    for(int i=1;i<nIFO;i++) sprintf(selection_wrc_den,"%s+oSNR[%d]",selection_wrc_den,i);

    sprintf(selection,"%s:(%s)/(%s):sqrt(likelihood)",selection,selection_wrc_num,selection_wrc_den);
    cout << "selection : " << selection << endl;
    itree->Draw(selection,cut,"goff");
  } else {
    itree->Draw("theta[1]:phi[1]:erA[0]:sqrt(likelihood)",cut,"goff");  // MDC file injected coordinates
  }

  int inj_size = itree->GetSelectedRows();
  double* iinj_theta = itree->GetV1(); 
  double* iinj_phi   = itree->GetV2(); 
  double* iierA0     = itree->GetV3(); 
  double* insnr      = itree->GetV4(); 
  double* inj_theta  = new double[inj_size];
  double* inj_phi    = new double[inj_size];
  double* erA0       = new double[inj_size];
  double* nsnr       = new double[inj_size];
  for (int i=0;i<inj_size;i++) {
    inj_theta[i]=iinj_theta[i];
    inj_phi[i]=iinj_phi[i];
    erA0[i]=iierA0[i];
    nsnr[i]=insnr[i];
  }

  if (pType=="WRC50") {
    char selection[256]="";

    char selection_wrc_num[256]="iSNR[0]+oSNR[0]-2*ioSNR[0]";
    for(int i=1;i<nIFO;i++) sprintf(selection_wrc_num,"%s+iSNR[%d]+oSNR[%d]-2*ioSNR[%d]",selection_wrc_num,i,i,i);

    char selection_wrc_den[256]="oSNR[0]";
    for(int i=1;i<nIFO;i++) sprintf(selection_wrc_den,"%s+oSNR[%d]",selection_wrc_den,i);

    sprintf(selection,"(%s)/(%s):erA[1]:erA[5]:erA[9]",selection_wrc_num,selection_wrc_den);
    cout << "selection : " << selection << endl;
    itree->Draw(selection,cut,"goff");
  } else {
    itree->Draw("erA[0]:erA[1]:erA[5]:erA[9]",cut,"goff");  // MDC file injected coordinates
  }

  double* ierA0 = itree->GetV1(); 
  double* ierA1 = itree->GetV2(); 
  double* ierA5 = itree->GetV3(); 
  double* ierA9 = itree->GetV3(); 
  double* erA[10];
  for(int i=0;i<10;i++) erA[i] = new double[inj_size];
  for (int i=0;i<inj_size;i++) {
    erA[0][i]=ierA0[i];
    erA[1][i]=ierA1[i];
    erA[5][i]=ierA5[i];
    erA[9][i]=ierA9[i];
  }

  itree->Draw("theta[0]:phi[0]:inet:type[1]",cut,"goff");  // WAVE file reconstructed coordinates

  cout << cut << endl;

  int size = itree->GetSelectedRows();
  double* rec_theta = itree->GetV1(); 
  double* rec_phi  = itree->GetV2(); 
  double* inet = itree->GetV3(); 
  double* itype = itree->GetV4(); 

  cout << "Selected entries : " << size << endl;

  skymap* SM = new skymap(ANGULAR_RESOLUTION,0,180,0,360);
  int L = SM->size(); 

  double* coverage = new double[L]; for (int l=0;l<L;l++) coverage[l]=0;
  double* mean = new double[L];     for (int l=0;l<L;l++) mean[l]=0;
  double* nri = new double[L];      for (int l=0;l<L;l++) nri[l]=0.;
  double* snr = new double[L];      for (int l=0;l<L;l++) snr[l]=0.;

  int* nEA0 = new int[L];
  for (int l=0;l<L;l++) nEA0[l]=0;
  float** EA0 = new float*[L];
  for (int l=0;l<L;l++) EA0[l]=new float[MAX_ENTRIES_PER_PIXEL];
  for (int l=0;l<L;l++) {
    for (int n=0;n<MAX_ENTRIES_PER_PIXEL;n++) {
      EA0[l][n]=0.;
    }
  }

  for (int i=0;i<size;i++) {  // loop over events

    double inj_ph=inj_phi[i];
    double inj_th=inj_theta[i];

    int inj_ind = SM->getSkyIndex(inj_th,inj_ph);

    if (erA[0][i]<erA[5][i]) coverage[inj_ind]++;
    nri[inj_ind]+=inet[i];
    snr[inj_ind]+=nsnr[i];
    mean[inj_ind]+=erA0[i];
    EA0[inj_ind][nEA0[inj_ind]++]=erA0[i];

    if(nEA0[inj_ind]>MAX_ENTRIES_PER_PIXEL-1) {cout << "Error : entries per pixels to big" << endl;exit(1);}
  }

  int* nMDC = NULL;
  int* nMDC0 = NULL;
  if(pType=="EFFICIENCY") {
    TString mdcFileName = ifName;
    mdcFileName.ReplaceAll("wave_","mdc_");
    cout << "opening file " << mdcFileName.Data();
    TFile *fileMDC = TFile::Open(mdcFileName.Data());
    if (fileMDC==NULL) {cout << "Error Opening mdc file" << endl;exit(1);}
    TTree* itreeMDC = (TTree *) gROOT->FindObject("mdc");
    if (itreeMDC==NULL) {cout << "Error Opening mdc tree" << endl;exit(1);}
    char mdc_cut[256]="run>0";
    if(type>0) sprintf(mdc_cut,"%s && type[0]==%d",mdc_cut,type);
    cout << "mdc_cut : " << mdc_cut << endl;
    itreeMDC->Draw("theta[0]:phi[0]",mdc_cut,"goff");
    int mdc_size = itreeMDC->GetSelectedRows();
    double* mdc_theta = itreeMDC->GetV1();
    double* mdc_phi  = itreeMDC->GetV2();

    nMDC = new int[L];
    nMDC0 = new int[L];

    for (int l=0;l<L;l++) nMDC0[l]=0;
    for (int l=0;l<L;l++) nMDC[l]=0;
    for (int i=0;i<mdc_size;i++) {  // loop over mdc
      int mdc_ind = SM->getSkyIndex(mdc_theta[i],mdc_phi[i]);
      nMDC0[mdc_ind]++;
    }
  }

  int* nEA = new int[L];
  for (int l=0;l<L;l++) nEA[l]=0;
  double** EA = new double*[L];
  for (int l=0;l<L;l++) EA[l]=new double[MAX_ENTRIES_PER_PIXEL];
  for (int l=0;l<L;l++) {
    for (int n=0;n<MAX_ENTRIES_PER_PIXEL;n++) {
      EA[l][n]=0.;
    }
  }

  double* snr_smooth = new double[L];
  for (int l=0;l<L;l++) snr_smooth[l]=0.;

  // skymap smoothing

  for(int l=0; l<L; l++){
    double th = SM->getTheta(l);
    double ph = SM->getPhi(l);
    double st = SM->getThetaStep(l);

    for(int in=-NPIX_SMOOTH; in<=NPIX_SMOOTH; in++) {
      double ith = th+in*st; if(ith<0) ith+=180; if(ith>180) ith-=180;  // FIX
      double sp = SM->getPhiStep(SM->getSkyIndex(ith,ph));  // FIX
      //double sp = SM->getPhiStep(SM->getSkyIndex(th+in*st,ph));
      for(int im=-NPIX_SMOOTH; im<=NPIX_SMOOTH; im++) {
        double iph = ph+im*sp; if(iph<0) iph+=360; if(iph>360) iph-=360;  // FIX
        int i = SM->getSkyIndex(ith,iph);       // neighbour sky index  // FIX
        //int i = SM->getSkyIndex(th+in*st,ph+im*sp);       // neighbour sky index
        if((nEA0[l]+nEA0[i])>MAX_ENTRIES_PER_PIXEL-1) {cout << "Error : entries per pixels to big" << endl;exit(1);}
        for (int n=0;n<nEA0[i];n++) EA[l][nEA[l]++]=EA0[i][n];
        snr_smooth[l]+=snr[i];
        if(pType=="EFFICIENCY") nMDC[l]+=nMDC0[i];
      }
    }
  }

  int nEA_max=0;
  for (int l=0;l<L;l++) if(nEA[l]>nEA_max) nEA_max=nEA[l];
  for (int l=0;l<L;l++) if(nEA[l]>0) coverage[l]/=(double)nEA[l];
  for (int l=0;l<L;l++) if(nEA[l]>0) nri[l]/=(double)nEA[l];
  for (int l=0;l<L;l++) if(nEA[l]>0) mean[l]/=(double)nEA[l];
  for (int l=0;l<L;l++) if(nEA[l]>0) snr_smooth[l]/=(double)nEA[l];

  double* tmp = new double[nEA_max];
  int* index = new int[nEA_max];
  for (int l=0;l<L;l++) if(nEA[l]>0) {
    for (int n=0;n<nEA[l];n++) tmp[n]=EA[l][n];
    TMath::Sort(nEA[l],tmp,index,false);
    for (int n=0;n<nEA[l];n++) EA[l][n]=tmp[index[n]];
  }  
  delete [] tmp;
  delete [] index;

  double* EA50 = new double[L];
  double* EA90 = new double[L];
  for (int l=0;l<L;l++) if(nEA[l]>MIN_TRIGGERS_PER_PIXEL) {
    EA50[l]=EA[l][nEA[l]/2];
    EA90[l]=EA[l][(nEA[l]*9)/10];
  } else {EA50[l]=0.; EA90[l]=0.;}

  // get label from file title
  TObjArray* token = TString(ifName).Tokenize(TString("/"));
  TObjString* otoken = (TObjString*)token->At(token->GetEntries()-1);
  TString ifile_label = otoken->GetString();
  ifile_label.ReplaceAll(".root","");
  cout << "ifile_label " << ifile_label.Data() << endl;

  char title_label[256]="";
  sprintf(title_label,"%s : ",pType.Data());
  char cut_label[256]="";
  sprintf(cut_label,"EC");
#ifdef CC_THR
  sprintf(cut_label,"%s_cc_%.2f",cut_label,CC_THR);
#endif
#ifdef RHO_THR
  sprintf(cut_label,"%s_rho_%.2f",cut_label,RHO_THR);
#endif
#ifdef ED_THR
  sprintf(cut_label,"%s_ed_%.2f",cut_label,ED_THR);
#endif
  sprintf(cut_label,"%s_R%dx%d",cut_label,int(ANGULAR_RESOLUTION),int(ANGULAR_RESOLUTION));
  cout << "ANGULAR_RESOLUTION : " << ANGULAR_RESOLUTION << endl;
  cout << "cut_label : " << cut_label << endl;

  gnetwork* gNET = new gnetwork;
  for (int n=0;n<list->GetSize();n++) gNET->add(pDetector[n]);

  gNET->setSkyMaps(ANGULAR_RESOLUTION,0,180,0,360);
  gNET->setAntenna();
  gNET->setDelay(const_cast<char*>(pDetector[0]->Name));

  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

#ifdef DISPLAY_WORLD_MAP
  gSM->SetWorldMap();
#endif

  if (pType=="MEDIAN50") for (int l=0;l<L;l++) gSM->set(l,EA50[l]);
  if (pType=="MEDIAN90") for (int l=0;l<L;l++) gSM->set(l,EA90[l]);
  if (pType=="WRC50")    for (int l=0;l<L;l++) gSM->set(l,EA50[l]);

  if (pType==("EFFICIENCY")) {
    for (int l=0;l<L;l++) if (nMDC[l]>0) gSM->set(l,double(nEA[l])/nMDC[l]); else gSM->set(l,0.0);
  }

  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  if (pType=="MEDIAN50")   h2->GetZaxis()->SetRangeUser(0.5,50);
  if (pType=="MEDIAN90")   h2->GetZaxis()->SetRangeUser(0.5,50);
  if (pType=="EFFICIENCY") h2->GetZaxis()->SetRangeUser(0.,1);
  if (pType=="WRC50")      h2->GetZaxis()->SetRangeUser(0.05,0.5);

  // set title
  TString sTITLE=Title;
  if(sTITLE=="") {
    if (pType=="MEDIAN50")   sTITLE = "Median Error Region 50%";
    if (pType=="MEDIAN90")   sTITLE = "Median Error Region 90%";
    if (pType=="EFFICIENCY") sTITLE = "Efficiency";
    if (pType=="WRC50")      sTITLE = "Median Normalized Residual Energy 50%";
  }
  if(Tag!="") sTITLE=sTITLE+" "+Tag; 
  gSM->SetTitle(sTITLE);

  if (pType=="MEDIAN50") gSM->SetLogz(true);
  if (pType=="MEDIAN90") gSM->SetLogz(true);
  if (pType=="WRC50")    gSM->SetLogz(true);
  gSM->Draw(0);

#ifdef DRAW_SITES
  char ifostr[64]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,pDetector[n]->Name);
  }
  cout << "oNetwork : " << ifostr << endl;

  // Set Network (only to print sites)
  gNET->DrawSitesShortLabel(kBlack);
  gNET->DrawSites(kBlack,2.0);
  gNET->DrawSitesArms(1000000,kWhite,2.0);
#endif

  // compute median50% average over the sky
  double MeanEA50=0;
  int nMeanEA50=0;
  for (int l=0;l<L;l++) if(nEA[l]>0) {MeanEA50+=EA50[l];nMeanEA50++;}
  MeanEA50/=nMeanEA50;
  cout << "MeanEA50 : " << MeanEA50 << endl;

  char skymap_file_name[256]="";
#ifdef WRITE_SKYMAP
  if(type>0) {
    sprintf(skymap_file_name,"SKYMAP_%s_%s_RES%d%s_%d.root",
            pType.Data(),ifile_label.Data(),int(ANGULAR_RESOLUTION*10.),DATA_LABEL,type);
  } else {
    sprintf(skymap_file_name,"SKYMAP_%s_%s_RES%d%s.root",
            pType.Data(),ifile_label.Data(),int(ANGULAR_RESOLUTION*10.),DATA_LABEL);
  }
  cout << skymap_file_name << endl;
/*
  gskymap oSM(ANGULAR_RESOLUTION,0,180,0,360);
  if (pType=="MEDIAN50") for (int l=0;l<L;l++) oSM.set(l,EA50[l]);
  if (pType=="MEDIAN90") for (int l=0;l<L;l++) oSM.set(l,EA90[l]);
  oSM.DumpObject(skymap_file_name);
*/
  gNET->DumpObject(skymap_file_name);
#endif

#ifdef WRITE_PLOT
  sprintf(skymap_file_name,"%s_%s_RES%d%s.png",
          pType.Data(),ifile_label.Data(),int(ANGULAR_RESOLUTION*10.),DATA_LABEL);
  cout << skymap_file_name << endl;
  gSM->Print(skymap_file_name);
  exit(0);
#endif

}

