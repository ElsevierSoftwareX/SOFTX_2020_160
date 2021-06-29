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


// this macro compute the FAD statistic 
// it is used in post production 
// can be run in standalone or it is called from the sim report

#include <vector>

#define FAR_FILE_NAME "rate_threshold_veto.txt"
#define LIV_FILE_NAME "live.txt"

// -----------------------------------------
// plot type : ptype
// -----------------------------------------
// "FARvsRHO"
// "REFFvsRHO"
// "REFFvsIRHO"
// "VISvsRHO"
// "PRODvsFAD"
// "EVTvsRHO"
// "FADvsRHO"
// "FAPvsRHO"
// "RECvsINJ"
// "OBSTIMEvsFAR"
// -----------------------------------------

#define DENSITY_FORMULA
#define APPLY_ENORM

void MakePlot(TString ptype, bool pAll); 
TGraphErrors* Plot(int mtype, TString ptitle, TString xtitle, TString ytitle, bool save=true);
void RECvsINJ(int mtype, TString ptitle, TString xtitle, TString ytitle);
void MakeHtml(TString ptitle);
void MakeHtmlIndex(TString* ptype, int psize);
void AddHtmlHeader(TString odir); 
void MakeMultiPlotFAD(TString mdc_inj_file);

// --------------------------------------------------------
// Global variables
// --------------------------------------------------------
TString   gODIR   = pp_dir;		// output dir
TString   gPTYPE  = "VISvsRHO";		// plot type
int       gMTYPE  = -1;                 // mdc type
bool      gFIT    = false;		// enable fit
TCanvas*  gCANVAS = NULL;		// canvas object
CWB::mdc* gMDC    = NULL;		// MDC object
TTree*    gTRWAVE = NULL;		// wave tree
TTree*    gTRMDC  = NULL;		// mdc tree
bool      gINSPIRAL = false;		// inspiral/burst mdc
double    gOBSTIME= 0.0;		// observation livetime
double    gBKGTIME= 0.0;		// background  livetime
double    gRHOMIN = 5.;			// minimum rho value selected for plots
int       gNBINS  = 100;                // number of bins used in mdc histo 
int       gNZBINS = -1; 		// le_0 -> standard FAD statistic
                                        // ge_0 -> modified FAD statistic  

vector<double> RHO;			// rho
vector<double> eRHO;			// RHO error
vector<double> VIS;			// visible volume`
vector<double> eVIS;			// VIS error
vector<double> FAR;			// FAR
vector<double> eFAR;			// FAR error
vector<double> dFAR;			// differendial FAR
vector<double> normHRSS;		// normalization hrss

// --bkgrep  : directory of the background report (used to read rate_threshold_veto.txt and live.txt)
// --hrss    : if it is a number then it is used as normalization constant for all MDC types (def=0)
//             if =0 then hrss rescale is not applied
//             if it is a file then it is the list of hrss used for each MDC type 
//             (format : for each line -> hrss)
// --gfit    : true/false -> enable/disable fit in the output plots (default=false)
// --rhomin  : minimum rho value selected for plots (default = 5)
// --units   : K -> Kpc,Kyr : M -> Mpc,Myr  (def=M)
// --distr   : formula/mdc -> radial distribution is computed from formula or from mdc injections (def=MDC)
// --nbins   : number of bins used in hist to computed the radial distribution from the mdc injections
// --header  : if true -> add cwb header to fad html file (def=false)
// --multi   : if true -> FAD multi plot for each mdc set are created and substituted in the sim report 
//             page to the eff_freq plots (def=false)
// --title   : title of the html page (default=FAD) : spaces must be filled with *
// --subtitle: subtitle of the html page (default="") : spaces must be filled with *


TString pp_fad_bkgrep="";
TString pp_fad_hrss="";
TString pp_fad_rhomin="";
TString pp_fad_gfit="";
TString pp_fad_units="";
TString pp_fad_distr="";
TString pp_fad_nbins="";
TString pp_fad_nzbins="";
TString pp_fad_title="";
TString pp_fad_subtitle="";

#define nPLOT 10

TString ptype[nPLOT] = {
		    "FARvsRHO",
		    "RECvsINJ",
		    "VISvsRHO",
		    "REFFvsRHO",
		    "REFFvsIRHO",
		    "FADvsRHO",
		    "PRODvsFAD",
		    "EVTvsRHO",
		    "FAPvsRHO",
		    "OBSTIMEvsFAR"
                   };


void cwb_mkfad(TString odir="fad", int nzbins=-1, double obstime=0, bool bexit=true) {
//
// odir    : output fad directory (default="fad")
// nzbins  : (default is -1)
//           - NOTE : if it is defined in pp_fad the input value is overwritten
//           - if nzbin=0 the FAD statistic is computed with classical FAD 
//           - if nzbin>0 the FAD statistic is computed until there are nzbins
//             consecutive bins with zero events inside
//           - if nzbin<0 the FAD statistic is computed with classical FAD and min-hold 
// obstime : zero lag observation time, if 0 then it is read from live.txt 
// bexit   : if true (def) then the macro exit at the end of the execution
//
// NOTE    : pp_fad must be defined in user_pparameters.C file or in the input 
//           fad config file (fad.in) provided to the line command cwb_mkfad   
//           Ex: cwb_mkfad M1 fad.in fad  
// 

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  // ------------------------------------------------
  // cwb_mkfad works only with simulation=4 mode
  // ------------------------------------------------
/*
  if(simulation!=4) {
    cout << "cwb_mkfad.C - Error : cwb_mkfad works only with simulation=4 mode" << endl;
    exit(1);
  }
*/

  // ------------------------------------------------
  // read pp_fad
  // ------------------------------------------------
  if(pp_fad=="") { 	// nothing to be done
    cout<<endl<<"cwb_mkfad.C : Warning : pp_fad is not defined, macro is not executed !!!"<<endl<<endl;
    exit(0);		
  }

  // get --bkgrep
  pp_fad_bkgrep = CWB::Toolbox::getParameter(pp_fad,"--bkgrep");
  if(pp_fad_bkgrep=="") {
    cout<<"cwb_mkfad.C : Error : pp_fad --bkgrep not defined"<<endl;exit(1);}

  TB.checkFile(pp_fad_bkgrep);

  // get --hrss
  pp_fad_hrss = CWB::Toolbox::getParameter(pp_fad,"--hrss");
  if(pp_fad_hrss=="") {pp_fad_hrss="0";
    cout<<"cwb_mkfad.C : Info : pp_fad --hrss not defined : set hrss=1"<<endl;}

  if(!pp_fad_hrss.IsFloat()) { 
    // hrss normalization values are define in a file
    TB.checkFile(pp_fad_hrss);
  }

  // get --rhomin
  pp_fad_rhomin = CWB::Toolbox::getParameter(pp_fad,"--rhomin");
  if(pp_fad_rhomin!="") {
    if(pp_fad_rhomin.IsFloat()) { 
      gRHOMIN = pp_fad_rhomin.Atof();
    } else {
      cout<<"cwb_mkfad.C : Error : pp_fad --rhomin must be a number"<<endl;
      exit(1);
    }
  }

  // get --gfit
  pp_fad_gfit = CWB::Toolbox::getParameter(pp_fad,"--gfit");
  if(pp_fad_gfit!="") {
    if((pp_fad_gfit!="true")&&(pp_fad_gfit!="false")) {
      cout<<"cwb_mkfad.C : Error : pp_fad --gfit bad value -> must be true/false"<<endl;
      exit(1);
    } else {
      if(pp_fad_gfit=="true")  gFIT=true;
      if(pp_fad_gfit=="false") gFIT=false;
    }
  }

  // get --units
  pp_fad_units = CWB::Toolbox::getParameter(pp_fad,"--units");
  if(pp_fad_units=="") pp_fad_units="M";
  pp_fad_units.ToUpper();
  if((pp_fad_units!="K")&&(pp_fad_units!="M")) {
    cout<<"cwb_mkfad.C : Error : pp_fad --units bad value -> must be K/M"<<endl;
    exit(1);
  } 

  // get --distr
  pp_fad_distr = CWB::Toolbox::getParameter(pp_fad,"--distr");
  if(pp_fad_distr=="") pp_fad_distr="mdc";
  pp_fad_distr.ToUpper();
  if((pp_fad_distr!="MDC")&&(pp_fad_distr!="FORMULA")) {
    cout<<"cwb_mkfad.C : Error : pp_fad --distr bad value -> must be MDC/FORMULA"<<endl;
    exit(1);
  } 

  // get --nzbins (if it is defined in pp_fad the input value is overwritten)
  gNZBINS = nzbins; 
  pp_fad_nzbins = CWB::Toolbox::getParameter(pp_fad,"--nzbins");
  if(pp_fad_nzbins!="") {
    TString _pp_fad_nzbins = pp_fad_nzbins;
    if((_pp_fad_nzbins[0]=='+')||(_pp_fad_nzbins[0]=='-')) _pp_fad_nzbins.Remove(0,1);
    if(_pp_fad_nzbins.IsDigit()) { 
      gNZBINS = pp_fad_nzbins.Atoi();
    } else {
      cout<<"cwb_mkfad.C : Error : pp_fad --nzbins must be an integer number"<<endl;
      exit(1);
    }
  } 

  // get --nbins
  pp_fad_nbins = CWB::Toolbox::getParameter(pp_fad,"--nbins");
  if(pp_fad_nbins!="") {
    if(pp_fad_nbins.IsDigit()) { 
      gNBINS = pp_fad_nbins.Atoi();
    } else {
      cout<<"cwb_mkfad.C : Error : pp_fad --nbins must be an integer number"<<endl;
      exit(1);
    }
  } 

  // get --header
  TString pp_fad_header = CWB::Toolbox::getParameter(pp_fad,"--header");
  if(pp_fad_header=="") pp_fad_header="false";
  pp_fad_header.ToUpper();

  // get --multi
  TString pp_fad_multi = CWB::Toolbox::getParameter(pp_fad,"--multi");
  if(pp_fad_multi=="") pp_fad_multi="false";
  pp_fad_multi.ToUpper();

  // get --title
  TString pp_fad_title = CWB::Toolbox::getParameter(pp_fad,"--title");
  pp_fad_title.ReplaceAll("*"," ");
  if(pp_fad_title=="") pp_fad_title="FAD";

  // get --subtitle
  TString pp_fad_subtitle = CWB::Toolbox::getParameter(pp_fad,"--subtitle");
  pp_fad_subtitle.ReplaceAll("*"," ");

  // ------------------------------------------------
  // create canvas
  // ------------------------------------------------
  gCANVAS = new TCanvas("canvas", "canvas",16,30,825,546);
  gCANVAS->Range(-19.4801,-9.25,-17.4775,83.25);
  gCANVAS->SetBorderSize(2);
  gCANVAS->SetFrameFillColor(0);
  gCANVAS->SetGridx();
  gCANVAS->SetGridy();

  gStyle->SetOptFit(kTRUE);

  // ------------------------------------------------
  // define output dir
  // ------------------------------------------------
  gODIR = TString(pp_dir)+TString("/")+odir;

  // ------------------------------------------------
  // export to CINT net,cfg
  // exec config plugin
  // if Burst MDC then extract mdcHRSS which is defined in config plugin
  // ------------------------------------------------
  char cmd[128];
  //sprintf(cmd,"network* net = (network*)%p;",net);
  sprintf(cmd,"network* net = new network;");
  gROOT->ProcessLine(cmd);
  CWB::config* cfg = new CWB::config;
  cfg->Import();
  sprintf(cmd,"CWB::config* cfg = (CWB::config*)%p;",cfg);
  gROOT->ProcessLine(cmd);
  sprintf(cmd,"int gIFACTOR=1;");
  gROOT->ProcessLine(cmd);
  configPlugin.Exec();
  //MDC.Print();
  gMDC = &MDC;

  gINSPIRAL = gMDC->GetInspiral()!="" ? true : false;         // Inspiral/Burst MDC

  if(!gINSPIRAL) {	// if burst then mdcHRSS must be declaren in the config plugin
    // verify if mdcHRSS exist
    TGlobal* global = (TGlobal*)gROOT->GetGlobal("mdcHRSS",true);
    if(global!=NULL) {
      cout << global->GetTypeName() << endl;
      if(TString(global->GetTypeName())!="vector<double>") {
        cout << "cwb_mkfad.C - Error : 'vector<double> mdcHRSS' do not exist in CINT" << endl;
        gSystem->Exit(1);
      }
    }
  }

  // ------------------------------------------------
  // open wave/mdc trees
  // ------------------------------------------------
  // open wave file
  TFile* fwave = TFile::Open(sim_file_name);
  gTRWAVE = (TTree*) gROOT->FindObject("waveburst");
  // open mdc file
  TFile *fmdc = TFile::Open(mdc_file_name);
  gTRMDC = (TTree*) gROOT->FindObject("mdc");

  // ------------------------------------------------
  // get observation time
  // ------------------------------------------------
  char liv_file_path[1024];
  sprintf(liv_file_path,"%s/%s",pp_fad_bkgrep.Data(),LIV_FILE_NAME);
  int countlag=0;
  gBKGTIME=GetLiveTime(liv_file_path,-1,-1,countlag);   // get full  livetime
  gOBSTIME=GetLiveTime(liv_file_path,0,0,countlag);   	// get zero lag livetime
  gBKGTIME-=gOBSTIME;		      			// background livetime
  if(obstime>0)  gOBSTIME=obstime;			// user provided zero lag livetime
  cout << "Zero Lag Observation Time : " << gOBSTIME << " sec" << endl;
  if(gOBSTIME==0) {
    cout << "cwb_mkfad.C - Error : Observation Time is zero" << endl;
    cout << "                      probably your set do not includes the zero lag" << endl << endl;
    exit(1);
  }
  cout << "Background Observation Time : " << gBKGTIME << " sec" << endl;

  // ------------------------------------------------
  // read cumulative rho far from rate_threshold_veto.txt file
  // ------------------------------------------------

  char far_file_path[1024];
  sprintf(far_file_path,"%s/%s",pp_fad_bkgrep.Data(),FAR_FILE_NAME);

  ifstream infar;
  infar.open(far_file_path, ios::in);
  if (!infar.good()) {cout << "Error Opening File : " << far_file_path << endl;exit(1);}
  while (1) {
    float irho,ifar;
    infar >> irho >> ifar;
    if (!infar.good()) break;
    if(irho<gRHOMIN) continue;
    if(ifar==0) break;
    RHO.push_back(irho);
    FAR.push_back(ifar);
    eFAR.push_back(sqrt(gBKGTIME*ifar)/gBKGTIME);
    dFAR.push_back(ifar);
  }
  infar.close();
  // add one more value with FAR=0
  float drho=RHO[RHO.size()-1]-RHO[RHO.size()-2];
  RHO.push_back(RHO[RHO.size()-1]+drho);
  FAR.push_back(0);
  eFAR.push_back(0);
  dFAR.push_back(0);

  // add eRHO
  eRHO.resize(RHO.size());
  for(int i=0;i<eRHO.size();i++) eRHO[i]=0;

  // compute the differential FAR
  for(int i=0;i<RHO.size()-1;i++) {
    dFAR[i]=dFAR[i]-dFAR[i+1];
  }

  for(int i=0;i<RHO.size()-1;i++) {
    //cout << i << " " <<RHO[i] << " " << dFAR[i] << endl;
  }

  // ------------------------------------------------
  // read normalization hrss
  // ------------------------------------------------
  int N = gMDC->wfList.size(); // get number of MDC types

  if(pp_fad_hrss.IsFloat()) { // apply the same valye to all mtypes

    for(int i=0;i<N;i++) normHRSS.push_back(pp_fad_hrss.Atof());

  } else { // hrss normalization values are define in a file

    ifstream in_hrss_norm;
    in_hrss_norm.open(pp_fad_hrss.Data());
    int count=0;
    double hrss_temp;
    while(1) {
      in_hrss_norm >> hrss_temp;
      if (!in_hrss_norm.good()) break;
      //cout << hrss_temp << endl;
      normHRSS.push_back(hrss_temp);
      count++;
      if(count>N) {
        cout << "cwb_mkfad.C - Errors : too many waveforms on " << pp_fad_hrss << " file " << endl; 
        exit(1);
      }
    }
    in_hrss_norm.close();
  }

  // ------------------------------------------------
  // create output dir
  // ------------------------------------------------
  TB.mkDir(gODIR,!pp_batch);

  // ------------------------------------------------
  // create multiPlot - used in sim report page
  // ------------------------------------------------
  if(pp_fad_multi=="TRUE") MakeMultiPlotFAD(mdc_inj_file);

  // ------------------------------------------------
  // create plots
  // ------------------------------------------------
  gMTYPE  = -1;   // reset initial value of mdc type
  gCANVAS->cd();  // must be done because MakeMultiPlotFAD set a different canvas
  for(int i=0;i<nPLOT;i++) MakePlot(ptype[i], false);

  // ------------------------------------------------
  // make fad.html
  // ------------------------------------------------
  MakeHtmlIndex(ptype,nPLOT);
  if(pp_fad_header=="TRUE") AddHtmlHeader(gODIR); 

  if(bexit) gSystem->Exit(0); else return;
}

void MakePlot(TString ptype, bool pAll) {

  gPTYPE = ptype;

  char rho_type[32];sprintf(rho_type,"rho[%d]",pp_irho);

  TString ptitle;
  TString xtitle;
  TString ytitle;

  if(gPTYPE=="REFFvsRHO")  {
     gCANVAS->SetLogx(false);
     gCANVAS->SetLogy(false);
     ptitle=TString("Effective Radius vs ")+rho_type;
     gStyle->SetStatY(0.3);
     xtitle=rho_type;
     ytitle="Reff (Kpc)";

  } else if(gPTYPE=="REFFvsIRHO") {
     gCANVAS->SetLogx(false);
     gCANVAS->SetLogy(true);
     ptitle=TString("Effective Radius vs inverse ")+rho_type;
     xtitle=TString(rho_type)+"^{-1}";
     ytitle="Reff (Kpc)";

  } else if(gPTYPE=="VISvsRHO") {
     gCANVAS->SetLogx(false);
     gCANVAS->SetLogy(true);
     ptitle=TString("Visible Volume vs ")+rho_type;
     xtitle=rho_type;
     ytitle="Vvis (Kpc^{3})";

  } else if(gPTYPE=="PRODvsFAD") {
     gCANVAS->SetLogx(true);
     gCANVAS->SetLogy(false);
     ptitle=TString("Productivity vs False Alarm Density");
     xtitle="FAD (Kpc^{-3} Kyr^{-1})";
     ytitle="Productivity (Kpc^{3} Kyr)";

  } else if(gPTYPE=="FADvsRHO") {
     gCANVAS->SetLogx(false);
     gCANVAS->SetLogy(true);
     ptitle=TString("False Alarm Density vs ")+rho_type;
     xtitle=rho_type;
     ytitle="FAD (Kpc^{-3} Kyr^{-1})";

  } else if(gPTYPE=="FARvsRHO") {
     gCANVAS->SetLogx(false);
     gCANVAS->SetLogy(true);
     ptitle=TString("False Alarm Rate vs ")+rho_type;
     xtitle=rho_type;
     ytitle="False Alarm Rate (sec^{-1})";

  } else if(gPTYPE=="FAPvsRHO") {
     gCANVAS->SetLogx(false);
     gCANVAS->SetLogy(true);
     ptitle=TString("False Alarm Probability vs ")+rho_type;
     xtitle=rho_type;
     ytitle="False Alarm Probability";

  } else if(gPTYPE=="EVTvsRHO") {
     gCANVAS->SetLogx(false);
     gCANVAS->SetLogy(true);
     ptitle=TString("Expected background events per Observation time vs ")+rho_type;
     xtitle=rho_type;
     ytitle="Expected Events";

  } else if(gPTYPE=="OBSTIMEvsFAR") {
     gCANVAS->SetLogx(true);
     gCANVAS->SetLogy(false);
     ptitle=TString("Observation Time vs False Alarm Rate");
     xtitle="FAR (sec^{-1})";
     ytitle="Observation Time (sec)";

  } else if(gPTYPE=="RECvsINJ") {
     gCANVAS->SetLogx(false);
     gCANVAS->SetLogy(true);
     ptitle="Reconstructed events vs Injected events";
     gStyle->SetOptStat(0);
     xtitle = "distance (Kpc)";
     ytitle = "counts";
  }

  // get number of MDC types
  int N = gMDC->wfList.size();

  if(gPTYPE=="RECvsINJ") {

    if(pAll) for(int i=0;i<=N;i++) RECvsINJ(i, ptitle, xtitle, ytitle);
    else     RECvsINJ(0, ptitle, xtitle, ytitle);

  } else { 	

    if(gPTYPE=="EVTvsRHO") {
      Plot(0, ptitle, xtitle, ytitle);
    } else {
  
#ifdef APPLY_ENORM
      if(pAll) for(int i=0;i<=N;i++) Plot(i, ptitle, xtitle, ytitle);
      else     Plot(0, ptitle, xtitle, ytitle);
#else
      for(int i=1;i<=N;i++) Plot(i, ptitle, xtitle, ytitle);
      //Plot(1, ptitle);
#endif

    } 
  }

  return;
}

void getVisibleVolume(int mtype) {

  char sel[1024];
  char cut[1024];
  TF1 *vdens = NULL;
//  TF1 *hdist = NULL;

  // INJECTED
  if(mtype==0) {
    sprintf(cut,"1==1");
  } else {
    //sprintf(cut,"type[0]==%d",mtype);
    sprintf(cut,"type==%d",mtype);
  }

  double min_dist = gTRMDC->GetMinimum("distance");
  double max_dist = gTRMDC->GetMaximum("distance");
#ifdef APPLY_ENORM
  double enorm_max=0;
  for(int i=0;i<normHRSS.size();i++) {
    double enorm=normHRSS[i]/mdcHRSS[i];
    if(enorm>enorm_max) enorm_max=enorm;
  }
  min_dist *= enorm_max;
  max_dist *= enorm_max;
#endif

  TH1F* hdist = new TH1F("hdist","hdist",gNBINS,1000*min_dist,1000*max_dist); // Mpc -> Kpc

  TTreeFormula mcut("mcut", cut, gTRMDC);
  gTRMDC->SetBranchStatus("*",false);
  gTRMDC->SetBranchStatus("distance",true);
  gTRMDC->SetBranchStatus("type",true);
  int type;
  float distance;
  gTRMDC->SetBranchAddress("distance",&distance);
  gTRMDC->SetBranchAddress("type",&type);
  int msize = gTRMDC->GetEntries();
  int nmdc = 0;
  for(int i=0;i<msize;i++) {
    gTRMDC->GetEntry(i);
    if(mcut.EvalInstance()==0) continue;
    distance*=1000; 
#ifdef APPLY_ENORM
    // compute Normalization to hrss @ 10 Kpc
    if(normHRSS.size()&&(normHRSS[type-1]>0)) {
      distance *= normHRSS[type-1]/mdcHRSS[type-1];
    }
#endif
    hdist->Fill(distance);
    nmdc++;
  }
  cout << "nmdc : " << nmdc << endl; 
  gTRMDC->SetBranchStatus("*",true);

/*
  double min_dist = gTRMDC->GetMinimum("distance");
  double max_dist = gTRMDC->GetMaximum("distance");

  // get total number of inject MDC with TYPE=type
  sprintf(sel,"1000*distance>>hdist1(%d,%f,%f)",gNBINS,1000*min_dist,1000*max_dist);	// Mpc -> Kpc
  if(mtype==0) {
    sprintf(cut,"");
  } else {
    sprintf(cut,"type[0]==%d",mtype);
  } 
  gTRMDC->Draw(sel,cut,"goff");  
  int nmdc = gTRMDC->GetSelectedRows();
  cout << "nmdc : " << nmdc << endl;
*/

  // get Sky Distribution parameters 
  TString sky_dist_formula = gMDC->GetSkyFile();		// get distribution formula
  vector<mdcpar> sky_parms = gMDC->GetSkyParms();
  bool error=false;
  int entries = gMDC->GetPar("entries",sky_parms,error);
  double min_dist = gMDC->GetPar("rho_min",sky_parms,error);	// get rho min
  if(error) min_dist=0.;
  double dist_max = gMDC->GetPar("rho_max",sky_parms,error);	// get rho max
  if(error) dist_max=0.;

  bool formula = ((sky_dist_formula!="")&&(dist_max>0)) ? true : false;

  formula&=(pp_fad_distr=="FORMULA");

  if(formula) { // radial distribution is extract from the config plugin

    cout << "sky_dist_formula : " << sky_dist_formula << endl;
    cout << "min_dist : " << min_dist << " dist_max : " << dist_max << endl;
  
    // compute the formula normalization
    // the integral between min_dist and dist_max must be the total number 
    // of the injected mdc = nmdc
    TF1 *fdist = new TF1("fdist",sky_dist_formula,min_dist,dist_max);
    double norm=fdist->Integral(min_dist,dist_max);
    norm=nmdc/norm;
    cout << norm << endl;
  
    // visible volume density formula
    char vdens_expression[256];
    sprintf(vdens_expression,"(4*TMath::Pi()*x*x)/(%f*(%s))",norm,sky_dist_formula.Data());
    cout << "vdens_expression : " << vdens_expression << endl;
    vdens = new TF1("vdens",vdens_expression,min_dist,dist_max);

  } else { 	// radial distribution is extract from the mdc injection root file

//hdist = (TH1F*)gDirectory->Get("hdist1");
    if (hdist!=NULL) {
      // compute radial distribution density histogram
      for(int i=1;i<=hdist->GetNbinsX();i++) {
        double value = hdist->GetBinContent(i);
        value /= hdist->GetBinWidth(i);	
        hdist->SetBinContent(i,value);
      }
      //hdist->Draw("HIST");
      //gCANVAS->Print("hdist.png");
    } else {
      cout << "cwb_mkfad.C - Error : hdist is NULL" << endl; 
      gSystem->Exit(1);
    }
  }  

  int nRHO = RHO.size();

  // compute visible volume vs rho
  VIS.resize(nRHO);
  eVIS.resize(nRHO);
  VIS[nRHO-1] = 0;
  eVIS[nRHO-1] = 0;
  for(int n=nRHO-1;n>=0;n--) {

    double rho_min = RHO[n];
    double rho_max = n==nRHO-1 ? 1.e80 : RHO[n+1];    // when n=nRHO-1 the range is [RHO[n],1e80]

    sprintf(sel,"range[1]:type[1]");
    if(mtype==0) {
      sprintf(cut,"abs(time[0]-time[%d])<%f && netcc[%d]>%f && rho[%d]>%f && rho[%d]<=%f",
              nIFO,T_win,pp_inetcc,T_cor,pp_irho,rho_min,pp_irho,rho_max);
    } else {
      sprintf(cut,"type[1]==%d && abs(time[0]-time[%d])<%f && netcc[%d]>%f && rho[%d]>%f && rho[%d]<=%f",
              mtype,nIFO,T_win,pp_inetcc,T_cor,pp_irho,rho_min,pp_irho,rho_max);
    }
    if(T_vED>0) sprintf(cut,"%s && neted[0]/ecor<%f",cut,T_vED);
    if(T_pen>0) sprintf(cut,"%s && penalty>%f",cut,T_pen);
    gTRWAVE->Draw(sel,cut,"goff");  
    int nwave = gTRWAVE->GetSelectedRows();
    //cout << "rho_min " << rho_min << "\trho_max " << rho_max << "\tnwave : " << nwave << endl; 
    double* range = gTRWAVE->GetV1();
    double* itype = gTRWAVE->GetV2();

    double Pi = TMath::Pi();

    for(int i=0;i<nwave;i++) {
      double dist=1000*range[i];	// Mpc -> Kpc
#ifdef APPLY_ENORM
      // compute Normalization to hrss @ 10 Kpc
      if(normHRSS.size()&&(normHRSS[itype[i]-1]>0)) {
        dist *= normHRSS[itype[i]-1]/mdcHRSS[itype[i]-1];
      }
#endif

      double ivdens;
      if(formula) {
        ivdens = vdens->Eval(dist);
      } else {
        ivdens = (4*TMath::Pi()*dist*dist)/hdist->Interpolate(dist);
      }
      VIS[n] +=ivdens;
      eVIS[n]+=pow(ivdens,2);
    }
    eVIS[n]=sqrt(eVIS[n]);

    if(n>0) VIS[n-1]=VIS[n];

    if(pp_fad_units=="M") {VIS[n]*=1.e-9; eVIS[n]*=1.e-9;}	// Kpc^3 -> Mpc^3

    //cout << n << " " << RHO[n] << " " << VIS[n] << endl;
    //cout << n << " " << eRHO[n] << " " << eVIS[n] << endl;
  }

  if(hdist) delete hdist;

  gMTYPE = mtype;	// update global value
  return;
}

TGraphErrors* Plot(int mtype, TString ptitle, TString xtitle, TString ytitle, bool save) {

  // set RHO,eRHO,VIS,eVIS
  // visible volume is computed only if mtype is different 
  // from the previous mtype value (stored in gMTYPE)
  if(mtype!=gMTYPE) getVisibleVolume(mtype);

  if(pp_fad_units=="M") {	// Kpc -> Mpc : Kyr -> Myr
    ptitle.ReplaceAll("Kpc","Mpc");
    xtitle.ReplaceAll("Kpc","Mpc");
    ytitle.ReplaceAll("Kpc","Mpc");
    ptitle.ReplaceAll("Kyr","Myr");
    xtitle.ReplaceAll("Kyr","Myr");
    ytitle.ReplaceAll("Kyr","Myr");
  }

  int nRHO = RHO.size();
  wavearray<double> x(nRHO);
  wavearray<double> ex(nRHO);
  wavearray<double> y(nRHO);
  wavearray<double> ey(nRHO);

  if(gPTYPE=="FARvsRHO") {

    for(int n=0;n<nRHO;n++) {
      x[n]=RHO[n];
      ex[n]=eRHO[n];
      y[n]=FAR[n];
      ey[n]=eFAR[n];
    }
  }

  if((gPTYPE=="FADvsRHO")||(gPTYPE=="EVTvsRHO")||(gPTYPE=="FAPvsRHO")||
     (gPTYPE=="PRODvsFAD")||(gPTYPE=="OBSTIMEvsFAR")) {

    double Kyr = 86400.*365.*1000.;
    double Myr = 86400.*365.*1000000.;
    wavearray<double> FAD(nRHO);
    wavearray<double> eFAD(nRHO);

    if(gNZBINS<0) {		// compute classic FAD statistic + min-hold
      for(int n=0;n<nRHO;n++) {
        FAD[n]=FAR[n]/VIS[n];	
        eFAD[n]=sqrt(pow(eFAR[n]/FAR[n],2)+pow(eVIS[n]/VIS[n],2))*FAD[n];
      } 

      // apply min hold to make FAD monothonic  
      wavearray<double> mFAD(nRHO);
      wavearray<double> meFAD(nRHO);
      mFAD=0;mFAD[0]=FAD[0];meFAD[0]=eFAD[0];
      for(int n=1;n<nRHO;n++) {
        if(FAD[n]<mFAD[n-1]) {mFAD[n]=FAD[n];meFAD[n]=eFAD[n];} 
        else                 {mFAD[n]=mFAD[n-1];meFAD[n]=meFAD[n-1];}
      }
      for(int n=0;n<nRHO;n++) {FAD[n]=mFAD[n];eFAD[n]=meFAD[n];}
    } else if(gNZBINS==0) {	// compute standard FAD statistic
      for(int n=0;n<nRHO;n++) {
        FAD[n]=dFAR[n]/VIS[n];	
        eFAD[n]=sqrt(pow(eFAR[n]/FAR[n],2)+pow(eVIS[n]/VIS[n],2))*FAD[n];
      } 
      for(int n=nRHO-1;n>0;n--) FAD[n-1]+=FAD[n];
    } else if(gNZBINS>0) {	// compute classic FAD statistic + max hold
      for(int n=0;n<nRHO;n++) {
        FAD[n]=FAR[n]/VIS[n];	
        eFAD[n]=sqrt(pow(eFAR[n]/FAR[n],2)+pow(eVIS[n]/VIS[n],2))*FAD[n];
      } 
      // find first gNZBINS consecutive zero FAR bin = zbin : start from n=1 because dFAR[0]=0
      int nzbins=0;
      for(int n=1;n<nRHO;n++) {
        if(dFAR[n]==0) {
          nzbins++;
          if(nzbins==gNZBINS) {
            nRHO=n-(gNZBINS-1);break;
          }
        } else nzbins=0;
      }
      x.resize(nRHO);
      // apply max hold to make FAD monothonic  
      // max hold FAD is applied only below bin<zbin
      wavearray<double> mFAD(nRHO);
      mFAD=0;mFAD[nRHO-1]=FAD[nRHO-1];
      for(int n=nRHO-2;n>=0;n--) {
        if(FAD[n]>mFAD[n+1]) mFAD[n]=FAD[n]; else mFAD[n]=mFAD[n+1];
      }
      for(int n=0;n<nRHO;n++) FAD[n]=mFAD[n];
    }

    double Xyr = (pp_fad_units=="M") ? Myr : Kyr;
    // sec -> Kyr/Myr
    for(int n=0;n<nRHO;n++) {FAD[n]*=Xyr;eFAD[n]*=Xyr;}
 
    for(int n=0;n<nRHO;n++) {
      x[n]=RHO[n];
      ex[n]=eRHO[n];
      if(gPTYPE=="FADvsRHO") { y[n]=FAD[n]; ey[n]=eFAD[n]; }
      if(gPTYPE=="PRODvsFAD") { 
        x[n]=FAD[n]; ex[n]=eFAD[n]; 
        // if FAD not change the productivity not change  
        if(n && (FAD[n]==FAD[n-1])) {
          y[n]=y[n-1]; ey[n]=ey[n-1]; 
        } else {
          y[n]=(gOBSTIME/Xyr)*VIS[n]; ey[n]=(gOBSTIME/Xyr)*eVIS[n]; 
        }
      }
      if(gPTYPE=="EVTvsRHO") { 
        // if FAD not change the expected events not change  
        if(n && (FAD[n]==FAD[n-1])) {
          y[n]=y[n-1]; ey[n]=ey[n-1]; 
        } else {
          y[n]=FAD[n]*(gOBSTIME/Xyr)*VIS[n]; ey[n]=0; 
        }
      }
      if(gPTYPE=="FAPvsRHO") { 
        // if FAD not change the FAP not change  
        if(n && (FAD[n]==FAD[n-1])) {
          y[n]=y[n-1]; ey[n]=ey[n-1]; 
        } else {
          y[n]=FAD[n]*(gOBSTIME/Xyr)*VIS[n];
          y[n]=1.-exp(-y[n]); 
        }
        ey[n]=0; 
      }
    }
  }

  if(gPTYPE=="OBSTIMEvsFAR") { 
    for(int n=0;n<nRHO;n++) {
      x[n]=FAR[n]; ex[n]=eFAR[n]; 
      y[n]=gOBSTIME; ey[n]=0; 
    }
  }

  if(gPTYPE=="VISvsRHO") {
 
    for(int n=0;n<nRHO;n++) {
      x[n]=RHO[n];
      ex[n]=eRHO[n];
      y[n]=VIS[n];
      ey[n]=eVIS[n];
    }
  }

  if(gPTYPE=="REFFvsRHO") {

    for(int n=0;n<nRHO;n++) {
      x[n]=RHO[n];
      ex[n]=eRHO[n];
      y[n]=pow(3*VIS[n]/(4.*TMath::Pi()),1./3.);
      ey[n]=y[n]/3.*(eVIS[n]/VIS[n]);
    }
  }

  if(gPTYPE=="REFFvsIRHO") {

    for(int n=0;n<nRHO;n++) {
      x[n]=1./RHO[n];
      ex[n]=eRHO[n] ? 1./eRHO[n] : 0.;
      y[n]=pow(3*VIS[n]/(4.*TMath::Pi()),1./3.);
      ey[n]=y[n]/3.*(eVIS[n]/VIS[n]);
    }
  }

  for(int n=0;n<nRHO;n++) {
    //cout << n << " " << x[n] << " " << ex[n] << " " << y[n] << " " << ey[n] << endl;
  }

  char title[256];

#ifdef APPLY_ENORM
  if(mtype==0) {
    sprintf(title,"%s : %s", ptitle.Data(),"ALL");
  } else {
    if(normHRSS.size()) {
      sprintf(title,"%s : %s : Strain @ 10Kpc = %g",
              ptitle.Data(),gMDC->wfList[mtype-1].name.Data(),normHRSS[mtype-1]);
    } else {
      sprintf(title,"%s : %s",
              ptitle.Data(),gMDC->wfList[mtype-1].name.Data());
    }
  }
#else
  if(mtype==0) {
    cout << "MakeFAD.C : Error - mtype=0 non allowed when data are not normalized" << endl;
    gSystem->Exit(1);
  }
  sprintf(title,"%s : %s : Strain @ 10Kpc = %g",
          ptitle.Data(),gMDC->wfList[mtype-1].name.Data(),mdcHRSS[mtype-1]);
#endif

  if(gPTYPE=="EVTvsRHO") sprintf(title,"%s : %s", ptitle.Data(),"ALL");

  TGraphErrors* gr = new TGraphErrors;
  int cnt=0;
  gr->SetPoint(cnt,x[0],y[0]);
  gr->SetPointError(cnt++,0,0);
  for(int i=1;i<nRHO;i++) {
    double dx = (x[i]-x[i-1])/10000.;
    gr->SetPoint(cnt,x[i],y[i-1]);
    //gr->SetPointError(cnt++,0,ey[i-1]);
    gr->SetPointError(cnt++,ex[i],ey[i-1]);
    gr->SetPoint(cnt,x[i]+dx,y[i]);
    //gr->SetPointError(cnt++,0,ey[i]);
    gr->SetPointError(cnt++,ex[i],ey[i]);
  }
  gr->GetHistogram()->GetXaxis()->SetTitle(xtitle.Data());
  gr->GetHistogram()->GetYaxis()->SetTitle(ytitle.Data());
  gr->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  gr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);

  // convert logx,logy to linx,liny when x,y ranges are less than a 1 order of magnitude
  double xmax=0;double xmin=1e80;
  double ymax=0;double ymin=1e80;
  for(int i=0;i<nRHO;i++) {if(x[i]>xmax) xmax=x[i]; if(x[i]<xmin && x[i]!=0) xmin=x[i];}
  for(int i=0;i<nRHO;i++) {if(y[i]>ymax) ymax=y[i]; if(y[i]<ymin && y[i]!=0) ymin=y[i];}
  if(xmax/xmin<10) gCANVAS->SetLogx(false);
  if(ymax/ymin<10) gCANVAS->SetLogy(false);
  gr->GetHistogram()->GetXaxis()->SetRangeUser(xmin,xmax);
  gr->GetHistogram()->GetYaxis()->SetRangeUser(0.95*ymin,1.05*ymax);
 
  gr->SetTitle(title);
  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  gr->SetFillColor(2);
  gr->SetFillStyle(3003);
  gr->Draw("3AL");

  if(gFIT) {
    TF1 *gfit=NULL;

    if(gPTYPE=="VISvsRHO")   gfit = new TF1("gfit","[0]+[1]/pow(x,3)",x[0],x[nRHO-1]);
    if(gPTYPE=="REFFvsRHO")  gfit = new TF1("gfit","[0]/x",x[0],x[nRHO-1]);
    if(gPTYPE=="REFFvsIRHO") gfit = new TF1("gfit","[0]*x",x[0],x[nRHO-1]);

    if((gPTYPE=="VISvsRHO") || (gPTYPE=="REFFvsRHO") || (gPTYPE=="REFFvsIRHO")) {
      gr->Fit(gfit,"R");
      gfit=gr->GetFunction("gfit");
      gfit->SetLineColor(kBlue);
      gfit->SetLineWidth(1);
    }
  }

  if(!save) return gr;

  // print plot
  char ofname[256];
  if(mtype==0) {
    sprintf(ofname,"%s/%s_%s.png",
            gODIR.Data(),gPTYPE.Data(),"ALL");
  } else {
    sprintf(ofname,"%s/%s_%s.png",
            gODIR.Data(),gMDC->wfList[mtype-1].name.Data());
  }
  cout << ofname << endl;
  gCANVAS->Print(ofname);
  if(mtype==0) {
     sprintf(ofname,"%s/%s_%s.txt",
             gODIR.Data(),gPTYPE.Data(),"ALL");
     cout << "Write values on: " << ofname << endl;
     ofstream outf;
     outf.open(ofname);
     //for (int index=0; index<x.size(); index++) {out << x.data[index] << "\t" << y.data[index] << "\t" << ex.data[index] << " " << ey.data[index] << endl;}
     for (int index=0; index<x.size(); index++) {outf << x.data[index] << "\t" << y.data[index] << endl;}
     outf.close();
  }

  return gr;
}

void RECvsINJ(int itype, TString ptitle, TString xtitle, TString ytitle) {

  char cut[256];
  char title[256];

  if(pp_fad_units=="M")  xtitle.ReplaceAll("Kpc","Mpc");

  // INJECTED
  float ufactor;
  if(pp_fad_units=="M") {	// Kyr -> Myr
    ufactor=1;
  } else {			// Kpc -> Mpc 
    ufactor=1000;
  }
  if(itype==0) {
    sprintf(cut,"1==1");
  } else {
    sprintf(cut,"type==%d",itype);
  }

  float dmin = gTRMDC->GetMinimum("distance");
  float dmax = gTRMDC->GetMaximum("distance");
#ifdef APPLY_ENORM
  double enorm_max=0;
  for(int i=0;i<normHRSS.size();i++) {
    double enorm=normHRSS[i]/mdcHRSS[i];
    if(enorm>enorm_max) enorm_max=enorm;
  }
#endif
  if(enorm_max==0) enorm_max=1;
  dmin*=ufactor*enorm_max;
  dmax*=ufactor*enorm_max;

  TH1F* mhist = new TH1F("mhist","mhist",100,dmin,dmax);

  TTreeFormula mcut("mcut", cut, gTRMDC);
  gTRMDC->SetBranchStatus("*",false);
  gTRMDC->SetBranchStatus("distance",true);
  gTRMDC->SetBranchStatus("type",true);
  int mtype;
  float distance;
  gTRMDC->SetBranchAddress("distance",&distance);
  gTRMDC->SetBranchAddress("type",&mtype);
  int msize = gTRMDC->GetEntries();
  int nmdc = 0;
  for(i=0;i<msize;i++) {
    gTRMDC->GetEntry(i);
    if(mcut.EvalInstance()==0) continue;
    distance*=ufactor; 
#ifdef APPLY_ENORM
    // compute Normalization to hrss @ 10 Kpc
    if(normHRSS.size()&&(normHRSS[mtype-1]>0)) {
      distance *= normHRSS[mtype-1]/mdcHRSS[mtype-1];
    }
#endif
    mhist->Fill(distance);
    nmdc++;
  }
  gTRMDC->SetBranchStatus("*",true);
  cout << "nmdc : " << nmdc << endl; 

  mhist->Draw("HIST");
  mhist->GetXaxis()->SetTitle(xtitle);
  mhist->GetYaxis()->SetTitle(ytitle);
  mhist->GetYaxis()->SetRangeUser(0.1,pow(10.,TMath::Ceil(TMath::Log10(mhist->GetMaximum()))));

  // DETECTED
  if(itype==0) {
    sprintf(cut,"abs(time[0]-time[%d])<%g && netcc[%d]>%g && rho[%d]>%g", 
            nIFO,T_win,pp_inetcc,T_cor,pp_irho,T_cut);
  } else {
    sprintf(cut,"type[1]==%d && abs(time[0]-time[%d])<%g && netcc[%d]>%g && rho[%d]>%g", 
            itype,nIFO,T_win,pp_inetcc,T_cor,pp_irho,T_cut);
  }
  if(T_vED>0) sprintf(cut,"%s && neted[0]/ecor<%f",cut,T_vED);
  if(T_pen>0) sprintf(cut,"%s && penalty>%f",cut,T_pen);

  TH1F* whist = new TH1F("whist","whist",100,dmin,dmax);

  TTreeFormula wcut("wcut", cut, gTRWAVE);
  gTRWAVE->SetBranchStatus("*",false);
  gTRWAVE->SetBranchStatus("range",true);
  gTRWAVE->SetBranchStatus("type",true);
  gTRWAVE->SetBranchStatus("time",true);
  gTRWAVE->SetBranchStatus("neted",true);
  gTRWAVE->SetBranchStatus("ecor",true);
  gTRWAVE->SetBranchStatus("netcc",true);
  gTRWAVE->SetBranchStatus("rho",true);
  gTRWAVE->SetBranchStatus("penalty",true);
  int wtype[2];
  float range[2];
  gTRWAVE->SetBranchAddress("range",range);
  gTRWAVE->SetBranchAddress("type",wtype);
  int wsize = gTRWAVE->GetEntries();
  int nwave = 0;
  for(i=0;i<wsize;i++) {
    gTRWAVE->GetEntry(i);
    if(wcut.EvalInstance()==0) continue;
    distance=ufactor*range[1]; 
#ifdef APPLY_ENORM
    // compute Normalization to hrss @ 10 Kpc
    if(normHRSS.size()&&(normHRSS[wtype[1]-1]>0)) {
      distance *= normHRSS[wtype[1]-1]/mdcHRSS[wtype[1]-1];
    }
#endif
    whist->Fill(distance);
    nwave++;
  }
  gTRWAVE->SetBranchStatus("*",true);
  cout << "nwave : " << nwave << endl; 
  whist->SetFillColor(kRed);
  whist->Draw("SAME");

  sprintf(title,"%s",cut);

  sprintf(title,"%s : inj = %d : rec : %d",ptitle.Data(),nmdc,nwave);
  mhist->SetTitle(title);

  // print plot
  char ofname[1024];
  if(itype==0) {
    sprintf(ofname,"%s/%s_%s.png",
            gODIR.Data(),gPTYPE.Data(),"ALL");
  } else {
    sprintf(ofname,"%s/%s_%s.png",
            gODIR.Data(),gPTYPE.Data(),gMDC->wfList[itype-1].name.Data());
  }
  cout << ofname << endl;
  gCANVAS->Print(ofname);

  // compute efficiency vs distance
  TH1F* ehist = new TH1F("ehist","ehist",100,dmin,dmax);
  for(int i=1;i<=mhist->GetNbinsX();i++) {
    double ndet = whist->GetBinContent(i);
    double ninj = mhist->GetBinContent(i);
    double eff = ninj ? ndet/ninj : 0.;
    ehist->SetBinContent(i,eff);
  }
  ehist->Draw("HIST");
  ehist->GetXaxis()->SetTitle(xtitle);
  ehist->GetYaxis()->SetTitle("efficiency");
  //ehist->GetYaxis()->SetRangeUser(0.1,pow(10.,TMath::Ceil(TMath::Log10(ehist->GetMaximum()))));
  ehist->SetTitle("Efficiency vs Distance");
  // print efficiency plot
  if(itype==0) {
    sprintf(ofname,"%s/%s_%s.png",
            gODIR.Data(),"EFFvsDIST","ALL");
  } else {
    sprintf(ofname,"%s/%s_%s.png",
            gODIR.Data(),"EFFvsDIST",gMDC->wfList[itype-1].name.Data());
  }
  cout << ofname << endl;
  gCANVAS->Print(ofname);

  delete mhist;
  delete whist;
  delete ehist;

  return;
}

void 
MakeHtml(TString ptitle) {

  char index_file[1024];

  sprintf(index_file,"%s/%s/fad.html",gODIR.Data(), gPTYPE.Data());

  ofstream out;
  out.open(index_file,ios::out);             // create index.html
  char oline[1024];
  out << "<html>" << endl;
  out << "<font color=\"black\" style=\"font-weight:bold;\"><center><p><h2>"
      << ptitle.Data() << "</h2><p><center></font>" << endl;
  out << "<hr><br>" << endl;
  int offset=1;
#ifdef APPLY_ENORM
  offset=0;
#endif
  int size=gMDC->wfList.size();
  if(gPTYPE=="EVTvsRHO") {offset=0; size=0;}

  for(int i=offset;i<=size;i++) {
    out << "<table border=0 cellpadding=2 align=\"center\">" << endl;
    out << "<tr align=\"center\">" << endl;
    if(i==0) {
      out << "<td><font style=\"font-weight:bold;\"><center><p><h2>"
          << "ALL" << "</h2><p><center></font></td>" << endl;
    } else {
      out << "<td><font style=\"font-weight:bold;\"><center><p><h2>"
          << gMDC->wfList[i-1].name.Data() << "</h2><p><center></font></td>" << endl;
    }
    out << "</tr>" << endl;
    out << "<tr align=\"center\">" << endl;
    char fname[1024];
    if(i==0) {
      sprintf(fname,"%s_%s_%s.png",gPTYPE.Data(),data_label,"ALL");
    } else {
      sprintf(fname,"%s_%s_%s.png",gPTYPE.Data(),data_label,gMDC->wfList[i-1].name.Data());
    }
    sprintf(oline,"<td><a href=\"%s\"><img src=\"%s\" width=650></a></td>",fname,fname);
    out << oline << endl;
    out << "</tr>" << endl;
    out << "</table>" << endl;
  }
  out << "</html>" << endl;
  out.close();

  return;
}

void 
MakeHtmlIndex(TString* pTYPE, int pSIZE) {

  // add EFFvsDIST plot to the pTYPE plot list
  int psize=pSIZE+1;
  TString* ptype = new TString[psize+1]; 
  int I=0;
  for(int i=0;i<pSIZE;i++) {
    ptype[I++]=pTYPE[i];
    if(pTYPE[i]=="RECvsINJ") ptype[I++]="EFFvsDIST";
  }

  TString info="";

  char index_file[1024];
  sprintf(index_file,"%s/fad.html",gODIR.Data());

  ofstream out;
  out.open(index_file,ios::out);             // create index.html
  char oline[1024];
  out << "<html>" << endl;
  out << "<font color=\"black\" style=\"font-weight:bold;\"><center><p><h1>"
      << pp_fad_title << "</h1><p></center></font>" << endl;
  if(pp_fad_subtitle!="") {
    out << "<font color=\"red\" style=\"font-weight:bold;\"><center><p><h3>"
        << pp_fad_subtitle << "</h3><p></center></font>" << endl;
  }
  out << "<hr>" << endl;
  out << "<h3>FAD Options</h3>" << endl;

  TString odir = TString(gSystem->WorkingDirectory())+"/"+gODIR;
  out << "<ul>" << endl;
  out << "<li> odir : <font color=\"red\"> " << odir << "</font>" << endl; 
  out << "</ul>" << endl;

  out << "<ul>" << endl;
  out << "<li> --bkgrep : <font color=\"red\"> " << pp_fad_bkgrep << "</font>" << endl; 
  out << "</ul>" << endl;

  if(pp_fad_hrss.IsFloat()) {
     if(pp_fad_hrss=="0")	info = " (hrss rescale is not applied)";
     else 			info = " (apply a constant hrss rescale)";
  } else 			info = " (hrss is rescaled using custom hrss factors)";

  out << "<ul>" << endl;
  out << "<li> --hrss : <font color=\"red\"> " << pp_fad_hrss << "</font>" << info << endl; 
  out << "</ul>" << endl;

  if(gFIT)	info = " (apply fit)";
  else		info = " (fit not applied)";

  out << "<ul>" << endl;
  out << "<li> --gfit : <font color=\"red\"> " << gFIT << "</font>" << info << endl; 
  out << "</ul>" << endl;

  out << "<ul>" << endl;
  out << "<li> --rhomin : <font color=\"red\"> " << gRHOMIN << "</font>" 
      << " (rho minimum used to compute FAD plots)" << endl; 
  out << "</ul>" << endl;

  if(gNZBINS<0)       info = " (un-biased FAD)";
  else if(gNZBINS==0) info = " (standard FAD)";
  else if(gNZBINS>0)  info = " (test FAD)";

  out << "<ul>" << endl;
  out << "<li> --nzbins : <font color=\"red\"> " << gNZBINS << "</font>" << info << endl; 
  out << "</ul>" << endl;

  if(pp_fad_units=="K")	info = " (Kpc,Kyr)";
  else 			info = " (Mpc,Myr)";

  out << "<ul>" << endl;
  out << "<li> --units : <font color=\"red\"> " << pp_fad_units << "</font>" << info << endl; 
  out << "</ul>" << endl;

  if(pp_fad_distr=="FORMULA")	info = " (radial distribution is computed from formula)";
  else				info = " (radial distribution is computed from mdc injections)";

  out << "<ul>" << endl;
  out << "<li> --distr : <font color=\"red\"> " << pp_fad_distr << "</font>" << info << endl; 
  out << "</ul>" << endl;

  out << "<ul>" << endl;
  out << "<li> --nbins : <font color=\"red\"> " << gNBINS << "</font>"
      << " (number of bins used in hist to computed the radial distribution from the MDC injections root file)" << endl; 
  out << "</ul>" << endl;

  out << "<hr><br>" << endl;

  TString _ptype[nPLOT+1];
  int index=0;
  _ptype[index++]="FAPvsRHO";
  for(int i=0;i<nPLOT+1;i++) if(ptype[i]!="FAPvsRHO") _ptype[index++]=ptype[i];

  for(int i=0;i<psize;i++) {
    out << "<table border=0 cellpadding=2 align=\"center\">" << endl;
    out << "<tr align=\"center\">" << endl;
    out << "</tr>" << endl;
    out << "<tr align=\"center\">" << endl;
    char fname[1024];
    if(i==0) {
      sprintf(fname,"%s_%s.png",_ptype[i].Data(),"ALL");
      sprintf(oline,"<td><a href=\"%s\"><img src=\"%s\" width=650></a></td>",fname,fname);
      out << oline << endl;
    } else {
      sprintf(fname,"%s_%s.png",_ptype[i].Data(),"ALL");
      sprintf(oline,"<td><a href=\"%s\"><img src=\"%s\" width=520></a></td>",fname,fname);
      out << oline << endl;
      sprintf(fname,"%s_%s.png",_ptype[i+1].Data(),"ALL");
      sprintf(oline,"<td><a href=\"%s\"><img src=\"%s\" width=520></a></td>",fname,fname);
      out << oline << endl;
    } 
    out << "</tr>" << endl;
    out << "</table>" << endl;
    if(i==0) out << "<hr><br>" << endl; else i++;
  }
  out << "</html>" << endl;
  out.close();

  return;
}

void AddHtmlHeader(TString odir) {

  char cmd[1024];

  sprintf(cmd,"root -n -l -b ${CWB_ROOTLOGON_FILE} '${HOME_CWB}/info/macros/MakeHTML.C\(\"%s\",false\)'",odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp ${HOME_WAT}//html/etc/html/ROOT.css %s",odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp ${HOME_WAT}//html/etc/html/ROOT.js %s",odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp ${HOME_WAT}//html/etc/html/tabber.css %s",odir.Data());
  gSystem->Exec(cmd);
  sprintf(cmd,"cp ${HOME_WAT}//html/etc/html/tabber.js %s",odir.Data());
  gSystem->Exec(cmd);
}

void MakeMultiPlotFAD(TString mdc_inj_file) {

  #define NTYPE_MAX 20
  #define NMDC_MAX  64

  // read injection file types
  char    imdc_set[NMDC_MAX][128];    // injection set
  size_t  imdc_type[NMDC_MAX];        // injection type
  char    imdc_name[NMDC_MAX][128];   // injection name
  double  imdc_fcentral[NMDC_MAX];    // injection central frequencies
  double  imdc_fbandwidth[NMDC_MAX];  // injection bandwidth frequencies
  size_t  imdc_index[NMDC_MAX];       // type reference array           
  size_t  imdc_iset[NMDC_MAX];        // injection set index            

  int ninj=ReadInjType(mdc_inj_file.Data(),NMDC_MAX,imdc_set,imdc_type,
                       imdc_name,imdc_fcentral,imdc_fbandwidth);
  if(ninj==0) {cout << "Error - no injection - terminated" << endl;exit(1);}

  TString* imdc_set_name = new TString[ninj];
  int nset=0;
  for(int i=0;i<ninj;i++) {
    bool bnew=true;
    for(int j=0;j<nset;j++) if(imdc_set[i]==imdc_set_name[j]) bnew=false;
    if(bnew) imdc_set_name[nset++]=imdc_set[i];
  }
  cout << "nset : " << nset << endl;

  for(int i=0;i<nset;i++) {
    for(int j=0;j<ninj;j++) if(imdc_set[j]==imdc_set_name[i]) imdc_iset[j]=i;
  }
  for (int iset=0;iset<nset;iset++) cout << iset << " " << imdc_set_name[iset].Data() << endl;

  Color_t colors[NMDC_MAX] = { 6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7,
                               6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7,
                               6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7,
                               6, 3, 2, 8,43, 7, 8, 4, 4, 2,43, 1, 3, 1, 6, 7 };
  Style_t markers[NMDC_MAX]= {20,21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,
                              21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,20,
                              21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,20,
                              21,22,23,29,29,21,22,20,21,22,23,29,20,21,22,20 };

  char title[256];
  char ofile[256];
  TCanvas*      canvas[NTYPE_MAX];
  TMultiGraph*  mg[NTYPE_MAX];
  TGraphErrors* gr[NMDC_MAX];
  TLegend*      legend[NTYPE_MAX];

  char rho_type[32];sprintf(rho_type,"rho[%d]",pp_irho);

  for(int iset=0;iset<nset;iset++) {

    int iset_size=0;    // iset size
    for(i=0; i<ninj; i++) if(imdc_iset[i]==iset) iset_size++;
    double hleg = 0.88-iset_size*0.08;
    legend[iset] = new TLegend(0.58,hleg,0.99,0.88,NULL,"brNDC");
    legend[iset]->AddEntry((TObject*)0, "", "");
 
    mg[iset] = new TMultiGraph();

    // set plot type 
    gPTYPE = "FADvsRHO";
    TString ptitle=imdc_set_name[iset];
    TString xtitle=rho_type;
    TString ytitle="FAD (Kpc^{-3} Kyr^{-1})";
    if(pp_fad_units=="M") ytitle.ReplaceAll("Kpc","Mpc");
    if(pp_fad_units=="M") ytitle.ReplaceAll("Kyr","Myr");

    double ymin=1e40;
    double ymax=0;
    for(int i=0; i<ninj; i++) {
      if(imdc_iset[i]!=iset) continue;  // skip unwanted injection types
      cout << "MakeMultiPlotFAD : Processing -> " << imdc_name[i] << endl;
      gr[i] = Plot(i+1, ptitle, xtitle, ytitle, false);
      gr[i]->SetLineColor(colors[i]);
      int N = gr[i]->GetN(); 
      double* Y = gr[i]->GetY(); 
      for(int k=0;k<N;k++) if(Y[k]>0 && Y[k]<ymin) ymin=Y[k];
      for(int k=0;k<N;k++) if(Y[k]>0 && Y[k]>ymax) ymax=Y[k];
      double fad_rho8 = gr[i]->Eval(8);	// get FAD @ rho=8
      char legname[32];sprintf(legname,"%s, %5.1E",imdc_name[i], fad_rho8);
      legend[iset]->AddEntry(gr[i],legname,"lp");
      mg[iset]->Add(gr[i]);
    }

    char namecanvas[8];
    sprintf(namecanvas,"c%i",iset);
    canvas[iset] = new TCanvas(namecanvas,namecanvas,125+iset*10,82,976,576);
    canvas[iset]->Clear();
    canvas[iset]->ToggleEventStatus();
    canvas[iset]->SetLogy(true);
    canvas[iset]->SetLogx(true);
    canvas[iset]->SetGridx();
    canvas[iset]->SetGridy();
    canvas[iset]->SetFillColor(kWhite);
    canvas[iset]->cd();

    mg[iset]->SetTitle(ptitle);
    mg[iset]->Paint("alp");
    mg[iset]->GetHistogram()->GetXaxis()->SetTitle(xtitle);
    mg[iset]->GetHistogram()->GetYaxis()->SetTitle(ytitle);
    mg[iset]->GetHistogram()->GetXaxis()->CenterTitle(true);
    mg[iset]->GetHistogram()->GetXaxis()->SetLabelFont(42);
    mg[iset]->GetHistogram()->GetXaxis()->SetTitleFont(42);
    mg[iset]->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
    mg[iset]->GetHistogram()->GetYaxis()->SetLabelFont(42);
    mg[iset]->GetHistogram()->GetYaxis()->SetTitleFont(42);
    mg[iset]->GetHistogram()->GetXaxis()->SetRangeUser(RHO[0],RHO[RHO.size()-2]);
    mg[iset]->GetHistogram()->GetYaxis()->SetRangeUser(0.90*ymin,1.1*ymax);
    mg[iset]->GetHistogram()->GetXaxis()->SetMoreLogLabels();

    mg[iset]->Draw("alp");

    char leg_header[32];sprintf(leg_header,"FAD @ rho[%d]=8",pp_irho);
    legend[iset]->SetHeader(leg_header);
    legend[iset]->SetBorderSize(1);
    legend[iset]->SetTextAlign(22);
    legend[iset]->SetTextFont(12);
    legend[iset]->SetLineColor(1);
    legend[iset]->SetLineStyle(1);
    legend[iset]->SetLineWidth(1);
    legend[iset]->SetFillColor(0);
    legend[iset]->SetFillStyle(1001);
    legend[iset]->SetTextSize(0.04);
    legend[iset]->SetLineColor(kBlack);
    legend[iset]->SetFillColor(kWhite);

    legend[iset]->Draw();

    sprintf(ofile,"%s/fad_rho_%s.gif",gODIR.Data(),imdc_set_name[iset].Data());
    cout << ofile << endl;
    canvas[iset]->SaveAs(ofile);
  }

  return;
}
