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


// ---------------------------------------------------------------------------------
// INFOS
// ---------------------------------------------------------------------------------

// This macro is used to tune the PP cut thresholds
//
// how to use - Ex: root -l -b TunePPCutsTool_Config.C 'TunePPCutsTool.C(100,2000,100,30,"output_rnd_walk_cuts.out")'
//
// where output_rnd_walk_cuts.out is the output random walk cut values
// 
// where TunePPCutsTool_Config.C is the configuration file. For example:
//

/*
{
  #include <O3/SEARCHES/OFFLINE/BBH/LH/PP_Cuts.hh>

  #define IFILE_BKG     "wave_O3_K15_C01_LH_BBH_BKG_xrun2.M1.V_hvetoLH.C_rho1_gt_6.root"
  #define IFILE_SIM     "wave_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.V_hvetoLH.root"  

  #define BKG_LIVE_TIME 35266962627.00  // non-zero lags : 77934 lags - 35266962627.00 sec = 1118.3 years

  // definition of the basic selection cuts used in O3b run
  TCut bin_cut_standard = bin1_cut;
  TCut bin_cut_basic = TCut("bin_cut_basic",(dqveto+lveto_cut).GetTitle());

  #define PP_CUTS "                                                       \n\
  #                                                                       \n\
  # name          cmp     mean    sigma   min     max     enabled_tune    \n\
  #                                                                       \n\
  netcc[0]        >       0.7     0.1     0.5     0.9     1               \n\
  netcc[2]        >       0.7     0.1     0.5     0.9     1               \n\
  norm            >       2.5     1.0     1.0     5.0     1               \n\
  Qveto[0]        >       0.25    0.1     0.0     0.5     1               \n\
  log10(penalty)  <       0.2     0.1     0.1     0.5     1               \n\
  frequency[0]    >       60.0    10.0    40.0    80.0    1               \n\
  frequency[0]    <       300.0   50.0    150.0   350.0   1               \n\
  chirp[1]        >       1.0     3.0     0.5     5.5     1               \n\
  chirp[1]        <       70.0    10.0    50.0    90.0    1               \n\
  #Qveto[2]       >       0.0     1.0     0.0     5.0     1               \n\
  "
}
*/

// ----------------------------------------------------------------------------
// includes
// ----------------------------------------------------------------------------

#include <boost/any.hpp>

// ----------------------------------------------------------------------------
// defines
// ----------------------------------------------------------------------------

#define YEAR    	(24*3600*365.)

#define nCUT_MAX	20	// number of PP cut max defined in the input config file

#define nBKG_MAX	10	// number of loudest background listed 

#define nIFAR		14	// number of ifar factor used to print inital/final detections

#define SIGMA_FACTOR	1.	// factor used to rescale input sigma

#define CHECK_DNEVT_VS_IFAR	// enable dnevt vs IFAR check -> the increment of detected events must be >=0 for each IFAR
				// see CheckIncrementalDetectionsVsIFAR function

// ----------------------------------------------------------------------------
// global variables
// ----------------------------------------------------------------------------

int 	gNBKG;			// number of background events for the correnspondig IFAR
int 	gNCTRY;			// number of random trials for central tuning
int 	gNIFTRY;		// number of random trials for inital/final tuning
int 	gNCUT;			// number of PP cuts
int	gNMDC;			// readed from the input mdc files

TString gCUT_NAME[nCUT_MAX];	// PP cut name
TString gCUT_CMP[nCUT_MAX];	// compare symbol '< , >, ...'
double  gCUT_MEA[nCUT_MAX];	// tuning mean     threshold (= standard PP cut threshold)
double  gCUT_SIG[nCUT_MAX];	// tuning sigma    threshold
double  gCUT_MIN[nCUT_MAX];	// tuning min      threshold
double  gCUT_MAX[nCUT_MAX];	// tuning max      threshold

TCut 	gZL_CUT;		// zero lag selection cut

TString gOFNAME;		// file name used to output the random walk of cut theresholds

// ----------------------------------------------------------------------------
// functions
// ----------------------------------------------------------------------------

int  GetSelectedEvents(TTree* tree_BKG, TTree* tree_SIM, TString user_cut, int ncut, double* cut_thr, double& rho_thr);
void PrintLoudest(TTree* tree_BKG, TString user_cut, int ncut, double* cut_thr, double rho_thr, TString label="", int nLoudest=nBKG_MAX);
int  ReadCutList(std::stringstream* in, TString* name, TString* cmp, double* mean, double* sigma, double* min, double* max);
void PrintPPCuts(int nevt_new, int nevt_old, int ncut, double* cut_thr, double rho_thri, TString user_cut, bool bppcuts, bool beff, bool bbcut);
void DumpIncrementalDetectionsVsIFAR(TTree* tree_BKG, TTree* tree_SIM, TString user_cut, double* cut_thr, TString ofname, bool app);
bool CheckIncrementalDetectionsVsIFAR(TTree* tree_BKG, TTree* tree_SIM, TString user_cut, double* cut_thr);
void PrintCovariance(TTree* tree, bool norm=true);
void DumpPPCuts(TString ofname, int ncut, double* cut_thr, bool app);

void TunePPCutsTool(double IFAR=10.0, int nCTRY=1000, int nIFTRY=200, int seed=0, TString ofname="") {
//
// IFAR		: input IFAR used to tunong the PP cut thresholds
// nCTRY	: input trials for central tuning
// nIFTRY	: input trials for initial and final tuning
// seed		: seed used for random number generator (optional)
// ofname       : file name used to output the random walk of cut theresholds, id ="" the output file is disabled
//		  the extension must be '.out'
//

  // check ofname extension
  if(ofname!="") {
    if(ofname.EndsWith(".out")) {
    } else {
      cout << endl << "TunePPCutsTool.C - Error : ofname extension must be '.out' !!!" << endl << endl;exit(1);
    }
  }

  gOFNAME  = ofname.ReplaceAll(".out",TString::Format("_IFAR_%.2f.out",IFAR));
  gNBKG    = BKG_LIVE_TIME/IFAR/YEAR;
  gNCTRY   = nCTRY;
  gNIFTRY  = nIFTRY;

  if(gNBKG==0) {cout << endl << "IFAR too high -> nBKG = 0, macro aborted" << endl << endl;exit(1);}

  cout << endl;
  cout << "--------------------------------------------------------------------------------------------" << endl;
  cout << " Input Parameters                                                                           " << endl;
  cout << "--------------------------------------------------------------------------------------------" << endl;
  cout << endl;
  cout << " IFAR = " << IFAR << " (years) ->" << "\tnBKG = " << gNBKG << "\tnCTRY = " 
       << gNCTRY << "\tnIFTRY = " << gNIFTRY << "\tseed = " << seed << endl;
  cout << endl;
  cout << " gOFNAME = " << gOFNAME << endl;
  cout << endl;

  gErrorIgnoreLevel = kError;  // disable tree selection warning messages

  // ---------------------------------------------------------------
  // Open input background ROOT file
  // ---------------------------------------------------------------

  TString ifile_BKG = IFILE_BKG;
  cout << " " << ifile_BKG << endl; 

  TFile *rfile_BKG = TFile::Open(ifile_BKG);
  if(rfile_BKG==NULL) {cout << "Error opening file: " << ifile_BKG << endl;exit(1);} 

  TTree* tree_BKG = (TTree *) gROOT->FindObject("waveburst");
  if(tree_BKG==NULL) {cout << "waveburst tree not found !!!" << endl;exit(1);}

  // Get detector names
  TString net="";
  vector<TString> ifo;
  TList* ifoList = tree_BKG->GetUserInfo();
  detectorParams dParams[NIFO_MAX];
  for (int n=0;n<ifoList->GetSize();n++) {
    detector* pDetector = (detector*)ifoList->At(n);
    dParams[n] = pDetector->getDetectorParams();
    ifo.push_back(dParams[n].name);
    net+=TString(dParams[n].name);
    //pDetector->print();
  }
  if(ifo.size()==0) {
    cout << "ifo names are not contained in the input wave root file: " << IFILE_BKG << endl;
    gSystem->Exit(1);
  }

  int nIFO = ifoList->GetSize();	// get number of detectors in the network

  //cout << "network = " << net << endl;

  // define zero lag cut
  gZL_CUT = TCut("gZL_CUT",TString::Format("(lag[%d]==0) && (slag[%d]==0)",nIFO,nIFO).Data());

  // ---------------------------------------------------------------
  // Open input simulation ROOT file
  // ---------------------------------------------------------------

  TString ifile_SIM = IFILE_SIM;
  cout << " " << ifile_SIM << endl; 

  TFile *rfile_SIM = TFile::Open(ifile_SIM);
  if(rfile_SIM==NULL) {cout << "Error opening file: " << ifile_SIM << endl;exit(1);} 

  TTree* tree_SIM = (TTree *) gROOT->FindObject("waveburst");
  if(tree_SIM==NULL) {cout << "waveburst tree not found !!!" << endl;exit(1);}

  // ---------------------------------------------------------------
  // Open input mdc ROOT file
  // ---------------------------------------------------------------

  TString ifile_MDC = IFILE_SIM;
  ifile_MDC.ReplaceAll("wave_","mdc_");
  cout << " " << ifile_MDC << endl; 

  TFile *rfile_MDC = TFile::Open(ifile_MDC);
  if(rfile_MDC==NULL) {cout << "Error opening file: " << ifile_MDC << endl;exit(1);} 

  TTree* tree_MDC = (TTree *) gROOT->FindObject("mdc");
  if(tree_MDC==NULL) {cout << "mdc tree not found !!!" << endl;exit(1);}

  // retrive the number of injections
  gNMDC = tree_MDC->GetEntries();
  cout << endl << " gNMDC = " << gNMDC << " (injections)" << endl << endl;

  // ---------------------------------------------------------------
  // clone tree -> memory resident -> faster
  // ---------------------------------------------------------------

  gROOT->cd(0);
  TTree* mtree_BKG = tree_BKG->CloneTree();
  TTree* mtree_SIM = tree_SIM->CloneTree();

  // ---------------------------------------------------------------
  // read cut list
  // ---------------------------------------------------------------

  std::stringstream pp_cuts(PP_CUTS);
  gNCUT = ReadCutList(&pp_cuts, gCUT_NAME, gCUT_CMP, gCUT_MEA, gCUT_SIG, gCUT_MIN, gCUT_MAX);

  TString user_cut = bin_cut_basic.GetTitle();

  // print standard PP cut string
  cout << endl << " bin1_cut_standard ..." << endl << endl; 
  cout << " -> " << bin_cut_standard.GetTitle() << endl << endl; 
  cout << endl << " bin1_cut_basic ..." << endl << endl; 
  cout << " -> " << bin_cut_basic.GetTitle() << endl << endl; 

  bool answer = CWB::Toolbox::question("do you want to continue ? ");
  if(!answer) exit(1);

  double cut_thr[nCUT_MAX];
  double cut_mean[nCUT_MAX];
  double cut_sigma[nCUT_MAX];
  double rho_thr;

  // init cut_mean, cut_sigma
  for(int i=0;i<gNCUT;i++) {
    cut_mean[i]  = gCUT_MEA[i];
    cut_sigma[i] = gCUT_SIG[i]*SIGMA_FACTOR;
  }

  // ---------------------------------------------------------------
  // print covariance tree_BKG 
  // ---------------------------------------------------------------

  cout << endl << " background covariance matrix ..." << endl << endl;
  PrintCovariance(tree_BKG);
  cout << endl << " simulation covariance matrix ..." << endl << endl;
  PrintCovariance(tree_SIM);

  answer = CWB::Toolbox::question("do you want to continue ? ");
  if(!answer) exit(1);

  // ---------------------------------------------------------------
  // print detections/loudest for standard PP cuts
  // ---------------------------------------------------------------

  int nevt_std = GetSelectedEvents(mtree_BKG, mtree_SIM, bin_cut_standard.GetTitle(), 0, NULL, rho_thr);
  cout << endl << " nevt standard = " << nevt_std << endl << endl;
  if(nevt_std==0) {cout << endl << " Error - input IFAR is too low ..." << endl << endl;exit(1);}
  PrintPPCuts(nevt_std, nevt_std, gNCUT, gCUT_MEA, rho_thr, user_cut, false, true, false);

  // ---------------------------------------------------------------
  // tune PP cuts
  // ---------------------------------------------------------------

  if(gOFNAME!="") DumpPPCuts(gOFNAME, gNCUT, cut_mean, false);
  if(gOFNAME!="") DumpPPCuts(gOFNAME, gNCUT, cut_mean, true);

  if(seed>=0) gRandom->SetSeed(seed);

  double cut_thr_max[nCUT_MAX];
  double rho_thr_max=0.0;
  int nevt_max=nevt_std;

  // INITIAL TUNING

  // generate random threshold values around the cut_thr_max (in the min:max range) individually for each cut
  cout << endl << "Initial tuning ..." << endl << endl;
  if(gOFNAME!="") DumpPPCuts(gOFNAME, 0, NULL, true);	// dump blank line
  for(int k=0;k<gNCUT;k++) cut_thr_max[k]=cut_mean[k];
  for(int i=0;i<gNCUT;i++) {
    if(gCUT_SIG[i]==0) continue;	// skip fixed cuts			
    cout << " " << i << " -> " << gCUT_NAME[i] << " ..." << endl;
    for(int k=0;k<gNCUT;k++) cut_thr[k]=cut_thr_max[k];
    for(int j=0;j<gNIFTRY;j++) {
      cut_thr[i]=gRandom->Uniform(gCUT_MIN[i],gCUT_MAX[i]);

      int nevt_new = GetSelectedEvents(mtree_BKG, mtree_SIM, user_cut, gNCUT, cut_thr, rho_thr);
#ifdef CHECK_DNEVT_VS_IFAR
      bool check_dnevt = (nevt_new>nevt_max) ? CheckIncrementalDetectionsVsIFAR(tree_BKG, tree_SIM, user_cut, cut_thr) : true;
#else
      bool check_dnevt = true;
#endif
      if(nevt_new>nevt_max && check_dnevt==true) {

	cout << endl << "  " << gCUT_NAME[i] << " : standard = " << gCUT_MEA[i] << " -> tuned = " << cut_thr[i] << endl;
        PrintPPCuts(nevt_new, nevt_std, gNCUT, cut_thr, rho_thr, user_cut, false, true, false);

        // update cut_mean
        for(int i=0;i<gNCUT;i++) cut_mean[i]=cut_thr[i];

        // update max values
        for(int i=0;i<gNCUT;i++) cut_thr_max[i]=cut_thr[i];
        rho_thr_max=rho_thr;
        nevt_max=nevt_new;

        if(gOFNAME!="") DumpPPCuts(gOFNAME, gNCUT, cut_mean, true);
      }
    }
  }
  if(gOFNAME!="") DumpPPCuts(gOFNAME, 0, NULL, true);	// dump blank line

  // CENTRAL TUNING

  // generate random threshold values around the cut_thr_max (in the min:max range) globaly for all cuts
  cout << endl << "Central tuning ..." << endl << endl;
  for(int i=0;i<gNCTRY;i++) {
    if(i%1000==0) cout << "Loop ->\t" << i << " / " << gNCTRY << endl;

    for(int i=0;i<gNCUT;i++) {
      cut_thr[i]=gRandom->Gaus(cut_mean[i],cut_sigma[i]);
      if(cut_thr[i]<gCUT_MIN[i]) cut_thr[i]=gCUT_MIN[i];
      if(cut_thr[i]>gCUT_MAX[i]) cut_thr[i]=gCUT_MAX[i];
    }

    int nevt_new = GetSelectedEvents(mtree_BKG, mtree_SIM, user_cut, gNCUT, cut_thr, rho_thr);
#ifdef CHECK_DNEVT_VS_IFAR
    bool check_dnevt = (nevt_new>nevt_max) ? CheckIncrementalDetectionsVsIFAR(tree_BKG, tree_SIM, user_cut, cut_thr) : true;
#else
    bool check_dnevt = true;
#endif
    if(nevt_new>nevt_max && check_dnevt==true) {

      PrintPPCuts(nevt_new, nevt_std, gNCUT, cut_thr, rho_thr, user_cut, false, true, false);

      // update cut_mean
      for(int i=0;i<gNCUT;i++) cut_mean[i]=cut_thr[i];

      // update max values
      for(int i=0;i<gNCUT;i++) cut_thr_max[i]=cut_thr[i];
      rho_thr_max=rho_thr;
      nevt_max=nevt_new;

      if(gOFNAME!="") DumpPPCuts(gOFNAME, gNCUT, cut_mean, true);
    }
  }

  // FINAL TUNING

  // generate random threshold values around the cut_thr_max (in the min:max range) individually for each cut
  cout << endl << "Final tuning ..." << endl << endl;
  if(gOFNAME!="") DumpPPCuts(gOFNAME, 0, NULL, true);	// dump blank line
  for(int i=0;i<gNCUT;i++) {
    if(gCUT_SIG[i]==0) continue;	// skip fixed cuts			
    cout << " " << i << " -> " << gCUT_NAME[i] << " ..." << endl;
    for(int k=0;k<gNCUT;k++) cut_thr[k]=cut_thr_max[k];
    for(int j=0;j<gNIFTRY;j++) {
      cut_thr[i]=gRandom->Uniform(gCUT_MIN[i],gCUT_MAX[i]);

      int nevt_new = GetSelectedEvents(mtree_BKG, mtree_SIM, user_cut, gNCUT, cut_thr, rho_thr);
#ifdef CHECK_DNEVT_VS_IFAR
      bool check_dnevt = (nevt_new>nevt_max) ? CheckIncrementalDetectionsVsIFAR(tree_BKG, tree_SIM, user_cut, cut_thr) : true;
#else
      bool check_dnevt = true;
#endif
      if(nevt_new>nevt_max && check_dnevt==true) {

	cout << endl << "  " << gCUT_NAME[i] << " : standard = " << gCUT_MEA[i] << " -> tuned = " << cut_thr[i] << endl;
        PrintPPCuts(nevt_new, nevt_std, gNCUT, cut_thr, rho_thr, user_cut, false, true, false);

        // update cut_mean
        for(int i=0;i<gNCUT;i++) cut_mean[i]=cut_thr[i];

        // update max values
        for(int i=0;i<gNCUT;i++) cut_thr_max[i]=cut_thr[i];
        rho_thr_max=rho_thr;
        nevt_max=nevt_new;

        if(gOFNAME!="") DumpPPCuts(gOFNAME, gNCUT, cut_mean, true);
      }
    }
  }
  if(gOFNAME!="") DumpPPCuts(gOFNAME, 0, NULL, true);	// dump blank line

  // if no improvements we use the input mean PP cut thresholds
  if(rho_thr_max==0.0) {
    rho_thr_max=rho_thr;
    for(int i=0;i<gNCUT;i++) cut_thr_max[i]=gCUT_MEA[i];
  }

  // ---------------------------------------------------------------
  // print detections/loudest for tuned PP cuts
  // ---------------------------------------------------------------

  cout << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "IFAR  = " << IFAR << " (years) " << "\tnBKG = " << gNBKG << endl;
  cout << "nCTRY = " << gNCTRY << "\tnIFTRY = " << gNIFTRY << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << endl;

  PrintLoudest(tree_BKG, bin_cut_standard.GetTitle(), 0, NULL, rho_thr, "standard", nBKG_MAX);		// loudest background with standard PP cuts
  PrintLoudest(tree_BKG, user_cut, gNCUT, cut_thr_max, rho_thr_max, "tuned", nBKG_MAX);			// loudest background with tuned PP cuts
  PrintPPCuts(nevt_max, nevt_std, gNCUT, cut_thr_max, rho_thr_max, user_cut, true, true, false);	// list tuned PP cuts

  // ---------------------------------------------------------------
  // Dump incremental detections: tuned/standard vs IFAR
  // ---------------------------------------------------------------

  DumpIncrementalDetectionsVsIFAR(tree_BKG, tree_SIM, user_cut, cut_thr_max, gOFNAME, true);

  exit(0);
}

int GetSelectedEvents(TTree* tree_BKG, TTree* tree_SIM, TString user_cut, int ncut, double* cut_thr, double& rho_thr) {

  TString sel="rho[1]";

  // make cut string
  TString cut="";
  if(user_cut!="") cut = TString::Format("(%s)",user_cut.Data());
  if(cut!="" && ncut>0) cut = cut + " && ";
  for(int i=0;i<ncut;i++) {
    TString sAND = i==ncut-1 ? "" : "&&"; 
    cut = cut + TString::Format("(%s%s%f)%s",gCUT_NAME[i].Data(),gCUT_CMP[i].Data(),cut_thr[i],sAND.Data());
  }

  TString cut_bkg = cut + TString::Format(" && !(%s)",gZL_CUT.GetTitle());
  tree_BKG->Draw(sel, cut_bkg, "goff");
  int nBKG = (Int_t)tree_BKG->GetSelectedRows();

  int ndet_SIM=0;
  if(nBKG>gNBKG) {

    // get rho[1] @ selected IFAR -> @ selected gNBKG 
    int *index = new Int_t[nBKG];
    double* value = tree_BKG->GetV1();
    TMath::Sort(nBKG,value,index,true);
    rho_thr = value[index[gNBKG]];
    delete [] index;
  
    cut = cut + TString::Format(" && rho[1]>%f",rho_thr);
    tree_SIM->Draw(sel,cut,"goff");
    ndet_SIM = (Int_t)tree_SIM->GetSelectedRows();
  }

  return ndet_SIM;
}

void PrintLoudest(TTree* tree_BKG, TString user_cut, int ncut, double* cut_thr, double rho_thr, TString label, int nLoudest) {

  // make cut string
  TString cut="";
  if(user_cut!="") cut = TString::Format("(%s)",user_cut.Data());
  if(cut!="" && ncut>0) cut = cut + " && ";
  for(int i=0;i<ncut;i++) {
    TString sAND = i==ncut-1 ? "" : "&&"; 
    cut = cut + TString::Format("(%s%s%f)%s",gCUT_NAME[i].Data(),gCUT_CMP[i].Data(),cut_thr[i],sAND.Data());
  }

  // print PP cuts string
  cout << endl << label << " PP cuts string ..." << endl;
  cout << endl << cut << endl << endl;

  // print Loudest background list
  cout << endl << label << " Loudest background list ..." << endl;
  TString cut_bkg = cut + TString::Format(" && !(%s)",gZL_CUT.GetTitle());
  tree_BKG->Draw("rho[1]",cut_bkg.Data(),"goff");
  int ndet_BKG = (Int_t)tree_BKG->GetSelectedRows();
  int *index = new Int_t[ndet_BKG];
  double* value = tree_BKG->GetV1();
  TMath::Sort(ndet_BKG,value,index,true);
  double rho_sup = value[index[0]];
  cout << endl;
  for(int i=0;i<nLoudest;i++) cout << "\t" << i << " -> \t" << value[index[i]] << endl;
  cout << endl;
  delete [] index;
}

void PrintPPCuts(int nevt_new, int nevt_old, int ncut, double* cut_thr, double rho_thr, TString user_cut, bool bppcuts, bool beff, bool bppcut) {

  // print PP cuts
  if(bppcuts) {
    char sout[1024]; 
    cout << endl << " PP cuts list ..." << endl << endl;
    // print header 
    sprintf(sout,"%*s %*s %*s %*s %*s %*s", 25, "name", 5, "cmp", 10, "tuned", 10, "mean", 14, "mean-tuned", 12, "increment");
    cout << endl << sout << endl << endl;
    for(int i=0;i<ncut;i++) {
      float dthr = 100*(cut_thr[i]-gCUT_MEA[i])/gCUT_MEA[i];
      sprintf(sout,"%*s %*s %*.2f %*.2f %*.2f %*.1f %s", 25, gCUT_NAME[i].Data(), 5, gCUT_CMP[i].Data(), 10, cut_thr[i], 10, gCUT_MEA[i], 14, cut_thr[i]-gCUT_MEA[i], 10, dthr, "%");
      cout << sout << endl;
    }
  }

  // print eff inprovement respect to the old PP cuts
  if(beff) {
    double dnevt = (nevt_new-nevt_old)/double(nevt_old);
    cout << endl;
    cout << "  @ rho[1] > " << rho_thr << "\t -> nevt_new = " << nevt_new << "\t -> dnevt = " << 100.*dnevt << " %" << endl;
    cout << endl;
  }

  // print PP cut string
  if(bppcut) {
    // make cut string
    TString cut="";
    if(user_cut!="") cut = TString::Format(" (%s)",user_cut.Data());
    if(cut!="" && ncut>0) cut = cut + " && ";
    for(int i=0;i<ncut;i++) {
      TString sAND = i==ncut-1 ? "" : "&&"; 
      cut = cut + TString::Format("(%s%s%f)%s",gCUT_NAME[i].Data(),gCUT_CMP[i].Data(),cut_thr[i],sAND.Data());
    }
    cout << endl << cut << endl << endl;
  }
}

void DumpPPCuts(TString ofname, int ncut, double* cut_thr, bool app) {

  FILE *fp;
  char mode[3];
  if(app) strcpy(mode, "a"); else strcpy(mode, "w");
  
  if((fp = fopen(ofname, mode)) == NULL ) {
     cout << " DumpPPCuts error: cannot open file " << ofname <<". \n";
     return;
  };

  char sout[1024];

  // dump header 
  sprintf(sout,"");
  for(int i=0;i<ncut;i++) sprintf(sout,"%s %*s",sout,15,gCUT_NAME[i].Data());
  //if(!app) cout << endl << sout << endl << endl;
  if(!app) fprintf(fp,"\n%s\n\n",sout);

  // dump values 
  sprintf(sout,"");
  for(int i=0;i<ncut;i++) sprintf(sout,"%s %*.4f",sout,15,cut_thr[i]); 
  //if(app) cout << sout << endl;
  if(app) fprintf(fp,"%s\n",sout);

  // dump blank line
  if(app && ncut==0) fprintf(fp,"\n");

  fclose(fp);
}

void DumpIncrementalDetectionsVsIFAR(TTree* tree_BKG, TTree* tree_SIM, TString user_cut, double* cut_thr, TString ofname, bool app) {

  FILE *fp=NULL;
  char mode[3];
  if(app) strcpy(mode, "a"); else strcpy(mode, "w");
  char sout[1024];

  if(ofname!="") {
    if((fp = fopen(ofname, mode)) == NULL) {
      cout << " DumpIncrementalDetectionsVsIFAR error: cannot open file " << ofname <<". \n";
      return;
    }
    if(fp!=NULL) fprintf(fp,"\n\nDump incremental detections: tuned/standard vs IFAR\n\n"); 
  }

  double ifar_min = 1;
  double ifar_max = BKG_LIVE_TIME/YEAR/4;

  double difar = (ifar_max-ifar_min)/nIFAR;

  cout << endl;
  cout << "Dump incremental detections: tuned/standard vs IFAR" << endl;
  cout << endl;
  double ifar = ifar_min;
  int gNBKG_save = gNBKG;
  for(int i=0;i<nIFAR;i++) {

    //double ifar = ifar_min + i*difar;
    ifar *= 2;

    gNBKG = BKG_LIVE_TIME/ifar/YEAR;

    double rho_thr_std, rho_thr_tuned;

    int nevt_std = GetSelectedEvents(tree_BKG, tree_SIM, bin_cut_standard.GetTitle(), 0, NULL, rho_thr_std);
    int nevt_tuned = GetSelectedEvents(tree_BKG, tree_SIM, user_cut, gNCUT, cut_thr, rho_thr_tuned);

    if(nevt_std!=0 && nevt_tuned!=0) {	// if (nevt_std==0 || nevt_tuned==0) the minimum rho threshold is not sufficient

      // print eff inprovement respect to the standard PP cuts
      double dnevt = (nevt_tuned-nevt_std)/double(nevt_std);
      sprintf(sout," @ ifar = %.2f (@ rho(std)=%.2f -> %d bkg)	 -> nevt(std/tuned) = %d/%d (@ rho=%.2f)	 -> dnevt = %.2f %s", \
              ifar, rho_thr_std, gNBKG, nevt_std, nevt_tuned, rho_thr_tuned, 100.*dnevt, "%");
      cout << sout << endl;
      if(fp!=NULL) fprintf(fp,"%s\n",sout);
    }
    if(gNBKG==1) break;
  }
  gNBKG = gNBKG_save;	// restore gNBKG
  cout << endl;

  if(fp!=NULL) {fprintf(fp,"\n");fclose(fp);} 
}

bool CheckIncrementalDetectionsVsIFAR(TTree* tree_BKG, TTree* tree_SIM, TString user_cut, double* cut_thr) {
//
// return false if not all detections (vs IFAR) dnevt are >=0
//

  double ifar_min = 1;
  double ifar_max = BKG_LIVE_TIME/YEAR/4;

  double difar = (ifar_max-ifar_min)/nIFAR;

  bool check_dnevt=true;
  double ifar = ifar_min;
  int gNBKG_save = gNBKG;
  for(int i=0;i<nIFAR;i++) {

    //double ifar = ifar_min + i*difar;
    ifar *= 2;

    gNBKG = BKG_LIVE_TIME/ifar/YEAR;

    double rho_thr_std, rho_thr_tuned;

    int nevt_std = GetSelectedEvents(tree_BKG, tree_SIM, bin_cut_standard.GetTitle(), 0, NULL, rho_thr_std);
    int nevt_tuned = GetSelectedEvents(tree_BKG, tree_SIM, user_cut, gNCUT, cut_thr, rho_thr_tuned);

    if(nevt_std!=0 && nevt_tuned!=0) {	// if (nevt_std==0 || nevt_tuned==0) the minimum rho threshold is not sufficient
      // check eff inprovement respect to the standard PP cuts
      double dnevt = (nevt_tuned-nevt_std)/double(nevt_std);
      if(dnevt<0) check_dnevt=false;
    }
    if(gNBKG==1 || check_dnevt==false) break;
  }
  gNBKG = gNBKG_save;	// restore gNBKG

  return check_dnevt;
}

int ReadCutList(std::stringstream* in, TString* name, TString* cmp, double* mean, double* sigma, double* min, double* max) {

  // Read Cuts file list

  int size=0;
  char str[1024];
  int fpos=0;
  while(true) {
    in->getline(str,1024);
    if (!in->good()) break;
    if(str[0] != '#' && str[0] != ' ') size++;
  }
  in->clear(ios::goodbit);
  in->seekg(0, ios::beg);
  if (size>nCUT_MAX) {cout << "ReadCutList Error - PP cuts > " << nCUT_MAX << endl;exit(1);}

  char sout[1024]; 
  char xname[256], xcmp[256];
  double xmean, xsigma, xmin, xmax;
  int xenabled; 
  int ncut=0;
  cout << endl << " Input cut list ..." << endl;
  // print header 
  sprintf(sout,"%*s %*s %*s %*s %*s %*s %*s %*s", 3, "id" , 25, "name", 10, "cmp", 10, "mean", 10, "sigma", 10, "min", 10, "max", 15, "tune-enabled");
  cout << endl << sout << endl << endl;
  while(true) {
    fpos=in->tellg();
    in->getline(str,1024);
    if (!in->good()) break;
    in->seekg(fpos, ios::beg);
    *in >> xname >> xcmp >> xmean >> xsigma >> xmin >> xmax >> xenabled;
    if(in->good() && !(xname[0] == '#' || xname[0]=='\0' || xname[0] == ' ')) {
      name[ncut] = xname;
      cmp[ncut]  = xcmp;
      mean[ncut] = xmean;
      sigma[ncut]= xenabled==0 ? 0 : xsigma;      // if xenabled=0 than sigma=0 and parameter is fixed to xmean;
      min[ncut]  = xmin;
      max[ncut]  = xmax;
      sprintf(sout,"%*d %*s %*s %*.2f %*.2f %*.2f %*.2f %*d", 3, ncut+1 , 25, xname, 10, xcmp, 10, xmean, 10, xsigma, 10, xmin, 10, xmax, 10, xenabled);
      cout << sout << endl;
      ncut++;
    } else in->clear();
  }
  cout << endl;

  return ncut;
}

void PrintCovariance(TTree* tree, bool norm) {

  int entries = tree->GetEntries();
  double* value[nCUT_MAX];

  for(int i=0;i<gNCUT;i++) {
    cout << i << " " << gCUT_NAME[i].Data() << endl;
    tree->Draw(gCUT_NAME[i], "", "goff");
    value[i] = new double[entries];
    double* val = tree->GetV1();
    for(int j=0;j<entries;j++) value[i][j]=val[j];
  }

  TPrincipal pca(gNCUT, "ND");
  double data[nCUT_MAX];

  for(int i=0;i<entries;i++) {
    for(int j=0;j<gNCUT;j++) data[j]=value[j][i];
    pca.AddRow(data);
  }

  // get covariance matrix
  TMatrixD* matrix = (TMatrixD*)pca.GetCovarianceMatrix();
  //matrix->Print("%11.3g");

  // do normalization
  TMatrixD nmatrix = *matrix;
  for(int i=0;i<matrix->GetNrows();i++) {
    for(int j=0;j<matrix->GetNcols();j++) {
      nmatrix(i,j)/=sqrt((*matrix)(i,i) * (*matrix)(j,j));
      nmatrix(i,j)*=100.;
    }
  }
  //nmatrix.Print();

  if(norm) {
    *matrix=nmatrix; 
    matrix->Print("f= %11.0f");
  } else {
    matrix->Print("f= %11.3g");
  }

  cout << endl;
  cout << " list of loudest correlation (>30%) ... " << endl;
  cout << endl;
  char sout[1024]; 
  for(int i=0;i<matrix->GetNrows();i++) {
    for(int j=0;j<matrix->GetNcols();j++) {
      if(fabs((*matrix)(i,j))>30 && i!=j && gCUT_NAME[i]!=gCUT_NAME[j]) {
        sprintf(sout," %*s vs %*s %*s = %*.2f %s", 35, gCUT_NAME[i].Data(), 15, gCUT_NAME[j].Data(), 15, " -> correlation", 7, (*matrix)(i,j), "%");
        cout << sout << endl;
        //cout << " " << gCUT_NAME[i] << "/" << gCUT_NAME[j] << "\t\tCorrelation = " << matrix(i,j) << " %" << endl;
      }
    }
  }
  cout << endl;

  for(int i=0;i<gNCUT;i++) delete [] value[i];
}

