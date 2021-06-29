#define XIFO 4

#pragma GCC system_header

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TMath.h"
#include "TCut.h"
#include <fstream>
#include <vector>
#include "Toolbox.hh"
#include "wavearray.hh"

// used to generate "far vr rho" file 
void 
MakeFARvsRHO(TString bin, TString wave_file, double bkg_livetime, TString ofName, bool ocolumns2=true) {

  #include "../../cwb/postproduction/O1/GW150914_search_bins.hh"

  bin.ToUpper();
  if(bin!="UNMODELED" && bin!="CONSTRAINED" && bin!="CHIRP" && bin!="EUNMODELED" && bin!="ECONSTRAINED") {
    cout << "MakeFAR - Error : bin not defined !!!" << endl;
    gSystem->Exit(1);
  }

  if(bkg_livetime<=0) {
    cout << "MakeFAR - Error : bkg_livetime <= 0 !!!" << endl;
    gSystem->Exit(1);
  }

  CWB::Toolbox::checkFile(wave_file);

  TString selection = "";
  if(bin=="UNMODELED")    selection = unmodeled.GetTitle();
  if(bin=="CONSTRAINED")  selection = constrained.GetTitle();
  if(bin=="CHIRP")        selection = chirp.GetTitle();
  if(bin=="EUNMODELED")   selection = eunmodeled.GetTitle();
  if(bin=="ECONSTRAINED") selection = econstrained.GetTitle();

  TFile *wfile = TFile::Open(wave_file.Data());
  TTree* twave = (TTree *) gROOT->FindObject("waveburst");

  char sel[1024];
  sprintf(sel,selection.Data());
  sprintf(sel,"%s && (lag[2]!=0 || slag[2]!=0)", sel);
  cout << sel << endl; 

  twave->Draw("rho[0]",sel,"goff");
  double* trho = twave->GetV1();
  int N = twave->GetSelectedRows();
  cout << "N " << N << endl;

  int pp_rho_min = 7;
  int pp_rho_max = 500;
  float pp_rho_wbin = 0.01;
  int pp_rho_bin = float(pp_rho_max-pp_rho_min)/pp_rho_wbin;
  double drho = double(pp_rho_max-pp_rho_min)/double(pp_rho_bin);

  wavearray<double> Xcut(pp_rho_bin);
  wavearray<double> Ycut(pp_rho_bin); Ycut = 0.;
  for(int j=0; j<Xcut.size(); j++) {
    Xcut.data[j] = pp_rho_min+j*drho;
    for(int i=0; i<N; i++) {
      if(trho[i] > Xcut.data[j]) {Ycut.data[j] += 1.;}
    }
  }

  vector<double> rho,far;
  for(int j=Xcut.size()-1; j>=0; j--) {
    if (Ycut[j]<=0) continue;
    rho.push_back(Xcut[j]);
    far.push_back(Ycut[j]);
  }  

  // write output file
  ofstream outa;
  outa.open(ofName.Data(),ios::out);
  if(!outa.good()) {cout << "Error Opening File " << ofName << endl;exit(1);}
  //outa << "RHO" << "\t\t" << "1/years" << "\t\t" << "e_RHO" << "\t" << "e_1/years" << endl << endl;
  for(int i=0;i<Xcut.size(); i++) if (Ycut[i]>0) {
    double x  = Xcut[i];			
    double ex = 0;
    double y  = Ycut[i]/bkg_livetime;
    double ey = sqrt(Ycut[i])/bkg_livetime;
    if(ocolumns2) {	// 2 columns
      if(y>=1) outa << x << "\t\t" << y << endl;
      if(y<1)  outa << x << "\t\t" << y << endl;
    } else { 		// 4 columns
      //cout << i << " " << x << " " << y << " " << ex << " " << ey << endl;
      if(y>=1) outa << x << "\t\t" << y << "\t\t" << ex << "\t" << ey << endl;
      if(y<1)  outa << x << "\t\t" << y << "\t" << ex << "\t" << ey << endl;
    }
  }
  outa.close();

  return;
}

