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

// This macro is used to add to the wave root file the redshift value reported in the mdc file
//
// how to use 
//            Ex: root -l -b 'AddRedshift2WAVE.C("merge/wave_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.root")'
//
// outputs:
//            new output wave file with leaf redshift
//              merge/wave_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.U_redshift.root
//
//            symbolic link to the merge/mdc_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.root
//              merge/mdc_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.U_redshift.root
//
//            symbolic link to the merge/merge_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.lst
//              merge/merge_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.U_redshift.lst
//


#define MAX_TREE_SIZE	100000000000LL

int binarySearch(double array[], int start, int end, double key);
int binarySearch(double array[], int size, double key);

void AddRedshift2WAVE(TString ifname) {
//
//  ifname	= is the name of the input wave root simulation file 
//		  Ex: "merge/wave_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.root"
//

  // ---------------------------------------------------------------
  // Open input WAVE ROOT file
  // ---------------------------------------------------------------

  TString iwfile = ifname;
  cout << endl << " " << iwfile << endl;

  TFile *iwroot = TFile::Open(iwfile);
  if(iwroot==NULL) {cout << "Error opening file: " << iwfile << endl;exit(1);}

  TTree* iwtree = (TTree *) gROOT->FindObject("waveburst");
  if(iwtree==NULL) {cout << "waveburst tree not found !!!" << endl;exit(1);}
  int wentries = (int)iwtree->GetEntries();
  cout << endl << " wntries = " << wentries << endl;

  // get number of detectors in the network
  TList* ifoList = iwtree->GetUserInfo();
  int nIFO = ifoList->GetSize();        
  cout << endl << " nIFO = " << nIFO << endl << endl;

  // ---------------------------------------------------------------
  // Open input MDC ROOT file
  // ---------------------------------------------------------------

  TString imfile = ifname;
  imfile.ReplaceAll("wave_","mdc_");
  cout << " " << imfile << endl;

  TFile *imroot = TFile::Open(imfile);
  if(imroot==NULL) {cout << "Error opening file: " << imfile << endl;exit(1);}

  TTree* imtree = (TTree *) gROOT->FindObject("mdc");
  if(imtree==NULL) {cout << "mdc tree not found !!!" << endl;exit(1);}

  // sort mdc injection wrt time of the first detector
  imtree->Draw("time[0]:redshift","","goff");
  int imentries = (int)imtree->GetSelectedRows();
  int *mindex = new int[imentries];
  double* imtime = imtree->GetV1(); 
  double* imredshift = imtree->GetV2();
  TMath::Sort(imentries,imtime,mindex,false);

  // create sorted mdc time array
  double* stime = new double[imentries];
  for (Int_t i=0;i<imentries;i++) stime[i] = imtime[mindex[i]];

  // ---------------------------------------------------------------
  // create updated output WAVE ROOT file with new redshift leaf
  // ---------------------------------------------------------------

  TString owfile = ifname;
  owfile.ReplaceAll(".root",".U_redshift.root");
  cout << endl << " new output tree file = " << owfile << endl << endl;

  //open new file to store the updated tree
  TFile owroot(owfile,"recreate");
  //Create an empty clone of the original input iwtree
  TTree *owtree = (TTree*)iwtree->CloneTree(0);
  owtree->SetMaxTreeSize(MAX_TREE_SIZE);

  // add redshift leaf to owtree 
  TBranch* branch;
  bool redshift_exists=false;
  TIter next(iwtree->GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if(TString("redshift").CompareTo(branch->GetName())==0) redshift_exists=true;
  }
  next.Reset();
  float wredshift;
  if (redshift_exists) iwtree->SetBranchAddress("redshift",&wredshift);
  else                 owtree->Branch("redshift",&wredshift,"redshift/F");

  double* iwtime = new double[2*nIFO];
  iwtree->SetBranchAddress("time",iwtime);

  cout << endl;
  for (Int_t i=0;i<wentries;i++) {
     if (i%1000==0 && i>0) cout << " write entry : " << i << endl;
     iwtree->GetEntry(i);
     int sindex = binarySearch(stime, 0, imentries, iwtime[nIFO]);
     wredshift = imredshift[mindex[sindex]];
     owtree->Fill();
  }
  owtree->Write();
  delete [] mindex;

  // ---------------------------------------------------------------
  // create symbolic link to input mdc root file & merge lst
  // ---------------------------------------------------------------

  TString cmd;

  TString omfile = ifname;
  omfile.ReplaceAll("wave_","mdc_");
  omfile.ReplaceAll(".root",".U_redshift.root");
  TString imfile_base = gSystem->BaseName(imfile.Data());	// base name
  cmd = TString::Format("ln -sf %s %s",imfile_base.Data(),omfile.Data());
  cout << endl << " " << cmd << endl << endl;
  gSystem->Exec(cmd.Data());

  TString ilfile = ifname;
  ilfile.ReplaceAll("wave_","merge_");
  ilfile.ReplaceAll(".root",".lst");
  TString olfile = ifname;
  olfile.ReplaceAll("wave_","merge_");
  olfile.ReplaceAll(".root",".U_redshift.lst");
  TString ilfile_base = gSystem->BaseName(ilfile.Data());	// base name
  cmd = TString::Format("ln -sf %s %s",ilfile_base.Data(),olfile.Data());
  cout << endl << " " << cmd << endl << endl;
  gSystem->Exec(cmd.Data());

  exit(0);
}

//______________________________________________________________________________
int binarySearch(double array[], int start, int end, double key) {
   // Determine the search point.
   // int searchPos = end + ((start - end) >> 1);
   // int searchPos = (start + end) / 2;        
   int searchPos = (start + end) >> 1;
   // If we crossed over our bounds or met in the middle, then it is not here.
   if (start > end)
      return -1;
   // Search the bottom half of the array if the query is smaller.
   if (array[searchPos] > key)
      return binarySearch (array, start, searchPos - 1, key);
   // Search the top half of the array if the query is larger.
   if (array[searchPos] < key)
      return binarySearch (array, searchPos + 1, end, key);
   // If we found it then we are done.
   if (array[searchPos] == key)
      return searchPos;

   return -1;
}

//______________________________________________________________________________
int binarySearch(double array[], int size, double key) {
   return binarySearch(array, 0, size - 1, key);
}

