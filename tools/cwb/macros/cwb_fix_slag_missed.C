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


{
  // this fix must be used if live tree has been produced without slag leaf

  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  TString cwb_merge_label;
  if(gSystem->Getenv("CWB_MERGE_LABEL")==NULL) {
    cout << "Error : environment CWB_MERGE_LABEL is not defined!!!" << endl;exit(1);
  } else {
    cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  }

  char net_file_name[256];
  char liv_file_name[256];

  if(simulation) {
    cout << "Error : cwb_fix_slag_missed.C must be applied only to production files !!!" << endl;exit(1);
  } else {
    sprintf(net_file_name,"%s/wave_%s.%s.root",merge_dir,data_label,cwb_merge_label.Data());
    sprintf(liv_file_name,"%s/live_%s.%s.root",merge_dir,data_label,cwb_merge_label.Data());
    cout << net_file_name << endl;
    cout << liv_file_name << endl;

    // Check if files exist
    TB.checkFile(net_file_name);
    TB.checkFile(liv_file_name);
  }

  // open wave tree
  TFile *fwave = TFile::Open(net_file_name);
  if(fwave==NULL) {cout << "Error opening file " << net_file_name << endl;exit(1);}
  TTree* twave = (TTree *) gROOT->FindObject("waveburst");    
  if(twave==NULL) {cout << "Error opening tree wave" << endl;exit(1);}

  // check if slag is presents in wave tree
  TBranch* branch;
  bool net_slag=false;
  TIter next(twave->GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if (TString("slag").CompareTo(branch->GetName())==0) net_slag=true;
  }
  next.Reset();

  // open live tree
  TFile *flive = TFile::Open(liv_file_name);
  if(flive==NULL) {cout << "Error opening file " << liv_file_name << endl;exit(1);}
  TTree* tlive = (TTree *) gROOT->FindObject("liveTime");    
  if(tlive==NULL) {cout << "Error opening tree livetime" << endl;exit(1);}

  // check if slag is presents in live tree
  TBranch* branch;
  bool liv_slag=false;
  TIter next(tlive->GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if (TString("slag").CompareTo(branch->GetName())==0) liv_slag=true;
  }
  next.Reset();
  if(liv_slag&&net_slag) {
    cout << "Files " << endl;
    cout << liv_file_name << endl;
    cout << net_file_name << endl;
    cout << "already contains slag leaf, exit" << endl;
    exit(1);
  }

  if(!liv_slag&&net_slag) {
    cout << "File " << net_file_name << " aleady contains slag leaf" << endl;
    cout << "only " << liv_file_name << " needs to be fixed" << endl;
    cout << "apply cwb_fix_live_slag_missed.C " << endl;
    exit(1);
  }

  if(liv_slag&&!net_slag) {
    cout << "File " << liv_file_name << " aleady contains slag leaf" << endl;
    cout << "only " << net_file_name << " needs to be fixed" << endl;
    cout << "something wrong during production, please check" << endl;
    exit(1);
  }

  // insert slag=0 in wave and live trees 
  // tree has been produced without slag
  bool fixit=false;
  if (!liv_slag) {
    char answer[256];
    strcpy(answer,"");
    do {
      cout << "File " << net_file_name << " do not contains slag leaf" << endl;
      cout << "File " << liv_file_name << " do not contains slag leaf" << endl;
      cout << "Do you want to fix them ? (y/n) ";
      cin >> answer;
      cout << endl << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"y")==0) fixit=true;
  }
  if(!fixit) exit(0);

  // -------------------------------------------------------------------------
  // create slag fix live root file
  // -------------------------------------------------------------------------

  TString fix_liv_file_name(liv_file_name);
  fix_liv_file_name.ReplaceAll(".root",".FIX.root");
  cout << "Fixed live file : " << fix_liv_file_name.Data() << endl;

  TFile* flive_fix = new TFile(fix_liv_file_name,"RECREATE");
  if (flive_fix->IsZombie()) {
    cout << "CWB::Toolbox::setVeto - Error opening file " << fix_liv_file_name.Data() << endl;
    exit(1);
  }
  TTree *tlive_fix = (TTree*)tlive->CloneTree(0);
  tlive_fix->SetMaxTreeSize(5000000000);
  // add slag leaf
  float  xslag[6];
  tlive_fix->Branch("slag",xslag,"slag[6]/F");
  tlive_fix->Write();

  Int_t  xrun;
  tlive->SetBranchAddress("run",&xrun);
  int ntrg_live = tlive->GetEntries();
  for(int i=0; i<ntrg_live; i++) {
    if(i%100000==0) cout << i << "/" << ntrg_live << endl;
    tlive->GetEntry(i);
    for (int n=0;n<nIFO;n++) xslag[n] = 0.;
    xslag[nIFO] = 0.;
    for (int n=nIFO+1;n<6;n++) xslag[n] = -1;
    tlive_fix->Fill();
  }

  tlive_fix->Write();
  flive_fix->Close();


  // -------------------------------------------------------------------------
  // create slag fix wave root file
  // -------------------------------------------------------------------------

  TString fix_net_file_name(net_file_name);
  fix_net_file_name.ReplaceAll(".root",".FIX.root");
  cout << "Fixed wave file : " << fix_net_file_name.Data() << endl;

  TFile* fwave_fix = new TFile(fix_net_file_name,"RECREATE");
  if (fwave_fix->IsZombie()) {
    cout << "CWB::Toolbox::setVeto - Error opening file " << fix_net_file_name.Data() << endl;
    exit(1);
  }
  TTree *twave_fix = (TTree*)twave->CloneTree(0);
  twave_fix->SetMaxTreeSize(5000000000);
  // add slag leaf
  float  xslag[6];
  twave_fix->Branch("slag",xslag,"slag[6]/F");
  twave_fix->Write();

  int ntrg_wave = twave->GetEntries();
  for(int i=0; i<ntrg_wave; i++) {
    if(i%100000==0) cout << i << "/" << ntrg_wave << endl;
    twave->GetEntry(i);
    for (int n=0;n<nIFO;n++) xslag[n] = 0.;
    xslag[nIFO] = 0.;
    for (int n=nIFO+1;n<6;n++) xslag[n] = -1;
    twave_fix->Fill();
  }

  twave_fix->Write();
  fwave_fix->Close();

  // -------------------------------------------------------------------------
  // create a symbolic link to merge*.lst file name
  // -------------------------------------------------------------------------
  vector<int> jobList;
  char ilstfname[1024];
  sprintf(ilstfname,"merge_%s.%s.lst",data_label,cwb_merge_label.Data());
  TString olstfname = fix_net_file_name;
  olstfname.ReplaceAll("wave_","merge_");
  olstfname.ReplaceAll(".root",".lst");
  olstfname.Remove(0,olstfname.Last('/')+1);    // strip path
  cout << olstfname << endl;
  int estat;
  Long_t id,size,flags,mt;
  estat = gSystem->GetPathInfo(TString(merge_dir)+"/"+ilstfname,&id,&size,&flags,&mt);
  if (estat==0) {
    char cmd[1024];
    sprintf(cmd,"cd %s;ln -sf %s %s",merge_dir,ilstfname,olstfname.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  exit(0);
}
