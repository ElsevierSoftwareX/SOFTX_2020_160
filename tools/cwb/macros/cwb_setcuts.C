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


// apply selection on leaves to merged root file : used by the cwb_setcuts command

{

  CWB::Toolbox TB;

  int estat;
  Long_t id,size,flags,mt;
  char cmd[1024];

  	  cwb_merge_label     = TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  TString cwb_jcuts_tree      = TString(gSystem->Getenv("CWB_JCUTS_TREE"));
  TString cwb_wcuts_tree      = TString(gSystem->Getenv("CWB_WCUTS_TREE"));
  TString cwb_cuts_label      = TString(gSystem->Getenv("CWB_CUTS_LABEL"));
  TString cwb_mcuts_tree      = TString(gSystem->Getenv("CWB_MCUTS_TREE"));
  TString cwb_lcuts_tree      = TString(gSystem->Getenv("CWB_LCUTS_TREE"));
  TString cwb_setcuts_options = TString(gSystem->Getenv("CWB_SETCUTS_OPTIONS"));

  TString mdir = merge_dir;

  // check if cwb_setcuts_options contains '--' than cwb_setcuts_options
  // is used to extract all setcuts parameters
  TString cwb_tcuts_tree = "";
  TString cwb_unique_evt = "false";
  if(cwb_setcuts_options.Contains("--")) {
    TString option="";
    // get the cwb_tcuts_tree
    cwb_tcuts_tree = TB.getParameter(cwb_setcuts_options,"--tcuts");
    // get the cwb_wcuts_tree
    cwb_wcuts_tree = TB.getParameter(cwb_setcuts_options,"--wcuts");
    // get the cwb_jcuts_tree
    cwb_jcuts_tree = TB.getParameter(cwb_setcuts_options,"--jcuts");
    // get the cwb_mcuts_tree
    cwb_mcuts_tree = TB.getParameter(cwb_setcuts_options,"--mcuts");
    // get the cwb_lcuts_tree
    cwb_lcuts_tree = TB.getParameter(cwb_setcuts_options,"--lcuts");
    // get the cwb_cuts_label
    cwb_cuts_label = TB.getParameter(cwb_setcuts_options,"--label");
    // get the cwb_unique_evt
    cwb_unique_evt = TB.getParameter(cwb_setcuts_options,"--unique");
  }
  // if cwb_tcuts_tree!="" than cwb_tcuts_tree is added to cwb_wcuts_tree
  if(cwb_tcuts_tree!="") {
     // check if tcut is defined
     TGlobal *global=(TGlobal*)gROOT->GetListOfGlobals()->FindObject(cwb_tcuts_tree.Data());
     if(global==NULL) {
        cout << "cwb_setcuts.C : Error - tcut " << cwb_tcuts_tree << " is not defined " << endl; 
        cout << "                must be included in the user_pparameters.C file" << endl;
        gSystem->Exit(1);
     }
     TCut *tcut = (TCut*)global->GetAddress();
     if(TString(tcut->GetName())!=cwb_tcuts_tree) {
        cout << "cwb_setcuts.C : Error - tcut " << cwb_tcuts_tree << " not correspond to global Tcut name" << endl; 
        cout << "                must be included in the user_pparameters.C file" << endl;
        gSystem->Exit(1);
     }
     cwb_tcuts_tree = tcut->GetTitle();	// extract TCut 
     cout << "cwb_tcuts_tree : " << cwb_tcuts_tree << endl;
     if(cwb_wcuts_tree=="") {
        cwb_wcuts_tree = cwb_tcuts_tree;
     } else {
       TString cwb_wcuts = cwb_tcuts_tree;
       cwb_wcuts = TString("(")+cwb_wcuts_tree+")&&("+cwb_tcuts_tree+")";
       cwb_wcuts_tree = cwb_wcuts;
     }
     cout << "cwb_wcuts_tree : " << cwb_wcuts_tree << endl;
  }
  // check if cwb_unique_evt can be applied 
  if((simulation==0)&&(cwb_unique_evt=="true")) {
    cout << "cwb_setcuts.C : Error - unique option can be applied only to simulated data" << endl; 
    gSystem->Exit(1);
  }
  // check if cwb_cuts_label is defined
  if((cwb_wcuts_tree!="")&&(cwb_cuts_label=="")) {
    cout << "cwb_setcuts.C : Error - cuts label not defined" << endl; 
    gSystem->Exit(1);
  }
  // check if cwb_cuts_label contains '.'
  if(cwb_cuts_label.Contains(".")) {
    cout << "cwb_setcuts.C : Error - cwb_cuts_label " << cwb_cuts_label 
         << " can not contains '.'" << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_wcuts_tree contains run leaf
  if(cwb_wcuts_tree.Contains("run")) {
    cout << "cwb_setcuts.C : Error - cwb_wcuts_tree " << cwb_jcuts_tree 
         << " can not contains leaf : run" << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_wcuts_tree contains slag leaf
  if(cwb_wcuts_tree.Contains("slag")) {
    cout << "cwb_setcuts.C : Error - cwb_wcuts_tree " << cwb_jcuts_tree 
         << " can not contains leaf : slag" << endl << endl;
    gSystem->Exit(1);
  }
  // check if cwb_wcuts_tree contains lag leaf
  if(cwb_wcuts_tree.Contains("lag")) {
    cout << "cwb_setcuts.C : Error - cwb_wcuts_tree " << cwb_jcuts_tree 
         << " can not contains leaf : lag" << endl << endl;
    gSystem->Exit(1);
  }
  if(cwb_wcuts_tree=="") cwb_wcuts_tree="run>0";
  if(cwb_jcuts_tree=="") cwb_jcuts_tree="run>0";

  // create wcuts
  TString cwb_wcuts = cwb_wcuts_tree;
  cwb_wcuts = TString("(")+cwb_wcuts_tree+")&&("+cwb_jcuts_tree+")";
  if(cwb_lcuts_tree!="") cwb_wcuts = TString("(")+cwb_wcuts+")&&("+cwb_lcuts_tree+")";

  // create mcuts
  TString cwb_mcuts = cwb_mcuts_tree;
  if(cwb_mcuts!="") 
    cwb_mcuts = TString("(")+cwb_mcuts_tree+")&&("+cwb_jcuts_tree+")";
  else
    cwb_mcuts = cwb_jcuts_tree;

  // create input wave root cuts file name
  char iwfname[1024];  
  sprintf(iwfname,"wave_%s.%s.root",data_label,cwb_merge_label.Data());

  // create output wave root cuts file name
  TString uwfname = mdir+"/"+iwfname;
  uwfname.ReplaceAll(".root",TString(".C_U.root"));

  // if cwb_unique_evt=true then create a wave root file with unique
  // reconstructed events
  if(cwb_unique_evt=="true") {
    TB.setUniqueEvents(mdir+"/"+iwfname, uwfname, nIFO, pp_irho);
    // create symbolic link of mdc file
    char imfname[1024];  
    sprintf(imfname,"mdc_%s.%s.root",data_label,cwb_merge_label.Data());
    TString omfname = uwfname;
    omfname.ReplaceAll("wave_","mdc_");
    omfname.Remove(0,omfname.Last('/')+1);	// strip path
    cout << omfname << endl;
    estat = gSystem->GetPathInfo(mdir+"/"+imfname,&id,&size,&flags,&mt);
    if(estat==0) {
      sprintf(cmd,"cd %s;ln -sf %s %s",mdir.Data(),imfname,omfname.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
    }
    // create symbolic link of merge list file
    char ilstfname[1024];  
    sprintf(ilstfname,"merge_%s.%s.lst",data_label,cwb_merge_label.Data());
    TString olstfname = uwfname;
    olstfname.ReplaceAll("wave_","merge_");
    olstfname.ReplaceAll(".root",".lst");
    olstfname.Remove(0,olstfname.Last('/')+1);	// strip path
    cout << olstfname << endl;
    estat = gSystem->GetPathInfo(mdir+"/"+ilstfname,&id,&size,&flags,&mt);
    if (estat==0) {
      sprintf(cmd,"cd %s;ln -sf %s %s",mdir.Data(),ilstfname,olstfname.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
    }
  }
  if((cwb_wcuts_tree="")&&(cwb_cuts_label=="")) exit(0);

  if(cwb_unique_evt=="true") {
    sprintf(iwfname,"wave_%s.%s.C_U.root",data_label,cwb_merge_label.Data());
    cwb_merge_label+=".C_U";
  }

  // create output wave root cuts file name
  TString owfname = mdir+"/"+iwfname;
  owfname.ReplaceAll(".root",TString(".C_")+cwb_cuts_label+".root");

  // create a merge*.lst file name & run selection cuts
  vector<int> jobList;
  char ilstfname[1024];  
  sprintf(ilstfname,"merge_%s.%s.lst",data_label,cwb_merge_label.Data());
  TString olstfname = owfname;
  olstfname.ReplaceAll("wave_","merge_");
  olstfname.ReplaceAll(".root",".lst");
  olstfname.Remove(0,olstfname.Last('/')+1);	// strip path
  cout << olstfname << endl;
  if(cwb_jcuts_tree=="") {	// no jobs cuts -> create a symbolic link
    estat = gSystem->GetPathInfo(mdir+"/"+ilstfname,&id,&size,&flags,&mt);
    if (estat==0) {
      sprintf(cmd,"cd %s;ln -sf %s %s",mdir.Data(),ilstfname,olstfname.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
    }
  } else {
    // build list of selected jobs
    vector<TString> jobFileList;
    vector<int> jobList = CWB::Toolbox::getMergeJobList(mdir+"/"+ilstfname,jobFileList);
    cout << "Number of jobs in the merged lists : " << jobList.size() << endl;
    // create a dummy tree, fill with merged jobs, apply cuts, get selected jobs
    TTree jtree("jtree","jtree");
    int run; jtree.Branch("run",&run,"run/I");
    TString *jfname = new TString();
    jtree.Branch("jfname", &jfname);
    // fill with merged jobs
    for(int i=0;i<jobList.size();i++) {run=jobList[i]; *jfname=jobFileList[i]; jtree.Fill();}
    // check cuts
    TTreeFormula formula("trCuts", cwb_jcuts_tree.Data(), &jtree);
    int err = formula.Compile(cwb_jcuts_tree.Data());
    if(err) {
      cout << "cwb_setcuts.C - wrong input cuts " << cwb_jcuts_tree << endl << endl;
      gSystem->Exit(1);
    }
    // get selected jobs 
    jtree.Draw("Entry$",cwb_jcuts_tree.Data(),"goff");
    int jsel = jtree.GetSelectedRows();
    double* entry = jtree.GetV1();			// get selected jobs
    jobList.clear();
    jobFileList.clear();
    for(int i=0;i<jsel;i++) {
      jtree.GetEntry(entry[i]);
      //cout << i << " " << run << " " << jfname->Data() << endl;
      jobList.push_back(run);				// list of selected jobs
      jobFileList.push_back(jfname->Data());		// list of selected file jobs
    }
    CWB::Toolbox::dumpFileList(jobFileList,mdir+"/"+olstfname);
  }

  // apply wcuts to wave file
  int nsel = CWB::Toolbox::setCuts(iwfname,mdir,mdir,"waveburst",cwb_wcuts,cwb_cuts_label);
  if(nsel<=0) {
    cout << "cwb_setcuts.C : Error - Number of selected waveburst entries = " << nsel << endl << endl;
    gSystem->Exit(1);
  }
  cout << "cwb_setcuts.C : Number of selected waveburst entries = " << nsel << endl;

  // apply mcuts to mdc file
  if(simulation) {
    char imfname[1024];  
    sprintf(imfname,"mdc_%s.%s.root",data_label,cwb_merge_label.Data());
    int nsel = CWB::Toolbox::setCuts(imfname,mdir,mdir,"mdc",cwb_mcuts,cwb_cuts_label);
    cout << "cwb_setcuts.C : Number of selected mdc entries = " << nsel << endl << endl;
    if(nsel<0) { // create a symbolic link to mdc*.root file name
      TString omfname = owfname;
      omfname.ReplaceAll("wave_","mdc_");
      omfname.Remove(0,omfname.Last('/')+1);	// strip path
      cout << omfname << endl;
      estat = gSystem->GetPathInfo(mdir+"/"+imfname,&id,&size,&flags,&mt);
      if(estat==0) {
        sprintf(cmd,"cd %s;ln -sf %s %s",mdir.Data(),imfname,omfname.Data());
        cout << cmd << endl;
        gSystem->Exec(cmd);
      }
    }
  }

  // apply cuts to live file
  if(!simulation) {
    char ilfname[1024];  
    sprintf(ilfname,"live_%s.%s.root",data_label,cwb_merge_label.Data());
    TString cwb_live_cuts = cwb_jcuts_tree;
    if(cwb_lcuts_tree!="") cwb_live_cuts = TString("(")+cwb_live_cuts+")&&("+cwb_lcuts_tree+")";
    int nsel = -1;
    if(cwb_live_cuts!="run>0") {
      nsel = CWB::Toolbox::setCuts(ilfname,mdir,mdir,"liveTime",cwb_live_cuts,cwb_cuts_label);
      cout << "cwb_setcuts.C : Number of selected live entries = " << nsel << endl;
    }
    if(nsel<0) { // create a symbolic link to live*.root file name
      TString olfname = owfname;
      olfname.ReplaceAll("wave_","live_");
      olfname.Remove(0,olfname.Last('/')+1);	// strip path
      cout << olfname << endl;
      estat = gSystem->GetPathInfo(mdir+"/"+ilfname,&id,&size,&flags,&mt);
      if(estat==0) {
        sprintf(cmd,"cd %s;ln -sf %s %s",mdir.Data(),ilfname,olfname.Data());
        cout << cmd << endl;
        gSystem->Exec(cmd);
      }
    }
  }

  exit(0);
}
