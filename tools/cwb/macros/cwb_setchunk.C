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


// adds a new chunk (chunk of data id) leaf to the selected entries in the merged wave root file 

{
  int estat;
  Long_t id,size,flags,mt;
  char cmd[1024];

  TString mdir = merge_dir;

  TString cwb_merge_label;
  if(gSystem->Getenv("CWB_MERGE_LABEL")==NULL) {
    cout << "cwb_setchunk Error : environment CWB_MERGE_LABEL is not defined!!!" << endl;exit(1);
  } else {
    cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  }

  // check if label has the correct format (M#)
  if(cwb_merge_label[0]!='M') {
    cout << endl;
    cout << "cwb_setchunk Error : label " << cwb_merge_label.Data() << " has bad format (M#)" << endl;
    cout << "Only merge label can be used for cwb_setchunk" << endl << endl;exit(1);
  } else {
    TString lcheck=cwb_merge_label;
    lcheck.Remove(0,1);
    if(!lcheck.IsDigit()) {
      cout << endl;
      cout << "cwb_setchunk Error : label " << cwb_merge_label.Data() << " has bad format (M#)" << endl;
      cout << "Only merge label can be used for cwb_setchunk" << endl << endl;exit(1);
    }
  }

  // get chunk num
  int cwb_setchunk_id=-1;
  if(gSystem->Getenv("CWB_SETCHUNK_ID")==NULL) {
    cout << endl << "cwb_setchunk Error : environment CWB_SETCHUNK_ID is not defined!!!" << endl << endl;
    gSystem->Exit(1);
  } else {
    if(TString(gSystem->Getenv("CWB_SETCHUNK_ID")).IsDigit()) {
      cwb_setchunk_id=TString(gSystem->Getenv("CWB_SETCHUNK_ID")).Atoi();
    } else {
      cout << endl << "cwb_setchunk Error : environment CWB_SETCHUNK_ID is not a positive integer number!!!" << endl << endl;
      gSystem->Exit(1);
    }
  }
  if(cwb_setchunk_id<1) {
    cout << endl << "cwb_setchunk Error : environment CWB_SETCHUNK_ID is not an integer number >0 !!!" << endl << endl;
    gSystem->Exit(1);
  }

  // create input wave root wave file name
  char iwfname[1024];  
  sprintf(iwfname,"wave_%s.%s.root",data_label,cwb_merge_label.Data());

  // create output wave root cuts file name
  TString owfname = mdir+"/"+iwfname;
  char schunk[32];sprintf(schunk,"chunk%d",cwb_setchunk_id);
  owfname.ReplaceAll(".root",TString(".K_")+schunk+".root");

  // apply setChunk to the wave file
  CWB::Toolbox::setChunk(iwfname,mdir,mdir,"waveburst",cwb_setchunk_id);

  if(simulation==0) {

    // create input wave root live file name
    char iwfname[1024];  
    sprintf(iwfname,"live_%s.%s.root",data_label,cwb_merge_label.Data());

    // create output wave root cuts file name
    TString owfname = mdir+"/"+iwfname;
    char schunk[32];sprintf(schunk,"chunk%d",cwb_setchunk_id);
    owfname.ReplaceAll(".root",TString(".K_")+schunk+".root");

    // apply setChunk to the wave file
    CWB::Toolbox::setChunk(iwfname,mdir,mdir,"liveTime",cwb_setchunk_id);

  } else {

    // create input wave root mdc file name
    char iwfname[1024];  
    sprintf(iwfname,"mdc_%s.%s.root",data_label,cwb_merge_label.Data());

    // create output wave root cuts file name
    TString owfname = mdir+"/"+iwfname;
    char schunk[32];sprintf(schunk,"chunk%d",cwb_setchunk_id);
    owfname.ReplaceAll(".root",TString(".K_")+schunk+".root");

    // apply setChunk to the wave file
    CWB::Toolbox::setChunk(iwfname,mdir,mdir,"mdc",cwb_setchunk_id);

  }

  // create a merge*.lst file name & run selection cuts
  vector<int> jobList;
  char ilstfname[1024];  
  sprintf(ilstfname,"merge_%s.%s.lst",data_label,cwb_merge_label.Data());
  TString olstfname = owfname;
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

  exit(0);
}
