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


// obsolete

void 
cwb_jnet(TString jName="", TString uName="") {

  bool jshell = jName!="" ? true : false;

  // if jName!="" jfile is opened from disk otherwise it is opened within TBrowser  
  TFile* jfile = jName=="" ? gROOT->GetFile() : jfile = new TFile(jName);
  if(jfile==NULL) {
    cout << "cwb_jnet - Error opening root file : " << jName << endl;
    if(jshell) exit(1); else return;
  }
  jfile->ls();

  // read config object
  CWB::config* cfg = (CWB::config*)jfile->Get("config");
  if(cfg==NULL) {
    cout << "cwb_jnet - config is not present in : " << jName << endl;
    if(jshell) exit(1); else return;
  }
  // export config to cint
  cfg->Export();
  // set auxiliary configuration
  if(uName!="") cfg.Import(uName);  // auxiliary config file
  //cfg->Print();

  // get analysis setup
  TString ANALYSYS = cfg->analysis;

  // get jName ...
  if(jName=="") {
    jName = jfile->GetPath();
    jName.ReplaceAll(":/","");
    cout << jName.Data() << endl;
  }

  jfile->Close();

  // create working directory jdir
  TObjArray* token = TString(jName).Tokenize(TString('/'));
  TObjString* stoken =(TObjString*)token->At(token->GetEntries()-1);
  int jobID=TString(stoken->GetString().ReplaceAll("job","").ReplaceAll(".root","")).Atoi();
  TString jdir = stoken->GetString().ReplaceAll(".root","");

  // -------------------------------------------------------------------------
  // Check if output file already exists and asks if you want to overwrite it
  // -------------------------------------------------------------------------
  bool overwrite=true;
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(jdir.Data(),&id,&size,&flags,&mt);
  if (estat==0) {
    char answer[256];
    strcpy(answer,"");
    do {
      cout << "Dir " << jdir.Data() << " already exist" << endl;
      cout << "Do you want to overwrite it ? (y/n) ";
      cin >> answer;
      cout << endl << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"n")==0) overwrite=false;
  } 
  if(overwrite) {
    // clean working directories
    gSystem->Exec(TString("rm -rf ")+jdir+"/"+config_dir+"/*");
    gSystem->Exec(TString("rm -rf ")+jdir+"/"+data_dir+"/*");
    gSystem->Exec(TString("rm -rf ")+jdir+"/"+tmp_dir+"/*");
  } else return;

  gSystem->Exec(TString("mkdir -p ")+jdir+"/"+config_dir);
  gSystem->Exec(TString("mkdir -p ")+jdir+"/"+data_dir);
  gSystem->Exec(TString("mkdir -p ")+jdir+"/"+tmp_dir);
  TString juser_parameters = "config/juser_parameters.C";
  if(uName!="") gSystem->Exec(TString("cp ")+uName+" "+jdir+"/"+juser_parameters);
  gSystem->cd(jdir);
//  else          juser_parameters="";
//  else          gSystem->Exec(TString("touch ")+jdir+"/"+juser_parameters);
  if(uName=="") {
    gSystem->Exec(TString("echo \"{\n\" >> ")+juser_parameters);  
    gSystem->Exec(TString("echo \"jobfOptions=CWB_JOBF_SAVE_ALL;\n\" >> ")+juser_parameters);  
    gSystem->Exec(TString("echo \"}\n\" >> ")+juser_parameters);  
  }
//  gSystem->Exec(TString("echo \"strcpy(config_dir,\"config\");\n\" >> ")+juser_parameters);  
  jName="../"+jName;

  // delete temp cfg
  delete cfg;

  int runID=1;
  if(ANALYSYS=="1G") {
    cwb1G CWB(jName,juser_parameters);
    CWB.run(runID);
  } else
  if(ANALYSYS=="2G") {
    cwb2G CWB(jName,juser_parameters);
    CWB.run(runID);
  } else {
    cout << "cwb_jnet - Error : analysis must be 1G or 2G" << endl;
    if(jshell) exit(1); else return;
  }

  jshell ? exit(0) : return;
}
