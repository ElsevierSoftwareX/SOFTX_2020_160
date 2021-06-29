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


// creates cwb working directory and subdirectories

//#define SYMBOLIC_CONDOR_LOG_DIR // used in ATLAS cluster to fix issue with HSM

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));

  bool cwb_mkdir_batch=false;
  if(TString(gSystem->Getenv("CWB_MKDIR_OPTION")).CompareTo("batch")==0) cwb_mkdir_batch=true;
  TString cwb_mkdir_wrkdir="";;
  if(gSystem->Getenv("CWB_MKDIR_WRKDIR")!=NULL) {
    // get working dir
    // if cwb_mkdir_wrkdir=""  : the working dir is the current
    // if cwb_mkdir_wrkdir!="" : the working dir is created
    cwb_mkdir_wrkdir=TString(gSystem->Getenv("CWB_MKDIR_WRKDIR"));
    if(cwb_mkdir_wrkdir!="") {
      // if present then strip last '/'
      if(cwb_mkdir_wrkdir[cwb_mkdir_wrkdir.Sizeof()-2]=='/') {
        cwb_mkdir_wrkdir.Resize(cwb_mkdir_wrkdir.Sizeof()-2); 
      }
    }
  }

  // check if wrkdir name is valid
  if(cwb_mkdir_wrkdir.Contains(".M")) {
    cout << "cwb_mkdir - Error : working directory name is not a valid name ..." << endl;
    cout << "\"" << cwb_mkdir_wrkdir << "\"" << endl;
    cout << "'.M' can not be used : it is used to label the merge data" << endl;
    cout << "cwb_mkdir aborted" << endl << endl; 
    gSystem->Exit(1);
  }
  if(cwb_mkdir_wrkdir.Contains(".V")) {
    cout << "cwb_mkdir - Error : working directory name is not a valid name ..." << endl;
    cout << "\"" << cwb_mkdir_wrkdir << "\"" << endl;
    cout << "'.V' can not be used : it is used to label the vetoed data" << endl;
    cout << "cwb_mkdir aborted" << endl << endl; 
    gSystem->Exit(1);
  }
  if(cwb_mkdir_wrkdir.Contains(".C")) {
    cout << "cwb_mkdir - Error : working directory name is not a valid name ..." << endl;
    cout << "\"" << cwb_mkdir_wrkdir << "\"" << endl;
    cout << "'.C' can not be used : it is used to label the cutted data" << endl;
    cout << "cwb_mkdir aborted" << endl << endl; 
    gSystem->Exit(1);
  }

  // get cluster site
  TString site_cluster = "";;
  if(gSystem->Getenv("SITE_CLUSTER")!=NULL) {
    site_cluster = TString(gSystem->Getenv("SITE_CLUSTER"));
  }

  #define NDIR 20

  TString dir_name[NDIR];
  int nDIR=0;
  if(cwb_mkdir_wrkdir!="") {
    dir_name[nDIR++] = cwb_mkdir_wrkdir;
    cout << endl;
    bool overwrite=TB.checkFile(cwb_mkdir_wrkdir,true,"working dir already exist !!!");   
    if(!overwrite) {cout << "cwb_mkdir aborted" << endl << endl; gSystem->Exit(1);}
  }
  dir_name[nDIR++] = tmp_dir;
  dir_name[nDIR++] = data_dir;
  dir_name[nDIR++] = config_dir;
  dir_name[nDIR++] = input_dir;
  dir_name[nDIR++] = output_dir;
  dir_name[nDIR++] = macro_dir;
  dir_name[nDIR++] = report_dir;
  dir_name[nDIR++] = condor_log;
#ifdef SYMBOLIC_CONDOR_LOG_DIR
  if(site_cluster=="ATLAS") {
    // used in ATLAS cluster to fix issue with HSM
    // condor log dir is in the local header nodes
    if(cwb_mkdir_wrkdir!="") {
      dir_name[nDIR++] = TString(condor_log_dir)+TString("/")+gSystem->BaseName(cwb_mkdir_wrkdir);
    } else {
      dir_name[nDIR++] = TString(condor_log_dir)+TString("/")++TString(data_label);
    }
    Long_t id,size=0,flags,mt;
    int estat = gSystem->GetPathInfo(dir_name[nDIR-1].Data(),&id,&size,&flags,&mt);
    if(estat==0) {
      cout << endl;
      cout << "cwb_mkdir - Error : condor log directory ..." << endl;
      cout << "\"" << dir_name[nDIR-1] << "\"" << endl;
      cout << "already exist, select a unique name for the working directory" << endl;
      cout << "cwb_mkdir aborted" << endl << endl; 
      gSystem->Exit(1);
    }
  } else {
    dir_name[nDIR++] = log_dir;
  }
#else
  dir_name[nDIR++] = log_dir;
#endif
  dir_name[nDIR++] = condor_dir;
  dir_name[nDIR++] = merge_dir;
  dir_name[nDIR++] = ced_dir;
  dir_name[nDIR++] = pp_dir;
  dir_name[nDIR++] = dump_dir;
  // www report dir is created only if the env WWW_PUBLIC_DIR  != ""
  if(TString(www_dir)!="") {
    if(cwb_mkdir_wrkdir!="") dir_name[nDIR] = TString(www_dir)+TString("/")+gSystem->BaseName(cwb_mkdir_wrkdir);
    else                     dir_name[nDIR] = TString(www_dir)+TString("/")+TString(data_label);
    dir_name[nDIR].ReplaceAll(report_dir+TString("/"),""); 
    // check if www dir exist, must be unique !!!
    cout << endl;
    bool overwrite=TB.checkFile(dir_name[nDIR],true,"www report directory must be unique !!!");   
    if(!overwrite) {cout << "cwb_mkdir aborted" << endl << endl; gSystem->Exit(1);}
    nDIR++;
  }
/*
  dir_name[nDIR]   = TString(www_dir)+TString("/")+TString(data_label)+TString("/")+TString(pp_dir);
  dir_name[nDIR].ReplaceAll(report_dir+TString("/"),""); 
  nDIR++;
  dir_name[nDIR]   = TString(www_dir)+TString("/")+TString(data_label)+TString("/")+TString(ced_dir);
  dir_name[nDIR].ReplaceAll(report_dir+TString("/"),""); 
  nDIR++;
  dir_name[nDIR]   = TString(www_dir)+TString("/")+TString(data_label)+TString("/")+TString(dump_dir);
  dir_name[nDIR].ReplaceAll(report_dir+TString("/"),""); 
  nDIR++;
*/

  char cmd[1024];  
  char ldir[1024];  
  for(int i=0;i<nDIR;i++) {

    if(dir_name[i]=="") continue;

    sprintf(ldir,"%s",dir_name[i].Data());

    // ----------------------------------------------------------------
    // Check if dir lists exist
    // ----------------------------------------------------------------
    Long_t id,size,flags,mt;
    int estat = gSystem->GetPathInfo(ldir,&id,&size,&flags,&mt);
    if (estat==0) {
      char answer[256];
      strcpy(answer,"");
      do {
        cout << endl;
        cout << "dir \"" << ldir << "\" already exist" << endl;
        cout << "Do you want to remove the files ? (y/n) ";
        if(cwb_mkdir_batch) strcpy(answer,"y"); else cin >> answer;
        cout << endl;
      } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
      if (strcmp(answer,"y")==0) {
        sprintf(cmd,"mkdir -p %s",ldir);
        cout << cmd << endl;
        gSystem->Exec(cmd);
      }
    } else {
      sprintf(cmd,"mkdir -p %s",ldir);
      cout << cmd << endl;
      gSystem->Exec(cmd);
//      sprintf(cmd,"rm %s/*",ldir);
//      cout << cmd << endl;
//      gSystem->Exec(cmd);
    }

    // if working dir is created then cd to wrkdir
    if(i==0 && cwb_mkdir_wrkdir!="") gSystem->ChangeDirectory(cwb_mkdir_wrkdir.Data());
  }

#ifdef SYMBOLIC_CONDOR_LOG_DIR
  if(site_cluster=="ATLAS") {
    // used in ATLAS cluster to fix issue with HSM
    // condor log dir is in the local header nodes
    if(cwb_mkdir_wrkdir!="") {
      sprintf(cmd,"ln -sf %s/%s %s",condor_log_dir.Data(),
              gSystem->BaseName(cwb_mkdir_wrkdir),log_dir);
    } else {
      sprintf(cmd,"ln -sf %s/%s %s",condor_log_dir.Data(),data_label,log_dir);
    }
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }
#endif

  TString cwb_rootlogon_file = TString(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  sprintf(cmd,"ln -sf %s rootlogon.C",cwb_rootlogon_file.Data());
  //cout << cmd << endl;
  //gSystem->Exec(cmd);

  TString cwb_dir_name = TString(gSystem->Getenv("CWB_MACROS"));
  sprintf(cmd,"ln -sf %s/README.cwb",cwb_dir_name.Data());
  //cout << cmd << endl;
  //gSystem->Exec(cmd);

  TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
  TString cwb_nets_name = TString(gSystem->Getenv("CWB_NETS_FILE"));
  if(gSystem->Getenv("_USE_PEGASUS")!=NULL) {
    // symbolic links for pegasus scripts
    sprintf(cmd,"ln -sf %s %s/cwb.sh",cwb_nets_name.Data(),condor_dir);
    cout << cmd << endl;
    gSystem->Exec(cmd);
  } else {
    // symbolic links for condor scripts
    sprintf(cmd,"ln -sf %s/cwb_loudest.sh %s/loudest.sh",cwb_scripts.Data(),condor_dir);
    cout << cmd << endl;
    gSystem->Exec(cmd);

    sprintf(cmd,"ln -sf %s/cwb_ced.sh %s/ced.sh",cwb_scripts.Data(),condor_dir);
    cout << cmd << endl;
    gSystem->Exec(cmd);

    sprintf(cmd,"ln -sf %s %s/cwb.sh",cwb_nets_name.Data(),condor_dir);
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  // www report links are created only if the env WWW_PUBLIC_DIR  != ""
  if(TString(www_dir)!="") {
    // data label
    char www_label[1024];
    TObjArray* token = TString(work_dir).Tokenize(TString("/"));
    sprintf(www_label,"%s",((TObjString*)token->At(token->GetEntries()-1))->GetString().Data());
    if(cwb_mkdir_wrkdir!="") {
      TObjArray* token2 = TString(cwb_mkdir_wrkdir).Tokenize(TString("/"));
      sprintf(www_label,"%s",((TObjString*)token2->At(token2->GetEntries()-1))->GetString().Data());
    }

    TString cwb_dir = cwb_mkdir_wrkdir.BeginsWith("/") ? cwb_mkdir_wrkdir : TString(work_dir)+"/"+cwb_mkdir_wrkdir;

    sprintf(cmd,"ln -sf %s/%s %s/%s",cwb_dir.Data(),pp_dir,www_dir,www_label);
    cout << cmd << endl;
    gSystem->Exec(cmd);

    sprintf(cmd,"ln -sf %s/%s %s/%s",cwb_dir.Data(),ced_dir,www_dir,www_label);
    cout << cmd << endl;
    gSystem->Exec(cmd);

    sprintf(cmd,"ln -sf %s/%s %s/%s",cwb_dir.Data(),dump_dir,www_dir,www_label);
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  exit(0);
}
