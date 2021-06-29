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


// cwb_clonedir.C macro, used to process cwb_clonedir options 

{
  #include <vector>

  CWB::Toolbox TB;

  // get cwb_clonedir options
  TString cwb_clonedir_options=TString(gSystem->Getenv("CWB_CLONEDIR_OPTIONS"));

  // check if cwb_clonedir options have a correct format
  if(gSystem->Getenv("CWB_CLONEDIR_CHECK")!=NULL) {
    bool check;
    CWB::Toolbox::getParameter(cwb_clonedir_options);
    // get output option
    TString cwb_clonedir_options_output = 
            CWB::Toolbox::getParameter(cwb_clonedir_options,"--output");
    cwb_clonedir_options_output.ToUpper();
    if(cwb_clonedir_options_output!="") { // check if the output option is correct
      check=false;
      if(cwb_clonedir_options_output=="LINKS") check=true;
      if(cwb_clonedir_options_output=="MERGE") check=true;
      if(!check) {
        cout << endl << "cwb_clonedir Error : --output available value are : links/merge!!!" << endl << endl;
        gSystem->Exit(1);
      }
    }

    // get config option
    TString cwb_clonedir_options_config = 
            CWB::Toolbox::getParameter(cwb_clonedir_options,"--config");
    cwb_clonedir_options_config.ToUpper();
    if(cwb_clonedir_options_config!="") { // check if the config option is correct
      check=false;
      if(cwb_clonedir_options_config=="CHECK") check=true;
      if(!check) {
        cout << endl << "cwb_clonedir Error : --check available values are : true/false!!!" << endl << endl;
        gSystem->Exit(1);
      }
    }

    // get simulation option
    TString cwb_clonedir_options_simulation = 
            CWB::Toolbox::getParameter(cwb_clonedir_options,"--simulation");
    cwb_clonedir_options_simulation.ToUpper();
    if(cwb_clonedir_options_simulation!="") { // check if the simulation option is correct
      check=false;
      if(cwb_clonedir_options_simulation=="TRUE")  check=true;
      if(cwb_clonedir_options_simulation=="FALSE") check=true;
      if(!check) {
        cout << endl << "cwb_clonedir Error : --simulation available values are : true/false!!!" << endl << endl;
        gSystem->Exit(1);
      }
    }
  
    // get chunk option
    TString cwb_clonedir_options_chunk = 
      CWB::Toolbox::getParameter(cwb_clonedir_options,"--chunk");
    if(cwb_clonedir_options_chunk!="") {
      // get label option
      TString cwb_clonedir_options_label = 
        CWB::Toolbox::getParameter(cwb_clonedir_options,"--label");
      if(cwb_clonedir_options_label!="") {
        cout << endl << "cwb_clonedir Error : --chunk option can not be used when --label option is used!!!" << endl << endl;
        gSystem->Exit(1);
      }
      if(!cwb_clonedir_options_chunk.IsDigit()) {
        cout << endl << "cwb_clonedir Error : ckunk option is not a positive integer number!!!" << endl << endl;
        gSystem->Exit(1);
      }
    }

    gSystem->Exit(0);
  }

  if(cwb_clonedir_options=="") {
    cout << endl;
    cout << "Destination directory already exist ..." << endl;
    cout << endl;
    gSystem->Exit(0);
  }

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  // get directory where cwb_clonedir has been executed
  TString cwb_clonedir_pwd="";
  if(gSystem->Getenv("CWB_CLONEDIR_PWD")!=NULL) {
    cwb_clonedir_pwd=TString(gSystem->Getenv("CWB_CLONEDIR_PWD"));
  } else {
    cout << "cwb_clonedir.C - Error : CWB_CLONEDIR_PWD not declared" << endl;
    exit(1);
  }

  // get cwb_clonedir_src
  TString cwb_clonedir_src="";
  if(gSystem->Getenv("CWB_CLONEDIR_SRC")!=NULL) {
    cwb_clonedir_src=TString(gSystem->Getenv("CWB_CLONEDIR_SRC"));
    if(!cwb_clonedir_src.BeginsWith("/")) {
      cwb_clonedir_src=cwb_clonedir_pwd+"/"+cwb_clonedir_src;
    }
  } else {
    cout << "cwb_clonedir.C - Error : CWB_CLONEDIR_SRC not declared" << endl;
    exit(1);
  }
  TB.checkFile(cwb_clonedir_src);

  // get cwb_clonedir_dest
  TString cwb_clonedir_dest="";
  if(gSystem->Getenv("CWB_CLONEDIR_DEST")!=NULL) {
    cwb_clonedir_dest=TString(gSystem->Getenv("CWB_CLONEDIR_DEST"));
  } else {
    cout << "cwb_clonedir.C - Error : CWB_CLONEDIR_DEST not declared" << endl;
    exit(1);
  }
  
  // get output option
  TString cwb_clonedir_options_output = 
    CWB::Toolbox::getParameter(cwb_clonedir_options,"--output");
  cwb_clonedir_options_output.ToUpper();
  
  // get config option
  TString cwb_clonedir_options_config = 
    CWB::Toolbox::getParameter(cwb_clonedir_options,"--config");
  cwb_clonedir_options_config.ToUpper();
  
  // get simulation option
  TString cwb_clonedir_options_simulation = 
    CWB::Toolbox::getParameter(cwb_clonedir_options,"--simulation");
  cwb_clonedir_options_simulation.ToUpper();
  
  // get chunk option
  TString cwb_clonedir_options_chunk = 
    CWB::Toolbox::getParameter(cwb_clonedir_options,"--chunk");
  
  // get label option
  TString cwb_clonedir_options_label = 
    CWB::Toolbox::getParameter(cwb_clonedir_options,"--label");

  // get stage option
  TString cwb_clonedir_options_jstage = 
    CWB::Toolbox::getParameter(cwb_clonedir_options,"--jstage");
  cwb_clonedir_options_jstage.ToUpper();
  // convert stage name to value
  TString cwb_stage_label="supercluster_";
  if(cwb_clonedir_options_jstage=="FULL")         cwb_stage_label="wave_";
  if(cwb_clonedir_options_jstage=="INIT")         cwb_stage_label="init_";
  if(cwb_clonedir_options_jstage=="STRAIN")       cwb_stage_label="strain_";
  if(cwb_clonedir_options_jstage=="CSTRAIN")      cwb_stage_label="cstrain_";
  if(cwb_clonedir_options_jstage=="COHERENCE")    cwb_stage_label="coherence_";
  if(cwb_clonedir_options_jstage=="SUPERCLUSTER") cwb_stage_label="supercluster_";
  if(cwb_clonedir_options_jstage=="LIKELIHOOD")   cwb_stage_label="wave_";

  char src_output_dir[1024];
  sprintf(src_output_dir,"%s/%s",cwb_clonedir_src.Data(),output_dir);

  // get rid of final "/" to makes gSystem->BaseName working properly
  if(cwb_clonedir_src.EndsWith("/")) cwb_clonedir_src.Resize(cwb_clonedir_src.Sizeof()-2);
  // extract the source data label
  const char* src_data_label = gSystem->BaseName(cwb_clonedir_src.Data());
 
  // check if source config files are compatible with the dest config files 
  if(cwb_clonedir_options_config=="CHECK") {
    char cmd[1024];
    char fparms[1024];

    sprintf(fparms,"%s/%s",cwb_clonedir_src.Data(),gSystem->ExpandPathName("$CWB_UPARAMETERS_FILE"));

    // check if files exist
    TB.checkFile(fparms);

    cout << "fparms   : " << fparms << endl;

    // diff config files
    sprintf(cmd,"diff %s %s | less",fparms,gSystem->ExpandPathName("$CWB_UPARAMETERS_FILE"));
    gSystem->Exec(cmd);

    char answer[256];
    strcpy(answer,"");
    do {
      cout << "Are src & dest config files compatibles ? (y/n) " << endl;
      cin >> answer;
      cout << endl << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"n")==0) gSystem->Exit(1);
  }
 
  // Create symbolic links under the dest output dir of the source output dir files
  // The link data_label is the dest data_label
  if(cwb_clonedir_options_output=="LINKS") {
    vector<TString> fileList = TB.getFileListFromDir(src_output_dir, ".root", 
                                                     cwb_stage_label, src_data_label);
    cout << "cwb_clonedir : create symbolic links in output dir ..." << endl;
    int newlinks=0;
    int estat;
    Long_t id,size,flags,mt;
    char cmd[1024];
    char link_to[1024];
    for(int i=0;i<fileList.size();i++) {
      if(i%100==0) cout << i << "/" << fileList.size() << endl;
      //cout << i << " " << fileList[i].Data() << endl;
      const char* file_name_base = gSystem->BaseName(fileList[i].Data());
      TString sfile_name_base = file_name_base;
      sfile_name_base.ReplaceAll(src_data_label,data_label);

      TString filePath = TString(output_dir)+"/"+sfile_name_base;
      estat = gSystem->GetPathInfo(filePath,&id,&size,&flags,&mt);
      if (estat!=0) {
        newlinks++;
        sprintf(link_to,"%s/%s",output_dir,sfile_name_base.Data());
        gSystem->Symlink(fileList[i].Data(),link_to);
      }
    }
    cout << "Number of new links " << newlinks << "/" << fileList.size() << endl;
  }

  char src_merge_dir[1024];
  sprintf(src_merge_dir,"%s/%s",cwb_clonedir_src.Data(),merge_dir);
 
  // merge wave&live under the dest output dir of the source merged files
  if(cwb_clonedir_options_output=="MERGE") {
    int estat;
    Long_t id,size,flags,mt;
    vector<TString> fileList = TB.getFileListFromDir(src_merge_dir, ".root", "", src_data_label);
    cout << "cwb_clonedir : merge src wave&live in output dir ..." << endl;
    TRegexp reg1(".*M[0-9]+.root");
    // compute max merge version
    int max_version=0;
    for(int i=0;i<fileList.size();i++) {
      // select merged file with *.M#.root format
      if(fileList[i].Contains(reg1)) {
        //cout << i << " " << fileList[i] << endl;
        TObjArray* token = TString(fileList[i]).Tokenize(TString("."));
        TString smergeID = ((TObjString*)token->At(token->GetEntries()-2))->GetString();
        smergeID.ReplaceAll("M","");
        if(smergeID.IsDigit()) {
          int mergeID = smergeID.Atoi();
          if(max_version<mergeID) max_version=mergeID;
        }
        delete token;
      }
    }
    // cout << "max_version : " << max_version << endl;
    TString chunk_label = cwb_clonedir_options_chunk!="" ? TString(".K_chunk")+cwb_clonedir_options_chunk : "";
    char fwave[1024];
    char flive[1024];
    char fmerge[1024];
    char flist[1024];
    char fmdc[1024];
    if(cwb_clonedir_options_label=="") {
      sprintf(fwave,"%s/wave_%s.M%d%s.root",src_merge_dir,src_data_label,max_version,chunk_label.Data());
      sprintf(flive,"%s/live_%s.M%d%s.root",src_merge_dir,src_data_label,max_version,chunk_label.Data());
      sprintf(fmerge,"%s/wave_%s.M%d%s.root",output_dir,src_data_label,max_version,chunk_label.Data());
      sprintf(flist,"%s/wave_%s.M%d%s.lst",output_dir,src_data_label,max_version,chunk_label.Data());
      sprintf(fmdc,"%s/mdc_%s.M%d%s.root",src_merge_dir,src_data_label,max_version,chunk_label.Data());
    } else {
      sprintf(fwave,"%s/wave_%s.%s.root",src_merge_dir,src_data_label,cwb_clonedir_options_label.Data());
      sprintf(flive,"%s/live_%s.%s.root",src_merge_dir,src_data_label,cwb_clonedir_options_label.Data());
      sprintf(fmerge,"%s/wave_%s.%s.root",output_dir,src_data_label,cwb_clonedir_options_label.Data());
      sprintf(flist,"%s/wave_%s.%s.lst",output_dir,src_data_label,cwb_clonedir_options_label.Data());
      sprintf(fmdc,"%s/mdc_%s.%s.root",src_merge_dir,src_data_label,cwb_clonedir_options_label.Data());
    }
    cout << "fwave  : " << fwave << endl;
    cout << "flive  : " << flive << endl;
    cout << "fmerge : " << fmerge << endl;
    cout << "flist  : " << flist << endl;
    cout << "fmdc   : " << fmdc << endl;
    // check if files exist
    estat = gSystem->GetPathInfo(fwave,&id,&size,&flags,&mt);
    if (estat!=0) {cout << "cwb_clonedir.C : Error - wave not exist : " << fwave << endl;exit(1);}
    if(cwb_clonedir_options_simulation=="TRUE") {
      estat = gSystem->GetPathInfo(fmdc,&id,&size,&flags,&mt);
      if (estat!=0) {cout << "cwb_clonedir.C : Error - mdc not exist : " << fmdc << endl;exit(1);}
    } else {
      estat = gSystem->GetPathInfo(flive,&id,&size,&flags,&mt);
      if (estat!=0) {cout << "cwb_clonedir.C : Error - live not exist : " << flive << endl;exit(1);}
    }
    // check if file with a different version .M# already exist in the output_dir
    vector<TString> checkList = TB.getFileListFromDir(output_dir, "", "", src_data_label);
    if(checkList.size()!=0) {
      cout << endl;
      cout << "cwb_clonedir.C : Error - files with source label : " << src_data_label << endl;
      cout << "already exist in the output directory : " << cwb_clonedir_dest << 
              "/" << output_dir << endl;
      cout << "Before to apply cwb_clonedir.C the following files must be removed !!!" << endl << endl;
      for(int i=0;i<checkList.size();i++) {
        cout << i << " " << checkList[i].Data() << endl;
      }
      gSystem->Exit(1);
    }
    if(chunk_label!="") {
      // check if file with same chunk label .K_chunk# already exist in the output_dir
      checkList = TB.getFileListFromDir(output_dir, "", "", chunk_label);
      if(checkList.size()!=0) {
        cout << endl;
        cout << "cwb_clonedir.C : Error - files with chunk label : " << chunk_label << endl;
        cout << "already exist in the output directory : " << cwb_clonedir_dest << 
                "/" << output_dir << endl;
        cout << "Before to apply cwb_clonedir.C the following files must be removed !!!" << endl << endl;
        for(int i=0;i<checkList.size();i++) {
          cout << i << " " << checkList[i].Data() << endl;
        }
        gSystem->Exit(1);
      }
    }
    // merge wave&live
    TFileMerger M;
    M.AddFile(fwave,true);
    if(cwb_clonedir_options_simulation=="TRUE") {
      M.AddFile(fmdc,true);
    } else {
      M.AddFile(flive,true);
    }
    M.OutputFile(fmerge);
    if(!M.Merge()) {
      char cmd[1024];
      sprintf(cmd,"rm %s",fmerge);
      //cout << cmd << endl;
      gSystem->Exec(cmd);
      cout << "cwb_clonedir.C : Error - Merge failed !!!" << endl;
      gSystem->Exit(1);
    }
    // write merge file list
    vector<TString> mergeList(2);
    mergeList[0]=fwave;
    if(cwb_clonedir_options_simulation=="TRUE") {
      mergeList[1]=fmdc;
    } else { 
      mergeList[1]=flive;
    }
    TB.dumpFileList(mergeList, flist);
  }

  gSystem->Exit(0);
}
