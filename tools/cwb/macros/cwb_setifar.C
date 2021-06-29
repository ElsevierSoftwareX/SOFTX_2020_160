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


// adds a new ifar (inverse alse alarm rate) leaf to the selected entries in the merged wave root file 
// for each selected entry a new leaf is created "ifar" which value is obtained from the CWB_SETFAR_FILE
// CWB_SETFAR_FILE is a text file which contains a corrispondence between rho and far (from background)

{
  int estat;
  Long_t id,size,flags,mt;
  char cmd[1024];

          cwb_merge_label      = TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  TString cwb_setifar_tsel     = TString(gSystem->Getenv("CWB_SETIFAR_TSEL"));
  TString cwb_setifar_file     = TString(gSystem->Getenv("CWB_SETIFAR_FILE"));
  TString cwb_setifar_label    = TString(gSystem->Getenv("CWB_SETIFAR_LABEL"));
  TString cwb_setifar_mode     = TString(gSystem->Getenv("CWB_SETIFAR_MODE"));
  TString cwb_setifar_options  = TString(gSystem->Getenv("CWB_SETIFAR_OPTIONS"));

  TString mdir = merge_dir;

  if(cwb_setifar_options.Contains("--")) {
    TString option="";
    // get the cwb_setifar_tsel
    cwb_setifar_tsel = CWB::Toolbox::getParameter(cwb_setifar_options,"--tsel");
    if(cwb_setifar_tsel=="") {
      // look if the explicit tsel definition is defined
      cwb_setifar_tsel = CWB::Toolbox::getParameter(cwb_setifar_options,"--xtsel");
      if(cwb_setifar_tsel=="") {
        cout << "cwb_setifar.C : Error - (--tsel) or (--xtsel) is not defined !!!" << endl << endl;
        gSystem->Exit(1);
      }
    } else {
      // search global definition (defined in *.hh pp cuts file)
      // check if cwb_setifar_tsel is defined

      TString TSEL = cwb_setifar_tsel;

      cwb_setifar_tsel.ReplaceAll("&&"," ");
      cwb_setifar_tsel.ReplaceAll("||"," ");
      cwb_setifar_tsel.ReplaceAll("("," ");
      cwb_setifar_tsel.ReplaceAll(")"," ");
      cwb_setifar_tsel.ReplaceAll("!"," ");

      TObjArray* token = TString(cwb_setifar_tsel).Tokenize(TString(' '));
      for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

         // check if tsel is defined
         TGlobal *global=(TGlobal*)gROOT->GetListOfGlobals()->FindObject(stok.Data());
         if(global==NULL) {
            cout << "cwb_setifar.C : Error - tsel \'" << stok << "\' -> wrong syntax or it is not defined " << endl;
            cout << "                must be included in the user_pparameters.C file" << endl << endl;
            gSystem->Exit(1);
         }
         TCut *tsel = (TCut*)global->GetAddress();
         if(TString(tsel->GetName())!=stok) {
            cout << "cwb_setifar.C : Error - tsel \'" << stok << "\' not correspond to global Tcut name" << endl;
            cout << "                must be included in the user_pparameters.C file" << endl << endl;
            gSystem->Exit(1);
         }

         TSEL.ReplaceAll(stok,TString::Format("(%s)",tsel->GetTitle()));
       }

       cwb_setifar_tsel = TSEL;    // extract TCut 
       cout << "cwb_setifar_tsel : " << cwb_setifar_tsel << endl;
    } 

    // get the cwb_setifar_file
    cwb_setifar_file = CWB::Toolbox::getParameter(cwb_setifar_options,"--file");
    if(cwb_setifar_file=="") {
      // look if the explicit file path definition is defined
      cwb_setifar_file = CWB::Toolbox::getParameter(cwb_setifar_options,"--xfile");
      if(cwb_setifar_file=="") {
        cout << "cwb_setifar.C : Error - (--file) or (--xfile) is not defined !!!" << endl << endl;
        gSystem->Exit(1);
      }
    } else {
      // search global definition far vs rho file
      // check if cwb_setifar_file is defined

      // check if cwb_setifar_file is an index of an array, ex: file_list[3] 
      int beg_bracket = cwb_setifar_file.Index("[");
      int end_bracket = cwb_setifar_file.Index("]");
      int len_bracket = end_bracket-beg_bracket+1;
      bool is_bracket=true;
      if(end_bracket!=cwb_setifar_file.Sizeof()-2) is_bracket=false;
      if(len_bracket<=0) is_bracket=false;
      TString sindx = cwb_setifar_file(beg_bracket+1,len_bracket-2);
      if(!sindx.IsDigit()) is_bracket=false;
      // extract array index
      int indx = is_bracket ? sindx.Atoi() : -1;
      if(indx>=0) cwb_setifar_file=cwb_setifar_file(0,beg_bracket);

      TGlobal *global=(TGlobal*)gROOT->GetListOfGlobals()->FindObject(cwb_setifar_file.Data());
      if(global==NULL) {
         cout << "cwb_setifar.C : Error - setifar_file \'" << cwb_setifar_file << "\' -> wrong syntax or is not defined " << endl;
         cout << "                must be included in the user_pparameters.C file" << endl << endl;
         gSystem->Exit(1);
      }
      int adim = global->GetArrayDim();
      if(adim==0) { 	// string
        cwb_setifar_file = *(TString*)global->GetAddress();
      }
      if(adim==1) { 	// array of strings dim=1
        int asize = global->GetMaxIndex(0); 
        if(indx<0 || indx>=asize) {
          cout << "cwb_setifar.C : Error - setifar_file \'" << cwb_setifar_file << "\' input array index not declared or not allowed, max=" << asize-1 << endl;
          cout << "                check user_pparameters.C file" << endl << endl;
          gSystem->Exit(1);
        }
        TString* afile = (TString*)global->GetAddress();
        cwb_setifar_file = afile[indx]; 
      }
      if(adim>1) {	// array of strings dim>1 (not allowed)
         cout << "cwb_setifar.C : Error - setifar_file \'" << cwb_setifar_file << "\' array dimension > 1 not allowed " << endl;
         cout << "                check user_pparameters.C file" << endl << endl;
         gSystem->Exit(1);
      } 
      cout << "cwb_setifar_file : " << cwb_setifar_file << endl;
    }

    // get the cwb_setifar_mode
    cwb_setifar_mode = CWB::Toolbox::getParameter(cwb_setifar_options,"--mode");
    if(cwb_setifar_mode=="") {
      cout << "cwb_setifar.C : Error - --mode is not defined !!!" << endl;
      gSystem->Exit(1);
    }

    // get the cwb_setifar_label
    cwb_setifar_label = CWB::Toolbox::getParameter(cwb_setifar_options,"--label");
    if(cwb_setifar_label=="") {
      cout << "cwb_setifar.C : Error - --label is not defined !!!" << endl;
      gSystem->Exit(1);
    }
  }

  // set inclusive mode 
  bool inclusive=true;
  cwb_setifar_mode.ToUpper();
  if(cwb_setifar_mode.BeginsWith("I")) inclusive=true;  
  if(cwb_setifar_mode.BeginsWith("E")) inclusive=false;  

  // check if cwb_ifar_label is defined
  if(cwb_setifar_label=="") {
    cout << "cwb_setifar.C : Error - ifar label not defined" << endl; 
    gSystem->Exit(1);
  }
  // check if cwb_setifar_label contains '.'
  if(cwb_setifar_label.Contains(".")) {
    cout << "cwb_setifar.C : Error - cwb_setifar_label " << cwb_setifar_label 
         << " can not contains '.'" << endl << endl;
    gSystem->Exit(1);
  }
  if(cwb_setifar_tsel=="") cwb_wcuts_tree="run>0";

  // create input wave root wave file name
  char iwfname[1024];  
  sprintf(iwfname,"wave_%s.%s.root",data_label,cwb_merge_label.Data());

  // create output wave root cuts file name
  TString owfname = mdir+"/"+iwfname;
  owfname.ReplaceAll(".root",TString(".S_")+cwb_setifar_label+".root");

  // apply setifar to the selected entries in the wave file
  int nsel = CWB::Toolbox::setIFAR(iwfname,mdir,mdir,"waveburst",cwb_setifar_tsel,
                                   cwb_setifar_file,pp_irho,cwb_setifar_label,inclusive);
  if(nsel<=0) {
    cout << "cwb_setifar.C : Warninig - Number of selected waveburst entries = " << nsel << endl << endl;
    //gSystem->Exit(1);
  }
  cout << "cwb_setifar.C : Number of selected waveburst entries = " << nsel << endl;

  // if cwb_setifar_label=same the original wave root file is overwrite
  if(cwb_setifar_label=="same") {
    sprintf(cmd,"/bin/mv %s %s/%s", owfname.Data(), mdir.Data(), iwfname);
    cout << cmd << endl;
    gSystem->Exec(cmd);
    exit(0);
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

  // for simulation!=0 create a symbolic link to mdc*.root file name
  if(simulation) {
    char ilfname[1024];
    sprintf(ilfname,"mdc_%s.%s.root",data_label,cwb_merge_label.Data());
    TString olfname = owfname;
    olfname.ReplaceAll("wave_","mdc_");
    olfname.Remove(0,olfname.Last('/')+1);	// strip path
    cout << olfname << endl;
    estat = gSystem->GetPathInfo(mdir+"/"+ilfname,&id,&size,&flags,&mt);
    if(estat==0) {
      sprintf(cmd,"cd %s;ln -sf %s %s",mdir.Data(),ilfname,olfname.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
    }
  }

  // for simulation=0 create a symbolic link to live*.root file name
  if(!simulation) {
    char ilfname[1024];
    sprintf(ilfname,"live_%s.%s.root",data_label,cwb_merge_label.Data());
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

  exit(0);
}
