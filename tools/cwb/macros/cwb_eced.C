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


// 1G easy ced : an interactive command to produce CED

{
  char cmd[1024];

  segLen=60;
  cedRHO=0.0;   // disable the cedRHO threshold
  simulation=0;
  factors[0]=1;

  // variables where to save the user parameters definitions
  int _nIFO=0;
  TString _ifo[NIFO_MAX];
  TString _channelNamesRaw[NIFO_MAX];
  TString _channelNamesMDC[NIFO_MAX];
  TString _frFiles[2*NIFO_MAX];
  dqfile _DQF[20];int _nDQF=0;

  char channelTypesRaw[NIFO_MAX][1024]; // detector data types - used only for cwb_eced by gw_data_find
  for(int i=0;i<NIFO_MAX;i++) strcpy(channelTypesRaw[i],"");

  unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
  char tag[256];sprintf(tag,"%d",Pid);

  // get cluster site
  TString site_cluster="";
  if(gSystem->Getenv("SITE_CLUSTER")!=NULL) {
    site_cluster=TString(gSystem->Getenv("SITE_CLUSTER"));
  }

  // get env options
  int xgps_event=0;
  TString cwb_eced_opts=TString(gSystem->Getenv("CWB_ECED_OPTS"));
  TString eced_cfg="";
  if(cwb_eced_opts!="") {
    if(cwb_eced_opts.IsFloat()) {    // env CWB_ECED_OPTS contains only the gps_event
      xgps_event=cwb_eced_opts.Atoi();
      cwb_eced_opts=""; 
    } else {                         // more options
      TString ECED_CFG = CWB::Toolbox::getParameter(cwb_eced_opts,"--cfg");
      if(ECED_CFG!="") {
        CWB::Toolbox::checkFile(gSystem->ExpandPathName(ECED_CFG.Data()));
        gROOT->Macro(gSystem->ExpandPathName(ECED_CFG.Data())); // set cfg user define cwb parameters
        eced_cfg=ECED_CFG;
      }
  
      TString ECED_GPS = CWB::Toolbox::getParameter(cwb_eced_opts,"--gps");          // set event gps time
      if(ECED_GPS.IsFloat()) xgps_event=ECED_GPS.Atoi(); 
  
      TString ECED_FACTOR = CWB::Toolbox::getParameter(cwb_eced_opts,"--factor");    // set factor
      if(ECED_FACTOR.IsFloat()) factors[0]=ECED_FACTOR.Atof(); 

      TString ECED_SEARCH = CWB::Toolbox::getParameter(cwb_eced_opts,"--search");    // set search type
      if(ECED_SEARCH.Sizeof()==2) SEARCH()=ECED_SEARCH[0];
  
      TString ECED_SIM = CWB::Toolbox::getParameter(cwb_eced_opts,"--sim");          // set simulation mode
      if(ECED_SIM.IsFloat()) { 
        int sim = ECED_SIM.Atoi(); 
        if((sim<0)&&(sim>2)) { 
          cout << "Error : bad --sim value [0/1/2] !!!" << endl;
          gSystem->Exit(1);
        } else simulation=sim;
      }

      TString ECED_TAG = CWB::Toolbox::getParameter(cwb_eced_opts,"--tag");          // set directory user TAG
      if(ECED_TAG!="") strcpy(tag,ECED_TAG.Data());
    }
  }
  if(xgps_event<=0) {
    cout << "cwb_eced - Error : gps_event must be > 0, check config/user_parameters.C" << endl;
    gSystem->Exit(1);
  } else {
    char sgps_event[16];sprintf(sgps_event,"%d",xgps_event);
    gSystem->Setenv("CWB_GPS_EVENT",sgps_event);
  }

  // get evn ndet parameters
  if(gSystem->Getenv("CWB_ECED_NDET")!=NULL) {
    if(TString(gSystem->Getenv("CWB_ECED_NDET")).IsDigit()) {
      nIFO=TString(gSystem->Getenv("CWB_ECED_NDET")).Atoi();
      if(nIFO==0) {
        cout << "cwb_eced - Error : ifos are not defined" << endl;
        gSystem->Exit(1);
      }
    } else {
      cout << "Error : environment CWB_ECED_NDET is not defined!!!" << endl;
      gSystem->Exit(1);
    }
  }

  nDQF=0;
  for(int n=0;n<nIFO;n++) {
    TString ecedOptions;
    char cwb_eced_det[64]; sprintf(cwb_eced_det,"CWB_ECED_DET%d",n+1); 
    if(gSystem->Getenv(cwb_eced_det)!=NULL) {
      ecedOptions=TString(gSystem->Getenv(cwb_eced_det));

      // if ecedOptions = L1,H1,H2,V1,G1
      if(ecedOptions.Sizeof()==3) {
        strcpy(ifo[n],ecedOptions.Data());
        ecedOptions=""; 
      } else {
        TString ECED_IFO = CWB::Toolbox::getParameter(ecedOptions,"--ifo");
        if(ECED_IFO!="") strcpy(ifo[n],ECED_IFO.Data());
        else {cout << "Error : Inline param --ifo is not defined!!!" << endl;exit(1);}
      }

      TString ECED_TYPE = CWB::Toolbox::getParameter(ecedOptions,"--type");
      if(ECED_TYPE!="") strcpy(channelTypesRaw[n],ECED_TYPE.Data());
      else { // search the definitions in the user configuration file
        for(int i=0;i<_nIFO;i++) if(_ifo[i]==ifo[n]) strcpy(channelTypesRaw[n],_frFiles[i].Data());
        for(int i=0;i<_nIFO;i++) if(_ifo[i]==ifo[n]) strcpy(channelTypesRaw[n+nIFO],_frFiles[i+_nIFO].Data());
      }

      TString ECED_CHRAW = CWB::Toolbox::getParameter(ecedOptions,"--chraw");
      if(ECED_CHRAW!="") strcpy(channelNamesRaw[n],ECED_CHRAW.Data());
      else { // search the definitions in the user configuration file
        for(int i=0;i<_nIFO;i++) if(_ifo[i]==ifo[n]) strcpy(channelNamesRaw[n],_channelNamesRaw[i].Data());
      }

      TString ECED_CHMDC = CWB::Toolbox::getParameter(ecedOptions,"--chmdc");
      if(ECED_CHMDC!="") strcpy(channelNamesMDC[n],ECED_CHMDC.Data());
      else { // search the definitions in the user configuration file
        for(int i=0;i<_nIFO;i++) if(_ifo[i]==ifo[n]) strcpy(channelNamesMDC[n],_channelNamesMDC[i].Data());
      }

      TString ECED_SHIFT = CWB::Toolbox::getParameter(ecedOptions,"--shift");
      if(ECED_SHIFT!="") {
        if(ECED_SHIFT.IsFloat()) {
          dataShift[n] = ECED_SHIFT.Atof();  
          if(fmod(dataShift[n],1)!=0)
            {cout << "Error : shift parameter is not an integer number!!! " << ECED_SHIFT << endl;gSystem->Exit(1);}
          for(int i=0;i<nDQF;i++) DQF[i].shift+=ECED_SHIFT.Atof();
        } else {cout << "Error : shift parameter is not a number!!! " << ECED_SHIFT << endl;gSystem->Exit(1);}
      } else { // search the definitions in the user configuration file
        //for(int i=0;i<_nIFO;i++) if(_ifo[i]==ifo[n]) strcpy(dataShift[n],_channelNamesMDC[i].Data());
      }

      // copy DQ definitions
      for(int i=0;i<_nDQF;i++) if(TString(_DQF[i].ifo)==ifo[n]) DQF[nDQF++]=_DQF[i];

      cout << n << " ifo " << ifo[n] << " type " << channelTypesRaw[n] 
           << " chraw " << channelNamesRaw[n] << " chmdc " << channelNamesMDC[n] << endl; 
    }
  }
  for(int i=0;i<nDQF;i++) cout << "DQF : " << i << " " << DQF[i].ifo << " " <<  DQF[i].file << endl;

  // set reference ifo 
  strcpy(refIFO,ifo[0]);

  // select default job segment parameters
  if(segLen<60) segLen=60;
  segMLS = segLen;
  segTHR = 0;

  // force nfactor=1
  nfactor = 1;     
  // for background force factors[0]=1
  if(simulation==0) factors[0]=1;

  // select a segment of 60 sec around the event gps time
  int gps_start = int(xgps_event)-(int(segLen/2.+1)+segEdge);
  int gps_stop  = int(xgps_event)+(int(segLen/2.+1)+segEdge);

  // if not defined then set default Raw channel names  & channel types
  for(int n=0;n<nIFO;n++) {

    TString ECED_IFO = ifo[n];

    if(ECED_IFO=="L1") {
      if(TString(channelNamesRaw[n])=="") sprintf(channelNamesRaw[n],"L1:LDAS-STRAIN");
      if(TString(channelTypesRaw[n])=="") sprintf(channelTypesRaw[n],"L1_LDAS_C02_L2");
    } else 
    if(ECED_IFO=="H1") {
      if(TString(channelNamesRaw[n])=="") sprintf(channelNamesRaw[n],"H1:LDAS-STRAIN");
      if(TString(channelTypesRaw[n])=="") sprintf(channelTypesRaw[n],"H1_LDAS_C02_L2");
    } else 
    if(ECED_IFO=="H2") {
      if(TString(channelNamesRaw[n])=="") sprintf(channelNamesRaw[n],"H2:LDAS-STRAIN");
      if(TString(channelTypesRaw[n])=="") sprintf(channelTypesRaw[n],"H2_LDAS_C02_L2");
    } else 
    if(ECED_IFO=="V1") {
      if(TString(channelNamesRaw[n])=="") sprintf(channelNamesRaw[n],"V1:h_16384Hz");
      if(TString(channelTypesRaw[n])=="") sprintf(channelTypesRaw[n],"HrecV2");
    } else 
    if(ECED_IFO=="G1") {
      if(TString(channelNamesRaw[n])=="") sprintf(channelNamesRaw[n],"G1:DER_DATA_H");
      if(TString(channelTypesRaw[n])=="") sprintf(channelTypesRaw[n],"G1_RDS_C01_L3");
    }  
  }

  TString cwb_inet_opts=TString(gSystem->Getenv("CWB_INET_OPTIONS"));
  // define eced directory
  char eced_dir[1024]="ECED_";
  for(int n=0;n<nIFO;n++) sprintf(eced_dir,"%s%s",eced_dir,ifo[n]);
  sprintf(eced_dir,"%s_TAG%s",eced_dir,tag);
  if(cwb_inet_opts=="ced") {
    sprintf(eced_dir,"%s_GPS%d",eced_dir,xgps_event);
    // Check if eced_dir exist
    Long_t id,size,flags,mt;
    int estat = gSystem->GetPathInfo(eced_dir,&id,&size,&flags,&mt);
    if(estat==0) {
      cout << "cwb_eced - Error : directory " << eced_dir << " already exist " << endl; 
      gSystem->Exit(1);
    }
  } 
  // redefine eced working directories
  sprintf(tmp_dir,   "%s/tmp",eced_dir);
  sprintf(data_dir,  "%s/tmp",eced_dir);
  sprintf(config_dir,"%s/tmp",eced_dir);
  sprintf(input_dir, "%s/tmp",eced_dir);
  sprintf(output_dir,"%s/tmp",eced_dir);
  sprintf(macro_dir, "%s/tmp",eced_dir);
  sprintf(report_dir,"%s/tmp",eced_dir);
  sprintf(log_dir,   "%s/tmp",eced_dir);
  sprintf(condor_dir,"%s/tmp",eced_dir);
  sprintf(merge_dir, "%s/tmp",eced_dir);
  sprintf(pp_dir,    "%s/tmp",eced_dir);
  sprintf(dump_dir,  "%s/tmp",eced_dir);
  if(cwb_inet_opts=="ced") {
    sprintf(ced_dir,   "%s/ced",eced_dir);
  } else {
    sprintf(ced_dir,   "%s/output",eced_dir);
  }

  // check if frame files are available with gw_data_find command 
  // fill ldf_command array to get frame files list for each detector 
  bool ldf_gcheck=true; 
  bool ldf_check[NIFO_MAX]; 
  char ldf_command[NIFO_MAX][1024];  
  for(int i=0;i<NIFO_MAX;i++) {strcpy(ldf_command[i],"");ldf_check[i]=true;}
  for(int n=0;n<nIFO;n++) {

    if(TString(channelTypesRaw[n]).Contains(".")) { // if type contains '.' then type is a frame files list
      // set input frame file
      strcpy(frFiles[n],channelTypesRaw[n]);
      strcpy(frFiles[n+nIFO],channelTypesRaw[n+nIFO]);
      continue;
    }

    // check cluster
    if((site_cluster!="CIT")&&(site_cluster!="ATLAS")) {
      cout << "cwb_eced error : gw_data_find can be used only in CIT or ATLAS cluster" << endl;
      gSystem->Exit(1);
    }

    // eced is enabled only for L1,H1,H2,V1,G1 detectors
    // gw_data_find is define only for such detectors
    if(strcmp(ifo[n],"L1")&&strcmp(ifo[n],"H1")&&strcmp(ifo[n],"H2")&&strcmp(ifo[n],"V1")&&strcmp(ifo[n],"G1")) {
      cout << "cwb_eced - Error : declared ifo " << ifo[n] 
           << " not enabled for eced (must be L1 or H1 or H2 or V1 or G1)" << endl;
      gSystem->Exit(1);
    }

    TString ECED_IFO = ifo[n];

    strcpy(ldf_command[n],"gw_data_find");

    // select the observatory
    char observatory[4]="";
    if(ECED_IFO=="L1") strcpy(observatory,"L");
    if(ECED_IFO=="H1") strcpy(observatory,"H");
    if(ECED_IFO=="H2") strcpy(observatory,"H");
    if(ECED_IFO=="V1") strcpy(observatory,"V");
    if(ECED_IFO=="G1") strcpy(observatory,"G");

    sprintf(ldf_command[n],"%s --observatory=%s",ldf_command[n],observatory);

    // select type
    if(TString(channelNamesRaw[n])=="") {
      cout << "cwb_eced - Error : strain channel name for ifo " << ifo[n] 
           << " not defined, check parameter channelNamesRaw[" << n << "] in config/user_parameters.C" << endl;
      gSystem->Exit(1);
    }

    // create temporary files to store gw_data_find output
    TString base = "eced";
    FILE* fp = gSystem->TempFileName(base);
    TString tmpFile =  CWB::Toolbox::getFileName(fp);
    fclose(fp);

    sprintf(ldf_command[n],"%s --type=%s",ldf_command[n],channelTypesRaw[n]);
    sprintf(ldf_command[n],"%s --gps-start-time=%f --gps-end-time=%f",
            ldf_command[n],gps_start-dataShift[n],gps_stop-dataShift[n]);
    sprintf(cmd,"%s --url-type=file > %s",ldf_command[n],tmpFile.Data());
    sprintf(frFiles[n],"%s/%s.frames",input_dir,ifo[n]); // set input frame file
    sprintf(ldf_command[n],"%s --url-type=file > %s",ldf_command[n],frFiles[n]);

    // -------------------
    // check if data exist 
    // -------------------

    // execute gw_data_find command
    cout << cmd << endl;
    int ret = gSystem->Exec(cmd);
    if(ret) {
      cout << "cwb_eced - gw_data_find error" << endl;
      gSystem->Exit(1);
    }

    // read size of frame file list 
    Long_t id, size, flags, modtime;
    gSystem->GetPathInfo(tmpFile, &id, &size, &flags, &modtime);

    // delete temporary file
    sprintf(cmd,"rm %s",tmpFile.Data());
    gSystem->Exec(cmd);

    // check size of frame file list 
    if(size==0) {ldf_check[n]=false;ldf_gcheck=false;} 
  }
  // check ligo_dta_find tests
  cout << endl;
  for(int n=0;n<nIFO;n++) {
    if(!ldf_check[n]) cout << ifo[n] << " cwb_eced - no data found !!!" << endl;
    else              cout << ifo[n] << " cwb_eced - data found" << endl;
  }
  if(!ldf_gcheck) {
    cout << endl << "check types : gw_data_find --show-types" << endl << endl;
    gSystem->Exit(1);
  }

  // create eced working directories
  bool check  = (cwb_inet_opts=="ced") ? true : false;
  bool bremove = (cwb_inet_opts=="ced") ? true : false;
  CWB::Toolbox::mkDir(eced_dir,check,bremove);
  CWB::Toolbox::mkDir(tmp_dir,check,bremove);
  CWB::Toolbox::mkDir(ced_dir,check,bremove);

  // create frame files
  for(int n=0;n<nIFO;n++) {
    // execute gw_data_find command
    cout << ldf_command[n] << endl;
    int ret = gSystem->Exec(ldf_command[n]);
    if(ret) {
      cout << "cwb_eced - gw_data_find error" << endl;
      gSystem->Exit(1);
    }
   
    // extract unique frame file list 
    //CWB::Toolbox::getUniqueFileList(frFiles[n],frFiles[n]);
  }

  // check size of frame file list 
  for(int n=0;n<nIFO;n++) {
    Long_t id, size, flags, modtime;
    gSystem->GetPathInfo(frFiles[n], &id, &size, &flags, &modtime);
    if(size==0) {
      cout << "cwb_eced - Error : frame file list \"" << frFiles[n] << "\" with size=0" << endl;
      gSystem->Exit(1);
    }
  }

  // create dummy user_parameters.C
  char uparameter[1024];sprintf(uparameter,"%s/user_parameters.C",tmp_dir);
  if(eced_cfg!="") {
    // used only for book-keeping
    gSystem->Exec(TString("cp ")+eced_cfg+TString(" ")+TString(uparameter));
  } else {
    // dummy user_parameters.C (only to avoid config::check error)
    gSystem->Exec(TString("touch ")+TString(uparameter));
  }
  gSystem->Setenv("CWB_UPARAMETERS_FILE",uparameter);

  // set ced output dir to be ced_dir
  gSystem->Setenv("CWB_CED_DIR",ced_dir);

  if(cwb_inet_opts=="ced") {
    // create www symbolic link
    sprintf(cmd,"mkdir -p %s/../ceds",www_dir);
    cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"ln -s %s/%s %s/../ceds/%s",gSystem->WorkingDirectory(),ced_dir,www_dir,eced_dir);
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  // create cat1 file
  for(int n=0;n<nIFO;n++) {

    ofstream out;
    char cat1_file[1024];sprintf(cat1_file,"%s/%s_cat1.in",input_dir,ifo[n]);
    out.open(cat1_file,ios::out);
    if (!out.good()) {cout << "cwb_eced - Error : Error Opening File : " << cat1_file << endl;exit(1);}
    out << gps_start << " " << gps_stop << endl;
    out.close();

    strcpy(DQF[nDQF].ifo, ifo[n]);
    sprintf(DQF[nDQF].file, "%s",cat1_file);
    DQF[nDQF].cat    = CWB_CAT1;
    DQF[nDQF].shift  = 0.;
    DQF[nDQF].invert = false;
    DQF[nDQF].c4     = false;
    nDQF++;
  }
}
