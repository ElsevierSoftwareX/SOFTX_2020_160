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


// interactive cwb pipeline (1G & 2G)

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  //TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_NETC_FILE"));

  TString cwb_batch="false";
  if(gSystem->Getenv("CWB_BATCH")!=NULL) {
    cwb_batch=TString(gSystem->Getenv("CWB_BATCH"));
  }

  if(cwb_batch=="false") strcpy(nodedir,tmp_dir); // tmp_dir is used in interactive mode
  strcpy(output_dir,data_dir);  // default

  int cwb_jobid=0;
  if(gSystem->Getenv("CWB_JOBID")==NULL) {
    cout << "Error : environment CWB_JOBID is not defined!!!" << endl;exit(1);
  } else {
    if(TString(gSystem->Getenv("CWB_JOBID")).IsDigit()) {
      cwb_jobid=TString(gSystem->Getenv("CWB_JOBID")).Atoi();
    } else {
      cout << "Error : environment CWB_JOBID is not defined!!!" << endl;exit(1);
    }
  }

  if(gSystem->Getenv("CWB_CED_DIR")!=NULL) {
    TString cwb_ced_dir=TString(gSystem->Getenv("CWB_CED_DIR"));
    if(cwb_ced_dir.Sizeof()>1) { 
      TB.checkFile(cwb_ced_dir);
      strcpy(output_dir,cwb_ced_dir.Data());  // user defined
    }
  }
 
  TString cwb_inet_options=TString(gSystem->Getenv("CWB_INET_OPTIONS"));
  if(cwb_inet_options.CompareTo("")!=0) {
    TString inet_options = cwb_inet_options; 
    //inet_options.ToUpper();
    if(!inet_options.Contains("--")) {	// parameters are used only by cwb_inet

      TObjArray* token = TString(inet_options).Tokenize(TString(' '));
      for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok=="true")   cedDump = true;
        if(stok=="ced")    cedDump = true;
        if(stok=="false")  cedDump = false;
        if(stok=="root")   cedDump = true;
        if(stok.Contains("cedDump=")) {
          TString cedDump_par=stok;
          cedDump_par.Remove(0,cedDump_par.Last('=')+1);
          if(cedDump_par=="true")  cedDump=true;
          if(cedDump_par=="false") cedDump=false;
        }

        if(stok.Contains("gps=")) {
          TString gps_par=stok;
          gps_par.Remove(0,gps_par.Last('=')+1);
          if(gps_par.IsFloat()) gSystem->Setenv("CWB_GPS_EVENT",gps_par);
        }
        if(stok.Contains("iwindow=")) {
          TString iwindow_par=stok;
          iwindow_par.Remove(0,iwindow_par.Last('=')+1);
          if(iwindow_par.IsFloat()) gap=iwindow_par.Atof();
        }
        if(stok.Contains("netCC=")) {
          TString netCC_par=stok;
          netCC_par.Remove(0,netCC_par.Last('=')+1);
          if(netCC_par.IsFloat()) netCC=netCC_par.Atof();
        }
        if(stok.Contains("Acore=")) {
          TString Acore_par=stok;
          Acore_par.Remove(0,Acore_par.Last('=')+1);
          if(Acore_par.IsFloat()) Acore=Acore_par.Atof();
        }
        if(stok.Contains("dump=")) {
          TString dump_par=stok;
          dump_par.Remove(0,dump_par.Last('=')+1);
          if(dump_par=="true")  dump=true;
          if(dump_par=="false") dump=false;
        }
        if(stok.Contains("plugin=")) {
          TString plugin_par=stok;
          plugin_par.Remove(0,plugin_par.Last('=')+1);
          if(plugin_par.EndsWith(".C"))  plugin=TMacro(plugin_par);
        }
        if(stok.Contains("search=")) {
          TString search_par=stok;
          search_par.Remove(0,search_par.Last('=')+1);
          if(search_par=="r")  SEARCH()='r';
          if(search_par=="i")  SEARCH()='i';
          if(search_par=="p")  SEARCH()='p';
          if(search_par=="l")  SEARCH()='l';
          if(search_par=="c")  SEARCH()='c';
          if(search_par=="e")  SEARCH()='e';
          if(search_par=="s")  SEARCH()='s';
          if(search_par=="g")  SEARCH()='g';
          if(search_par=="b")  SEARCH()='b';
        }
        if(stok.Contains("optim=")) {
          TString optim_par=stok;
          optim_par.Remove(0,optim_par.Last('=')+1);
          if(optim_par=="true")  optim=true;
          if(optim_par=="false") optim=false;
        }

      }

    } else {				// parameters are used only by CWB_Plugin_cwb_inet.C

      // get CWB_Plugin_cwb_inet.C
      TMacro cwb_inet_plugin = gSystem->ExpandPathName("$HOME_WAT/tools/cwb/plugins/CWB_Plugin_cwb_inet.C"); 
      // creates a temporary file name
      char fpluginName[1024]; 
      unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
      sprintf(fpluginName,"%s/CWB_Plugin_cwb_inet_%s_%d_job%d.XXXXXX",
              tmp_dir,data_label,Pid,cwb_jobid);
      CWB::Toolbox::mksTemp(fpluginName);  // create temporary file, must be deleted only at the end of run
      // write plugin, if plugin is defined then we merge cwb_inet_plugin & plugin
      ofstream out;
      out.open(fpluginName,ios::out);
      if(!out.good()) {cout << "cwb_inet.C - Error : Opening File : " << fpluginName << endl;gSystem->Exit(1);}
      if(TString(plugin.GetName())!="") {	
        TList* fLines = plugin.GetListOfLines();
        TObjString *obj;
        TIter next(fLines);
        out <<  "// --> BEGIN USER PLUGIN CODE" << endl;
        while ((obj = (TObjString*) next())) {
          TString line = obj->GetName();
          line.ReplaceAll("CWB_Plugin(","CWB_UserPlugin(");	// rename plugin name
          out <<  line.Data() << endl;
        }
        out <<  "// --> END USER PLUGIN CODE" << endl << endl;
      }
      // write CWB_Plugin_cwb_inet.C code
      TList* fLines = cwb_inet_plugin.GetListOfLines();
      TObjString *obj;
      TIter next(fLines);
      out <<  "// --> BEGIN CWB_INET PLUGIN CODE" << endl;
      while ((obj = (TObjString*) next())) {
        TString line = obj->GetName();
        out <<  line.Data() << endl;
        // add CWB_PLUGIN (user plugin) call after the CWB_Plugin line
        if((TString(plugin.GetName())!="")&&(line.Contains("CWB_Plugin("))) {
          out << endl;
          out << "  CWB_UserPlugin(jfile, cfg, net, x, ifo, type); // CALL USER PLUGIN CODE" << endl;
          out << endl;
        }
      }
      out <<  "// --> END CWB_INET PLUGIN CODE" << endl;
      out.close();

      if(cwb_inet_options.Contains("frdisplay")) plugin = fpluginName;
      if(cwb_inet_options.Contains("psd"))       plugin = fpluginName;
      if(cwb_inet_options.Contains("inj"))       plugin = fpluginName;
      if(cwb_inet_options.Contains("wdm"))       plugin = fpluginName;
      if(cwb_inet_options.Contains("nrms"))      plugin = fpluginName;
      if(cwb_inet_options.Contains("emax"))      plugin = fpluginName;
      if(cwb_inet_options.Contains("sparse"))    plugin = fpluginName;
      if(cwb_inet_options.Contains("ced"))       plugin = fpluginName;
      plugin.SetName(fpluginName);
    }
  }

  int gps_event=0;
  TString cwb_gps_event=TString(gSystem->Getenv("CWB_GPS_EVENT"));
  if(cwb_gps_event.CompareTo("")!=0) {
    if(!cwb_gps_event.IsFloat()) {cout<< "Error : CWB_GPS_EVENT is not a number" << endl;exit(1);}
    if(cwb_gps_event.Atoi()>0) gps_event=cwb_gps_event.Atoi();
  }

  network NET_CED;                          // network
  detector* pD_CED[NIFO_MAX];                   // pointers to detectors
  for(int n=0;n<NIFO_MAX;n++) pD_CED[n]=NULL;

  segTHR=0;				    // remove the cat2 check
  if(simulation) {
    float mdc_factor=0;
    TString cwb_mdc_factor=TString(gSystem->Getenv("CWB_MDC_FACTOR"));
    if(cwb_mdc_factor.CompareTo("")!=0) {
      if(!cwb_mdc_factor.IsFloat()) {cout<< "Error : CWB_MDC_FACTOR is not a number" << endl;exit(1);}
      if(cwb_mdc_factor.Atof()>0) mdc_factor=cwb_mdc_factor.Atof();
    }
    if(mdc_factor!=0) {
      nfactor=1;
      factors[0]=mdc_factor;
    }
  } else {
    int job_lag=-1;	  // if <0 all lags defined in user parameters are selected
    TString cwb_job_lag=TString(gSystem->Getenv("CWB_JOB_LAG"));
    if(cwb_job_lag.CompareTo("")!=0) {
      if(!cwb_job_lag.IsFloat()) {cout<< "Error : CWB_JOB_LAG is not a number" << endl;exit(1);}
      job_lag=cwb_job_lag.Atoi();
    }

    // select only one lag
    if(job_lag>=0) {
      mlagStep= 0;
      lagSize = 1;        // number of lags (simulation=1)
      lagOff  = job_lag;  // first lag id (lagOff=0 - include zero lag )
    } else job_lag=0;

    for(int i=0; i<nIFO; i++) {
      if(strlen(ifo[i])>0) pD_CED[i] = new detector(ifo[i]);        // built in detector
      else                 pD_CED[i] = new detector(detParms[i]);   // user define detector
    }

    for(int i=0; i<nIFO; i++) NET_CED.add(pD_CED[i]);

    // WARNING : segLen now is fixed, should be constrained by cat1  !!!

    wavearray<double> x(segLen*16384);
    x.rate(16384.);
    pD_CED[0]->TFmap=x;

    int lags_ced = NET_CED.setTimeShifts(lagSize,lagStep,lagOff,lagMax,lagFile,lagMode,lagSite);
    cout<<"lag step: "<<lagStep<<endl;
    cout<<"number of time lags_ced : " << lags_ced << endl;

    // print selected lag
    printf("%8s ","lag");
    for(int n=0; n<nIFO; n++) printf("%12.12s%2s","ifo",NET_CED.getifo(n)->Name);
    printf("\n");
    printf("%8d ",job_lag);
    for(int n=0; n<nIFO; n++) printf("%14.5f",pD_CED[n]->lagShift.data[0]);
    printf("\n");

    pD_CED[0]->TFmap.resize(0);
  }

  if(gps_event==0) return; 

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}

  if(gps_event>0) {
    for(int n=0; n<nIFO; n++) {

      strcpy(DQF[nDQF].ifo, ifo[n]);
      sprintf(DQF[nDQF].file, "%s/%s_%s.gps_%d",tmp_dir,ifo[n],data_label,int(gps_event));
      DQF[nDQF].cat    = CWB_CAT2;
      DQF[nDQF].shift  = 0.;
      DQF[nDQF].invert = false;
      DQF[nDQF].c4     = true;
      nDQF++;

      cout << DQF[nDQF-1].file << endl; 

      ofstream out;
      out.open(DQF[nDQF-1].file,ios::out);
      cout << "Write file : " << DQF[nDQF-1].file << endl;
      if (!out.good()) {cout << "Error Opening File : " << DQF[nDQF-1].file << endl;exit(1);}
      out.precision(14);
      if(simulation) {
        int istart = int(gps_event)-gap;
        int istop  = int(gps_event)+gap;
        out << "1 " << istart << " " << istop << " " << 2*gap << endl;
      } else {
        int jobSegLen = 0;
        if(slagSize>0) {   // SLAG Segments
          if(slagSize==1) {
            if((slagMin!=0)||(slagMax!=0)||(slagOff!=0)||(slagFile!=NULL)) {
              cout << "cwb_inet gps_event>1 is not implemented with : !!!" << endl;
              cout << "if slagSize=1 -> slagMin=slagMax=slagOff=0, slagFile=NULL "<<endl;
              exit(1);
            }
          } else {
            cout << "cwb_inet gps_event>1 is not implemented with slagSize>1 !!!" << endl;
            exit(1);
          }
        } else {           // Standard Segments 
          CWB_CAT dqcat = CWB_CAT1;
          vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, dqcat);
          vector<waveSegment> jobList=TB.getJobList(cat1List, segLen, segMLS, segEdge);        // get  standard job list
	  if(jobList.size()==0) {  
            cout << endl << "cwb_inet standard job list size = 0 !!!" << endl << endl;
            exit(1);
          }
          jobSegLen = jobList[cwb_jobid-1].stop-jobList[cwb_jobid-1].start; 
        }
  
        double tShift = pD_CED[n]->lagShift.data[0] - pD_CED[0]->lagShift.data[0]; 
        if(tShift<0) tShift=jobSegLen-fabs(tShift); 
  
        int istart = int(gps_event+tShift)-gap;
        int istop  = int(gps_event+tShift)+gap;
        out << n+1 << " " << istart << " " << istop << " " << 2*gap << endl;
      }
      out.close();
    }
  }

  for(int n=0; n<nIFO; n++) if(pD_CED[n]!=NULL) delete pD_CED[n];
}
