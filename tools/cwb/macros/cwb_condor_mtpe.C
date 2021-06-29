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

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  if(TString(condor_tag)=="") {
    cout << endl;
    cout << "cwb_condor_mtpe.C : Error - the accounting_group is not defined !!!" << endl;
    cout << "The accounting_group must be defined in the user_parameters.C file" << endl;
    cout << "See the following link:" << endl;
    cout <<" https://ldas-gridmon.ligo.caltech.edu/accounting/condor_groups/determine_condor_account_group.html" << endl;
    cout << "Examples : " << endl;
    cout << "strcpy(condor_tag,\"ligo.dev.o2.burst.allsky.cwboffline\");" << endl;
    cout << "strcpy(condor_tag,\"ligo.prod.o2.burst.allsky.cwboffline\");" << endl;
    cout << "If you don't need it set : strcpy(condor_tag,\"disabled\");" << endl << endl;
    exit(1);
  }
  if(TString(condor_tag)=="disabled") strcpy(condor_tag,"");

  if(!(simulation==4) && !(simulation==0)) {
    cout << endl << "cwb_condor_mtpe Error : allowed only with simulation=0/4 !!!" << endl << endl;
    gSystem->Exit(1);
  }

  bool singleDetector=false;
  if(nIFO==1) { // Single Detector Mode
    CWB::config config;
    config.Import();
    config.SetSingleDetectorMode();
    config.Export();
    singleDetector=true;
  }
  if(nIFO==2) { // 2 detectors with the same name "same detector -> nIFO=1"
    if(TString(ifo[0])==TString(ifo[1])) singleDetector=true;
  }

  // get PE jobid
  int cwb_condor_mtpe_jobid=-1; 
  if(gSystem->Getenv("CWB_CONDOR_MTPE_JOBID")==NULL) {
    cout << endl << "cwb_condor_mtpe Error : environment CWB_CONDOR_MTPE_JOBID is not defined!!!" << endl << endl;
    gSystem->Exit(1);
  } else {
    if(TString(gSystem->Getenv("CWB_CONDOR_MTPE_JOBID")).IsDigit()) {
      cwb_condor_mtpe_jobid=TString(gSystem->Getenv("CWB_CONDOR_MTPE_JOBID")).Atoi();
    } else {
      cout << endl << "cwb_condor_mtpe Error : environment CWB_CONDOR_MTPE_JOBID is not an integer number!!!" << endl << endl;
      gSystem->Exit(1);
    }
  }
  if(cwb_condor_mtpe_jobid<1) {
    cout << endl << "cwb_condor_mtpe Error : environment CWB_CONDOR_MTPE_JOBID is not an integer number >0 !!!" << endl << endl;
    gSystem->Exit(1);
  }

  // get PE trials
  int cwb_condor_mtpe_trials=-1; 
  int cwb_condor_mtpe_offset=1; 

  // for simulation=4 the trials are defined in the nfactor, factors parameters
  if(simulation==4 && nfactor>1) {
    cwb_condor_mtpe_trials = factors[0]+nfactor-1;
    cwb_condor_mtpe_offset = factors[0];
    if(cwb_condor_mtpe_trials<1) {
      cout << endl << "cwb_condor_mtpe Error : (simulation=4) nfactor is not an integer number >0 !!!" << endl << endl;
      gSystem->Exit(1);
    }
    if(cwb_condor_mtpe_offset<1) {
      cout << endl << "cwb_condor_mtpe Error : (simulation=4) factors[0] is not an integer number >0 !!!" << endl << endl;
      gSystem->Exit(1);
    }
  }

  // for simulation=0 the trials are defined in the PE parPlugin
  if((simulation==0)||(simulation==4 && nfactor==1)) {
    bool multitask=false;
    TString pe_options = parPlugin;
    if(pe_options.CompareTo("")!=0) {
      cout << pe_options << endl;
      if(!pe_options.Contains("--")) {  
  
        TObjArray* token = TString(pe_options).Tokenize(TString(' '));

        for(int j=0;j<token->GetEntries();j++) {

          TObjString* tok = (TObjString*)token->At(j);
          TString stok = tok->GetString();

          if(stok.Contains("pe_multitask=")) {
            TString pe_multitask=stok;
            pe_multitask.Remove(0,pe_multitask.Last('=')+1);
            if(pe_multitask=="true")  multitask=true;
            if(pe_multitask=="false") multitask=false;
          }

          if(stok.Contains("pe_trials=") || stok.Contains("pe_retry=")) {
            TString pe_trials=stok;
            pe_trials.Remove(0,pe_trials.Last('=')+1);
            if(pe_trials.IsDigit()) cwb_condor_mtpe_trials=pe_trials.Atoi();
          }
        }
      }
    }
    if(!multitask) {
      cout << endl << "cwb_condor_mtpe Error : (simulation=0) pe_multitask must be enabled !!!" << endl << endl;
      gSystem->Exit(1);
    }
    if(cwb_condor_mtpe_trials<1 || cwb_condor_mtpe_trials>99) {
      cout << endl << "cwb_condor_mtpe Error : (simulation=0) pe_trials in PE parPlugin is not an integer number >0 && <100 !!!" << endl << endl;
      gSystem->Exit(1);
    }
  }

  char condor_label[1024];
  sprintf(condor_label,"%s.mtpe",data_label);

  // create condor sub file

  char ofile_condor_sub[1024];
  sprintf(ofile_condor_sub,"%s/%s.sub",condor_dir,condor_label);

  FILE *fP=NULL;
  if((fP = fopen(ofile_condor_sub, "w")) == NULL) {
    cout << endl << "cwb_condor_mtpe.C : Error - cannot open file " << ofile_condor_sub << endl << endl;
    exit(1);
  }
  cout << endl << ofile_condor_sub << endl;

  char full_condor_dir[1024];
  char full_condor_out_dir[1024];
  char full_condor_err_dir[1024];

  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);
  sprintf(full_condor_out_dir,"%s/%s",work_dir,log_dir);
  sprintf(full_condor_err_dir,"%s/%s",work_dir,log_dir);

  fprintf(fP,"universe = vanilla\n");
  fprintf(fP,"getenv = true\n");
  fprintf(fP,"priority = $(PRI)\n");
  fprintf(fP,"on_exit_hold = ( ExitCode != 0 )\n");
  fprintf(fP,"request_memory = 2000\n");
  fprintf(fP,"executable = loudest.sh\n");
  fprintf(fP,"job_machine_attrs = Machine\n");
  fprintf(fP,"job_machine_attrs_history_length = 5\n");
  fprintf(fP,"requirements = target.machine =!= MachineAttrMachine1 && target.machine =!= MachineAttrMachine2 && target.machine =!= MachineAttrMachine3 && target.machine =!= MachineAttrMachine4 && target.machine =!= MachineAttrMachine5\n");
  fprintf(fP,"environment = CWB_JOBID=$(PID);CWB_UFILE=$(CWB_UFILE);CWB_STAGE=$(CWB_STAGE);CWB_MDC_FACTOR=$(CWB_MDC_FACTOR);CWB_CED_DIR=$(CWB_CED_DIR)\n");
  if(TString(condor_tag)!="") fprintf(fP,"accounting_group = %s\n",condor_tag);
  fprintf(fP,"output = %s/$(PID)_$(CWB_MDC_FACTOR)_%s.out\n",full_condor_out_dir,condor_label);
  fprintf(fP,"error = %s/$(PID)_$(CWB_MDC_FACTOR)_%s.err\n",full_condor_err_dir,condor_label);
  fprintf(fP,"log = %s/%s.log\n",condor_log,condor_label);
  fprintf(fP,"notification = never\n");
  fprintf(fP,"rank=memory\n");
  fprintf(fP,"queue\n");

  fclose(fP);

  // create condor dag file

  char ofile_condor_dag[1024];
  sprintf(ofile_condor_dag,"%s/%s.dag",condor_dir,condor_label);

  ofstream out;
  out.open(ofile_condor_dag,ios::out);
  if (!out.good()) {cout << endl << "cwb_condor_mtpe Error Opening File : " << ofile_condor_dag << endl << endl;exit(1);}
  cout << ofile_condor_dag << endl << endl;

  TString cwb_uparameters_file=TString(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  for(int n=cwb_condor_mtpe_offset; n<=cwb_condor_mtpe_trials; n++) {
    char ostring[1024*16];
    sprintf(ostring,"JOB A%i_%d %s/%s.sub",cwb_condor_mtpe_jobid,n,full_condor_dir,condor_label);
    out << ostring << endl;
    sprintf(ostring,"VARS A%i_%d PID=\"%i\" CWB_UFILE=\"%s\" CWB_STAGE=\"CWB_STAGE_FULL\" CWB_MDC_FACTOR=\"%d\" CWB_CED_DIR=\"output\"",
            cwb_condor_mtpe_jobid,n,cwb_condor_mtpe_jobid,cwb_uparameters_file.Data(),n);
    out << ostring << endl;
    sprintf(ostring,"RETRY A%i_%d 3000",cwb_condor_mtpe_jobid,n);
    out << ostring << endl;
  }

  if((simulation==0)||(simulation==4 && nfactor==1)) {
    char ostring[1024*16];
    out << "PARENT ";
    for(int n=1; n<cwb_condor_mtpe_trials; n++) {
      sprintf(ostring,"A%i_%d ",cwb_condor_mtpe_jobid,n);
      out << ostring;
    }
    sprintf(ostring,"CHILD A%i_%d",cwb_condor_mtpe_jobid,cwb_condor_mtpe_trials);
    out << ostring << endl;
    out.close();
  }

  exit(0);
}

