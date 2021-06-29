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
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  if(TString(condor_tag)=="") {
    cout << endl;
    cout << "cwb_condor_create_ced.C : Error - the accounting_group is not defined !!!" << endl;
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

  // creates ced dir for loudest events
  TString pp_ced_dir = TString(pp_dir)+TString("/ced");
  CWB::Toolbox::mkDir(pp_ced_dir,!pp_batch);

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

  char events_sorted[256];
  sprintf(events_sorted,"%s/events_sorted.txt",netdir);
  CWB::Toolbox::checkFile(events_sorted);

  char pm[8];
  char c3[8];
  float icc, icc2, icc3, irho, iacor, ilag, islag, ilik, ipen, icHH, ivHH, ivED;
  int ifreq, ilow, ihigh;
  float idur;
  int isize, irate, irun;
  float phi, theta, psi;

  double* itime = new double[nIFO];
  double* iSNR  = new double[nIFO];
  double* ihrss = new double[nIFO];

/*
  char ofile_name[256];
  sprintf(ofile_name,"%s/event_ced_parameters.txt",netdir);
  cout << ofile_name << endl;

  ofstream out;
  out.open(ofile_name,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << ofile_name << endl;exit(1);}

  int ievt=0;
  ifstream f_ev(events_sorted);
  //cout << "jobid gps_evt[0=(full job)/sec] ced_dump[false/true] mdc_factor lag[0=(full lags)/#lag]" << endl;
  while(ievt<pp_max_nloudest_list) {
    f_ev>>pm>>c3>>irho>>icc>>iacor>>ilag>>islag>>ilik>>ipen>>icHH>>ifreq>>ilow>>ihigh>>idur>>isize>>irate>>irun;
    for(int i=0;i<nIFO;i++) f_ev>>itime[i];
    for(int i=0;i<nIFO;i++) f_ev>>iSNR[i];
    for(int i=0;i<nIFO;i++) f_ev>>ihrss[i];
    f_ev>>phi>>theta>>psi;
    ievt++;
    if (!f_ev.good()) break;
    //cout << irun << " " << 0 << " " << "true" << " " << 0 << " " << ilag << " " << ced_dir << endl; 
    out << irun << " " << 0 << " " << "true" << " " << 0 << " " << ilag << " " << ced_dir << endl; 
  }
  cout << endl;
  f_ev.close();
  out.close();
*/

  char condor_label[1024];
  sprintf(condor_label,"%s.%s.R%s.ced",data_label,cwb_merge_label.Data(),pp_label.Data());

  // create condor sub file

  char ofile_condor_sub[1024];
  sprintf(ofile_condor_sub,"%s/%s.sub",condor_dir,condor_label);

  FILE *fP=NULL;
  if((fP = fopen(ofile_condor_sub, "w")) == NULL) {
    cout << "cwb_condor_create_ced.C : Error - cannot open file " << ofile_condor_sub << endl;
    exit(1);
  }
  cout << ofile_condor_sub << endl;

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
  fprintf(fP,"executable = ced.sh\n");
  fprintf(fP,"job_machine_attrs = Machine\n");
  fprintf(fP,"job_machine_attrs_history_length = 5\n");
  fprintf(fP,"requirements = target.machine =!= MachineAttrMachine1 && target.machine =!= MachineAttrMachine2 && target.machine =!= MachineAttrMachine3 && target.machine =!= MachineAttrMachine4 && target.machine =!= MachineAttrMachine5\n");
  fprintf(fP,"environment = CWB_JOBID=$(PID);CWB_GPS_EVENT=$(CWB_GPS_EVENT);CWB_INET_OPTIONS=$(CWB_INET_OPTIONS);CWB_MDC_FACTOR=$(CWB_MDC_FACTOR);CWB_JOB_LAG=$(CWB_JOB_LAG);CWB_CED_DIR=$(CWB_CED_DIR);CWB_BATCH=$(CWB_BATCH)\n");
  if(TString(condor_tag)!="") fprintf(fP,"accounting_group = %s\n",condor_tag);
  fprintf(fP,"output = %s/$(PID)_$(CWB_JOB_LAG)_%s.out\n",full_condor_out_dir,condor_label);
  fprintf(fP,"error = %s/$(PID)_$(CWB_JOB_LAG)_%s.err\n",full_condor_err_dir,condor_label);
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
  if (!out.good()) {cout << "Error Opening File : " << ofile_condor_dag << endl;exit(1);}
  cout << ofile_condor_dag << endl;

  int ievt=0;
  vector<TString> JTAG;
  ifstream f_ev(events_sorted);
  while(ievt<pp_max_nloudest_list) {
    f_ev>>pm>>c3>>irho>>icc>>icc2>>icc3>>iacor>>ilag>>islag>>ilik>>ipen>>icHH>>ifreq>>ilow>>ihigh>>idur>>isize>>irate>>irun;
    for(int i=0;i<nIFO;i++) f_ev>>itime[i];
    for(int i=0;i<nIFO;i++) f_ev>>iSNR[i];
    for(int i=0;i<nIFO;i++) f_ev>>ihrss[i];
    f_ev>>phi>>theta>>psi;
    if (!f_ev.good()) break;
    if (!TString(pm).Contains("+")) continue;
    ievt++;
    //cout << irun << " " << 0 << " " << "true" << " " << 0 << " " << ilag << " " << ced_dir << endl; 
    //out << irun << " " << 0 << " " << "true" << " " << 0 << " " << ilag << " " << ced_dir << endl; 
    //int jID = irun;
    int jID = ievt;
    char jtag[1024];sprintf(jtag,"%i_%i",irun,(int)ilag);
    bool bjtag=false;for(int j=0;j<JTAG.size();j++) if(JTAG[j]==jtag) bjtag=true;
    if(!bjtag) JTAG.push_back(jtag); else continue;	// skip job if already created
    char ostring[1024];
    sprintf(ostring,"JOB A%i_%s %s.sub",jID,jtag,condor_label);
    out << ostring << endl;
    int evt_gps_time = singleDetector ? itime[0] : 0;
    sprintf(ostring,"VARS A%i_%s PID=\"%i\" CWB_GPS_EVENT=\"%i\" CWB_INET_OPTIONS=\"ced\" CWB_MDC_FACTOR=\"0\" CWB_JOB_LAG=\"%i\" CWB_CED_DIR=\"%s\" CWB_BATCH=\"true\"",jID,jtag,irun,evt_gps_time,(int)ilag,pp_ced_dir.Data());
    out << ostring << endl;
    sprintf(ostring,"RETRY A%i_%s 3000",jID,jtag);
    out << ostring << endl;
  }
  out.close();
  f_ev.close();

  exit(0);
}

