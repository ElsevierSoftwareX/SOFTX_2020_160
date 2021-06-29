int
CreateSubWaveMDC(TString label, TString condor_dir, TString out_dir, TString err_dir, TString log_dir) {

  char ofile[256];
  sprintf(ofile,"%s/%s.sub",condor_dir.Data(),label.Data());

  FILE *fP=NULL;
  if((fP = fopen(ofile, "w")) == NULL) {
    cout << "CWB::Toolbox::createSubFile : Error - cannot open file " << ofile << endl;
    exit(1);
  }

  fprintf(fP,"universe = vanilla\n");
  fprintf(fP,"getenv = true\n");
  fprintf(fP,"priority = $(PRI)\n");
  fprintf(fP,"on_exit_hold = ( ExitCode != 0 )\n");
  fprintf(fP,"request_memory = 3000\n");
  fprintf(fP,"executable = WaveMDC.sh\n");
  fprintf(fP,"environment = WMDC_JOBID=$(PID);WMDC_SEG_START=$(WMDC_SEG_START);WMDC_SEG_STOP=$(WMDC_SEG_STOP)\n");
  fprintf(fP,"accounting_group = ligo.prod.o1.burst.allsky.cwboffline\n");
  fprintf(fP,"output = %s/$(PID)_%s.out\n",out_dir.Data(),label.Data());
  fprintf(fP,"error = %s/$(PID)_%s.err\n",err_dir.Data(),label.Data());
  fprintf(fP,"log = %s/%s.log\n",log_dir.Data(),label.Data());
  fprintf(fP,"notification = never\n");
  fprintf(fP,"rank=memory\n");
  fprintf(fP,"queue\n");

  fclose(fP);

  return 0;
}
