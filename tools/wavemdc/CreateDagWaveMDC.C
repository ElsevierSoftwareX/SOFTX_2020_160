int
CreateDagWaveMDC(TString condor_dir, TString label, int jobmin, int jobmax, int jobstep=1) {

  char ofile[256];
  sprintf(ofile,"%s/%s.dag",condor_dir.Data(),label.Data());
  ofstream out;
  out.open(ofile,ios::out);

  int jID=1;
  for(int n=jobmin;n<=jobmax;n+=jobstep) {
    char ostring[256];
    sprintf(ostring,"JOB A%i %s/%s.sub",jID,condor_dir.Data(),label.Data());
    out << ostring << endl;
    int start=n;
    int stop = (n+jobstep-1)<=jobmax ? n+jobstep-1 : jobmax;
    sprintf(ostring,"VARS A%i PID=\"%i\" WMDC_SEG_START=\"%i\" WMDC_SEG_STOP=\"%i\" ",jID,jID,start,stop);
    out << ostring << endl;
    sprintf(ostring,"RETRY A%i 3000",jID);
    out << ostring << endl;
    jID++;
  }
  out.close();

  return 0;
}
