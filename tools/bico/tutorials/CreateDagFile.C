//
// Create DAG File for bico condor jobs
// Author : Gabriele Vedovato

{

  #define FIRST_JOB_ID 1
  #define LAST_JOB_ID  242

  #define DATA_LABEL "JD3_SIM_HBRST1_G1V1_JD3_G1V1_run5r"


  int start=FIRST_JOB_ID;
  int end=LAST_JOB_ID;
  char ofile[256];
  sprintf(ofile,"%s/%s.dag",gSystem->WorkingDirectory(),DATA_LABEL);
  ofstream out;
  out.open(ofile,ios::out);
  for (int i=start;i<end+1;i++) {
    char ostring[256];
    sprintf(ostring,"JOB A%i %s.sub",i,DATA_LABEL);
    out << ostring << endl;
    sprintf(ostring,"VARS A%i PID=\"%i\"",i,i);
    out << ostring << endl;
    sprintf(ostring,"RETRY A%i 3000",i);
    out << ostring << endl;
  }
  out.close();
  exit(0);
}
