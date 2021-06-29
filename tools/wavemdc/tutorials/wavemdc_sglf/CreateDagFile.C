#define FIRST_JOB_ID 210
#define LAST_JOB_ID 310
#define FILE_PER_JOB 10

#define DATA_LABEL "TestBurstMDC"


{
  int start=FIRST_JOB_ID;
  int end=LAST_JOB_ID;
  char ofile[256];
  sprintf(ofile,"%s/%s.dag",gSystem->WorkingDirectory(),DATA_LABEL);
  ofstream out;
  out.open(ofile,ios::out);
  int jobId=1; 
  for (int i=start;i<end+1;i+=FILE_PER_JOB) {
    char ostring[256];
    sprintf(ostring,"JOB A%i %s.sub",jobId,DATA_LABEL);
    out << ostring << endl;
    sprintf(ostring,"VARS A%i PID=\"%i\" MDC_SEG_START=\"%i\" MDC_SEG_STOP=\"%i\" ",jobId,jobId,i,i+FILE_PER_JOB-1);
    out << ostring << endl;
    sprintf(ostring,"RETRY A%i 3000",jobId);
    out << ostring << endl;
    jobId++;
  }
  out.close();
  exit(0);
}
