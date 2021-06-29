//
// Create DAQ file for bico condor jobs
// Author : Gabriele Vedovato

{
  //#define CHANNEL_LIST "input/H1_PEM_Channels.txt"
  #define CHANNEL_LIST "input/H-H1_RDS_R_L1_coh.txt"
  #define DATA_LABEL "TestBico_H1_RDS_R_L1_4"


  ifstream in;
  in.open(CHANNEL_LIST);
  if(!in.good()) {cout << "Error Opening File : " << CHANNEL_LIST << endl;exit(1);}
  char ofile[256];
  sprintf(ofile,"%s/%s.dag",gSystem->WorkingDirectory(),DATA_LABEL);
  ofstream out;
  out.open(ofile,ios::out);

  int id=1;
  char chname[1024];
  while(1) {
    in >> chname;
    if (!in.good()) break;
    char ostring[256];
    sprintf(ostring,"JOB A%i %s.sub",id,DATA_LABEL);
    out << ostring << endl;
    sprintf(ostring,"VARS A%i PID=\"%i\" BICO_CHNAME=\"%s\"",id,id,chname);
    out << ostring << endl;
    sprintf(ostring,"RETRY A%i 3000",id);
    id++;
  }
  out.close();
  in.close();

  exit(0);
}
