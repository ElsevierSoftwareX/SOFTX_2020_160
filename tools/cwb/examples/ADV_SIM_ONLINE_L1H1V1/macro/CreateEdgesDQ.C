#define SCRATCH 2

{

  char jobListFile[256];
  sprintf(jobListFile,"%s/%s.sjob",dump_dir,data_label);

  char dqEdgesFile[256];
  sprintf(dqEdgesFile,"%s/edges.in",input_dir);

  ifstream in;
  in.open(jobListFile);
  if(!in.good()) {cout << "Error Opening File : " << jobListFile << endl;gSystem->Exit(1);}

  ofstream out;
  out.open(dqEdgesFile,ios::out);
  cout << "Create DQ edges -> " << dqEdgesFile << endl;

  int start,stop;
  while(1) {
    in >> start >> stop;
    if (!in.good()) break;
    //cout << start << " " << stop << endl;
    out << start+SCRATCH << "\t" << stop-SCRATCH << endl;
  }

  in.close();
  out.close();

  gSystem->Exit(0);
}
