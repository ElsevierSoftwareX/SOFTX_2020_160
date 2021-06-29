
#include <vector>

void RemoveSuperclusterFiles() {

  CWB::Toolbox TB;

  char output_dir[1024];
  sprintf(output_dir,"%s/output", gSystem->WorkingDirectory());

  cout << endl << "output dir : " << output_dir << endl;
  bool answer = CWB::Toolbox::question("do you want to remove all supercluster files from the output dir ? ");
  if(!answer) exit(1);

  cout << "Starting reading output directory ..." << endl;
  vector<TString> fileList = TB.getFileListFromDir(output_dir,".root","supercluster_","",true);
  int nfile = fileList.size();
  for(int n=0;n<nfile;n++) {
    if(n==0) cout << n << " " << fileList[n].Data()<< endl;
//    if(n<1000) cout << n << " " << fileList[n].Data()<< endl;

    if (n%1000==0) cout << "cwb_condor benchmark - " << n << "/" << fileList.size() << " files" << endl;
    gSystem->Unlink(fileList[n].Data());
  }

  exit(0);
}
