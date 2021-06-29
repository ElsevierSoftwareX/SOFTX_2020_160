{

  CWB::Toolbox TB;

  TString frDir = "frames";
  TString label = "L1H1V1-TestBurstMDC";
  TString logFile = label+"-Log.txt";

  bool first=true;

  vector<TString> logList;

  vector<TString> dirList = TB.getFileListFromDir(frDir, label);
  for(int i=0;i<dirList.size();i++) {
    cout << dirList[i] << endl;
    vector<TString> fileList = TB.getFileListFromDir(dirList[i], ".txt");
    for(int j=0;j<fileList.size();j++) {
      logList.push_back(fileList[j].Data());
      //cout << j << " " << fileList[j].Data() << endl;
    }
  }

  TString lstFile = frDir+"/"+label+"-Log.lst";
  ofstream out;
  out.open(lstFile.Data(),ios::out);
  if (!out.good()) {cout << "MergeLogFiles.C - Error Opening File : " << lstFile.Data() << endl;exit(1);}
  for(int i=0;i<logList.size();i++) out << logList[i].Data() << endl;
  out.close();

  TString sortFile = lstFile;
  sortFile.ReplaceAll(".lst",".sort");

  char cmd[256];
  sprintf(cmd,"sort %s > %s",lstFile.Data(),sortFile.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);


  ifstream in;
  in.open(sortFile.Data(),ios::in);
  if (!in.good()) {cout << "MergeLogFiles.C - Error Opening Sorted File : " << sortFile.Data() << endl;exit(1);}

  char log[1024];
  while(true) {
    in.getline(log,1024);
    if (!in.good()) break;
    char cmd[256];
    if(first) {
      sprintf(cmd,"cat %s >  %s/%s",log,frDir.Data(),logFile.Data());
      first=false;
    } else {
      sprintf(cmd,"cat %s >> %s/%s",log,frDir.Data(),logFile.Data());
    }
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  in.close();

  sprintf(cmd,"rm %s %s",lstFile.Data(),sortFile.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);

  TString frlFile = lstFile;
  frlFile.ReplaceAll(".lst",".frl");
  sprintf(cmd,"ls %s/%s/*/*.gwf > %s",gSystem->WorkingDirectory(),frDir.Data(),frlFile.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);

  exit(0);

