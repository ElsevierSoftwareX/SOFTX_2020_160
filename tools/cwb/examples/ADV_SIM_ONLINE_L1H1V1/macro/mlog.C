{
  #include <vector>

  CWB::Toolbox TB;

  network* net=MDC.GetNetwork();
  int nIFO=net->ifoListSize();
  char ifoLabel[64]="";
  for(int i=0;i<nIFO;i++) sprintf(ifoLabel,"%s%s",ifoLabel,net->ifoName[i]);

  TString logFile = TString("Log-")+TString(ifoLabel)+"-"+frLabel+".txt";

  bool first=true;

  vector<TString> logList;

  vector<TString> dirList = TB.getFileListFromDir(frDir, "logs");
  for(int i=0;i<dirList.size();i++) {
    cout << dirList[i] << endl;
    vector<TString> fileList = TB.getFileListFromDir(dirList[i], ".txt");
    for(int j=0;j<fileList.size();j++) {
      logList.push_back(fileList[j].Data());
      //cout << j << " " << fileList[j].Data() << endl;
    }
  }

  TString lstFile = frDir+"/Log-"+TString(ifoLabel)+"-"+frLabel+".lst";
  ofstream out;
  out.open(lstFile.Data(),ios::out);
  if (!out.good()) {cout << "LogWaveMDC.C - Error Opening File : " << lstFile.Data() << endl;gSystem->Exit(1);}
  for(int i=0;i<logList.size();i++) out << logList[i].Data() << endl;
  out.close();

  TString sortFile = lstFile;
  sortFile.ReplaceAll(".lst",".sort");

  char cmd[256];
  sprintf(cmd,"sort %s > %s",lstFile.Data(),sortFile.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);

  // write merged log
  ifstream in;
  in.open(sortFile.Data(),ios::in);
  if (!in.good()) {cout << "LogWaveMDC.C - Error Opening Sorted File : " << sortFile.Data() << endl;gSystem->Exit(1);}

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

  // write file frame list
  TString frlFile = lstFile;
  frlFile.ReplaceAll(".lst",".frl");
  sprintf(cmd,"ls %s/%s/*/*.gwf > %s",gSystem->WorkingDirectory(),frDir.Data(),frlFile.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);

  // write log root
  TString rootFile = logFile;
  rootFile.ReplaceAll(".txt",".root");
  cout << TString(frDir+"/"+rootFile).Data() << endl;
  MDC.SetSkyDistribution(MDC_LOGFILE,frDir+"/"+logFile,0);
//  MDC.DumpLog(frDir+"/"+rootFile,frLabel);   // commented because it crash. why?

  // write log header
  TString logHeaderFile = logFile+".header";
  MDC.DumpLogHeader(frDir+"/"+logHeaderFile,frLabel);

  // merge log and log header
  sprintf(cmd,"cat %s/%s >>  %s/%s",frDir.Data(),logFile.Data(),frDir.Data(),logHeaderFile.Data());
  //cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"rm %s/%s",frDir.Data(),logFile.Data());
  //cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"mv %s/%s %s/%s",frDir.Data(),logHeaderFile.Data(),frDir.Data(),logFile.Data());
  //cout << cmd << endl;
  gSystem->Exec(cmd);

  // dump waveforms
  sprintf(cmd,"mkdir -p %s/Waveforms",frDir.Data());
  gSystem->Exec(cmd);
  cout << "Dump waveforms ..." << endl;
  //MDC.PrintWaveformsList();
  MDC.Dump(frDir+"/Waveforms");

  gSystem->Exit(0);
}
