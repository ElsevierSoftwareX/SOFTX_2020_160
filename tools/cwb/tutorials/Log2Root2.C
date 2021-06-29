// 
// convert BurstMDC ascii log file to cWB root file
// Author : Gabriele Vedovato                     

{

  #define LOGFILE "BurstMDC-RD2_G1V1_AS-Log.txt"
  #define N_IFO 2


  TString ifo[N_IFO] = {"G1","V1"};
 
  CWB::mdc MDC(N_IFO,ifo);

  // write log root
  vector<mdcpar> par;
  TString rootFile = LOGFILE;
  rootFile.ReplaceAll(".txt",".root");
  cout << "Input  Log Root File : " << LOGFILE << endl;
  cout << "Output Log Root File : " << rootFile << endl;

  MDC.SetSkyDistribution(MDC_LOGFILE,LOGFILE,par);

  MDC.DumpLog(rootFile);

  exit(0);
}
