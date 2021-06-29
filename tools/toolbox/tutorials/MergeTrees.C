// Show how to merge two wave/mdc root files
// The input root files are in the output directory
//   the file merged are "wave/mdc*.M1.root"
// The output merged files are saved to merge directory
//   the output merged files are renamed to wave/mdc_OLABEL.root
// Author : Gabriele Vedovato

#define WAVE
#define OLABEL  "S6A-VSR2_HLV_SIM_EOBNRv2_NSNS_2G_merge_dbg.M1"

{
  CWB::Toolbox TB;

#ifdef WAVE
  TString treeName = "waveburst";
#else
  TString treeName = "mdc";
#endif
  TString dir_name = "output";
  TString merge_dir = "merge";
#ifdef WAVE
  TString beginsWith = "wave";
#else
  TString beginsWith = "mdc";
#endif
  TString endsWith = ".M1.root";
  TString label = OLABEL;

  vector<TString> fileList = TB.getFileListFromDir(dir_name,endsWith,beginsWith);
  for(int i=0;i<fileList.size();i++) {
    cout << i << " " << fileList[i].Data() << endl;
  }
  TB.mergeTrees(fileList, treeName, merge_dir, label);

  TString ifileName = merge_dir+"/merge_"+label+".root";
#ifdef WAVE
  TString ofileName = merge_dir+"/wave_"+label+".root";
#else
  TString ofileName = merge_dir+"/mdc_"+label+".root";
#endif

  char cmd[1024];
  sprintf(cmd,"mv %s %s",ifileName.Data(),ofileName.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);

  exit(0);
}
