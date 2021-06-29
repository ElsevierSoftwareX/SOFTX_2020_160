// merge root trees
// Author : Gabriele Vedovato

{
  #define INPUT_ROOT_FILE_DIR   "output"
  #define TREE_NAME        	"waveburst"
  #define OUTPUT_MERGE_DIR 	"merge"
  #define OFILE_NAME      	"merge.root"

  #include <vector>

  CWB::Toolbox TB;

  vector<TString> fileList = TB.getFileListFromDir(INPUT_ROOT_FILE_DIR, ".root");

  for(int i=0;i<fileList.size();i++) {
    cout << i << " " << fileList[i] << endl;
  }

  TString odir = ".";
  TString label = MERGE_LABEL;

  TB.mergeTrees(fileList, TREE_NAME, OUTPUT_MERGE_DIR, OFILE_NAME, false);

  exit(0);

}

