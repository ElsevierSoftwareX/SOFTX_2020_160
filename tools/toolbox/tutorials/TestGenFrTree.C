//
// Test Create/Sort FrTree
// Author : Gabriele Vedovato

#define LIST_FILE_NAME "SEGMENTS/lists/GHLTV-HBRST14_S6D_R1.frames"
#define ROOT_FILE_NAME "root/GHLTV-HBRST14_S6D_R1.root"
#define SORT_FILE_NAME "root/GHLTV-HBRST14_S6D_R1-Sorted.root"

{

  cwbtb tb;
  int nfiles=0;
  nfiles=tb.frl2FrTree(LIST_FILE_NAME,ROOT_FILE_NAME);
  cout << "nfiles : " << nfiles << endl;
  nfiles=tb.sortFrTree(ROOT_FILE_NAME,SORT_FILE_NAME);
  cout << "nfiles : " << nfiles << endl;

  exit(0);
}
