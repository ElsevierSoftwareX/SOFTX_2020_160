//
// List frame : start, stop, length, path 
// Author : Gabriele Vedovato


{
  #define FRLIST_FILE "L1_NINJA2_G1000176_EARLY_RECOLORED.list"

  CWB::frame fr(FRLIST_FILE);

  cout << fr.getNfiles() << endl;

  std::vector<frfile> frlist = fr.getFrList();  

  cout << frlist.size() << endl;

  for(int i=0;i<frlist.size();i++) {

    cout << i << " " << frlist[i].start << " " << frlist[i].stop << " " << frlist[i].length << endl; 
    cout << frlist[i].file[0].Data() << endl << endl;
  }
  exit(0);
}
