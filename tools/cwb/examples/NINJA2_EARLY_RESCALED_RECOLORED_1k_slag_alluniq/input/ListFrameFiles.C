//
// List frame : start, stop, length, path 
// Author : Gabriele Vedovato


#define FRLIST_FILE "L1_NINJA2_G1000176_EARLY_RECOLORED.list"
//#define FRLIST_FILE "L1_LDAS_C02_L2.list"
//#define FRLIST_FILE "L1_LDAS_C02_L2_big.list"

{

  #include <vector>

  CWB::frame fr(FRLIST_FILE);

  cout << fr.getNfiles() << endl;

  std::vector<frfile> frlist = fr.getFrList();  

  cout << frlist.size() << endl;

  for(int i=1;i<frlist.size();i++) {

  if ( frlist[i].start- frlist[i-1].stop > 0.) {cout << frlist[i].start << " " << frlist[i].stop << " " << frlist[i].length << " " << frlist[i].start- frlist[i-1].stop << endl; }
//    cout << frlist[i].file[0].Data() << endl << endl;
  }
  exit(0);
}
