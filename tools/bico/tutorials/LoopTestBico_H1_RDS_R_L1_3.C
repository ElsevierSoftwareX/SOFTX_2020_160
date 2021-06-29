//
// Loop bicoherence over channels
// Author : Gabriele Vedovato

{
  #define CHANNEL_LIST "input/H1_PEM_Channels.txt"


  gROOT->LoadMacro("TestBico_H1_RDS_R_L1_3.C");

  ifstream in;
  in.open(CHANNEL_LIST);
  if(!in.good()) {cout << "Error Opening File : " << CHANNEL_LIST << endl;exit(1);}

  char chname[1024];
  while(1) {
   in >> chname;
   if (!in.good()) break;
   cout << "---------------------------------------------------------------" << endl;
   cout << chname << endl;
   cout << "---------------------------------------------------------------" << endl;
   TestBico_H1_RDS_R_L1_3(chname);
  }

  in.close();
}
