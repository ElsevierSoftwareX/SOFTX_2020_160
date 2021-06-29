#define IDIR "/home/shubhanshu.tiwari/EBBH/catalogue/O1/O1_C01_SIM_eBBH_10_10_snr_20_e_0_500/output/"
#define OFILE "ebbh_10_10_snr_20_e_0_500.txt"

//#define IDIR "/home/shubhanshu.tiwari/EBBH/catalogue/O1/O1_C01_SIM_eBBH_10_10_snr_20_e_6_500/output/"
//#define OFILE "ebbh_10_10_snr_20_e_6_500.txt"

//#define IDIR "/home/shubhanshu.tiwari/EBBH/catalogue/O1/O1_C01_SIM_eBBH_10_10_snr_20_e_8_500/output/"
//#define OFILE "ebbh_10_10_snr_20_e_8_500.txt"

{  
/*
  TTree itree("ebbh","ebbh")
  itree.ReadFile("ebbh_10_10_snr_20_e_0_500.txt","m1/F:m2:s1/F:s2/F:e/F:rp/F:dist/F:redshift/F:mch/F:ech/F:ei/F")
  itree.Draw("mch")
  itree.Draw("ei")
  itree.Draw("ech")
  itree.Draw("e")
*/
//  TString ifile = TString::Format("%s/*_1_job4000.root",IDIR);
//  ChirpMass(ifile);

  gROOT->LoadMacro("/home/vedovato/WP/eBBH/MACROS/ChirpMass.C");

  EBBH* ebbh = new EBBH;

  TString TAG = ""; 

  vector<TString> fList = CWB::Toolbox::getFileListFromDir(IDIR, ".root", "wave_", TAG, true);

  cout << "Number of files : " << fList.size() << endl;

  ofstream out;
  out.open(OFILE);
  for(int n=0; n<fList.size(); n++) {
    cout << fList[n] << endl; 
    double status = ChirpMass(fList[n],0,0,0,ebbh);

    cout << n << " -> Eccentricity Index = " << ebbh->ei << " status : " << status << endl;

    out << ebbh->mass[0]    << " " << ebbh->mass[0]      << " " 
        << ebbh->spin[2]    << " " << ebbh->spin[5]      << " " 
        << ebbh->e0         << " " << ebbh->rp0          << " " << ebbh->dist << " "   << ebbh->redshift << " "
        << ebbh->ch_mass[0] << " " << ebbh->ch_energy[0] << " " << ebbh->ei   << endl;
  }
  out.close();

  exit(0);
}
