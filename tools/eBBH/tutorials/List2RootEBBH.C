// Convert eBBH parms file list to a root file which contains the eBBH waveforms.  

void 
List2RootEBBH(TString ifName) {

  if(!ifName.EndsWith(".lst")) {
     cout << "List2RootEBBH - bad file extension  : " << ifName << endl;
     cout << "Extention must be .lst" << endl;
     gSystem->Exit(1);
  }

  // the distance of source is G*M/c^2  meters
  double G  = watconstants::GravitationalConstant();
  double M  = watconstants::SolarMass();
  double c  = watconstants::SpeedOfLightInVacuo();
  double pc = watconstants::Parsec();
  double distance_source_Kpc = G*M/(c*c)/pc/1.e3;

  TString ofName = ifName;
  ofName.ReplaceAll(".lst",".root");

  TFile efile(ofName,"RECREATE");

  int id;             
  double m1,m2,rp0,e0;
  wavearray<double>* hp = new wavearray<double>;
  wavearray<double>* hx = new wavearray<double>;

  TTree etree("ebbh","ebbh");
  etree.Branch("id",&id,"id/I");
  etree.Branch("m1",&m1,"m1/D");
  etree.Branch("m2",&m2,"m2/D");
  etree.Branch("rp0",&rp0,"rp0/D");
  etree.Branch("e0",&e0,"e0/D");
  etree.Branch("hp","wavearray<double>",&hp,32000,0);
  etree.Branch("hx","wavearray<double>",&hx,32000,0);

  ifstream in;
  in.open(ifName,ios::in);
  if (!in.good()) {cout << "List2RootEBBH - Error Opening File : " << ifName << endl;exit(1);}

  // get number of entries
  int entries=0;          
  char str[1024];         
  while(true) {           
    in.getline(str,1024); 
    if (!in.good()) break;
    if(str[0] != '#') entries++;
  }                             
  cout << "entries " << entries << endl;
  in.clear(ios::goodbit);               
  in.seekg(0, ios::beg);                

  int fpos=0;         
  while (1) {         
    fpos=in.tellg();  
    in.getline(str,1024);
    if(str[0] == '#') continue;
    if (!in.good()) break;     

    std::stringstream linestream(str);
    if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0)) {
       cout << "List2RootEBBH - Wrong Format for File : " << ifName << endl;
       cout << "input line : " << endl;                                            
       cout << str << endl;                                                        
       cout << "must be : " << endl;                                               
       cout << "event# " << " m1 " << " m2 " << " rp0 " << " e0 " << endl;         
       exit(1);                                                                    
    }                                                                              

    cout << "Create eBBH with parms : " << id << " " << m1 << " " << m2 << " " << rp0 << " " << e0 << " " << endl;
    getEBBH(m1,m2,rp0,e0,*hp,*hx);

    // rescale hp,hx to 10 Kpc
    (*hp)*=distance_source_Kpc/10.;
    (*hx)*=distance_source_Kpc/10.;

    etree.Fill();
  }                          

  in.close();

  etree.Write();
  efile.Close();

  gSystem->Exit(0);
}
