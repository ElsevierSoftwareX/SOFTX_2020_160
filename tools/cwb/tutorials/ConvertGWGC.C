// 
// convert GWGCCatalog_Rev1d7.txt format to GWGC format for inspiral input sky distribution
// Author : Gabriele Vedovato

{
  #define GWGC_FILE "../data/GWGCCatalog_Rev1d7.txt"
  #define GWGC_FILE_INSP "../data/GWGCCatalog_Rev1d7_InspiralFormat.txt"

  ifstream in;
  in.open(GWGC_FILE,ios::in);
  if (!in.good()) {cout << "Error Opening Input File : " << GWGC_FILE << endl;exit(1);}
  cout << "Opening Input File : " << GWGC_FILE << endl;

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

  ofstream out;
  out.open(GWGC_FILE_INSP,ios::out);
  if (!out.good()) {cout << "Error Opening Output File : " << GWGC_FILE_INSP << endl;exit(1);}

int cnt=0;
  char iline[1024];
  in.getline(iline,1024);  // skip first line (header)
  while (1) {

    in.getline(iline,1024);
    if (!in.good()) break;
    TObjArray* tok = TString(iline).Tokenize(TString('|'));

    TObjString* tname    = (TObjString*)tok->At(1);
    TObjString* tra      = (TObjString*)tok->At(2);
    TObjString* tdec     = (TObjString*)tok->At(3);
    TObjString* tdist    = (TObjString*)tok->At(14);
    TObjString* tabs_mag = (TObjString*)tok->At(13);

    TString name   = tname->GetString();
    double ra      = tra->GetString().Atof();
    double DEC     = tdec->GetString().Atof();
    double dist    = tdist->GetString().Atof();     // Mpc
    double lum     = tabs_mag->GetString().Atof();

    DEC = DEC>0 ? 90-DEC : -DEC-90;
    ra-=12;
//    dist*=1000;  // Kpc

    int    dec_d   = int(DEC); 
    int    dec_m   = int((DEC-dec_d)*60.);  
    if(dec_d<0) dec_d=-dec_d; 
    if(dec_m<0) dec_m=-dec_m; 
    char   dec_sgn = DEC>=0 ? '+' : '-';

    int    ra_h    = int(ra);  
    int    ra_m    = int(60.*(ra-ra_h));  
    if(ra_h<0) ra_h=-ra_h; 
    if(ra_m<0) ra_m=-ra_m; 
    char   ra_sgn  = ra>=0 ? '+' : '-';

    out << name.Data() << "\t" << ra_sgn <<ra_h <<":"<<ra_m << "\t"
                               << dec_sgn<<dec_d<<":"<<dec_m << "\t"
                               << dist << "\t" << lum << "\t" << 1 << endl;
/*
if((cnt++)%1000==0) {
cout << "ra " << ra << " dec " << DEC << " dist " << endl;
    cout << name.Data() << "\t" << ra_sgn <<ra_h <<":"<<ra_m << "\t"
                              << dec_sgn<<dec_d<<":"<<dec_m << "\t"
                              << dist << "\t" << lum << "\t" << 1 << endl;
//exit(0);
}
*/
  }

  out.close();

  exit(0);
}
