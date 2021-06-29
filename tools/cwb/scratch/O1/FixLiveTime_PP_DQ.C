// Author : GV
// this procedure read background files (rho vs rate) computed with livetime before DQ
// and write new background files (rho vs rate) computed with livetime after DQ
//

#define LIVE_TIME_BEFORE_PP_DQ	67953.8		// years
#define LIVE_TIME_AFTER_PP_DQ	67406.8		// years

void 
FixLiveTime_PP_DQ(TString ifName, bool icol2=false, bool ocol2=false, bool force_exit=true, TString tag="") {

  if(!ifName.EndsWith(".txt")) {cout << "Error, file extention not '.txt'" << endl;exit(1);}

  TString ofName = ifName;
  if(tag=="") ofName.ReplaceAll(".txt","_LTDQ.txt");
  else        ofName.ReplaceAll(".txt",TString("_")+tag+"_LTDQ.txt");

  cout << "Input  File : " << ifName << endl;
  cout << "Output File : " << ofName << endl;

  double live_time_ratio = LIVE_TIME_BEFORE_PP_DQ/LIVE_TIME_AFTER_PP_DQ;

  double rho,erho;     
  double far,efar;     

  ifstream in;
  in.open(ifName.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << ifName.Data() << endl;exit(1);}

  ofstream out;
  out.open(ofName,ios::out);
  if(!out.good()) {cout << "Error : Opening File : " << ofName << endl;gSystem->Exit(1);}

  // fixed live time
  while (1) {
    if(icol2) in >> rho >> far;
    else      in >> rho >> far >> erho >> efar;
    if (!in.good()) break;
    double far_LTDQ = far*live_time_ratio;
    double efar_LTDQ = efar*live_time_ratio;
    //cout << rho << " " << far << " " << far_LTDQ << endl;
    if(ocol2) out << rho << "\t" << far_LTDQ << endl;
    else      out << rho << "\t" << far_LTDQ << "\t" << erho << "\t" << efar_LTDQ << endl;
  }

  in.close();
  out.close();

  if(force_exit) exit(0);
}
