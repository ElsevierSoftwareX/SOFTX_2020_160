// this macro produces an html page : 1G vs 2G efficiency 50% comparison
// Author : Gabriele Vedovato
//
// run ex : root 'Efficiency2Gvs1G.C("definitions.C",50)'
// definitions.C contains all definitions
//
// --------------------------------------------------------------
// 1G vs 2G pipeline definitions
// --------------------------------------------------------------
// #define TITLE_LABEL ""	// report title	
// #define OFILE_LABEL ""	// report/dump/s/Efficiency2Gvs1G_OFILE_LABEL
// --------------------------------------------------------------
// 1G pipeline definitions
// --------------------------------------------------------------
// #define PDIR_1G ""		// working dir
// #define SDIR_1G ""		// report sub directory (report/postprod/XXX/data)
// #define SWWW_1G ""		// simulation www link
// #define BWWW_1G ""		// background www link
// --------------------------------------------------------------
// 2G pipeline definitions
// --------------------------------------------------------------
// #define PDIR_2G ""		// working dir
// #define SDIR_2G ""		// report sub directory (report/postprod/XXX/data)
// #define SWWW_2G ""		// simulation www link
// #define BWWW_2G ""		// background www link
//

int GetEfficiency(double eff, TString directory, vector<TString>& vname, vector<double>& vhrss);

void Efficiency2Gvs1G(TString cfgFile, int efficiency=50) {

  if(efficiency<0 || efficiency>100) {
    cout << "Efficiency2Gvs1G : Error - efficiency must be [0,100]" << endl;
    exit(1);
  }

  gROOT->Macro(cfgFile);

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));

  gROOT->Macro(gSystem->ExpandPathName("$CWB_PARAMETERS_FILE"));
  gROOT->Macro(gSystem->ExpandPathName("$CWB_UPARAMETERS_FILE"));
//  gROOT->Macro(gSystem->ExpandPathName("$CWB_PPARAMETERS_FILE"));

  vector<TString> vname[2];
  vector<double>  vhrss[2];

  GetEfficiency(efficiency/100., TString(PDIR_1G)+"/"+SDIR_1G, vname[0], vhrss[0]); 
  GetEfficiency(efficiency/100., TString(PDIR_2G)+"/"+SDIR_2G, vname[1], vhrss[1]); 

  cout << "1G efficiencies" << endl;
  for(int i=0;i<vname[0].size();i++) {
   cout << vname[0][i] << " " << vhrss[0][i] << endl;
  } 
  cout << endl;

  cout << "2G efficiencies" << endl;
  for(int i=0;i<vname[1].size();i++) {
    cout << vname[1][i] << " " << vhrss[1][i] << endl;
  } 
  cout << endl;

  char title[1024];
  sprintf(title,"%s - Efficiency %d",TITLE_LABEL,efficiency);

  cout << "1G vs 2G efficiencies" << endl;
  char eff_file[1024];sprintf(eff_file,"Efficiency2Gvs1G_%s_eff%d.html",OFILE_LABEL,efficiency);
  ofstream out;
  out.open(eff_file,ios::out);              // create Efficiency2Gvs1G.html
  char oline[1024];
  out << "<html>" << endl;
  out << "<font color=\"red\" style=\"font-weight:bold;\"><center><p><h2>"
      << title << "</h2><p><center></font>" << endl;

  out << "<h4>1G : " << endl;
  sprintf(oline,"<a target=\"_parent\" href=\"%s\">bkg report page</a>",BWWW_1G);
  out << oline << endl;
  out << "---" << endl;
  sprintf(oline,"<a target=\"_parent\" href=\"%s\">sim report page</a>",SWWW_1G);
  out << oline << endl;
  out << "</h4>" << endl;
  out << "<h4>2G : " << endl;
  sprintf(oline,"<a target=\"_parent\" href=\"%s\">bkg report page</a>",BWWW_2G);
  out << oline << endl;
  out << "---" << endl;
  sprintf(oline,"<a target=\"_parent\" href=\"%s\">sim report page</a>",SWWW_2G);
  out << oline << endl;
  out << "</h4>" << endl;

  out << "<table border=1 cellpadding=12 align=\"center\">" << endl;
  out << "<br>" << endl;

  out << "<tr align=\"center\">" << endl;
  out << "<td align=\"center\"><font color=\"red\"> ID </font></td>" << endl;
  out << "<td align=\"center\"><font color=\"red\"> MDC </font></td>" << endl;
  out << "<td align=\"center\"><font color=\"red\"> 1G efficiency 50% </font></td>" << endl;
  out << "<td align=\"center\"><font color=\"red\"> 2G efficiency 50% </font></td>" << endl;
  out << "<td align=\"center\"><font color=\"red\"> 100*(2G-1G)/Min(1G,2G) </font></td>" << endl;
  out << "</tr>" << endl;

  for(int i=0;i<vname[1].size();i++) {
    for(int j=0;j<vname[0].size();j++) {
      if(vname[1][i]==vname[0][j]) {

        out << "<tr align=\"center\">" << endl;
        sprintf(oline,"<td align=\"center\"><font color=\"black\">%d</font></td>",i+1);
        out << oline << endl;
        sprintf(oline,"<td align=\"center\"><font color=\"blue\">%s</font></td>",vname[1][i].Data());
        out << oline << endl;
        sprintf(oline,"<td align=\"center\"><font color=\"black\">%2.2g</font></td>",vhrss[0][j]);
        out << oline << endl;
        sprintf(oline,"<td align=\"center\"><font color=\"black\">%2.2g</font></td>",vhrss[1][i]);
        out << oline << endl;
        double percentage = 100.*(vhrss[1][i]-vhrss[0][j])/TMath::Min(vhrss[0][j],vhrss[1][i]);
        if(percentage<=0) {
          sprintf(oline,"<td align=\"center\"><font color=\"blue\">%2.1f</font></td>",percentage);
        } else {
          sprintf(oline,"<td align=\"center\"><font color=\"red\">%2.1f</font></td>",percentage);
        }
        out << oline << endl;
        out << "</tr>" << endl;

        cout << i << " " << vname[1][i] << " " << vhrss[1][i] << " "  << vhrss[0][j] << endl;
      }
    } 
  } 

  out << "</table>" << endl;
  out << "</html>" << endl;
  out.close();

  cout << endl << eff_file << endl << endl;

  exit(0);
}

int
GetEfficiency(double eff, TString directory, vector<TString>& vname, vector<double>& hrss) {

  vector<TString> fileList = CWB::Toolbox::getFileListFromDir(directory, ".txt", "","fit_parameters");

  for(int i=0;i<fileList.size();i++) {

    cout << i << " " << fileList[i].Data() << endl;
    // 16 2.014E+01 2.695E-21 +- 7.012E-23 3.454E-01 5.000E-01 7.552E-01 SG1615Q100

    ifstream in(fileList[i]);
    if(!in.good()) {cout << "Error Opening File : " << fileList[i] << endl;exit(1);}

    char dummy[1024];
    char name[1024];
    double chi2, hrss50, sigma, betam, betap, hrssEr;
    while(1) {
      in >> dummy >> chi2 >> hrss50 >> dummy >> dummy >> sigma >> betam >> betap >> name;
      if (!in.good()) break;
      vname.push_back(name);

      double par0=TMath::Log10(hrss50);
      TF1 fitFunc("logNfit",logNfit,pow(10.0,-23.0),pow(10.0,-18.5),5);
      fitFunc.SetNpx(10000);
      fitFunc.SetParameters(par0,sigma,betam,betap,0);
      double fhrss=fitFunc.GetX(eff,pow(10.0,-25.0),pow(10.0,-18.5));
      hrss.push_back(fhrss);
    }
    in.close();
  }

  return 0;
}
