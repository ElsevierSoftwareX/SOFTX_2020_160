void MakeFooter() {

  ifstream in;
  in.open("etc/html/footer_template.html",ios::in);
  if (!in.good()) {cout << "Error Opening File : " << "etc/html/footer_template.html" << endl;exit(1);}

  ofstream out;
  out.open("etc/html/footer.html",ios::out);

  // get current date
  TDatime date;
  date.Set();

  // get current git version
  gSystem->Exec("git rev-parse --abbrev-ref HEAD | tr '\n' ' ' > gitversion.txt");
  gSystem->Exec("git -C ${HOME_WAT} tag -l --points-at HEAD | tr '\n' ' ' >> gitversion.txt");
  gSystem->Exec("git -C ${HOME_WAT} rev-parse HEAD | tr '\n' ' ' >> gitversion.txt");
  ifstream ingit;
  ingit.open("gitversion.txt",ios::in);
  if (!ingit.good()) {cout << "Error Opening File : " << "gitversion.txt" << endl;exit(1);}
  char git_version[2048];
  ingit.getline(git_version,2048);
  gSystem->Exec("rm gitversion.txt");
  ingit.close();

  char str[2048];
  while(true) {
    in.getline(str,2048);
    if (!in.good()) break;
    //cout << str << endl;
    TString line = str;
    // replace GENERATED with generation date
    line.ReplaceAll("DATE_GENERATED",date.AsString());
    line.ReplaceAll("GIT_VERSION",git_version);
    out << line.Data() << endl;
  }

  out.close();

  exit(0);
}

