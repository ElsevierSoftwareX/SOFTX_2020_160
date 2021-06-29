// ===============================================================================
// this macro makes ROOT5 user_parameters C++ compliant (ROOT6)
// ===============================================================================

void FixUserParametersROOT6(TString ifName, TString type) {

#ifdef _USE_ROOT6

  // working dir
  char work_dir[1024];
  sprintf(work_dir,"%s",gSystem->WorkingDirectory());

  // data label
  char data_label[1024];
  sprintf(data_label,"%s",gSystem->BaseName(work_dir));
  //cout << "--------------------------> data_label : " << data_label << endl;

  Long_t id,size,flags,mt;
  int estat = gSystem->GetPathInfo(ifName.Data(),&id,&size,&flags,&mt);
  if (estat==0) {

    ifstream in;
    in.open(ifName.Data(),ios::in);
    if (!in.good()) {cout << "Error Opening File : " << ifName.Data() << endl;exit(1);}

    ofstream out;
    TString ifName_tmp = ifName+".tmp";
    out.open(ifName_tmp.Data(),ios::out);
    if (!out.good()) {cout << "Error Opening File : " << ifName_tmp << endl;exit(1);}

    TString nDQF="";

    char str[1024];
    while(true) {
      in.getline(str,1024);
      if (!in.good()) break;
      TString ostr = str;
      ostr.ReplaceAll(" ","");
      ostr.ReplaceAll("\t","");
      if(type=="p") {				// config/user_parameters.C
        if(ostr.BeginsWith("search=")) {
          ostr = str;
          ostr="// "+ostr+"\t// #ROOT5";
          out << ostr.Data() << endl;
          ostr = str;
          ostr.ReplaceAll("search","cfg_search");
          ostr+="\t// #ROOT6";
        } else if(ostr.BeginsWith("gamma=")) {
          ostr = str;
          ostr="// "+ostr+"\t// #ROOT5";
          out << ostr.Data() << endl;
          ostr = str;
          ostr.ReplaceAll("gamma","cfg_gamma");
          ostr+="\t// #ROOT6";
        } else if(ostr.BeginsWith("nDQF=")) {
          ostr = str;
          ostr = ostr(ostr.First("=")+1,ostr.First(";")-ostr.First("=")-1);
          nDQF=ostr;
          ostr = str;
        } else if(!ostr.BeginsWith("//")&&ostr.Contains("nDQF")&&(nDQF!="")) {
          ostr = str;
          ostr="// "+ostr+"\t// #ROOT5";
          out << ostr.Data() << endl;
          ostr = str;
          ostr.ReplaceAll("nDQF",nDQF);
          ostr+="\t// #ROOT6";
        } else ostr = str;
      }
      if(type=="mp") {				// merge/cwb_user_parameters.C
        if(ostr.BeginsWith("chardata_label[1024]=")) {
          ostr = str;
          ostr="// "+ostr+"\t// #ROOT5";
          out << ostr.Data() << endl;
          ostr = str;
          TString old_label = ostr(ostr.First("=")+1,ostr.First(";")-ostr.First("=")-1);
          ostr.ReplaceAll(old_label,TString(" \"")+TString(data_label)+TString("\""));
          ostr+="\t// #ROOT6";
        } else if (ostr.BeginsWith("char*lagFile=\"")) {
          // fix wrong user_parameter file syntax produced by the config class 
          ostr = str;
          ostr="// "+ostr+"\t// #ROOT5";
          out << ostr.Data() << endl;
          ostr = str;
          ostr = ostr(0, ostr.First("=")+1)+TString(" NULL;");
          ostr+="\t// #ROOT6";
        } else if (ostr.BeginsWith("char*slagFile=\"")) {
          // fix wrong user_parameter file syntax produced by the config class 
          ostr = str;
          ostr="// "+ostr+"\t// #ROOT5";
          out << ostr.Data() << endl;
          ostr = str;
          ostr = ostr(0, ostr.First("=")+1)+TString(" NULL;");
          ostr+="\t// #ROOT6";
        } else ostr = str;
      }
      if(type=="pp") {				// config/user_pparameters.C
        if(ostr.BeginsWith("intchunkID=")) {
          ostr = str;
          ostr="// "+ostr+"\t// #ROOT5";
          out << ostr.Data() << endl;
          ostr = str;
          ostr.ReplaceAll("int ","");
          ostr+="\t// #ROOT6";
        } else if(ostr.BeginsWith("TStringcalibVer=")) {
          ostr = str;
          ostr="// "+ostr+"\t// #ROOT5";
          out << ostr.Data() << endl;
          ostr = str;
          ostr.ReplaceAll("TString ","");
          ostr+="\t// #ROOT6";
        } else ostr = str;
      }
      out << ostr.Data() << endl;
    }
    out.close();
    in.close();

    char cmd[1024];
    sprintf(cmd,"mv %s %s",ifName_tmp.Data(),ifName.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);
  }
#endif

  exit(0);
}

