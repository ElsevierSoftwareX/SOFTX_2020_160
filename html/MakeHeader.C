void MakeHeader() {

  ifstream in;
  in.open("etc/html/header_template.html",ios::in);
  if (!in.good()) {cout << "Error Opening File : " << "etc/html/header_template.html" << endl;exit(1);}

  // get SITE_CLUSTER
  TString site_cluster="VIRTUALBOX";                                            // default value
  if(gSystem->Getenv("SITE_CLUSTER")!=NULL) {
    site_cluster=TString(gSystem->Getenv("SITE_CLUSTER"));
  }
  TString cluster_site_logo  = "cluster_site_logo_modern.png";
  TString cluster_site_url1  = "http://www.ligo.caltech.edu/";
  TString cluster_site_name1 = "LIGO Homepage";
  TString cluster_site_url2  = "http://www.virgo-gw.eu/";
  TString cluster_site_name2 = "VIRGO Homepage";
  if(site_cluster=="ATLAS") {
    cluster_site_logo  = "atlas_logo_modern.png";
    cluster_site_url1  = "http://www.aei.mpg.de/14026/AEI_Hannover/";
    cluster_site_name1 = "AEI Hannover Homepage";
    cluster_site_url2  = "http://www.aei.mpg.de/14026/AEI_Hannover";
    cluster_site_name2 = "AEI Hannover Homepage";
  }
  if(site_cluster=="CIT") {
    cluster_site_logo  = "ligo_virgo_logo_modern.png";
    cluster_site_url1  = "http://www.ligo.caltech.edu/";
    cluster_site_name1 = "LIGO Homepage";
    cluster_site_url2  = "https://www.virgo-gw.eu/";
    cluster_site_name2 = "VIRGO Homepage";
  }
  if(site_cluster=="VIRTUALBOX") {
    cluster_site_logo  = "cluster_site_logo_modern.png";
    cluster_site_url1  = "https://www.gw-openscience.org/about/";
    cluster_site_name1 = "GWOSC Homepage";
    cluster_site_url2  = "https://www.gw-openscience.org/about/";
    cluster_site_name2 = "GWOSC Homepage";
  }

  // get HOME_WWW
  TString home_www="~waveburst/waveburst";  					// default value
  if(gSystem->Getenv("HOME_WWW")!=NULL) {
    home_www=TString(gSystem->Getenv("HOME_WWW"));
  }

  // get HOME_CED_WWW
  TString home_ced_www="~waveburst/waveburst/ced-1.0-modern";    		// default value
  if(gSystem->Getenv("HOME_CED_WWW")!=NULL) {
    home_ced_www=TString(gSystem->Getenv("HOME_CED_WWW"));
  }

  // get CWB_DOC_URL
  TString cwb_doc_url="https://ldas-jobs.ligo.caltech.edu/~waveburst/doc";    	// default value
  if(gSystem->Getenv("CWB_DOC_URL")!=NULL) {
    cwb_doc_url=TString(gSystem->Getenv("CWB_DOC_URL"));
  }

  // get CWB_GIT_URL
  TString cwb_git_url="https://git.ligo.org";    				// default value
  if(gSystem->Getenv("CWB_GIT_URL")!=NULL) {
    cwb_git_url=TString(gSystem->Getenv("CWB_GIT_URL"));
  }

  // get issues URL (Ex: https://gitlab.com/groups/gwburst/-/issues)
  TString cwb_git_issues_url=cwb_git_url;
  TString partA=cwb_git_issues_url;
  TString partB=cwb_git_issues_url;
  partA.Remove(partA.Last('/'),partA.Sizeof()-1);
  partB.Remove(0,partB.Last('/')+1);
  cwb_git_issues_url=partA+"/groups/"+partB+"/-/issues";

  ofstream out;
  out.open("etc/html/header.html",ios::out);

  char str[2048];
  while(true) {
    in.getline(str,2048);
    if (!in.good()) break;
    //cout << str << endl;
    TString line = str;
    // replace HOME_CED_WWW, HOME_CED_WWW, CWB_DOC_URL with used defined values (watenv)
    line.ReplaceAll("HOME_WWW",home_www);
    line.ReplaceAll("HOME_CED_WWW",home_ced_www);
    if(cwb_doc_url.Contains("gwburst.gitlab.io")) {	// public cWB documentation	
      line.ReplaceAll("/CWB_DOC_URL/cwb/man",cwb_doc_url);
    } else {
      line.ReplaceAll("CWB_DOC_URL",cwb_doc_url);
    }
    line.ReplaceAll("CWB_GIT_URL",cwb_git_url);
    line.ReplaceAll("CWB_GIT_ISSUES_URL",cwb_git_issues_url);
    line.ReplaceAll("CLUSTER_SITE_LOGO",cluster_site_logo);
    line.ReplaceAll("CLUSTER_SITE_URL1",cluster_site_url1);
    line.ReplaceAll("CLUSTER_SITE_NAME1",cluster_site_name1);
    line.ReplaceAll("CLUSTER_SITE_URL2",cluster_site_url2);
    line.ReplaceAll("CLUSTER_SITE_NAME2",cluster_site_name2);

    line.ReplaceAll(TString("=\"/http"), TString("=\"http"));
    line.ReplaceAll(TString("=\"//"), TString("=\"/"));

    out << line.Data() << endl;
  }

  out.close();

  exit(0);
}

