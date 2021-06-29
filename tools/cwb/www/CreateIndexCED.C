#define CEDINDEX_HEADER 	"cedindex_header.html"
#define CEDINDEX_BODY 		"cedindex_body.html"
#define CEDINDEX_TFMAP    	"cedindex_tfmap.html"
#define CEDINDEX_LIKELIHOOD 	"cedindex_likelihood.html"
#define CEDINDEX_POLARGRAM 	"cedindex_polargram.html"
#define CEDINDEX_REC_SIGNAL 	"cedindex_rec_signal.html"
#define CEDINDEX_INJ_SIGNAL 	"cedindex_inj_signal.html"
#define CEDINDEX_SKYMAP 	"cedindex_skymap.html"

#define CEDINDEX_DIR		"index"

void CreateIndexCED(int nifo, TString* ifo, bool sim=false) {

  // get SITE_CLUSTER
  TString site_cluster="VIRTUALBOX";                                            // default value
  if(gSystem->Getenv("SITE_CLUSTER")!=NULL) {
    site_cluster=TString(gSystem->Getenv("SITE_CLUSTER"));
  }
  TString cluster_site_logo  = "cluster_site_logo_modern.png";
  TString cluster_site_url1  = "http://www.ligo.caltech.edu/";
  TString cluster_site_name1 = "LIGO Homepage";
  TString cluster_site_url2  = "https://www.virgo-gw.eu/";
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
  TString home_www="~waveburst/waveburst";                                      // default value
  if(gSystem->Getenv("HOME_WWW")!=NULL) {
    home_www=TString(gSystem->Getenv("HOME_WWW"));
  }

  // get HOME_CED_WWW
  TString home_ced_www="~waveburst/waveburst/ced-1.0-modern";                   // default value
  if(gSystem->Getenv("HOME_CED_WWW")!=NULL) {
    home_ced_www=TString(gSystem->Getenv("HOME_CED_WWW"));
  }

  // get CWB_DOC_URL
  TString cwb_doc_url="https://ldas-jobs.ligo.caltech.edu/~waveburst/doc";      // default value
  if(gSystem->Getenv("CWB_DOC_URL")!=NULL) {
    cwb_doc_url=TString(gSystem->Getenv("CWB_DOC_URL"));
  }

  // get CWB_GIT_URL
  TString cwb_git_url="https://git.ligo.org";                                   // default value
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

  // It creates CEDINDEX_DIR dir if not already there
  if (!gSystem->OpenDirectory(CEDINDEX_DIR)) gSystem->MakeDirectory(CEDINDEX_DIR);

  ifstream in;
  in.open(CEDINDEX_BODY,ios::in);
  if(!in.good()) {cout << "Error Opening File : " << CEDINDEX_BODY << endl;exit(1);}

  char ofileName[256]="";
  if(sim) sprintf(ofileName,"%s/cedindex_sim_",CEDINDEX_DIR); 
  else    sprintf(ofileName,"%s/cedindex_",CEDINDEX_DIR); 
  for(int n=0;n<nifo;n++) sprintf(ofileName,"%s%s",ofileName,ifo[n].Data());
  sprintf(ofileName,"%s.html",ofileName);
  cout << "ofileName : " <<  ofileName << endl;

  ofstream out;
  out.open(ofileName,ios::out);
  if(!out.good()) {cout << "Error Opening File : " << ofileName << endl;exit(1);}

  char line[2048];

  ifstream inx;
  inx.open(CEDINDEX_HEADER,ios::in);
  if(!inx.good()) {cout << "Error Opening File : " << CEDINDEX_HEADER << endl;exit(1);}
  while(1) {
    inx.getline(line,2048);
    if(!inx.good()) break;
    TString sline = line;
    // replace HOME_CED_WWW, HOME_CED_WWW, CWB_DOC_URL with used defined values (watenv)
    sline.ReplaceAll("HOME_WWW",home_www);
    sline.ReplaceAll("HOME_CED_WWW",home_ced_www);
    if(cwb_doc_url.Contains("gwburst.gitlab.io")) {     // public cWB documentation     
      sline.ReplaceAll("/CWB_DOC_URL/cwb/man",cwb_doc_url);
    } else {
      sline.ReplaceAll("CWB_DOC_URL",cwb_doc_url);
    }
    sline.ReplaceAll("CWB_GIT_URL",cwb_git_url);
    sline.ReplaceAll("CWB_GIT_ISSUES_URL",cwb_git_issues_url);
    sline.ReplaceAll("CLUSTER_SITE_LOGO",cluster_site_logo);
    sline.ReplaceAll("CLUSTER_SITE_URL1",cluster_site_url1);
    sline.ReplaceAll("CLUSTER_SITE_NAME1",cluster_site_name1);
    sline.ReplaceAll("CLUSTER_SITE_URL2",cluster_site_url2);
    sline.ReplaceAll("CLUSTER_SITE_NAME2",cluster_site_name2);

    sline.ReplaceAll(TString("=\"/http"), TString("=\"http"));
    sline.ReplaceAll(TString("=\"//"), TString("=\"/"));

    out << sline.Data() << endl;
  }
  inx.close();

  while(1) {
    in.getline(line,2048);
    if(!in.good()) break;
    //cout << line << endl;
    TString sline = line;

    if(cwb_doc_url!="") {
      sline.ReplaceAll("<!--CWB_DOC_URL","");
      sline.ReplaceAll("CWB_DOC_URL-->","");
      sline.ReplaceAll("XCWB_DOC_URL",cwb_doc_url.Data());
    }

    sline.ReplaceAll("HOME_CED_WWW",home_ced_www.Data());
    sline.ReplaceAll(TString("=\"/http"), TString("=\"http"));
    sline.ReplaceAll(TString("=\"//"), TString("=\"/"));
    out << sline.Data() << endl;

    if(sline.Contains("<!--TFMAP-->")) {
      for(int n=0;n<nifo;n++) {
        ifstream inx;
        inx.open(CEDINDEX_TFMAP,ios::in);
        if(!inx.good()) {cout << "Error Opening File : " << CEDINDEX_TFMAP << endl;exit(1);}
        while(1) {
          inx.getline(line,1024);
          if(!inx.good()) break;
          TString sline2 = line;
          sline2.ReplaceAll("IFO",ifo[n]);
          out << sline2.Data() << endl;
        }
        inx.close();
      }
    }

    if(sline.Contains("<!--LIKELIHOOD-->")) {
      ifstream inx;
      inx.open(CEDINDEX_LIKELIHOOD,ios::in);
      if(!inx.good()) {cout << "Error Opening File : " << CEDINDEX_LIKELIHOOD << endl;exit(1);}
      while(1) {
        inx.getline(line,1024);
        if(!inx.good()) break;
        out << line << endl;
      }
      inx.close();
    }

    if(sline.Contains("<!--POLARGRAM-->")) {
      ifstream inx;
      inx.open(CEDINDEX_POLARGRAM,ios::in);
      if(!inx.good()) {cout << "Error Opening File : " << CEDINDEX_POLARGRAM << endl;exit(1);}
      while(1) {
        inx.getline(line,1024);
        if(!inx.good()) break;
        out << line << endl;
      }
      inx.close();
    }

    if(sline.Contains("<!--SIGNAL-->")) {
      for(int n=0;n<nifo;n++) {
        ifstream inx;
        inx.open(CEDINDEX_REC_SIGNAL,ios::in);
        if(!inx.good()) {cout << "Error Opening File : " << CEDINDEX_REC_SIGNAL << endl;exit(1);}
        while(1) {
          inx.getline(line,1024);
          if(!inx.good()) break;
          TString sline2 = line;
          sline2.ReplaceAll("IFO",ifo[n]);
          out << sline2.Data() << endl;
        }
        if(!sim) continue;
        ifstream iny;
        iny.open(CEDINDEX_INJ_SIGNAL,ios::in);
        if(!iny.good()) {cout << "Error Opening File : " << CEDINDEX_INJ_SIGNAL << endl;exit(1);}
        while(1) {
          iny.getline(line,1024);
          if(!iny.good()) break;
          TString sline2 = line;
          sline2.ReplaceAll("IFO",ifo[n]);
          out << sline2.Data() << endl;
        }
        iny.close();
      }
    }

    if(sline.Contains("<!--SKYMAP-->")&&(nifo>1)) {
      ifstream inx;
      inx.open(CEDINDEX_SKYMAP,ios::in);
      if(!inx.good()) {cout << "Error Opening File : " << CEDINDEX_SKYMAP << endl;exit(1);}
      while(1) {
        inx.getline(line,1024);
        if(!inx.good()) break;
        out << line << endl;
      }
      inx.close();
    }
  }

  return;
}

