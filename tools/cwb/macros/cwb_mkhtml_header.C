/*
# Copyright (C) 2019 Gabriele Vedovato, Marco Drago
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


// make the html header for simulation/production reports : used by the cwb_report command

{
  #include <unistd.h>
  #include <stdio.h>

  // check if files exist
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  char html_footer[1024];
  sprintf(html_footer,"%s/html/etc/html/footer.html",gSystem->Getenv("HOME_WAT"));
  if(gSystem->AccessPathName(html_footer)) {
    cout << endl;
    cout << "Please do : " << endl;
    cout << "cd $HOME_WAT" << endl;
    cout << "make html_footer" << endl;
    cout << endl;
    exit(1);
  }

  // if CWB_DOC_URL is define then man infos are added to web pages 
  TString cwb_doc_url="";
  if(gSystem->Getenv("CWB_DOC_URL")!=NULL) {
    cwb_doc_url=TString(gSystem->Getenv("CWB_DOC_URL"));
  }

  // get site cluster name
  TString site_cluster="";
  if(gSystem->Getenv("SITE_CLUSTER")==NULL) {
    cout << "Error : environment SITE_CLUSTER is not defined!!!" << endl;exit(1);
  } else {
    site_cluster=TString(gSystem->Getenv("SITE_CLUSTER"));
  }

  // get the number of jobs in the merged list
  int iversion;
  cwb_merge_label;
  if(gSystem->Getenv("CWB_MERGE_LABEL")==NULL) {
    cout << "Error : environment CWB_MERGE_LABEL is not defined!!!" << endl;exit(1);
  } else {
    cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  }
  TObjArray* token = TString(cwb_merge_label).Tokenize(TString("."));
  TString merge_label=((TObjString*)token->At(0))->GetString();
  // check if label has the correct format (M#)
  if(merge_label[0]!='M') {
    cout << "Error : label " << merge_label.Data() << " has bad format (M#)" << endl;exit(1);
  } else {
    TString lcheck=merge_label;
    lcheck.Remove(0,1);
    if(!lcheck.IsDigit()) {
      cout << "Error : label " << merge_label.Data() << " has bad format (M#)" << endl;exit(1);
    } else {
      iversion=lcheck.Atoi();
    }
  }
  char ilstfname[1024];
  sprintf(ilstfname,"%s/merge_%s.%s.lst",merge_dir,data_label,cwb_merge_label.Data());
  vector<TString> merge_jobFileList;
  vector<int> merge_jobList = CWB::Toolbox::getMergeJobList(ilstfname,merge_jobFileList);
  cout << "Number of jobs in the merged lists : " << merge_jobList.size() << endl;

  if(nfactor<=0) nfactor=1;	// fix nfactor if not initialized

  strcpy(work_dir,gSystem->WorkingDirectory());

  // get the number of job submit by condor
  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);

  char condor_dag_file[1024];
  sprintf(condor_dag_file,"%s/%s.dag",full_condor_dir,data_label);
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(condor_dag_file,&id,&size,&flags,&mt);
  vector<int> condor_jobList;
  if (estat==0) {
    condor_jobList=CWB::Toolbox::getCondorJobList(full_condor_dir, data_label);
    cout << "Number of jobs/factors in the condor lists : " << condor_jobList.size() << "/" << nfactor << endl;
  }

  char merge_job_list_str[64];
  if(merge_jobList.size()<nfactor*condor_jobList.size()) {
    sprintf(merge_job_list_str,"<font color=\"red\">%d</font>",(int)merge_jobList.size());
  } else {
    sprintf(merge_job_list_str,"%d",(int)merge_jobList.size());
  }
  char condor_job_list_str[64];
  char percentage_job_list_str[64];
  if(condor_jobList.size()>0 && merge_jobList.size()>0) {
    if(simulation) {
      sprintf(condor_job_list_str,"(%d*%d)", nfactor,(int)condor_jobList.size());
      sprintf(percentage_job_list_str,"<font color=\"blue\">%2.1f%%</font>",
              100.*merge_jobList.size()/(nfactor*condor_jobList.size()));
    } else {
      sprintf(condor_job_list_str,"%d", (int)condor_jobList.size());
      sprintf(percentage_job_list_str,"<font color=\"blue\">%2.1f%%</font>",
              100.*merge_jobList.size()/condor_jobList.size());
    }
  } else {	// super run
    sprintf(condor_job_list_str,"0");
    sprintf(percentage_job_list_str,"<font color=\"blue\">%2.1f%%</font>",0.);
  }

  // save configuration files to report directory
  TString home_wat;
  if(gSystem->Getenv("HOME_WAT")==NULL) {
    cout << "Error : environment HOME_WAT is not defined!!!" << endl;exit(1);
  } else {
    home_wat=TString(gSystem->Getenv("HOME_WAT"));
  }
  TString cwb_netc_file;
  if(gSystem->Getenv("CWB_NETC_FILE")==NULL) {
    cout << "Error : environment CWB_NETC_FILE is not defined!!!" << endl;exit(1);
  } else {
    cwb_netc_file=TString(gSystem->Getenv("CWB_NETC_FILE"));
  }
  TString cwb_parameters_file;
  if(gSystem->Getenv("CWB_PARAMETERS_FILE")==NULL) {
    cout << "Error : environment CWB_PARAMETERS_FILE is not defined!!!" << endl;exit(1);
  } else {
    cwb_parameters_file=TString(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  }
  TString cwb_uparameters_file;
  if(gSystem->Getenv("CWB_UPARAMETERS_FILE")==NULL) {
    cout << "Error : environment CWB_UPARAMETERS_FILE is not defined!!!" << endl;exit(1);
  } else {
    cwb_uparameters_file=TString(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  }
  TString cwb_pparameters_file;
  if(gSystem->Getenv("CWB_PPARAMETERS_FILE")==NULL) {
    cout << "Error : environment CWB_PPARAMETERS_FILE is not defined!!!" << endl;exit(1);
  } else {
    cwb_pparameters_file=TString(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  }
  TString cwb_upparameters_file;
  if(gSystem->Getenv("CWB_UPPARAMETERS_FILE")==NULL) {
    cout << "Error : environment CWB_UPPARAMETERS_FILE is not defined!!!" << endl;exit(1);
  } else {
    cwb_upparameters_file=TString(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  }

  char cmd[1024];
  char odir[1024];
  sprintf(odir,"%s/%s",pp_dir,pp_data_dir);	// output report directory

  // get production lib versions from merged root file  
  TString merged_file_name = simulation ? sim_file_name : net_file_name;
  TFile* wfile = new TFile(merged_file_name);
  if(wfile==NULL) {cout << "Error : file " << merged_file_name.Data() << " not found" <<  endl;exit(1);}
  // read history object & build a string with the multistage git versions
  CWB::History* hst = (CWB::History*)wfile->Get("history");
  TString merge_cuts="";
  TString pr_wat_ver="",pr_wat_git="",pr_cfg_ver="",pr_framelib_ver="",pr_lal_ver="",pr_root_ver="",full_parameters_file="";         
  if(hst) {                                                                                             
    TString log="";                                                                                     
    TList* stageList = hst->GetStageNames();    // get stage list                                       
    TList* typeList  = hst->GetTypeNames();	// get type list
    bool hwat=false,hgit=false,hcfg_branch=false,hcfg_tag=false,hcfg_diff=false,hfrlib=false,hroot=false,hlal=false;
    for(int j=0;j<typeList->GetSize();j++) {	// check if history types are defined
      TObjString* typeObjString = (TObjString*)typeList->At(j);
      TString typeName = typeObjString->GetString();
      if(typeName=="WATVERSION")   hwat=true;
      if(typeName=="GITVERSION")   hgit=true;
      if(typeName=="FRLIBVERSION") hfrlib=true;
      if(typeName=="ROOTVERSION")  hroot=true;
      if(typeName=="LALVERSION")   hlal=true;
      if(typeName=="CWB_CONFIG_BRANCH") hcfg_branch=true;
      if(typeName=="CWB_CONFIG_TAG")    hcfg_tag=true;
      if(typeName=="CWB_CONFIG_DIFF")   hcfg_diff=true;
    }
    delete typeList;
    for(int i=0;i<stageList->GetSize();i++) {   // get lib versions                                     
      TObjString* stageObjString = (TObjString*)stageList->At(i);                                       
      TString stageName = stageObjString->GetString();                                                  
      char* stage = const_cast<char*>(stageName.Data());                                                
      TString stageLabel = stageName; stageLabel.Resize(2); stageLabel+=TString("-");                   
      log = hwat   ? hst->GetHistory(stage,const_cast<char*>("WATVERSION"))   : "";   
      if(log!="") pr_wat_ver += log+" ";              
      log = hgit   ? hst->GetHistory(stage,const_cast<char*>("GITVERSION"))   : "";   
      if(log!="") pr_wat_git += stageLabel+log+" ";   
      log = hfrlib ? hst->GetHistory(stage,const_cast<char*>("FRLIBVERSION")) : ""; 
      if(log!="") pr_framelib_ver = log;                   
      log = hroot  ? hst->GetHistory(stage,const_cast<char*>("ROOTVERSION"))  : "";  
      if(log!="") pr_root_ver = log;                   
      log = hlal   ? hst->GetHistory(stage,const_cast<char*>("LALVERSION"))   : "";   
      if(log!="") pr_lal_ver = log;                   
      log = hcfg_branch ? hst->GetHistory(stage,const_cast<char*>("CWB_CONFIG_BRANCH"))   : "";   
      if(log!="") pr_cfg_ver = log;                   
      log = hcfg_tag    ? hst->GetHistory(stage,const_cast<char*>("CWB_CONFIG_TAG"))      : "";   
      if(log!="") pr_cfg_ver = log;                   
      log = hcfg_diff   ? hst->GetHistory(stage,const_cast<char*>("CWB_CONFIG_DIFF"))     : "";   
      if(log!="") pr_cfg_ver += "/"+log;                   
      // get cuts parameters used in setCuts
      if(stageName=="CUTS") merge_cuts=hst->GetHistory(stage,const_cast<char*>("PARAMETERS"));
    }
    if(pr_wat_ver=="") pr_wat_ver="n/a";
    if(pr_wat_git=="") pr_wat_git="n/a";
    if(pr_framelib_ver=="") pr_framelib_ver="n/a";
    if(pr_root_ver=="") pr_root_ver="n/a";
    if(pr_lal_ver=="") pr_lal_ver="n/a";
    if(pr_cfg_ver=="") pr_cfg_ver="n/a";
    delete stageList;
    delete hst;
  }
  wfile->Close();

  if(pp_jet_benckmark && merge_jobList.size()>0) {
    // Create jet_benchmark.png (estimation time plot)
    char max_opt[64]="";
    if(pp_jet_benckmark>0) sprintf(max_opt,"--max %d",pp_jet_benckmark);
    sprintf(cmd,"export CWB_BENCH_OPTS=");
    sprintf(cmd,"%s'--bench jet --plot hist %s --save %s/jet_benchmark.png';",cmd,max_opt,odir);
    sprintf(cmd,"%s root -n -l -b ",cmd);
    sprintf(cmd,"%s ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE}",cmd);
    sprintf(cmd,"%s ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE}",cmd);
    sprintf(cmd,"%s ${CWB_MACROS}/cwb_condor_benchmark.C",cmd);
    cout << cmd << endl;
    gSystem->Exec(cmd);	
  }
  if(pp_mem_benckmark && merge_jobList.size()>0) {
    // Create mem_benchmark.png (maximum memory)
    char max_opt[64]="";
    if(pp_mem_benckmark>0) sprintf(max_opt,"--max %d",pp_mem_benckmark);
    sprintf(cmd,"export CWB_BENCH_OPTS=");
    sprintf(cmd,"%s'--bench mem --plot hist %s --save %s/mem_benchmark.png';",cmd,max_opt,odir);
    sprintf(cmd,"%s root -n -l -b ",cmd);
    sprintf(cmd,"%s ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE}",cmd);
    sprintf(cmd,"%s ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE}",cmd);
    sprintf(cmd,"%s ${CWB_MACROS}/cwb_condor_benchmark.C",cmd);
    cout << cmd << endl;
    gSystem->Exec(cmd);	
  }

  // Create the job_status.html file
  int max_jobs = 0; 
  for(int i=0;i<condor_jobList.size();i++)      // find max job id in the condor dag list
      if(condor_jobList[i]>max_jobs) max_jobs=condor_jobList[i];
  if(max_jobs==0) max_jobs = merge_jobList.size();	// online mode
  vector<int> jobStatus(max_jobs);
  vector<TString> jobFile(max_jobs);
  for (int i=0;i<max_jobs;i++) jobStatus[i]=-1; // excluded jobs           tag =-1
  for (int i=0;i<condor_jobList.size();i++) 	// jobs in the dag file    tag = 0  
      jobStatus[condor_jobList[i]-1]=0;  
  for (int i=0;i<merge_jobList.size();i++) {	// jobs in the merge file  tag >=1 
      if(merge_jobList[i]>max_jobs) {
        cout << "cwb_mkhtml_header.C - Error : number of jobs in merged list is > job in condor list" << endl;
        gSystem->Exit(1);
      }
      jobStatus[merge_jobList[i]-1]+=1;  
      jobFile[merge_jobList[i]-1]=merge_jobFileList[i];  
  }
  char jstatus_file[1024];sprintf(jstatus_file,"%s/job_status.html",odir); 
  ofstream out;
  out.open(jstatus_file,ios::out);		// create job_status.html
  char oline[1024];
  out << "<html>" << endl;
  char jstatus[1024];
  sprintf(jstatus,"<font color=\"black\">Job Status : </font>");
  sprintf(jstatus,"%s<font color=\"red\">%s</font>",jstatus,merge_job_list_str);
  sprintf(jstatus,"%s<font color=\"black\">/%s - </font>",jstatus,condor_job_list_str);
  sprintf(jstatus,"%s<font color=\"blue\">%s</font>",jstatus,percentage_job_list_str);
  out << "<font style=\"font-weight:bold;\"><center><p><h2>"
      << jstatus << "</h2><p><center></font>" << endl;
  if(((pp_jet_benckmark)||(pp_mem_benckmark)) && (merge_jobList.size()>0)) {
    out << "<hr><br>" << endl;
    out << "<table border=0 cellpadding=2 align=\"center\">" << endl;
    out << "<tr align=\"center\">" << endl;
    if(pp_jet_benckmark) 
      out << "<td><a href=\"jet_benchmark.png\"><img src=\"jet_benchmark.png\" width=450></a></td>" << endl;
    if(pp_mem_benckmark) 
      out << "<td><a href=\"mem_benchmark.png\"><img src=\"mem_benchmark.png\" width=450></a></td>" << endl;
    out << "</tr>" << endl;
    out << "<tr align=\"center\">" << endl;
    if(pp_jet_benckmark>0) 
      out << "<td>max time < " << pp_jet_benckmark << " hours </td>" << endl;
    else
      out << "<td>max hour = auto </td>" << endl;
    if(pp_mem_benckmark>0) 
      out << "<td>max mem < " << pp_mem_benckmark << " MB </td>" << endl;
    else
      out << "<td>max mem = auto </td>" << endl;
    out << "</tr>" << endl;
    out << "</table>" << endl;
  }
  out << "<hr><br>" << endl;
  //out << "<a>(<td><font color=\"blue\"> blue </font></td> : used in the report )</a>" << endl;
  out << "<a><td><font color=\"black\"> List of unfinished jobs </font></td></a>" << endl;
  out << "<br><br>" << endl;
  out << "<table border=0 cellpadding=2 align=\"center\">" << endl;
/*
  out << "<tr align=\"center\">" << endl;
  out << "<td>jobID</td>" << endl;
  if(simulation) out << "<td>factors</td>" << endl;
  out << "<td>File Name</td>" << endl;
  out << "</tr>" << endl;
*/
  if(pp_job_status) {				// write job_status.html
    for (int i=0;i<max_jobs;i++) {		
      if(jobStatus[i]==0) {			// jobs not done
        int ib=i; int ie=i;
        // extract range of jobs in the same status
        for(int j=ib;j<max_jobs;j++) if(jobStatus[j]==0) ie=j; else break;
        i=ie;
        out << "<tr align=\"center\">" << endl;
        if(ib==ie)
          sprintf(oline,"<td><font color=\"blue\">%d</font></td>",ib+1);
        else 
          sprintf(oline,"<td><font color=\"blue\">%d:%d</font></td>",ib+1,ie+1);
        out << oline << endl;
        if(simulation) {
          sprintf(oline,"<td></td>");
          out << oline << endl;
        }
        sprintf(oline,"<td></td>");
        out << oline << endl;
        sprintf(oline,"<td align=\"left\"><font color=\"red\">%s</font></td>","unfinished/excluded");
        out << oline << endl;
        out << "</tr>" << endl;
      }
      if(simulation&&(jobStatus[i]>=1)&&(jobStatus[i]<nfactor)) {	// sim jobs not completed
        TString status = "unfinished";
        out << "<tr align=\"center\">" << endl;
        sprintf(oline,"<td><font color=\"black\">%d</font></td>",i+1);
        out << oline << endl;
        sprintf(oline,"<td><font color=\"red\">%d</font>/%d</td>",jobStatus[i],nfactor);
        out << oline << endl;
        TObjArray* token = jobFile[i].Tokenize(TString('_'));
        TObjString* stoken =(TObjString*)token->At(token->GetEntries()-2);
        jobFile[i].ReplaceAll("_"+stoken->GetString()+"_job","_*_job"); 
        if(token) delete token;
        sprintf(oline,"<td align=\"left\"><font color=\"blue\">%s</font></td>",jobFile[i].Data());
        out << oline << endl;
        sprintf(oline,"<td align=\"left\"><font color=\"red\">%s</font></td>",status.Data());
        out << oline << endl;
        out << "</tr>" << endl;
      }
    }
  }
  out << "</table>" << endl;
  out << "</html>" << endl;
  out.close();

  // convert user macro into html and copy to report dir
  THtml html;
  html.SetEtcDir(gSystem->ExpandPathName("$HOME_WAT/html/etc/html"));
  html.SetProductName("CWB");
  TString html_input_dir="$CWB_TOOLBOX:$CWB_GWAT:$CWB_HISTORY:$HOME_CWB";
  html_input_dir+=":$CWB_STFT:$CWB_BICO:$HOME_WAT/wat:$ROOTSYS/include";
  html_input_dir+=":"+TString(odir);
  html.SetInputDir(html_input_dir.Data());

  sprintf(cmd,"cp %s/tools/config.csh %s/%s/config.csh",home_wat.Data(),pp_dir,pp_data_dir);
  gSystem->Exec(cmd);

  // redirect stderr to /dev/null to getrid of messages produced by html.Convert
  fpos_t pos;
  fflush(stderr); fgetpos(stderr, &pos);
  int fd = dup(fileno(stderr));
  freopen("/dev/null", "w", stderr);

  char ftitle[1024];

  sprintf(ftitle,"<h2 class=\"convert\" align=\"center\"> %s </h2>","cwb_parameters.C");
  html.Convert(cwb_parameters_file.Data(),ftitle,odir);
  TString cwb_uparameters_path = TString(work_dir)+"/"+cwb_uparameters_file;
  sprintf(ftitle,"<h2 class=\"convert\" align=\"center\"> %s </h2>","user_parameters.C");
  html.Convert(cwb_uparameters_path.Data(),ftitle,odir);
  sprintf(ftitle,"<h2 class=\"convert\" align=\"center\"> %s </h2>","cwb_pparameters.C");
  html.Convert(cwb_pparameters_file.Data(),ftitle,odir);
  TString cwb_upparameters_path = TString(work_dir)+"/"+cwb_upparameters_file;
  sprintf(ftitle,"<h2 class=\"convert\" align=\"center\"> %s </h2>","user_pparameters.C");
  html.Convert(cwb_upparameters_path.Data(),ftitle,odir);

  // restore the stderr output
  fflush(stderr); dup2(fd, fileno(stderr)); close(fd);
  clearerr(stderr); fsetpos(stderr, &pos);

  // create plugin html
  gRandom->SetSeed(0);
  int rnID = int(gRandom->Rndm(13)*1.e9);
  UserGroup_t* _uinfo = gSystem->GetUserInfo();
  TString _uname = _uinfo->fUser;
  char tdir[1024];
  sprintf(tdir,"/dev/shm/%s/%d",_uname.Data(),rnID);
  gSystem->Exec(TString("mkdir -p ")+TString(tdir));

  if(TString(plugin.GetName()).Sizeof()>1) {
    char pluginName[1024]="";
    sprintf(pluginName,"%s/CWB_Plugin.C",tdir);
    plugin.SaveSource(pluginName);
    sprintf(ftitle,"<h2 class=\"convert\" align=\"center\"> %s </h2>","CWB_Plugin.C");
    html.Convert(pluginName,ftitle,odir);
    gSystem->Exec(TString("rm -f ")+TString(pluginName));
  }
  if(TString(configPlugin.GetName()).Sizeof()>1) {
    char configPluginName[1024]="";
    sprintf(configPluginName,"%s/CWB_configPlugin.C",tdir);
    configPlugin.SaveSource(configPluginName);
    sprintf(ftitle,"<h2 class=\"convert\" align=\"center\"> %s </h2>","CWB_configPlugin.C");
    html.Convert(configPluginName,ftitle,odir);
    gSystem->Exec(TString("rm -f ")+TString(configPluginName));
  }

  ofstream hout;
  char ofile[1024];
  sprintf(ofile,"%s/header.html", pp_dir);
  cout << "make header file : " << ofile << endl;
  hout.open(ofile,ios::out);
  if (!hout.good()) {cout << "Error Opening File : " << ofile << endl;exit(1);}

  TString Rho_LF,RhOut,RhAcor,svED,sPEN,sIFAR,sWIN,sHRSS,sHRSS_LF,sfLow,sfHigh,sfResample;
  TString sdelta,sgamma,ssubnet,sbpp,snetRHO,snetCC;
  TString sT_cor,sT_cut,sT_scc,spp_vetoes;
  TString spp_irho,spp_inetcc;
  TString ssimulation;

  char s[128];
  sprintf(s,"%g",delta);  sdelta=s;
  sprintf(s,"%g",GAMMA());sgamma=s;
  if (TString(analysis)=="2G")
   {
     sprintf(s,"%.2f",subnet);
     ssubnet = s;
   }
  else ssubnet = "-";
  //sprintf(s,"%.2f",subnet); ssubnet = TString(analysis)=="2G" ? s : "-";
  sprintf(s,"%.6f",bpp);    sbpp=s;
  sprintf(s,"%.2f",netRHO); snetRHO=s;
  sprintf(s,"%.2f",netCC);  snetCC=s;
  sprintf(s,"%.2f",T_cor);  sT_cor=s;
  sprintf(s,"%.2f",T_cut);  sT_cut=s;
  sprintf(s,"%.2f",T_out);  RhOut=s;
  sprintf(s,"%.2f",T_acor); RhAcor=s;
  sprintf(s,"%.2f",T_vED);  svED=s;
  sprintf(s,"%.2f",T_scc);  sT_scc=s;
  sprintf(s,"%.2f",T_pen);  sPEN=s;
  sprintf(s,"%.2f",T_ifar); sIFAR=s;
  sprintf(s,"%.2f",T_win);  sWIN=s;
  sprintf(s,"%.2f",fLow);   sfLow=s;
  sprintf(s,"%.2f",fHigh);  sfHigh=s;
  sprintf(s,"%d",pp_irho);  spp_irho=s;
  sprintf(s,"%d",pp_inetcc); spp_inetcc=s;
  sprintf(s,"%d",simulation); ssimulation=s;

  // search type
  char ssearch[32]="";
  if(optim) sprintf(ssearch,"SRA(%c)",SEARCH());
  else      sprintf(ssearch,"MRA(%c)",SEARCH());
  if(pattern==0)                sprintf(ssearch,"%s:Packet(0)",ssearch);
  if((pattern!=0 && pattern<0)) sprintf(ssearch,"%s:Packet(%d)",ssearch,pattern);
  if((pattern!=0 && pattern>0)) sprintf(ssearch,"%s:Packet(+%d)",ssearch,pattern);

  // user frequencies band cuts
  char pp_fbandcuts[1024]="";
  for(int j=0;j<nFCUT;j++) {
    sprintf(pp_fbandcuts,"%s [%3.1f:%3.1f]",pp_fbandcuts,lowFCUT[j],highFCUT[j]);
  }

  // post production vetoes
  spp_vetoes="";
#ifdef CAT2_VETO
  spp_vetoes+="cat2 ";
#endif
#ifdef HVETO_VETO
  spp_vetoes+="hveto ";
#endif
#ifdef CAT3_VETO
  spp_vetoes+="cat3 ";
#endif
#ifdef PEM_VETO
  spp_vetoes+="pem ";
#endif

  // data rate after resample
  double sRate = fResample>0 ? fResample>>levelR : inRate>>levelR;	
  sprintf(s,"%.2f",sRate);  sfResample=s;

  // compute sky resolution (deg^2)
  char sSkyMapRes[64];
  if(healpix) {	// healpix sky fragmentation
    int npix = 12*pow(4,healpix);
    double sphere_solid_angle = 4*TMath::Pi()*pow(180/TMath::Pi(),2);
    double skyres = sphere_solid_angle/npix; 
    sprintf(sSkyMapRes,"healpix - order = %d - skyres = %2.3f (deg^2)",(int)healpix,skyres);
  } else {	// cwb built in sky fragmentation
    double skyres = angle*angle;
    sprintf(sSkyMapRes,"built-in - skyres = %2.3f (deg^2)",angle*angle);
  }

  char sTFres[128]="";
  for(int n=l_high;n>=l_low;n--) 
    sprintf(sTFres,"%s %gx(1/%g)",sTFres,(sRate/2)/TMath::Power(2,n),sRate/TMath::Power(2,n));
 
  ifstream in;
  in.open(html_header_template,ios::in);
  if (!in.good()) {cout << "Error Opening File : " << html_header_template << endl;exit(1);}

  char pp_framelib_ver[32]; sprintf(pp_framelib_ver,"%f",FRAMELIB_VERSION);
  const char *pp_root_ver = gROOT->GetVersion();
  CWB::mdc MDC;

  TString pp_cfg_ver = GetGitInfos("branch","$CWB_CONFIG");
  pp_cfg_ver += GetGitInfos("tag","$CWB_CONFIG");
  if(GetGitInfos("diff","$CWB_CONFIG")!="") pp_cfg_ver += "/M";  

  char istring[1024];
  while (1) {
    in.getline(istring,1024);
    if (!in.good()) break;
    TString ostring(istring);
    ostring.ReplaceAll("FLOW",sfLow);
    ostring.ReplaceAll("FHIGH",sfHigh);
    ostring.ReplaceAll("SKYMAPRES",sSkyMapRes);
    ostring.ReplaceAll("TFRES",sTFres);
    ostring.ReplaceAll("FRESAMPLE",sfResample);
    ostring.ReplaceAll("RUN_LABEL",RunLabel);
    ostring.ReplaceAll("ST_COR",sT_cor);
    ostring.ReplaceAll("RHO_LF",Rho_LF);
    ostring.ReplaceAll("ST_CUT",sT_cut);
    ostring.ReplaceAll("RH_OUT",RhOut);
    ostring.ReplaceAll("RH_ACOR",RhAcor);
    ostring.ReplaceAll("svED",svED);
    ostring.ReplaceAll("ST_SCC",sT_scc);
    ostring.ReplaceAll("sPEN",sPEN);
    ostring.ReplaceAll("sIFAR",sIFAR);
    ostring.ReplaceAll("sWIN",sWIN);
    ostring.ReplaceAll("sHRSS_LF",sHRSS_LF);
    ostring.ReplaceAll("sHRSS",sHRSS);
    ostring.ReplaceAll("sPP_VETOES",spp_vetoes);
//    ostring.ReplaceAll("SUBTITLE",subtitle);
    ostring.ReplaceAll("BUILD_DATE",wat::Time("now").GetDateString());
    ostring.ReplaceAll("PR_WAT_VER",pr_wat_ver);
    ostring.ReplaceAll("PR_WAT_GIT",pr_wat_git);
    ostring.ReplaceAll("PR_ROOT_VER",pr_root_ver);
    ostring.ReplaceAll("PR_FRAMELIB_VER",pr_framelib_ver);
    ostring.ReplaceAll("PR_LAL_VER",pr_lal_ver);
    ostring.ReplaceAll("PP_WAT_VER",watversion('s'));
    ostring.ReplaceAll("PP_WAT_GIT",watversion('r'));
    ostring.ReplaceAll("PP_ROOT_VER",pp_root_ver);
    ostring.ReplaceAll("PP_FRAMELIB_VER",pp_framelib_ver);
    ostring.ReplaceAll("PR_CFG_VER",pr_cfg_ver);
    ostring.ReplaceAll("PP_CFG_VER",pp_cfg_ver);
    if(GetLALVersion()!="") {
      ostring.ReplaceAll("PP_LAL_VER",GetLALVersion());
    } else {
      ostring.ReplaceAll("PP_LAL_VER","");
    }
    ostring.ReplaceAll("PIPELINE",analysis);
    ostring.ReplaceAll("SEARCH",ssearch);
    ostring.ReplaceAll("SIMULATION",ssimulation);
    ostring.ReplaceAll("DELTA",sdelta);
    ostring.ReplaceAll("GAMMA",sgamma);
    ostring.ReplaceAll("SUBNET",ssubnet);
    ostring.ReplaceAll("BPP",sbpp);
    ostring.ReplaceAll("PP_IRHO",spp_irho);
    ostring.ReplaceAll("PP_INETCC",spp_inetcc);
    ostring.ReplaceAll("NETRHO",snetRHO);
    ostring.ReplaceAll("NETCC",snetCC);
    ostring.ReplaceAll("PP_FREQ_BAND_CUTS",pp_fbandcuts);
    ostring.ReplaceAll("MERGE_CUTS",merge_cuts);
    ostring.ReplaceAll("TITLE",title);
    ostring.ReplaceAll("SITE_CLUSTER",site_cluster);
    ostring.ReplaceAll("WORK_DIR",work_dir);
    ostring.ReplaceAll("PP_DATA_DIR",pp_data_dir);
    ostring.ReplaceAll("JOB_MERGE",merge_job_list_str);
    ostring.ReplaceAll("JOB_CONDOR",condor_job_list_str);
    ostring.ReplaceAll("JOB_PERCENTAGE",percentage_job_list_str);
    if (T_scc>0) {
      ostring.ReplaceAll("<!--CUT_SCC","");
      ostring.ReplaceAll("CUT_SCC-->","");
    }
    if (T_pen>0) {
      ostring.ReplaceAll("<!--CUT_PEN","");
      ostring.ReplaceAll("CUT_PEN-->","");
    }
    if (T_ifar>0) {
      ostring.ReplaceAll("<!--CUT_IFAR","");
      ostring.ReplaceAll("CUT_IFAR-->","");
    }
    if(simulation) {
      ostring.ReplaceAll("<!--CUT_WIN","");
      ostring.ReplaceAll("CUT_WIN-->","");
    }
    ostring.ReplaceAll("</html>","");
    if (T_vED>0) {
      ostring.ReplaceAll("<!--CUT_NED","");
      ostring.ReplaceAll("CUT_NED-->","");
    }
    if (spp_vetoes!="") {
      ostring.ReplaceAll("<!--CUT_PP_VETOES","");
      ostring.ReplaceAll("CUT_PP_VETOES-->","");
    }
    if(TString(plugin.GetName()).Sizeof()>1) {
      ostring.ReplaceAll("<!--CWB_PLUGIN","");
      ostring.ReplaceAll("CWB_PLUGIN-->","");
    }
    if(TString(configPlugin.GetName()).Sizeof()>1) {
      ostring.ReplaceAll("<!--CWB_CONF_PLUGIN","");
      ostring.ReplaceAll("CWB_CONF_PLUGIN-->","");
    }
    if (TString(analysis)=="1G") {
      ostring.ReplaceAll("<!--CWB1G_PARAMETERS","");
      ostring.ReplaceAll("CWB1G_PARAMETERS-->","");
    } else {
      ostring.ReplaceAll("<!--CWB2G_PARAMETERS","");
      ostring.ReplaceAll("CWB2G_PARAMETERS-->","");
    }
    if(cwb_doc_url!="") {
      ostring.ReplaceAll("<!--CWB_DOC_URL","");
      ostring.ReplaceAll("CWB_DOC_URL-->","");
      ostring.ReplaceAll("XCWB_DOC_URL",cwb_doc_url.Data());
    }
    ostring.ReplaceAll("</html>","");
    hout << ostring.Data() << endl;
  }
  in.close();
  hout.close();

  exit(0);
}
