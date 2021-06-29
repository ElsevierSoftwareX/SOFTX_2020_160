/*
# Copyright (C) 2019 Gabriele Vedovato
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

#ifndef TOOLFUN_HH
#define TOOLFUN_HH

#include "TGlobal.h"
 
#ifdef _USE_LAL
#include <lal/LALConfig.h>
#include <lal/LALVCSInfoHeader.h>
#endif

#ifdef _USE_EBBH
#include "cvode/cvode.h"
#endif

//_______________________________________________________________________________________
inline Double_t PoissonIFunction(Double_t* x, Double_t* par) {
//
// Int Poisson Function
//

  return par[0]*TMath::PoissonI(x[0],par[1]);
}

//_______________________________________________________________________________________
inline Double_t logNfit(Double_t *x, Double_t *par) {
//
// log normal function used for efficiency curves fits
//

  double y = (log10(x[0])-par[0]);
  if(par[4]) y=-y;
  double s = y<0 ? par[1]*exp(y*par[2]) : par[1]*exp(y*par[3]);

  if(y>0) {
    if(par[3]>1./y) {s = par[1]*par[3]*exp(1.); y = 1.;}
    y = s>0 ? fabs(y/s) : 100.;
    return 1-TMath::Erfc(y)/2;
  }

  if(y<0) {
    y = s>0 ? fabs(y/s) : 100.;
    return TMath::Erfc(y)/2;
  }
  return 0.5;
}

//_______________________________________________________________________________________
inline int DrawMDC(TString sel="", TString cut="", TString opt="") {
//
// Draw mdc wavetree 
// Input parameters are the same as TTree::Draw
//

  TTree* tree = (TTree *) gROOT->FindObject("mdc"); 
  if(tree) {
     if(sel!="") tree->Draw(sel,cut,opt);
     else tree->StartViewer();
  } else {cout << "DrawMDC : mdc tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetSelectedRows();
}

//_______________________________________________________________________________________
inline int DrawWAVE(TString sel="", TString cut="", TString opt="") {
//
// Draw wave wavetree 
// Input parameters are the same as TTree::Draw
//

  TTree* tree = (TTree *) gROOT->FindObject("waveburst"); 
  if(tree) {
     if(sel!="") tree->Draw(sel,cut,opt);
     else tree->StartViewer();
  } else {cout << "DrawWAVE : waveburst tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetSelectedRows();
}

//_______________________________________________________________________________________
inline int DrawLIVE(TString sel="", TString cut="", TString opt="") {
//
// Draw liveTime 
// Input parameters are the same as TTree::Draw
//

  TTree* tree = (TTree *) gROOT->FindObject("liveTime"); 
  if(tree) {
     if(sel!="") tree->Draw(sel,cut,opt);
     else tree->StartViewer();
  } else {cout << "DrawLIVE : liveTime tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetSelectedRows();
}

//_______________________________________________________________________________________
inline int ScanMDC(TString sel="", TString cut="", TString opt="") {
//
// Scan mdc wavetree 
// Input parameters are the same as TTree::Draw
//

  TTree* tree = (TTree *) gROOT->FindObject("mdc"); 
  if(tree) {
     if(sel!="") tree->Scan(sel,cut,opt);
     else tree->StartViewer();
  } else {cout << "ScanMDC : mdc tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetSelectedRows();
}

//_______________________________________________________________________________________
inline int ScanWAVE(TString sel="", TString cut="", TString opt="") {
//
// Scan wave wavetree 
// Input parameters are the same as TTree::Draw
//

  TTree* tree = (TTree *) gROOT->FindObject("waveburst"); 
  if(tree) {
     if(sel!="") tree->Scan(sel,cut,opt);
     else tree->StartViewer();
  } else {cout << "ScanWAVE : waveburst tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetSelectedRows();
}

//_______________________________________________________________________________________
inline int ScanLIVE(TString sel="", TString cut="", TString opt="") {
//
// Scan liveTime 
// Input parameters are the same as TTree::Draw
//

  TTree* tree = (TTree *) gROOT->FindObject("liveTime"); 
  if(tree) {
     if(sel!="") tree->Scan(sel,cut,opt);
     else tree->StartViewer();
  } else {cout << "ScanLIVE : liveTime tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetSelectedRows();
}

//_______________________________________________________________________________________
inline int PrintMDC() {
//
// Print mdc tree 
//

  TTree* tree = (TTree *) gROOT->FindObject("mdc"); 
  if(tree) tree->Print();
  else {cout << "PrintMDC : mdc tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetEntries();
}

//_______________________________________________________________________________________
inline int PrintWAVE() {
//
// Print wave tree 
//

  TTree* tree = (TTree *) gROOT->FindObject("waveburst"); 
  if(tree) tree->Print();
  else {cout << "PrintWAVE : waveburst tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetEntries();
}

//_______________________________________________________________________________________
inline int PrintLIVE() {
//
// Print live tree 
//

  TTree* tree = (TTree *) gROOT->FindObject("liveTime"); 
  if(tree) tree->Print();
  else {cout << "PrintLIVE : liveTime tree not found !!!" << endl;return 0;}
  return (Int_t)tree->GetEntries();
}

//_______________________________________________________________________________________
inline TTree* GetMDC() {
//
// Get mdc tree pointer
//

  return (TTree *) gROOT->FindObject("mdc"); 
}

//_______________________________________________________________________________________
inline TTree* GetWAVE() {
//
// Get wave tree pointer
//

  return (TTree *) gROOT->FindObject("waveburst"); 
}

//_______________________________________________________________________________________
inline TTree* GetLIVE() {
//
// Get live tree pointer
//

  return (TTree *) gROOT->FindObject("liveTime"); 
}

//_______________________________________________________________________________________
inline TString GetFileLabel(TTree* tree, int run, int lag, int slag, double segEdge, TString psfix) {
//
// return a label 
// label = "START_STOP_psfix_slagSLAG_lagLAG_1_jobRUN"
//

  char cut[1024];
  sprintf(cut,"run==%i",run);
  char draw[256];
  sprintf(draw,"start[0]-left[0]:stop[0]+right[0]");
  tree->Draw(draw,cut,"goff");
  if (tree->GetSelectedRows()==0) {cout <<"GetFileLabel : No events" << endl; exit(1);}
  double* segstart=tree->GetV1();
  double* segstop=tree->GetV2();
  char label[1024];
  int start = TMath::Nint(segstart[0])+segEdge;
  int stop  = TMath::Nint(segstop[0])-TMath::Nint(segstart[0])-2*segEdge;
  sprintf(label,"%i_%i_%s_slag%i_lag%i_%i_job%i",
          start,stop, psfix.Data(),slag,lag,1,run);
  return label;
}

//_______________________________________________________________________________________
inline void GetProcInfo(TString str="") {
//
// return pocess infos
// Ex : virtual : 279 (mb)  rss  : 44 (mb)
// 

   gSystem->Exec("date");
   TString s;
   FILE *f = fopen(Form("/proc/%d/statm", gSystem->GetPid()), "r");
   s.Gets(f);
   Long_t total, rss;
   sscanf(s.Data(), "%ld %ld", &total, &rss);
   cout << str.Data() << " virtual : " <<  total * 4 / 1024 
        << " (mb)  rss  : " <<  rss * 4 / 1024 << " (mb)" << endl;
   fclose(f);
   return;
}

//_______________________________________________________________________________________
inline double GetStart(TTree* tree, int nifo, int run, float rho, 
                       double time0, double time1, TString analysis, int irho) {
//
// 
//

  char cut[1024];
  sprintf(cut,"run==%i && abs(rho[%d]-%f)<0.001 && abs(time[0]-%f)<0.0001 && abs(time[1]-%f)<0.0001",
               run,irho,rho,time0,time1);
  char draw[256];
  sprintf(draw,"start[%i]",nifo);
  tree->Draw(draw,cut,"goff");
  if (tree->GetSelectedRows()>1)  {
    cout << "GetStart : Too many events @ the same time : " << tree->GetSelectedRows() << endl; 
    cout << cut << endl;
    exit(1);
  }
  if (tree->GetSelectedRows()==0) {cout <<"GetStart : No events" << endl; exit(1);}
  double* start=tree->GetV1();
  return start[0];
}

//_______________________________________________________________________________________
inline TString GetString(TTree* tree, int run, int lag, TString psfix) {
//
//
//

  char cut[1024];
  sprintf(cut,"run==%i",run);
  char draw[256];
  sprintf(draw,"start[0]-left[0]:stop[0]+right[0]");
  tree->Draw(draw,cut,"goff");
  double* segstart=tree->GetV1();
  double* segstop=tree->GetV2();
  char s[1024];
  sprintf(s,"%i_%i_%s_%i_id%i",TMath::Nint(segstart[0])+8,
          TMath::Nint(segstop[0])-TMath::Nint(segstart[0])-16,psfix.Data(),lag,run);
  TString st(s);
  return st;
}

//_______________________________________________________________________________________
inline void CheckAnalysis() {
//
// check consistency of analysis CINT variable with the environment CWB_ANALYSIS
//

  // import analysis value from CINT
  char analysis[8];
  TGlobal* global = gROOT->GetGlobal("analysis",true);                               
  //if(global!=NULL) analysis = *(char*)global->GetAddress();                         
  if(global!=NULL) {                         
    memcpy((void*)&analysis,(void*)global->GetAddress(),sizeof(analysis)*sizeof(char));
  } else {cout << "CheckAnalysis : analysis is not defined !!!" << endl;exit(1);}

  if(gSystem->Getenv("CWB_ANALYSIS")!=NULL) {
    if(TString(analysis)!=TString(gSystem->Getenv("CWB_ANALYSIS"))) {
      cout << endl;
      cout << "Error : analysis=" << analysis;
      cout << " is inconsistent with the environment CWB_ANALYSIS="
           << gSystem->Getenv("CWB_ANALYSIS") << endl;
      cout << "        check analysis parameter in user_parameters.C" << endl;
      cout << "        check CWB_ANALYSIS env in watenv setup" << endl;
      cout << "        use 'cwb_setpipe [1G/1g/1,2G/2g/2]' command to switch analysis type [1G/2G]" << endl;
      cout << endl;
      gSystem->Exit(1);
    }
  } else {
    cout << "" << endl;
    cout << "CWB_ANALYSIS env not defined in watenv setup" << endl;
    cout << "" << endl;
    cout << " - add in watenv setup the following statement" << endl;
    cout << "" << endl;
    cout << "   setenv CWB_ANALYSIS        '1G'       # 1G  analysis " << endl;
    cout << "   or" << endl;
    cout << "   setenv CWB_ANALYSIS        '2G'       # 2G  analysis " << endl;
    cout << "" << endl;
    cout << " - use 'cwb_setpipe [1G/1g/1,2G/2g/2]' command to switch analysis type [1G/2G]" << endl;
    cout << "" << endl;
    gSystem->Exit(1);
  }
}

//_______________________________________________________________________________________
inline TString GetGitInfos(TString option="path", TString igit_path="$CWB_CONFIG") {
//
// Get git version infos
//
//
// input  option        : diff,tag,branch,tag
//        igit_path   	: input git path, default env variable=$CWB_CONFIG
//
// output s             : return git version infos
//

  char git_path[1024] = "";
  if(igit_path.BeginsWith("$")) {	// get enviromantal variable
    igit_path.ReplaceAll("$","");
    if(gSystem->Getenv(igit_path.Data())!=NULL) {
      strcpy(git_path,gSystem->Getenv(igit_path.Data()));
    }
  } else {
    strcpy(git_path,igit_path.Data());
  }

  if(option=="path") return git_path;

  if(TString(git_path)=="") return "";

  // Check if path exists
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(git_path,&id,&size,&flags,&mt);
  if(estat!=0) return "";

  char cmd[1024];

  if(option=="diff") {
    sprintf(cmd,"git -C %s  diff HEAD", git_path);
  }

  if(option=="hash") {
    sprintf(cmd,"git -C %s  rev-parse HEAD", git_path);
  }

  if(option=="branch") {
    sprintf(cmd,"git -C %s branch -l | grep '*'", git_path);
  }

  if(option=="tag") {
    sprintf(cmd,"git -C %s tag -l --points-at HEAD", git_path);
  }

  if(option=="url") {
    sprintf(cmd,"git -C %s config --get remote.origin.url", git_path);
  }

  // redirect stderr to /dev/null to getrid of messages produced by cmd
  fpos_t poserr; fflush(stderr); fgetpos(stderr, &poserr);
  int fderr = dup(fileno(stderr)); freopen("/dev/null", "w", stderr);
  // redirect stdout to /dev/null to getrid of messages produced by cmd
  fpos_t posout; fflush(stdout); fgetpos(stdout, &posout);
  int fdout = dup(fileno(stdout)); freopen("/dev/null", "w", stdout);

  int  the_err = gSystem->Exec(cmd);

  // restore the stderr output
  fflush(stderr); dup2(fderr, fileno(stderr)); close(fderr);
  clearerr(stderr); fsetpos(stderr, &poserr);
  // restore the stdout output
  fflush(stdout); dup2(fdout, fileno(stdout)); close(fdout);
  clearerr(stdout); fsetpos(stdout, &posout);


  if(the_err==0) { 
    TString the_output = gSystem->GetFromPipe(cmd);
    if(option=="branch") {
      the_output.ReplaceAll("*","");
      the_output.ReplaceAll(" ","");
      if(the_output=="(nobranch)") the_output="";
    }
    return the_output;
  }

  return "";
}

//_______________________________________________________________________________________
inline void PrintLogoCWB(TString LALVersion="", TString cwb_library_path="$HOME_WAT", TString cwb_config_path="$CWB_CONFIG") {
//
// Print CWB logo
//

  // check operative system
  // ------------------------------------------------------------------------------------
  TString OS = "";
  if(TString(gSystem->GetBuildArch()).Contains("linux")) OS="Linux";
  if(TString(gSystem->GetBuildArch()).Contains("macos")) OS="Darwin";

  char line[256];                                       
  TString sline;                                        
  const char *root_version = gROOT->GetVersion();       
  cout << endl;           
  cout << "     ****************************************************************************" << endl;
  cout << "     *                                                                          *" << endl;
  cout << "     *             W E L C O M E  to  C W B                                     *" << endl;
  cout << "     *                                                                          *" << endl;
  // print WAT version                                                        
  sprintf(line,"     *             cWB    Version  %s (XIFO=%s)",watversion('s'),watversion('i'));      
  sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   

  // print library GIT branch                                                       
  //TString wat_branch = GetGitInfos("branch",cwb_library_path);
  //TString wat_tag    = GetGitInfos("tag",cwb_library_path);
  //TString wat_diff   = GetGitInfos("diff",cwb_library_path);
  TString wat_branch = watversion('b');
  TString wat_tag    = watversion('g');
  TString wat_diff   = "";
  if(wat_branch!="") {
    if(wat_diff!="") wat_branch=wat_branch+"/M";
    sprintf(line,"     *                     Branch  %s",wat_branch.Data());         
  } else if(wat_tag!="") { 
    if(wat_diff!="") wat_branch=wat_branch+"/M";
    sprintf(line,"     *                        Tag  %s",wat_tag.Data());         
  } else {
    sprintf(line,"     *                             %s","Undefined");         
  }
  sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   

  // print library GIT hash                                                       
  //sprintf(line,"     *                       Hash  %s",GetGitInfos("hash",cwb_library_path).Data());      
  sprintf(line,"     *                       Hash  %s",watversion('r'));      
  sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   
  sprintf(line,"     *                 Short Hash  %s",watversion('R'));      
  sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   
  cout << "     *                                                                          *" << endl;

  // print LAL version                                                        
  if(LALVersion!="") {                                         
    sprintf(line,"     *             LAL     Branch  %s",LALVersion.Data());                                
    sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                 
  }                                                                           

  // FRAMELIB_VERSION is not defined in Darwin OS (to be checked)             
  if(OS=="Linux") {                                                           
    sprintf(line,"     *             FRLIB  Version  %2.2f",FRAMELIB_VERSION);
    sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                 
  }                                                                           

  // print CVODE version
  if(WAT::USE_EBBH()) {     
#ifdef SUNDIALS_PACKAGE_VERSION 	// SUNDIALS_PACKAGE_VERSION = 2.7.0
    sprintf(line,"     *             CVODE  Version  %s",SUNDIALS_PACKAGE_VERSION);                                
    sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                 
#endif
#ifdef SUNDIALS_VERSION 		// SUNDIALS_VERSION >= 3.2.1
    sprintf(line,"     *             CVODE  Version  %s",SUNDIALS_VERSION);                                
    sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                 
#endif
  }
  cout << "     *                                                                          *" << endl;

  // print ROOT version                                                       
  sprintf(line,"     *             Based  on ROOT  %s",root_version);         
  sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   
  
  cout << "     *                                                                          *" << endl;
  if(WAT::USE_ROOT6())    cout << "     *             ROOT6           ENABLED                                      *" << endl;         
  if(WAT::USE_CPP11())    cout << "     *             CP11            ENABLED                                      *" << endl;         
  if(WAT::USE_ICC())      cout << "     *             ICC             ENABLED                                      *" << endl;         
  if(!WAT::USE_LAL())     cout << "     *             LAL             DISABLED                                     *" << endl;         
  if(!WAT::USE_HEALPIX()) cout << "     *             HEALPix         DISABLED                                     *" << endl;         
  if(!WAT::USE_EBBH())    cout << "     *             EBBH            DISABLED                                     *" << endl;

  if(gSystem->Getenv("_USE_PEGASUS")!=NULL) 
                          cout << "     *             PEGASUS         ENABLED                                      *" << endl;         
  if(gSystem->Getenv("_USE_LSF")!=NULL) 
                          cout << "     *             LSF             ENABLED                                      *" << endl;         
  if(gSystem->Getenv("_USE_OSG")!=NULL) 
                          cout << "     *             OSG             ENABLED                                      *" << endl;         
#if (WAT_VERSION_DEVEL != 0) 
  cout << "     *                                                                          *" << endl;
  sprintf(line,"     *             DEVELOPMENT VERSION : %d",WAT_VERSION_DEVEL);         
  sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   
#endif 
  cout << "     *                                                                          *" << endl;

  // print config GIT branch
  TString cfg_branch = GetGitInfos("branch",cwb_config_path);
  TString cfg_tag    = GetGitInfos("tag",cwb_config_path);
  TString cfg_diff   = GetGitInfos("diff",cwb_config_path);
  bool cfg_status = true;
  if(cfg_branch!="") {
    if(cfg_diff!="") cfg_branch=cfg_branch+"/M";
    sprintf(line,"     *             CONFIG  Branch  %s",cfg_branch.Data());         
  } else if(cfg_tag!="") { 
    if(cfg_diff!="") cfg_branch=cfg_branch+"/M";
    sprintf(line,"     *             CONFIG     Tag  %s",cfg_tag.Data());         
  } else {
    sprintf(line,"     *             CONFIG          %s","Undefined");         
    cfg_status = false;
  }
  sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   

  // print config GIT hash                                                       
  if(cfg_status) {
    sprintf(line,"     *                       Hash  %s",GetGitInfos("hash",cwb_config_path).Data());      
    sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   
  }

  // print data config GIT branch
  TString cwb_data_path = gSystem->ExpandPathName(cwb_config_path+"/DATA");
  TString data_branch = GetGitInfos("branch",cwb_data_path);
  TString data_tag    = GetGitInfos("tag",cwb_data_path);
  TString data_diff   = GetGitInfos("diff",cwb_data_path);
  bool data_status = true;

  char data_name[1024];
  int len;
  if ((len = readlink(cwb_data_path.Data(), data_name, sizeof(data_name)-1)) != -1) data_name[len] = '\0';

  if(data_branch!="") {
    if(data_diff!="") data_branch=data_branch+"/M";
    if(data_branch==cfg_branch) {					// default config DATA
      sprintf(line,"     *             DATA    Branch  %s",data_name);         
    } else {
      sprintf(line,"     *             DATA    Branch  %s",data_branch.Data());         
    }
  } else if(data_tag!="") { 
    if(data_diff!="") data_branch=data_branch+"/M";
    if(data_tag==cfg_tag) {						// default config DATA
      sprintf(line,"     *             DATA       Tag  %s",data_name);         
    } else { 
      sprintf(line,"     *             DATA       Tag  %s",data_tag.Data());         
    }
  } else {
    sprintf(line,"     *             DATA            %s","UNDEFINED");         
    data_status = false;
  }
  sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   

  // print data config GIT hash                                                       
  if(data_status) {
    sprintf(line,"     *                       Hash  %s",GetGitInfos("hash",cwb_data_path).Data());      
    sline=line;sline.Resize(80);sline+="*";cout<<sline<<endl;                   
  }

  cout << "     *                                                                          *" << endl;
  cout << "     ****************************************************************************" << endl;
  cout << endl;                                                               
  cout << watversion('x') << endl;                                            
  cout << "   Compiled on " << watversion('k') << " "                         
       << watversion('m') << " " << watversion('n') << endl;                  
  cout << "                 " << watversion('t') << endl;                     
  cout << endl;                                                               

  return;
}

//_______________________________________________________________________________________
inline int Draw(TChain& rf, string par, string cut, string opt)
{ return rf.Draw(par.c_str(),cut.c_str(),opt.c_str()); }

//_______________________________________________________________________________________
inline int Draw(TTree& rf, string par, string cut, string opt)
{ return rf.Draw(par.c_str(),cut.c_str(),opt.c_str()); }

//_______________________________________________________________________________________
inline TCanvas* DRAW(TChain& rf, string par, string cut, string opt, int iX=600, int iY=600){ 
  TCanvas *c1 = new TCanvas("c","C",0,0,iX,iY);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);
  c1->SetBottomMargin(0.103939);
  rf.Draw(par.c_str(),cut.c_str(),opt.c_str()); 
  return c1;
}

//_______________________________________________________________________________________
inline TCanvas* Draw(TH1F* h, int iX=600, int iY=600){

  TCanvas *c1 = new TCanvas("c","C",0,0,iX,iY);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);
  c1->SetBottomMargin(0.103939);
  h->Draw();

  return c1;
}  

//_______________________________________________________________________________________
inline TCanvas* Draw(TH2F* h, char* opt, int iX=600, int iY=600){

  TCanvas *c1 = new TCanvas("c","C",0,0,iX,iY);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetLogx(kFALSE);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);
  c1->SetBottomMargin(0.103939);
  h->Draw(opt);

  return c1;
}  

//_______________________________________________________________________________________
inline TCanvas* Draw(TH1F** h, int n=1, int m=0, int iX=600, int iY=600){

  TCanvas *c1 = new TCanvas("c","C",0,0,iX,iY);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetLogx(kFALSE);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);
  c1->SetBottomMargin(0.103939);

  int i=0;
  if(!m) {
    (*(h+i))->Draw();
    for(i=1; i<n; i++) (*(h+i))->Draw("same");
  }
  else {
    c1->Divide(m,n);
    for(int i=0; i<n*m; i++) { c1->cd(i+1); (*(h+i))->Draw(); }
  }
  return c1;
}  

//_______________________________________________________________________________________
inline TCanvas* Draw(TH2F** h, char* opt, int n=1, int m=0, int iX=600, int iY=600){

  TCanvas *c1 = new TCanvas("c","C",0,0,iX,iY);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetLogx(kFALSE);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);
  c1->SetBottomMargin(0.103939);

  if(!m) {
    (*h)->Draw(opt);
    for(int i=1; i<n; i++) (*(h+i))->Draw("same");
  }
  else {
    c1->Divide(m,n);
    for(int i=0; i<n*m; i++) { c1->cd(i+1); (*(h+i))->Draw(opt); }
  }
  return c1;
}  

//_______________________________________________________________________________________
inline TCanvas* Draw(TGraphErrors** h, char* opt, int n=1, int m=0, int iX=600, int iY=600){

  TCanvas *c1 = new TCanvas("c","C",0,0,iX,iY);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetLogx(kFALSE);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);
  c1->SetBottomMargin(0.103939);

  if(!m) {
    (*h)->Draw(opt);
    for(int i=1; i<n; i++) (*(h+i))->Draw("P");
  }
  else {
    c1->Divide(m,n);
    for(int i=0; i<n*m; i++) { c1->cd(i+1); (*(h+i))->Draw(opt); }
  }
  return c1;
}  

//_______________________________________________________________________________________
inline TCanvas* Draw(TGraph** h, char* opt, int n=1, int m=0, int iX=600, int iY=600){

  TCanvas *c1 = new TCanvas("c","C",0,0,iX,iY);
  c1->SetBorderMode(0);
  c1->SetFillColor(0);
  c1->SetBorderSize(2);
  c1->SetLogx(kFALSE);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.1517039);
  c1->SetTopMargin(0.0772727);
  c1->SetBottomMargin(0.103939);

  if(!m) {
    (*h)->Draw(opt);
    for(int i=1; i<n; i++) (*(h+i))->Draw("P");
  }
  else {
    c1->Divide(m,n);
    for(int i=0; i<n*m; i++) { c1->cd(i+1); (*(h+i))->Draw(opt); }
  }
  return c1;
}

//_______________________________________________________________________________________
inline double GetLiveTime(TString fliveName, int SLAG, int LAG, int& nlags) {
//
// return live time in sec read from live.txt file produced in the background report
//
// input  fliveName : livetime file name
//        SLAG      : select LAG, SLAG<0 -> any slag
//        LAG       : select LAG, LAG<0  -> any lag
//
// output nlags     : return the selected lags
//
// examples         : LAG=0,SLAG=0   -> select zero LAG live time
//                  : LAG=-1,SLAG=-1 -> select all LAG,SLAG
//                  : LAG=1,SLAG=-1  -> select LAG=-1 and all SLAG
//
// livetime format  : slag=SLAG   lag=LAG   L1:slag= 0.00  lag= 0.00  ...  live=2145137.00
//

  double livetime=0.;
  nlags=0;

  ifstream inliv(fliveName);
  if (!inliv.good()) {cout << "GetLiveTime : Error Opening File : " << fliveName << endl;exit(1);}

  char line[1024];
  while(1) {
    // read line
    inliv.getline(line,256);
    if (!inliv.good()) break;
    TString Oliv(line);
    TObjArray* token = TString(Oliv).Tokenize(TString('='));
    // skip nonzero livetime line
    if(Oliv.Contains("nonzero")) continue;
    // check line format
    if(!TString(line).Contains("slag=")) {cout<<"GetLiveTime : Bad line format : "<<line<<endl<<endl;exit(1);}
    if(!TString(line).Contains("lag="))  {cout<<"GetLiveTime : Bad line format : "<<line<<endl<<endl;exit(1);}
    if(!TString(line).Contains("live=")) {cout<<"GetLiveTime : Bad line format : "<<line<<endl<<endl;exit(1);}
    // last token contains livetime
    TObjString* tok=(TObjString*)token->At(token->GetEntries()-1);
    token = TString(Oliv).Tokenize(TString(' '));
    // get slag number
    TObjString* tokslag=(TObjString*)token->At(0);
    int slag=TString(tokslag->GetString()).ReplaceAll("slag=","").Atoi();
    // get lag number
    TObjString* toklag=(TObjString*)token->At(1);
    int lag=TString(toklag->GetString()).ReplaceAll("lag=","").Atoi();
    // get livetime with lag=LAG & slag=SLAG
    if((lag==LAG||LAG<0) && (slag==SLAG||SLAG<0)) {
      livetime+=TString(tok->GetString()).Atof();
      nlags++;
    }
  }
  inliv.close();

  return livetime;
}

//_______________________________________________________________________________________
inline int ReadInjType(TString ifName, int ntype_max, char set[][128], size_t type[], 
                       char name[][128], double fcentral[], double fbandwidth[]) {
//
// Read the user defined injection types file : ifName
//
// ifName                 // list mdc types file name (ascii)
//                           format : set type name fcentral fbandwidth
// ntype_max              // max number of mdc types
// set[NMDC_MAX][128];    // injection set
// type[NMDC_MAX];        // injection type
// name[NMDC_MAX][128];   // injection name
// fcentral[NMDC_MAX];    // injection central frequencies
// fbandwidth[NMDC_MAX];  // injection bandwidth frequencies

  // open injection file list mdc types
  ifstream in;
  cout << "inj file: " << ifName << endl;
  in.open(ifName.Data(),ios::in);
  if (!in.good()) {cout << "ReadMdcType : Error Opening File : " << ifName << endl;exit(1);}

  int ninj=0;
  char str[1024];
  int fpos=0;
  while (1) {
    fpos=in.tellg();
    in.getline(str,1024);
    if (!in.good() && strlen(str)!=0) {
      cout << endl;
      cout << "ReadMdcType : Error Reading File : " << ifName << endl;
      cout << "at line : \"" << str << "\"" << endl;
      cout << "check if new line is defined at the end of the string" << endl;
      exit(1);
    }
    if (!in.good()) break;
    if((str[0] == '#')||(str[0] == ' ')||(str[0] == 0)) continue;
    in.seekg(fpos, ios::beg);

    in >> set[ninj] >> type[ninj] >> name[ninj]
       >> fcentral[ninj] >> fbandwidth[ninj];

    fpos=in.tellg();
    in.seekg(fpos+1, ios::beg);

    if (!in.good()) break;
    cout << " " << set[ninj] << " " << type[ninj] << " " << name[ninj]
         << " " << fcentral[ninj] << " " << fbandwidth[ninj] << endl;
    ninj++;
    if(ninj>=ntype_max) {
      cout << "ReadMdcType : Error - max allowed injection types must be < "
           << ntype_max << endl;
      exit(1);
    }
  }
  in.close();

  return ninj;
}

//_______________________________________________________________________________________
inline void MakePlotsHtmlTable(ofstream* out, TString title, TString png1, TString png2="") { 
//
// macro used to generate a table with 2 plots in html code
//
// out	 : pinter to the output streamer
// title : table title
// png1  : path of first png file name 
// png2  : path of second png file name 
//
// if png2="" then only one plot is showed 
//

  int psize = 960;
  if(png2.IsDigit()) {
    psize=png2.Atoi();
    png2="";
  } else {
    psize = png2=="" ? 960 : 480;
  }

  *out << "<br></td></tr></table>" << endl;                                                      
  *out << endl;                                                                                  
  *out << "<div align=\"center\"><font color=\"red\"><h2>"<<title<<"</h2></font></div>" << endl;

  *out << "<p><table summary=\"\"><tr align=\"left\"><td valign=\"top\" width=\"50%\"> 	\
           <div align=\"center\"><ul><br/>" << endl;                                               
  *out << "<a class=\"image\" title=\""<<png1<<"\">" << endl;                                    
  *out << "<img src=\""<<png1<<"\" onerror=\"this.style.visibility = 'hidden'\" width=\""<<psize<<"\"> </a>" << endl;                             
  *out << "</br></ul><br><br>" << endl;                                                          
  *out << endl; 
  if(png2=="") return;                                                                                 
  *out << "</div>" << endl;                                                                      
  *out << "<p></td><td valign=\"top\" width=\"50%\"><div align=\"center\"><ul><br/>" << endl;    
  *out << "<a class=\"image\" title=\""<<png2<<"\">" << endl;                                    
  *out << "<img src=\""<<png2<<"\" onerror=\"this.style.visibility = 'hidden'\" width=\""<<psize<<"\"> </a>" << endl;                             
  *out << "</br></ul><br><br>" << endl;                                                          
  *out << endl;                                                                                  
  *out << "</div>" << endl;
}

//_______________________________________________________________________________________
inline void MakePlotsHtmlCellTable(ofstream* out, TString title, TString png1, TString subtitle1="", TString png2="", TString subtitle2="") { 
//
// macro used to generate a cell with a table with 2 plots in html code
//
// out	 : pointer to the output streamer
// title : table title
// subtitle1 : title for cell 1
// png1  : path of first png file name 
// subtitle2 : title for cell 2 
// png2  : path of second png file name 
//
// if png2="" then only one plot is showed 
//

  int h = 380;
  *out << "<font color=\"red\" style=\"font-weight:bold;\"><center><p><h3>"<< title << "</h3><p><center></font>" << endl;
  *out << "<p>" << endl;
  *out << endl;
  *out << "<table cellspacing=\"0\" cellpadding=\"6\" border=\"0\" align=\"center\" width=\"1100\">" << endl;
  *out << "<tr>" << endl;
  *out << "<th>" << subtitle1 << "</th>" << endl;
  *out << "<th>" << subtitle2 << "</th>" << endl;
  *out << "</tr>" << endl;
  *out << "<tr>" << endl;
  *out << "<td align=\"center\"><a target=\"_parent\" href=\""<< png1 <<"\"><img src=\""<< png1 << "\" height=" << h <<"></a></td>";
  if(png2=="") return; 
  *out << "<td align=\"center\"><a target=\"_parent\" href=\""<< png2 <<"\"><img src=\""<< png2 << "\" height=" << h <<"></a></td>";
  *out << "</tr>" << endl;
  *out << "</table>" << endl;
  *out << "<br><br>" << endl;
}
//_______________________________________________________________________________________
inline double GetPrecision(int cluster_size_threshold=0, int healpix_order=0) {
//
// set $CWB_PARAMETERS_FILE precision variable
// is used in network::likelihoodWP to reduce the computational time of big clusters
// cluster_size_threshold is the threshold 
// healpix_order is the healpix_orser skymap used when cluster_size is greater 
// of cluster_size_threshold*resolution levels

  if(cluster_size_threshold>=65536) {
     cout << "GetPrecision Error : cluster_size_threshold must be < 65536" << endl;
     exit(-1);
  };
  if(healpix_order<0) {
     cout << "GetPrecision Error : healpix_order must be >= 0" << endl;
     exit(-1);
  };
  double precision = cluster_size_threshold+65536*healpix_order;

  return precision;
}

//_______________________________________________________________________________________
inline void AddRho2FAR(double rho, wavearray<double>* Xfar, wavearray<double>* Yfar, double far_rho_min, double far_drho) {
//
// Add entry rho to the ifar statistic
// is used in a compiled version for the pp report to speed up the computation

  for(int j=0; j<Xfar->size(); j++) {
    Xfar->data[j] = far_rho_min+j*far_drho;
    if(rho > Xfar->data[j]) {Yfar->data[j] += 1.;}
  }
}

//_______________________________________________________________________________________
#ifdef _USE_LAL
inline TString GetLALVersion(TString option="") {
//                                                                                         
// return LAL version (if LAL is not enabled in mdc an empty string is returned)                                                  
//                                                                                         
//                                                                                         
// Input:  option  - lalVersion,lalVersionMajor,lalVersionMinor,lalVersionMicro,
//                   lalVersionDevel,lalBuildDate,lalConfigureArgs,lalConfigureDate

  if(option=="lalVersion")       return LAL_VERSION;
  if(option=="lalVersionMajor")  {char s[32];sprintf(s,"%d",LAL_VERSION_MAJOR); return s;}
  if(option=="lalVersionMinor")  {char s[32];sprintf(s,"%d",LAL_VERSION_MINOR); return s;}
  if(option=="lalVersionMicro")  {char s[32];sprintf(s,"%d",LAL_VERSION_MICRO); return s;}
  if(option=="lalVersionDevel")  {char s[32];sprintf(s,"%d",LAL_VERSION_DEVEL); return s;}

  TString date = LAL_VCS_DATE;
  date.Remove(date.Index(' '),date.Sizeof());
  TString vcs_branch = TString(LAL_VCS_BRANCH)+" ("+date+")"; 

  return vcs_branch;
}
#else
inline TString GetLALVersion(TString options="") {return "";}
#endif

#endif
