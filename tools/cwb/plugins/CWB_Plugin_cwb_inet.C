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


#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "gwavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TMath.h"
#include "mdc.hh"
#include <vector>

void CWB_Plugin_RemoveTemporaryFiles(TFile* jfile, CWB::config* cfg);
TString CWB_Plugin_CheckIFO(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                            network* net, TString ifo);
TString CWB_Plugin_CheckTYPE(TString cwb_inet_options,TFile* jfile, CWB::config* cfg, TString type);
void CWB_Plugin_frdisplay(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                          network* net, WSeries<double>* x, TString ifo, int type);
void CWB_Plugin_psd(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                          network* net, WSeries<double>* x, TString ifo, int type);
void CWB_Plugin_inj(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                          network* net, WSeries<double>* x, TString ifo, int type);
void CWB_Plugin_nRMS(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                          network* net, WSeries<double>* x, TString ifo, int type);
void CWB_Plugin_wdm(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                          network* net, WSeries<double>* x, TString ifo, int type);
void CWB_Plugin_emax(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                          network* net, WSeries<double>* xdummy, TString sdummy, int type);
void CWB_Plugin_sparse(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                          network* net, TString ifo, int type);
void CWB_Plugin_ced(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                          network* net, TString ifo, int type);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// Plugin to show data with frdisplay & psd with the cwb_inet command

  if(TString(cfg->analysis)!="2G") {
    cout << "CWB_Plugin_cwb_inet.C : allowed only for 2G analysis" << endl;
    gSystem->Exit(1);
  }

  // get cwb_inet options
  TString cwb_inet_options=TString(gSystem->Getenv("CWB_INET_OPTIONS"));

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_cwb_inet.C : ifo : " << ifo.Data() << " - type : " << type << endl;
  cout << "       options : " << cwb_inet_options << endl;
  cout << endl;

  // get ifo : return if not equal to input plugin parameter (ifo)
  TString cwb_inet_tool = CWB::Toolbox::getParameter(cwb_inet_options,"--tool");
  cwb_inet_tool.ToUpper();
  TString cwb_inet_type = CWB::Toolbox::getParameter(cwb_inet_options,"--type");
  cwb_inet_type.ToUpper();

  // exit if cwb_inet_tool can not be applied 
  bool error = true;
  // import istage,jstage
  //CWB_STAGE gISTAGE=-1; IMPORT(CWB_STAGE,gISTAGE)
  int gISTAGE=-1; IMPORT(int,gISTAGE)
  int gJSTAGE=-1; IMPORT(int,gJSTAGE)

  if((gISTAGE<CWB_STAGE_STRAIN)&&(gJSTAGE>=CWB_STAGE_STRAIN)||(gJSTAGE==CWB_STAGE_FULL)) {
    if((cwb_inet_tool=="PSD")&&((cwb_inet_type=="STRAIN")||(cwb_inet_type=="MDC"))) error=false;
    if((cwb_inet_tool=="WDM")&&((cwb_inet_type=="STRAIN")||(cwb_inet_type=="MDC"))) error=false;
    if((cwb_inet_tool=="FRDISPLAY")&&((cwb_inet_type=="STRAIN")||(cwb_inet_type=="MDC"))) error=false;
    if(cwb_inet_tool=="INJ") error=false;
  }
  if((gISTAGE<CWB_STAGE_CSTRAIN)&&(gJSTAGE>=CWB_STAGE_CSTRAIN)||(gJSTAGE==CWB_STAGE_FULL)) {
    if((cwb_inet_tool=="WDM")&&(cwb_inet_type=="STRAIN")) error=false;
    if((cwb_inet_tool=="PSD")&&(cwb_inet_type=="WHITE")) error=false;
    if((cwb_inet_tool=="WDM")&&(cwb_inet_type=="WHITE")) error=false;
    if((cwb_inet_tool=="FRDISPLAY")&&(cwb_inet_type=="WHITE")) error=false;
    if(cwb_inet_tool=="INJ") error=false;
    if(cwb_inet_tool=="NRMS") error=false;
  }
  if((gISTAGE<CWB_STAGE_COHERENCE)&&(gJSTAGE>=CWB_STAGE_COHERENCE)||(gJSTAGE==CWB_STAGE_FULL)) {
    if(cwb_inet_tool=="EMAX") error=false;
  }
  if((gISTAGE<CWB_STAGE_SUPERCLUSTER)&&(gJSTAGE>=CWB_STAGE_SUPERCLUSTER)||(gJSTAGE==CWB_STAGE_FULL)) {
    if((cwb_inet_tool=="SPARSE")&&(cwb_inet_type=="SUPERCLUSTER")) error=false;
  }
  if((gISTAGE<CWB_STAGE_LIKELIHOOD)&&(gJSTAGE>=CWB_STAGE_LIKELIHOOD)||(gJSTAGE==CWB_STAGE_FULL)) {
    if((cwb_inet_tool=="SPARSE")&&(cwb_inet_type=="LIKELIHOOD")) error=false;
    if(cwb_inet_tool=="CED") error=false;
  }

  if(error) {
    cout << "CWB_Plugin_cwb_inet.C : cwb_inet option can not applied for this stage range , process terminated!!!" << endl;
    cout << endl;
    cout << "OPTIONS     : " << cwb_inet_options << endl;
    cout << "BEGIN STAGE : " << cwb::GetStageString(CWB_STAGE(gISTAGE)) << endl;
    cout << "END STAGE   : " << cwb::GetStageString(CWB_STAGE(gJSTAGE)) << endl;
    cout << endl;
    CWB_Plugin_RemoveTemporaryFiles(NULL, cfg);
    gSystem->Exit(1);
  }

  // execute plugin
       if(cwb_inet_tool=="FRDISPLAY") CWB_Plugin_frdisplay(cwb_inet_options, jfile, cfg, net, x, ifo, type);
  else if(cwb_inet_tool=="PSD")       CWB_Plugin_psd(cwb_inet_options, jfile, cfg, net, x, ifo, type);
  else if(cwb_inet_tool=="INJ")       CWB_Plugin_inj(cwb_inet_options, jfile, cfg, net, x, ifo, type);
  else if(cwb_inet_tool=="NRMS")      CWB_Plugin_nRMS(cwb_inet_options, jfile, cfg, net, x, ifo, type);
  else if(cwb_inet_tool=="WDM")       CWB_Plugin_wdm(cwb_inet_options, jfile, cfg, net, x, ifo, type);
  else if(cwb_inet_tool=="EMAX")      CWB_Plugin_emax(cwb_inet_options, jfile, cfg, net, x, ifo, type);
  else if(cwb_inet_tool=="SPARSE")    CWB_Plugin_sparse(cwb_inet_options, jfile, cfg, net, ifo, type);
  else if(cwb_inet_tool=="CED")       CWB_Plugin_ced(cwb_inet_options, jfile, cfg, net, ifo, type);
  else {
    cout << "CWB_Plugin_cwb_inet.C : Bad inet option --tool " << cwb_inet_tool << endl;
    cout << "Select : --tool FRDISPLAY/PSD/NRMS/WDM/EMAX/SPARSE/CED" << endl << endl;
    CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
    gSystem->Exit(1);
  }
}

void 
CWB_Plugin_ced(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                  network* net, TString ifo, int type) {

  if(type==CWB_PLUGIN_ILIKELIHOOD) { 	// save ced to jfile
    cfg->cedDump = true;
    unsigned int jobfOptions = cfg->jobfOptions;
    jobfOptions |= CWB_JOBF_SAVE_CED;
    cfg->jobfOptions = (CWB_JOBF_OPTIONS)jobfOptions;
    return;
  }	
  if(type!=CWB_PLUGIN_CLOSE_JOB) return;

  TString cwb_inet_draw = CWB::Toolbox::getParameter(cwb_inet_options,"--draw");
  cwb_inet_draw.ToUpper();

  // get job file name
  TString jname = jfile->GetPath();
  jname.ReplaceAll(":/","");
  cout << jname.Data() << endl;

  if(cwb_inet_draw=="TRUE") { 
    // creates a temporary file name
    char fmacroName[1024];
    unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
    sprintf(fmacroName,"%s/CWB_Plugin_cwb_inet_ced_macro_%s_%d_job%d.XXXXXX",
            cfg->tmp_dir,cfg->data_label,Pid,net->nRun);
    CWB::Toolbox::mksTemp(fmacroName);  // create temporary file, must be deleted only at the end of run
    ofstream out;
    out.open(fmacroName,ios::out);
    if(!out.good()) {
      cout << "CWB_Plugin_cwb_inet.C - Error : Opening File : " << fmacroName << endl;
      gSystem->Exit(1);
    }
    char line[1024]; 
    sprintf(line,"{TFile* ifile = new TFile(\"%s\");TBrowser br(\"CWB\",ifile);}",jname.Data()); 
    out << line << endl; 
    out.close();

    char cmd[1024]; 
    // start root browser
    sprintf(cmd,"root -l %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
    // remove temporary files
    sprintf(cmd,"rm %s",jname.Data());
    cout << cmd << endl; gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
  }

  CWB_Plugin_RemoveTemporaryFiles(NULL, cfg);
  gSystem->Exit(0);
}

void 
CWB_Plugin_sparse(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                  network* net, TString ifo, int type) {

  // get type : return if not equal to input plugin parameter (type)
  TString cwb_inet_type = CWB_Plugin_CheckTYPE(cwb_inet_options, jfile, cfg, type);
  if(cwb_inet_type=="") return;
  cwb_inet_type.ToLower();

  bool check=false;
  if((type==CWB_PLUGIN_ISUPERCLUSTER)&&(cwb_inet_type=="supercluster")) check=true;
  if((type==CWB_PLUGIN_ILIKELIHOOD)&&(cwb_inet_type=="likelihood")) check=true;
  if(!check) return;

  // get ifo : return if not equal to input plugin parameter (ifo)
  TString cwb_inet_ifo = CWB::Toolbox::getParameter(cwb_inet_options,"--ifo");
  // check if ifo is in the network
  if(cwb_inet_ifo!="") CWB_Plugin_CheckIFO(cwb_inet_options, jfile, cfg, net, ifo);

  watplot* WTS = new watplot(const_cast<char*>("WTS"));

  TString cwb_inet_draw = CWB::Toolbox::getParameter(cwb_inet_options,"--draw");
  cwb_inet_draw.ToUpper();
  TString odir = cfg->dump_dir; 

  // get draw mode 
  TString cwb_inet_mode = CWB::Toolbox::getParameter(cwb_inet_options,"--mode");
  int mode = (cwb_inet_mode!="") ? cwb_inet_mode.Atoi() : 1;
  TString gmode="";
  if(mode==0) mode=1;
  if(mode==1) gmode="(E00+E90)/2";
  if(mode==2) gmode="sqrt((E00+E90)/2)";
  if(mode==3) gmode="amplitude:00";
  if(mode==4) gmode="energy:00";
  if(mode==5) gmode="|amplitude:00|";
  if(mode==6) gmode="amplitude:90";
  if(mode==7) gmode="energy:90";
  if(mode==8) gmode="|amplitude:90|";

  // get core pixels selection
  // core=true  -> rebuild TF map using only core pixels
  // core=false -> rebuild TF map using all pixels
  TString cwb_inet_core = CWB::Toolbox::getParameter(cwb_inet_options,"--core");
  cwb_inet_core = ((cwb_inet_core=="true")||(cwb_inet_core=="false")) ? cwb_inet_core : "true";	

  // open ioutput root file
  char fName[1024];
  sprintf(fName,"%s/sparse_%s_%s_%d_%s_job%lu.root",
          odir.Data(),cwb_inet_type.Data(),cwb_inet_ifo.Data(),mode,cfg->data_label,net->nRun);
  TFile* ofile = new TFile(fName,"RECREATE");

  int nIFO = net->ifoListSize();
  detector* pD[NIFO_MAX];                       //! pointers to detectors
  for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);
  //for(int n=0;n<nIFO;n++) cout << n << " detector Name : " << pD[n]->Name << endl;

  bool singleDetector = false;
  if(nIFO==2 && TString(pD[1]->Name)==pD[0]->Name) singleDetector=true;

  for(int n=0;n<nIFO;n++) {    // loop over detectors
    if(cwb_inet_ifo!="") if(net->getifo(n)->Name!=cwb_inet_ifo) continue;

    // create ifo root directory to store sparse TF
    TDirectory* dsparse = ofile->mkdir(pD[n]->Name);
    dsparse->cd();

    int nSS = pD[n]->vSS.size();
    for(int i=0;i<nSS;i++) {    // loop over levels

      int level;
      if(cwb_inet_type=="supercluster") level = i+cfg->l_low;
      if(cwb_inet_type=="likelihood")   level = i+cfg->l_low;

      // create ifo root directory to store sparse TF
      char dname[32]; sprintf(dname,"%s:level-%d",pD[n]->Name,level);
      TDirectory* dlevel = dsparse->mkdir(dname);
      dlevel->cd();

      size_t rateANA=cfg->inRate>>cfg->levelR;
      int rate  = rateANA>>level;
      double dt = 1000./rate;
      double df = rateANA/2./double(1<<level);
      cout << "level : " << level << "\t rate : " << rate 
           << "\t df(hz) : " << df << "\t dt(ms) : " << dt << endl;

      bool core = (cwb_inet_core=="true") ? true : false;
      pD[n]->vSS[i].Expand(core);    // rebuild wseries from sparse table
      //cout << "REBUILD SIZE " << pD[n]->vSS[i].size() << endl;

      WSeries<double>* WS = (WSeries<double>*)(&(pD[n]->vSS[i]));

      // Plot WDM Scalogram
      double start = WS->start();
      double stop  = WS->start()+WS->size()/WS->rate();
      double flow  = cfg->fLow;
      double fhigh = cfg->fHigh;

      WTS->plot(WS, mode, start, stop,const_cast<char*>("COLZ"));
      WTS->hist2D->GetYaxis()->SetRangeUser(flow, fhigh);
      char title[64]; 
      sprintf(title,"SPARSE TF (%s) - ifo:%s - type:sparse - level:%d - dt:%g (ms) - df=%g (Hz)",
              gmode.Data(),net->getifo(n)->Name,level,dt,df);
      WTS->hist2D->SetTitle(title);
      char name[64]; sprintf(name,"%s:tf:sparse:%d",net->getifo(n)->Name,level);
      WTS->canvas->Write(name);

      // project hist2D to time axis
      sprintf(name,"%s:t:sparse:%d",net->getifo(n)->Name,level);
      sprintf(title,"SPARSE:TIME (%s) - ifo:%s - type:sparse - level:%d - dt:%g (ms) - df=%g (Hz)",
              gmode.Data(),net->getifo(n)->Name,level,dt,df);
      TH1D *tProjH2 = WTS->hist2D->ProjectionX(name);
      WTS->canvas->Clear();
      tProjH2->SetTitle(title);
      tProjH2->SetStats(kFALSE);
      tProjH2->Draw();
      WTS->canvas->Write(name);

      // project hist2D to frequency axis
      int hist2D_size = WTS->hist2D->GetNbinsX();
      int edgeSize = (int)(hist2D_size*cfg->segEdge/(2*cfg->segEdge+cfg->segLen));
      sprintf(name,"%s:f:sparse:%d",net->getifo(n)->Name,level);
      sprintf(title,"SPARSE:FREQ (%s) - ifo:%s - type:sparse - level:%d - dt:%g (ms) - df=%g (Hz)",
              gmode.Data(),net->getifo(n)->Name,level,dt,df);
      TH1D *fProjH2 = WTS->hist2D->ProjectionY(name,edgeSize,hist2D_size-edgeSize);
      WTS->canvas->Clear();
      fProjH2->SetTitle(title);
      fProjH2->SetStats(kFALSE);
      fProjH2->Draw();
      WTS->canvas->Write(name);

    }
    if(singleDetector) break;
  }

  delete WTS;
  delete ofile;

  if(cwb_inet_draw=="TRUE") { 
    // creates a temporary file name
    char fmacroName[1024];
    unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
    sprintf(fmacroName,"%s/CWB_Plugin_cwb_inet_sparse_macro_%s_%d_job%d.XXXXXX",
            cfg->tmp_dir,cfg->data_label,Pid,net->nRun);
    CWB::Toolbox::mksTemp(fmacroName);  // create temporary file, must be deleted only at the end of run
    ofstream out;
    out.open(fmacroName,ios::out);
    if(!out.good()) {cout << "CWB_Plugin_cwb_inet.C - Error : Opening File : " << fmacroName << endl;gSystem->Exit(1);}
    char line[1024]; 
    sprintf(line,"{TFile* ifile = new TFile(\"%s\");TBrowser br(\"CWB\",ifile);}",fName); 
    out << line << endl; 
    out.close();

    char cmd[1024]; 
    // start root browser
    sprintf(cmd,"root -l %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
    // remove temporary files
    sprintf(cmd,"rm %s",fName);
    cout << cmd << endl; gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
  } else {
    cout << endl << "Dump Sparse : " << fName << endl << endl;
  }

  CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
  gSystem->Exit(0);
}

void 
CWB_Plugin_nRMS(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                network* net, WSeries<double>* x, TString ifo, int type) {

  if(type!=CWB_PLUGIN_ODATA_CONDITIONING) return;

  // get type : return if not equal to input plugin parameter (type)
  TString cwb_inet_type = CWB_Plugin_CheckTYPE(cwb_inet_options, jfile, cfg, type);
  if(cwb_inet_type=="") return;

  // get ifo : return if not equal to input plugin parameter (ifo)
  TString cwb_inet_ifo = CWB_Plugin_CheckIFO(cwb_inet_options, jfile, cfg, net, ifo);
  if(cwb_inet_ifo=="") return;

  if(x==NULL) {
    cout << "CWB_Plugin_cwb_inet.C : Error : NULL data " << endl;
    gSystem->Exit(1);
  }

  TString cwb_inet_gps = CWB::Toolbox::getParameter(cwb_inet_options,"--gps");
  double gps = 0;
  if(cwb_inet_gps!="") {
    if(!cwb_inet_gps.IsFloat()) {
      cout<< "CWB_Plugin_cwb_inet.C : Error : --gps value is not a number" << endl;
      exit(1);
    }
    if(cwb_inet_gps.Atof()>0) gps=cwb_inet_gps.Atof();
    cout.precision(14);
    cout << "CWB_Plugin_cwb_inet.C : gps = " << gps << endl;
  }

  TString cwb_inet_draw = CWB::Toolbox::getParameter(cwb_inet_options,"--draw");
  cwb_inet_draw.ToUpper();
  TString odir = cfg->dump_dir; 

  // get oneside option	// default is TRUE
  TString cwb_inet_oneside = CWB::Toolbox::getParameter(cwb_inet_options,"--oneside");
  cwb_inet_oneside.ToUpper();
  if(cwb_inet_oneside=="") cwb_inet_oneside="TRUE";
  if((cwb_inet_oneside!="TRUE")&&(cwb_inet_oneside!="FALSE")) {
    cout << "CWB_Plugin_cwb_inet.C : Error : wrong --oneside option , must be true/false" << endl;
    gSystem->Exit(1);
  }
  TString psdtype = (cwb_inet_oneside=="TRUE") ? "oneside" : "doubleside";
  bool oneside = (cwb_inet_oneside=="TRUE") ? true : false;

  // get nRMS
  WSeries<double> nRMS;
  int nIFO = net->ifoListSize();
  for(int n=0;n<nIFO;n++) if(ifo==net->getifo(n)->Name) nRMS=net->getifo(n)->nRMS;

  int levels = nRMS.getLevel()+1;		// number of levels
  int slices = nRMS.size()/levels;		// number of nRMS samples 
  double length = slices*cfg->whiteStride;	// nRMS len in sec

  // nRMS do not come from a TF transform
  // nRMS params must be fixed before to pass to watplot
  double rate = nRMS.rate();
  double fNinquist = rate/2.;
  nRMS.rate(nRMS.size()/length);
  WDM<double>* wdm = (WDM<double>*) nRMS.pWavelet;
  wdm->nSTS=nRMS.size();

  // Plot nRMS Scalogram
  double start = nRMS.start();
  double stop  = nRMS.start()+nRMS.size()/nRMS.rate();
  double flow  = cfg->fLow;
  double fhigh = cfg->fHigh;

  if(gps>0) { 	// if options --gps value > 0 we extract rms at that time and produce PSD

    if(gps<start || gps>stop) {
      cout.precision(14);
      cout << "CWB_Plugin_cwb_inet.C : Error : --gps value is out of the data range ["  
           << start << ":" << stop << "]" << endl;
      gSystem->Exit(1);
    }

    int M = nRMS.getLevel();
    const double scale = 1./nRMS.wrate();
    int igps = (gps-nRMS.start())/scale + 1;	// time index

    // dump nRMS to file
    char file[1024];
    double chuncklen=8;	// default = 8 sec
    // get chuncklen options - Example --chuncklen 8
    TString cwb_inet_chuncklen = CWB::Toolbox::getParameter(cwb_inet_options,"--chuncklen");
    if(cwb_inet_chuncklen.IsFloat()) chuncklen = cwb_inet_chuncklen.Atof();
    if(chuncklen==0) chuncklen=x->size()/x->rate()-2*cfg->segEdge; // chuncklen = full buffer
    sprintf(file,"%s/nrms_%s_gps%10.2f_%s_%s_%d_%s_job%lu.txt",
            odir.Data(),psdtype.Data(),gps,cwb_inet_type.Data(),ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    cout << endl << "Dump nRMS " << "(chuncklen = " << chuncklen << ")"
         << " : " << file << endl << endl;
    ofstream out;
    out.open(file,ios::out);
    if(!out.good()) {cout << "CWB_Plugin_cwb_inet.C : Error Opening File : " << file << endl;exit(1);}

    double inRate = cfg->fResample>0 ? cfg->fResample : cfg->inRate;
    double rescale = oneside ? sqrt(2.)/sqrt(inRate) : 1./sqrt(inRate);

    double df=fNinquist/M;
    for(int j=1; j<M; ++j) {
      double f = j*df+df/2.;
      if(f<flow || f>fhigh) continue;
      out << j*df+df/2. << " " << nRMS.data[igps*(M+1)+j]*rescale << endl; 
    }
    out.close();

    if(cwb_inet_draw=="TRUE") { 
      TString cwb_inet_save = CWB::Toolbox::getParameter(cwb_inet_options,"--save");
      cwb_inet_save.ToUpper();
      int save  = cwb_inet_save=="TRUE" ? 1 : 0;
      TString cwb_inet_range = CWB::Toolbox::getParameter(cwb_inet_options,"--range");
      cwb_inet_range.ToUpper();
      int range = cwb_inet_range=="FIX" ? 1 : 0;
      char cmd[1024];
      // execute cwb_draw_sensitivity
      sprintf(cmd,"export CWB_SENSITIVITY_FILE_NAME=\"%s\";",file);
      sprintf(cmd,"%s export CWB_SENSITIVITY_SAVE_PLOT=%d;",cmd,save);
      sprintf(cmd,"%s export CWB_SENSITIVITY_RANGE_FIX=%d;",cmd,range);
      sprintf(cmd,"%s root -n -l ${CWB_MACROS}/cwb_draw_sensitivity.C",cmd);
      cout << cmd << endl; gSystem->Exec(cmd);
  
      // remove temporary files
      sprintf(cmd,"rm %s",file);
      cout << cmd << endl; gSystem->Exec(cmd);
    }

  } else { 	// produce the T/F plot of nRMS

    // open ioutput root file
    char fName[1024];
    sprintf(fName,"%s/nRMS_%s_%d_%s_job%lu.root",
            odir.Data(),ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
    TFile* ofile = new TFile(fName,"RECREATE");
    cout << endl << "Dump nRMS : " << fName << endl << endl;

    watplot* WTS = new watplot(const_cast<char*>("nRMS"));
    WTS->plot(nRMS, 3, start, stop,const_cast<char*>("COLZ"));
    WTS->hist2D->GetYaxis()->Set(WTS->hist2D->GetNbinsY(),0, fNinquist); // correct freq range
    WTS->hist2D->GetYaxis()->SetRangeUser(flow, fhigh);
    char title[64]; 
    sprintf(title,"nRMS TF - ifo:%s", ifo.Data());
    WTS->hist2D->SetTitle(title);
    char name[64]; sprintf(name,"%s",ifo.Data());
    WTS->canvas->SetLogz(true);
    WTS->canvas->Write(name);
//    delete WTS;		// this instruction is commented because it produces a crash (why?)

    delete ofile;

    if(cwb_inet_draw=="TRUE") { 
      // creates a temporary file name
      char fmacroName[1024];
      unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
      sprintf(fmacroName,"%s/CWB_Plugin_cwb_inet_nRMS_macro_%s_%d_job%d.XXXXXX",
              cfg->tmp_dir,cfg->data_label,Pid,net->nRun);
      CWB::Toolbox::mksTemp(fmacroName);  // create temporary file, must be deleted only at the end of run
      ofstream out;
      out.open(fmacroName,ios::out);
      if(!out.good()) {cout << "CWB_Plugin_cwb_inet.C - Error : Opening File : " << fmacroName << endl;gSystem->Exit(1);}
      char line[1024]; 
      sprintf(line,"{TFile* ifile = new TFile(\"%s\");TBrowser br(\"CWB\",ifile);}",fName); 
      out << line << endl; 
      out.close();
  
      char cmd[1024]; 
      // start root browser
      sprintf(cmd,"root -l %s",fmacroName);
      cout << cmd << endl; gSystem->Exec(cmd);
      // remove temporary files
      sprintf(cmd,"rm %s",fName);
      cout << cmd << endl; gSystem->Exec(cmd);
      sprintf(cmd,"rm %s",fmacroName);
      cout << cmd << endl; gSystem->Exec(cmd);
    }
  }

  CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
  gSystem->Exit(0);
}

void 
CWB_Plugin_wdm(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
               network* net, WSeries<double>* x, TString ifo, int type) {

  // get type : return if not equal to input plugin parameter (type)
  TString cwb_inet_type = CWB_Plugin_CheckTYPE(cwb_inet_options, jfile, cfg, type);
  if(cwb_inet_type=="") return;
  cwb_inet_type.ToLower();

  // get ifo : return if not equal to input plugin parameter (ifo)
  TString cwb_inet_ifo = CWB_Plugin_CheckIFO(cwb_inet_options, jfile, cfg, net, ifo);
  if(cwb_inet_ifo=="") return;

  if(x==NULL) {
    cout << "CWB_Plugin_cwb_inet.C : Error : NULL data " << endl;
    gSystem->Exit(1);
  }

  watplot* WTS = new watplot(const_cast<char*>("WTS"));

  wavearray<double> y = *x;                   // copy original data to y (contains the final handled data)
  WSeries<double> w;                          // temporary array for data manipulation
  Meyer<double> B(1024);           	      // set wavelet for resampling
  WDM<double>* pwdm = NULL;

  if(type==CWB_PLUGIN_STRAIN) {		      // we must apply resampling
    if(cfg->fResample>0) {                    // RESAMPLING
      y.FFTW(1);
      y.resize(cfg->fResample/y.rate()*y.size());
      y.FFTW(-1);
      y.rate(cfg->fResample);
    }
    w.Forward(y,B,cfg->levelR);
    w.getLayer(y,0);
  }

  TString cwb_inet_draw = CWB::Toolbox::getParameter(cwb_inet_options,"--draw");
  cwb_inet_draw.ToUpper();
  TString odir = cfg->dump_dir; 

  // get draw mode 
  TString cwb_inet_mode = CWB::Toolbox::getParameter(cwb_inet_options,"--mode");
  int mode = (cwb_inet_mode!="") ? cwb_inet_mode.Atoi() : 2;
  TString gmode="";
  if(mode==0) mode=2;
  if(mode==1) gmode="(E00+E90)/2";
  if(mode==2) gmode="sqrt((E00+E90)/2)";
  if(mode==3) gmode="amplitude:00";
  if(mode==4) gmode="energy:00";
  if(mode==5) gmode="|amplitude:00|";
  if(mode==6) gmode="amplitude:90";
  if(mode==7) gmode="energy:90";
  if(mode==8) gmode="|amplitude:90|";

  // open ioutput root file
  char fName[1024];
  sprintf(fName,"%s/wdm_%s_%s_%d_%s_job%lu.root",
          odir.Data(),cwb_inet_type.Data(),ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
  TFile* ofile = new TFile(fName,"RECREATE");
  cout << endl << "Dump WDM : " << fName << endl << endl;

  for(int level=cfg->l_high; level>=cfg->l_low; level--) {    // loop over levels

    // create ifo root directory to store wdm TF
    char dname[32]; sprintf(dname,"%s:level-%d",ifo.Data(),level);
    TDirectory* dwdm = ofile->mkdir(dname);
    dwdm->cd();

    int layers = level>0 ? 1<<level : 0;      // get layers
    cout << "level : " << level << " layers : " << layers << endl;

    size_t rateANA=cfg->inRate>>cfg->levelR;
    int rate  = rateANA>>level;
    double dt = 1000./rate;
    double df = rateANA/2./double(1<<level);
    cout << "level : " << level << "\t rate : " << rate 
         << "\t df(hz) : " << df << "\t dt(ms) : " << dt << endl;

    pwdm = net->getwdm(layers+1);             // get pointer to wdm transform
    if(pwdm==NULL) {
      cout << "CWB_Plugin_cwb_inet.C : Error - WDM not defined !!!" << endl;
      gSystem->Exit(1);
    }

    w.Forward(y,*pwdm);                       // apply wdm transformation
    WSeries<double>* WS = (WSeries<double>*)(&w);

    // Plot WDM Scalogram
    double start = WS->start();
    double stop  = WS->start()+WS->size()/WS->rate();
    double flow  = cfg->fLow;
    double fhigh = cfg->fHigh;

    WTS->plot(WS, mode, start, stop,const_cast<char*>("COLZ"));
    WTS->hist2D->GetYaxis()->SetRangeUser(flow, fhigh);
    char title[64]; 
    sprintf(title,"WDM:TF (%s) - ifo:%s - type:%s - level:%d - dt:%g (ms) - df=%g (Hz)",
            gmode.Data(),ifo.Data(),cwb_inet_type.Data(),level,dt,df);
    WTS->hist2D->SetTitle(title);
    char name[64]; sprintf(name,"%s:tf:%s-%d",ifo.Data(),cwb_inet_type.Data(),level);
    WTS->canvas->Write(name);

    // Project hist2D to time axis
    sprintf(name,"%s:t:%s-%d",ifo.Data(),cwb_inet_type.Data(),level);
    sprintf(title,"WDM:TIME (%s) - ifo:%s - type:%s - level:%d - dt:%g (ms) - df=%g (Hz)",
            gmode.Data(),ifo.Data(),cwb_inet_type.Data(),level,dt,df);
    TH1D *tProjH2 = WTS->hist2D->ProjectionX(name);
    WTS->canvas->Clear();
    tProjH2->SetTitle(title);
    tProjH2->SetStats(kFALSE);
    tProjH2->Draw();
    WTS->canvas->Write(name);

    // Project hist2D to frequency axis
    int hist2D_size = WTS->hist2D->GetNbinsX();
    int edgeSize = (int)(hist2D_size*cfg->segEdge/(2*cfg->segEdge+cfg->segLen));
    sprintf(name,"%s:f:%s-%d",ifo.Data(),cwb_inet_type.Data(),level);
    sprintf(title,"WDM:FREQ (%s) - ifo:%s - type:%s - level:%d - dt:%g (ms) - df=%g (Hz)",
            gmode.Data(),ifo.Data(),cwb_inet_type.Data(),level,dt,df);
    TH1D *fProjH2 = WTS->hist2D->ProjectionY(name,edgeSize,hist2D_size-edgeSize);
    WTS->canvas->Clear();
    fProjH2->SetTitle(title);
    fProjH2->SetStats(kFALSE);
    fProjH2->Draw();
    WTS->canvas->Write(name);

    w.Inverse();                              // inverse transform
    y = w;                                    // copy manipulated data to y
  }

  delete ofile;
//  delete WTS;		// this instruction is commented because it produces a crash (why?)

  if(cwb_inet_draw=="TRUE") { 
    // creates a temporary file name
    char fmacroName[1024];
    unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
    sprintf(fmacroName,"%s/CWB_Plugin_cwb_inet_wdm_macro_%s_%d_job%d.XXXXXX",
            cfg->tmp_dir,cfg->data_label,Pid,net->nRun);
    CWB::Toolbox::mksTemp(fmacroName);  // create temporary file, must be deleted only at the end of run
    ofstream out;
    out.open(fmacroName,ios::out);
    if(!out.good()) {cout << "CWB_Plugin_cwb_inet.C - Error : Opening File : " << fmacroName << endl;gSystem->Exit(1);}
    char line[1024]; 
    sprintf(line,"{TFile* ifile = new TFile(\"%s\");TBrowser br(\"CWB\",ifile);}",fName); 
    out << line << endl; 
    out.close();

    char cmd[1024]; 
    // start root browser
    sprintf(cmd,"root -l %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
    // remove temporary files
    sprintf(cmd,"rm %s",fName);
    cout << cmd << endl; gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
  }

  CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
  gSystem->Exit(0);
}

void 
CWB_Plugin_emax(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                network* net, WSeries<double>* xdummy, TString sdummy, int type) {

  if(type!=CWB_PLUGIN_XCOHERENCE) return;

  // import factor
  int gIFACTOR=-1; IMPORT(int,gIFACTOR)
  // import resolution
  size_t gILEVEL=-1; IMPORT(size_t,gILEVEL)

  // get level : if level>0 then select only one level
  TString cwb_inet_level = CWB::Toolbox::getParameter(cwb_inet_options,"--level");
  int ilevel = cwb_inet_level.Atoi();
  if(ilevel>0 && ilevel!=gILEVEL) return;

  TString cwb_inet_draw = CWB::Toolbox::getParameter(cwb_inet_options,"--draw");
  cwb_inet_draw.ToUpper();
  TString odir = cfg->dump_dir; 

  int level = gILEVEL;
  int layers = level>0 ? 1<<level : 0;      // get layers

  size_t rateANA=cfg->inRate>>cfg->levelR;
  int rate  = rateANA>>level;
  double dt = 1000./rate;
  double df = rateANA/2./double(1<<level);
  cout << "level : " << level << "\t rate(hz) : " << rate << "\t layers : " << layers
       << "\t df(hz) : " << df << "\t dt(ms) : " << dt << endl;

  int nIFO = net->ifoListSize();
  detector* pD[NIFO_MAX];                            // pointers to detectors
  for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);

  // open ioutput root file
  char fName[1024];
  sprintf(fName,"%s/emax_%d_%s_job%lu.root",
          odir.Data(),int(pD[0]->getTFmap()->start()),cfg->data_label,net->nRun);
  TFile* ofile = (ilevel||level==cfg->l_high) ? new TFile(fName,"RECREATE") : new TFile(fName,"UPDATE");
  if(ofile==NULL||!ofile->IsOpen())
    {cout << "CWB_Plugin_cwb_inet.C - Error opening : " << fName <<  endl;exit(1);}
  cout << endl << "Dump EMAX : " << fName << endl << endl;

  watplot WTS(const_cast<char*>("WTS"));

  bool singleDetector = false;
  if(nIFO==2 && TString(pD[1]->Name)==pD[0]->Name) singleDetector=true;

  WSeries<double> WSE = *pD[0]->getTFmap();	     
  for(int n=0; n<=nIFO; n++) {                       // produce TF maps with max over the sky energy

    if(n<nIFO) WSE.add(*pD[n]->getTFmap());	     // produce TF maps with sum of energy over the ifo

    TString ifo = n==nIFO ? "NET" : net->ifoName[n];

    // create ifo root directory to store emax TF
    TDirectory* demax = ofile->mkdir(ifo.Data());
    demax->cd();

    WSeries<double>* WS = n==nIFO ? &WSE : pD[n]->getTFmap();

    //scalogram maps
    double start = WS->start()+cfg->segEdge;
    double stop  = WS->start()+WS->size()/WS->rate()-cfg->segEdge;
    double flow  = cfg->fLow;
    double fhigh = cfg->fHigh;
    WTS.plot(*WS, 0, start, stop,const_cast<char*>("COLZ"));
    WTS.hist2D->GetYaxis()->SetRangeUser(flow, fhigh);

    // Dump EMAX Scalogram
    char title[64]; 
    sprintf(title,"EMAX:TF (energy) - ifo:%s - factor:%d - level:%d - dt:%g (ms) - df=%g (Hz)",
            ifo.Data(),gIFACTOR,level,dt,df);
    WTS.hist2D->SetTitle(title);
    char name[64]; sprintf(name,"%s:tf:%d-%d",ifo.Data(),gIFACTOR,level);
    WTS.canvas->Write(name);

    // Project hist2D to time axis
    sprintf(name,"%s:t:%d-%d",ifo.Data(),gIFACTOR,level);
    sprintf(title,"EMAX:TIME (energy) - ifo:%s - factor:%d - level:%d - dt:%g (ms) - df=%g (Hz)",
            ifo.Data(),gIFACTOR,level,dt,df);
    TH1D *tProjH2 = WTS.hist2D->ProjectionX(name);
    WTS.canvas->Clear();
    tProjH2->SetTitle(title);
    tProjH2->SetStats(kFALSE);
    tProjH2->Draw();
    WTS.canvas->Write(name);

    // Project hist2D to frequency axis
    int hist2D_size = WTS.hist2D->GetNbinsX();
    int edgeSize = (int)(hist2D_size*cfg->segEdge/(2*cfg->segEdge+cfg->segLen));
    sprintf(name,"%s:f:%d-%d",ifo.Data(),gIFACTOR,level);
    sprintf(title,"EMAX:FREQ (energy) - ifo:%s - factor:%d - level:%d - dt:%g (ms) - df=%g (Hz)",
            ifo.Data(),gIFACTOR,level,dt,df);
    TH1D *fProjH2 = WTS.hist2D->ProjectionY(name,edgeSize,hist2D_size-edgeSize);
    WTS.canvas->Clear();
    fProjH2->SetTitle(title);
    fProjH2->SetStats(kFALSE);
    fProjH2->Draw();
    WTS.canvas->Write(name);

    if(singleDetector) break;
  }

  delete ofile;

  if(!ilevel && level!=cfg->l_low) return;

  if(cwb_inet_draw=="TRUE") { 
    // creates a temporary file name
    char fmacroName[1024];
    unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
    sprintf(fmacroName,"%s/CWB_Plugin_cwb_inet_wdm_macro_%s_%d_job%d.XXXXXX",
            cfg->tmp_dir,cfg->data_label,Pid,net->nRun);
    CWB::Toolbox::mksTemp(fmacroName);  // create temporary file, must be deleted only at the end of run
    ofstream out;
    out.open(fmacroName,ios::out);
    if(!out.good()) {cout << "CWB_Plugin_cwb_inet.C - Error : Opening File : " << fmacroName << endl;gSystem->Exit(1);}
    char line[1024]; 
    sprintf(line,"{TFile* ifile = new TFile(\"%s\");TBrowser br(\"CWB\",ifile);}",fName); 
    out << line << endl; 
    out.close();

    char cmd[1024]; 
    // start root browser
    sprintf(cmd,"root -l %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
    // remove temporary files
    sprintf(cmd,"rm %s",fName);
    cout << cmd << endl; gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
  }

  CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
  gSystem->Exit(0);
}

void 
CWB_Plugin_psd(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
               network* net, WSeries<double>* x, TString ifo, int type) {

  // get type : return if not equal to input plugin parameter (type)
  TString cwb_inet_type = CWB_Plugin_CheckTYPE(cwb_inet_options, jfile, cfg, type);
  if(cwb_inet_type=="") return;

  // get ifo : return if not equal to input plugin parameter (ifo)
  TString cwb_inet_ifo = CWB_Plugin_CheckIFO(cwb_inet_options, jfile, cfg, net, ifo);
  if(cwb_inet_ifo=="") return;

  if(x==NULL) {
    cout << "CWB_Plugin_cwb_inet.C : Error : NULL data " << endl;
    gSystem->Exit(1);
  }

  TString cwb_inet_draw = CWB::Toolbox::getParameter(cwb_inet_options,"--draw");
  cwb_inet_draw.ToUpper();
  TString odir = cfg->dump_dir; 
  //TString odir = (cwb_inet_draw=="TRUE") ? cfg->tmp_dir : cfg->dump_dir; 

  double chuncklen=8;	// default = 8 sec
  // get chuncklen options - Example --chuncklen 8
  TString cwb_inet_chuncklen = CWB::Toolbox::getParameter(cwb_inet_options,"--chuncklen");
  if(cwb_inet_chuncklen.IsFloat()) chuncklen = cwb_inet_chuncklen.Atof();
  if(chuncklen==0) chuncklen=x->size()/x->rate()-2*cfg->segEdge; // chuncklen = full buffer

  // get oneside option	// default is TRUE
  TString cwb_inet_oneside = CWB::Toolbox::getParameter(cwb_inet_options,"--oneside");
  cwb_inet_oneside.ToUpper();
  if(cwb_inet_oneside=="") cwb_inet_oneside="TRUE";
  if((cwb_inet_oneside!="TRUE")&&(cwb_inet_oneside!="FALSE")) {
    cout << "CWB_Plugin_cwb_inet.C : Error : wrong --oneside option , must be true/false" << endl;
    gSystem->Exit(1);
  }
  TString psdtype = (cwb_inet_oneside=="TRUE") ? "oneside" : "doubleside";
  bool oneside = (cwb_inet_oneside=="TRUE") ? true : false;

  char file[1024];
  sprintf(file,"%s/psd_%s_%s_%s_%d_%s_job%lu.txt",
          odir.Data(),psdtype.Data(),cwb_inet_type.Data(),ifo.Data(),int(x->start()),cfg->data_label,net->nRun);
  cout << endl << "Dump PSD " << "(chuncklen = " << chuncklen << ")"
       << " : " << file << endl << endl;
  CWB::Toolbox::makeSpectrum(file, *x, chuncklen, cfg->segEdge, oneside); // dump spectrum

  if(cwb_inet_draw=="TRUE") { 
    TString cwb_inet_save = CWB::Toolbox::getParameter(cwb_inet_options,"--save");
    cwb_inet_save.ToUpper();
    int save  = cwb_inet_save=="TRUE" ? 1 : 0;
    TString cwb_inet_range = CWB::Toolbox::getParameter(cwb_inet_options,"--range");
    cwb_inet_range.ToUpper();
    int range = cwb_inet_range=="FIX" ? 1 : 0;
    char cmd[1024];
    // execute cwb_draw_sensitivity
    sprintf(cmd,"export CWB_SENSITIVITY_FILE_NAME=\"%s\";",file);
    sprintf(cmd,"%s export CWB_SENSITIVITY_SAVE_PLOT=%d;",cmd,save);
    sprintf(cmd,"%s export CWB_SENSITIVITY_RANGE_FIX=%d;",cmd,range);
    sprintf(cmd,"%s root -n -l ${CWB_MACROS}/cwb_draw_sensitivity.C",cmd);
    cout << cmd << endl; gSystem->Exec(cmd);

    // remove temporary files
    sprintf(cmd,"rm %s",file);
    cout << cmd << endl; gSystem->Exec(cmd);
  }

  CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
  gSystem->Exit(0);
}

void 
CWB_Plugin_inj(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
               network* net, WSeries<double>* x, TString ifo, int type) {

  if(type!=CWB_PLUGIN_OREADDATA) return;

  TString cwb_inet_draw = CWB::Toolbox::getParameter(cwb_inet_options,"--draw");
  cwb_inet_draw.ToUpper();
  TString odir = cfg->dump_dir; 

  int nIFO = net->ifoListSize();
  detector* pD[NIFO_MAX];                            // pointers to detectors
  for(int n=0;n<nIFO;n++) pD[n] = net->getifo(n);

  // open ioutput root file
  char fName[1024];
  sprintf(fName,"%s/inj_%d_%s_job%lu.root",
          odir.Data(),int(pD[0]->getTFmap()->start()),cfg->data_label,net->nRun);
  TFile* ofile = new TFile(fName,"RECREATE");
  if(ofile==NULL||!ofile->IsOpen())
    {cout << "CWB_Plugin_cwb_inet.C - Error opening : " << fName <<  endl;exit(1);}
  cout << endl << "Dump INJ : " << fName << endl << endl;

  for(int n=0;n<nIFO;n++) { 

    // create ifo root directory to store injections
    TDirectory* dinj = ofile->mkdir(pD[n]->Name);
    dinj->cd();

    wavearray<double> w = pD[n]->HoT;
    std::vector<double>* pT = net->getmdcTime();
    double dT = cfg->iwindow/2.;
    double offset = cfg->segEdge;

    int j,nstop,nstrt;
    size_t k;
    double F,T;
    size_t K    = pT->size();
    size_t N    = w.size();
    double rate = w.rate();                        // simulation rate

    dT = fabs(dT);

    // isolate injections
    for(k=0; k<K; k++) {

      T = (*pT)[k] - w.start();
 
      nstrt = int((T - dT)*rate);
      nstop = int((T + dT)*rate);
      if(nstrt<=0) nstrt = 0;
      if(nstop>=int(N)) nstop = N;
      if(nstop<=0) continue;                     // outside of the segment
      if(nstrt>=int(N)) continue;                // outside of the segment

      gwavearray<double> gx(nstop-nstrt);
      gx.rate(rate);
      gx.start(0);

      for(j=nstrt; j<nstop; j++) gx.data[j-nstrt]=w.data[j];
      gx.SetComment(net->mdcList[k]);
      char name[256];sprintf(name,"%s:%2.2f",pD[n]->Name,nstrt/rate); 
      gx.Write(name);
    }
  }
  delete ofile;

  if(cwb_inet_draw=="TRUE") { 
    // creates a temporary file name
    char fmacroName[1024];
    unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
    sprintf(fmacroName,"%s/CWB_Plugin_cwb_inet_wdm_macro_%s_%d_job%d.XXXXXX",
            cfg->tmp_dir,cfg->data_label,Pid,net->nRun);
    CWB::Toolbox::mksTemp(fmacroName);  // create temporary file, must be deleted only at the end of run
    ofstream out;
    out.open(fmacroName,ios::out);
    if(!out.good()) {cout << "CWB_Plugin_cwb_inet.C - Error : Opening File : " << fmacroName << endl;gSystem->Exit(1);}
    char line[1024]; 
    sprintf(line,"{TFile* ifile = new TFile(\"%s\");TBrowser br(\"CWB\",ifile);}",fName); 
    out << line << endl; 
    out.close();

    char cmd[1024]; 
    // start root browser
    sprintf(cmd,"root -l %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
    // remove temporary files
    sprintf(cmd,"rm %s",fName);
    cout << cmd << endl; gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",fmacroName);
    cout << cmd << endl; gSystem->Exec(cmd);
  }

  CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
  gSystem->Exit(0);
}

void 
CWB_Plugin_frdisplay(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, 
                     network* net, WSeries<double>* x, TString ifo, int type) {

  // get type : return if not equal to input plugin parameter (type)
  TString cwb_inet_type = CWB_Plugin_CheckTYPE(cwb_inet_options, jfile, cfg, type);
  if(cwb_inet_type=="") return;

  // get ifo : return if not equal to input plugin parameter (ifo)
  TString cwb_inet_ifo = CWB_Plugin_CheckIFO(cwb_inet_options, jfile, cfg, net, ifo);
  if(cwb_inet_ifo=="") return;

  if(x==NULL) {
    cout << "CWB_Plugin_cwb_inet.C : Error : NULL data " << endl;
    gSystem->Exit(1);
  }

  // get ifo index
  int ifoID =0; for(int n=0;n<cfg->nIFO;n++) if(ifo==net->getifo(n)->Name) {ifoID=n;break;}

  // get frdisplay directory
  TString home_frdisplay="";
  if(gSystem->Getenv("HOME_FRDISPLAY")==NULL) {
    cout << "Error : environment HOME_FRDISPLAY is not defined!!!" << endl;
    gSystem->Exit(1);
  } else {
    home_frdisplay=TString(gSystem->Getenv("HOME_FRDISPLAY"));
  }
  cout << "home_frdisplay : " << home_frdisplay.Data() << endl;

  // get data buffer start,stop   
  int job_start = x->start();
  int job_stop  = x->stop(); 
  cout << "job_start : " << job_start << " job_stop : " << job_stop << endl;
                                                                          
  // prepare frame list format for frdisplay                                
  gRandom->SetSeed(0);                                                      
  int rnID = gRandom->Uniform(0,10000000);   // random name ID     

  // write temporary frame file
  char chName[64];  sprintf(chName,"%s",cfg->channelNamesRaw[ifoID]);
  char frLabel[64]; sprintf(frLabel,"CWB_Plugin_cwb_inet_%s_%s_%d",ifo.Data(),cwb_inet_type.Data(),rnID);
  TString frFile = CWB::Toolbox::WriteFrameFile(*x, chName, "FRDISPLAY", frLabel, cfg->tmp_dir);

  // write temporary ffl file
  char ffl[1024];
  sprintf(ffl,"%s/CWB_Plugin_cwb_inet_%s_%s_%d.ffl",cfg->tmp_dir,ifo.Data(),cwb_inet_type.Data(),rnID);
  cout << "ffl : " << ffl << endl;                                         

  ofstream out;
  out.open(ffl,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << ffl << endl;gSystem->Exit(1);}
  out << frFile.Data() << " " << job_start << " " << job_stop-job_start << " " << 0 << " " << 0 << endl;
  cout << frFile.Data() << endl;
  out.close();

  // build cmd frdisplay command
  char cmd[1024];
  sprintf(cmd,"%s/FrDisplay -d 5 -proc -t %s -i %s",home_frdisplay.Data(),chName,ffl);

  // get filter options - Example --filter 50,
  TString cwb_inet_hpf = CWB::Toolbox::getParameter(cwb_inet_options,"--hpf");
  if(cwb_inet_hpf!="") sprintf(cmd,"%s -k \"-Bu -Hp -o 6 -a %s\"",cmd,cwb_inet_hpf.Data());

  char xoptions[256]="";
  // get decimateby options - Example --decimateby 16
  TString cwb_inet_decimateby = CWB::Toolbox::getParameter(cwb_inet_options,"--decimateby");
  if(cwb_inet_decimateby!="") sprintf(xoptions,"%s -decimateby %s",xoptions,cwb_inet_decimateby.Data());
  // get -downmix options - Example --downmix 400
  TString cwb_inet_downmix = CWB::Toolbox::getParameter(cwb_inet_options,"--downmix");
  if(cwb_inet_downmix!="") sprintf(xoptions,"%s -downmix %s",xoptions,cwb_inet_downmix.Data());
  // get -uscaleby options - Example --uscaleby 0.5
  TString cwb_inet_uscaleby = CWB::Toolbox::getParameter(cwb_inet_options,"--uscaleby");
  if(cwb_inet_uscaleby!="") sprintf(xoptions,"%s -uscaleby %s",xoptions,cwb_inet_uscaleby.Data());
  // get -decimategain options - Example --decimategain 48
  TString cwb_inet_decimategain = CWB::Toolbox::getParameter(cwb_inet_options,"--decimategain");
  if(cwb_inet_decimategain!="") sprintf(xoptions,"%s -decimategain %s",xoptions,cwb_inet_decimategain.Data());

  // execute frdisplay command
  if(TString(xoptions)!="") sprintf(cmd,"%s -x \"%s\"",cmd,xoptions);
  cout << cmd << endl; gSystem->Exec(cmd);

  // remove temporary files
  sprintf(cmd,"rm %s %s",ffl,frFile.Data());
  cout << cmd << endl; gSystem->Exec(cmd);

  CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
  gSystem->Exit(0);
}

TString 
CWB_Plugin_CheckTYPE(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, TString type)  {

  // get type : return if not equal to input plugin parameter (type)
  TString cwb_inet_type = CWB::Toolbox::getParameter(cwb_inet_options,"--type");
  cwb_inet_type.ToUpper();
  if((cwb_inet_type=="STRAIN")||(cwb_inet_type=="MDC")||(cwb_inet_type=="WHITE")||
     (cwb_inet_type=="SUPERCLUSTER")||(cwb_inet_type=="LIKELIHOOD")) {
    bool selected=false;
    if((type==CWB_PLUGIN_DATA)         &&(cwb_inet_type=="STRAIN"))       selected=true;
    if((type==CWB_PLUGIN_MDC)          &&(cwb_inet_type=="MDC"))          selected=true;
    if((type==CWB_PLUGIN_WHITE)        &&(cwb_inet_type=="WHITE"))        selected=true;
    if((type==CWB_PLUGIN_IDATA_CONDITIONING)&&(cwb_inet_type=="STRAIN"))  selected=true;
    if((type==CWB_PLUGIN_ODATA_CONDITIONING)&&(cwb_inet_type=="STRAIN"))  selected=true;
    if((type==CWB_PLUGIN_ISUPERCLUSTER)&&(cwb_inet_type=="SUPERCLUSTER")) selected=true;
    if((type==CWB_PLUGIN_ILIKELIHOOD)  &&(cwb_inet_type=="LIKELIHOOD"))   selected=true;
    if(!selected) return "";
  } else {
    cout << "CWB_Plugin_cwb_inet.C : Option error --type " << cwb_inet_type << endl << endl;
    cout << "Select : --type strain/mdc/white/supercluster/likelihood" << endl << endl;
    CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
    gSystem->Exit(1);
  }

  return cwb_inet_type;
}

TString 
CWB_Plugin_CheckIFO(TString cwb_inet_options, TFile* jfile, CWB::config* cfg, network* net, TString ifo)  {

  // get ifo : return if not equal to input plugin parameter (ifo)
  TString cwb_inet_ifo = CWB::Toolbox::getParameter(cwb_inet_options,"--ifo");
  cwb_inet_ifo.ToUpper();
  // check if ifo is in the network
  bool isPresentIFO=false;
  for(int n=0;n<cfg->nIFO;n++) if(cwb_inet_ifo==net->getifo(n)->Name) isPresentIFO=true;
  if(isPresentIFO) {
    if(ifo!="0") if(cwb_inet_ifo!=ifo) return "";
  } else {
    cout << "CWB_Plugin_cwb_inet.C : Option error --ifo " << cwb_inet_ifo << endl;
    char ifo_opts[32];strcpy(ifo_opts,net->getifo(0)->Name);
    for(int n=1;n<cfg->nIFO;n++) sprintf(ifo_opts,"%s/%s",ifo_opts,net->getifo(n)->Name);
    cout << "Select : --ifo " << ifo_opts << endl << endl;
    CWB_Plugin_RemoveTemporaryFiles(jfile, cfg);
    gSystem->Exit(1);
  }

  return cwb_inet_ifo;
}

void 
CWB_Plugin_RemoveTemporaryFiles(TFile* jfile, CWB::config* cfg)  {

  // remove temporary job file 
  if(jfile) {
    TString jname = jfile->GetPath();
    bool remove = true;
    if(jname.Contains("init_")) remove=false;
    if(jname.Contains("strain_")) remove=false;
    if(jname.Contains("cstrain_")) remove=false;
    if(jname.Contains("coherence_")) remove=false;
    if(jname.Contains("supercluster_")) remove=false;
    if(jname.Contains("likelihood_")) remove=false;
    if(remove) {
      jname.ReplaceAll(":/","");
      //cout << jname.Data() << endl;
      gSystem->Exec(TString("rm "+jname).Data());
    }
  }

  // remove temporary plugin file
  gSystem->Exec((TString("rm ")+cfg->plugin.GetTitle()).Data());

  // remove temporary plugin config file
  TString configPluginName = cfg->configPlugin.GetTitle();
  if(configPluginName!="") 
    gSystem->Exec((TString("rm ")+configPluginName).Data());

  // search output root file in the system list
  // & remove output root file 
  TList *files = (TList*)gROOT->GetListOfFiles();
  if(files) {
    TIter next(files);
    TSystemFile *file;
    while ((file=(TSystemFile*)next())) {
      TString fname = file->GetName();
      if(fname.Contains("wave_")) gSystem->Exec((TString("rm ")+fname).Data());
    }
  }
}
