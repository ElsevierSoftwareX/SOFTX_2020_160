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


// make the html production report : used by the cwb_report command
{

  #define NEPS 28 
  #define NEPS2 19 

  #include "wat.hh"

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  bool singleDetector=false; 
  if(nIFO==1) { // Single Detector Mode
    CWB::config config;
    config.Import();
    config.SetSingleDetectorMode();
    config.Export();
    singleDetector=true; 
  }
  if(nIFO==2) { // 2 detectors with the same name "same detector -> nIFO=1" 
    if(TString(ifo[0])==TString(ifo[1])) nIFO=1;   
  }

  // if CWB_DOC_URL is define then man infos are added to web pages
  TString cwb_doc_url="";
  if(gSystem->Getenv("CWB_DOC_URL")!=NULL) {
    cwb_doc_url=TString(gSystem->Getenv("CWB_DOC_URL"));
  }

  // check vetoes
  bool bcat2  = false;
  bool bhveto = false;
  bool bcat3  = false;

  for(int i=0;i<nVDQF;i++) if(VDQF[i].cat==CWB_CAT2)  bcat2=true;
  for(int i=0;i<nVDQF;i++) if(VDQF[i].cat==CWB_HVETO) bhveto=true;
  for(int i=0;i<nVDQF;i++) if(VDQF[i].cat==CWB_CAT3)  bcat3=true;

  TString NetCC,Rho,RhOut,RhAcor,svED,sPEN,sHRSS,sfLow,sfHigh,sfResample,spp_rho_min,spp_rho_max;

  char s[16];
  sprintf(s,"%.2f",T_cor);      NetCC=s;
  sprintf(s,"%.2f",T_cut);      Rho=s;
  sprintf(s,"%.2f",T_out);      RhOut=s;
  sprintf(s,"%.2f",T_acor);     RhAcor=s;
  sprintf(s,"%.2f",T_vED);      svED=s;
  sprintf(s,"%.2f",T_pen);      sPEN=s;
  sprintf(s,"%.2f",fLow);       sfLow=s;
  sprintf(s,"%.2f",fHigh);      sfHigh=s;
  sprintf(s,"%.2f",pp_rho_min); spp_rho_min=s;
  sprintf(s,"%.2f",pp_rho_max); spp_rho_max=s;

  TString psfix(data_label);
  TFile* psp = TFile::Open(net_file_name);
  cout<<"Opening "<<net_file_name<<endl;
  TTree* wave_tree = (TTree*)psp->Get("waveburst");

  // create body.html file

  ofstream out;
  char fileout[1024];
  sprintf(fileout,"%s/body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << fileout << endl;exit(1);}

  TString PP_DATA_PATH, PP_DATA_DIR;
  sprintf(netdir,"%s/%s",pp_dir,pp_data_dir);
  PP_DATA_PATH=netdir;
  sprintf(netdir,"%s",pp_data_dir);
  PP_DATA_DIR=netdir;

  TString RATEBACKGROUND;

  // convert eps plots to gif plots
  vector<TString> epsList = CWB::Toolbox::getFileListFromDir(PP_DATA_PATH, ".eps");
  for(int i=0;i<epsList.size();i++) {
    TString ofile=epsList[i];
    ofile.ReplaceAll(".eps",".gif");
    char cmd[2048];
    sprintf(cmd,"convert %s %s",epsList[i].Data(),ofile.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",epsList[i].Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  // get livetime   
  char fileliv[1024];
  sprintf(fileliv,"%s/live.txt",PP_DATA_PATH.Data());
  cout << fileliv << endl;
  int countlag=0;
  double OLIVETIME = GetLiveTime(fileliv,0,0,countlag);	  // get zero  livetime
  double LIVETIME  = GetLiveTime(fileliv,-1,-1,countlag); // get total livetime
  LIVETIME-=OLIVETIME;	// non zero live time
  countlag-=1;		// subtract the zero lag 
  cout.precision(14);
  cout << OLIVETIME << " " << LIVETIME << endl;

  ifstream bkgin;
  bkgin.open(html_body_prod_template,ios::in);
  if (!bkgin.good()) {cout << "Error Opening File : " << html_body_prod_template << endl;exit(1);}

  char exec[1024];
  char istring[1024];

  while (1) {
    bkgin.getline(istring,256);
    if (!bkgin.good()) break;
    TString ostringa(istring);
    if(bcat3||bhveto) {
      ostringa.ReplaceAll("<!--VETO","");
      ostringa.ReplaceAll("VETO-->","");
    } else {
      ostringa.ReplaceAll("<!--NOVETO","");
      ostringa.ReplaceAll("VETONO-->","");
    }
    if(TString(analysis)=="1G") {
      ostringa.ReplaceAll("<!--RHO_PF","");
      ostringa.ReplaceAll("RHO_PF-->","");
    } else {
      ostringa.ReplaceAll("<!--RHO_SUBNET","");
      ostringa.ReplaceAll("RHO_SUBNET-->","");
    }
    if(pp_rho_vs_rate_no_multiplicity) {
      ostringa.ReplaceAll("<!--MULTI","");
      ostringa.ReplaceAll("MULTI-->","");
    }
    ostringa.ReplaceAll("FLOW",sfLow);
    ostringa.ReplaceAll("FHIGH",sfHigh);
    ostringa.ReplaceAll("PP_RHO_MIN",spp_rho_min);
    ostringa.ReplaceAll("PP_RHO_MAX",spp_rho_max);
    ostringa.ReplaceAll("PP_DATA_DIR",PP_DATA_DIR);
    ostringa.ReplaceAll("NETCC",NetCC);
    ostringa.ReplaceAll("RHO",Rho);
    ostringa.ReplaceAll("RH_OUT",RhOut);
    ostringa.ReplaceAll("RH_ACOR",RhAcor);
    ostringa.ReplaceAll("svED",svED);
    TString Rate= "rate = " + RATEBACKGROUND; 
    ostringa.ReplaceAll("RATEBACKGROUND",Rate);
    ostringa.ReplaceAll("SUBTITLE",subtitle);
    ostringa.ReplaceAll("TITLE",title);
    char ltime[1024];
    sprintf(ltime, "zero lag : <font color=\"red\"> %.2f sec = %.2f days</font>",
            OLIVETIME,OLIVETIME/86400.);
    ostringa.ReplaceAll("LIVETIME",ltime);
    sprintf(ltime, "non-zero lags : <font color=\"red\"> %i lags - %.2f sec = %.2f days = %.1f years </fonts>",
            countlag,LIVETIME,LIVETIME/(float)86400.,LIVETIME/86400./365.);
    ostringa.ReplaceAll("LIVETIM2",ltime);
    ostringa.ReplaceAll("PLOT_LIST",pp_data_dir);
    if(cwb_doc_url!="") {
      ostringa.ReplaceAll("<!--CWB_DOC_URL","");
      ostringa.ReplaceAll("CWB_DOC_URL-->","");
      ostringa.ReplaceAll("XCWB_DOC_URL",cwb_doc_url.Data());
    }
	
    ostringa.ReplaceAll("</html>","");
    out << ostringa.Data() << endl;
  }
  bkgin.close();

  char _ifo[NIFO_MAX][8];
  for(int n=0;n<nIFO;n++) {
    if(strlen(ifo[n])>0) strcpy(_ifo[n],ifo[n]);           // built in detector
    else                 strcpy(_ifo[n],detParms[n].name); // user define detector
  }

  char file_ev[1024];
  sprintf(file_ev,"%s/events.txt",PP_DATA_PATH.Data());
  ifstream for_ev(file_ev);
  char file_header[1024];
  sprintf(file_header,"%s/events_header.txt",PP_DATA_PATH.Data());
  ofstream for_header(file_header);
  char file_notsorted[1024];
  sprintf(file_notsorted,"%s/events_notsorted.txt",PP_DATA_PATH.Data());
  ofstream for_notsorted(file_notsorted);
  char in_ev[100000];
  while(1) {
    for_ev.getline(in_ev,100000);
    if (!for_ev.good()) break;
    TString sin_ev(in_ev);
    if (sin_ev.Contains("#")) for_header << sin_ev.Data() << endl;
    else for_notsorted << sin_ev.Data() << endl;
  }
  for_header.close();
  for_notsorted.close();
  char file_sorted[1024];
  sprintf(file_sorted,"%s/events_sorted.txt",PP_DATA_PATH.Data());
  sprintf(exec,"sort -g -r --key=3 %s > %s",file_notsorted, file_sorted); 
  gSystem->Exec(exec);
  sprintf(exec,"cat %s %s > %s/EVENTS.txt",file_header, file_sorted, PP_DATA_PATH.Data());
  gSystem->Exec(exec);

  sprintf(file_notsorted,"%s/events_notsorted.txt",PP_DATA_PATH.Data());
  ifstream f_ev(file_sorted);
  char pm[2];
  char c3[2];
  float icc, irho, iacor, ilag, islag, ilik, ipen, icHH, ivHH, ivED, icc2, icc3; //SLAG
  int ifreq, iband;
  float idur;
  int isize, irate, irun;
  float SNR0, SNR1, SNR2, SNR3, hrss0, hrss1, hrss2, hrss3;
  double time0, time1, time2, time3;
  float phi, theta, psi;
  int ecount=0;

  char os[1024];

  out << "<table border=0 cellpadding=2 class=\"datagrid\">" << endl;
  out << "<tr align=\"center\">"<< endl;
  out << "<td>ID</td>"<< endl;
  if(bhveto) out << "<td>KW</td>"<< endl;
  if(bcat3)  out << "<td>cat3  </td>"<< endl;
  out << "<td>rho["<<pp_irho<<"]</td>"<< endl;
  out << "<td>cc["<<pp_inetcc<<"]</td>"<< endl;
  if(TString(analysis)=="2G") out << "<td>subnet</td>"<< endl;
//  out << "<td>ecor</td>"<< endl;
  out << "<td>lag</td>"<< endl;
  out << "<td>slag</td>"<< endl;
  out << "<td>SNRnet</td>"<< endl;
  if(TString(analysis)=="1G") out << "<td>pf</td>"<< endl;
  if(TString(analysis)=="1G") out << "<td>vED</td>"<< endl;
  out << "<td>freq</td>"<< endl;
  out << "<td>bw</td>"<< endl;
  out << "<td>dur</td>"<< endl;
  out << "<td>size</td>"<< endl;
  if(optim==true) out << "<td>res</td>"<< endl;
  out << "<td>run</td>"<< endl;
  for (int nn=0;nn<nIFO;nn++) out << "<td>GPS " << _ifo[nn] << "</td>"<< endl;
  for (int nn=0;nn<nIFO;nn++) out << "<td>SNR " << _ifo[nn] << "</td>"<< endl;
/*
  for (int nn=0;nn<nIFO;nn++) out << "<td>hrss " << _ifo[nn] << "</td>"<< endl;
  out << "<td>phi</td>"<< endl;
  out << "<td>theta</td>"<< endl;
  out << "<td>psi</td>"<< endl;
*/
  out << "</tr>"<< endl;

  double start[5];
  while(ecount<pp_max_nloudest_list) {
    if (nIFO>3) {
       f_ev>>pm>>c3>>irho>>icc>>icc2>>icc3>>iacor>>ilag>>islag>>ilik>>ipen>>icHH>>ifreq>>iband>>idur>>isize>>irate>>irun>>time0>>time1>>time2>>time3>>SNR0>>SNR1>>SNR2>>SNR3>>hrss0>>hrss1>>hrss2>>hrss3>>phi>>theta>>psi;
    } else if (nIFO>2) { 
       f_ev>>pm>>c3>>irho>>icc>>icc2>>icc3>>iacor>>ilag>>islag>>ilik>>ipen>>icHH>>ifreq>>iband>>idur>>isize>>irate>>irun>>time0>>time1>>time2>>SNR0>>SNR1>>SNR2>>hrss0>>hrss1>>hrss2>>phi>>theta>>psi;
    } else {
       f_ev>>pm>>c3>>irho>>icc>>icc2>>icc3>>iacor>>ilag>>islag>>ilik>>ipen>>icHH>>ifreq>>iband>>idur>>isize>>irate>>irun>>time0>>time1>>SNR0>>SNR1>>hrss0>>hrss1>>phi>>theta>>psi;
    }
    if (!f_ev.good()) break;
    TString ppm(pm);
    TString cc3(c3);
    if (ppm.Contains("+")) {
      ecount++;
      cout << ecount << " ";cout.flush();
      out << "<tr align=\"center\">"<< endl;
      sprintf(os,"run==%i",irun);
      TString s_ced=GetFileLabel(wave_tree,irun,ilag,islag,segEdge,psfix);
      char namedir[1024];
      sprintf(namedir,"");
      if(singleDetector&&strlen(channelNamesRaw[0])>1) {
        start[0]=GetStart(wave_tree, 0, irun, irho, time0, time1, analysis, pp_irho);
        sprintf(namedir, "%s_%s_%.3f", _ifo[0],channelNamesRaw[0],start[0]);
      } else {
        for (int nn=0; nn<nIFO; nn++) sprintf(namedir,"%s%s",namedir,_ifo[nn]);
        for (int nn=0; nn<nIFO; nn++) {
          start[nn]=GetStart(wave_tree, nn, irun, irho, time0, time1, analysis, pp_irho);
          sprintf(namedir,"%s_%.3f",namedir,start[nn]);
        }
      }
      sprintf(os,"<td><a href=\"ced/ced_%s/%s\" target=\"_blank\">%i</a></td>",s_ced.Data(),namedir,ecount);
      out << os << endl;
      if(bhveto) {
        sprintf(os,"<td>%s</td>",ppm.Remove(0,1).Data());	
        out << os << endl;
      }
      if(bcat3) {
        sprintf(os,"<td>%s</td>",cc3.Remove(0,1).Data());	
        out << os << endl;
      }
      sprintf(os,"<td>%.2f</td>",irho);
      out << os << endl;
      sprintf(os,"<td>%.2f</td>",icc);
      out << os << endl;
      if(TString(analysis)=="2G") {
        sprintf(os,"<td>%.2f</td>",icc2);
        out << os << endl;
      }
/*
      sprintf(os,"<td>%.2f</td>",iacor);
      out << os << endl;
*/
      sprintf(os,"<td>%i</td>",(int)ilag);
      out << os << endl;
      sprintf(os,"<td>%i</td>",(int)islag);
      out << os << endl;
      sprintf(os,"<td>%3.1f</td>",sqrt(ilik));
      out << os << endl;
      sprintf(os,"<td>%.3f</td>",ipen);
      if(TString(analysis)=="1G") out << os << endl;
      sprintf(os,"<td>%.3f</td>",icHH);
      if(TString(analysis)=="1G") out << os << endl;
      sprintf(os,"<td><font color=\"black\">%i</font></td>",ifreq);
      out << os << endl;
      sprintf(os,"<td>%i</td>",iband);
      out << os << endl;
      sprintf(os,"<td>%.3f</td>",idur);
      out << os << endl;
      sprintf(os,"<td>%i</td>",isize);
      out << os << endl;
      sprintf(os,"<td>%i</td>",irate);
      if(optim==true) out << os << endl;
      sprintf(os,"<td>%i</td>",irun);
      out << os << endl;
      sprintf(os,"<td>%.2f</td>",time0);
      out << os << endl;
      if (nIFO>1){
        sprintf(os,"<td>%.2f</td>",time1);
        out << os << endl;
      }
      if (nIFO>2){
        sprintf(os,"<td>%.2f</td>",time2);
        out << os << endl;
      }
      if (nIFO>3){
        sprintf(os,"<td>%.2f</td>",time3);
        out << os << endl;
      }
      sprintf(os,"<td>%3.1f</td>",sqrt(SNR0));
      out << os << endl;
      if (nIFO>1){
        sprintf(os,"<td>%3.1f</td>",sqrt(SNR1));
        out << os << endl;
      }
      if (nIFO>2){
        sprintf(os,"<td>%3.1f</td>",sqrt(SNR2));
        out << os << endl;
      }
      if (nIFO>3){
        sprintf(os,"<td>%3.1f</td>",sqrt(SNR3));
        out << os << endl;
      }
/*
      sprintf(os,"<td>%.1e</td>",hrss0);
      out << os << endl;
      sprintf(os,"<td>%.1e</td>",hrss1);
      out << os << endl;
      if (nIFO>2){
        sprintf(os,"<td>%.1e</td>",hrss2);
        out << os << endl;
      }
      if (nIFO>3){
        sprintf(os,"<td>%.1e</td>",hrss3);
        out << os << endl;
      }
      sprintf(os,"<td>%.2f</td>",phi);
      out << os << endl;
      sprintf(os,"<td>%.2f</td>",theta);
      out << os << endl;
      sprintf(os,"<td>%.2f</td>",psi);
      out << os << endl;
*/
      out << "</tr>" << endl;
    }

  }
  cout << endl;
  f_ev.close();
  out << "</table>" << endl;
  out << "<p>" << endl;
  out << endl;

  cout << exec << endl;
  gSystem->Exec(exec);

  sprintf(file_notsorted,"%s/events_notsorted.txt",PP_DATA_PATH.Data());
  sprintf(exec,"rm %s",file_notsorted);
  cout << exec << endl;
  gSystem->Exec(exec);
  sprintf(exec,"rm %s",file_header);
  cout << exec << endl;
  gSystem->Exec(exec);

  out << "</table>" << endl;
  out << "<p>" << endl;
  out << endl;
  out.close();
  exit(0);
}

