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


// make the html simulation report : used by the cwb_report command
{

  #define NMDC_MAX 64

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  // initialize the ifo name (takes into account the non built-in detectors)
  char IFO[NIFO_MAX][NMDC_MAX];
  for(int n=0;n<nIFO;n++) {
    if(strlen(ifo[n])!=0) strcpy(IFO[n],ifo[n]); else strcpy(IFO[n],detParms[n].name);
  }

  // If fad is enable & simulation==4 & multi=true 
  // FAD multi plot for each mdc set are substituted in the sim report
  // page to the eff_freq plots
  TString pp_fad_multi = "FALSE";
  if(pp_fad!="" && simulation==4) {
    // get multi parameter from pp_fad
    pp_fad_multi = CWB::Toolbox::getParameter(pp_fad,"--multi");
    if(pp_fad_multi=="") pp_fad_multi="false";
    pp_fad_multi.ToUpper();
  }

  // read injection file types
  char    imdc_set[NMDC_MAX][128];    // injection set
  size_t  imdc_type[NMDC_MAX];        // injection type
  char    imdc_name[NMDC_MAX][128];   // injection name
  double  imdc_fcentral[NMDC_MAX];    // injection central frequencies
  double  imdc_fbandwidth[NMDC_MAX];  // injection bandwidth frequencies
  size_t  imdc_index[NMDC_MAX];       // type reference array
  size_t  imdc_iset[NMDC_MAX];        // injection set index

  int ninj=ReadInjType(mdc_inj_file,NMDC_MAX,imdc_set,imdc_type,
                       imdc_name,imdc_fcentral,imdc_fbandwidth);
  if(ninj==0) {cout << "cwb_mkhtml_sim.C : Error - no injection - terminated" << endl;exit(1);}

  TString* imdc_set_name = new TString[ninj];
  int nset=0;
  for(int i=0;i<ninj;i++) {
    bool bnew=true;
    for(int j=0;j<nset;j++) if(imdc_set[i]==imdc_set_name[j]) bnew=false;
    if(bnew) imdc_set_name[nset++]=imdc_set[i];
  }
  cout << "nset : " << nset << endl;

  // copy injection file list mdc types to report data dir

  char cmd[1024]; 
  sprintf(cmd,"cp %s %s/%s/injectionList.txt",mdc_inj_file,pp_dir,pp_data_dir); 
  cout << cmd << endl;
  gSystem->Exec(cmd);

  // create body.html file

  sprintf(netdir,"%s/%s",pp_dir,pp_data_dir);
  TString PP_DATA_PATH=netdir;
  sprintf(netdir,"%s",pp_data_dir);
  TString PP_DATA_DIR=netdir;

  ofstream out;
  char fileout[256];
  sprintf(fileout,"%s/body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout,ios::out);
  if (!out.good()) {cout << "cwb_mkhtml_sim.C : Error Opening File : " << fileout << endl;exit(1);}

  // if CWB_DOC_URL is define then man infos are added to web pages
  TString cwb_doc_url="";
  if(gSystem->Getenv("CWB_DOC_URL")!=NULL) {
    cwb_doc_url=TString(gSystem->Getenv("CWB_DOC_URL"));
  }

  out << "<hr>" << endl;
  char plot_list[1024];
  sprintf(plot_list,"<a target=\"_parent\" name=\"Full Plot List\" href=\"%s\"><h3>Full Plot List</h3></a>",pp_data_dir);
  if(cwb_doc_url!="") out<<"<table width=100%> <tr> <td align=\"left\">"<<endl;
  out << plot_list << endl;
  if(cwb_doc_url!="") {
    out<<"</td>"<<endl;
    out<<"<td align=\"right\">"<<endl;
    out<<"<a target=\"_parent\" href=\""<<cwb_doc_url.Data()
       <<"/cwb/man/Simulation-part.html#Simulation-part\">infos</a>"<<endl;
    out<<"</td> </tr> </table>"<<endl;
  }

  // read evt_parameters_ALL.txt 
  char evt_par_fname[2048];
  ifstream in_evt_par;
  sprintf(evt_par_fname,"%s/evt_parameters_ALL.txt",PP_DATA_PATH.Data());
  in_evt_par.open(evt_par_fname,ios::in);
  if(!in_evt_par.good()) 
    {cout << "cwb_mkhtml_sim.C : Error Opening File : " << evt_par_fname << endl;exit(1);}
  int epar_id;
  TString epar_mdc_name;
  float f_mean,f_rms,t_mean,t_rms;
  TString EPAR_MDC_NAME[NMDC_MAX];
  float F_MEAN[NMDC_MAX],T_RMS[NMDC_MAX];
  int ecnt=0;
  while (1) {
    in_evt_par >> epar_id >> epar_mdc_name >> f_mean >> f_rms >> t_mean >> t_rms;
    if (!in_evt_par.good()) break;
    EPAR_MDC_NAME[ecnt]=epar_mdc_name;
    F_MEAN[ecnt]=f_mean;
    T_RMS[ecnt]=t_rms;
    cout << "EVT_PAR : " << EPAR_MDC_NAME[ecnt] << " " <<  F_MEAN[ecnt] << " " << T_RMS[ecnt] << endl;
    ecnt++;
  }
  in_evt_par.close();

  char outwrite[2048];
  char file[2048];
  for (int iset=0;iset<nset;iset++) {

    cout << "MDC Set: " << iset << " " << endl;
    char th_psfix[256];
    sprintf(th_psfix,"_%d_%d",int(T_cor*100),int(T_cut*10));
    //TString threshold(th_psfix);
    TString threshold("");
    cout << "DIR: " << PP_DATA_PATH.Data() << endl;

/* ?
    char file_fit[2048];
    ofstream out_fit;
    sprintf(file_fit,"%s/fit_parameters_%s.txt",PP_DATA_PATH.Data(),threshold.Data());
    out_fit.open(file_fit,ios::out);
    int set_count=0;
*/

    cout << imdc_set_name[iset].Data() << endl;
    char file[256];
    sprintf(file,"%s/fit_parameters_%s%s.txt",PP_DATA_PATH.Data(),imdc_set_name[iset].Data(),threshold.Data());
    cout << file << endl;
    ifstream simin_hb;
    simin_hb.open(file,ios::in);
    if (!simin_hb.good()) {cout << "cwb_mkhtml_sim.C : Error Opening File : " << file << endl;exit(1);}
    int hb_count;
    TString hb_piumeno, hb_temp_wave;
    float hb_chi2, hb_err, hb_par1, hb_par2, hb_par3;
    double hb_hrss10[NMDC_MAX], hb_hrss50[NMDC_MAX], hb_hrss90[NMDC_MAX];
    TF1* gfit;
    TString hb_waveform[NMDC_MAX];
    for (int i=0; i<NMDC_MAX; i++) {
      hb_hrss10[i]=0;
      hb_hrss50[i]=0;
      hb_hrss90[i]=0;
      hb_waveform[i]="";
    }
    int maxcount=0;
    while (1) {
      simin_hb >> hb_count >> hb_chi2 >> hb_hrss50[maxcount] >> hb_piumeno 
               >> hb_err >> hb_par1 >> hb_par2 >> hb_par3 >> hb_temp_wave;
      if (!simin_hb.good()) break;
      hb_waveform[maxcount]=(TString)hb_temp_wave;
/* ?
      out_fit << " " << set_count << " " << hb_chi2 << " " << hb_hrss50[maxcount] 
              << " " << hb_piumeno << " " << hb_err << " " << hb_par1 << " " 
              << hb_par2 << " " << hb_par3 << " " << hb_temp_wave << endl;
      set_count++;
*/
      double inf = simulation==2 ? log10(factors[0]/2.) : -23;
      double sup = simulation==2 ? log10(2.*factors[nfactor-1]) : -18.5;

      if(simulation==1 && pp_factor2distance) {
        inf = log10(pp_factor2distance/factors[nfactor-1]);
        sup = log10(pp_factor2distance/factors[0]);
      }

      double par0=TMath::Log10(hb_hrss50[maxcount]);
      gfit = new TF1("logNfit",logNfit,pow(10.0,inf),pow(10.0,sup),5);
      gfit->SetNpx(100000);
      gfit->SetParameters(par0,hb_par1,hb_par2,hb_par3,pp_factor2distance);
      hb_hrss10[maxcount]=gfit->GetX(.1,pow(10.0,inf),pow(10.0,sup));
      hb_hrss90[maxcount]=gfit->GetX(.9,pow(10.0,inf),pow(10.0,sup));
      if(gfit->Eval(hb_hrss90[maxcount])<0.89) hb_hrss90[maxcount]=-1;
      maxcount++;
    }

    // table title
    out << "<font color=\"red\" style=\"font-weight:bold;\">"" <center>\t<p><h2>Detection efficiency for injections "
        << imdc_set_name[iset].Data() << "</h2>\t<p><center></font>" << endl;
    out <<"<table cellspacing=\"0\" cellpadding=\"6\" border=\"1\" align=\"center\">" << endl;

    // table row 1 (description)
    out <<"<tr align=\"center\">"<<endl;
/* ?
    sprintf(file,"%s/fit_parameters_%s.txt",PP_DATA_DIR.Data(),threshold.Data());
*/
    sprintf(file,"%s/fit_parameters_%s.txt",PP_DATA_DIR.Data(),imdc_set_name[iset].Data());
    out <<"<th colspan=\"4\"><a target=\"_parent\" href="<<file<<">fit pararameter</a></th>"<<endl;
    if(pp_fad_multi=="FALSE") {
      if(pp_factor2distance) {
        out <<"<th colspan=\""<<nIFO+4<<"\">distance vs frequency</th>"<<endl;
      } else {
        out <<"<th colspan=\""<<nIFO+4<<"\">efficiency vs frequency</th>"<<endl;
      }
    } else {
      out <<"<th colspan=\""<<nIFO+4<<"\">False Alarm Density vs rho</th>"<<endl;
    }
    out <<"</tr>"<<endl;

    // table row 2 (figures)
    out<<"<tr align=\"center\">"<<endl;
    sprintf(outwrite,
      "<td colspan=\"4\" align=\"center\"><a target=\"_parent\" href=\"%s/eff_%s%s.gif\"><img src=\"%s/eff_%s%s.gif\" height=300></a></td>",
      PP_DATA_DIR.Data(),imdc_set_name[iset].Data(),threshold.Data(),PP_DATA_DIR.Data(),imdc_set_name[iset].Data(),threshold.Data());
    out << outwrite << endl;
    if(pp_fad_multi=="FALSE") {
      sprintf(outwrite,
        "<td colspan=\"%d\" align=\"center\"><a target=\"_parent\" href=\"%s/eff_freq_%s.gif\"><img src=\"%s/eff_freq_%s.gif\" height=300></a></td>",
        nIFO+4,PP_DATA_DIR.Data(),imdc_set_name[iset].Data(),PP_DATA_DIR.Data(),imdc_set_name[iset].Data());
    } else {
      sprintf(outwrite,
        "<td colspan=\"%d\" align=\"center\"><a target=\"_parent\" href=\"%s/../fad/fad_rho_%s.gif\"><img src=\"%s/../fad/fad_rho_%s.gif\" height=300></a></td>",
        nIFO+4,PP_DATA_DIR.Data(),imdc_set_name[iset].Data(),PP_DATA_DIR.Data(),imdc_set_name[iset].Data());
    }
    out << outwrite << endl;
    out<<"</tr>"<<endl;

    // table row 3 (mdc_name/hrss10%/hrss50%/hrss90%/ifo_hrss_nre/time/freq/eff titles)
    out<<"<tr align=\"center\">"<<endl;
    out <<"<td>  waveform </td>"<<endl;
    if(simulation==1 && pp_factor2distance) {
      out <<"<td>  distance@10% </td>"<<endl;
      out <<"<td>  distance@50% </td>"<<endl;
      out <<"<td>  distance@90% </td>"<<endl;
    } else if(simulation==2) {
      out <<"<td>  snr@10% </td>"<<endl;
      out <<"<td>  snr@50% </td>"<<endl;
      out <<"<td>  snr@90% </td>"<<endl;
    } else {
      out <<"<td>  hrss@10% </td>"<<endl;
      out <<"<td>  hrss@50% </td>"<<endl;
      out <<"<td>  hrss@90% </td>"<<endl;
    }
    for (int nn=0;nn<nIFO;nn++) out << "<th>"<<IFO[nn]<<"</th>"<<endl;
    out<<"<th>time (rms)</th>"<<endl;
    out<<"<th>freq (mean)</th>"<<endl;
    out<<"<th>eff</th>"<<endl;
    out<<"</tr>"<<endl;

    // table row 4 (mdc_name/hrss10%/hrss50%/hrss90%  values)
    out<<"<tr align=\"center\">"<<endl;
 
    sprintf(outwrite,"<td>%s",hb_waveform[0].Data()); out << outwrite << endl; 
    for (int i=1; i<maxcount; i++) { 
      sprintf(outwrite,"<br>%s",hb_waveform[i].Data()); out << outwrite << endl; 
    }
    out<<"</td>"<<endl; 

    if(pp_show_eff_fit_curve) {	// write hrss@10/50/90 from efficiency fits

      sprintf(outwrite,"<td>  %.2e",hb_hrss10[0]); out << outwrite << endl; 
      for (int i=1; i<maxcount; i++) { 
        sprintf(outwrite,"<br>%.2e",hb_hrss10[i]); out << outwrite << endl; 
      }
      out<<"</td>"<<endl; 

      sprintf(outwrite,"<td>  %.2e",hb_hrss50[0]); out << outwrite << endl; 
      for (int i=1; i<maxcount; i++) { 
        sprintf(outwrite,"<br>%.2e",hb_hrss50[i]); out << outwrite << endl; 
      }
      out<<"</td>"<<endl; 

      if(hb_hrss90[0]>0) {
        sprintf(outwrite,"<td>  %.2e",hb_hrss90[0]); out << outwrite << endl; 
      } else {
        sprintf(outwrite,"<td>  NA"); out << outwrite << endl; 
      }
      for (int i=1; i<maxcount; i++) { 
        if(hb_hrss90[i]>0) {
          sprintf(outwrite,"<br>%.2e",hb_hrss90[i]); out << outwrite << endl; 
        } else {
          sprintf(outwrite,"<br>NA"); out << outwrite << endl; 
        }
      }
      out<<"</td>"<<endl; 
    
    } else {			// write blanck spaces

      for(int k=0;k<3;k++) {
        out << "<td>  " << endl; 
        for (int i=1; i<maxcount; i++) out << "<br>" << endl; 
        out<<"</td>"<<endl; 
      }
    }

    // table row 4 (ifo_hrss_nre/time/freq/eff links)
    for (int nn=0;nn<nIFO;nn++) { 
      out << "<td>" << endl;
      for (int i=0; i<maxcount; i++) {
        sprintf(outwrite,"<a target=\"_parent\" href=\"%s/hrss_%s_%s.gif\">hrss</a> / ",
                PP_DATA_DIR.Data(),IFO[nn],hb_waveform[i].Data());
        sprintf(outwrite,"%s<a target=\"_parent\" href=\"%s/nre_%s_%s.gif\">nre</a><br>",
                outwrite,PP_DATA_DIR.Data(),IFO[nn],hb_waveform[i].Data());
        out << outwrite << endl;
      }
      out << "</td>" << endl;
    }

    out << "<td align=\'right\'>" << endl;
    for (int i=0; i<maxcount; i++) {
      int ik=0;for(int k=0;k<ecnt;k++) if(EPAR_MDC_NAME[k]==hb_waveform[i]) ik=k;
      sprintf(outwrite,"<a target=\"_parent\" href=\"%s/t_%s.gif\">%3.1f  ms</a><br>",
              PP_DATA_DIR.Data(),hb_waveform[i].Data(),1000*T_RMS[ik]);
      out << outwrite << endl;
    }
    out << "</td>" << endl;

    out << "<td align=\'right\'>" << endl;
    for (int i=0; i<maxcount; i++) {
      int ik=0;for(int k=0;k<ecnt;k++) if(EPAR_MDC_NAME[k]==hb_waveform[i]) ik=k;
      sprintf(outwrite,"<a target=\"_parent\" href=\"%s/f_%s.gif\">%3.2f  Hz</a><br>",
              PP_DATA_DIR.Data(),hb_waveform[i].Data(),F_MEAN[ik]);
      out << outwrite << endl;
    }
    out << "</td>" << endl;

    out << "<td>" << endl;
    for (int i=0; i<maxcount; i++) {
      sprintf(outwrite,"<a target=\"_parent\" href=\"%s/%s.gif\">plot</a> / ",
              PP_DATA_DIR.Data(),hb_waveform[i].Data());
      out << outwrite << endl;

      sprintf(outwrite,"<a target=\"_parent\" href=\"%s/eff_%s.txt\">txt</a><br>",
              PP_DATA_DIR.Data(),hb_waveform[i].Data());
      out << outwrite << endl;
    }
    out << "</td>" << endl;

    out << "</tr>" << endl;
  
    out << "</table>" << endl;
/* ?
    out_fit.close();
*/
 
    out << "<p>" << endl;
    out << endl;
  }
  out.close();
  exit(0);
}

