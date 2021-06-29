/*
# Copyright (C) 2019 Francesco Salemi
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


// This macro makes the html page for the cbc report :
// used by the cbc_plot_ifar command
// Author : Francesco Salemi
{
#ifdef REDSHIFT
  bool Redshift = 1;
#endif

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  // create body.html file

  // if CWB_DOC_URL is define then man infos are added to web pages
  TString cwb_doc_url = "";
  if (gSystem->Getenv("CWB_DOC_URL") != NULL) {
    cwb_doc_url = TString(gSystem->Getenv("CWB_DOC_URL"));
  }

  ofstream out;
  char fileout[256];

  sprintf(netdir, "%s/%s", pp_dir, pp_data_dir);
  TString PP_DATA_PATH = netdir;
  sprintf(netdir, "%s", pp_data_dir);
  TString PP_DATA_DIR = netdir;
  TString HTML_DATA_DIR = "html"; // is it needed?
  // char file[2048];

  // convert eps plots to gif plots
  vector<TString> epsList =
      CWB::Toolbox::getFileListFromDir(PP_DATA_PATH, ".eps");
  for (int i = 0; i < epsList.size(); i++) {
    TString oFile = epsList[i];
    oFile.ReplaceAll(".eps", ".png");
    char cmd[2048];
    sprintf(cmd, "convert %s -resize 996x774 %s", epsList[i].Data(),
            oFile.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    cout << "Removing eps files...." << endl;
    sprintf(cmd, "rm %s", epsList[i].Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  // SORTING of ascii files

#ifdef WRITE_ASCII
  char file_notsorted[1024];
  sprintf(file_notsorted, "%s/recovered_signals.txt", PP_DATA_PATH.Data());
  char file_sorted[1024];
  sprintf(file_sorted, "%s/recovered_signals_GPSsorted.txt",
          PP_DATA_PATH.Data());
  char exec[1024];
  sprintf(exec, "sort -g --key=1 %s > %s", file_notsorted, file_sorted);
  gSystem->Exec(exec);
  sprintf(file_notsorted, "%s/injected_signals.txt", PP_DATA_PATH.Data());
  sprintf(file_sorted, "%s/injected_signals_GPSsorted.txt",
          PP_DATA_PATH.Data());
  sprintf(exec, "sort -g --key=1 %s > %s", file_notsorted, file_sorted);
  gSystem->Exec(exec);

  TString xml;
// cout << XML[0] << endl;
/*       for (int i = 1; i < nfactor + 1; i++) {
       // if (Redshift) {
                xml = XML[i - 1];
                xml.ReplaceAll(".xml", ".found.txt");
                sprintf(file_notsorted, "%s/%s", PP_DATA_PATH.Data(),
                        xml.Data());
                xml.ReplaceAll(".found.txt", "_sorted.found.txt");
                sprintf(file_sorted, "%s/%s", PP_DATA_PATH.Data(),
                        xml.Data());

        //} else {
                sprintf(file_notsorted, "%s/recovered_signals_%d.txt",
                        PP_DATA_PATH.Data(), i);
                sprintf(file_sorted,
                        "%s/recovered_signals_%d_GPSsorted.txt",
                        PP_DATA_PATH.Data(), i);
        //}
        sprintf(exec, "sort -g --key=1 %s > %s", file_notsorted,
                file_sorted);
        gSystem->Exec(exec);
        if (Redshift) {
                xml = XML[i - 1];
                xml.ReplaceAll(".xml", ".inj.txt");
                sprintf(file_notsorted, "%s/%s", PP_DATA_PATH.Data(),
                        xml.Data());
                xml.ReplaceAll(".inj.txt", "_sorted.inj.txt");
                sprintf(file_sorted, "%s/%s", PP_DATA_PATH.Data(),
                        xml.Data());

        } else {
                sprintf(file_notsorted, "%s/injected_signals_%d.txt",
                        PP_DATA_PATH.Data(), i);
                sprintf(file_sorted,
                        "%s/injected_signals_%d_GPSsorted.txt",
                        PP_DATA_PATH.Data(), i);
        //}
        sprintf(exec, "sort -g --key=1 %s > %s", file_notsorted,
                file_sorted);
        gSystem->Exec(exec);
} */
#endif

  network *net = NULL;
  CWB::config *cfg = new CWB::config;
  strcpy(cfg->tmp_dir, "tmp");
  CWB::mdc MDC(net);
  // load MDC setup  (define MDC set)
  // nfactor=4;
  CWB_PLUGIN_EXPORT(MDC)

  int gIFACTOR = -1;
  IMPORT(int, gIFACTOR)
  CWB_PLUGIN_EXPORT(net)
  CWB_PLUGIN_EXPORT(cfg)
  CWB_PLUGIN_EXPORT(gIFACTOR)

  // --------------------------------------------------------------------
  // create main html body
  // --------------------------------------------------------------------

  sprintf(fileout, "%s/main_body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout, ios::out);
  if (!out.good()) {
    cout << "cwb_mkhtml_cbc.C : Error Opening File : " << fileout << endl;
    exit(1);
  }

  out << "<br><br>" << endl;
  MakePlotsHtmlCellTable(
      &out, "Effective Radii", "data/Effective_radius.png",
      "Effective radius: m1 vs m2", "data/Effective_radius_chi.png",
      "Effective radius: M<sub>total</sub> vs &chi;<sub>eff</sub>");
  MakePlotsHtmlCellTable(&out, "Sensitive Distance", "data/ROC_IFAR_Mtot.png",
                         "Receiver Operating Curves (ROCs)",
                         "data/Distance_vs_total_mass_ifar.png",
                         "Recovered events vs FAR");
  MakePlotsHtmlCellTable(&out, "Injected distance distribution",
                     "data/distance_distribution.png", "","", "");

  out.close();

  // --------------------------------------------------------------------
  // create ROC html body
  // --------------------------------------------------------------------

  sprintf(fileout, "%s/ROC_body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout, ios::out);
  if (!out.good()) {
    cout << "cwb_mkhtml_cbc.C : Error Opening File : " << fileout << endl;
    exit(1);
  }

  out << "<br><br>" << endl;
  MakePlotsHtmlCellTable(&out, "Sensitive distance for M<sub>total</sub> bins ",
                         "data/ROC_IFAR_Mtot.png", "Inverse False Alarm Rate",
                         "data/ROC_rho1_Mtot.png", "Magnitude Test Statistic");
  MakePlotsHtmlCellTable(&out, "Sensitive distance for M<sub>chirp</sub> bins ",
                         "data/ROC_IFAR_chirp.png", "Inverse False Alarm Rate",
                         "data/ROC_rho1_chirp.png", "Magnitude Test Statistic");
  MakePlotsHtmlCellTable(&out, "Sensitive distance for &eta; bins ",
                         "data/ROC_IFAR_eta.png", "Inverse False Alarm Rate",
                         "data/ROC_rho1_eta.png", "Magnitude Test Statistic");
  MakePlotsHtmlCellTable(&out, "Sensitive distance for Cos(&iota;) bins ",
                         "data/ROC_IFAR_iota.png", "Inverse False Alarm Rate",
                         "data/ROC_rho1_iota.png", "Magnitude Test Statistic");
  MakePlotsHtmlCellTable(
      &out, "Sensitive distance for &chi;<sub>eff</sub> bins ",
      "data/ROC_IFAR_chieff.png", "Inverse False Alarm Rate",
      "data/ROC_rho1_chieff.png", "Magnitude Test Statistic");
  MakePlotsHtmlCellTable(&out, "Sensitive distance for &chi;<sub>p</sub> bins ",
                         "data/ROC_IFAR_chip.png", "Inverse False Alarm Rate",
                         "data/ROC_rho1_chip.png", "Magnitude Test Statistic");
  // MakePlotsHtmlTable(&out,"Injected distance
  // distribution","data/Injected_distances_distribution.png","580");

  out.close();

  // --------------------------------------------------------------------
  // create parspace html body
  // --------------------------------------------------------------------

  sprintf(fileout, "%s/parspace_body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout, ios::out);
  if (!out.good()) {
    cout << "cwb_mkhtml_cbc.C : Error Opening File : " << fileout << endl;
    exit(1);
  }

  out << "<br><br>" << endl;
  MakePlotsHtmlCellTable(
      &out, "Injected distributions for various signal parameters",
      "data/mtot_distribution.png", "Total mass distribution",
      "data/mchirp_distribution.png", "Chirp mass distribution");
  MakePlotsHtmlCellTable(&out, "", "data/iota_distribution.png",
                         "Inclination distribution", "data/eta_distribution.png",
                         "Symmetric Mass ratio distribution");
                         
  MakePlotsHtmlCellTable(&out, "", "data/chieff_distribution.png",
                         "&chi;<sub>eff</sub>  distribution",
                         "data/chip_distribution.png",
                         "&chi;<sub>p</sub> distribution");
  // MakePlotsHtmlCellTable(&out,"SNR
  // plots","data/Injected_snr_distributions.png","","data/Estimated_snr_vs_Injected_snr.png","");
  // MakePlotsHtmlTable(&out,"Relative snr
  // loss","data/Relative_snr_Loss.png","580");

  out.close();

  // --------------------------------------------------------------------
  // create distance html body
  // --------------------------------------------------------------------

  sprintf(fileout, "%s/distance_body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout, ios::out);
  if (!out.good()) {
    cout << "cwb_mkhtml_cbc.C : Error Opening File : " << fileout << endl;
    exit(1);
  }

  out << "<br><br>" << endl;
  MakePlotsHtmlCellTable(&out, "Distance plots", "data/Distance_vs_mtot.png",
                         "Distance vs Total mass",
                         "data/Distance_vs_mchirp.png",
                         "Distance vs Chirp mass");
  MakePlotsHtmlCellTable(&out, "", "data/Distance_vs_eta.png",
                         "Distance vs Mass Ratio", 
                         "data/Distance_vs_iota.png",
                         "Distance vs Cos(&iota;) Inclination");
  MakePlotsHtmlCellTable(&out, "", "data/Distance_vs_chieff.png",
                         "Distance vs &chi;<sub>eff</sub>",
                         "data/Distance_vs_chip.png",
                         "Distance vs &chi;<sub>p</sub>");
                         
  MakePlotsHtmlTable(&out, "Range vs rho", "data/Range.png", "580");

  out.close();
  // --------------------------------------------------------------------
  // create SNR html body
  // --------------------------------------------------------------------

  sprintf(fileout, "%s/snr_body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout, ios::out);
  if (!out.good()) {
    cout << "cwb_mkhtml_cbc.C : Error Opening File : " << fileout << endl;
    exit(1);
  }

  out << "<br><br>" << endl;
  MakePlotsHtmlCellTable(&out, "SNR plots",
                         "data/Injected_snr_distributions.png", "",
                         "data/Estimated_snr_vs_Injected_snr.png", "");
  MakePlotsHtmlTable(&out, "Relative snr loss", "data/Relative_snr_Loss.png",
                     "580");

  out.close();

  exit(0);
}
