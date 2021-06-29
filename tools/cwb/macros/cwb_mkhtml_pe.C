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


// make the html pe report : used by the cwb_report command
{

  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
  CWB::Toolbox::checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

  // create body.html file

  sprintf(netdir,"%s/%s",pp_dir,pp_data_dir);
  TString PP_DATA_PATH=netdir;
  sprintf(netdir,"%s",pp_data_dir);
  TString PP_DATA_DIR=netdir;

  // if CWB_DOC_URL is define then man infos are added to web pages
  TString cwb_doc_url="";
  if(gSystem->Getenv("CWB_DOC_URL")!=NULL) {
    cwb_doc_url=TString(gSystem->Getenv("CWB_DOC_URL"));
  }

  ofstream out;
  char fileout[256];

  // --------------------------------------------------------------------
  // create prc html body
  // --------------------------------------------------------------------

  sprintf(fileout,"%s/prc_body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout,ios::out);
  if (!out.good()) {cout << "cwb_mkhtml_pe.C : Error Opening File : " << fileout << endl;exit(1);}

  MakePlotsHtmlTable(&out,"antenna pattern (same psd) : sensitivity & alignment",
                     "data/antpat_sensitivity.png","data/antpat_alignment.png");
  MakePlotsHtmlTable(&out,"injected vs reconstructed locations of the detected events",
                     "data/inj_detected_antpat.png","data/rec_detected_antpat.png");
  MakePlotsHtmlTable(&out,"injected vs detected events distributions : SNRnet & distance",
                     "data/inj_vs_rec_snr.png","data/inj_vs_rec_distance.png");
  MakePlotsHtmlTable(&out,"search area & pp plot","data/search_area.png","data/pp_plot_prc.png");
  MakePlotsHtmlTable(&out,"cos_theta distribution","data/cos_theta.png","580");
  MakePlotsHtmlTable(&out,"Sky Localization Median vs Network SNR","data/median50_vs_snr.png","data/median90_vs_snr.png");

  out.close();

  // --------------------------------------------------------------------
  // create wrc html body
  // --------------------------------------------------------------------

  sprintf(fileout,"%s/wrc_body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout,ios::out);
  if (!out.good()) {cout << "cwb_mkhtml_pe.C : Error Opening File : " << fileout << endl;exit(1);}

  MakePlotsHtmlTable(&out,"Reconstructed SNRnet vs Injected SNRnet","data/osnr_vs_isnr.png","580");
  MakePlotsHtmlTable(&out,"Residual Energy (Distrbution/PP-plot)","data/nre_vs_isnr.png","data/pp_plot_nre.png");
  MakePlotsHtmlTable(&out,"Overlap/Fitting Factors","data/of_vs_isnr.png","data/ff_vs_isnr.png");
  MakePlotsHtmlTable(&out,"","data/of_vs_ff.png","data/of_vs_nre.png");
  MakePlotsHtmlTable(&out,"Frequency Residual Energy (Distrbution/PP-plot)","data/pp_plot_fre.png","580");

  out.close();

  // --------------------------------------------------------------------
  // create mchirp html body
  // --------------------------------------------------------------------

  sprintf(fileout,"%s/mchirp_body.html", pp_dir);
  cout << fileout << endl;
  out.open(fileout,ios::out);
  if (!out.good()) {cout << "cwb_mkhtml_pe.C : Error Opening File : " << fileout << endl;exit(1);}

  MakePlotsHtmlTable(&out,"Chirp Mass","data/omch_vs_imch.png","data/dmch_vs_isnr.png");
  MakePlotsHtmlTable(&out,"Chirp Mass Error","data/emch_vs_imch.png","data/emch_vs_isnr.png");
  MakePlotsHtmlTable(&out,"Chirp Ellipticity","data/elch_vs_imch.png","data/elch_vs_isnr.png");
  MakePlotsHtmlTable(&out,"Chirp Pixel Fraction","data/pfch_vs_imch.png","data/pfch_vs_isnr.png");
  MakePlotsHtmlTable(&out,"Chirp Energy Fraction","data/efch_vs_imch.png","data/efch_vs_isnr.png");
  MakePlotsHtmlTable(&out,"PE Chirp Mass","data/opemch_vs_isnr.png","data/epemch_vs_isnr.png");

  out.close();


  exit(0);
}

