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

/**********************************************************
 * Package:      CBCTool Class Library
 * File name:    CBCTool.hh
 * Author:       Francesco Salemi (francesco.salemi@aei.mpg.de)
 **********************************************************/

#ifndef CBCTOOL_HH
#define CBCTOOL_HH

#include "Math/BrentRootFinder.h"
#include "Math/WrappedTF1.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraphSmooth.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TRatioPlot.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TTree.h"
#include "TTreeIndex.h"

#include "TGraph2DErrors.h"
#include "TPaveText.h"
#include "TText.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Meyer.hh"
#include "netevent.hh"
#include "network.hh"
#include "wavearray.hh"
#include "wavecomplex.hh"

#include "FrameL.h"

#include "History.hh"
#include "Toolbox.hh"
#include "Toolfun.hh"

#define LST_TREE_NAME "frl"

using namespace std;

namespace CWB {

class CBCTool {
   public:
    /* ******************************** */
    /* * post processing methods      * */
    /* ******************************** */

    static Long64_t GetTreeIndex(TTree* tree, const char* selection);

    static void doROCPlot(int bkg_entries,
                          double* rho_bkg,
                          int* index,
                          float RHO_BIN,
                          double liveTot,
                          TF1* AverageRad,
                          TCanvas* c1,
                          TString odir,
                          bool write_ascii);

    static void doROCPlot(int bkg_entries,
                          double* rho_bkg,
                          int* index,
                          float RHO_BIN,
                          float RHO_NBINS,
                          float RHO_MIN,
                          double liveTot,
                          double* Rrho,
                          double* eRrho,
                          double* Trho,
                          TCanvas* c1,
                          TString odir,
                          bool write_ascii);

    static TF1* doRangePlot(int RHO_NBINS,
                            double* Trho,
                            double* Rrho,
                            double* eRrho,
                            float RHO_MIN,
                            float T_out,
                            TCanvas* c1,
                            TString networkname,
                            TString odir,
                            bool write_ascii);

    static double getLiveTime(int nIFO,
                              TChain& liv,
                              int lag_number,
                              int slag_number,
                              int dummy);

    static double getLiveTime2(TChain& liv);

    static double getZeroLiveTime(int nIFO, TChain& liv);

    static double getZeroLiveTime2(int nIFO, TChain& liv);

    static double getFAR(float rho, TH1* hc, double liveTot);

    static TH2F* FillSLagHist(int NIFO_MAX,
                              TChain& live,
                              int NSlag,
                              double SlagMin,
                              double SlagMax);

    static void doChirpFARPlot(int sel_events,
                               double* recMchirp,
                               double* injMchirp,
                               double* far,
                               TCanvas* c1,
                               TString odir);

    static double calc_isco_radius(double a);

    static double calc_isco_freq(double a);

    static double _final_spin_diff(double a_f,
                                   double eta,
                                   double delta_m,
                                   double S,
                                   double Delta);

    static double _final_spin_diff(Double_t* x, Double_t* par);

    static double bbh_final_mass_and_spin_non_precessing(double m1,
                                                         double m2,
                                                         double chi1,
                                                         double chi2);

    static double chip(double m1,
                       double m2,
                       double s1x,
                       double s1y,
                       double s1z,
                       double s2x,
                       double s2y,
                       double s2z);

    static void AddChip(TString filein, TString treename);

    static void CreateDistanceParplots(char* sim_file_name,
                                       char* mdc_file_name,
                                       char* netdir,
                                       TString opt,
                                       double MINX,
                                       double MAXX,
                                       double MAXDISTANCE,
                                       int NBIN_DIST,
                                       float T_ifar,
                                       float T_win,
                                       int nIFO);

    static TGraphErrors* CreateGraphRadiusIFAR(char* sim_file_name,
                                               char* mdc_file_name,
                                               TString SEL,
                                               float shell_volume,
                                               Color_t color,
                                               TString opt,
                                               double liveTot,
                                               float T_ifar,
                                               float T_win,
                                               int TRIALS,
                                               int nIFO,
                                               float VT,
                                               float Tscale);

   private:
    ClassDef(CBCTool, 1)
};

}  // namespace CWB

#endif
