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

#include "CBCTool.hh"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "watversion.hh"
//#include "GCBCTool.hh"
#include <Riostream.h>
#include "Math/BrentRootFinder.h"
#include "Math/WrappedTF1.h"
#include "TF1.h"
#include "TThread.h"

#define CCAST(PAR) const_cast<char*>(PAR)

// Without this macro the THtml doc for CWB::CBCTool can not be generated
ClassImp(CWB::CBCTool)
// int compareSegments(waveSegment a, waveSegment b) {return a.start < b.start;}

// definitions used by CWB::CBCTool::mergeCWBTrees with threads

#define MAX_THREADS 8

    ////______________________________________________________________________________

    Long64_t CWB::CBCTool::GetTreeIndex(TTree* tree, const char* selection) {

    tree->Draw("Entry$>>hist(Entries$,0,Entries$)", selection, "goff");
    TH1I* hist = (TH1I*)gDirectory->Get("hist");
    Long64_t iEntry = hist->GetBinLowEdge(hist->FindFirstBinAbove(0));
    delete hist;
    return iEntry;
}

//______________________________________________________________________________
/*
TH1* CWB::CBCTool::GetCumulative(TH1* in, Bool_t forward)
{
   //  Return a pointer to an histogram containing the cumulative The
   //  cumulative can be computed both in the forward (default) or backward
   //  direction; the name of the new histogram is constructed from
   //  the name of this histogram with the suffix suffix appended.
   //
   // The cumulative distribution is formed by filling each bin of the
   // resulting histogram with the sum of that bin and all previous
   // (forward == kTRUE) or following (forward = kFALSE) bins.
   //
   // note: while cumulative distributions make sense in one dimension, you
   // may not be getting what you expect in more than 1D because the concept
   // of a cumulative distribution is much trickier to define; make sure you
   // understand the order of summation before you use this method with
   // histograms of dimension >= 2.

   const Int_t nbinsx = in->GetNbinsX();
   const Int_t nbinsy = in->GetNbinsY();
   const Int_t nbinsz = in->GetNbinsZ();
   TH1* hintegrated = (TH1*) in->Clone();
   hintegrated->Reset();
   if (forward) { // Forward computation
      Double_t sum = 0.;
      for (Int_t binz = 1; binz <= nbinsz; ++binz) {
         for (Int_t biny = 1; biny <= nbinsy; ++biny) {
            for (Int_t binx = 1; binx <= nbinsx; ++binx) {
               const Int_t bin = hintegrated->GetBin(binx, biny, binz);
               sum += in->GetBinContent(bin);
               hintegrated->SetBinContent(bin, sum);
            }
         }
      }
   } else { // Backward computation
      Double_t sum = 0.;
      for (Int_t binz = nbinsz; binz >= 1; --binz) {
         for (Int_t biny = nbinsy; biny >= 1; --biny) {
            for (Int_t binx = nbinsx; binx >= 1; --binx) {
               const Int_t bin = hintegrated->GetBin(binx, biny, binz);
               sum += in->GetBinContent(bin);
               hintegrated->SetBinContent(bin, sum);
            }
         }
      }
   }
   return hintegrated;
}

*/

//______________________________________________________________________________

void CWB::CBCTool::doROCPlot(int bkg_entries,
                             double* rho_bkg,
                             int* index,
                             float RHO_BIN,
                             double liveTot,
                             TF1* AverageRad,
                             TCanvas* c1,
                             TString odir,
                             bool write_ascii) {
    //
    // It produces ROC Plot from background and simulation ntuples
    //
    // Input
    //       bkg_entries : number of background events
    //       rho_bkg : array of doubles containing the events rho
    //       index : array of integers containing the sorted entries in
    //       ascending order  // MY: too much useless stuff, need some
    //       clening and rearrangement! RHO_BIN : float containing the rho
    //       bin  //MY: remove it? it should be defined in user-pp conf file
    //       liveTot : double containing the total background live time in
    //       seconds as calculated from the liveTime ntuple AverageRad :
    //       fitting function of the measured average radius as a function
    //       of rho (a power-law with some exponential, A*rho^B*exp(-C))
    //	 calculated by the following doRangePlot method
    //	 BEWARE, the fit does a decent job on all sofar tested
    // waveforms, but no check is done! Look at the Range.png. 	 The reasoning
    // in favour of using the radius from fitting in place of the measured
    // one (which produces some overestimate bias at low rho)
    //       is to have a better estimate at high rho, where the statistics
    //       is poor.
    //	 c1 : the formatted canvas used for all plots by cbc_plots.C
    //	 odir : output dir
    //	 write_ascii : boolean to toggle on/off the writing of ascii
    // output

    int nbins =
        -TMath::Nint(rho_bkg[index[bkg_entries - 1]] - rho_bkg[index[0]]) /
        RHO_BIN;
    cout << "Rho max : " << rho_bkg[index[0]]
         << " Rho min : " << rho_bkg[index[bkg_entries - 1]]
         << " nbins : " << nbins << endl;

    // definition and filling of the events rho histogram
    TH1D* h = new TH1D("h", "h", nbins, rho_bkg[index[bkg_entries - 1]],
                       rho_bkg[index[0]] * 1.05);
    for (int i = 0; i < bkg_entries; i++) {
        h->Fill(rho_bkg[i]);
    }
    // hc is the cumulative histogram of h
    // TH1* hc = h->GetCumulative(kFALSE);
    const Int_t nbinsx = h->GetNbinsX();
    const Int_t nbinsy = h->GetNbinsY();
    const Int_t nbinsz = h->GetNbinsZ();
    TH1* hc = (TH1*)h->Clone();
    hc->Reset();
    Double_t sum = 0.;
    for (Int_t binz = nbinsz; binz >= 1; --binz) {
        for (Int_t biny = nbinsy; biny >= 1; --biny) {
            for (Int_t binx = nbinsx; binx >= 1; --binx) {
                const Int_t bin = hc->GetBin(binx, biny, binz);
                sum += h->GetBinContent(bin);
                hc->SetBinContent(bin, sum);
            }
        }
    }

    double far[nbins], efar[nbins], rad[nbins], rhobin[nbins], erad[nbins];
    for (int i = 0; i < nbins; i++) {
        far[i] = (double)hc->GetBinContent(i) / liveTot;
        efar[i] = (double)hc->GetBinError(i) / liveTot;
        rhobin[i] = hc->GetBinCenter(
            i);  // MY : should I consider the left bound of the bin?
                 // this choice should be more conservative
        rad[i] = AverageRad->Eval(
            rhobin[i]);  // BEWARE: no guarantee that the fitting
                         // function returns correct estimates
        erad[i] = 0.0;
        //		cout << "(far, rho, radius) -
        //"<<far[i]<<","<<hc->GetBinCenter(i)<<","<<rad[i]<<" ";
    }
    // cout<<endl;
    // cout << "FAR max : "<< far[0] << " FAR min : "<<far[nbins-1] <<endl;
    c1->Clear();
    c1->SetLogy(1);
    c1->SetLogx(1);

    hc->Scale(1. / liveTot);
    hc->GetYaxis()->SetTitle("FAR [Hz]");
    hc->GetXaxis()->SetTitle("#rho");
    hc->Draw("LP");
    char fname[1024];
    sprintf(fname, "%s/FAR.png", odir.Data());
    c1->Update();
    c1->SaveAs(fname);  // A dump of the FAR used for the ROC

    c1->Clear();
    c1->SetLogx(1);
    c1->SetLogy(0);        // MY: log?
    c1->SetTheta(89.999);  // the setting of the theta and phi is for the
                           // viewing angle: there is no TGraphError with "pcol"
                           // option, so the trick is to use the 2D version and
                           // use a viewing angle the closest as possible to the
                           // vertical (Drawbacks: some artefacts due to that
                           // can be noticed; no gridlines)
    c1->SetPhi(0.0001);

    TGraph2DErrors* gr2 =
        new TGraph2DErrors(nbins, far, rad, rhobin, efar, erad, erad);
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(20);
    gr2->SetMinimum(0.0);
    double RHO_MAX = ceil(rho_bkg[index[0]] * 1.05);
    gr2->SetMaximum(RHO_MAX);
    // gr1->Draw("ALP");
    // TMultiGraph *multi = new TMultiGraph();
    // multi->Add(gr2);
    // multi->Paint("AP");
    gr2->SetTitle("ROC Curve: Average Range vs FAR @rho threshold");
    gr2->GetHistogram()->GetXaxis()->SetTitle("FAR [Hz]");
    gr2->GetHistogram()->GetZaxis()->SetTitle("#rho");
    //  multi->GetHistogram()->GetXaxis()->SetRangeUser(RHO_MIN,10.);
    //  multi->GetHistogram()->GetXaxis()->SetRangeUser(1e-12, 1e-7);
    gr2->GetHistogram()->GetXaxis()->SetLimits(
        1e-12, 1e-7);  // MY: Add similar line for the y axis? declare the X
                       // and Y bounds in the definition
    gr2->GetHistogram()->GetYaxis()->SetTitle("Average Range [Mpc]");
    gr2->GetXaxis()->SetTitleOffset(1.3);
    gr2->GetYaxis()->SetTitleOffset(1.25);
    gr2->GetXaxis()->SetTickLength(0.02);
    gr2->GetYaxis()->SetTickLength(0.01);
    gr2->GetXaxis()->CenterTitle(kTRUE);
    gr2->GetYaxis()->CenterTitle(kTRUE);
    gr2->GetZaxis()->CenterTitle(kTRUE);
    gr2->GetXaxis()->SetTitleFont(42);
    gr2->GetXaxis()->SetLabelFont(42);
    gr2->GetYaxis()->SetTitleFont(42);
    gr2->GetYaxis()->SetLabelFont(42);
    //  gr2->GetYaxis()->SetMoreLogLabels(kTRUE);
    gr2->GetYaxis()->SetNoExponent(kTRUE);

    gr2->Draw("err zcolpcol");
    sprintf(fname, "%s/ROC.png", odir.Data());
    c1->Update();
    c1->SaveAs(fname);

    if (write_ascii) {
        sprintf(fname, "%s/ROC.txt", odir.Data());
        FILE* frange = fopen(fname, "w");
        fprintf(frange, "#rho	FAR[Hz] eFAR range[Mpc] erange\n");
        for (int i = 0; i < nbins; i++) {
            fprintf(frange, "%3.3g %4.3g %4.3g %4.4g %4.4g\n",
                    hc->GetBinCenter(i), far[i], efar[i], rad[i], erad[i]);
        }
        fclose(frange);
    }

    delete h, hc, gr2;

    return;
}

//______________________________________________________________________________

void CWB::CBCTool::doROCPlot(int bkg_entries,
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
                             bool write_ascii) {
    //
    // It produces ROC Plot from background and simulation ntuples
    //
    // Input
    //       bkg_entries : number of background events
    //       rho_bkg : array of doubles containing the events rho
    //       index : array of integers containing the sorted entries in
    //       ascending order  // MY: too much useless stuff, need some
    //       clening and rearrangement! RHO_BIN : float containing the rho
    //       bin  //MY: remove it? it should be defined in user-pp conf file
    //       liveTot : double containing the total background live time in
    //       seconds as calculated from the liveTime ntuple Rrho and eRrho :
    //       arrays of doubles containing the Radius and error radius
    //       estimates at fixed rho thresholds Trho
    //
    //	 c1 : the formatted canvas used for all plots by cbc_plots.C
    //	 odir : output dir
    //	 write_ascii : boolean to toggle on/off the writing of ascii
    // output

    // int nbins =
    // -TMath::Nint(rho_bkg[index[bkg_entries-1]]-rho_bkg[index[0]])/RHO_BIN;
    int nbins = (int)RHO_NBINS;
    cout << "Rho max : " << rho_bkg[index[0]]
         << " Rho min : " << rho_bkg[index[bkg_entries - 1]]
         << " nbins : " << nbins << endl;

    // definition and filling of the events rho histogram
    TH1D* h = new TH1D("h", "h", nbins, RHO_MIN, RHO_NBINS * RHO_BIN);
    for (int i = 0; i < bkg_entries; i++) {
        h->Fill(rho_bkg[i]);
    }
    // hc is the cumulative histogram of h
    // TH1* hc = h->GetCumulative(kFALSE);

    const Int_t nbinsx = h->GetNbinsX();
    const Int_t nbinsy = h->GetNbinsY();
    const Int_t nbinsz = h->GetNbinsZ();
    TH1* hc = (TH1*)h->Clone();
    hc->Reset();
    Double_t sum = 0.;
    for (Int_t binz = nbinsz; binz >= 1; --binz) {
        for (Int_t biny = nbinsy; biny >= 1; --biny) {
            for (Int_t binx = nbinsx; binx >= 1; --binx) {
                const Int_t bin = hc->GetBin(binx, biny, binz);
                sum += h->GetBinContent(bin);
                hc->SetBinContent(bin, sum);
            }
        }
    }

    double far[nbins], efar[nbins], rad[nbins], rhobin[nbins], erad[nbins],
        vol[nbins], evol[nbins], erhobin[nbins];
    int j = 0;
    for (int i = 0; i < nbins; i++) {
        far[i] = (double)hc->GetBinContent(i + 1) / liveTot * 365. * 86400.;
        efar[i] = (double)hc->GetBinError(i + 1) / liveTot * 365. * 86400.;
        rhobin[i] = hc->GetBinCenter(
            i + 1);  // MY : should I consider the left bound of the
                     // bin? this choice should be more conservative
        while (rhobin[i] > Trho[j]) {
            j++;
        }
        erhobin[i] = RHO_BIN;
        rad[i] = Rrho[j];
        erad[i] = eRrho[j];
        vol[i] = 4. / 3. * TMath::Pi() * pow(Rrho[j], 3.);
        evol[i] = 4. * TMath::Pi() * pow(Rrho[j], 2.) * eRrho[j];
        //		cout << "(far, rho, radius) -
        //"<<far[i]<<","<<hc->GetBinCenter(i)<<","<<rad[i]<<" ";
    }
    // cout<<endl;
    // cout << "FAR max : "<< far[0] << " FAR min : "<<far[nbins-1] <<endl;
    c1->Clear();
    c1->SetLogy(1);
    c1->SetLogx(1);

    hc->Scale(365. * 86400. / liveTot);
    hc->GetYaxis()->SetTitle("FAR [year^{-1}]");
    hc->GetXaxis()->SetTitle("#rho");
    hc->Draw("LP");
    char fname[1024];
    sprintf(fname, "%s/FAR.png", odir.Data());
    c1->Update();
    c1->SaveAs(fname);  // A dump of the FAR used for the ROC

    c1->Clear();
    c1->SetLogx(1);
    c1->SetLogy(0);        // MY: log?
    c1->SetTheta(89.999);  // the setting of the theta and phi is for the
                           // viewing angle: there is no TGraphError with "pcol"
                           // option, so the trick is to use the 2D version and
                           // use a viewing angle the closest as possible to the
                           // vertical (Drawbacks: some artefacts due to that
                           // can be noticed; no gridlines)
    c1->SetPhi(0.0001);

    TGraph2DErrors* gr2 =
        new TGraph2DErrors(nbins, far, rad, rhobin, efar, erad, erhobin);
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(20);
    gr2->SetMinimum(0.0);
    double RHO_MAX = ceil(rho_bkg[index[0]] * 1.05);
    gr2->SetMaximum(RHO_MAX);
    // gr1->Draw("ALP");
    // TMultiGraph *multi = new TMultiGraph();
    // multi->Add(gr2);
    // multi->Paint("AP");
    gr2->SetTitle("ROC Curve: Average Range vs FAR @rho threshold");
    gr2->GetHistogram()->GetXaxis()->SetTitle("FAR [year^{-1}]");
    gr2->GetHistogram()->GetZaxis()->SetTitle("#rho");
    //  multi->GetHistogram()->GetXaxis()->SetRangeUser(RHO_MIN,10.);
    //  multi->GetHistogram()->GetXaxis()->SetRangeUser(1e-12, 1e-7);
    gr2->GetHistogram()->GetXaxis()->SetLimits(
        1e-5, 1e1);  // MY: Add similar line for the y axis? declare the X
                     // and Y bounds in the definition
    gr2->GetHistogram()->GetYaxis()->SetTitle("Average Range [Mpc]");
    gr2->GetXaxis()->SetTitleOffset(1.3);
    gr2->GetYaxis()->SetTitleOffset(1.25);
    gr2->GetXaxis()->SetTickLength(0.02);
    gr2->GetYaxis()->SetTickLength(0.01);
    gr2->GetXaxis()->CenterTitle(kTRUE);
    gr2->GetYaxis()->CenterTitle(kTRUE);
    gr2->GetZaxis()->CenterTitle(kTRUE);
    gr2->GetXaxis()->SetTitleFont(42);
    gr2->GetXaxis()->SetLabelFont(42);
    gr2->GetYaxis()->SetTitleFont(42);
    gr2->GetYaxis()->SetLabelFont(42);
    //  gr2->GetYaxis()->SetMoreLogLabels(kTRUE);
    gr2->GetYaxis()->SetNoExponent(kTRUE);

    gr2->Draw("err zcolpcol");
    sprintf(fname, "%s/ROC.png", odir.Data());
    c1->Update();
    c1->SaveAs(fname);

    if (write_ascii) {
        sprintf(fname, "%s/ROC.txt", odir.Data());
        FILE* frange = fopen(fname, "w");
        fprintf(frange,
                "#rho	FAR[year^{-1}] eFAR range[Mpc] "
                "erange\n");
        for (int i = 0; i < nbins; i++) {
            fprintf(frange, "%3.3g %4.3g %4.3g %4.4g %4.4g\n",
                    hc->GetBinCenter(i), far[i], efar[i], rad[i], erad[i]);
        }
        fclose(frange);
    }

    TGraph2DErrors* gr3 =
        new TGraph2DErrors(nbins, far, vol, rhobin, efar, erad, erhobin);
    gr3->SetTitle("ROC Curve: Average Volume vs FAR @rho threshold");
    gr3->GetHistogram()->GetYaxis()->SetTitle("Average volume [Mpc^3]");
    gr3->GetHistogram()->GetXaxis()->SetTitle("FAR [year^{-1}]");
    gr3->GetHistogram()->GetZaxis()->SetTitle("#rho");
    gr3->GetHistogram()->GetXaxis()->SetLimits(
        1e-5, 1e2);  // MY: Add similar line for the y axis? declare the X
                     // and Y bounds in the definition

    gr3->Draw("err zcolpcol");
    sprintf(fname, "%s/ROCV.png", odir.Data());
    c1->Update();
    c1->SaveAs(fname);

    delete h, hc, gr2;

    return;
}

//______________________________________________________________________________

TF1* CWB::CBCTool::doRangePlot(int RHO_NBINS,
                               double* Trho,
                               double* Rrho,
                               double* eRrho,
                               float RHO_MIN,
                               float T_out,
                               TCanvas* c1,
                               TString networkname,
                               TString odir,
                               bool write_ascii) {
    //
    // It produces Range Plot from simulation ntuple
    //	AverageRad : fitting function of the measured average radius as
    // a function of rho (a power-law with some exponential, A*rho^B*exp(-C))
    //	 BEWARE, the fit does a decent job on all sofar tested
    // waveforms, but no check is done! Look at the Range.png. 	 The reasoning
    // in favour of using the radius from fitting in place of the measured
    // one (which produces some overestimate bias at low rho)
    //       is to have a better estimate at high rho, where the statistics
    //       is poor.
    //
    // Input
    //       RHO_NBINS : number of rho bins  // MY: remove this in favour of
    //       the standard user_pp variables? Trho : array of doubles
    //       containg the rho thresholds at the left margin of the bins
    //	 Rrho : array of doubles containing the average radius as
    // calculated at the left margin of the bin 	 eRrho : array of
    // doubles containing the error on the average radius 	 RHO_MIN : float
    // with minimal rho // MY: do I really need this? 	 T_out : float with
    // output threashold 	 c1 : the formatted canvas used for all plots by
    // cbc_plots.C 	 networkname : TString with the network name //MY: read it from
    // the ntuple? 	 odir : output dir 	 write_ascii : boolean to toggle
    // on/off the writing of ascii output

    char fname[1024];

    if (write_ascii) {
        sprintf(fname, "%s/range_threshold.txt", odir.Data());
        FILE* frange = fopen(fname, "w");
        fprintf(frange, "#rho	range[Mpc] \n");
        for (int i = 0; i < RHO_NBINS; i++) {
            fprintf(frange, "%3.2f %4.3e\n", Trho[i], Rrho[i]);
        }
        fclose(frange);
    }

    // TF1 *f1 = new TF1("f1","[0]*pow(x,[1])",T_out,10.);
    TF1* f2 = new TF1("f2", "[0]*pow(x,[1])*TMath::Exp(-x*[2])", T_out,
                      Trho[RHO_NBINS - 1]);  // Empirical function to fit
                                             // the radius vs rho
    // TF1 *f3 = new TF1("f3","[0]+[1]/pow(x,1)",T_out,10.);
    f2->SetParameters(500., -1., 0.0);

    TGraphErrors* grtmp = new TGraphErrors(RHO_NBINS, Trho, Rrho, NULL, eRrho);
    grtmp->SetMarkerStyle(20);
    grtmp->SetMarkerSize(1.0);
    grtmp->GetHistogram()->GetXaxis()->SetTitle("#rho");
    grtmp->GetHistogram()->GetXaxis()->SetRangeUser(RHO_MIN,
                                                    Trho[RHO_NBINS - 1]);
    grtmp->GetHistogram()->GetYaxis()->SetRangeUser(
        f2->Eval(Trho[RHO_NBINS - 1]), f2->Eval(T_out));
    grtmp->GetHistogram()->GetYaxis()->SetTitle("Average Range [Mpc]");
    grtmp->GetXaxis()->SetTitleOffset(1.3);
    grtmp->GetYaxis()->SetTitleOffset(1.25);
    grtmp->GetXaxis()->SetTickLength(0.01);
    grtmp->GetYaxis()->SetTickLength(0.01);
    grtmp->GetXaxis()->CenterTitle(kTRUE);
    grtmp->GetYaxis()->CenterTitle(kTRUE);
    grtmp->GetXaxis()->SetTitleFont(42);
    grtmp->GetXaxis()->SetLabelFont(42);
    grtmp->GetYaxis()->SetTitleFont(42);
    grtmp->GetYaxis()->SetLabelFont(42);
    grtmp->GetYaxis()->SetMoreLogLabels(kTRUE);
    grtmp->GetYaxis()->SetNoExponent(kTRUE);

    grtmp->Fit("f2", "R");

    c1->Clear();
    c1->SetLogy();
    c1->SetLogx();
    gStyle->SetOptFit(kTRUE);
    sprintf(fname, "%s : Range (%s) ", networkname.Data(),
            f2->GetExpFormula().Data());
    grtmp->SetTitle(fname);
    grtmp->Draw("ALP");

    sprintf(fname, "%s/Range.png", odir.Data());
    c1->Update();
    c1->SaveAs(fname);

    delete grtmp;

    return f2;
}

//______________________________________________________________________________

double CWB::CBCTool::getLiveTime(int nIFO,
                                 TChain& liv,
                                 int lag_number,
                                 int slag_number,
                                 int dummy) {
    //
    // It calculates the total live time of all time shifts.
    // If defining a lag_number/slag_number, it calculates live time for
    // these only. This is a semplified version of the method (same name)
    // that can be found in the standard toolbox: it's much faster as it
    // does less things... Input: nIFO        : detector number
    //        liv         : tree containing live time information
    //        lag_number  : specific lag number
    //        slag_number : specific slag number
    //        dummy       : used to increase calculation speed
    //
    // Examples : (lag_number==-1) && (slag_number==-1) -> all non zero lags
    // : excluded (lag==0&&slag==0)
    //            (lag_number>=0)  && (slag_number>=0)  -> select
    //            (lag==lag_number) && (slag==slag_number) (lag_number>=0)
    //            && (slag_number==-1) -> select lag=lag_number   : excluded
    //            (lag==0&&slag==0) (lag_number==-1) && (slag_number>=0)  ->
    //            select slag=slag_number : excluded (lag==0&&slag==0)
    //
    // Note: the variable dummy has been introduced to optimize
    //       speed of the method (the reason is unknown!!!)

    float xlag[NIFO_MAX + 1];
    float xslag[NIFO_MAX + 1];
    double xlive;
    double liveTot = 0.;
    double Live = 0.;

    // check if slag is presents in liv tree
    TBranch* branch;
    bool slagFound = false;
    TIter next(liv.GetListOfBranches());
    while ((branch = (TBranch*)next())) {
        if (TString("slag").CompareTo(branch->GetName()) == 0)
            slagFound = true;
    }
    next.Reset();
    if (!slagFound) {
        cout << "CWB::Toolbox::getLiveTime : Error - live tree do not "
                "contains slag leaf"
             << endl;
        gSystem->Exit(1);
    }

    //  liv.SetCacheSize(10000000);
    //  liv.AddBranchToCache("slag");
    // liv.AddBranchToCache("lag");
    // liv.AddBranchToCache("live");
    liv.SetBranchAddress("slag", xslag);
    liv.SetBranchAddress("lag", xlag);

    liv.SetBranchAddress("live", &xlive);
    liv.SetBranchStatus("gps", false);

    Long64_t ntrg = liv.GetEntries();
    for (Long64_t i = 0; i < ntrg; i++) {
        liv.GetEntry(i);

        if ((lag_number >= 0) && (slag_number >= 0)) {
            Live = (xlag[nIFO] == lag_number && xslag[nIFO] == slag_number)
                       ? xlive
                       : 0.;  // lag/slag live time
        }
        if ((lag_number >= 0) && (slag_number == -1)) {
            Live = ((xlag[nIFO] == lag_number) &&
                    !((xlag[nIFO] == 0) && (xslag[nIFO] == 0)))
                       ? xlive
                       : 0.;  // lag live time
        }
        if ((lag_number == -1) && (slag_number >= 0)) {
            Live = ((xslag[nIFO] == slag_number) &&
                    !((xlag[nIFO] == 0) && (xslag[nIFO] == 0)))
                       ? xlive
                       : 0.;  // slag live time
        }
        if ((lag_number == -1) && (slag_number == -1)) {
            Live = !((xlag[nIFO] == 0) && (xslag[nIFO] == 0))
                       ? xlive
                       : 0.;  // non-zero live time
        }

        liveTot += Live;
    }
    return liveTot;
}

//______________________________________________________________________________

double CWB::CBCTool::getLiveTime2(TChain& liv) {
    // It calculates the total live time of all time shifts (both zero lag
    // and non-zero lags) : remember to subtract the zero-lag livetime if
    // you need precise times! This is a semplified version of the previous
    // method: it only reads&sums the live branch...faster than the speed of
    // light...:) (very useful for large liveTime ntuples) Input:
    //        liv         : tree containing live time information

    double xlive;
    double liveTot = 0.;
    TBranch* LIVE = liv.GetBranch("live");
    LIVE->SetAddress(&xlive);
    Long64_t nentries = LIVE->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        LIVE->GetEntry(i);
        liveTot += xlive;
    }
    return liveTot;
}
//______________________________________________________________________________

double CWB::CBCTool::getZeroLiveTime(int nIFO, TChain& liv) {
    // It calculates the total live time of zero lag
    // This is mostly needed to correct the background  estimate from the
    // previous method

    // Input: nIFO        : detector number
    //        liv         : tree containing live time information
    double liveTot = 0.;
    char cut[256];

    sprintf(cut, "(lag[%d]==0&&slag[%d]==0)", nIFO, nIFO);
    liv.SetEstimate(-1);
    Long64_t nentries = liv.Draw("live", cut, "goff");
    double* xlive = liv.GetV1();
    for (Long64_t i = 0; i < nentries; i++) {
        liveTot += xlive[i];
    }
    return liveTot;
}

//______________________________________________________________________________

//______________________________________________________________________________

double CWB::CBCTool::getZeroLiveTime2(int nIFO, TChain& liv) {
    // It calculates the total live time of zero lag
    // This is mostly needed to correct the background  estimate from the
    // previous method. Yet another emplementation to improve speed wrt the
    // previous method (best so far) : it reads from branches just when it
    // is really needed.
    //

    // Input: nIFO        : detector number
    //        liv         : tree containing live time information

    float xlag[NIFO_MAX + 1];
    float xslag[NIFO_MAX + 1];
    double xlive;
    double liveTot = 0.;

    Long64_t nentries = liv.GetEntries();
    TBranch* lag = liv.GetBranch("lag");
    lag->SetAddress(xlag);
    TBranch* slag = liv.GetBranch("slag");
    slag->SetAddress(xslag);
    TBranch* live = liv.GetBranch("live");
    live->SetAddress(&xlive);

    for (Long64_t i = 0; i < nentries; i++) {
        // T->LoadTree(i);
        slag->GetEntry(i);
        if (xslag[nIFO] != 0)
            continue;
        lag->GetEntry(i);
        if (xlag[nIFO] != 0)
            continue;
        live->GetEntry(i);
        liveTot += xlive;
    }
    return liveTot;
}

//______________________________________________________________________________

//______________________________________________________________________________

double CWB::CBCTool::getFAR(float rho, TH1* hc, double liveTot) {
    int w = 0;
    while ((rho > (float)hc->GetBinCenter(w)) && (w < hc->GetNbinsX() - 1)) {
        w++;
    }
    double far = (double)hc->GetBinContent(w + 1) / liveTot;
    double efar = (double)hc->GetBinError(w + 1) / liveTot;

    return far, efar;
}

//______________________________________________________________________________

TH2F* CWB::CBCTool::FillSLagHist(int NIFO_MAX,
                                 TChain& live,
                                 int NSlag,
                                 double SlagMin,
                                 double SlagMax) {
    TH2F* lSlag = new TH2F("LSLAG", "Live time distribution over slags", NSlag,
                           SlagMin / 86400., SlagMax / 86400., NSlag,
                           SlagMin / 86400., SlagMax / 86400.);
    cout << "Start filling lSlag histogram : " << endl;
    //  float  xlag[NIFO_MAX+1];
    float xslag[NIFO_MAX + 1];
    double xlive;
    live.SetBranchStatus("*", 0);  // disable all branches
    live.SetBranchStatus("slag", 1);
    live.SetBranchStatus("live", 1);
    live.SetBranchAddress("slag", xslag);
    // live.SetBranchAddress("lag",xlag);
    live.SetBranchAddress("live", &xlive);
    // live.SetBranchStatus("gps",false);

    long int ntrg = live.GetEntries();
    for (long int l = 0; l < ntrg; l++) {
        live.GetEntry(l);
        lSlag->Fill(xslag[0] / 86400., xslag[1] / 86400., xlive);
        if (l % 10000000 == 0) {
            cout << scientific << (double)l / ntrg << "% - ";
        }
    }

    return lSlag;
}
//______________________________________________________________________________

void CWB::CBCTool::doChirpFARPlot(int sel_events,
                                  double* recMchirp,
                                  double* injMchirp,
                                  double* far,
                                  TCanvas* c1,
                                  TString odir) {
    c1->Clear();
    c1->SetLogx(0);
    c1->SetLogy(0);
    c1->SetLogz(1);
    c1->SetTheta(89.999);  // the setting of the theta and phi is for the
                           // viewing angle: there is no TGraphError with "pcol"
                           // option, so the trick is to use the 2D version and
                           // use a viewing angle the closest as possible to the
                           // vertical (Drawbacks: some artefacts due to that
                           // can be noticed; no gridlines)
    c1->SetPhi(0.0001);
    double efar[sel_events];
    for (int i = 0; i < sel_events; i++) {
        efar[i] = 0.0;
    }

    TGraph2DErrors* gr2 = new TGraph2DErrors(sel_events, recMchirp, injMchirp,
                                             far, efar, efar, efar);
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(20);
    // gr2->SetMinimum(1e-4);
    // gr2->SetMaximum(1e1);
    gr2->SetTitle("");
    gr2->GetHistogram()->GetXaxis()->SetTitle("Injected mchirp");
    gr2->GetHistogram()->GetYaxis()->SetTitle("Recovered mchirp");
    gr2->GetHistogram()->GetXaxis()->SetLimits(
        0., 70.);  // MY: Add similar line for the y axis? declare the X and
                   // Y bounds in the definition
    gr2->GetHistogram()->GetYaxis()->SetLimits(
        0., 70.);  // MY: Add similar line for the y axis? declare the X and
                   // Y bounds in the definition
    gr2->GetHistogram()->GetZaxis()->SetTitle("False Alarm Rate [years^{-1}]");
    gr2->GetXaxis()->SetTitleOffset(1.3);
    gr2->GetYaxis()->SetTitleOffset(1.25);
    gr2->GetXaxis()->SetTickLength(0.02);
    gr2->GetYaxis()->SetTickLength(0.01);
    gr2->GetXaxis()->CenterTitle(kTRUE);
    gr2->GetYaxis()->CenterTitle(kTRUE);
    gr2->GetZaxis()->CenterTitle(kTRUE);
    gr2->GetXaxis()->SetTitleFont(42);
    gr2->GetXaxis()->SetLabelFont(42);
    gr2->GetYaxis()->SetTitleFont(42);
    gr2->GetYaxis()->SetLabelFont(42);
    //  gr2->GetYaxis()->SetMoreLogLabels(kTRUE);
    gr2->GetYaxis()->SetNoExponent(kTRUE);

    gr2->Draw("err zcolpcol");
    char fname[1024];
    sprintf(fname, "%s/ChirpFAR.png", odir.Data());
    c1->Update();
    c1->SaveAs(fname);

    return;
}

//______________________________________________________________________________
double CWB::CBCTool::calc_isco_radius(double a) {
    // Calculate the ISCO radius of a Kerr BH as a function of the Kerr
    // parameter a : Kerr parameter
    // Ref. Eq. (2.5) of Ori, Thorne Phys Rev D 62 124022 (2000)

    double z1 = 1. + pow((1. - pow(a, 2.)), 1. / 3) *
                         (pow(1. + a, 1. / 3) + pow(1. - a, 1. / 3));
    double z2 = pow(3. * pow(a, 2.) + pow(z1, 2.), 1. / 2);
    double a_sign = TMath::Sign(1., a);
    return 3 + z2 - pow((3. - z1) * (3. + z1 + 2. * z2), 1. / 2) * a_sign;
}

//______________________________________________________________________________

double CWB::CBCTool::calc_isco_freq(double a) {
    /*
    Calculate the ISCO frequency of a Kerr BH as a function of the Kerr
    parameter

    a : Kerr parameter (numpy array)
    */

    double r_isco = CWB::CBCTool::calc_isco_radius(a);
    double u_isco = pow(r_isco, -0.5);
    double v_isco =
        u_isco *
        pow(1. - a * pow(u_isco, 3.) + pow(a, 2.) * pow(u_isco, 6.), 1. / 3.);
    return pow(v_isco, 3.) / TMath::Pi();
}
//______________________________________________________________________________

double CWB::CBCTool::_final_spin_diff(double a_f,
                                      double eta,
                                      double delta_m,
                                      double S,
                                      double Delta) {
    // Internal function: the final spin is determined by minimizing this
    // function

    // calculate ISCO radius
    double r_isco = CWB::CBCTool::calc_isco_radius(a_f);

    // angular momentum at ISCO -- Eq.(2.8) of Ori, Thorne Phys Rev D 62
    // 124022 (2000)
    double J_isco =
        (3 * pow(r_isco, 1. / 2) - 2 * a_f) * 2. / pow(3 * r_isco, 1. / 2);

    // fitting coefficients - Table X1 of Healy et al Phys Rev D 90, 104004
    // (2014) [forth order fits]
    double L0 = 0.686710;
    double L1 = 0.613247;
    double L2a = -0.145427;
    double L2b = -0.115689;
    double L2c = -0.005254;
    double L2d = 0.801838;
    double L3a = -0.073839;
    double L3b = 0.004759;
    double L3c = -0.078377;
    double L3d = 1.585809;
    double L4a = -0.003050;
    double L4b = -0.002968;
    double L4c = 0.004364;
    double L4d = -0.047204;
    double L4e = -0.053099;
    double L4f = 0.953458;
    double L4g = -0.067998;
    double L4h = 0.001629;
    double L4i = -0.066693;

    double a_f_new =
        pow((4. * eta), 2.) *
            (L0 + L1 * S + L2a * Delta * delta_m + L2b * pow(S, 2.) +
             L2c * pow(Delta, 2.) + L2d * pow(delta_m, 2.) +
             L3a * Delta * S * delta_m + L3b * S * pow(Delta, 2.) +
             L3c * pow(S, 3.) + L3d * S * pow(delta_m, 2.) +
             L4a * Delta * pow(S, 2) * delta_m +
             L4b * pow(Delta, 3.) * delta_m + L4c * pow(Delta, 4.) +
             L4d * pow(S, 4.) + L4e * pow(Delta, 2.) * pow(S, 2.) +
             L4f * pow(delta_m, 4.) + L4g * Delta * pow(delta_m, 3.) +
             L4h * pow(Delta, 2.) * pow(delta_m, 2.) +
             L4i * pow(S, 2.) * pow(delta_m, 2.)) +
        S * (1. + 8. * eta) * pow(delta_m, 4.) +
        eta * J_isco * pow(delta_m, 6.);

    return TMath::Abs(a_f - a_f_new);
}

//______________________________________________________________________________

double CWB::CBCTool::_final_spin_diff(Double_t* x, Double_t* par) {
    double a_f = x[0];
    double eta = par[0];
    double delta_m = par[1];
    double S = par[2];
    double Delta = par[3];
    // Internal function: the final spin is determined by minimizing this
    // function

    // calculate ISCO radius
    double r_isco = CWB::CBCTool::calc_isco_radius(a_f);

    // angular momentum at ISCO -- Eq.(2.8) of Ori, Thorne Phys Rev D 62
    // 124022 (2000)
    double J_isco =
        (3 * pow(r_isco, 1. / 2) - 2 * a_f) * 2. / pow(3 * r_isco, 1. / 2);

    // fitting coefficients - Table X1 of Healy et al Phys Rev D 90, 104004
    // (2014) [forth order fits]
    double L0 = 0.686710;
    double L1 = 0.613247;
    double L2a = -0.145427;
    double L2b = -0.115689;
    double L2c = -0.005254;
    double L2d = 0.801838;
    double L3a = -0.073839;
    double L3b = 0.004759;
    double L3c = -0.078377;
    double L3d = 1.585809;
    double L4a = -0.003050;
    double L4b = -0.002968;
    double L4c = 0.004364;
    double L4d = -0.047204;
    double L4e = -0.053099;
    double L4f = 0.953458;
    double L4g = -0.067998;
    double L4h = 0.001629;
    double L4i = -0.066693;

    double a_f_new =
        pow((4. * eta), 2.) *
            (L0 + L1 * S + L2a * Delta * delta_m + L2b * pow(S, 2.) +
             L2c * pow(Delta, 2.) + L2d * pow(delta_m, 2.) +
             L3a * Delta * S * delta_m + L3b * S * pow(Delta, 2.) +
             L3c * pow(S, 3.) + L3d * S * pow(delta_m, 2.) +
             L4a * Delta * pow(S, 2) * delta_m +
             L4b * pow(Delta, 3.) * delta_m + L4c * pow(Delta, 4.) +
             L4d * pow(S, 4.) + L4e * pow(Delta, 2.) * pow(S, 2.) +
             L4f * pow(delta_m, 4.) + L4g * Delta * pow(delta_m, 3.) +
             L4h * pow(Delta, 2.) * pow(delta_m, 2.) +
             L4i * pow(S, 2.) * pow(delta_m, 2.)) +
        S * (1. + 8. * eta) * pow(delta_m, 4.) +
        eta * J_isco * pow(delta_m, 6.);

    return TMath::Abs(a_f - a_f_new);
}

//______________________________________________________________________________

double CWB::CBCTool::bbh_final_mass_and_spin_non_precessing(double m1,
                                                            double m2,
                                                            double chi1,
                                                            double chi2) {
    /*
    Calculate the mass and spin of the final BH resulting from the
    merger of two black holes with non-precessing spins

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    */

    // binary parameters
    double m = m1 + m2;
    double q = m1 / m2;
    double eta = q / pow((1. + q), 2.);
    double delta_m = (m1 - m2) / m;

    double S1 = chi1 * pow(m1, 2.);     // spin angular momentum 1
    double S2 = chi2 * pow(m2, 2.);     // spin angular momentum 2
    double S = (S1 + S2) / pow(m, 2.);  // symmetric spin (dimensionless --
                                        // called \tilde{S} in the paper)
    double Delta =
        (S2 / m2 - S1 / m1) / m;  // antisymmetric spin (dimensionless --
                                  // called tilde{Delta} in the paper

    //# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    //# compute the final spin
    //# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    //# res = so.minimize_scalar(_final_spin_diff, bounds=(-0.99999,
    // 0.99999), args=(eta, delta_m, S, Delta), method='Bounded', tol=1e-6,
    // options={'maxiter':100, 'disp':False})
    // a_f = res.x

    // x, cov_x = so.leastsq(_final_spin_diff, 0., args=(eta, delta_m, S,
    // Delta))
    TF1 f("spin_diff", CWB::CBCTool::_final_spin_diff, -1.0, 1.0, 4);
    f.SetParameters(eta, delta_m, S, Delta);
    ROOT::Math::WrappedTF1 wf1(f);

    // Create the Integrator
    ROOT::Math::BrentRootFinder brf;

    // Set parameters of the method
    brf.SetFunction(wf1, -1.0, 1.0);
    brf.Solve();

    // cout << brf.Root() << endl;
    double a_f = brf.Root();

    //# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    //# now compute the final mass
    //# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    double r_isco = CWB::CBCTool::calc_isco_radius(a_f);

    //# fitting coefficients - Table X1 of Healy et al Phys Rev D 90, 104004
    //(2014) # [forth order fits]
    double M0 = 0.951507;
    double K1 = -0.051379;
    double K2a = -0.004804;
    double K2b = -0.054522;
    double K2c = -0.000022;
    double K2d = 1.995246;
    double K3a = 0.007064;
    double K3b = -0.017599;
    double K3c = -0.119175;
    double K3d = 0.025000;
    double K4a = -0.068981;
    double K4b = -0.011383;
    double K4c = -0.002284;
    double K4d = -0.165658;
    double K4e = 0.019403;
    double K4f = 2.980990;
    double K4g = 0.020250;
    double K4h = -0.004091;
    double K4i = 0.078441;

    //# binding energy at ISCO -- Eq.(2.7) of Ori, Thorne Phys Rev D 62
    // 124022 (2000)
    double E_isco =
        (1. - 2. / r_isco + a_f / pow(r_isco, 1.5)) /
        pow(1. - 3. / r_isco + 2. * a_f / pow(r_isco, 1.5), 1. / 2.);

    //# final mass -- Eq. (14) of Healy et al Phys Rev D 90, 104004 (2014)
    double mf = pow(4. * eta, 2.) *
                    (M0 + K1 * S + K2a * Delta * delta_m + K2b * pow(S, 2.) +
                     K2c * pow(Delta, 2.) + K2d * pow(delta_m, 2.) +
                     K3a * Delta * S * delta_m + K3b * S * pow(Delta, 2.) +
                     K3c * pow(S, 3.) + K3d * S * pow(delta_m, 2.) +
                     K4a * Delta * pow(S, 2.) * delta_m +
                     K4b * pow(Delta, 3.) * delta_m + K4c * pow(Delta, 4.) +
                     K4d * pow(S, 4.) + K4e * pow(Delta, 2.) * pow(S, 2.) +
                     K4f * pow(delta_m, 4.) + K4g * Delta * pow(delta_m, 3.) +
                     K4h * pow(Delta, 2.) * pow(delta_m, 2.) +
                     K4i * pow(S, 2.) * pow(delta_m, 2.)) +
                (1 + eta * (E_isco + 11.)) * pow(delta_m, 6.);

    return mf * m;
}

double CWB::CBCTool::chip(double m1,
                          double m2,
                          double s1x,
                          double s1y,
                          double s1z,
                          double s2x,
                          double s2y,
                          double s2z) {
    /*
    Calculate the dimensionless precession spin parameter chi_p in source
    frame

    m1, m2: component masses in solar masses
    s{1,2}_{x,y,z}: dimensionless spins of two BHs
    */

    double M = m1 + m2;
    double m1_2 = m1 * m1;
    double m2_2 = m2 * m2;
    double eta = m1 * m2 / (M * M); /* Symmetric mass-ratio */

    /* Aligned spins */
    double chi1_l = s1z; /* Dimensionless aligned spin on BH 1 */
    double chi2_l = s2z; /* Dimensionless aligned spin on BH 2 */

    /* Magnitude of the spin projections in the orbital plane */
    double S1_perp = m1_2 * sqrt(s1x * s1x + s1y * s1y);
    double S2_perp = m2_2 * sqrt(s2x * s2x + s2y * s2y);

    /* From this we can compute chip*/
    double A1 = 2 + (3 * m2) / (2 * m1);
    double A2 = 2 + (3 * m1) / (2 * m2);
    double ASp1 = A1 * S1_perp;
    double ASp2 = A2 * S2_perp;
    double num = (ASp2 > ASp1) ? ASp2 : ASp1;
    double den = (m2 > m1) ? A2 * m2_2 : A1 * m1_2;

    double chip =
        num / den; /*  chip = max(A1 Sp1, A2 Sp2) / (A_i m_i^2) for i index
                      of larger BH (See Eqn. 32 in technical document) */

    return chip;
}
//______________________________________________________________________________

void CWB::CBCTool::AddChip(TString filein, TString treename) {
    // CWB::CBCTool cbcTool;
    TFile* f = new TFile(filein, "update");
    TTree* T = (TTree*)f->Get(treename);
    float mass[2];
    float spin[6];
    float chip;
    TBranch* bchip = T->Branch("chip", &chip, "chip/F");
    T->SetBranchAddress("mass", mass);
    T->SetBranchAddress("spin", spin);
    Long64_t nentries = T->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        T->GetEntry(i);
        chip = CWB::CBCTool::chip(mass[0], mass[1], spin[0], spin[1], spin[2],
                                  spin[3], spin[4], spin[5]);
        // cout <<  mass[0] << " " << mass[1] << " " << spin[0] << " "
        // << spin[1] << " " <<  spin[2] << " " << spin[3] << " " <<
        // spin[4] << " " <<spin[5] << endl;
        bchip->Fill();
    }
    T->Print();
    T->Write();
    delete f;
}

//______________________________________________________________________________

// Creates various distance vs pars plots
// Note : this method is used to generate the CBC report

#define MAXY 50000.0
#define MAXSNR 150.0

void CWB::CBCTool::CreateDistanceParplots(char* sim_file_name,
                                          char* mdc_file_name,
                                          char* netdir,
                                          TString opt = "",
                                          double MINX = 0.0,
                                          double MAXX = 1.0,
                                          double MAXDISTANCE = 5000.,
                                          int NBIN_DIST = 10,
                                          float T_ifar = 0.0,
                                          float T_win = 0.2,
                                          int nIFO = 2) {
    TCanvas* co_canvas3 = new TCanvas("sd3", "SD3", 3, 47, 1000, 802);
    co_canvas3->SetGridx();
    co_canvas3->SetGridy();
    co_canvas3->SetLogy();

    float myifar, netcc[3];
    float rho[2];
    double mytime[6];
    float factor, mydistance, mchirp;
    float mass[2];
    float spin[6];
    float chip;
    float iSNR[nIFO];
    float snr[nIFO];
    float chirp[6];
    float range[2];
    float iota[2];

    TFile* filein = new TFile(sim_file_name);
    TTree* sim = nullptr;
    filein->GetObject("waveburst", sim);
    if (!sim->GetListOfBranches()->FindObject("chip")) {
        cout << "Adding Chi_p branch to wave tree" << endl;
        CWB::CBCTool::AddChip(sim_file_name, "waveburst");
    }
    sim->SetBranchAddress("mass", mass);
    sim->SetBranchAddress("factor", &factor);
    sim->SetBranchAddress("range", range);
    sim->SetBranchAddress("chirp", chirp);
    sim->SetBranchAddress("rho", rho);
    sim->SetBranchAddress("netcc", netcc);
    sim->SetBranchAddress("ifar", &myifar);
    sim->SetBranchAddress("time", mytime);
    sim->SetBranchAddress("spin", spin);
    sim->SetBranchAddress("chip", &chip);
    sim->SetBranchAddress("iSNR", iSNR);
    sim->SetBranchAddress("iota", iota);

    TFile* filein2 = new TFile(mdc_file_name);
    TTree* mdc = nullptr;
    filein2->GetObject("mdc", mdc);
    gROOT->cd();
    if (!mdc->GetListOfBranches()->FindObject("chip")) {
        cout << "Adding Chi_p branch to mdc tree" << endl;
        CWB::CBCTool::AddChip(mdc_file_name, "mdc");
    }

    mdc->SetBranchAddress("time", mytime);
    mdc->SetBranchAddress("mass", mass);
    mdc->SetBranchAddress("factor", &factor);
    mdc->SetBranchAddress("distance", &mydistance);
    mdc->SetBranchAddress("mchirp", &mchirp);
    mdc->SetBranchAddress("spin", spin);
    mdc->SetBranchAddress("chip", &chip);
    mdc->SetBranchAddress("snr", snr);
    mdc->SetBranchAddress("iota", iota);

    int nevts = (int)mdc->GetEntries();
    float CYS = 86400. * 365.25;

    cout << nevts << " injected signals " << sim->GetEntries("ifar>0")
         << " recovered signals" << endl;
    int countv = 0;
    int countvifar = 0;
    int cnt = 0;
    std::vector<double> xi, xr, yi, yr;

    auto inj = new TH2F("Injected snr vs stat inj", "", NBIN_DIST, MINX, MAXX,
                        100, 0.0, MAXSNR);
    // auto SNR= new TH2F("iSNR vs stat", "", 10, MINX, MAXX, 100,0.0,
    // 100.0);
    TMultiGraph* mg = new TMultiGraph();

    if (opt.Contains("chieff")) {
        inj->GetXaxis()->SetTitle("#chi_{eff}");
        mg->GetXaxis()->SetTitle("#chi_{eff}");
    } else if (opt.Contains("chip")) {
        inj->GetXaxis()->SetTitle("#chi_{p}");
        mg->GetXaxis()->SetTitle("#chi_{p}");
    } else if (opt.Contains("chirp")) {
        inj->GetXaxis()->SetTitle("Chirp Mass (M_{#odot})");
        mg->GetXaxis()->SetTitle("Chirp Mass (M_{#odot})");
    } else if (opt.Contains("mtot")) {
        inj->GetXaxis()->SetTitle("Total Mass (M_{#odot})");
        mg->GetXaxis()->SetTitle("Total Mass (M_{#odot})");
    } else if (opt.Contains("eta")) {
        inj->GetXaxis()->SetTitle("#eta, Symmetric Mass Ratio");
        mg->GetXaxis()->SetTitle("#eta, Symmetric Mass Ratio");
    } else if (opt.Contains("iota")) {
        inj->GetXaxis()->SetTitle("Cos(#iota), Inclination");
        mg->GetXaxis()->SetTitle("Cos(#iota), Inclination");
    } else if (opt.Contains("distance")) {
        inj->GetXaxis()->SetTitle("Sensitive Distance (Mpc)");
        mg->GetXaxis()->SetTitle("Injected snr");
        // cout <<"Distance plots" << endl;
    } else {
        cout << "Not a valid option! "
                "opt={\"chip\",\"chirp\",\"eta\",\"mtot\", \"chieff\", "
                " \"iota\"}"
             << endl;
        exit(1);
    }

    // inj->GetYaxis()->SetRangeUser(10., MAXDISTANCE);
    inj->GetYaxis()->SetTitle("Injected snr");
    inj->GetXaxis()->SetTitleOffset(1.3);
    inj->GetYaxis()->SetTitleOffset(1.3);
    inj->GetXaxis()->SetTickLength(0.01);
    inj->GetYaxis()->SetTickLength(0.01);
    inj->GetXaxis()->CenterTitle(kTRUE);
    inj->GetYaxis()->CenterTitle(kTRUE);
    inj->GetXaxis()->SetTitleFont(42);
    inj->GetXaxis()->SetLabelFont(42);
    inj->GetYaxis()->SetTitleFont(42);
    inj->GetYaxis()->SetLabelFont(42);
    inj->SetMarkerStyle(20);
    inj->SetMarkerSize(0.5);
    inj->SetMarkerColor(2);
    inj->SetTitle("");
    TH2F* rec = (TH2F*)inj->Clone("Injected snr vs stat rec");
    // rec->GetYaxis()->SetRangeUser(10., MAXDISTANCE/1000);
    // TH2F *SNRrec = (TH2F *)SNR->Clone("iSNR vs stat rec");
    rec->SetMarkerColor(4);

    mg->GetYaxis()->SetTitle("Sensitive Distance (Mpc)");
    // mg->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
    mg->GetXaxis()->SetTitleOffset(1.3);
    mg->GetYaxis()->SetTitleOffset(1.3);
    mg->GetXaxis()->SetTickLength(0.01);
    mg->GetYaxis()->SetTickLength(0.01);
    mg->GetXaxis()->CenterTitle(kTRUE);
    mg->GetYaxis()->CenterTitle(kTRUE);
    mg->GetXaxis()->SetTitleFont(42);
    mg->GetXaxis()->SetLabelFont(42);
    mg->GetYaxis()->SetTitleFont(42);
    mg->GetYaxis()->SetLabelFont(42);
    mg->SetTitle("");

    // Loop over mdc TTree
    float SNR2, mSNR;
    for (int g = 0; g < (int)mdc->GetEntries(); g++) {
        mdc->GetEntry(g);
        SNR2 = pow(snr[0], 2.0) + pow(snr[1], 2.0);
        for (int i = 2; i < nIFO; i++) {
            SNR2 += pow(snr[i], 2.0);
        }
        mSNR = TMath::Sqrt(SNR2);
        yi.push_back(mydistance);
        if (opt.Contains("chieff")) {
            xi.push_back((spin[2] * mass[0] + spin[5] * mass[1]) /
                         (mass[1] + mass[0]));
        } else if (opt.Contains("chip")) {
            xi.push_back(chip);
        } else if (opt.Contains("chirp")) {
            xi.push_back(mchirp);
        } else if (opt.Contains("mtot")) {
            xi.push_back(mass[0] + mass[1]);
        } else if (opt.Contains("eta")) {
            xi.push_back(mass[0] * mass[1] / pow(mass[0] + mass[1], 2.0));
        } else if (opt.Contains("iota")) {
            xi.push_back(iota[1]);
        } else if (opt.Contains("distance")) {
            xi.push_back(mydistance);
        }

        inj->Fill(xi[xi.size() - 1], mSNR);
        // SNR->Fill(xi[xi.size()-1], mSNR);
    }
    // cout << MINX << " " << MAXX << " " << NBIN_DIST << " " << MAXDISTANCE
    // << endl;
    inj->Draw("p");
    // Loop over sim TTree
    for (int g = 0; g < (int)sim->GetEntries(); g++) {
        sim->GetEntry(g);
        // ifactor = (int)factor - 1;
        // cout << "g=" << g << " trueindex=" <<index1[g] <<" IFAR=" <<
        // myifar/CYS << " RHO1=" << rho[1] <<endl;
        if (myifar <= T_ifar * CYS) {
            countvifar++;
            // cout << g << " " <<index1[g] <<" " << myifar/CYS << "
            // ";
            continue;
        }

        if ((mytime[0] - mytime[nIFO]) < -T_win ||
            (mytime[0] - mytime[nIFO]) > 2 * T_win) {
            countv++;
            continue;
        }  // NOT checking for detector 1 and 2: very small bias...
        SNR2 = 0.0;
        for (int i = 0; i < nIFO; i++) {
            SNR2 += iSNR[i];
        }
        mSNR = TMath::Sqrt(SNR2);
        yr.push_back(range[1]);
        if (opt.Contains("chieff")) {
            xr.push_back((spin[2] * mass[0] + spin[5] * mass[1]) /
                         (mass[1] + mass[0]));
            cnt++;
        } else if (opt.Contains("chip")) {
            xr.push_back(chip);
            cnt++;
        } else if (opt.Contains("chirp")) {
            xr.push_back(chirp[0]);
            cnt++;
        } else if (opt.Contains("mtot")) {
            xr.push_back(mass[1] + mass[0]);
            cnt++;
        } else if (opt.Contains("eta")) {
            xr.push_back(mass[0] * mass[1] / pow(mass[0] + mass[1], 2.0));
            cnt++;
        } else if (opt.Contains("iota")) {
            xr.push_back(iota[1]);
            cnt++;
        } else if (opt.Contains("distance")) {
            xr.push_back(range[1]);
            cnt++;
        }
        rec->Fill(xr[xr.size() - 1], mSNR);
    }

    cout << endl;
    // cout << countvifar << " events vetoed by T_ifar : " << T_ifar <<
    // endl;
    // cout << countv << " events vetoed by T_win" << endl;
    // cout << rec->GetEntries() << " events selected" << endl;

    char lab[1024];
    char fname[1024];
    char fname2[1024];
    char fname3[1024];

    sprintf(fname, "%s/iSNR_vs_%s.eps", netdir, opt.Data());
    sprintf(fname3, "%s/Distance_vs_%s.eps", netdir, opt.Data());
    sprintf(fname2, "%s/%s_distribution.png", netdir, opt.Data());

    // D_Chi_rec->GetYaxis()->SetRangeUser(10.,3*MAXDISTANCE);
    // co_canvas3->SetLogx();
    inj->Draw("p");
    rec->Draw("p same");
    auto leg_D = new TLegend(0.6, 0.1, 0.9, 0.25, "", "brNDC");
    sprintf(lab, "Injections: %i", (int)mdc->GetEntries());
    leg_D->AddEntry("", lab, "a");
    sprintf(lab, "found: %i", cnt);
    leg_D->AddEntry(rec, lab, "p");
    sprintf(lab, "missed: %i", (int)mdc->GetEntries() - cnt);
    leg_D->AddEntry(inj, lab, "p");
    leg_D->SetFillColor(0);
    leg_D->SetFillColorAlpha(0, 0.9);
    leg_D->Draw();
    co_canvas3->SaveAs(fname);
    co_canvas3->SetLogx(0);
    // Test new scatter plot
    TGraph* rec_gr = new TGraph(xr.size(), &xr[0], &yr[0]);
    TGraph* inj_gr = new TGraph(xi.size(), &xi[0], &yi[0]);
    TH1D* injx;
    TH1D* recx;
    TH1D* snrx;
    if (!opt.Contains("distance")) {
        inj_gr->SetMarkerColor(2);
        inj_gr->SetMarkerStyle(20);
        inj_gr->SetMarkerSize(0.5);

        rec_gr->SetMarkerColor(4);
        rec_gr->SetMarkerStyle(20);
        rec_gr->SetMarkerSize(0.5);

        mg->Add(inj_gr);
        mg->Add(rec_gr);
        // mg->Draw("aple3");
        co_canvas3->Clear();
        mg->GetXaxis()->SetLimits(MINX, MAXX);
        mg->GetYaxis()->SetRangeUser(10., MAXDISTANCE);
        mg->Draw("ap");
        // rec_gr->Draw("ap");
        leg_D->Draw();
        co_canvas3->SaveAs(fname3);
        injx = (TH1D*)inj->ProjectionX();
        recx = (TH1D*)rec->ProjectionX();
        snrx = (TH1D*)inj->ProfileX();
    } else {
        injx = (TH1D*)inj->ProjectionX();
        recx = (TH1D*)rec->ProjectionX();
        snrx = (TH1D*)inj->ProfileX();
    }

    injx->Sumw2();
    // injx->SetFillColor(kRed);
    injx->SetFillColorAlpha(kRed, 0.3);
    injx->SetFillStyle(3004);
    recx->Sumw2();
    recx->SetFillColorAlpha(kBlue, 0.99);
    recx->SetFillStyle(3001);

    co_canvas3->Clear();
    co_canvas3->SetLogx(0);

    // Ratio plot
    auto rp1 = new TRatioPlot(recx, injx, "divsym");
    rp1->SetH1DrawOpt("PEBAR");
    rp1->SetH2DrawOpt("PEBAR");
    rp1->SetGraphDrawOpt("PB");
    rp1->SetSeparationMargin(0.01);

    // rp1->SetConfidenceIntervalColors();
    rp1->Draw("PEBAR");
    rp1->SetSplitFraction(0.5);
    // if (!opt.Contains("distance")){
    rp1->GetUpperRefYaxis()->SetRangeUser(1., MAXY);
    // rp1->GetLowerRefGraph()->SetMinimum(0.0);
    // rp1->GetLowerRefGraph()->SetMaximum(1.0);
    /* } else {
                  //co_canvas3->SetLog();
                  rp1->GetUpperRefYaxis()->SetRangeUser(1.,MAXY);
                  //rp1->GetUpperRefXaxis()->SetRangeUser(1.,2000.);
                  rp1->GetLowerRefGraph()->SetMinimum(0.0);
        rp1->GetLowerRefGraph()->SetMaximum(1.0);
                  //rp1->GetLowYaxis()->SetNdivisions(510);
    }*/

    rp1->GetLowerRefYaxis()->SetTitle("ratio");
    rp1->GetLowerRefYaxis()->CenterTitle(kTRUE);
    rp1->GetUpperRefYaxis()->SetTitle("entries");
    rp1->GetUpperRefYaxis()->CenterTitle(kTRUE);
    // rp1->GetUpperRefYaxis()->SetLineColor(kBlue);
    rp1->GetLowerRefYaxis()->SetLabelColor(kBlue);
    rp1->GetLowerRefXaxis()->CenterTitle(kTRUE);
    rp1->GetLowerRefXaxis()->SetTitleOffset(1.3);
    rp1->GetLowerRefGraph()->SetFillColor(kBlue);
    rp1->GetLowerRefGraph()->SetFillStyle(3001);
    // leg_D2->Draw();

    auto leg_D3 = new TLegend(0.6, 0.8455, 0.8965, 0.94555, "", "brNDC");
    // sprintf(lab, "Injections: %i", (int)mdc->GetEntries());
    // leg_D3->AddEntry("", lab, "a");
    sprintf(lab, "Injected: %i", (int)mdc->GetEntries());
    leg_D3->AddEntry(inj, lab, "p");
    sprintf(lab, "Found: %i", cnt);
    leg_D3->AddEntry(rec, lab, "p");
    // leg_D3->SetFillColor(0);
    leg_D3->SetFillColorAlpha(0, 0.9);
    leg_D3->Draw();
    // TGraph *g = rp1->GetLowerRefGraph();
    // rp11 = TRatioPlot(h1,h0)
    // rp1->GetUpperPad()->cd()
    // h2->Draw("same")
    rp1->GetLowerPad()->cd();
    // TPad* mp = (TPad*)rp1->GetLowerPad();

    Float_t rightmax = 1.2 * snrx->GetMaximum();
    Float_t rightmin = 0.8 * snrx->GetMinimum();
    double ymax = rp1->GetLowerRefYaxis()->GetXmax();
    double ymin = rp1->GetLowerRefYaxis()->GetXmin();
    double xmax = rp1->GetLowerRefXaxis()->GetXmax();
    Float_t scale = (ymax) / (rightmax);
    // cout << "rightmax: " << rightmax << " Uymax: "<<gPad->GetUymax() << "
    // scale: " << scale << endl; g->Draw(); rp1->GetLowerPad()->cd();
    snrx->Scale(scale);
    snrx->SetLineColor(kGreen + 2);
    snrx->SetMarkerColor(kGreen + 2);
    // snrx->SetFillColor(kGreen+2);
    snrx->SetFillColorAlpha(kGreen + 2, 0.4);
    snrx->SetFillStyle(3001);
    snrx->SetMarkerStyle(20);
    snrx->SetMarkerSize(0.5);

    snrx->Draw("SAME PEBAR");
    //	TGaxis*myaxis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
    // gPad->GetUxmax(),gPad->GetUymax(),	0,rightmax,510,"+L");
    TGaxis* myaxis =
        new TGaxis(xmax, ymin, xmax, ymax, rightmin, rightmax, 510, "+SL");
    myaxis->SetLineColor(kGreen + 2);
    myaxis->SetLabelColor(kGreen + 2);
    myaxis->SetTitle("Average Injected SNR");
    myaxis->SetTitleOffset(0.7);
    myaxis->SetTickLength(0.01);
    myaxis->CenterTitle(kTRUE);
    myaxis->SetTitleFont(42);
    myaxis->SetLabelFont(42);
    myaxis->SetLabelSize(0.06);
    myaxis->SetTitleSize(0.06);

    myaxis->Draw();

    auto leg_D4 = new TLegend(0.6, 0.75, 0.9, 0.975, "", "brNDC");
    // sprintf(lab, "Injections: %i", (int)mdc->GetEntries());
    // leg_D3->AddEntry("", lab, "a");
    sprintf(lab, "Ratio (i.e. Found / Injected)");
    leg_D4->AddEntry(recx, lab, "p");
    sprintf(lab, "Average Injected SNR");
    leg_D4->AddEntry(snrx, lab, "p");
    // leg_D3->SetFillColor(0);
    leg_D4->SetFillColorAlpha(0, 0.9);
    // leg_D4->Draw();
    // rp1->GetLowYaxis()->SetNdivisions(505);
    co_canvas3->Update();

    if (!opt.Contains("I")) {
        rp1->GetLowerPad()->SetEditable(kFALSE);
        co_canvas3->SaveAs(fname2);
        filein->Close();
        filein2->Close();
        delete filein, filein2;
        delete co_canvas3;
        delete inj, rec, injx, recx, snrx, leg_D, leg_D3, leg_D4, rp1, myaxis;
        xi.clear(), xr.clear(), yi.clear(), yr.clear();
        delete inj_gr, rec_gr, mg;
    }

    // exit(0);
}

//______________________________________________________________________________

// Create Graphs of Radius vs IFAR
// Note : this method is used to generate the CBC report

#define MINMAXIFAR 4278.

TGraphErrors* CWB::CBCTool::CreateGraphRadiusIFAR(char* sim_file_name,
                                                  char* mdc_file_name,
                                                  TString SEL,
                                                  float shell_volume,
                                                  Color_t color = kBlue,
                                                  TString opt = "default",
                                                  double liveTot = 1e6,
                                                  float T_ifar = 0.0,
                                                  float T_win = 0.2,
                                                  int TRIALS = 1,
                                                  int nIFO = 2,
                                                  float VT = 1.0,
                                                  float Tscale = 1.0) {
    float myifar, netcc[3];
    float rho[2];
    double mytime[6];
    float factor, distance, mchirp;
    float mass[2];
    float spin[6];
    // float chi[3];

    float chirp[6];
    float range[2];

    CWB::CBCTool cbcTool;

    TFile* filein = new TFile(sim_file_name);
    TTree* sim_org = nullptr;
    filein->GetObject("waveburst", sim_org);
    gROOT->cd();
    if (!sim_org->GetListOfBranches()->FindObject("chip")) {
        cout << "Adding Chi_p branch to wave tree" << endl;
        cbcTool.AddChip(sim_file_name, "waveburst");
    }
    TTree* sim = sim_org->CopyTree(SEL);
    sim->SetBranchAddress("mass", mass);
    sim->SetBranchAddress("factor", &factor);
    sim->SetBranchAddress("range", range);
    sim->SetBranchAddress("chirp", chirp);
    sim->SetBranchAddress("rho", rho);
    sim->SetBranchAddress("netcc", netcc);
    sim->SetBranchAddress("ifar", &myifar);
    sim->SetBranchAddress("time", mytime);
    sim->SetBranchAddress("spin", spin);

    TFile* filein2 = new TFile(mdc_file_name);
    TTree* mdc_org = nullptr;
    filein2->GetObject("mdc", mdc_org);
    gROOT->cd();
    SEL.ReplaceAll("chirp[0]", "mchirp");
    if (!mdc_org->GetListOfBranches()->FindObject("chip")) {
        cout << "Adding Chi_p branch to mdc tree" << endl;
        cbcTool.AddChip(mdc_file_name, "mdc");
    }
    TTree* mdc = mdc_org->CopyTree(SEL);
    mdc->SetBranchAddress("time", mytime);
    mdc->SetBranchAddress("mass", mass);
    mdc->SetBranchAddress("factor", &factor);
    mdc->SetBranchAddress("distance", &distance);
    mdc->SetBranchAddress("mchirp", &mchirp);
    mdc->SetBranchAddress("spin", spin);

    double dV = 0.0;
    std::vector<double> vdv;
    std::vector<double> vifar;
    std::vector<double> vrho;
    std::vector<double> vV;
    std::vector<double> veV;
    std::vector<double> vR;
    std::vector<double> veR;
    std::vector<double> vsifar;
    std::vector<double> vseifar;
    std::vector<double> vVr;
    std::vector<double> veVr;
    std::vector<double> vRr;
    std::vector<double> veRr;
    std::vector<double> vsrho;
    TGraphErrors* co_gr = nullptr;

    int nevts = (int)mdc->GetEntries();
    if (nevts == 0) {
        cout << "No events left after cut!" << endl;
        vR.push_back(0.0);
        vsifar.push_back(0.0);
        co_gr = new TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, 0);
        return co_gr;
    }
    float shell_volume_per_injection = shell_volume / nevts;
    int ifactor;
    float maxIFAR = 0.0;
    float CYS = 86400. * 365.25;
    if (sim->GetListOfBranches()->FindObject("ifar")) {
        maxIFAR = TMath::CeilNint(sim->GetMaximum("ifar") / CYS);
        // cout << "Maximum empirically estimated IFAR : " << maxIFAR <<
        // " [years]" << endl;
    } else {
        cout << "Missing ifar branch: either use cbc_plots or add it "
                "to wave tree."
             << endl;
        exit(1);
    }
    if (opt.Contains("rho")) {
        sim->BuildIndex("0", "rho[1]*10000");  // BEWARE rho[1] could be
                                               // search depandent!
    } else {
        sim->BuildIndex(
            "ifar",
            "rho[1]*1000");  // BEWARE rho[1] could be search depandent!
    }
    TTreeIndex* I1 = (TTreeIndex*)sim->GetTreeIndex();  // get the tree index
    Long64_t* index1 =
        I1->GetIndex();  // create an array of entries in sorted order

    // cout << nevts << " injected signals " << sim->GetEntries("ifar>0") <<
    // " recovered signals" <<endl;
    int countv = 0;
    int countvifar = 0;
    bool Error = false;
    double mpcTscale = 1.0;
    for (int g = 0; g < (int)sim->GetEntries(); g++) {
        sim->GetEntry(index1[g]);
        ifactor = (int)factor - 1;
        // cout << "g=" << g << " trueindex=" <<index1[g] <<" IFAR=" <<
        // myifar/CYS << " RHO1=" << rho[1] <<endl;
        if (myifar <= T_ifar * CYS) {
            countvifar++;
            // cout << g << " " <<index1[g] <<" " << myifar/CYS << "
            // ";
            continue;
        }

        if ((mytime[0] - mytime[nIFO]) < -T_win ||
            (mytime[0] - mytime[nIFO]) > 2 * T_win) {
            countv++;
            continue;
        }  // NOT checking for detector 1 and 2: very small bias...

        if (opt.Contains("DDistrVolume")) {
            dV = shell_volume_per_injection;
        } else if (opt.Contains("DDistrUniform")) {
            dV = pow(range[1], 2) * shell_volume_per_injection;
        } else if (opt.Contains("DDistrChirpMass")) {
            dV = pow(range[1], 2) * shell_volume_per_injection *
                 pow(chirp[0] / 1.22, 5. / 6.);
        } else if (opt.Contains("Redshift")) {
            dV = VT / (double)nevts;
            mpcTscale = Tscale * 1e-9;  // Conversion from Gpc^3 to Mpc^3
        } else {
            if (!Error) {
                cout << "No defined distance distribution? "
                        "WARNING: Assuming uniform in volume"
                     << endl;
                Error = 1;
            }
            dV = shell_volume_per_injection;
        }
        // vdv.push_back(dV+internal_volume);
        vdv.push_back(dV);
        vifar.push_back(myifar);
        vrho.push_back(rho[1]);
    }

    // cout << endl;
    // cout << countvifar << " events vetoed by T_ifar : " << T_ifar <<
    // endl; cout << countv << " events vetoed by T_win" << endl; cout <<
    // vdv.size() << " events selected" << endl;

    // cout << "DEB2" << endl;
    vV.push_back(vdv[vdv.size() - 1]);
    veV.push_back(pow(vdv[vdv.size() - 1], 2));
    if (vifar[vdv.size() - 1] < MINMAXIFAR * CYS) {
        vsifar.push_back(vifar[vdv.size() - 1]);
    } else {
        vsifar.push_back(MINMAXIFAR * CYS);
    }
    vVr.push_back(vdv[vdv.size() - 1]);
    veVr.push_back(pow(vdv[vdv.size() - 1], 2));
    vsrho.push_back(vrho[vdv.size() - 1]);
    // cout << "vifar[vdv.size()-1] = " << vifar[vdv.size()-1] << endl;
    // cout << "vsifar[0] = " << vsifar[0] << endl;
    // break;
    int mcount_ifar = 0;
    int mcount_rho = 0;
    for (int i = vdv.size() - 1; i >= 0; i--) {
        if (vifar[i] == 0) {
            continue;
        }
        if (vifar[i] > vsifar[0]) {
            vV[0] += vdv[i];
            veV[0] += pow(vdv[i], 2);
        } else if (vifar[i] == vsifar[mcount_ifar]) {
            vV[mcount_ifar] += vdv[i];
            veV[mcount_ifar] += pow(vdv[i], 2);
        } else {
            vsifar.push_back(vifar[i]);
            vseifar.push_back(TMath::Sqrt(TMath::Nint(liveTot * vifar[i])));
            vV.push_back(vV[mcount_ifar] + vdv[i]);
            veV.push_back(veV[mcount_ifar] + pow(vdv[i], 2));
            mcount_ifar++;
        }
        // cout << i << " " << vsifar[mcount_ifar] << " " <<
        // vV[mcount_ifar] << endl;
        if (vrho[i] == vsrho[mcount_rho]) {
            vVr[mcount_rho] += vdv[i];
            veVr[mcount_rho] += pow(vdv[i], 2);
        } else {
            vsrho.push_back(vrho[i]);
            vVr.push_back(vVr[mcount_rho] + vdv[i]);
            veVr.push_back(veVr[mcount_rho] + pow(vdv[i], 2));
            mcount_rho++;
        }
    }
    // cout << "Length of ifar/volume vector: " << vV.size() << endl;
    // cout << "Length of rho/volume vector: " << vVr.size() << endl;
    //	float
    for (int i = 0; i < (int)vV.size(); i++) {
        veV[i] = TMath::Sqrt(veV[i]);
        vsifar[i] /= (TRIALS * CYS);
        vR.push_back(pow(3. / (4. * TMath::Pi() * mpcTscale) * vV[i], 1. / 3.));
        veR.push_back(pow(3. / (4 * TMath::Pi() * mpcTscale), 1. / 3.) * 1. /
                      3 * pow(vV[i], -2. / 3.) * veV[i]);
        //	cout << i << " " <<  vsifar[i] << " " << vR[i] << " "
        //<<vV[i] <<endl;
    }

    for (int i = 0; i < (int)vVr.size(); i++) {
        veVr[i] = TMath::Sqrt(veVr[i]);
        vRr.push_back(
            pow(3. / (4. * TMath::Pi() * mpcTscale) * vVr[i], 1. / 3.));
        veRr.push_back(pow(3. / (4 * TMath::Pi() * mpcTscale), 1. / 3.) * 1. /
                       3 * pow(vVr[i], -2. / 3.) * veVr[i]);
        // cout <<  vsrho[i] << " ";
    }
    // cout << endl;

    // std::vector<TGraphErrors> v_gr(2);
    // TGraphErrors *co_gr = new TGraphErrors[2];
    // TGraphErrors co_gr;
    // TGraphErrors* co_gr = new TGraphErrors(vR.size(), &vsifar[0], &vR[0],
    // 0, &veR[0]);
    if (opt.Contains("rho")) {
        co_gr = new TGraphErrors(vRr.size(), &vsrho[0], &vRr[0], 0, &veRr[0]);
    } else {
        co_gr = new TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, &veR[0]);
    }
    // co_gr[0] = TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, &veR[0]);
    // co_gr->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    // co_gr->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
    // co_gr->GetXaxis()->SetLimits(MINIFAR,MAXIFAR);
    // co_gr->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    co_gr->GetXaxis()->SetTitleOffset(1.3);
    co_gr->GetYaxis()->SetTitleOffset(1.25);
    co_gr->GetXaxis()->SetTickLength(0.01);
    co_gr->GetYaxis()->SetTickLength(0.01);
    co_gr->GetXaxis()->CenterTitle(kTRUE);
    co_gr->GetYaxis()->CenterTitle(kTRUE);
    co_gr->GetXaxis()->SetTitleFont(42);
    co_gr->GetXaxis()->SetLabelFont(42);
    co_gr->GetYaxis()->SetTitleFont(42);
    co_gr->GetYaxis()->SetLabelFont(42);
    co_gr->SetMarkerStyle(20);
    co_gr->SetMarkerSize(1.0);
    co_gr->SetMarkerColor(1);
    co_gr->SetLineColor(color);
    co_gr->SetLineWidth(3);
    co_gr->SetTitle("");
    co_gr->SetFillColor(color);
    co_gr->SetFillStyle(3001);
    // co_gr->Draw("aple3");

    filein->Close();
    filein2->Close();
    delete filein, filein2;
    delete mdc, mdc_org, sim, sim_org;
    vifar.clear();
    vdv.clear();
    vrho.clear();
    vsifar.clear();
    vseifar.clear();
    vsrho.clear();
    vR.clear();
    veR.clear();
    veRr.clear();
    vV.clear();
    veV.clear();
    vRr.clear();
    vVr.clear();
    veVr.clear();
    return co_gr;

    // exit(0);
}

//______________________________________________________________________________

/*
void CWB::CBCTool::doEffectiveRadiusChiPlot(){


          int spin_mtot_bins = 0;
          double V0_spin_mtot = 0.0;
         for(int j=0; j<NBINS_MTOT+1; j++){
            for(int k=0; k<NBINS_SPIN+1; k++){

             //  volume_first_shell[j][k] =
efficiency_first_shell->GetBinContent(j+1,k+1);
              // if(factor_events_rec->GetBinContent(j+1,k+1) != 0.) {
                //	    error_volume_first_shell[j][k]
= 1./TMath::Sqrt(factor_events_rec->GetBinContent(j+1,k+1));
                //	    massbins++;
                //}
                   if(spin_mtot_volume[j][k] != 0.) {
                          spin_mtot_volume[j][k] =
shell_volume*spin_mtot_volume[j][k]; ///  Warning: neglecting the internal
volume...  + volume_internal_sphere*volume_first_shell[j][k]; V0_spin_mtot +=
spin_mtot_volume[j][k]; error_spin_mtot_volume[j][k] =
shell_volume*TMath::Sqrt(error_spin_mtot_volume[j][k]); ///Warning: neglecting
the internal volume...+
volume_internal_sphere*volume_first_shell[j][k]*error_volume_first_shell[j][k];

                          spin_mtot_radius[j][k] =
pow(3.*spin_mtot_volume[j][k]/(4*TMath::Pi()),1./3);

                         error_spin_mtot_radius[j][k] =
(1./3)*pow(3./(4*TMath::Pi()),1./3)*pow(1./pow(spin_mtot_volume[j][k],2),1./3)*error_spin_mtot_volume[j][k];
                         spin_mtot_bins++;
                   }
                   cout<<j<< " "<<k<< " "<<spin_mtot_volume[j][k]<<"
"<<spin_mtot_radius[j][k]<<endl;
            }
         }
         cout << "Average Spin-Mtot Volume at threshold V0 =
"<<V0_spin_mtot/spin_mtot_bins <<endl; c1->Clear();

         TH2F* h_spin_mtot_radius = new
TH2F("h_spin_mtot_radius","",NBINS_SPIN,MINCHI,MAXCHI,NBINS_MTOT,minMtot,maxMtot);
        //h_spin_mtot_radius->GetXaxis()->SetRangeUser(MIN_plot_mass1,MAX_plot_mass1);
        //h_spin_mtot_radius->GetYaxis()->SetRangeUser(MIN_plot_mass2,MAX_plot_mass2);
        h_spin_mtot_radius->GetXaxis()->SetTitle("#chi_{z}");
        h_spin_mtot_radius->GetYaxis()->SetTitle("Total Mass (M_{#odot})");
        h_spin_mtot_radius->GetXaxis()->SetTitleOffset(1.3);
        h_spin_mtot_radius->GetYaxis()->SetTitleOffset(1.25);
        h_spin_mtot_radius->GetXaxis()->CenterTitle(kTRUE);
        h_spin_mtot_radius->GetYaxis()->CenterTitle(kTRUE);
        h_spin_mtot_radius->GetXaxis()->SetNdivisions(410);
        h_spin_mtot_radius->GetYaxis()->SetNdivisions(410);
        h_spin_mtot_radius->GetXaxis()->SetTickLength(0.01);
        h_spin_mtot_radius->GetYaxis()->SetTickLength(0.01);
        h_spin_mtot_radius->GetZaxis()->SetTickLength(0.01);
        h_spin_mtot_radius->GetXaxis()->SetTitleFont(42);
        h_spin_mtot_radius->GetXaxis()->SetLabelFont(42);
        h_spin_mtot_radius->GetYaxis()->SetTitleFont(42);
        h_spin_mtot_radius->GetYaxis()->SetLabelFont(42);
        h_spin_mtot_radius->GetZaxis()->SetLabelFont(42);
        h_spin_mtot_radius->GetZaxis()->SetLabelSize(0.03);
        h_spin_mtot_radius->SetTitle("");


         for(int i=1; i<=NBINS_MTOT+1; i++){
            for(int j=1; j<=NBINS_SPIN+1; j++){
              h_spin_mtot_radius->SetBinContent(j,i,spin_mtot_radius[i-1][j-1]);
              h_spin_mtot_radius->SetBinError(j,i,error_spin_mtot_radius[i-1][j-1]);
                  cout<<j<< " "<<i<< "
"<<h_spin_mtot_radius->GetBinContent(j,i)<<endl;
            }
         }


        h_spin_mtot_radius->SetContour(NCont);
        h_spin_mtot_radius->SetEntries(1);  // This option needs to be enabled
when filling 2D histogram with SetBinContent h_spin_mtot_radius->Draw("colz text
colsize=2");  // Option to write error associated to the bin content
        //h_spin_mtot_radius->Draw("colz text");
        h_spin_mtot_radius->GetZaxis()->SetRangeUser(0,MAX_EFFECTIVE_RADIUS/2.);

         TExec *ex2 = new TExec("ex2","gStyle->SetPaintTextFormat(\".0f\");");
         ex2->Draw();


         //char radius_title[256];
         sprintf(radius_title,"Effective radius Mtot vs #chi_{z} (Mpc, %.1f <
#chi_{z} < %.1f)",MINCHI,MAXCHI);



         TPaveText *p_radius = new
TPaveText(0.325301,0.926166,0.767068,0.997409,"blNDC");
         p_radius->SetBorderSize(0);
         p_radius->SetFillColor(0);
         p_radius->SetTextColor(1);
         p_radius->SetTextFont(32);
         p_radius->SetTextSize(0.045);
         TText *text = p_radius->AddText(radius_title);
         p_radius->Draw();

         sprintf(fname,"%s/Effective_radius_chi.png",netdir);

         c1->Update();
         c1->SaveAs(fname);


}

*/
