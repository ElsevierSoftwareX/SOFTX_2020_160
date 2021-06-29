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


//
// Creates various distance vs pars plots
// Note : this macro is used to generate the CBC report
// Author : Francesco Salemi

#include <TComplex.h>
#include <TRandom.h>
#include <TStyle.h>
#include <fstream>
#include <iostream>
#include "Math/Rotation3D.h"
#include "Math/Vector3Dfwd.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TRatioPlot.h"
#include "TRotation.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TVector3.h"

#define MAXY 50000.0
#define MAXSNR 150.0

// TGraphErrors* CreateGraphRadiusIFAR (char* sim_file_name, char*
// mdc_file_name, float shell_volume, double liveTot=1e6, float T_ifar=0.0, float
// T_win=0.2, int TRIALS = 1, TString SEL="", Color_t color=kBlue, TString
// opt="default"){
void CreateDistanceParplots(char* sim_file_name,
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

		CWB::CBCTool cbcTool;

        TFile* filein = new TFile(sim_file_name);
        TTree* sim = nullptr;
        filein->GetObject("waveburst", sim);
        if (!sim->GetListOfBranches()->FindObject("chip")) {
                cout << "Adding Chi_p branch to wave tree" << endl;
                cbcTool.AddChip(sim_file_name, "waveburst");
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
                cbcTool.AddChip(mdc_file_name, "mdc");
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

        auto inj = new TH2F("Injected snr vs stat inj", "", NBIN_DIST, MINX,
                            MAXX, 100, 0.0, MAXSNR);
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
                        "opt={\"chip\",\"chirp\",\"eta\",\"mtot\", \"chieff\",  \"iota\"}"
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
        
        //Loop over mdc TTree
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
                        xi.push_back(mass[0] * mass[1] /
                                     pow(mass[0] + mass[1], 2.0));
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
                        xr.push_back(mass[0] * mass[1] /
                                     pow(mass[0] + mass[1], 2.0));
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
       // cout << countvifar << " events vetoed by T_ifar : " << T_ifar << endl;
        //cout << countv << " events vetoed by T_win" << endl;
        //cout << rec->GetEntries() << " events selected" << endl;

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
        
        //Ratio plot
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
        //rp1->GetLowerRefGraph()->SetMinimum(0.0);
        //rp1->GetLowerRefGraph()->SetMaximum(1.0);
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
        Float_t scale = ( ymax ) / (rightmax );
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
        //gPad->GetUxmax(),gPad->GetUymax(),	0,rightmax,510,"+L");
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
                delete inj, rec, injx, recx, snrx, leg_D, leg_D3, leg_D4, rp1,
                    myaxis;
                xi.clear(), xr.clear(), yi.clear(), yr.clear();
                delete inj_gr, rec_gr, mg;
        }

        // exit(0);
}
