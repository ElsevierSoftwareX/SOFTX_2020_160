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
// Draw Radius vs IFAR
// Note : this macro is used to generate the CBC report
// Author : Francesco Salemi

#define MINIFAR 1
#define MAXIFAR 10000.
#define MINRADIUS 0
#define MAXRADIUS 1000
#define MINRHO 5
#define MAXRHO 100

// This line is to load the library
// R__LOAD_LIBRARY($HOME_CWB/macros/cbcTool.CreateGraphRadiusIFAR_C.so)

void DrawRadiusIFARplots(char* sim_file_name,
                         char* mdc_file_name,
                         float shell_volume,
                         TString opt) {
    TCanvas* co_canvas = new TCanvas("sd2", "SD2", 3, 47, 1000, 802);
    co_canvas->SetGridx();
    co_canvas->SetGridy();
    co_canvas->SetLogx();

    std::array<float, 5> BBH_MTotBins = {4.0, 27.25, 51.50, 75.75, 100.0};
    std::array<float, 5> IMBHB_MTotBins = {100.0, 200.0, 400.0, 600.0, 800.0};
    std::array<float, 5> MTotBins = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<float, 5> MChirpBins = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<float, 5> BBH_MChirpBins = {1.74, 8.07, 14.92, 21.77, 100.0};
    std::array<float, 5> IMBHB_MChirpBins = {
        1.74, 100.0, 149.2, 217.7, 500.0};  // TODO find more rational binning
    std::array<float, 4> ChieffBins = {-1.0, -0.4, 0.4, 1.0};
    std::array<float, 4> ChipBins = {0.0, 0.3333, 0.6666, 1.0};
    std::array<float, 4> CosiBins = {-1.0, -0.86602540, 0.86602540, 1.0};
    TString CosiBinsT[4] = {"-1", "-#sqrt{3}/2", "#sqrt{3}/2", "1"};
    std::array<float, 5> EtaBins = {0.25, 0.24, 0.22, 0.18, 0.0};

    // Colors array
    Color_t col_arr[4] = {kBlue, kCyan, kGreen - 9, kOrange};
    Color_t col_arr2[3] = {kViolet - 1, kGreen - 9, kRed};

    if (opt.Contains("IMBHB")) {
        MTotBins = IMBHB_MTotBins;
        MChirpBins = IMBHB_MChirpBins;
    } else {
        MTotBins = BBH_MTotBins;
        MChirpBins = BBH_MChirpBins;
    }
    CWB::Toolbox TB;
    CWB::CBCTool cbcTool;

    //========== Sensitive Distance vs IFAR for varius Mtot bins
    //================
    cout << "Sensitive Distance vs IFAR ROC plot: Mtot" << endl;
    TMultiGraph* mg = new TMultiGraph();
    TString sel;
    int bin_count = 0;
    // TGraphErrors *co_gr = new TGraphErrors[4];
    std::vector<TGraphErrors*> co_gr(4);
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9, "", "brNDC");

    while (bin_count < 4) {
        sel.Form("(mass[0]+mass[1] >= %f) && (mass[0]+mass[1] < %f)",
                 MTotBins[bin_count], MTotBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume, col_arr[bin_count],
            opt, liveTot, T_ifar, T_win, TRIALS, nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg->Add(co_gr[bin_count]);
            sel.Form("M_{total} #in [%3.2f, %3.2f] M_{#odot}",
                     MTotBins[bin_count], MTotBins[bin_count + 1]);
            leg->AddEntry(co_gr[bin_count], sel, "l");
        }
        bin_count++;
    }
    mg->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
    mg->GetXaxis()->SetLimits(MINIFAR, MAXIFAR);
    // mg->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mg->GetXaxis()->SetTitleOffset(1.3);
    mg->GetYaxis()->SetTitleOffset(1.25);
    mg->GetXaxis()->SetTickLength(0.01);
    mg->GetYaxis()->SetTickLength(0.01);
    mg->GetXaxis()->CenterTitle(kTRUE);
    mg->GetYaxis()->CenterTitle(kTRUE);
    mg->GetXaxis()->SetTitleFont(42);
    mg->GetXaxis()->SetLabelFont(42);
    mg->GetYaxis()->SetTitleFont(42);
    mg->GetYaxis()->SetLabelFont(42);

    mg->Draw("aple3");

    leg->SetFillColorAlpha(0, 0.7);
    leg->Draw();

    char fname[1024];
    sprintf(fname, "%s/ROC_IFAR_Mtot.png", netdir);
    co_canvas->Update();
    co_canvas->SaveAs(fname);

    //========== Sensitive Distance vs rho[1] for varius Mtot bins
    //================
    cout << "Sensitive Distance vs rho[1] ROC plot: Mtot" << endl;
    TMultiGraph* mgr = new TMultiGraph();
    bin_count = 0;
    while (bin_count < 4) {
        sel.Form("(mass[0]+mass[1] >= %f) && (mass[0]+mass[1] < %f)",
                 MTotBins[bin_count], MTotBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume, col_arr[bin_count],
            opt + " rho", liveTot, T_ifar, T_win, TRIALS, nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mgr->Add(co_gr[bin_count]);
        }
        // sel.Form("M_{total} #in [%3.2f, %3.2f]
        // M_{#odot}",MTotBins[bin_count],MTotBins[bin_count+1]);
        // leg->AddEntry(co_gr[bin_count], sel, "l");
        bin_count++;
    }

    mgr->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mgr->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
    mgr->GetXaxis()->SetLimits(MINRHO, MAXRHO);
    // mgr->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mgr->GetXaxis()->SetTitleOffset(1.3);
    mgr->GetYaxis()->SetTitleOffset(1.25);
    mgr->GetXaxis()->SetTickLength(0.01);
    mgr->GetYaxis()->SetTickLength(0.01);
    mgr->GetXaxis()->CenterTitle(kTRUE);
    mgr->GetYaxis()->CenterTitle(kTRUE);
    mgr->GetXaxis()->SetTitleFont(42);
    mgr->GetXaxis()->SetLabelFont(42);
    mgr->GetYaxis()->SetTitleFont(42);
    mgr->GetYaxis()->SetLabelFont(42);

    mgr->Draw("aple3");
    leg->SetFillColorAlpha(0, 1.0);
    leg->Draw();

    sprintf(fname, "%s/ROC_rho1_Mtot.png", netdir);
    co_canvas->SetLogy();
    co_canvas->Update();
    co_canvas->SaveAs(fname);
    co_canvas->SetLogy(0);

    //========== Sensitive Distance vs IFAR for varius Xeff bins
    //================
    cout << "Sensitive Distance vs IFAR ROC plot: Chieff" << endl;
    TMultiGraph* mg2 = new TMultiGraph();
    TLegend* leg2 = new TLegend(0.6, 0.8, 0.9, 0.9, "", "tlNDC");
    bin_count = 0;
    while (bin_count < 3) {
        sel.Form(
            "(spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])>= %f && "
            "(spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])< %f",
            ChieffBins[bin_count], ChieffBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume,
            col_arr2[bin_count], opt, liveTot, T_ifar, T_win, TRIALS, nIFO, VT,
            Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg2->Add(co_gr[bin_count]);
            sel.Form("#chi_{eff} #in [%3.1f, %3.1f]", ChieffBins[bin_count],
                     ChieffBins[bin_count + 1]);
            leg2->AddEntry(co_gr[bin_count], sel, "l");
        }
        bin_count++;
    }

    mg2->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg2->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
    mg2->GetXaxis()->SetLimits(MINIFAR, MAXIFAR);
    // mg2->GetYaxis()->SetRangeUser(200,MAXRADIUS);
    mg2->GetXaxis()->SetTitleOffset(1.3);
    mg2->GetYaxis()->SetTitleOffset(1.25);
    mg2->GetXaxis()->SetTickLength(0.01);
    mg2->GetYaxis()->SetTickLength(0.01);
    mg2->GetXaxis()->CenterTitle(kTRUE);
    mg2->GetYaxis()->CenterTitle(kTRUE);
    mg2->GetXaxis()->SetTitleFont(42);
    mg2->GetXaxis()->SetLabelFont(42);
    mg2->GetYaxis()->SetTitleFont(42);
    mg2->GetYaxis()->SetLabelFont(42);

    mg2->Draw("aple3");

    // char lab[256];
    leg2->SetFillColor(0);
    leg2->SetFillColorAlpha(0, 0.7);
    leg2->Draw();

    sprintf(fname, "%s/ROC_IFAR_chieff.png", netdir);
    co_canvas->Update();
    co_canvas->SaveAs(fname);

    //========== Sensitive Distance vs IFAR for varius Xp bins
    //================
    cout << "Sensitive Distance vs IFAR ROC plot: Chip" << endl;
    TMultiGraph* mg4 = new TMultiGraph();
    TLegend* leg4 = new TLegend(0.6, 0.8, 0.9, 0.9, "", "tlNDC");
    bin_count = 0;
    while (bin_count < 3) {
        sel.Form("chip >= %f && chip < %f", ChipBins[bin_count],
                 ChipBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume,
            col_arr2[bin_count], opt, liveTot, T_ifar, T_win, TRIALS, nIFO, VT,
            Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg4->Add(co_gr[bin_count]);
            sel.Form("#chi_{p} #in [%3.1f, %3.1f]", ChipBins[bin_count],
                     ChipBins[bin_count + 1]);
            leg4->AddEntry(co_gr[bin_count], sel, "l");
        }
        bin_count++;
    }

    mg4->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg4->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
    mg4->GetXaxis()->SetLimits(MINIFAR, MAXIFAR);
    // mg2->GetYaxis()->SetRangeUser(200,MAXRADIUS);
    mg4->GetXaxis()->SetTitleOffset(1.3);
    mg4->GetYaxis()->SetTitleOffset(1.25);
    mg4->GetXaxis()->SetTickLength(0.01);
    mg4->GetYaxis()->SetTickLength(0.01);
    mg4->GetXaxis()->CenterTitle(kTRUE);
    mg4->GetYaxis()->CenterTitle(kTRUE);
    mg4->GetXaxis()->SetTitleFont(42);
    mg4->GetXaxis()->SetLabelFont(42);
    mg4->GetYaxis()->SetTitleFont(42);
    mg4->GetYaxis()->SetLabelFont(42);

    mg4->Draw("aple3");

    // char lab[256];
    leg4->SetFillColor(0);
    leg4->SetFillColorAlpha(0, 0.7);
    leg4->Draw();

    sprintf(fname, "%s/ROC_IFAR_chip.png", netdir);
    co_canvas->Update();
    co_canvas->SaveAs(fname);

    //========== Sensitive Distance vs IFAR for varius Cos(iota) bins
    //================
    cout << "Sensitive Distance vs IFAR ROC plot: Inclination" << endl;
    TMultiGraph* mg6 = new TMultiGraph();
    TLegend* leg6 = new TLegend(0.6, 0.8, 0.9, 0.9, "", "tlNDC");
    bin_count = 0;
    while (bin_count < 3) {
        sel.Form("iota[1] >= %f && iota[1] < %f", CosiBins[bin_count],
                 CosiBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume,
            col_arr2[bin_count], opt, liveTot, T_ifar, T_win, TRIALS, nIFO, VT,
            Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg6->Add(co_gr[bin_count]);
            sel.Form("#Cos(#iota) #in [%s, %s]", CosiBinsT[bin_count].Data(),
                     CosiBinsT[bin_count + 1].Data());
            leg6->AddEntry(co_gr[bin_count], sel, "l");
        }
        bin_count++;
    }

    mg6->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg6->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
    mg6->GetXaxis()->SetLimits(MINIFAR, MAXIFAR);
    // mg2->GetYaxis()->SetRangeUser(200,MAXRADIUS);
    mg6->GetXaxis()->SetTitleOffset(1.3);
    mg6->GetYaxis()->SetTitleOffset(1.25);
    mg6->GetXaxis()->SetTickLength(0.01);
    mg6->GetYaxis()->SetTickLength(0.01);
    mg6->GetXaxis()->CenterTitle(kTRUE);
    mg6->GetYaxis()->CenterTitle(kTRUE);
    mg6->GetXaxis()->SetTitleFont(42);
    mg6->GetXaxis()->SetLabelFont(42);
    mg6->GetYaxis()->SetTitleFont(42);
    mg6->GetYaxis()->SetLabelFont(42);

    mg6->Draw("aple3");

    // char lab[256];
    leg6->SetFillColor(0);
    leg6->SetFillColorAlpha(0, 0.7);
    leg6->Draw();

    sprintf(fname, "%s/ROC_IFAR_iota.png", netdir);
    co_canvas->Update();
    co_canvas->SaveAs(fname);

    //========== Sensitive Distance vs rho[1] for varius Xeff bins
    //================
    cout << "Sensitive Distance vs rho[1] ROC plot: Chieff" << endl;
    TMultiGraph* mg2r = new TMultiGraph();
    bin_count = 0;
    while (bin_count < 3) {
        sel.Form(
            "(spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])>= %f && "
            "(spin[2]*mass[0]+spin[5]*mass[1])/(mass[1]+mass[0])< %f",
            ChieffBins[bin_count], ChieffBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume,
            col_arr2[bin_count], opt + " rho", liveTot, T_ifar, T_win, TRIALS,
            nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg2r->Add(co_gr[bin_count]);
            // sel.Form("#chi_{eff} #in [%3.1f, %3.1f]", ChieffBins[bin_count],
            // ChieffBins[bin_count+1]); leg2->AddEntry(co_gr[bin_count], sel,
            // "l");
        }
        bin_count++;
    }
    mg2r->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg2r->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
    mg2r->GetXaxis()->SetLimits(MINRHO, MAXRHO);
    // mg2r->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mg2r->GetXaxis()->SetTitleOffset(1.3);
    mg2r->GetYaxis()->SetTitleOffset(1.25);
    mg2r->GetXaxis()->SetTickLength(0.01);
    mg2r->GetYaxis()->SetTickLength(0.01);
    mg2r->GetXaxis()->CenterTitle(kTRUE);
    mg2r->GetYaxis()->CenterTitle(kTRUE);
    mg2r->GetXaxis()->SetTitleFont(42);
    mg2r->GetXaxis()->SetLabelFont(42);
    mg2r->GetYaxis()->SetTitleFont(42);
    mg2r->GetYaxis()->SetLabelFont(42);

    mg2r->Draw("aple3");
    leg2->SetFillColorAlpha(0, 1.0);
    leg2->Draw();

    sprintf(fname, "%s/ROC_rho1_chieff.png", netdir);
    co_canvas->SetLogy();
    co_canvas->Update();
    co_canvas->SaveAs(fname);
    co_canvas->SetLogy(0);

    //========== Sensitive Distance vs rho[1] for varius Xp bins
    //================
    cout << "Sensitive Distance vs rho[1] ROC plot: Chip" << endl;
    TMultiGraph* mg4r = new TMultiGraph();
    bin_count = 0;
    while (bin_count < 3) {
        sel.Form("chip >= %f && chip < %f", ChipBins[bin_count],
                 ChipBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume,
            col_arr2[bin_count], opt + " rho", liveTot, T_ifar, T_win, TRIALS,
            nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg4r->Add(co_gr[bin_count]);
            // sel.Form("#chi_{p} #in [%3.1f, %3.1f]", ChipBins[bin_count],
            // ChipBins[bin_count+1]); leg4->AddEntry(co_gr[bin_count], sel,
            // "l");
        }
        bin_count++;
    }

    mg4r->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg4r->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
    mg4r->GetXaxis()->SetLimits(MINRHO, MAXRHO);
    // mg2r->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mg4r->GetXaxis()->SetTitleOffset(1.3);
    mg4r->GetYaxis()->SetTitleOffset(1.25);
    mg4r->GetXaxis()->SetTickLength(0.01);
    mg4r->GetYaxis()->SetTickLength(0.01);
    mg4r->GetXaxis()->CenterTitle(kTRUE);
    mg4r->GetYaxis()->CenterTitle(kTRUE);
    mg4r->GetXaxis()->SetTitleFont(42);
    mg4r->GetXaxis()->SetLabelFont(42);
    mg4r->GetYaxis()->SetTitleFont(42);
    mg4r->GetYaxis()->SetLabelFont(42);

    mg4r->Draw("aple3");
    leg4->SetFillColorAlpha(0, 1.0);
    leg4->Draw();

    sprintf(fname, "%s/ROC_rho1_chip.png", netdir);
    co_canvas->SetLogy();
    co_canvas->Update();
    co_canvas->SaveAs(fname);
    co_canvas->SetLogy(0);

    //========== Sensitive Distance vs rho[1] for varius Cos iota  bins
    //================
    cout << "Sensitive Distance vs rho[1] ROC plot: Inclination" << endl;
    TMultiGraph* mg6r = new TMultiGraph();
    bin_count = 0;
    while (bin_count < 3) {
        sel.Form("iota[1] >= %f && iota[1] < %f", CosiBins[bin_count],
                 CosiBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume,
            col_arr2[bin_count], opt + " rho", liveTot, T_ifar, T_win, TRIALS,
            nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg6r->Add(co_gr[bin_count]);
            // sel.Form("#Cos(#iota) #in [%s, %s]", CosiBinsT[bin_count],
            // CosiBinsT[bin_count+1]); leg6->AddEntry(co_gr[bin_count], sel,
            // "l");
        }
        bin_count++;
    }

    mg6r->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg6r->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
    mg6r->GetXaxis()->SetLimits(MINRHO, MAXRHO);
    // mg2r->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mg6r->GetXaxis()->SetTitleOffset(1.3);
    mg6r->GetYaxis()->SetTitleOffset(1.25);
    mg6r->GetXaxis()->SetTickLength(0.01);
    mg6r->GetYaxis()->SetTickLength(0.01);
    mg6r->GetXaxis()->CenterTitle(kTRUE);
    mg6r->GetYaxis()->CenterTitle(kTRUE);
    mg6r->GetXaxis()->SetTitleFont(42);
    mg6r->GetXaxis()->SetLabelFont(42);
    mg6r->GetYaxis()->SetTitleFont(42);
    mg6r->GetYaxis()->SetLabelFont(42);

    mg6r->Draw("aple3");
    leg6->SetFillColorAlpha(0, 1.0);
    leg6->Draw();

    sprintf(fname, "%s/ROC_rho1_iota.png", netdir);
    co_canvas->SetLogy();
    co_canvas->Update();
    co_canvas->SaveAs(fname);
    co_canvas->SetLogy(0);

    //========== Sensitive Distance vs IFAR for varius Mchirp bins
    //================
    cout << "Sensitive Distance vs IFAR ROC plot: Mchirp" << endl;
    TMultiGraph* mg3 = new TMultiGraph();
    TLegend* leg3 = new TLegend(0.6, 0.7, 0.9, 0.9, "", "brNDC");
    bin_count = 0;
    while (bin_count < 4) {
        sel.Form("(chirp[0] >= %f) && (chirp[0] < %f)", MChirpBins[bin_count],
                 MChirpBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume, col_arr[bin_count],
            opt, liveTot, T_ifar, T_win, TRIALS, nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg3->Add(co_gr[bin_count]);
            sel.Form("M_{chirp} #in [%3.2f, %3.2f] M_{#odot}",
                     MTotBins[bin_count], MTotBins[bin_count + 1]);
            leg3->AddEntry(co_gr[bin_count], sel, "l");
        }
        bin_count++;
    }

    mg3->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg3->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
    mg3->GetXaxis()->SetLimits(MINIFAR, MAXIFAR);
    // mg3->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mg3->GetXaxis()->SetTitleOffset(1.3);
    mg3->GetYaxis()->SetTitleOffset(1.25);
    mg3->GetXaxis()->SetTickLength(0.01);
    mg3->GetYaxis()->SetTickLength(0.01);
    mg3->GetXaxis()->CenterTitle(kTRUE);
    mg3->GetYaxis()->CenterTitle(kTRUE);
    mg3->GetXaxis()->SetTitleFont(42);
    mg3->GetXaxis()->SetLabelFont(42);
    mg3->GetYaxis()->SetTitleFont(42);
    mg3->GetYaxis()->SetLabelFont(42);

    mg3->Draw("aple3");

    // char lab[256];
    leg3->SetFillColor(0);
    leg3->SetFillColorAlpha(0, 0.7);
    leg3->Draw();

    sprintf(fname, "%s/ROC_IFAR_chirp.png", netdir);
    co_canvas->Update();
    co_canvas->SaveAs(fname);

    //========== Sensitive Distance vs rho[1] for varius Mchirp bins
    //================
    cout << "Sensitive Distance vs rho[1] ROC plot: Mchirp" << endl;
    TMultiGraph* mg3r = new TMultiGraph();
    bin_count = 0;
    while (bin_count < 4) {
        sel.Form("(chirp[0] >= %f) && (chirp[0] < %f)", MChirpBins[bin_count],
                 MChirpBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume, col_arr[bin_count],
            opt + " rho", liveTot, T_ifar, T_win, TRIALS, nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg3r->Add(co_gr[bin_count]);
            // sel.Form("M_{chirp} #in [%3.2f, %3.2f]
            // M_{#odot}",MTotBins[bin_count],MTotBins[bin_count+1]);
            // leg->AddEntry(co_gr[bin_count], sel, "l");
        }
        bin_count++;
    }

    mg3r->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg3r->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
    mg3r->GetXaxis()->SetLimits(MINRHO, MAXRHO);
    // mg3r->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mg3r->GetXaxis()->SetTitleOffset(1.3);
    mg3r->GetYaxis()->SetTitleOffset(1.25);
    mg3r->GetXaxis()->SetTickLength(0.01);
    mg3r->GetYaxis()->SetTickLength(0.01);
    mg3r->GetXaxis()->CenterTitle(kTRUE);
    mg3r->GetYaxis()->CenterTitle(kTRUE);
    mg3r->GetXaxis()->SetTitleFont(42);
    mg3r->GetXaxis()->SetLabelFont(42);
    mg3r->GetYaxis()->SetTitleFont(42);
    mg3r->GetYaxis()->SetLabelFont(42);

    mg3r->Draw("aple3");
    leg3->SetFillColorAlpha(0, 1.0);
    leg3->Draw();

    sprintf(fname, "%s/ROC_rho1_chirp.png", netdir);
    co_canvas->SetLogy();
    co_canvas->Update();
    co_canvas->SaveAs(fname);

    //========== Sensitive Distance vs IFAR for varius Mass Ratio bins
    //================
    cout << "Sensitive Distance vs IFAR ROC plot: Symmetric Mass Ratio" << endl;
    TMultiGraph* mg5 = new TMultiGraph();
    TLegend* leg5 = new TLegend(0.6, 0.7, 0.9, 0.9, "", "brNDC");
    bin_count = 0;
    while (bin_count < 4) {
        sel.Form(
            "(mass[0]*mass[1]/pow(mass[0]+mass[1],2.0) <= %f) && "
            "(mass[0]*mass[1]/pow(mass[0]+mass[1],2.0) > %f)",
            EtaBins[bin_count], EtaBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume, col_arr[bin_count],
            opt, liveTot, T_ifar, T_win, TRIALS, nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg5->Add(co_gr[bin_count]);
            sel.Form("#eta #in [%3.3f, %3.3f]", EtaBins[bin_count + 1],
                     EtaBins[bin_count]);
            leg5->AddEntry(co_gr[bin_count], sel, "l");
        }
        bin_count++;
    }

    mg5->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg5->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
    mg5->GetXaxis()->SetLimits(MINIFAR, MAXIFAR);
    // mg->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mg5->GetXaxis()->SetTitleOffset(1.3);
    mg5->GetYaxis()->SetTitleOffset(1.25);
    mg5->GetXaxis()->SetTickLength(0.01);
    mg5->GetYaxis()->SetTickLength(0.01);
    mg5->GetXaxis()->CenterTitle(kTRUE);
    mg5->GetYaxis()->CenterTitle(kTRUE);
    mg5->GetXaxis()->SetTitleFont(42);
    mg5->GetXaxis()->SetLabelFont(42);
    mg5->GetYaxis()->SetTitleFont(42);
    mg5->GetYaxis()->SetLabelFont(42);

    mg5->Draw("aple3");

    // char lab[256];
    // leg5->SetFillColor(0);
    leg5->SetFillColorAlpha(0, 0.9);
    leg5->Draw();

    sprintf(fname, "%s/ROC_IFAR_eta.png", netdir);
    co_canvas->SetLogy(0);
    co_canvas->Update();
    co_canvas->SaveAs(fname);

    //========== Sensitive Distance vs rho[1] for varius mass ratio bins
    //================
    cout << "Sensitive Distance vs rho[1] ROC plot: Symmetric Mass Ratio"
         << endl;
    TMultiGraph* mg5r = new TMultiGraph();
    bin_count = 0;
    while (bin_count < 4) {
        sel.Form(
            "(mass[0]*mass[1]/pow(mass[0]+mass[1],2.0) <= %f) && "
            "(mass[0]*mass[1]/pow(mass[0]+mass[1],2.0) > %f)",
            EtaBins[bin_count], EtaBins[bin_count + 1]);
        co_gr[bin_count] = cbcTool.CreateGraphRadiusIFAR(
            sim_file_name, mdc_file_name, sel, shell_volume, col_arr[bin_count],
            opt + " rho", liveTot, T_ifar, T_win, TRIALS, nIFO, VT, Tscale);
        if (co_gr[bin_count]->GetN() > 1) {
            mg5r->Add(co_gr[bin_count]);
            // sel.Form("#eta #in [%3.2f,
            // %3.2f]",EtaBins[bin_count],EtaBins[bin_count+1]);
            // leg5->AddEntry(co_gr[bin_count], sel, "l");
        }
        bin_count++;
    }

    mg5r->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
    mg5r->GetXaxis()->SetTitle("Magnitude Test Statistic (rho[1])");
    mg5r->GetXaxis()->SetLimits(MINRHO, MAXRHO);
    // mg3r->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
    mg5r->GetXaxis()->SetTitleOffset(1.3);
    mg5r->GetYaxis()->SetTitleOffset(1.25);
    mg5r->GetXaxis()->SetTickLength(0.01);
    mg5r->GetYaxis()->SetTickLength(0.01);
    mg5r->GetXaxis()->CenterTitle(kTRUE);
    mg5r->GetYaxis()->CenterTitle(kTRUE);
    mg5r->GetXaxis()->SetTitleFont(42);
    mg5r->GetXaxis()->SetLabelFont(42);
    mg5r->GetYaxis()->SetTitleFont(42);
    mg5r->GetYaxis()->SetLabelFont(42);

    mg5r->Draw("aple3");
    leg5->SetFillColorAlpha(0, 1.0);
    leg5->Draw();

    sprintf(fname, "%s/ROC_rho1_eta.png", netdir);
    co_canvas->SetLogy(1);
    co_canvas->Update();
    co_canvas->SaveAs(fname);

    delete co_canvas;
    // delete mg,mgr,mg2,mg2r,mg3,mg3r;
    delete mg, mgr, mg2, mg2r, mg3, mg3r, mg4, mg4r, mg5, mg5r, mg6, mg6r;
    delete leg, leg2, leg3, leg4, leg5, leg6;
    // delete co_gr1, co_gr2, co_gr3, co_gr4;
    // delete[] co_gr;
    // exit(0);
}
