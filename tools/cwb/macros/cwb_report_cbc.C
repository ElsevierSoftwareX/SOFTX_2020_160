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



#define RHO_MIN 5.0
#define RHO_BIN 0.1
#define RHO_NBINS 5000
#define CONTOURS 7

{
    //###############################################################################################################################
    // Palette
    //###############################################################################################################################
#define NRGBs 6
#define NCont 99
gStyle->SetPalette(57);
/*
    gStyle->SetNumberContours(NCont);
    double stops[NRGBs] = {0.10, 0.25, 0.45, 0.60, 0.75, 1.00};
    double red[NRGBs] = {0.00, 0.00, 0.00, 0.97, 0.97, 0.10};
    double green[NRGBs] = {0.97, 0.30, 0.40, 0.97, 0.00, 0.00};
    double blue[NRGBs] = {0.97, 0.97, 0.00, 0.00, 0.00, 0.00};
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
*/
    // -----------------------------------------------------------
    // CBC
    // -----------------------------------------------------------

    int NBINS_mass = (int)((MAX_MASS - MIN_MASS) / MASS_BIN);

#ifndef MIN_plot_mass1
    float MIN_plot_mass1 = 0.0;
#endif
#ifndef MIN_plot_mass2
    float MIN_plot_mass2 = 0.0;
#endif
#ifndef MAX_plot_mass1
    float MAX_plot_mass1 = 150.0;
#endif
#ifndef MAX_plot_mass2
    float MAX_plot_mass2 = 100.0;
#endif
#ifndef MAX_EFFECTIVE_RADIUS
    float MAX_EFFECTIVE_RADIUS = 3000.;
#endif
#ifndef MASS_BIN
    float MASS_BIN = 15.0;
#endif
#ifndef MIN_MASS
    float MIN_MASS = 0.0;
#endif
#ifndef MAX_MASS
    float MAX_MASS = 150.0;
#endif

#ifdef MAX_MASS2
    float max_mass2 = MAX_MASS2;
#else
    float max_mass2 = MAX_MASS;
#endif

#ifdef MIN_MASS1
    float min_mass1 = MIN_MASS1;
#else
    float min_mass1 = MIN_MASS;
#endif

#ifdef MAX_MASS1
    float max_mass1 = MAX_MASS1;
#else
    float max_mass1 = MAX_MASS;
#endif

#ifdef MIN_MASS2
    float min_mass2 = MIN_MASS2;
#else
    float min_mass2 = MIN_MASS;
#endif

    int NBINS_mass1 = (int)((max_mass1 - min_mass1) / MASS_BIN);
    int NBINS_mass2 = (int)((max_mass2 - min_mass2) / MASS_BIN);

#ifndef CHI_BIN
    float CHI_BIN = 0.66666;
#endif
#ifndef MINCHI
    float MINCHI = -1.0;
#endif
#ifndef MAXCHI
    float MAXCHI = 1.0;
#endif
    cout << "cwb_report_cbc starts..." << endl;
    
    //Search types: BBH or IMBHB
    TString SearchType = gSystem->Getenv("CBC_SEARCH_TYPE");
    cout << "CBC Search type: " <<  gSystem->Getenv("CBC_SEARCH_TYPE") <<endl;
    //  cout << "Mass1 : ["<<MIN_MASS<<","<<MAX_MASS<<"] with
    //  "<<NBINS_mass<<" bins"<<endl;
    cout << "Mass1 : [" << min_mass1 << "," << max_mass1 << "] with "
         << NBINS_mass1 << " bins" << endl;
    cout << "Mass2 : [" << min_mass2 << "," << max_mass2 << "] with "
         << NBINS_mass2 << " bins" << endl;
         
    CWB::Toolbox TB;
    CWB::CBCTool cbcTool;

    TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
    TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
    TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
    TB.checkFile(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
    TB.checkFile(gSystem->Getenv("CWB_UPPARAMETERS_FILE"));
    TB.checkFile(gSystem->Getenv("CWB_EPPARAMETERS_FILE"));

    TB.mkDir(netdir, true);

#ifndef _USE_ROOT6
// the CWB_CAT_NAME declared in CWB_EPPARAMETERS_FILE is not visible. why?
// the include is defined again
#undef GTOOLBOX_HH
#endif
#include "GToolbox.hh"

    gStyle->SetTitleFillColor(kWhite);
    // FgStyle->SetLineColor(kWhite);
    gStyle->SetNumberContours(256);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetStatBorderSize(1);
    gStyle->SetOptStat(kFALSE);

    // remove the red box around canvas
    gStyle->SetFrameBorderMode(0);
    gROOT->ForceStyle();

    char networkname[256];
    if (strlen(ifo[0]) > 0)
        strcpy(networkname, ifo[0]);
    else
        strcpy(networkname, detParms[0].name);
    for (int i = 1; i < nIFO; i++) {
        if (strlen(ifo[i]) > 0)
            sprintf(networkname, "%s-%s", networkname,
                    ifo[i]);  // built in detector
        else
            sprintf(networkname, "%s-%s", networkname,
                    detParms[i].name);  // user define detector
    }
    // int gIFACTOR=1;
    //  double FACTORS[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // float CYS = 31556952;  //~86400.*365.25;  //Exact conversion from
    // year to seconds
    float CYS =
        86400. * 365.25;  // Exact conversion from Julian year to seconds
    network* net = NULL;
    CWB::config* cfg = new CWB::config;
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
    float* minMtot = new float[nfactor];
    float* maxMtot = new float[nfactor];
    float* minMChirp = new float[nfactor];
    float* maxMChirp = new float[nfactor];
    float* minDistanceXML = new float[nfactor];
    float* maxDistanceXML = new float[nfactor];
    float* minDistance = new float[nfactor];
    float* maxDistance = new float[nfactor];
    float* minRatio = new float[nfactor];
    float* maxRatio = new float[nfactor];
    float* shell_volume = new float[nfactor];
    bool DDistrUniform, DDistrVolume, FixedFiducialVolume, DDistrChirpMass,
        bminMtot, bmaxMtot, bminDistance, bmaxDistance, bminRatio, bmaxRatio,
        Redshift, PointMasses, minchi, bminMass1, bminMass2, bmaxMass1,
        bmaxMass2;
    TString* waveforms = new TString[nfactor];
    float* FACTORS = new float[nfactor];
    float ShellminDistance = 9999999.0;
    float ShellmaxDistance = 0.0;
    TString opt = SearchType;

    // break;
    cout << "Number of Factors:" << nfactor << endl;
    for (int l = 0; l < nfactor; l++) {
        gIFACTOR = l + 1;
        FACTORS[l] = gIFACTOR;
        gROOT->Macro(configPlugin.GetTitle());
        TString Insp = MDC.GetInspiral();
        if (MDC.GetInspiralOption("--waveform") != "") {
            waveforms[gIFACTOR - 1] = MDC.GetInspiralOption("--waveform");
        }
        if (MDC.GetInspiralOption("--min-mtotal") != "") {
            minMtot[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--min-mtotal").Atof();
            bminMtot = 1;
        }
        if (MDC.GetInspiralOption("--max-mtotal") != "") {
            maxMtot[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--max-mtotal").Atof();
            bmaxMtot = 1;
        }

        if (MDC.GetInspiralOption("--min-distance") != "") {
            minDistanceXML[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--min-distance").Atof();
            bminDistance = 1;
        }
        if (ShellminDistance > minDistanceXML[gIFACTOR - 1]) {
            ShellminDistance = minDistanceXML[gIFACTOR - 1];
        }
        if (MDC.GetInspiralOption("--max-distance") != "") {
            maxDistanceXML[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--max-distance").Atof();
            bmaxDistance = 1;
        }
        if (ShellmaxDistance < maxDistanceXML[gIFACTOR - 1]) {
            ShellmaxDistance = maxDistanceXML[gIFACTOR - 1];
        }

        if (MDC.GetInspiralOption("--min-mratio") != "") {
            minRatio[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--min-mratio").Atof();
            bminRatio = 1;
        }
        if ((MDC.GetInspiralOption("--min-mass1") != "") &&
            (MDC.GetInspiralOption("--min-mass2") != "")) {
            minRatio[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--min-mass1").Atof() /
                MDC.GetInspiralOption("--min-mass2").Atof();
            bminRatio = 1;
        }
        if (MDC.GetInspiralOption("--max-mratio") != "") {
            maxRatio[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--max-mratio").Atof();
            bmaxRatio = 1;
        }
        if ((MDC.GetInspiralOption("--min-mass1") != "") &&
            (MDC.GetInspiralOption("--max-mass2") != "")) {
            minRatio[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--min-mass1").Atof() /
                (float)MDC.GetInspiralOption("--max-mass2").Atof();
            bminRatio = 1;
            bminMass1 = 1;
            bmaxMass2 = 1;
        }
        if ((MDC.GetInspiralOption("--max-mass1") != "") &&
            (MDC.GetInspiralOption("--min-mass2") != "")) {
            maxRatio[gIFACTOR - 1] =
                (float)MDC.GetInspiralOption("--max-mass1").Atof() /
                (float)MDC.GetInspiralOption("--min-mass2").Atof();
            bmaxRatio = 1;
            bminMass2 = 1;
            bmaxMass1 = 1;
        }
        if (MDC.GetInspiralOption("--d-distr") != "") {
            if (MDC.GetInspiralOption("--d-distr").CompareTo("uniform") == 0) {
                DDistrUniform = 1;
                opt = "DDistrUniform";
                cout << "Uniform distribution in linear distance" << endl;
                shell_volume[gIFACTOR - 1] =
                    4. * TMath::Pi() *
                    (maxDistanceXML[gIFACTOR - 1] / 1000. -
                     minDistanceXML[gIFACTOR - 1] / 1000.);
            } else if (MDC.GetInspiralOption("--d-distr").CompareTo("volume") ==
                       0) {
                DDistrVolume = 1;
                opt = "DDistrVolume";
                cout << "Uniform distribution in volume" << endl;
                shell_volume[gIFACTOR - 1] =
                    4. * TMath::Pi() *
                    (pow(maxDistanceXML[gIFACTOR - 1] / 1000., 3) -
                     pow(minDistanceXML[gIFACTOR - 1] / 1000., 3)) /
                    3;
            } else {
                cout << "No defined distance distribution? "
                        "Or different from uniform and volume?"
                     << endl;
                exit(1);
            }

            cout << "Shell volume: " << shell_volume[gIFACTOR - 1] << endl;
        }
        if (MDC.GetInspiralOption("--dchirp-distr") != "") {
            if (MDC.GetInspiralOption("--dchirp-distr").CompareTo("uniform") ==
                0) {
                DDistrChirpMass = 1;
                opt = "DDistrChirpMass";
                shell_volume[gIFACTOR - 1] =
                    4. * TMath::Pi() *
                    (maxDistanceXML[gIFACTOR - 1] / 1000. -
                     minDistanceXML[gIFACTOR - 1] / 1000.);
                maxMChirp[gIFACTOR - 1] =
                    maxMtot[gIFACTOR - 1] *
                    pow(minRatio[gIFACTOR - 1], 3. / 5.) /
                    pow(1 + minRatio[gIFACTOR - 1], 6. / 5.);
                minMChirp[gIFACTOR - 1] =
                    minMtot[gIFACTOR - 1] *
                    pow(maxRatio[gIFACTOR - 1], 3. / 5.) /
                    pow(1 + maxRatio[gIFACTOR - 1], 6. / 5.);
                maxDistance[gIFACTOR - 1] =
                    maxDistanceXML[gIFACTOR - 1] *
                    pow(maxMChirp[gIFACTOR - 1] / 1.22,
                        5. / 6.);  // Rescale the distances to
                                   // get the extremes for chirp
                                   // distances
                if (ShellmaxDistance < maxDistance[gIFACTOR - 1]) {
                    ShellmaxDistance = maxDistance[gIFACTOR - 1];
                }
                minDistance[gIFACTOR - 1] =
                    minDistanceXML[gIFACTOR - 1] *
                    pow(minMChirp[gIFACTOR - 1] / 1.22, 5. / 6.);
                cout << "Uniform distribution in Chirp Mass "
                        "distance"
                     << endl;
                cout << "MaxDistance: " << maxDistance[gIFACTOR - 1] << endl;
            } else {
                cout << "No defined distance "
                        "distribution? Or different from "
                        "uniform in distance?"
                     << endl;
                exit(1);
            }
        }
        cout << "MDC set: " << gIFACTOR << endl;
        cout << "xml conf: waveform=" << waveforms[gIFACTOR - 1]
             << " minMtot=" << minMtot[gIFACTOR - 1]
             << " maxMtot=" << maxMtot[gIFACTOR - 1]
             << " minDistance=" << minDistance[gIFACTOR - 1]
             << " maxDistance=" << maxDistance[gIFACTOR - 1]
             << " minRatio=" << minRatio[gIFACTOR - 1]
             << " maxRatio=" << maxRatio[gIFACTOR - 1] << endl;
    }

    // break;

#ifdef FIXMAXDISTANCE
    bmaxDistance = 1;
    maxDistance[gIFACTOR - 1] = FIXMAXDISTANCE;
    shell_volume[gIFACTOR - 1] =
        4. / 3. * TMath::Pi() * pow(maxDistance[gIFACTOR - 1] / 1000., 3.);
#endif

#ifdef FIXMINRATIO
    bminRatio = 1;
    minRatio[gIFACTOR - 1] = FIXMINRATIO;

#endif
/*	if(bminRatio){
            float MINRATIO = minRatio[gIFACTOR - 1];
        } else {
            float MINRATIO = 0.02;
        }
            float MINRATIO = 0.02;
        cout << "MINRATIO = " << MINRATIO << endl;*/
#ifdef FIXMAXRATIO
    bmaxRatio = 1;
    maxRatio[gIFACTOR - 1] = FIXMAXRATIO;

#endif
/*	if(bmaxRatio) {
            float MAXRATIO = maxRatio[gIFACTOR - 1];
        } else {
            float MAXRATIO = 50.0;
        }
        cout << "MAXRATIO = " << MAXRATIO << endl;*/
#ifdef FIXMINTOT
    bminMtot = 1;
    minMtot[gIFACTOR - 1] = FIXMINTOT;
#endif

#ifdef FIXMAXTOT
    bmaxMtot = 1;
    maxMtot[gIFACTOR - 1] = FIXMAXTOT;
#endif

#ifdef REDSHIFT
    Redshift = 1;
    opt += " Redshift";
    // check if VT and TMC have been properly initialized in the user pp file
    if (VT == -1.0) {
        cout << "Undefined VT! When using the \"REDSHIFT\" option, VT has "
                "to be initialized to the proper fiducial volume times time"
                " for the simulation which has been injected."
             << endl;
        cout << " Add something like: VT = XY; in your user pp file,"
                " where XY is a double in [Gpc^3*years]."
             << endl;
        exit(1);
    }
    if (TMC == -1.0) {
        cout << "Undefined TMC! When using the \"REDSHIFT\" option, TMC has "
                "to be initialized to the proper MonteCarlo time "
                "for the simulation which has been injected."
             << endl;
        cout << "Add something like: TMC = XY; in your user pp file, "
                "where XY is an int in [s]."
             << endl;
        exit(1);
    }
#endif
#ifdef POINTMASSES
    PointMasses = 1;
#endif
    //#ifdef MINCHI
    // Setting default to filling all spin related histograms... TODO:
    // review and eventually remove all ""if" statements on minchi
    minchi = 1;
    //#endif

    float MINMtot = 0.0;
    if (bminMtot) {
        float MINMtot = 0.99 * minMtot[gIFACTOR - 1];
    }

    float MAXMtot = 100.0;
    int NBINS_MTOT = 0;
    // if (bmaxMtot) {
    // float MAXMtot = 1.01*maxMtot[gIFACTOR-1];
    // int NBINS_MTOT =
    // TMath::FloorNint((maxMtot[gIFACTOR-1]-minMtot[gIFACTOR-1])/MASS_BIN/2.);
    MAXMtot = MAX_MASS;
    NBINS_MTOT = TMath::FloorNint((MAX_MASS - MIN_MASS) / MASS_BIN / 2.);
    cout << "NBINS_MTOT: " << NBINS_MTOT << endl;

    //} else {
    if (!bmaxMtot) {
        cout << "Undefined maximal total mass!! Using default, i.e. "
                "100.0"
             << endl;
        // exit(1);
    }

    float MINDISTANCE = 0.0;
    if (bminDistance) {
        MINDISTANCE = 0.9 * ShellminDistance;
        cout << "MINDISTANCE = " << MINDISTANCE << endl;
    } else {
        cout << "Undefined minimal distance!! Using default, i.e. 0" << endl;
    }

    float MAXDISTANCE = 5000000;  // 5 Gpc
    if (bmaxDistance) {
        MAXDISTANCE = maxDistance[gIFACTOR - 1];
        //                    TMath::Max(ShellmaxDistance,
        //                    maxDistance[gIFACTOR - 1]);
        cout << "MAXDISTANCE = " << MAXDISTANCE << endl;
    } else {
        cout << "Undefined maximal distance !! Using default, i.e. 5 "
                "Gpc."
                "You can define a MAXDISTANCE in pp par file, e.g. "
                "#define FIXMINDISTANCE 5000000"
             << endl;
    }

    float MINCHIRP = 100.0;
    float MAXCHIRP = 0.0;
    float MINRATIO = 1.0;
    float MAXRATIO = 10.0;
    if ((bminRatio) && (bmaxRatio)) {
        MAXRATIO = maxRatio[gIFACTOR - 1];
        MINRATIO = minRatio[gIFACTOR - 1];  /// Temporary fix
    } else {
        cout << "Undefined min/max Ratio.. Using default [1; 10] " << endl;
    }
    MINCHIRP = MINMtot * pow(MAXRATIO, 3. / 5.) / pow(1 + MAXRATIO, 6. / 5.);
    // MAXCHIRP = MAXMtot * pow(MINRATIO, 3. / 5.) / pow(1 + MINRATIO, 6.
    // / 5.);
    MAXCHIRP = MAXMtot / pow(2, 6. / 5.);
    /*for (int l = 0; l < nfactor; l++) {
            if (MINRATIO > minRatio[l]) {
                    MINRATIO = minRatio[l];
            }
            if (MAXRATIO < maxRatio[l]) {
                    MAXRATIO = maxRatio[l];



            if (MINCHIRP > minMtot[l] * pow(maxRatio[l], 3. / 5.) /
                               pow(1 + minRatio[l], 6. / 5.)) {
                    MINCHIRP = minMtot[l] *
                               pow(maxRatio[l], 3. / 5.) /
                               pow(1 + minRatio[l], 6. / 5.);
            };
            if (MAXCHIRP < maxMtot[l] * pow(minRatio[l], 3. / 5.) /
                               pow(1 + minRatio[l], 6. / 5.)) {
                    MAXCHIRP = maxMtot[l] *
                               pow(minRatio[l], 3. / 5.) /
                               pow(1 + minRatio[l], 6. / 5.);
            };
    }*/

    for (int l = 0; l < nfactor - 1; l++) {
        if ((minDistanceXML[l] == minDistanceXML[l + 1]) &&
            (maxDistanceXML[l] == maxDistanceXML[l + 1])) {
            FixedFiducialVolume = 1;
        } else {
            // FixedFiducialVolume=0;
            FixedFiducialVolume = 1;
            cout << "Beware: different fiducial volumes for "
                    "different factors!!"
                 << endl;
            //   exit(1);
        }
    }
    //  cout << "Plotting bounds: MINMtot=" << MINMtot << " MAXMtot=" <<
    //  MAXMtot
    //      << " MINRATIO=" << MINRATIO << " MAXRATIO=" << MAXRATIO
    //     << " MINDISTANCE=" << MINDISTANCE << " MAXDISTANCE=" <<
    //     MAXDISTANCE
    //    << " MINCHIRP=" << MINCHIRP << " MAXCHIRP=" << MAXCHIRP << endl;
    // If not vetoes are defined, then the char string veto_not_vetoed is
    // forced to be equal to ch2
    if (strlen(veto_not_vetoed) == 0) {
        sprintf(veto_not_vetoed, "%s", ch2);
    }

    // break;
    //###############################################################################################################################
    // Definitions of Canvas and histograms
    //###############################################################################################################################

    TCanvas* c1 = new TCanvas("c1", "c1", 3, 47, 1000, 802);
    // c1->Range(-1.216392, -477.6306, 508.8988, 2814.609);
    // c1->SetFillColor(0);
    // c1->SetBorderMode(0);
    // c1->SetBorderSize(2);
    c1->SetGridx();
    c1->SetGridy();
    // c1->SetRightMargin(0.1154618);
    // c1->SetTopMargin(0.07642487);
    // c1->SetBottomMargin(0.1450777);

    TCanvas* c2 = new TCanvas("c2", "c2", 3, 47, 1000, 802);
    c2->Range(-1.216392, -477.6306, 508.8988, 2814.609);
    // c1->SetFillColor(0);
    // c1->SetBorderMode(0);
    // c1->SetBorderSize(2);
    //c2->SetGridx();
    //c2->SetGridy();
    c2->SetRightMargin(0.154618);
    c2->SetTopMargin(0.07642487);
    c2->SetBottomMargin(0.1450777);

    TH2F* inj_events = new TH2F("inj_events", "D_Minj", NBINS_mass, MIN_MASS,
                                MAX_MASS, NBINS_mass2, min_mass2, max_mass2);
    inj_events->GetXaxis()->SetRangeUser(MIN_plot_mass1, MAX_plot_mass1);
    inj_events->GetYaxis()->SetRangeUser(MIN_plot_mass2, MAX_plot_mass2);
    inj_events->GetXaxis()->SetTitle("Mass 1 (M_{#odot})");
    inj_events->GetYaxis()->SetTitle("Mass 2 (M_{#odot})");
    inj_events->GetXaxis()->SetTitleOffset(1.3);
    inj_events->GetYaxis()->SetTitleOffset(1.25);
    inj_events->GetXaxis()->CenterTitle(kTRUE);
    inj_events->GetYaxis()->CenterTitle(kTRUE);
    inj_events->GetXaxis()->SetNdivisions(410);
    inj_events->GetYaxis()->SetNdivisions(410);
    inj_events->GetXaxis()->SetTickLength(0.01);
    inj_events->GetYaxis()->SetTickLength(0.01);
    inj_events->GetZaxis()->SetTickLength(0.01);
    inj_events->GetXaxis()->SetTitleFont(42);
    inj_events->GetXaxis()->SetLabelFont(42);
    inj_events->GetYaxis()->SetTitleFont(42);
    inj_events->GetYaxis()->SetLabelFont(42);
    inj_events->GetZaxis()->SetLabelFont(42);
    inj_events->GetZaxis()->SetLabelSize(0.03);
    // inj_events->GetXaxis()->SetNdivisions(10,kFALSE);
    // inj_events->GetYaxis()->SetNdivisions(10,kFALSE);

    inj_events->SetTitle("");

    inj_events->SetContour(NCont);

    TH2F* rec_events = new TH2F("rec_events", "D_Mrec", NBINS_mass, MIN_MASS,
                                MAX_MASS, NBINS_mass2, min_mass2, max_mass2);
    rec_events->GetXaxis()->SetRangeUser(MIN_plot_mass1, MAX_plot_mass1);
    rec_events->GetYaxis()->SetRangeUser(MIN_plot_mass2, MAX_plot_mass2);
    rec_events->GetXaxis()->SetTitle("Mass 1 (M_{#odot})");
    rec_events->GetYaxis()->SetTitle("Mass 2 (M_{#odot})");
    rec_events->GetXaxis()->SetTitleOffset(1.3);
    rec_events->GetYaxis()->SetTitleOffset(1.25);
    rec_events->GetXaxis()->CenterTitle(kTRUE);
    rec_events->GetYaxis()->CenterTitle(kTRUE);
    rec_events->GetXaxis()->SetNdivisions(410);
    rec_events->GetYaxis()->SetNdivisions(410);
    rec_events->GetXaxis()->SetTickLength(0.01);
    rec_events->GetYaxis()->SetTickLength(0.01);
    rec_events->GetZaxis()->SetTickLength(0.01);
    rec_events->GetXaxis()->SetTitleFont(42);
    rec_events->GetXaxis()->SetLabelFont(42);
    rec_events->GetYaxis()->SetTitleFont(42);
    rec_events->GetYaxis()->SetLabelFont(42);
    rec_events->GetZaxis()->SetLabelFont(42);
    rec_events->GetZaxis()->SetLabelSize(0.03);
    rec_events->SetTitle("");

    rec_events->SetContour(NCont);
    TH2F* factor_events_inj = new TH2F[nfactor];
    // TH1* events_inj[nfactor];
    for (int i = 0; i < nfactor; i++) {
        factor_events_inj[i] =
            // new TH2F("factor_events_inj", "", NBINS_mass, MIN_MASS,
            TH2F("factor_events_inj", "", NBINS_mass, MIN_MASS, MAX_MASS,
                 NBINS_mass2, min_mass2, max_mass2);
        // events_inj[i] = new TH1("events_inj","",
        // int(maxDistance[i]/1000.), 0.0, maxDistance[i]/1000.);
    }
    TH2F* factor_events_rec =
        new TH2F("factor_events_rec", "", NBINS_mass, MIN_MASS, MAX_MASS,
                 NBINS_mass2, min_mass2, max_mass2);

    TH2F* D_Mtot_inj =
        new TH2F("Distance vs Mtot inj.", "", 1000, MINMtot, MAXMtot, 5000,
                 MINDISTANCE / 1000., MAXDISTANCE / 1000.);

    D_Mtot_inj->GetXaxis()->SetTitle("Total mass (M_{#odot})");
    D_Mtot_inj->GetYaxis()->SetTitle("Distance (Mpc)");
    D_Mtot_inj->GetXaxis()->SetTitleOffset(1.3);
    D_Mtot_inj->GetYaxis()->SetTitleOffset(1.3);
    D_Mtot_inj->GetXaxis()->SetTickLength(0.01);
    D_Mtot_inj->GetYaxis()->SetTickLength(0.01);
    D_Mtot_inj->GetXaxis()->CenterTitle(kTRUE);
    D_Mtot_inj->GetYaxis()->CenterTitle(kTRUE);
    D_Mtot_inj->GetXaxis()->SetTitleFont(42);
    D_Mtot_inj->GetXaxis()->SetLabelFont(42);
    D_Mtot_inj->GetYaxis()->SetTitleFont(42);
    D_Mtot_inj->GetYaxis()->SetLabelFont(42);
    D_Mtot_inj->SetMarkerStyle(20);
    D_Mtot_inj->SetMarkerSize(0.5);
    D_Mtot_inj->SetMarkerColor(2);
    D_Mtot_inj->SetTitle("");
    // D_Mtot_inj->Draw("ap");
    D_Mtot_inj->SetName("D_Mtotinj");

    int NBINS_SPIN = (int)((MAXCHI - MINCHI) / CHI_BIN);
    cout << "NBINS_SPIN: " << NBINS_SPIN << endl;
    TH2F* inj_events_spin_mtot =
        new TH2F("inj_spin_Mtot", "", NBINS_SPIN, MINCHI, MAXCHI, NBINS_MTOT,
                 MIN_MASS, MAX_MASS);
    TH2F* rec_events_spin_mtot =
        new TH2F("rec_spin_Mtot", "", NBINS_SPIN, MINCHI, MAXCHI, NBINS_MTOT,
                 MIN_MASS, MAX_MASS);
    TH2F* factor_events_spin_mtot_inj = new TH2F[nfactor];
    for (int i = 0; i < nfactor; i++) {
        factor_events_spin_mtot_inj[i] =
            // new TH2F("factor_events_spin_mtot_inj", "", NBINS_SPIN,
            TH2F("factor_events_spin_mtot_inj", "", NBINS_SPIN, MINCHI, MAXCHI,
                 NBINS_MTOT, MIN_MASS, MAX_MASS);
    }

    TH1F* Dt = new TH1F("Dt", "", 1000, -0.5, 0.5);
    Dt->SetMarkerStyle(20);
    Dt->SetMarkerSize(0.5);
    Dt->SetMarkerColor(4);

    TH2F* rhocc =
        new TH2F("rhocc", "", 100, 0., 1., 100, pp_rho_min, pp_rho_max);
    rhocc->SetTitle("0 < cc < 1");
    rhocc->SetTitleOffset(1.3, "Y");
    rhocc->GetXaxis()->SetTitle("network correlation");
    rhocc->GetYaxis()->SetTitle("#rho");
    rhocc->SetStats(kFALSE);
    rhocc->SetMarkerStyle(20);
    rhocc->SetMarkerSize(0.5);
    rhocc->SetMarkerColor(1);

    TH2F* rho_pf =
        new TH2F("rho_pf", "", 100, -1., 2., 100, pp_rho_min, pp_rho_max);
    rho_pf->SetTitle("chi2");
    rho_pf->GetXaxis()->SetTitle("log10(chi2)");
    rho_pf->SetTitleOffset(1.3, "Y");
    rho_pf->GetYaxis()->SetTitle("#rho");
    rho_pf->SetMarkerStyle(20);
    rho_pf->SetMarkerColor(1);
    rho_pf->SetMarkerSize(0.5);
    rho_pf->SetStats(kFALSE);

    Float_t xq[6] = {8.0, 15.0, 25.0, 35.0, 45.0, 55.0};
    Float_t* yq = new Float_t[101];
    for (int i = 0; i < 101; i++) {
        yq[i] = -50.0 + i;
    }
    TH2F* dchirp_rec =
        new TH2F("dchirp rec.", "Chirp Mass estimate", 5, xq, 100, yq);
    dchirp_rec->GetXaxis()->SetTitle("Chirp mass (M_{#odot})");
    dchirp_rec->GetYaxis()->SetTitle(
        "Chirp mass: injected-estimated (M_{#odot})");
    dchirp_rec->GetXaxis()->SetTitleOffset(1.3);
    dchirp_rec->GetYaxis()->SetTitleOffset(1.3);
    dchirp_rec->GetXaxis()->SetTickLength(0.01);
    dchirp_rec->GetYaxis()->SetTickLength(0.01);
    dchirp_rec->GetXaxis()->CenterTitle(kTRUE);
    dchirp_rec->GetYaxis()->CenterTitle(kTRUE);
    dchirp_rec->GetXaxis()->SetTitleFont(42);
    dchirp_rec->GetXaxis()->SetLabelFont(42);
    dchirp_rec->GetYaxis()->SetTitleFont(42);
    dchirp_rec->GetYaxis()->SetLabelFont(42);
    dchirp_rec->SetMarkerStyle(20);
    dchirp_rec->SetMarkerSize(0.5);
    dchirp_rec->SetMarkerColor(2);

    TH3F* D_dchirp_rec =
        new TH3F("Distance vs dchirp rec.", "", 5, 10.0, 50., 100, -50, 50.,
                 5000, MINDISTANCE / 1000., MAXDISTANCE / 1000);
    // TH3F *D_dchirp_rec = new TH3F("Distance vs dchirp
    // rec.","",5,xq,100,-50,50.,.5000,MINDISTANCE/1000.,MAXDISTANCE/1000);
    D_dchirp_rec->GetXaxis()->SetTitle("Chirp mass (M_{#odot})");
    D_dchirp_rec->GetYaxis()->SetTitle(
        "Chirp mass: injected-estimated (M_{#odot})");
    D_dchirp_rec->GetZaxis()->SetTitle("Distance (Mpc)");
    D_dchirp_rec->GetXaxis()->SetTitleOffset(1.3);
    D_dchirp_rec->GetYaxis()->SetTitleOffset(1.3);
    D_dchirp_rec->GetXaxis()->SetTickLength(0.01);
    D_dchirp_rec->GetYaxis()->SetTickLength(0.01);
    D_dchirp_rec->GetXaxis()->CenterTitle(kTRUE);
    D_dchirp_rec->GetYaxis()->CenterTitle(kTRUE);
    D_dchirp_rec->GetXaxis()->SetTitleFont(42);
    D_dchirp_rec->GetXaxis()->SetLabelFont(42);
    D_dchirp_rec->GetYaxis()->SetTitleFont(42);
    D_dchirp_rec->GetYaxis()->SetLabelFont(42);
    D_dchirp_rec->SetMarkerStyle(20);
    D_dchirp_rec->SetMarkerSize(0.5);
    D_dchirp_rec->SetMarkerColor(2);
    D_dchirp_rec->SetTitle("");

    // break;
    //
    //###############################################################################################################################
    // Loop over the injected and recovered events
    //###############################################################################################################################
    TChain sim("waveburst");
    TChain mdc("mdc");
    sim.Add(sim_file_name);
    mdc.Add(mdc_file_name);

    double MAXIFAR = 0.0;
    if (sim.GetListOfBranches()->FindObject("ifar")) {
        MAXIFAR = TMath::CeilNint(sim.GetMaximum("ifar") / CYS / TRIALS);
        cout << "Maximum empirically estimated IFAR : " << MAXIFAR << " [years]"
             << endl;
    } else {
        cout << "Missing ifar branch: either use cbc_plots or add it "
                "to wave tree..."
             << endl;
        exit(1);
    }

    TH3F* D_Mtot_rec3 = new TH3F("Distance vs Mtot rec. vs ifar", "", 100,
                                 MINMtot, MAXMtot, 1000, MINDISTANCE / 1000.,
                                 MAXDISTANCE / 1000., 1000, 0., MAXIFAR);
    // D_Mtot_rec->SetName("D_rec");
    D_Mtot_rec3->SetMarkerStyle(20);
    D_Mtot_rec3->SetMarkerSize(0.5);
    D_Mtot_rec3->SetMarkerColor(4);

    int l, m, mtt, cz;
    double mytime[6];
    float factor, distance, mchirp;
    float mass[2];
    float spin[6];
    float chi[3];
    int ifactor;
    double NEVTS;
    mdc.SetBranchAddress("time", mytime);
    mdc.SetBranchAddress("mass", mass);
    mdc.SetBranchAddress("factor", &factor);
    mdc.SetBranchAddress("distance", &distance);
    mdc.SetBranchAddress("mchirp", &mchirp);
    mdc.SetBranchAddress("spin", spin);

    bool write_ascii = false;
    char fname3[2048];
    char line[1024];
    double liveZero;

#ifdef WRITE_ASCII
    write_ascii = true;
#ifdef LIVE_ZERO
    liveZero = LIVE_ZERO;
#else
    cout << "ERROR! Missing zero lag livetime...Please add:" << endl;
    cout << "\"#define LIVE_ZERO XXXXXXXXXX\" where XXXXXXXXXX is the "
            "livetime in seconds to your config/user_pparameters.C file...."
         << endl;
    exit(1);
#endif
#endif
    //	if (write_ascii) {
    sprintf(fname3, "%s/injected_signals.txt", netdir);
    ofstream finj;
    finj.open(fname3, std::ofstream::out);
    sprintf(line,
            "#GPS@L1 mass1 mass2 distance spinx1 spiny1 spinz1 "
            "spinx2 spiny2 spinz2 \n");
    finj << line << endl;
    ofstream* finj_single = new ofstream[nfactor];
    TString xml;
    // cout << fname3 <<endl;
    for (int l = 0; l < nfactor; l++) {
        /* if (Redshift) {
                 cout << fname3 <<endl;
                 xml = XML[l];
                 xml.ReplaceAll(".xml", ".inj.txt");
                 sprintf(fname3, "%s/%s", netdir, xml.Data());
                 cout << fname3 <<endl;
      //   } else {
     //	    cout << fname3 <<endl;
       */
        sprintf(fname3, "%s/injected_signals_%d.txt", netdir, l + 1);
        //	    cout << fname3 <<endl;
        //	    }
        finj_single[l].open(fname3, std::ofstream::out);
        finj_single[l] << line << endl;
    }
    //	}

    for (int g = 0; g < (int)mdc.GetEntries(); g++) {
        mdc.GetEntry(g);
        ifactor = (int)factor - 1;
        if (Redshift) {
            if (PointMasses) {
                mass[0] = SFmasses[ifactor][0];
                mass[1] = SFmasses[ifactor][1];
            }
            mchirp = pow(mass[0] * mass[1], 3. / 5.) /
                     pow(mass[0] + mass[1], 1. / 5.);
        }

        inj_events->Fill(mass[0], mass[1]);
        D_Mtot_inj->Fill(mass[1] + mass[0], distance);
        // D_Mchirp_inj->Fill(mchirp, distance);
        // if(minchi){
        for (int i = 0; i < 3; i++) {
            chi[i] = (spin[i] * mass[0] + spin[i + 3] * mass[1]) /
                     (mass[1] + mass[0]);
        }

        inj_events_spin_mtot->Fill(chi[2], mass[1] + mass[0]);
        //}
        // D_q_inj->Fill(mass[0] / mass[1], distance);

        //		for(int i = 0; i<nfactor; i++){
        //			if(fabs(factor-FACTORS[i])<TOLERANCE){

        if ((ifactor > nfactor - 1) || (ifactor < 0)) {
            cout << "ifactor: " << ifactor << endl;
            exit(1);
        }
        factor_events_inj[ifactor].Fill(mass[0], mass[1]);
        if (minchi) {
            factor_events_spin_mtot_inj[ifactor].Fill(chi[2],
                                                      mass[1] + mass[0]);
        }

        //					break;
        //			}
        //		}
        if (write_ascii) {
            sprintf(line,
                    "%10.3f %3.2f %3.2f %3.2f %3.2f %3.2f "
                    "%3.2f %3.2f %3.2f %3.2f\n",
                    mytime[0], mass[0], mass[1], distance, spin[0], spin[1],
                    spin[2], spin[3], spin[4], spin[5]);
            finj << line << endl;
            finj_single[ifactor] << line << endl;
        }
    }
    if (write_ascii) {
        finj.close();
    }
    NEVTS = 0.0;
    for (int l = 0; l < nfactor; l++) {
        if (write_ascii) {
            finj_single[l].close();
        }
        // Total number of events over all mass bins and factors, since
        // all factors share the same Fiducial volume
        NEVTS += factor_events_inj[l].GetEntries();
    }

    cout << "Injected signals: " << mdc.GetEntries() << endl;
    cout << "Injected signals in histogram factor_events_inj: " << NEVTS
         << endl;
    // char cut[512];
    // sprintf(cut,"rho[%d]>%f && ifar>%f && %s",pp_irho,T_cut,T_ifar,ch2);
    float myifar, ecor, m1, m2, netcc[3], neted, penalty;
    float rho[2];

    float chirp[6];
    float range[2];
    float frequency[2];
    float iSNR[3], sSNR[3];
    sim.SetBranchAddress("mass", mass);
    sim.SetBranchAddress("factor", &factor);
#ifdef OLD_TREE
    sim.SetBranchAddress("distance", &distance);
    sim.SetBranchAddress("mchirp", &mchirp);
#else
    sim.SetBranchAddress("range", range);
    sim.SetBranchAddress("chirp", chirp);
#endif
    sim.SetBranchAddress("rho", rho);
    sim.SetBranchAddress("netcc", netcc);
    sim.SetBranchAddress("neted", &neted);
    sim.SetBranchAddress("ifar", &myifar);
    sim.SetBranchAddress("ecor", &ecor);
    sim.SetBranchAddress("penalty", &penalty);
    sim.SetBranchAddress("time", mytime);
    sim.SetBranchAddress("iSNR", iSNR);
    sim.SetBranchAddress("sSNR", sSNR);
    sim.SetBranchAddress("spin", spin);
    sim.SetBranchAddress("frequency", frequency);

#ifdef HVETO_VETO
    cout << "Adding hveto flags.." << endl;
    UChar_t veto_hveto_L1, veto_hveto_H1, veto_hveto_V1;
    sim.SetBranchAddress("veto_hveto_L1", &veto_hveto_L1);
    sim.SetBranchAddress("veto_hveto_H1", &veto_hveto_H1);
// sim.SetBranchAddress("veto_hveto_V1",&veto_hveto_V1);
#endif

#ifdef CAT3_VETO
    cout << "Adding cat3 flags.." << endl;
    UChar_t veto_cat3_L1, veto_cat3_H1, veto_cat3_V1;
    sim.SetBranchAddress("veto_cat3_L1", &veto_cat3_L1);
    sim.SetBranchAddress("veto_cat3_H1", &veto_cat3_H1);
    // sim.SetBranchAddress("veto_cat3_V1",&veto_cat3_V1);

#endif

    float** volume = new float*[NBINS_mass1];
    float** volume_first_shell = new float*[NBINS_mass1];
    float** radius = new float*[NBINS_mass1];
    float** error_volume = new float*[NBINS_mass1];
    float** error_volume_first_shell = new float*[NBINS_mass1];
    float** error_radius = new float*[NBINS_mass1];

    for (int i = 0; i < NBINS_mass1; i++) {
        volume[i] = new float[NBINS_mass2];
        volume_first_shell[i] = new float[NBINS_mass2];
        radius[i] = new float[NBINS_mass2];
        error_volume[i] = new float[NBINS_mass2];
        error_volume_first_shell[i] = new float[NBINS_mass2];
        error_radius[i] = new float[NBINS_mass2];
        for (int j = 0; j < NBINS_mass2; j++) {
            volume[i][j] = 0.;
            volume_first_shell[i][j] = 0.;
            radius[i][j] = 0.;
            error_volume[i][j] = 0.;
            error_volume_first_shell[i][j] = 0.;
            error_radius[i][j] = 0.;
        }
    }

    // if(minchi){
    float** spin_mtot_volume =
        new float*[NBINS_MTOT +
                   1];  /// WHY " + 1" ? Probably I never manage to make it
                        /// correct and there were overflows...FS
    float** spin_mtot_radius = new float*[NBINS_MTOT + 1];
    float** error_spin_mtot_volume = new float*[NBINS_MTOT + 1];
    float** error_spin_mtot_radius = new float*[NBINS_MTOT + 1];

    for (int i = 0; i < NBINS_MTOT + 1; i++) {
        spin_mtot_volume[i] = new float[NBINS_SPIN + 1];
        spin_mtot_radius[i] = new float[NBINS_SPIN + 1];
        error_spin_mtot_volume[i] = new float[NBINS_SPIN + 1];
        error_spin_mtot_radius[i] = new float[NBINS_SPIN + 1];

        for (int j = 0; j < NBINS_SPIN + 1; j++) {
            spin_mtot_volume[i][j] = 0.;
            error_spin_mtot_volume[i][j] = 0.;
            // volume_first_shell[i][j] = 0.;
            // error_volume_first_shell[i][j] = 0.;
            spin_mtot_radius[i][j] = 0.;
            error_spin_mtot_radius[i][j] = 0.;
        }
    }
    //}
    char fname[1024];
    // if (write_ascii) {
    sprintf(fname, "%s/recovered_signals.txt", netdir);
    ofstream fev;
    fev.open(fname, std::ofstream::out);
    //#ifdef BKG_NTUPLE
    sprintf(line,
            "#GPS@L1			  FAR[Hz] eFAR[Hz] Pval "
            "ePval   factor rho frequency  iSNR  sSNR \n");
    fev << line << endl;
    //#else
    //  fprintf(fev,"#GPS@L1       factor rho frequency iSNR  sSNR
    //  \n");
    //#endif
    ofstream* fev_single = new ofstream[nfactor];
    for (int l = 1; l < nfactor + 1; l++) {
        // if (Redshift) {
        //        xml = XML[l - 1];
        //        xml.ReplaceAll(".xml", ".found.txt");
        //        sprintf(fname, "%s/%s", netdir, xml.Data());
        //} else {
        sprintf(fname, "%s/recovered_signals_%d.txt", netdir, l);
        //}
        fev_single[l - 1].open(fname, std::ofstream::out);
        //#ifdef BKG_NTUPLE
        fev_single[l - 1] << line << endl;
        //#else
        //		    fprintf(fev_single[l-1],"#GPS@L1
        // factor rho  frequency iSNR  sSNR \n");
        //#endif
    }
    // }

    double Vrho[RHO_NBINS], eVrho[RHO_NBINS], Rrho[RHO_NBINS], eRrho[RHO_NBINS],
        Trho[RHO_NBINS];
    for (int i = 0; i < RHO_NBINS; i++) {
        Vrho[i] = 0.;
        eVrho[i] = 0.;
        Rrho[i] = 0.;
        eRrho[i] = 0.;
        Trho[i] = RHO_MIN + i * RHO_BIN;
    }
    double dV, dV1, dV_spin_mtot, nevts, internal_volume;
    int nT;
    int countv = 0;
    int cnt = 0;
    int cnt2 = 0;
    int cntfreq = 0;
    bool bcut = false;

    double liveTot = sim.GetMaximum("ifar");

    double BKG_LIVETIME_yr = liveTot / CYS;
    double BKG_LIVETIME_Myr = BKG_LIVETIME_yr / (1.e6);

    cout.precision(14);
    cout << "Total live time ---> background: " << liveTot << " s" << endl;
    cout << "Total live time ---> background: " << BKG_LIVETIME_yr << " yr"
         << endl;
    cout << "Total live time ---> background: " << BKG_LIVETIME_Myr << " Myr"
         << endl;

    std::vector<double> vdv;
    std::vector<double> vifar;
    // break;
    for (int g = 0; g < (int)sim.GetEntries(); g++) {
        sim.GetEntry(g);
        ifactor = (int)factor - 1;
        if (myifar <= T_ifar * CYS * TRIALS) {
            countv++;
            // cout << countv << " Vetoed for ifar : "<< myifar << "
            // with rho[1]=" <<rho[1]<< endl; cout << countv << "
            // Vetoed for netcc[0] : "<< netcc[0] << endl; cout <<
            // countv << " Vetoed for neted : "<< neted / ecor <<
            // endl; cout << countv << " Vetoed for chi2 : "<<
            // penalty << endl;
            continue;
        }
        // Replaced RHO_MIN with T_cut as in cbc_plots_sim4.C
        if (rho[pp_irho] <= T_cut) {
            countv++;
            continue;
        }
        if (netcc[0] <= T_cor) {
            countv++;

            continue;
        }
        if ((mytime[0] - mytime[nIFO]) < -T_win ||
            (mytime[0] - mytime[nIFO]) > 2 * T_win) {
            countv++;
            continue;
        }  // NOT checking for detector 1 and 2: very small bias...
        if (T_vED > 0) {
            if (neted / ecor >= T_vED) {
                countv++;
                continue;
            }
        }
        if (T_pen > 0) {
            if (penalty <= T_pen) {
                countv++;
                continue;
            }
        }
        bcut = false;
        for (int j = 0; j < nFCUT; j++) {
            if ((frequency[0] > lowFCUT[j]) && (frequency[0] <= highFCUT[j]))
                bcut = true;
        }
        if (bcut) {
            countv++;
            cntfreq++;
            continue;
        }
        // if(factor<11.){continue;}
        if (++cnt % 1000 == 0) {
            cout << cnt << " - ";
        }
        Dt->Fill(mytime[0] - mytime[nIFO]);

        if (Redshift) {
            if (PointMasses) {
                mass[0] = SFmasses[ifactor][0];
                mass[1] = SFmasses[ifactor][1];
            }
            chirp[0] = pow(mass[0] * mass[1], 3. / 5.) /
                       pow(mass[0] + mass[1], 1. / 5.);
        }

#ifdef OLD_TREE

        range[1] = distance;
        chirp[0] = mchirp;

#endif
        if (range[1] == 0.0) {
            continue;
        }
        // D_Mtot_rec->Fill(mass[1] + mass[0], range[1]);
        D_Mtot_rec3->Fill(mass[1] + mass[0], range[1], myifar / TRIALS / CYS);
        // D_Mchirp_rec->Fill(chirp[0], range[1]);
        // D_q_rec->Fill(mass[0] / mass[1], range[1]);
        // D_rec->Fill(range[1]);
        rhocc->Fill(netcc[0], rho[pp_irho]);

        float chi2 = penalty > 0 ? log10(penalty) : 0;
        rho_pf->Fill(chi2, rho[pp_irho]);

        m1 = mass[0];
        m2 = mass[1];
        l = TMath::FloorNint((m2 - min_mass2) / MASS_BIN);
        if ((l >= NBINS_mass2) || (l < 0)) {
            cout << "Underflow/Overflow => Enlarge mass range! l=" << l
                 << " NBINS_mass=" << NBINS_mass2 << " m2=" << m2
                 << " min_mass2=" << min_mass2 << endl;
            exit(1);
        }
        m = TMath::FloorNint((m1 - MIN_MASS) / MASS_BIN);
        if ((m >= NBINS_mass) || (m < 0)) {
            cout << "Underflow/Overflow => Enlarge mass range! m=" << m
                 << " NBINS_mass=" << NBINS_mass << " m1=" << m1
                 << " MIN_MASS=" << MIN_MASS << endl;
            exit(1);
        }
        rec_events->Fill(mass[0], mass[1]);
        D_dchirp_rec->Fill(chirp[0], chirp[0] - chirp[1],
                           range[1]);  // In case of Redshift, this is wrong.
        dchirp_rec->Fill(chirp[0], chirp[0] - chirp[1]);  // same as
                                                          // above
                                                          // if(minchi){
        for (int i = 0; i < 3; i++) {
            chi[i] = (spin[i] * mass[0] + spin[i + 3] * mass[1]) /
                     (mass[1] + mass[0]);
        }

        mtt = TMath::FloorNint((m1 + m2 - MIN_MASS) / MASS_BIN / 2.);
        if ((mtt > NBINS_MTOT) || (mtt < 0)) {
            cout << "mt=" << mtt << " is larger than NBINS_MTOT=" << NBINS_MTOT
                 << " Mtot=" << m1 + m2 << " minMtot=" << MIN_MASS << endl;
            continue;
        }
        cz = TMath::FloorNint((chi[2] - MINCHI) / CHI_BIN);
        if ((cz > NBINS_SPIN) || (cz < 0)) {
            cout << "cz=" << cz << " is larger than NBINS_SPIN=" << NBINS_SPIN
                 << " chi[2]=" << chi[2] << " MINCHI=" << MINCHI << endl;
            continue;
        }
        rec_events_spin_mtot->Fill(chi[2], mass[1] + mass[0]);
        //}
        //	if(factor==FACTORS[nfactor]){factor_events_rec->Fill(mass[1],mass[0]);dV
        //= 1./(pow(factor,3)*factor_events_inj[i]->GetBinContent(m+1,l+1));}break;}
        /////BEWARE! Sorted factors!

        //	for(int i = 0; i<nfactor; i++){
        //					if(fabs(factor-FACTORS[i])<TOLERANCE){

        if (DDistrVolume) {
            dV = shell_volume[ifactor];
            //	    dV1 =
            // 1./(pow(factor,3)*factor_events_inj[i]->GetEntries());
            cnt2++;
        } else if (DDistrUniform) {
            dV = pow(range[1], 2) * shell_volume[ifactor];

        } else if (DDistrChirpMass) {
            dV =
                pow(range[1], 2) * shell_volume[ifactor] *
                pow(chirp[0] / 1.22,
                    5. /
                        6.);  // 1.22=1.4*2^-1./5  as in
                              // https://dcc.ligo.org/DocDB/0119/P1500086/008/cbc_ahope_paper.pdf

        }

        else if (Redshift) {
            dV = VT;  // BEWARE: here it's actually dV*T
            nevts = (double)NINJ[ifactor];
        }
        //	else {cout << "No defined distance distribution?????! Or
        // different from uniform and volume???"<<endl;exit(1);}
        else {
            cout << "No defined distance distribution? "
                    "WARNING: Assuming uniform in volume"
                 << endl;
            DDistrVolume = 1;
            dV = shell_volume[ifactor];
        }
        if (FixedFiducialVolume) {
            // Total number of events over all mass bins and
            // factors, since all factors share the same Fiducial
            // volume
            nevts = NEVTS;
        }

        else {
            nevts = factor_events_inj[ifactor].GetEntries();
        }
        //    internal_volume =
        //    4./3.*TMath::Pi()*pow(minDistance[0]/1000.,3.);
        //   if(INCLUDE_INTERNAL_VOLUME){dV + = internal_volume;}

#ifdef POINTMASSES
        nevts = (double)NINJ[ifactor];  // HERE a control on NINJ
#endif                                  // vector is needed....

        dV1 = dV / nevts;
        if (minchi) {
            nevts = 0.0;
            if (FixedFiducialVolume) {
                for (int i = 0; i < nfactor; i++) {
                    nevts += factor_events_spin_mtot_inj[i].GetBinContent(
                        cz + 1, mtt + 1);
                }
            } else {
                nevts = factor_events_spin_mtot_inj[ifactor].GetBinContent(
                    cz + 1, mtt + 1);
            }
#ifdef POINTMASSES
            nevts = (double)NINJ[ifactor];
#endif
            /// Temporay patch for redshifted mass distributions,
            /// i.e.pow(3. / (4. * TMath::Pi() * Tscale) * vV[i], 1.
            /// / 3.) point-like in source frame and spread over
            /// multiple bins in the detector frame

            dV_spin_mtot = dV / nevts;
        }
        nevts = 0;
        for (int i = 0; i < nfactor; i++) {
            nevts += factor_events_inj[i].GetBinContent(m + 1, l + 1);
        }
#ifdef POINTMASSES
        nevts = (double)NINJ[ifactor];
#endif
        /// Temporay patch for redshifted mass distributions, i.e.
        /// point-like in source frame and spread over multiple bins
        /// in the detector frame
        dV /= nevts;
        //	if(i==nfactor-1){factor_events_rec->Fill(mass[1],mass[0]);
        // if(fabs(factor-INTERNAL_FACTOR)<TOLERANCE){factor_events_rec->Fill(mass[1],mass[0]);
        //	if(fabs(factor-FACTORS[0])<TOLERANCE){
        factor_events_rec->Fill(mass[0], mass[1]);
        //	}
        //		break;
        //	}
        //	}

        nT = TMath::Min(TMath::Floor((rho[pp_irho] - RHO_MIN) / RHO_BIN),
                        (double)RHO_NBINS) +
             1;
        for (int i = 0; i < nT; i++) {
#ifdef VOL_M1M2

            Vrho[i] += dV;
            eVrho[i] += pow(dV, 2);
#else
            Vrho[i] += dV1;
            eVrho[i] += pow(dV1, 2);
#endif
        }
        //	  cout << "rho[1]="<<rho[1]<<"  nT="<<nT<<"
        // Vrho[0]="<<Vrho[0] <<" Vrho[1]="<<Vrho[1]<<"
        // Vrho[2]="<<Vrho[2]<<endl;
        // This is useless now.... since I replace dRHO_MIN with T_cut
        // if(rho[pp_irho]<=T_cut){continue;}
        vdv.push_back(dV1);
        vifar.push_back(myifar);
        if (write_ascii) {
            sprintf(line,
                    "%10.3f  %4.3e %4.3e %4.3e %4.3e  %3.2f	%3.2f "
                    "%3.2f %3.2f %3.2f\n",
                    mytime[2], 1. / myifar,
                    TMath::Sqrt(TMath::Nint(liveTot / myifar)) / liveTot,
                    1.0 - TMath::Exp(-liveZero / myifar),
                    TMath::Sqrt(TMath::Nint(liveTot / myifar)) / liveTot *
                        liveZero * TMath::Exp(-liveZero / myifar),
                    factor, rho[pp_irho], frequency[0], sqrt(iSNR[0] + iSNR[1]),
                    sqrt(sSNR[0] + sSNR[1]));  // BEWARE!!! SNR
                                               // calculation for
                                               // 2-fold network...
            fev << line << endl;
            fev_single[ifactor] << line << endl;
        }

        volume[m][l] += dV;
        error_volume[m][l] += pow(dV, 2);
        if (minchi) {
            spin_mtot_volume[mtt][cz] += dV_spin_mtot;
            // error_spin_mtot_volume[mtt][cz] += pow(dV_spin_mtot,
            // 2);
            error_spin_mtot_volume[mtt][cz] += TMath::Power(dV_spin_mtot, 2.0);
        }
    }
    cout << endl;
    if (write_ascii) {
        fev.close();
        for (int l = 0; l < nfactor; l++) {
            fev_single[l].close();
        }
    }
    cout << "Recovered entries: " << cnt << endl;
    cout << "Recovered entries: " << cnt2 << endl;
    cout << "Recovered entries cut by frequency: " << cntfreq << endl;
    cout << "Recovered entries vetoed: " << countv << endl;
    cout << "dV : " << dV << " dV1 : " << dV1 << endl;

    // break;
    internal_volume = 4. / 3. * TMath::Pi() * pow(minDistance[0] / 1000., 3.);
    if (INCLUDE_INTERNAL_VOLUME) {
        for (int i = 0; i < vdv.size(); i++) {
            if (vdv[i] > 0.0) {
                vdv[i] += internal_volume;
            }
        }
        for (int i = 0; i < RHO_NBINS; i++) {
            if (Vrho[i] > 0.0) {
                Vrho[i] += internal_volume;
            }
        }

        for (int i = 0; i < NBINS_MTOT + 1; i++) {
            for (int j = 0; j < NBINS_SPIN + 1; j++) {
                if (spin_mtot_volume[i][j] > 0.0) {
                    spin_mtot_volume[i][j] += internal_volume;
                }
            }
        }

        for (int i = 0; i < NBINS_mass; i++) {
            for (int j = 0; j < NBINS_mass2; j++) {
                if (volume[i][j] > 0.0) {
                    volume[i][j] += internal_volume;
                }
            }
        }
    }

    Int_t* mindex = new Int_t[vdv.size()];
    TMath::Sort((int)vdv.size(), &vifar[0], mindex, true);
    std::vector<double> vV;
    std::vector<double> veV;
    std::vector<double> vsifar;
    std::vector<double> vseifar;
    std::vector<double> vfar;
    std::vector<double> vefar;

    vV.push_back(vdv[mindex[0]]);
    veV.push_back(pow(vdv[mindex[0]], 2));
    vsifar.push_back(vifar[mindex[0]]);
    vfar.push_back(1. / vifar[mindex[0]]);
    vefar.push_back(TMath::Sqrt(TMath::Nint(liveTot / vifar[mindex[0]])) /
                    liveTot);
    // break;
    int mcount = 0;
    for (int i = 1; i < vdv.size(); i++) {
        if (vifar[mindex[i]] == 0) {
            continue;
        }
        if (vifar[mindex[i]] == vsifar[mcount]) {
            vV[mcount] += vdv[mindex[i]];
            veV[mcount] += pow(vdv[mindex[i]], 2);
        } else {
            vsifar.push_back(vifar[mindex[i]]);
            vfar.push_back(1. / vifar[mindex[i]]);
            vefar.push_back(
                TMath::Sqrt(TMath::Nint(liveTot / vifar[mindex[i]])) / liveTot);
            vseifar.push_back(
                TMath::Sqrt(TMath::Nint(liveTot * vifar[mindex[i]])));
            vV.push_back(vV[mcount] + vdv[mindex[i]]);
            veV.push_back(veV[mcount] + pow(vdv[mindex[i]], 2));
            mcount++;
        }
    }
    cout << "Length of ifar/volume vector: " << vV.size() << endl;
    for (int i = 0; i < vV.size(); i++) {
        veV[i] = TMath::Sqrt(veV[i]);
        vfar[i] *= TRIALS * CYS;
        vefar[i] *= TRIALS * CYS;
        vsifar[i] /= (TRIALS * CYS);
    }

    // break;

    c1->Clear();

    TGraphErrors* gr =
        new TGraphErrors(vV.size(), &vfar[0], &vV[0], &vefar[0], &veV[0]);
    gr->GetYaxis()->SetTitle("Volume [Mpc^{3}]");
    if (Redshift) {
        gr->GetYaxis()->SetTitle("Volume*Time  [Gpc^{3}*yr]");
    }
    gr->GetXaxis()->SetTitle("FAR [yr^{-1}]");
    gr->GetXaxis()->SetTitleOffset(1.3);
    gr->GetYaxis()->SetTitleOffset(1.25);
    gr->GetXaxis()->SetTickLength(0.01);
    gr->GetYaxis()->SetTickLength(0.01);
    gr->GetXaxis()->CenterTitle(kTRUE);
    gr->GetYaxis()->CenterTitle(kTRUE);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetXaxis()->SetLabelFont(42);
    gr->GetYaxis()->SetTitleFont(42);
    gr->GetYaxis()->SetLabelFont(42);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    gr->SetMarkerColor(1);
    gr->SetLineColor(kBlue);
    gr->SetTitle("");

    gr->Draw("ap");

    c1->SetLogx(1);
    sprintf(fname, "%s/ROCV.png", netdir);
    c1->Update();
    c1->SaveAs(fname);

    std::vector<double> vR;
    std::vector<double> veR;
    float Tscale = 1.0;
    int actual_nfactor = 0;
    float FMC = 0.0;
    if (Redshift) {
        for (int i = 0; i < nfactor; i++) {
            FMC += factor_events_inj[i].GetEntries() / NINJ[i];
            if (factor_events_inj[i].GetEntries() > 0) {
                actual_nfactor++;
            }
        }
        Tscale = TMC * FMC / actual_nfactor /
                 CYS;  /// Scaling factor from T MonteCarlo and the
                       /// averaged (over the various MCs) time covered
                       /// by the analysis calculated in Years.
    }
    for (int i = 0; i < vV.size(); i++) {
        vR.push_back(pow(3. / (4. * TMath::Pi() * Tscale) * vV[i], 1. / 3.));
        veR.push_back(pow(3. / (4 * TMath::Pi() * Tscale), 1. / 3.) * 1. / 3 *
                      pow(vV[i], -2. / 3.) * veV[i]);
    }
    cout << "Zero-lag livetime: " << Tscale * 365.25 << " [days]" << endl;
    c1->Clear();

    TGraphErrors* gr2 =
        new TGraphErrors(vR.size(), &vfar[0], &vR[0], &vefar[0], &veR[0]);
    gr2->GetYaxis()->SetTitle("Sensitive Range [Mpc]");
    if (Redshift) {
        gr2->GetYaxis()->SetTitle("Sensitive Range [Gpc]");
    }
    gr2->GetXaxis()->SetTitle("FAR [yr^{-1}]");
    gr2->GetXaxis()->SetTitleOffset(1.3);
    gr2->GetYaxis()->SetTitleOffset(1.25);
    gr2->GetXaxis()->SetTickLength(0.01);
    gr2->GetYaxis()->SetTickLength(0.01);
    gr2->GetXaxis()->CenterTitle(kTRUE);
    gr2->GetYaxis()->CenterTitle(kTRUE);
    gr2->GetXaxis()->SetTitleFont(42);
    gr2->GetXaxis()->SetLabelFont(42);
    gr2->GetYaxis()->SetTitleFont(42);
    gr2->GetYaxis()->SetLabelFont(42);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(1.0);
    gr2->SetMarkerColor(1);
    gr2->SetLineColor(kBlue);
    gr2->SetTitle("");
    gr2->Draw("ap");

    sprintf(fname, "%s/ROC.png", netdir);
    c1->Update();
    c1->SaveAs(fname);

    if (write_ascii) {
        sprintf(fname, "%s/ROC.txt", netdir);
        FILE* frange = fopen(fname, "w");
        fprintf(frange, "#rho FAR[year^{-1}] eFAR range[Gpc] erange\n");
        for (int i = 0; i < vR.size(); i++) {
            fprintf(frange, "%3.3g %4.3g %4.3g %4.4g %4.4g\n", 0.0, vfar[i],
                    vefar[i], vR[i], veR[i]);
        }
        fclose(frange);
    }

    //	gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawRadiusIFAR.C"));
    // gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/AddChip.C"));
    // gROOT->LoadMacro(gSystem->ExpandPathName(
    //   "$HOME_CWB/macros/CreateGraphRadiusIFAR.C"));
    gROOT->LoadMacro(
        gSystem->ExpandPathName("$HOME_CWB/macros/DrawRadiusIFARplots.C"));
    // gROOT->LoadMacro(gSystem->ExpandPathName("$HOME_CWB/macros/DrawRadiusRhoplots.C"));
    // gROOT->LoadMacro(gSystem->ExpandPathName(
    //  "$HOME_CWB/macros/CreateDistanceParplots.C"));
    // DrawRadiusIFAR(vR, veR, vsifar, vseifar, netdir);
    cout << "Volume Density distribution passed as an option: " << opt.Data()
         << endl;
    DrawRadiusIFARplots(sim_file_name, mdc_file_name, shell_volume[0], opt);
    cbcTool.CreateDistanceParplots(sim_file_name, mdc_file_name, netdir, "mtot",
                                   MINMtot, MAXMtot, MAXDISTANCE / 1000, 10,
                                   T_ifar, T_win, nIFO);
    cbcTool.CreateDistanceParplots(sim_file_name, mdc_file_name, netdir,
                                   "mchirp", MINCHIRP, MAXCHIRP,
                                   MAXDISTANCE / 1000, 10, T_ifar, T_win, nIFO);
    cbcTool.CreateDistanceParplots(sim_file_name, mdc_file_name, netdir, "eta",
                                   0.1, 0.25, MAXDISTANCE / 1000, 10, T_ifar,
                                   T_win, nIFO);
    // CreateDistanceParplots(sim_file_name,mdc_file_name,netdir,"eccentricity",-1.,1.,
    // MAXDISTANCE/1000, 10, T_ifar, T_win, nIFO);
    cbcTool.CreateDistanceParplots(sim_file_name, mdc_file_name, netdir, "iota",
                                   -1., 1., MAXDISTANCE / 1000, 10, T_ifar,
                                   T_win, nIFO);
    cbcTool.CreateDistanceParplots(sim_file_name, mdc_file_name, netdir,
                                   "chieff", -1.0, 1.0, MAXDISTANCE / 1000, 10,
                                   T_ifar, T_win, nIFO);
    cbcTool.CreateDistanceParplots(sim_file_name, mdc_file_name, netdir, "chip",
                                   0.0, 1.0, MAXDISTANCE / 1000, 10, T_ifar,
                                   T_win, nIFO);
    cbcTool.CreateDistanceParplots(
        sim_file_name, mdc_file_name, netdir, "distance", MINDISTANCE / 1000,
        MAXDISTANCE / 1000, MAXDISTANCE / 1000, 30, T_ifar, T_win, nIFO);

    // DrawRadiusRhoplots(sim_file_name, mdc_file_name,
    // shell_volume[0],opt);
    //  break;
    for (int i = 0; i < RHO_NBINS; i++) {
        eVrho[i] = TMath::Sqrt(eVrho[i]);
    }
    cout << "Vrho[0] = " << Vrho[0] << " +/- " << eVrho[0] << endl;
    cout << "Vrho[RHO_NBINS-1] = " << Vrho[RHO_NBINS - 1] << " +/- "
         << eVrho[RHO_NBINS - 1] << endl;
    // for(int
    // i=0;i<mdc_entries;i++){inj_events->Fill(im1[i],im2[i]);iMtot[i]=im1[i]+im2[i];}
    inj_events->Draw("colz");

    int MAX_AXIS_Z = inj_events->GetBinContent(inj_events->GetMaximumBin()) + 1;
    inj_events->GetZaxis()->SetRangeUser(0, MAX_AXIS_Z);

    char inj_title[256];
    sprintf(inj_title, "Injected events");

    TPaveText* p_inj =
        new TPaveText(0.325301, 0.926166, 0.767068, 0.997409, "blNDC");
    p_inj->SetBorderSize(0);
    p_inj->SetFillColor(0);
    p_inj->SetTextColor(1);
    p_inj->SetTextFont(32);
    p_inj->SetTextSize(0.045);
    TText* text = p_inj->AddText(inj_title);
    p_inj->Draw();

    // #ifdef MIN_CHI
    // sprintf(fname,"%s/Injected_mass1_mass2_chi_%f_%f.eps",netdir,MIN_CHI,MAX_CHI);
    // #else
    sprintf(fname, "%s/Injected_mass1_mass2.eps", netdir);
    // #endif
    c1->Update();
    c1->SaveAs(fname);
    //###############################################################################################################################
    // dchirp candle plots
    //###############################################################################################################################
    c1->Clear();
    c1->SetLogx(0);
    TH1* hcandle = D_dchirp_rec->Project3D("yx");
    hcandle->SetMarkerSize(0.5);
    hcandle->Draw("CANDLE");
    sprintf(fname, "%s/Dchirp_candle.png", netdir);
    c1->Update();
    c1->SaveAs(fname);
    c1->Clear();
    dchirp_rec->SetMarkerSize(0.5);
    dchirp_rec->Draw("CANDLE");
    sprintf(fname, "%s/Dchirp_candle2.png", netdir);
    c1->Update();
    c1->SaveAs(fname);
    //###############################################################################################################################
    // Delta time
    //###############################################################################################################################
    c1->Clear();
    c1->SetLogy(1);
    Dt->Draw();
    sprintf(fname, "%s/Delta_t.png", netdir);
    c1->Update();
    c1->SaveAs(fname);

    //###############################################################################################################################
    // RHO CC/chi2 PLOTS
    //###############################################################################################################################

    c1->Clear();
    if (pp_rho_log)
        c1->SetLogy(kTRUE);
    else
        c1->SetLogy(kFALSE);
    rhocc->Draw("colz");
    sprintf(fname, "%s/rho_cc.eps", netdir);
    c1->Update();
    c1->SaveAs(fname);

    c1->Clear();
    c1->SetLogx(kFALSE);
    if (pp_rho_log)
        c1->SetLogy(kTRUE);
    else
        c1->SetLogy(kFALSE);
    rho_pf->Draw("colz");
    sprintf(fname, "%s/rho_pf.eps", netdir);
    c1->Update();
    c1->SaveAs(fname);

    //###############################################################################################################################
    // SNR PLOTS
    //###############################################################################################################################
    char sel[1024] = "";
    char newcut[2048] = "";
    char newcut2[2048] = "";
    // char title[1024] = "";
    TString myptitle;
    TString myxtitle;
    TString myytitle;
    c1->Clear();
    c1->SetLogx(true);
    c1->SetLogy(true);
    myptitle = "Number of reconstructed events as a function of injected SNR";
    gStyle->SetOptStat(0);

    myxtitle = "Injected network snr";
    // xtitle = "sqrt(snr[0]^2";
    // for(int i=1;i<nIFO;i++){xtitle += " + snr["+xtitle.Itoa(i,10)+"]^2";}
    // xtitle +=")";
    myytitle = "counts";

    // INJECTED
    sprintf(newcut, "");
    sprintf(sel, "sqrt(pow(snr[%d],2)", 0);
    for (int i = 1; i < nIFO; i++) {
        sprintf(sel, "%s + pow(snr[%d],2)", sel, i);
    }

    sprintf(sel, "%s)>>hist(500)", sel);
    mdc.Draw(sel, "");

    int nmdc = mdc.GetSelectedRows();
    cout << "nmdc : " << nmdc << endl;

    TH2F* htemp = (TH2F*)gPad->GetPrimitive("hist");
    htemp->GetXaxis()->SetTitleOffset(1.35);
    htemp->GetYaxis()->SetTitleOffset(1.50);
    htemp->GetXaxis()->CenterTitle(true);
    htemp->GetYaxis()->CenterTitle(true);
    htemp->GetXaxis()->SetLabelFont(42);
    htemp->GetXaxis()->SetTitleFont(42);
    htemp->GetYaxis()->SetLabelFont(42);
    htemp->GetYaxis()->SetTitleFont(42);
    htemp->GetXaxis()->SetTitle(myxtitle);
    htemp->GetYaxis()->SetTitle(myytitle);
    // htemp->GetXaxis()->SetRangeUser(4e-25,1e-21);
    // htemp->GetXaxis()->SetRangeUser(gTRMDC->GetMinimum("snr"),gTRMDC->GetMaximum("snr"));
    htemp->GetXaxis()->SetRangeUser(
        1, pow(10., TMath::Ceil(TMath::Log10(htemp->GetMaximum()))));
    htemp->GetYaxis()->SetRangeUser(
        1, pow(10., TMath::Ceil(TMath::Log10(htemp->GetMaximum()))));
    htemp->SetLineColor(kBlack);
    htemp->SetLineWidth(3);
    sprintf(newcut,
            "(((time[0]-time[%d])>-%g) || (time[0]-time[%d])<%g) "
            "&& rho[%d]> %g",
            nIFO, T_win, nIFO, 2 * T_win, pp_irho, T_cut);

    // DETECTED iSNR
    sprintf(sel, "sqrt(iSNR[%d]", 0);
    for (int i = 1; i < nIFO; i++) {
        sprintf(sel, "%s + iSNR[%d]", sel, i);
    }
    sprintf(sel, "%s)>>hist2(500)", sel);
    cout << "cut " << newcut << endl;

    // sim.SetFillColor(kBlue);
    // htemp->SetLineColor(kBlue);
    sim.Draw(sel, newcut, "same");
    TH2F* htemp2 = (TH2F*)gPad->GetPrimitive("hist2");
    htemp2->SetFillColor(kRed);
    htemp2->SetFillStyle(3017);
    htemp2->SetLineColor(kRed);
    htemp2->SetLineWidth(2);

    int nwave = sim.GetSelectedRows();
    cout << "nwave : " << nwave << endl;
    sprintf(title, "%s", newcut);

    // DETECTED iSNR after final cuts
    sprintf(newcut2, "%s && %s", newcut, veto_not_vetoed);
    cout << "final cut " << newcut2 << endl;
    TString sel_fin = sel;
    sel_fin.ReplaceAll("hist2", "hist3");

    // sim.SetFillColor(kRed);
    // htemp->SetLineColor(kRed);
    sim.Draw(sel_fin, newcut2, "same");
    TH2F* htemp3 = (TH2F*)gPad->GetPrimitive("hist3");
    htemp3->SetFillColor(kBlue);
    htemp3->SetFillStyle(3017);
    htemp3->SetLineColor(kBlue);
    htemp3->SetLineWidth(2);
    int nwave_final = sim.GetSelectedRows();
    cout << "nwave_final : " << nwave_final << endl;
    sprintf(title, "%s", newcut);

    sprintf(title, "%s", myptitle.Data());
    htemp->SetTitle(title);
    char lab[256];
    leg_snr = new TLegend(0.53,0.70,0.95,0.92, "", "brNDC");
    sprintf(lab, "Injections Average SNR: %g", htemp->GetMean());
    leg_snr->AddEntry("", lab, "a");
    sprintf(lab, "Injected: %i", nmdc);
    leg_snr->AddEntry(htemp, lab, "l");
    sprintf(lab, "Found(minimal cuts): %i", nwave);
    leg_snr->AddEntry(htemp2, lab, "l");
    sprintf(lab, "Found(final cuts): %i", nwave_final);
    leg_snr->AddEntry(htemp3, lab, "l");
    leg_snr->SetFillColor(0);
    leg_snr->Draw();

    sprintf(fname, "%s/Injected_snr_distributions.png", netdir);
    c1->Update();
    c1->SaveAs(fname);

    // Plot estimated network snr vs injected network snr
    sim.SetMarkerStyle(20);
    sim.SetMarkerSize(0.5);
    sim.SetMarkerColor(kRed);

    sprintf(sel, "sqrt(sSNR[%d]", 0);
    for (int i = 1; i < nIFO; i++) {
        sprintf(sel, "%s + sSNR[%d]", sel, i);
    }
    sprintf(sel, "%s):sqrt(iSNR[%d]", sel, 0);
    for (int i = 1; i < nIFO; i++) {
        sprintf(sel, "%s + iSNR[%d]", sel, i);
    }
    sprintf(sel, "%s)>>hist4(500)", sel);
    sim.Draw(sel, newcut, "");
    TH2F* htemp4 = (TH2F*)gPad->GetPrimitive("hist4");
    htemp4->SetTitle("Estimated vs Injected network SNR");
    htemp4->GetXaxis()->SetTitleOffset(1.35);
    htemp4->GetYaxis()->SetTitleOffset(1.50);
    htemp4->GetXaxis()->CenterTitle(true);
    htemp4->GetYaxis()->CenterTitle(true);
    htemp4->GetXaxis()->SetLabelFont(42);
    htemp4->GetXaxis()->SetTitleFont(42);
    htemp4->GetYaxis()->SetLabelFont(42);
    htemp4->GetYaxis()->SetTitleFont(42);
    htemp4->GetXaxis()->SetTitle("Injected network SNR");
    htemp4->GetYaxis()->SetTitle("Estimated network SNR");
    htemp4->GetXaxis()->SetRangeUser(
        5, pow(10., TMath::Ceil(TMath::Log10(htemp4->GetMaximum()))));
    htemp4->GetYaxis()->SetRangeUser(
        5, pow(10., TMath::Ceil(TMath::Log10(htemp4->GetMaximum()))));

    // htemp4->SetMarkerStyle(20);
    // htemp4->SetMarkerSize(0.5);
    // htemp4->SetMarkerColor(kRed);
    sim.SetMarkerColor(kBlue);
    sim.Draw(sel, newcut2, "same");
    // TH2F *htemp5 = (TH2F*)gPad->GetPrimitive("hist4");
    // htemp5->SetMarkerStyle(20);
    // htemp5->SetMarkerSize(0.5);
    // htemp5->SetMarkerColor(kBlue);

    sprintf(fname, "%s/Estimated_snr_vs_Injected_snr.eps", netdir);
    c1->Update();
    c1->SaveAs(fname);

    // Plot relative snr loss

    sprintf(sel, "(sqrt(sSNR[%d]", 0);
    for (int i = 1; i < nIFO; i++) {
        sprintf(sel, "%s + sSNR[%d]", sel, i);
    }
    sprintf(sel, "%s)-sqrt(iSNR[%d]", sel, 0);
    for (int i = 1; i < nIFO; i++) {
        sprintf(sel, "%s + iSNR[%d]", sel, i);
    }
    sprintf(sel, "%s))/sqrt(iSNR[%d]", sel, 0);
    for (int i = 1; i < nIFO; i++) {
        sprintf(sel, "%s + iSNR[%d]", sel, i);
    }
    sprintf(sel, "%s)>>hist5", sel);
    cout << "Selection: " << sel << endl;

    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    c1->SetLogx(false);

    sim.Draw(sel, newcut2);
    TH2F* htemp5 = (TH2F*)gPad->GetPrimitive("hist5");
    htemp5->SetTitle("");
    htemp5->GetXaxis()->SetTitleOffset(1.35);
    htemp5->GetYaxis()->SetTitleOffset(1.50);
    htemp5->GetXaxis()->CenterTitle(true);
    htemp5->GetYaxis()->CenterTitle(true);
    htemp5->GetXaxis()->SetLabelFont(42);
    htemp5->GetXaxis()->SetTitleFont(42);
    htemp5->GetYaxis()->SetLabelFont(42);
    htemp5->GetYaxis()->SetTitleFont(42);
    htemp5->GetXaxis()->SetTitle("(Estimated snr -Injected snr)/Injected snr");
    htemp5->GetYaxis()->SetTitle("Counts");
    sim.GetHistogram()->Fit("gaus");

    sprintf(fname, "%s/Relative_snr_Loss.png", netdir);
    c1->Update();
    c1->SaveAs(fname);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    //###############################################################################################################################
    //  Chirp FAR PLOT
    //###############################################################################################################################
    /*
     c1->Clear();
     c1->SetLogx(0);
     c1->SetLogy(0);
     c1->SetLogz(1);

       sprintf(sel,"chirp[0]:chirp[1]:rho[%d]",pp_irho);

       sim.Draw(sel,newcut2,"goff"); //select 3 columns
       int sel_events = sim.GetSelectedRows();
       double sfar[sel_events];
               double sfar2[sel_events];
           double recMchirp2[sel_events];
               double injMchirp2[sel_events];
       double* srho=sim.GetV3();
       double* recMchirp=sim.GetV1();
       double* injMchirp=sim.GetV2();
       double tmp;
       for(int i=0;i<sel_events;i++){
                                                                         tmp
   = cbcTool.getFAR(srho[i],hc,liveTot)*86400.*365.; if(tmp != 0.0) {sfar[i]
   = tmp;} else{sfar[i] = 1./liveTot;}

       }

       Int_t  *index2 = new Int_t[sel_events];
               TMath::Sort(sel_events,sfar,index2,kFALSE);

               for(int i=0;i<sel_events;i++){

                                                                         sfar2[i]
   = sfar[index2[i]];
                                                                         recMchirp2[i]
   = recMchirp[index2[i]];
                                                                         injMchirp2[i]
   = injMchirp[index2[i]];
               }
               cbcTool.doChirpFARPlot(sel_events, recMchirp2, injMchirp2,
   sfar2, c1, netdir);

   c1->SetLogz(0);
  */
    //###############################################################################################################################
    //  DISTANCE PLOTS
    //###############################################################################################################################
    c1->Clear();

    // break;
    D_Mtot_inj->GetYaxis()->SetRangeUser(10., 3 * MAXDISTANCE);
    D_Mtot_inj->Draw("p");
    D_Mtot_rec3->GetZaxis()->SetRange(0, 1.);
    TH1* t0 = D_Mtot_rec3->Project3D("yx_0");
    t0->SetMarkerColor(kBlue);
    t0->Draw("p same");
    D_Mtot_rec3->GetZaxis()->SetRange(1., 10.);
    TH1* t1 = D_Mtot_rec3->Project3D("yx_1");
    t1->SetMarkerColor(kGreen + 2);
    t1->Draw("p same");
    D_Mtot_rec3->GetZaxis()->SetRange(10., 100.);
    TH1* t10 = D_Mtot_rec3->Project3D("yx_10");
    t10->SetMarkerColor(kCyan);
    t10->Draw("p same");
    D_Mtot_rec3->GetZaxis()->SetRange(100., MAXIFAR);
    TH1* t100 = D_Mtot_rec3->Project3D("yx_100");
    t100->SetMarkerColor(kOrange);
    t100->Draw("p same");
    TLegend* leg_D2 =
        new TLegend(0.579719, 0.156425, 0.85, 0.327409, "", "brNDC");
    leg_D2->AddEntry(t0, "IFAR<1 yr", "p");
    leg_D2->AddEntry(t1, "1 yr < IFAR < 10 yr", "p");
    leg_D2->AddEntry(t10, "10 yr < IFAR < 100 yr", "p");
    leg_D2->AddEntry(t100, "IFAR > 100 yr", "p");

    leg_D2->SetFillColor(0);
    leg_D2->Draw();

    sprintf(fname, "%s/Distance_vs_total_mass_ifar.eps", netdir);
    c1->Update();
    c1->SaveAs(fname);

    //###############################################################################################################################
    // Efficiency
    //###############################################################################################################################

    c1->Clear();
    c1->SetLogy(false);
    TH2F* efficiency = (TH2F*)rec_events->Clone();
    efficiency->SetName("efficiency");
    efficiency->Divide(inj_events);
    efficiency->GetZaxis()->SetRangeUser(0, 1.0);
    efficiency->SetTitle("");

    efficiency->Draw("colz text");

    TExec* ex1 = new TExec("ex1", "gStyle->SetPaintTextFormat(\".2f\");");
    ex1->Draw();

    char eff_title[256];
    sprintf(eff_title, "Efficiency");
    if (minchi) {
        sprintf(eff_title, "%s (%.1f < #chi < %.1f)", eff_title, MINCHI,
                MAXCHI);
    }

    TPaveText* p_eff =
        new TPaveText(0.325301, 0.926166, 0.767068, 0.997409, "blNDC");
    p_eff->SetBorderSize(0);
    p_eff->SetFillColor(0);
    p_eff->SetTextColor(1);
    p_eff->SetTextFont(32);
    p_eff->SetTextSize(0.045);
    text = p_eff->AddText(eff_title);
    p_eff->Draw();

    if (minchi) {
        sprintf(fname, "%s/Efficiency_mass1_mass2_chi_%f_%f.png", netdir,
                MINCHI, MAXCHI);
    } else {
        sprintf(fname, "%s/Efficiency_mass1_mass2.png", netdir);
    }
    c1->Update();
    c1->SaveAs(fname);

    TH2F* efficiency_first_shell = (TH2F*)factor_events_rec->Clone();
    efficiency_first_shell->Divide(&factor_events_inj[nfactor - 1]);

    //###############################################################################################################################
    // Effective radius calculation
    //###############################################################################################################################
    int massbins = 0;
    double V0 = 0.0;
    if (Redshift) {
        Tscale *= 1e-9;
    }  // Normalization to Mpc^3
    for (int j = 0; j < NBINS_mass; j++) {
        for (int k = 0; k < NBINS_mass2; k++) {
            // volume_first_shell[j][k] =
            // efficiency_first_shell->GetBinContent(j+1,k+1);
            if (factor_events_rec->GetBinContent(j + 1, k + 1) != 0.) {
                //  error_volume_first_shell[j][k] =
                //  1./TMath::Sqrt(factor_events_rec->GetBinContent(j+1,k+1));
                massbins++;
            }

            // volume[j][k] = shell_volume*volume[j][k] +
            // volume_internal_sphere*volume_first_shell[j][k];
            // volume[j][k] = shell_volume*volume[j][k];
            V0 += volume[j][k];
            if (error_volume[j][k] != 0.) {
                error_volume[j][k] = TMath::Sqrt(error_volume[j][k]);
            }

            radius[j][k] = pow(3. * volume[j][k] / (4 * TMath::Pi() * Tscale),
                               1. / 3);  // Tscale calculated above: it's either
                                         // 1.0 or if Redshift it's the ratio
                                         // between injected and tot MC signals
            if (volume[j][k] != 0.)
                error_radius[j][k] =
                    (1. / 3) * pow(3. / (4 * TMath::Pi() * Tscale), 1. / 3) *
                    pow(1. / pow(volume[j][k], 2), 1. / 3) * error_volume[j][k];
        }
    }
    cout << "Average Volume at threshold V0 = " << V0 / massbins << "on "
         << massbins << " mass bins" << endl;
    c2->Clear();
    //	break;
    TH2F* h_radius = new TH2F("h_radius", "", NBINS_mass, MIN_MASS, MAX_MASS,
                              NBINS_mass2, min_mass2, max_mass2);
    h_radius->GetXaxis()->SetRangeUser(MIN_plot_mass1, MAX_plot_mass1);
    h_radius->GetYaxis()->SetRangeUser(MIN_plot_mass2, MAX_plot_mass2);
    h_radius->GetXaxis()->SetTitle("Mass 1 (M_{#odot})");
    h_radius->GetYaxis()->SetTitle("Mass 2 (M_{#odot})");
    h_radius->GetZaxis()->SetTitle("Sensitive Distance [Gpc]");
    h_radius->GetXaxis()->SetTitleOffset(1.3);
    h_radius->GetYaxis()->SetTitleOffset(1.25);
    h_radius->GetZaxis()->SetTitleOffset(1.5);
    h_radius->GetXaxis()->CenterTitle(kTRUE);
    h_radius->GetYaxis()->CenterTitle(kTRUE);
    h_radius->GetZaxis()->CenterTitle(kTRUE);
    h_radius->GetXaxis()->SetNdivisions(410);
    h_radius->GetYaxis()->SetNdivisions(410);
    h_radius->GetXaxis()->SetTickLength(0.01);
    h_radius->GetYaxis()->SetTickLength(0.01);
    h_radius->GetZaxis()->SetTickLength(0.01);
    h_radius->GetXaxis()->SetTitleFont(42);
    h_radius->GetXaxis()->SetLabelFont(42);
    h_radius->GetYaxis()->SetTitleFont(42);
    h_radius->GetYaxis()->SetLabelFont(42);
    h_radius->GetZaxis()->SetLabelFont(42);
    h_radius->GetZaxis()->SetLabelSize(0.03);
    h_radius->SetTitle("");
    h_radius->SetMarkerSize(1.5);
    h_radius->SetMarkerColor(kWhite);
    /*TH2F *h_radius_text = new TH2F("h_radius_text", "", NBINS_mass,
    MIN_MASS, MAX_MASS, NBINS_mass2, min_mass2, max_mass2); TH2F
    *h_radius_text2 = new TH2F("h_radius_text2", "", NBINS_mass, MIN_MASS,
                                      MAX_MASS, NBINS_mass2, min_mass2,
    max_mass2);
    */
    for (int i = 1; i <= NBINS_mass; i++) {
        for (int j = 1; j <= NBINS_mass2; j++) {
            h_radius->SetBinContent(i, j, radius[i-1][j-1]*0.001);
            // h_radius_text->SetBinContent(i, j, radius[i - 1][j -
            // 1]);
            //   h_radius_text2->SetBinContent(i, j, error_radius[i
            //   - 1][j - 1]);
            h_radius->SetBinError(i, j, error_radius[i-1][j-1]*0.001);
        }
    }

    h_radius->SetContour(NCont);
    h_radius->SetEntries(1);  // This option needs to be enabled when
                              // filling 2D histogram with SetBinContent
   // h_radius->Draw("colz text colsize=2");  // Option to write error
    h_radius->Draw("colz TEXT");  // Option to write error
                                            // associated to the bin content
                                            // h_radius->Draw("colz text");

    h_radius->GetZaxis()->SetRangeUser(0, MAX_EFFECTIVE_RADIUS / 2./1000);

     gStyle->SetPaintTextFormat("4.2f");
    //TExec* ex2 = new TExec("ex2", "gStyle->SetPaintTextFormat(\".0f\");");
    //ex2->Draw();

    char radius_title[256];

    sprintf(radius_title, "%s : Effective radius (Mpc)", networkname);

    p_radius = new TPaveText(0.325301, 0.926166, 0.767068, 0.997409, "blNDC");
    p_radius->SetBorderSize(0);
    p_radius->SetFillColor(0);
    p_radius->SetTextColor(1);
    p_radius->SetTextFont(32);
    p_radius->SetTextSize(0.045);
    text = p_radius->AddText(radius_title);
    p_radius->Draw();
    float M1, M2;
    Double_t x0, Y0, z0;
    char val[20];
    Int_t nGraphs = 0;
    Int_t TotalConts = 0;
    TLatex Tl;
    TLatex Tl2;
    int point;

#ifdef PLOT_CHIRP

    TH2F* chirp_mass = new TH2F(
        "Chirp_Mass", "", NBINS_mass * 10, MIN_plot_mass1, MAX_plot_mass1 * 1.1,
        NBINS_mass2 * 10, MIN_plot_mass2, MAX_plot_mass2 * 1.1);
    for (Int_t i = 0; i < NBINS_mass * 10; i++) {
        for (Int_t j = 0; j < NBINS_mass2 * 10; j++) {
            M1 = MIN_plot_mass1 +
                 i * (MAX_plot_mass1 - MIN_plot_mass1 * 1.1) / NBINS_mass / 10.;
            M2 = MIN_plot_mass2 +
                 j * (MAX_plot_mass2 - MIN_plot_mass2 * 1.1) / NBINS_mass / 10.;
            chirp_mass->SetBinContent(
                i, j,
                (float)TMath::Power(pow(M1 * M2, 3.) / (M1 + M2), 1. / 5));
        }
    }
    Double_t contours[CONTOURS];
    // contours[0] = 35.;
    for (int i = 0; i < CONTOURS; i++) {
        contours[i] = TMath::Nint((i + 1) * MAXCHIRP / CONTOURS);
    }
    chirp_mass->GetZaxis()->SetRangeUser(0, MAXCHIRP);

    chirp_mass->SetContour(CONTOURS, contours);

    TCanvas* c1t = new TCanvas("c1t", "Contour List", 610, 0, 600, 600);
    c1t->SetTopMargin(0.15);

    // Draw contours as filled regions, and Save points
    chirp_mass->Draw("CONT Z LIST");
    c1t->Update();

    // Get Contours
    TObjArray* conts =
        (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = NULL;
    TGraph* curv = NULL;
    TGraph* gc = NULL;
    nGraphs = 0;
    TotalConts = 0;

    if (conts == NULL) {
        printf("*** No Contours Were Extracted!\n");
        TotalConts = 0;
        return;
    } else {
        TotalConts = conts->GetSize();
    }

    printf("TotalConts = %d\n", TotalConts);

    for (int i = 0; i < TotalConts; i++) {
        contLevel = (TList*)conts->At(i);
        printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
        nGraphs += contLevel->GetSize();
    }

    nGraphs = 0;

    Tl.SetTextSize(0.02);
    // break;
    for (int i = 0; i < TotalConts; i++) {
        contLevel = (TList*)conts->At(i);
        z0 = contours[i];

        printf("Z-Level Passed in as:  Z = %f\n", z0);

        // Get first graph from list on curves on this level
        curv = (TGraph*)contLevel->First();
        for (int j = 0; j < contLevel->GetSize(); j++) {
            point = curv->GetN() - 1;
            curv->GetPoint(point, x0, Y0);
            point--;
            // printf("\tCoordinate Point # %d : x0 = %f  y0 = %f
            // \n",point,x0,y0);

            while (((x0 < MIN_plot_mass1) || (x0 > MAX_plot_mass1) ||
                    (Y0 < MIN_plot_mass2) || (Y0 > MAX_plot_mass2)) &&
                   (point > 0)) {
                curv->GetPoint(point, x0, Y0);
                // printf("\tCoordinate Point # %d : x0 = %f  y0
                // = %f \n",point,x0,y0);
                point--;
            }
            curv->SetLineWidth(1);
            curv->SetLineStyle(3);
            // curv->SetLineColor(2);
            nGraphs++;
            printf("\tGraph: %d  -- %d Elements\n", nGraphs, curv->GetN());
            printf("\tCoordinate Point # %d : x0 = %f  y0 = %f \n", point, x0,
                   Y0);
            // Draw clones of the graphs to avoid deletions in case
            // the 1st
            // pad is redrawn.
            c1->cd();
            gc = (TGraph*)curv->Clone();
            gc->Draw("C");
            // sprintf(val,"#color[2]{#scale[0.9]{#it{M_{chirp}}=%g}}",z0);
            // sprintf(val,"#scale[0.9]{#it{M_{chirp}}=%g}",z0);
            sprintf(val, "#it{M_{c}}=%g", z0);
            // if(i==0){x0 = 110.;y0=10v.;}
            Tl.DrawLatex(x0 - 1., Y0 + 0.5, val);

            // x0=0.;
            // y0=0.;
            // curv = (TGraph*)contLevel->After(curv); // Get Next
            // graph
        }
    }
    c1->Update();

#endif

#ifdef PLOT_MASS_RATIO
    TH2F* mass_ratio = new TH2F(
        "Mass_Ratio", "", NBINS_mass * 10, MIN_plot_mass1, MAX_plot_mass1 * 1.1,
        NBINS_mass2 * 10, MIN_plot_mass2, MAX_plot_mass2 * 1.1);

    for (int i = 0; i < NBINS_mass * 10; i++) {
        for (int j = 0; j < NBINS_mass * 10; j++) {
            M1 = MIN_MASS +
                 i * (MAX_plot_mass1 - MIN_plot_mass1 * 1.1) / NBINS_mass / 10.;
            M2 = MIN_MASS +
                 j * (MAX_plot_mass2 - MIN_plot_mass2 * 1.1) / NBINS_mass / 10.;

            mass_ratio->SetBinContent(i, j, (float)M1 / M2);
        }
    }

    Double_t contours2[5];
    contours2[0] = 1.;
    contours2[1] = 1.5;
    contours2[2] = 2.;
    contours2[3] = 3.;
    contours2[4] = 4.;

    mass_ratio->SetContour(5, contours2);

    TCanvas* c1t2 = new TCanvas("c1t2", "Contour List2", 610, 0, 600, 600);
    c1t2->SetTopMargin(0.15);

    // Draw contours as filled regions, and Save points
    mass_ratio->Draw("CONT Z LIST");
    c1t2->Update();

    // Get Contours
    TObjArray* conts2 =
        (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel2 = NULL;
    TGraph* curv2 = NULL;
    TGraph* gc2 = NULL;

    nGraphs = 0;
    TotalConts = 0;

    if (conts2 == NULL) {
        printf("*** No Contours Were Extracted!\n");
        TotalConts = 0;
        return;
    } else {
        TotalConts = conts2->GetSize();
    }
    printf("TotalConts = %d\n", TotalConts);

    for (int i = 0; i < TotalConts; i++) {
        contLevel2 = (TList*)conts2->At(i);
        printf("Contour %d has %d Graphs\n", i, contLevel2->GetSize());
        nGraphs += contLevel2->GetSize();
    }

    nGraphs = 0;
    Tl2.SetTextSize(0.02);

    for (int i = 0; i < TotalConts; i++) {
        contLevel2 = (TList*)conts2->At(i);
        z0 = contours2[i];

        printf("Z-Level Passed in as:  Z = %f\n", z0);

        // Get first graph from list on curves on this level
        curv2 = (TGraph*)contLevel2->First();

        for (int j = 0; j < contLevel2->GetSize(); j++) {
            point = curv2->GetN() - 1;
            curv2->GetPoint(point, x0, Y0);
            point--;
            printf("\tCoordinate Point # %d : x0 = %f  y0 = %f \n", point, x0,
                   Y0);

            while (((x0 < MIN_plot_mass1) || (x0 > MAX_plot_mass1 - 4.) ||
                    (Y0 < MIN_plot_mass2) || (Y0 > MAX_plot_mass2 - 2)) &&
                   (point > 0)) {
                curv2->GetPoint(point, x0, Y0);
                point--;
            }

            curv2->SetLineWidth(1);
            curv2->SetLineStyle(3);

            nGraphs++;
            printf("\tGraph: %d  -- %d Elements\n", nGraphs, curv2->GetN());
            printf("\tCoordinate Point # %d : x0 = %f  y0 = %f \n", point, x0,
                   Y0);

            // Draw clones of the graphs to avoid deletions in case
            // the 1st pad is redrawn.
            c1->cd();
            gc2 = (TGraph*)curv2->Clone();
            gc2->Draw("C");

            sprintf(val, "#it{q}=%.2g", z0);
            Tl2.DrawLatex(x0, Y0, val);

            // curv2 = (TGraph*)contLevel2->After(curv2); // Get
            // Next graph
        }
    }
    c1->Update();
#endif
    //  #ifdef MIN_CHI
    // sprintf(fname,"%s/Effective_radius_chi_%f_%f.png",netdir,MIN_CHI,MAX_CHI);
    // #else
    sprintf(fname, "%s/Effective_radius.png", netdir);
    // sprintf(fname,"%s/Effective_radius.png",netdir);
    // #endif
    c2->Update();
    c2->SaveAs(fname);

    //###############################################################################################################################
    // Effective radius spin-mtot calculation
    //###############################################################################################################################
    if (minchi) {
        int spin_mtot_bins = 0;
        double V0_spin_mtot = 0.0;
        for (int j = 0; j < NBINS_MTOT; j++) {
            for (int k = 0; k < NBINS_SPIN; k++) {
                //  volume_first_shell[j][k] =
                //  efficiency_first_shell->GetBinContent(j+1,k+1);
                // if(factor_events_rec->GetBinContent(j+1,k+1)
                // != 0.) {
                //	    error_volume_first_shell[j][k] =
                // 1./TMath::Sqrt(factor_events_rec->GetBinContent(j+1,k+1));
                //	    massbins++;
                //}
                if (spin_mtot_volume[j][k] != 0.) {
                    // spin_mtot_volume[j][k] =
                    // shell_volspin_mtot_volume[j][k]; ///
                    // Warning: neglecting the internal
                    // volume...  +
                    // volume_internal_sphere*volume_first_shell[j][k];
                    V0_spin_mtot += spin_mtot_volume[j][k];
                    error_spin_mtot_volume[j][k] = TMath::Sqrt(
                        error_spin_mtot_volume
                            [j]
                            [k]);  /// Warning:
                                   /// neglecting the
                                   /// internal
                                   /// volume...+
                                   /// volume_internal_sphere*volume_first_shell[j][k]*error_volume_first_shell[j][k];

                    spin_mtot_radius[j][k] = pow(3. * spin_mtot_volume[j][k] /
                                                     (4 * TMath::Pi() * Tscale),
                                                 1. / 3);

                    error_spin_mtot_radius[j][k] =
                        (1. / 3) *
                        pow(3. / (4 * TMath::Pi() * Tscale), 1. / 3) *
                        pow(1. / pow(spin_mtot_volume[j][k], 2), 1. / 3) *
                        error_spin_mtot_volume[j][k];
                    spin_mtot_bins++;
                }
                //  cout << j << " " << k << " " <<
                //  spin_mtot_volume[j][k]
                // << " " << spin_mtot_radius[j][k] << endl;
            }
        }
        cout << "Average Spin-Mtot Volume at threshold V0 = "
             << V0_spin_mtot / spin_mtot_bins << endl;
        c2->Clear();

        TH2F* h_spin_mtot_radius =
            new TH2F("h_spin_mtot_radius", "", NBINS_SPIN, MINCHI, MAXCHI,
                     NBINS_MTOT, MIN_MASS, MAX_MASS);
        // h_spin_mtot_radius->GetXaxis()->SetRangeUser(MIN_plot_mass1,MAX_plot_mass1);
        // h_spin_mtot_radius->GetYaxis()->SetRangeUser(MIN_plot_mass2,MAX_plot_mass2);
        h_spin_mtot_radius->GetXaxis()->SetTitle("#chi_{z}");
        h_spin_mtot_radius->GetYaxis()->SetTitle("Total Mass (M_{#odot})");
        h_spin_mtot_radius->GetZaxis()->SetTitle("Sensitive Distance [Gpc]");
        h_spin_mtot_radius->GetXaxis()->SetTitleOffset(1.3);
        h_spin_mtot_radius->GetYaxis()->SetTitleOffset(1.25);
        h_spin_mtot_radius->GetZaxis()->SetTitleOffset(1.5);
        h_spin_mtot_radius->GetXaxis()->CenterTitle(kTRUE);
        h_spin_mtot_radius->GetYaxis()->CenterTitle(kTRUE);
        h_spin_mtot_radius->GetZaxis()->CenterTitle(kTRUE);
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
        h_spin_mtot_radius->SetMarkerColor(kWhite);
        h_spin_mtot_radius->SetMarkerSize(
            1.5);  // to scale the numbers in each cell

        for (int i = 1; i <= NBINS_MTOT; i++) {
            for (int j = 1; j <= NBINS_SPIN; j++) {
                h_spin_mtot_radius->SetBinContent(
                    j, i, spin_mtot_radius[i-1][j-1]*0.001);
                h_spin_mtot_radius->SetBinError(
                    j, i, error_spin_mtot_radius[i-1][j-1]*0.001);
                cout << j << " " << i << " "
                     << h_spin_mtot_radius->GetBinContent(j, i) << endl;
            }
        }

        h_spin_mtot_radius->SetContour(NCont);
        h_spin_mtot_radius->SetEntries(1);  // This option needs to be enabled
                                            // when filling 2D histogram with
                                            // SetBinContent
        //	h_spin_mtot_radius->Draw("colz text0e colsize=1");  //
        // Option to write error associated to the bin content
        h_spin_mtot_radius->Draw(
            "colz text");  // Option to write error associated to the
                           // bin content
                           //	h_spin_mtot_radius->Draw("colz text");
        h_spin_mtot_radius->GetZaxis()->SetRangeUser(0,
                                                     MAX_EFFECTIVE_RADIUS / 2./1000);

     //   TExec* ex3 = new TExec("ex2", "gStyle->SetPaintTextFormat(\".0f\");");
        // TExec *ex2 = new
        // TExec("ex2","gStyle->SetPaintTextFormat(\".2f\");");
       // ex3->Draw();

        sprintf(fname, "%s/Effective_radius_chi.png", netdir);

        c2->Update();
        c2->SaveAs(fname);
    }

    //###############################################################################################################################
    // Calculating Radius vs rho
    //###############################################################################################################################
    for (int i = 0; i < RHO_NBINS; i++) {
//  Vrho[i]*=shell_volume;
// eVrho[i]*=shell_volume;
#ifndef VOL_M1M2
        Vrho[i] /= massbins;
        eVrho[i] /= massbins;
#endif
        Rrho[i] = pow(3. * Vrho[i] / (4 * TMath::Pi() * Tscale), 1. / 3);
        eRrho[i] = pow(3. / (4 * TMath::Pi() * Tscale), 1. / 3) * 1. / 3 *
                   pow(Vrho[i], -2. / 3) * eVrho[i];

        /*if (i % 100 == 0) {
                cout << "Rho bin: " << Trho[i] << " Radius: " << Rrho[i]
                     << " +/- " << eRrho[i] << endl;
        }*/
    }

    TF1* f2 = cbcTool.doRangePlot(RHO_NBINS, Trho, Rrho, eRrho, RHO_MIN, T_cut,
                                  c1, networkname, netdir, write_ascii);

    //###############################################################################################################################
    // Deletions
    //###############################################################################################################################

    cout << "Deletions..." << endl;

    delete[] minMtot, maxMtot, minMChirp, maxMChirp, minDistanceXML,
        maxDistanceXML, minDistance, maxDistance, minRatio, maxRatio,
        shell_volume, FACTORS;
    delete[] waveforms, factor_events_inj, factor_events_spin_mtot_inj;
    /* for (int l = 0; l < nfactor; l++) {
            delete[] finj_single[l];
            delete[] fev_single[l];
     }*/
    // delete[] finj, fev;
    // delete[] finj_single;
    // delete[] fev_single;
    delete c1, c2;  // TCanvas
    delete inj_events, rec_events, factor_events_rec, D_Mtot_inj,
        inj_events_spin_mtot, rec_events_spin_mtot, rhocc, rho_pf, dchirp_rec,
        D_dchirp_rec;  // TH2F
    delete gr, gr2;    // TGraphErrors
    // delete p_inj;

    for (int i = 0; i < NBINS_mass1; i++) {
        delete[] volume[i];
        delete[] volume_first_shell[i];
        delete[] radius[i];
        delete[] error_volume[i];
        delete[] error_volume_first_shell[i];
        delete[] error_radius[i];
    }
    delete[] volume;
    delete[] volume_first_shell;
    delete[] radius;
    delete[] error_volume;
    delete[] error_volume_first_shell;
    delete[] error_radius;

    for (int i = 0; i < NBINS_MTOT + 1; i++) {
        delete[] spin_mtot_volume[i];
        delete[] spin_mtot_radius[i];
        delete[] error_spin_mtot_volume[i];
        delete[] error_spin_mtot_radius[i];
    }
    delete[] spin_mtot_volume;
    delete[] spin_mtot_radius;
    delete[] error_spin_mtot_volume;
    delete[] error_spin_mtot_radius;

    //#ifdef BKG_NTUPLE
    //	cout<<"Mass averaged Visible volume @rho="<<RHO_MIN<<" :
    //"<<Vrho[0]<<" on "<<massbins<< " mass bins" <<endl;
    //	cout<<"OLIVETIME_Myr : "<<OLIVETIME_Myr<<" shell_volume :
    //"<<shell_volume[0]<<endl;
    // break;
    /// ROC Curve
    // cbcTool.doROCPlot(bkg_entries, rho_bkg, index, RHO_BIN, liveTot, f2,
    // c1, netdir, write_ascii);

    //  float rmin = TMath::Max(RHO_MIN,rho_bkg[index[bkg_entries-1]]);

    //	cbcTool.doROCPlotIFAR(sim,final_cut,c1, netdir, write_ascii);
    //#endif

    /*
    #ifdef FAD

    TGraph* grtmp = new TGraph(RHO_NBINS,Trho,Vrho);
      //      f1->SetParameters(500.,-2.3);
            f2->SetParameters(500.,-2.5,0.0);
            grtmp->Fit("f2","R");

            grtmp->SetMarkerStyle(20);
            grtmp->SetMarkerSize(1.0);
            //grtmp->Draw("alp");
     TMultiGraph *multigraph = new TMultiGraph();
      multigraph->Add(grtmp);
      multigraph->Paint("AP");
      multigraph->GetHistogram()->GetXaxis()->SetTitle("#rho");
      multigraph->GetHistogram()->GetXaxis()->SetRangeUser(RHO_MIN,RHO_MAX);
      multigraph->GetHistogram()->GetYaxis()->SetRangeUser(f2->Eval(RHO_MAX),f2->Eval(RHO_MIN));
      multigraph->GetHistogram()->GetYaxis()->SetTitle("Productivity(Mpc^{3}
    Myr)");
      multigraph->GetXaxis()->SetTitleOffset(1.3);
      multigraph->GetYaxis()->SetTitleOffset(1.25);
      multigraph->GetXaxis()->SetTickLength(0.01);
      multigraph->GetYaxis()->SetTickLength(0.01);
      multigraph->GetXaxis()->CenterTitle(kTRUE);
      multigraph->GetYaxis()->CenterTitle(kTRUE);
      multigraph->GetXaxis()->SetTitleFont(42);
      multigraph->GetXaxis()->SetLabelFont(42);
      multigraph->GetYaxis()->SetTitleFont(42);
      multigraph->GetYaxis()->SetLabelFont(42);
      multigraph->GetYaxis()->SetMoreLogLabels(kTRUE);
      multigraph->GetYaxis()->SetNoExponent(kTRUE);

      c1->Clear();
       c1->SetLogy(1);
      gStyle->SetOptFit(kTRUE);

      multigraph->Draw("ALP");
        sprintf(radius_title,"%s : Productivity (%s) ", networkname,
    waveform.Data());
        p_radius->Clear();
      TText *text = p_radius->AddText(radius_title);
       p_radius->Draw();

    sprintf(fname,"%s/Productivity.png",netdir);
     c1->Update();
     c1->SaveAs(fname);



     double VFAD[bkg_entries];
     //double VRHO[bkg_entries];
     double MU[bkg_entries];
     double FAD = 0.0;
     double dFAD = 0.0;
     int cnt2 = 0;
     int lag_cnt = 0;

    double K = OLIVETIME_Myr/BKG_LIVETIME_Myr;
    for (int i = 0;i<bkg_entries;i++){

    //#ifdef OLD_FAD
          dFAD = 1.0/f2->Eval(VRHO[i]);
          FAD += dFAD;
          //FAD = (i+1)*dFAD;

          VFAD[i] = FAD*K;
          //VRHO[i] = rho_bkg[index[i]];
          MU[i] = VFAD[i]/dFAD ;
     //     if ((++cnt2%500==0)||fabs(FAD-0.001) < 0.001)){cout << i << "
    Rho : " << rho_bkg[index[i]] << " Productivity : " << MU[i]/FAD <<" FAD
    : "<<FAD<<" MU : "<<MU[i]<< endl;}
    #ifdef OPEN_LAG

          while((VRHO[i]<=LRHO[lag_cnt])&&lag_cnt < fl_entries){

            // for(int j = -1;j<1;j++)  { cout << i+j << " Rho : " <<
    VRHO[i+j] << " Productivity : " << MU[i+j]/VFAD[i+j] <<" FAD :
    "<<VFAD[i+j]<<" MU : "<<MU[i+j]<< endl;}
                 cout << i+j << " Rho : " << VRHO[i] << " Productivity : "
    << MU[i]/VFAD[i] <<" FAD : "<<VFAD[i]<<" MU : "<<MU[i]<< endl;
                if (i==0){cout << "BEWARE!!! Loudest lag event larger than
    loudest bkg event: "<< LRHO[lag_cnt] <<" "<<VRHO[0]
    <<endl;LFAD[lag_cnt]=FAD;LMU[lag_cnt] = MU[i];lag_cnt++;}
                else {
                      #ifdef OLD_FAD
                      LFAD[lag_cnt] =
    (VFAD[i]-VFAD[i-1])/(VRHO[i]-VRHO[i-1])*(LRHO[lag_cnt]-VRHO[i-1])+VFAD[i-1];
                      LMU[lag_cnt] =
    (MU[i]-MU[i-1])/(VRHO[i]-VRHO[i-1])*(LRHO[lag_cnt]-VRHO[i-1])+MU[i-1];
                      cout << "Lag["<< OPEN_LAG << "]: #"<< lag_cnt <<"
    rho=" << LRHO[lag_cnt] << " FAD="<<LFAD[lag_cnt] << "
    MU="<<LMU[lag_cnt]<<endl;
                      lag_cnt++;
                      #endif
                }
          }

    #endif

         //
          //RhoH->Fill(rho_bkg[index[i]],FAD);

    }

    cout << "Loop ends..."<<endl;

      c1->Clear();
      c1->SetLogy(1);
      gStyle->SetOptFit(kFALSE);
    TGraph* grFAD = new TGraph(bkg_entries,VRHO,VFAD);

    grFAD->GetHistogram()->GetXaxis()->SetRangeUser(RHO_MIN,RHO_MAX);
     // multigraph->GetHistogram()->GetYaxis()->SetRangeUser(0.01,10.0);
      grFAD->GetXaxis()->SetTitle("#rho");
     grFAD->GetYaxis()->SetTitle("FAD (Mpc^{-3} Myr^{-1})");
     grFAD->GetXaxis()->SetTitleOffset(1.3);
     grFAD->GetYaxis()->SetTitleOffset(1.25);
     grFAD->GetXaxis()->SetTickLength(0.01);
     grFAD->GetYaxis()->SetTickLength(0.01);
     grFAD->GetXaxis()->CenterTitle(kTRUE);
     grFAD->GetYaxis()->CenterTitle(kTRUE);
     grFAD->GetXaxis()->SetTitleFont(42);
     grFAD->GetXaxis()->SetLabelFont(42);
     grFAD->GetYaxis()->SetTitleFont(42);
     grFAD->GetYaxis()->SetLabelFont(42);
     grFAD->SetMarkerStyle(20);
     grFAD->SetMarkerSize(1.0);
     grFAD->SetMarkerColor(1);
     grFAD->SetLineColor(kBlue);
     grFAD->SetTitle("");
     grFAD->Draw("alp");

      #ifdef OPEN_LAG
     TGraph* gr_LAG = new TGraph(fl_entries,LRHO,LFAD);
     gr_LAG->SetMarkerStyle(20);
     gr_LAG->SetMarkerSize(1.0);
     gr_LAG->SetMarkerColor(2);
     gr_LAG->Draw("p,SAME");
     #endif

     char leg[128];
     sprintf(leg,"Background");
     leg_FAD = new TLegend(0.707831,0.735751,0.940763,0.897668,"","brNDC");
     leg_FAD->AddEntry(grFAD,leg,"p");

     #ifdef OPEN_LAG

     sprintf(leg,"Events from lag[%i]",OPEN_LAG);
     leg_FAD->AddEntry(gr_LAG,leg,"p");
     #endif
     leg_FAD->SetFillColor(0);
     leg_FAD->Draw();



     char FAD_title[256];
     sprintf(FAD_title,"%s : FAD distribution (%s) ", networkname,
    waveform.Data());



     TPaveText *title1 = new
    TPaveText(0.2941767,0.9274611,0.7359438,0.9974093,"blNDC");
     title1->SetBorderSize(0);
     title1->SetFillColor(0);
     title1->SetTextColor(1);
     title1->SetTextFont(32);
     title1->SetTextSize(0.045);
     TText *text = title1->AddText(FAD_title);
     title1->Draw();



      //RhoH->Draw("LP");
    sprintf(fname,"%s/FAD.eps",netdir);
     c1->Update();
     c1->SaveAs(fname);


     c1->Clear();
      c1->SetLogy(1);
      gStyle->SetOptFit(kFALSE);
    TGraph* grMU = new TGraph(bkg_entries,VRHO,MU);
    grMU->GetHistogram()->GetYaxis()->SetTitle("#mu");
    grMU->GetHistogram()->GetXaxis()->SetTitle("#rho");
    grMU->GetHistogram()->GetXaxis()->SetRangeUser(RHO_MIN,RHO_MAX);

    grMU->GetXaxis()->SetTitleOffset(1.3);
    grMU->GetYaxis()->SetTitleOffset(1.25);
    grMU->GetXaxis()->SetTickLength(0.01);
    grMU->GetYaxis()->SetTickLength(0.01);
    grMU->GetXaxis()->CenterTitle(kTRUE);
    grMU->GetYaxis()->CenterTitle(kTRUE);
    grMU->GetXaxis()->SetTitleFont(42);
    grMU->GetXaxis()->SetLabelFont(42);
    grMU->GetYaxis()->SetTitleFont(42);
    grMU->GetYaxis()->SetLabelFont(42);
    grMU->SetMarkerStyle(20);
    grMU->SetMarkerSize(1.0);
    grMU->SetMarkerColor(1);
    grMU->SetLineColor(kBlue);
    grMU->SetTitle("");
    grMU->Draw("alp");

     #ifdef OPEN_LAG
     TGraph* gr_LAG_MU = new TGraph(fl_entries,LRHO,LMU);
     gr_LAG_MU->SetMarkerStyle(20);
     gr_LAG_MU->SetMarkerSize(1.0);
     gr_LAG_MU->SetMarkerColor(2);
     gr_LAG_MU->Draw("p,SAME");
     #endif


     sprintf(FAD_title,"%s : Expected events per Observation time ",
    networkname);


     title1->Clear();
     TText *text = title1->AddText(FAD_title);
     title1->Draw();
     leg_FAD->Draw();
    sprintf(fname,"%s/MU.eps",netdir);
     c1->Update();
     c1->SaveAs(fname);

    #ifdef OPEN_LAG
     wave.Scan("rho[1]:time[0]:netcc[0]:run",lag_cut,"colsize=12");
     #endif

     #endif
    #endif
    */
    // return vR, veR, vsifar, netdir;
    gSystem->Exit(0);
}
