{

// --------------------------------------------------------------
// 1G vs 2G-ISRA 
// --------------------------------------------------------------

nROC=2;

REP_WWW = "https://www.atlas.aei.uni-hannover.de/~vedovato/LSC/reports/";

ROC_TITLE = "ROC <br> S6A-L1H1V1 <br> <font color=\"blue\">1G</font> vs <font color=\"red\">2G-ISRA</font>";

ROC_DIR = "report/dump/ROC_1G_2GISRA"; 

ROC_WWW = REP_WWW+"/S6A_R4_SIM_BRST_L1H1V1_2G_run43/dump/ROC_1G_2GISRA";

// --------------------------------------------------------------
// NAME reports 
// --------------------------------------------------------------

ROC_NAME[0] = "1G";
ROC_NAME[1] = "2G-ISRA";

// --------------------------------------------------------------
// SIM/BKG line type
// --------------------------------------------------------------

LINE_COLOR[0] = 4;
LINE_COLOR[1] = 2;

LINE_STYLE[0] = 1;
LINE_STYLE[1] = 1;

// --------------------------------------------------------------
// SIM/BKG reports dirs
// --------------------------------------------------------------

SIM_DIR[0] = "/home/vedovato/tst/Re_S6A_R3_LF_BRST_L1H1V1_1G_Manual/report/postprod/M1.R_rsearch_i0cc70_i1rho40_vED40_win01_freq64_2048/data/";
SIM_DIR[1] = "/home/vedovato/S6/coherent/offline/SIMULATION/S6A_R4_SIM_BRST_L1H1V1_2G_run43/report/postprod/M1.C_cc2_gt_0d45.V_hvetoLHV_cat3LHV.R_ISRA_hveto_cat3_i0cc70_i1rho56_win01_freq64_2048/data";

BKG_DIR[0] = "/home/waveburst/Official_results/coherent/offline/PRODUCTION/S6A_R3_BKG_LF_L1H1V1_1GM_wat539/report/postprod/M1.FIX.V_hvetoLHV_cat3LHV.R_rsearch_hveto_cat3_i0cc70_i1rho60_vED40_freq64_2048/data";
BKG_DIR[1] = "/home/vedovato/S6/coherent/offline/PRODUCTION/S6A_R4_BKG_L1H1V1_2G_run17/report/postprod/M1.C_cc2_gt_0d45.V_hvetoLHV_cat3LHV.R_ISRA_hveto_cat3_i0cc70_i1rho50_freq64_2048/data";

// --------------------------------------------------------------
// SIM/BKG reports : www links
// --------------------------------------------------------------

SIM_WWW[0] = REP_WWW+"/Re_S6A_R3_LF_BRST_L1H1V1_1G_Manual/postprod/M1.R_rsearch_i0cc70_i1rho40_vED40_win01_freq64_2048/";
SIM_WWW[1] = REP_WWW+"/S6A_R4_SIM_BRST_L1H1V1_2G_run43/postprod/M1.C_cc2_gt_0d45.V_hvetoLHV_cat3LHV.R_ISRA_hveto_cat3_i0cc70_i1rho56_win01_freq64_2048/";

BKG_WWW[0] = "https://www.atlas.aei.uni-hannover.de/~waveburst/LSC/reports/Reports_1G_wat539/S6A_R3_BKG_LF_L1H1V1_1GM_wat539/postprod/M1.FIX.V_hvetoLHV_cat3LHV.R_rsearch_hveto_cat3_i0cc70_i1rho60_vED40_freq64_2048/";
BKG_WWW[1] = REP_WWW+"/S6A_R4_BKG_L1H1V1_2G_run17/postprod/M1.C_cc2_gt_0d45.V_hvetoLHV_cat3LHV.R_ISRA_hveto_cat3_i0cc70_i1rho50_freq64_2048/";

}
