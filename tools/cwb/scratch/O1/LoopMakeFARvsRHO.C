{

  #define WAVE_FILE "merge/wave_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a.M1.C_rho_gt7.V_cat2LH_hvetoLH_cat3LH.root"
  #define BKG_LIVETIME  2145888576000.

  gROOT->LoadMacro("MakeFARvsRHO_C.so");

  TString odir = "report/dump/";
  TString ofile_unmodeled    = odir+"FARvsRHO_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a_unmodeled.txt";
  TString ofile_eunmodeled   = odir+"FARvsRHO_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a_eunmodeled.txt";
  TString ofile_constrained  = odir+"FARvsRHO_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a_constrained.txt";
  TString ofile_econstrained = odir+"FARvsRHO_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a_econstrained.txt";
  TString ofile_chirp        = odir+"FARvsRHO_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a_chirp.txt";

  MakeFARvsRHO("eunmodeled",   WAVE_FILE, BKG_LIVETIME, ofile_eunmodeled);
  MakeFARvsRHO("unmodeled",    WAVE_FILE, BKG_LIVETIME, ofile_unmodeled);
  MakeFARvsRHO("eunmodeled",   WAVE_FILE, BKG_LIVETIME, ofile_eunmodeled);
  MakeFARvsRHO("constrained",  WAVE_FILE, BKG_LIVETIME, ofile_constrained);
  MakeFARvsRHO("econstrained", WAVE_FILE, BKG_LIVETIME, ofile_econstrained);
  MakeFARvsRHO("chirp",        WAVE_FILE, BKG_LIVETIME, ofile_chirp);

  exit(0);
}
