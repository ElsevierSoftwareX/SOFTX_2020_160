{
  // Definitions of the search bins used for the GW150914 analysis
  #include "../../cwb/postproduction/O1/GW150914_search_bins.hh"

  // Merge Label
  TString merge_label = "M1.V_cat2LH_hvetoLH_cat3LH";

  // Definition of the FAR file names
  TString far_dir         = "../ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a/report/dump/";
  TString far_unmodeled   = far_dir+"FARvsRHO_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a_unmodeled.txt";
  TString far_constrained = far_dir+"FARvsRHO_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a_constrained.txt";
  TString far_chirp       = far_dir+"FARvsRHO_ER8b_12Sep20Oct_C01_BKG_LF_rMRA_run1ato10a_chirp.txt";

  char cmd[1024];

  // add ifar parameter to the wave root file for the unmodeled bin
  sprintf(cmd,"${CWB_SCRIPTS}/cwb_setifar.csh %s \"%s\" %s unmodeled",
          merge_label.Data(),unmodeled.GetTitle(),far_unmodeled.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);

  // add ifar parameter to the wave root file for the constrained bin
  sprintf(cmd,"${CWB_SCRIPTS}/cwb_setifar.csh %s.S_unmodeled \"%s\" %s constrained",
          merge_label.Data(),constrained.GetTitle(),far_constrained.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);

  // add ifar parameter to the wave root file for the chirp bin
  sprintf(cmd,"${CWB_SCRIPTS}/cwb_setifar.csh %s.S_unmodeled.S_constrained \"%s\" %s chirp",
          merge_label.Data(),chirp.GetTitle(),far_chirp.Data());
  cout << cmd << endl;
  gSystem->Exec(cmd);
  
  exit(0);
}
