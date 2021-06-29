{
//!DATA_CONDITIONING
// Implements the 2G frequency cuts in pixel selection stage : configPlugin File Example

  cout << "CWB_Plugin_fCuts_Config.C" << endl;

  // frequency user fcut definition

  int nFCUT=10;	// nFCUT is used in CWB_Plugin_fCuts.C plugin (max nFCUT is FCUT_MAX)
  // entry format : {"IFO", freq_low, freq_high, "list of levels used by fcut"}
  // freq range removed : [freq_low(Hz), freq_high(Hz)]
  // list of levels used by fcut : "#lev1,#lev2,...,#levN"
  fcut fCUT[nFCUT] = {
                    {"L1" , 58,  62, "8,9,10"},
                    {"L1" , 89,  96, "8,9,10"},
                    {"L1" ,102, 123, "8,9,10"},

                    {"H1" , 58,  62, "8,9,10"},
                    {"H1" , 83,  88, "8,9,10"},
                    {"H1" ,102, 122, "8,9,10"},

                    {"V1" , 48,  52, "8,9,10"},
                    {"V1" , 82,  89, "8,9,10"},
                    {"V1" , 91,  96, "8,9,10"},
                    {"V1" ,106, 112, "8,9,10"}
                  };

  // copy to FCUT array (used by CWB_Plugin_fCuts.C plugin)
  for(int j=0;j<nFCUT;j++) FCUT[j]=fCUT[j];

}
