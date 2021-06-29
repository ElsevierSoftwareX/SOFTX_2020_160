{
  //!REGRESSION
  // Config Plugin for linear/bilinear ranking analysis for one detector

  //---------------------------------------------------------------------
  // MAIN REGRESSION RANK
  //---------------------------------------------------------------------

  bool SAVE_TREE                        = true;                 // if true then results are saved to tree
  bool SAVE_ASCII                       = true;                 // if true then results are saved to ascii file

  TString REGRESSION_RANK               = "LINEAR";             // regression rank [LINEAR/BILINEAR]
  int  CUT_LOW_FREQ                     = 32;                   // if > 0 the data in [0:CUT_LOW_FREQ] are set to 0

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS FOR H1
  //---------------------------------------------------------------------


  TString FRLIST_WITNESS             	= "input/H1_RDS_R_L1_937473000_947017000.frl";
  TString CHLIST_LINEAR     		= "input/S6B_H1_Linear_Regression_Channel_List.txt";
  TString CHLIST_WITNESS	       	= "input/S6B_H1_Witness_Channel_List.txt";

  TString CHNAME_LINEAR                 = "H0:PEM-COIL_MAGX";

  int     RESAMPLING_INDEX              = 11;           // resample target/witness channels to pow(2,RESAMPLING_INDEX)

  double  WFLOW                         = 0;            // lower frequency used to made the bilinear channels
  double  WFHIGH                        = 10;           // high  frequency used to made the bilinear channels

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS
  //---------------------------------------------------------------------

  double LAYER_WIDTH                    = 1.0;          // frequency layer resolution used in the regression analysis
  double fPOWERLINE                     = 60.0;         // powerline frequency (Hz)
  int    lPOWERLINE                     = 2;            // low power line harmonic (lPOWERLINE*fPOWERLINE)
  int    hPOWERLINE                     = 2;            // high power line harmonic (hPOWERLINE*fPOWERLINE)
  double FWIDTH                         = 5;            // frequency width used by regression

  // PARAMETERS FOR LINEAR REGRESSION

  int    L_NFILTER                      = 5;            // half-size length of a unit filter (setFilter)
  double L_APPLY_THRESHOLD              = 0.2;          // threshold used in apply
  double L_SOLVE_THRESHOLD              = 0.0;          // eigenvalue threshold (solve)
  double L_SOLVE_NEIGEN_PER_LAYER       = 0;            // number of selected eigenvalues (solve)
  char   L_SOLVE_REGULATOR              = 'h';          // regulator (solve)

  // PARAMETERS FOR BILINEAR REGRESSION

  int    B_NFILTER                      = 5;            // half-size length of a unit filter (setFilter)
  double B_APPLY_THRESHOLD              = 0.2;          // threshold used in apply
  double B_SOLVE_THRESHOLD              = 0.0;          // eigenvalue threshold (solve)
  double B_SOLVE_NEIGEN_PER_LAYER       = 0;            // number of selected eigenvalues (solve)
  char   B_SOLVE_REGULATOR              = 'h';          // regulator (solve)

}
