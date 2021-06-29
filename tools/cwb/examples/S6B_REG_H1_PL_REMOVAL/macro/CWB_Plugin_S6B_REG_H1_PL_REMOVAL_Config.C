{

  //---------------------------------------------------------------------
  // MAIN
  //---------------------------------------------------------------------

  bool APPLY_REGRESSION			= true;
  bool SAVE_FRAME			= true;
  bool SAVE_PLOT			= true;
  bool SAVE_PLOT_INDATA			= false;
  bool CUT_0_32_HZ			= true;

  TString OPLOT_DIR       		= "oplots";
  TString OFRAME_DIR     		= "oframes";

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS FOR H1
  //---------------------------------------------------------------------

  TString H1_FRLIST_WITNESS		= "input/H1_RDS_R_L1_94240000_942500000.frl";
  TString H1_WITNESS_CHANNEL_LIST	= "input/WitnessChannels_H1.lst";

  //---------------------------------------------------------------------
  // REGRESSION PARAMETERS
  //---------------------------------------------------------------------

  int    nFILTER			= 5;
  double LAYER_WIDTH			= 1.0;
  int    nPL				= 5;       // number of processed power lines
  double FWIDTH				= 5;       // frequency width used by regression

  // parameters to clean power lines

  double APPLY_THRESHOLD		= 0.2;
  double SOLVE_THRESHOLD		= 0.0;
  double SOLVE_NEIGEN_PER_LAYER		= 0;
  char   SOLVE_REGULATOR		= 'h';

}
