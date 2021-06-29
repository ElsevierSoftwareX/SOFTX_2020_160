{
  //#define PP_64_200_HZ
  //#define PP_200_2048_HZ
  //#define PP_1600_5000_HZ

  #define RUN_LABEL "ADV 2G"

  // ----------------------------------------------------------
  // Declare standard cuts & standard pp_label
  // ----------------------------------------------------------

  user_pp_label = "";

  // ----------------------------------------------------------
  // thresholds declaration
  // ----------------------------------------------------------

  T_cor      = 0.5;        // cc cut
  //T_cut      = 6.0;        // rho high frequency cut
  T_cut      = 5.0;        // rho high frequency cut
  T_CUT      = 0.0;        // rho low frequency cut  (only for simulations)
  T_vED      = 0.0;        // vED threshold

  T_pen      = 0.0;        // penalty threshold
  T_out      = 0.0;        // output threshold

  T_cut_veto = 6.9;        // rho high frequency cut
  T_CUT_veto = 6.9;        // rho low frequency cut

  T_win      = 10.0;

  hours      = 24;         // bin size in hours for rate vs time plot

  pp_max_nloudest_list = 25; // maximum number of loudest event in the report list

  pp_irho = 1;

  lowFCUT[nFCUT] = 0.;  highFCUT[nFCUT] = 32.; nFCUT++;
  lowFCUT[nFCUT] = 1024.;  highFCUT[nFCUT] = 8196.; nFCUT++;

  pp_jet_benckmark = -1;
  pp_mem_benckmark = -1;

}
