#define nRHO 150

#define RUN_LABEL "ADV H1 VS H2H3"

#define CWB_SIMPLOT_USE_FACTOR

{
  // ----------------------------------------------------------
  // Declare standard cuts & standard pp_label
  // ----------------------------------------------------------

  user_pp_label = "";

  // ----------------------------------------------------------
  // thresholds declaration
  // ----------------------------------------------------------

  T_cor      = 0.50;       // cc cut
  T_cut      = 2.0;        // rho high frequency cut
  T_CUT      = 2.0;        // rho low frequency cut  (only for simulations)
  T_vED      = 0.0;        // vED threshold
  T_pen      = 0.0;        // penalty threshold
  T_out      = 3.5;        // output threshold

  T_cut_veto = 6.9;        // rho high frequency cut
  T_CUT_veto = 6.9;        // rho low frequency cut

  hours      = 24;         // bin size in hours for rate vs time plot 

  // obsolete parameters
  T_acor     = 0.0;        // output threshold
  T_hrss     = 0.0;        // penalty threshold
  i_hrss1    = 0;
  i_hrss2    = 1;
  T_rms      = 0.;         // penalty threshold

  ratio_hrss  = 1;          //CORRECTION FREQUENCY
  pp_freq_offset = 0.;      //CORRECTION FREQUENCY  (SIM)

}
