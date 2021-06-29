//#define PP_64_200_HZ
//#define PP_200_2048_HZ
//#define PP_1600_5000_HZ

//#define CAT2_VETO
//#define HVETO_VETO
//#define CAT3_VETO
//#define PEM_VETO

#define nRHO 150

#define RUN_LABEL "S6A-VSR2"

{
  // ----------------------------------------------------------
  // Declare standard cuts & standard pp_label
  // ----------------------------------------------------------

  user_pp_label = "";

  // ----------------------------------------------------------
  // thresholds declaration
  // ----------------------------------------------------------

  T_cor      = 0.60;       // cc cut
  T_cut      = 3.5;        // rho high frequency cut
  T_CUT      = 3.5;        // rho low frequency cut  (only for simulations)
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


  // ----------------------------------------------------------
  // VETO cuts
  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  // ----------------------------------------------------------

//  const int nvdqf=18;
  const int nvdqf=16;
  dqfile vdqf[nvdqf] = {
         {"L1", "input/S6A_OFFLINE_L1SCIENCE.txt", CWB_CAT0,  0., false, false},
         {"H1", "input/S6A_OFFLINE_H1SCIENCE.txt", CWB_CAT0,  0., false, false},
         {"V1", "input/S6A_OFFLINE_V1SCIENCE.txt", CWB_CAT0,  0., false, false},
         {"L1", "input/S6A_OFFLINE_L1_DQCAT1SEGMENTS.txt", CWB_CAT1,  0., true, false},
         {"H1", "input/S6A_OFFLINE_H1_DQCAT1SEGMENTS.txt", CWB_CAT1,  0., true, false},
         {"V1", "input/S6A_OFFLINE_V1_DQCAT1SEGMENTS.txt", CWB_CAT1,  0., true, false},
         {"L1", "input/S6A_OFFLINE_L1_DQCAT2SEGMENTS.txt", CWB_CAT2,  0., true, false},
         {"H1", "input/S6A_OFFLINE_H1_DQCAT2SEGMENTS.txt", CWB_CAT2,  0., true, false},
         {"V1", "input/S6A_OFFLINE_V1_DQCAT2SEGMENTS.txt", CWB_CAT2,  0., true, false},
         {"L1", "input/S6A_OFFLINE_L1_DQCAT4SEGMENTS.txt", CWB_CAT1,  0., true, false},
         {"H1", "input/S6A_OFFLINE_H1_DQCAT4SEGMENTS.txt", CWB_CAT1,  0., true, false},
         {"V1", "input/S6A_OFFLINE_V1_DQCAT4SEGMENTS.txt", CWB_CAT1,  0., true, false},
//         {"L1", "input/L1-HVETO_COMBINED_VETO_SEGS_BOTH-930960000-4233600_1.2.txt", CWB_HVETO, 0., false, false},
//         {"H1", "input/H1-HVETO_COMBINED_VETO_SEGS_BOTH-930960000-4838400.txt", CWB_HVETO, 0., false, false},
         {"V1", "input/V1_KW_HVETO_S6A.txt", CWB_HVETO, 0., false, false},
         {"L1", "input/S6A_OFFLINE_L1_DQCAT3SEGMENTS.txt", CWB_CAT3,  0., false, false},
         {"H1", "input/S6A_OFFLINE_H1_DQCAT3SEGMENTS.txt", CWB_CAT3,  0., false, false},
         {"V1", "input/S6A_OFFLINE_V1_DQCAT3SEGMENTS.txt", CWB_CAT3,  0., false, false},
        };

}
