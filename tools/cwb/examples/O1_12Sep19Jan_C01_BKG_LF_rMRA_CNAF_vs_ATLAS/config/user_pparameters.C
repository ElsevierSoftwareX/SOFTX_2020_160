#define RUN_LABEL "O1 LF UnModeled C01 : 2015-10-20 13:32:03 UTC Tue To 2016-01-01 03:00:00 UTC Fri"
#define CAT2_VETO
#define CAT3_VETO
//#define HVETO_VETO
{
  T_cor      = 0.70;       // cc cut
  T_cut      = 0.0;        // rho high frequency cut
  T_out	     = 0.0;
  hours      = 1.;         // bin size in hours for rate vs time plot

  pp_irho = 0;
  pp_inetcc =  0;
  pp_rho_min = 5;
  pp_rho_max = 20;
//  pp_rho_max = 250;
//  pp_rho_log = true;

  pp_max_nloudest_list = 100; 

  lowFCUT[nFCUT] = 0;    highFCUT[nFCUT] = 32; nFCUT++;
  lowFCUT[nFCUT] = 992;  highFCUT[nFCUT] = 1024; nFCUT++;
 
//  pp_jet_benckmark = -1;
//  pp_mem_benckmark = -1;


  // ----------------------------------------------------------
  // VETO cuts
  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  // ----------------------------------------------------------

  const int nvdqf=4;
  dqfile vdqf[nvdqf] = {
         {"L1" ,"/home/vedovato/O1//DQvetos/O1_20Oct01Jan_C01/L1Cat2.txt",   CWB_CAT2, 0., true,   false},
         {"H1" ,"/home/vedovato/O1//DQvetos/O1_20Oct01Jan_C01/H1Cat2.txt",   CWB_CAT2, 0., true,   false}
         {"L1" ,"/home/vedovato/O1//DQvetos/O1_20Oct01Jan_C01/L1Cat4.txt",   CWB_CAT3, 0., false,  false},
         {"H1" ,"/home/vedovato/O1//DQvetos/O1_20Oct01Jan_C01/H1Cat4.txt",   CWB_CAT3, 0., false,  false}
         };

}
