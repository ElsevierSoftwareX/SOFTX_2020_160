{

//#define BBH_LH
#define IMBHB_LH

//#define BBH_LV
//#define IMBHB_LV

//#define BBH_HV
//#define IMBHB_HV

// -------------------------------------------------------
// O3_K15_C01_LH_BBH_BKG_xrun2
// -------------------------------------------------------

#ifdef BBH_LH

  #include <O3/SEARCHES/OFFLINE/BBH/LH/PP_Cuts.hh>

  #define IFILE_BKG     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/BBH/LH/BKG/O3_K15_C01_LH_BBH_BKG_xrun2/merge/wave_O3_K15_C01_LH_BBH_BKG_xrun2.M1.V_hvetoLH.C_rho1_gt_6.root"
  #define IFILE_SIM     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/BBH/LH/SIM/O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1/merge/wave_O3_K15_C01_LH_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.V_hvetoLH.root"  

  #define BKG_LIVE_TIME	35266962627.00	// non-zero lags : 77934 lags - 35266962627.00 sec = 408182.44 days = 1118.3 years

  // definition of the basic selection cuts used in O3b run
  TCut bin_cut_standard = bin1_cut;
  TCut bin_cut_basic = TCut("bin_cut_basic",(dqveto+lveto_cut).GetTitle());

  #define PP_CUTS "                                                       \n\
  #                                                                       \n\
  # name          cmp     mean    sigma   min     max     enabled_tune    \n\
  #                                                                       \n\
  netcc[0]        >       0.7     0.1     0.5     0.9     1               \n\
  netcc[2]        >       0.7     0.1     0.5     0.9     1               \n\
  norm            >       2.5     1.0     1.0     5.0     1               \n\
  Qveto[0]        >       0.25    0.1     0.0     0.5     1               \n\
  log10(penalty)  <       0.2     0.1     0.1     0.4     1               \n\
  frequency[0]    >       60.0    10.0    40.0    80.0    1               \n\
  frequency[0]    <       300.0   50.0    150.0   350.0   1               \n\
  chirp[1]        >       1.0     3.0     0.5     5.5     1               \n\
  chirp[1]        <       70.0    10.0    50.0    90.0    1               \n\
  #Qveto[2]       >       0.0     1.0     0.0     5.0     1               \n\
  "

#endif

// -------------------------------------------------------
// O3_K15_C01_LH_IMBHB_BKG_xrun2
// -------------------------------------------------------

#ifdef IMBHB_LH

  #include <O3/SEARCHES/OFFLINE/IMBHB/LH/PP_Cuts.hh>

  #define IFILE_BKG     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/IMBHB/LH/BKG/O3_K15_C01_LH_IMBHB_BKG_xrun2/merge/wave_O3_K15_C01_LH_IMBHB_BKG_xrun2.M1.V_hvetoLH.C_rho1_gt_6.root"
  #define IFILE_SIM     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/IMBHB/LH/SIM/O3_K15_C01_LH_IMBHB_SIM_NR-MIX_19MBins_xrun1/merge/wave_O3_K15_C01_LH_IMBHB_SIM_NR-MIX_19MBins_xrun1.M1.V_hvetoLH.root"

  #define BKG_LIVE_TIME	35266962627.00	// non-zero lags : 77934 lags - 35266962627.00 sec = 408182.44 days = 1118.3 years

  // definition of the basic selection cuts used in O3b run
  TCut bin_cut_standard = bin1_cut;
  TCut bin_cut_basic = TCut("bin_cut_basic",(dqveto_cut+TF_cut+"chirp[1]>-100").GetTitle());

  #define PP_CUTS "									\n\
  #											\n\
  # name          		cmp     mean    sigma   min     max     enabled_tune	\n\
  #											\n\
  netcc[0] 			> 	0.8 	0.1 	0.5 	0.9 	1		\n\
  netcc[2] 			> 	0.7 	0.1 	0.5 	0.9 	1		\n\
  norm 				> 	4.0 	1.0 	2.0 	6.0 	1		\n\
  Qveto[0] 			> 	0.09 	0.1 	0.0 	0.3 	1		\n\
  log10(penalty) 		< 	0.4 	0.1 	0.1 	0.6 	1		\n\
  frequency[0] 			> 	24.0 	20.0 	16.0 	32.0 	1		\n\
  frequency[0] 			< 	100.0 	40.0 	50.0 	150.0 	1		\n\
  abs(chirp[1]) 		> 	10.0 	3.0 	5.0 	15.0 	1		\n\
  abs(chirp[1]/Qveto[0]) 	> 	15.0 	5.0 	10.0 	15.0	1		\n\
  #Qveto[2]                     >       0.0     1.0     0.0     5.0     1		\n\
  "

#endif


// -------------------------------------------------------
// O3_K15_C01_LV_BBH_BKG_xrun1
// -------------------------------------------------------

#ifdef BBH_LV

  #include <O3/SEARCHES/OFFLINE/BBH/LV/PP_Cuts.hh>

  #define IFILE_BKG     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/BBH/LV/BKG/O3_K15_C01_LV_BBH_BKG_xrun1/merge/wave_O3_K15_C01_LV_BBH_BKG_xrun1.M1.V_hvetoL_userLV.C_rho1_gt_6.root"
  #define IFILE_SIM     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/BBH/LV/SIM/O3_K15_C01_LV_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1/merge/wave_O3_K15_C01_LV_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.V_hvetoL.root"  

  #define BKG_LIVE_TIME	4101136694.00	// non-zero lags : 5994 lags - 4101136694.00 sec = 47466.86 days = 130.0 years

  // definition of the basic selection cuts used in O3b run
  TCut bin_cut_standard = bin1_cut;
  TCut bin_cut_basic = TCut("bin_cut_basic",(dqveto+lveto_cut).GetTitle());

  #define PP_CUTS "                                                       \n\
  #                                                                       \n\
  # name          cmp     mean    sigma   min     max     enabled_tune    \n\
  #                                                                       \n\
  netcc[0]        >       0.5     0.1     0.3     0.7     1               \n\
  netcc[2]        >       0.5     0.1     0.3     0.7     1               \n\
  norm            >       2.5     1.0     1.0     5.0     1               \n\
  Qveto[0]        >       0.25    0.1     0.0     0.5     1               \n\
  log10(penalty)  <       0.2     0.1     0.1     0.4     1               \n\
  frequency[0]    >       60.0    10.0    40.0    80.0    1               \n\
  frequency[0]    <       300.0   50.0    150.0   350.0   1               \n\
  chirp[1]        >       1.0     3.0     0.5     5.5     1               \n\
  chirp[1]        <       70.0    10.0    50.0    90.0    1               \n\
  #Qveto[2]       >       0.0     1.0     0.0     5.0     1               \n\
  "

#endif

// -------------------------------------------------------
// O3_K15_C01_LV_IMBHB_BKG_xrun1
// -------------------------------------------------------

#ifdef IMBHB_LV

  #include <O3/SEARCHES/OFFLINE/IMBHB/LV/PP_Cuts.hh>

  #define IFILE_BKG     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/IMBHB/LV/BKG/O3_K15_C01_LV_IMBHB_BKG_xrun1/merge/wave_O3_K15_C01_LV_IMBHB_BKG_xrun1.M1.V_hvetoL_userLV.C_rho1_gt_6.root"
  #define IFILE_SIM     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/IMBHB/LV/SIM/O3_K15_C01_LV_IMBHB_SIM_NR-MIX_19MBins_xrun1/merge/wave_O3_K15_C01_LV_IMBHB_SIM_NR-MIX_19MBins_xrun1.M1.V_hvetoL.root"

  #define BKG_LIVE_TIME	4101136694.00	// non-zero lags : 5994 lags - 4101136694.00 sec = 47466.86 days = 130.0 years

  // definition of the basic selection cuts used in O3b run
  TCut bin_cut_standard = bin1_cut;
  TCut bin_cut_basic = TCut("bin_cut_basic",(dqveto_cut+TF_cut+"chirp[1]>-100").GetTitle());

  #define PP_CUTS "									\n\
  #											\n\
  # name          		cmp     mean    sigma   min     max     enabled_tune	\n\
  #											\n\
  netcc[0] 			> 	0.5 	0.1 	0.3 	0.7 	1		\n\
  netcc[2] 			> 	0.5 	0.1 	0.3 	0.7 	1		\n\
  norm 				> 	4.0 	1.0 	2.0 	6.0 	1		\n\
  Qveto[0] 			> 	0.09 	0.1 	0.0 	0.3 	1		\n\
  log10(penalty) 		< 	0.4 	0.1 	0.1 	0.6 	1		\n\
  frequency[0] 			> 	24.0 	20.0 	16.0 	32.0 	1		\n\
  frequency[0] 			< 	100.0 	40.0 	50.0 	150.0 	1		\n\
  abs(chirp[1]) 		> 	10.0 	3.0 	5.0 	15.0 	1		\n\
  abs(chirp[1]/Qveto[0]) 	> 	15.0 	5.0 	10.0 	15.0	1		\n\
  #Qveto[2]                     >       0.0     1.0     0.0     5.0     1		\n\
  "

#endif


// -------------------------------------------------------
// O3_K15_C01_HV_BBH_BKG_xrun1
// -------------------------------------------------------

#ifdef BBH_HV

  #include <O3/SEARCHES/OFFLINE/BBH/HV/PP_Cuts.hh>

  #define IFILE_BKG     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/BBH/HV/BKG/O3_K15_C01_HV_BBH_BKG_xrun1/merge/wave_O3_K15_C01_HV_BBH_BKG_xrun1.M1.V_hvetoH_userHV.C_rho1_gt_6.root"
  #define IFILE_SIM     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/BBH/HV/SIM/O3_K15_C01_HV_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1/merge/wave_O3_K15_C01_HV_BBH_SIM_BBH_BROAD_ISOTROPIC_xrun1.M1.V_hvetoH.root"  

  #define BKG_LIVE_TIME	3795059499.00	// non-zero lags : 5994 lags - 3795059499.00 sec = 43924.30 days = 120.3 years

  // definition of the basic selection cuts used in O3b run
  TCut bin_cut_standard = bin1_cut;
  TCut bin_cut_basic = TCut("bin_cut_basic",(dqveto+lveto_cut).GetTitle());

  #define PP_CUTS "                                                       \n\
  #                                                                       \n\
  # name          cmp     mean    sigma   min     max     enabled_tune    \n\
  #                                                                       \n\
  netcc[0]        >       0.5     0.1     0.3     0.7     1               \n\
  netcc[2]        >       0.5     0.1     0.3     0.7     1               \n\
  norm            >       2.5     1.0     1.0     5.0     1               \n\
  Qveto[0]        >       0.25    0.1     0.0     0.5     1               \n\
  log10(penalty)  <       0.2     0.1     0.1     0.4     1               \n\
  frequency[0]    >       60.0    10.0    40.0    80.0    1               \n\
  frequency[0]    <       300.0   50.0    150.0   350.0   1               \n\
  chirp[1]        >       1.0     3.0     0.5     5.5     1               \n\
  chirp[1]        <       70.0    10.0    50.0    90.0    1               \n\
  #Qveto[2]       >       0.0     1.0     0.0     5.0     1               \n\
  "

#endif

// -------------------------------------------------------
// O3_K15_C01_HV_IMBHB_BKG_xrun1
// -------------------------------------------------------

#ifdef IMBHB_HV

  #include <O3/SEARCHES/OFFLINE/IMBHB/HV/PP_Cuts.hh>

  #define IFILE_BKG     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/IMBHB/HV/BKG/O3_K15_C01_HV_IMBHB_BKG_xrun1/merge/wave_O3_K15_C01_HV_IMBHB_BKG_xrun1.M1.V_hvetoH_userHV.C_rho1_gt_6.root"
  #define IFILE_SIM     "/home/gabriele.vedovato/O3/SEARCHES/OFFLINE/IMBHB/HV/SIM/O3_K15_C01_HV_IMBHB_SIM_NR-MIX_19MBins_xrun1/merge/wave_O3_K15_C01_HV_IMBHB_SIM_NR-MIX_19MBins_xrun1.M1.V_hvetoH.root"

  #define BKG_LIVE_TIME	3795059499.00	// non-zero lags : 5994 lags - 3795059499.00 sec = 43924.30 days = 120.3 years

  // definition of the basic selection cuts used in O3b run
  TCut bin_cut_standard = bin1_cut;
  TCut bin_cut_basic = TCut("bin_cut_basic",(dqveto_cut+TF_cut+"chirp[1]>-100").GetTitle());

  #define PP_CUTS "									\n\
  #											\n\
  # name          		cmp     mean    sigma   min     max     enabled_tune	\n\
  #											\n\
  netcc[0] 			> 	0.5 	0.1 	0.3 	0.7 	1		\n\
  netcc[2] 			> 	0.5 	0.1 	0.3 	0.7 	1		\n\
  norm 				> 	4.0 	1.0 	2.0 	6.0 	1		\n\
  Qveto[0] 			> 	0.09 	0.1 	0.0 	0.3 	1		\n\
  log10(penalty) 		< 	0.4 	0.1 	0.1 	0.6 	1		\n\
  frequency[0] 			> 	24.0 	20.0 	16.0 	32.0 	1		\n\
  frequency[0] 			< 	100.0 	40.0 	50.0 	150.0 	1		\n\
  abs(chirp[1]) 		> 	10.0 	3.0 	5.0 	15.0 	1		\n\
  abs(chirp[1]/Qveto[0]) 	> 	15.0 	5.0 	10.0 	15.0	1		\n\
  #Qveto[2]                     >       0.0     1.0     0.0     5.0     1		\n\
  "

#endif

}
