//
// Config File for Bico Jobs
// Author : Gabriele Vedovato


{

#define SRATE 16384.
//#define BSIZE (16384*20)
#define BSIZE (16384*80)

/*
//  
#define BCM_ORDER 1
//#define BCM_X_NUM_SLICES 2048
//#define BCM_Y_NUM_SLICES 2048
#define BCM_X_NUM_SLICES 256
#define BCM_Y_NUM_SLICES 256
#define BCM_X_SLICE_INDEX 0
#define BCM_Y_SLICE_INDEX 0
*/

/*
//  
#define BCM_ORDER 1
#define BCM_X_NUM_SLICES 256
#define BCM_Y_NUM_SLICES 256
#define BCM_X_SLICE_INDEX 5
#define BCM_Y_SLICE_INDEX 5
*/

/*
// ---------------------------
// line 60.0
// ---------------------------
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 512
#define BCM_Y_NUM_SLICES 2048
#define BCM_X_SLICE_INDEX 3
#define BCM_Y_SLICE_INDEX 0
*/

/*
// ---------------------------
// line 120.0
// ---------------------------
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 512
#define BCM_Y_NUM_SLICES 2048
#define BCM_X_SLICE_INDEX 7
#define BCM_Y_SLICE_INDEX 0
*/


// ---------------------------
// line 180
// ---------------------------
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 1024
#define BCM_Y_NUM_SLICES 2048
#define BCM_X_SLICE_INDEX 22  
#define BCM_Y_SLICE_INDEX 0


/*
// ---------------------------
// lines 180, 184.5, 185.385, 189.387
// ---------------------------
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 512
#define BCM_Y_NUM_SLICES 2048
#define BCM_X_SLICE_INDEX 11
#define BCM_Y_SLICE_INDEX 0
*/

/*
// ---------------------------
// lines 180, 184.5, 185.385, 189.387
// ---------------------------
#define BCM_ORDER 2
#define BCM_Y_NUM_SLICES 2048
#define BCM_Y_NUM_SLICES 512
#define BCM_X_SLICE_INDEX 11
#define BCM_Y_SLICE_INDEX 0
*/

/*
// ---------------------------
// line 199.523
// ---------------------------
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 512
#define BCM_Y_NUM_SLICES 2048
#define BCM_X_SLICE_INDEX 12
#define BCM_Y_SLICE_INDEX 0
*/

/*
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 512
#define BCM_Y_NUM_SLICES 1024
#define BCM_X_SLICE_INDEX 15
#define BCM_Y_SLICE_INDEX 2
*/

/*
#define BCM_ORDER 1
#define BCM_X_NUM_SLICES 1024
#define BCM_Y_NUM_SLICES 1024
#define BCM_X_SLICE_INDEX 4
#define BCM_Y_SLICE_INDEX 4
*/

/*
// ---------------------------
// line 210.0
// ---------------------------
#define BCM_ORDER 2
#define BCM_X_NUM_SLICES 512
#define BCM_Y_NUM_SLICES 2048
#define BCM_X_SLICE_INDEX 13
#define BCM_Y_SLICE_INDEX 0
*/



//#define FRLIST_NAME "input/S6A_R3_H1_LDAS_C02_L2.frames"
//#define CHNAME "H1:LDAS-STRAIN"
//#define CHNAME "L1:LDAS-STRAIN"
//#define CHNAME "V1:h_16384Hz"

//#define SEGLIST "input/S6B_H1SCIENCE_EXCLUDECAT1CAT4.2c.txt"

#define FRLIST_NAME_1 "input/H1_LDAS_C02_L2.frl"
#define FRLIST_NAME_2 "input/H1_RDS_R_L1.frl"

#define CHNAME_1 "H1:LDAS-STRAIN"
//#define CHNAME_1 "H1:LSC-DARM_CTRL"
//#define CHNAME_1 "H1:LSC-DARM_ERR"

//#define CHNAME_2 "H0:PEM-LVEA2_V1"
//#define CHNAME_2 "H0:PEM-LSC1_MAGZ"
//#define CHNAME_2 "H0:PEM-HAM6_ACCX"
//#define CHNAME_2 "H0:PEM-COIL_MAGX"
//#define CHNAME_2 "H0:PEM-ISCT10_ACCX"
//#define CHNAME_2 "H0:PEM-LVEA_SEISX"
//#define CHNAME_2 "H0:PEM-BSC1_MAG1X"
//#define CHNAME_2 "H0:PEM-BSC1_MAG1Y"
//#define CHNAME_2 "H0:PEM-BSC1_MAG1Z"

//#define CHNAME_2 "H1:SUS-ETMX_COIL_LL"
//#define CHNAME_2 "H1:SUS-ETMX_COIL_LR"
//#define CHNAME_2 "H1:SUS-ETMX_COIL_UL"
//#define CHNAME_2 "H1:SUS-ETMX_COIL_UR"

#define CHNAME_2 "H1:SUS-ITMX_COIL_LL"
//#define CHNAME_2 "H1:SUS-ITMX_COIL_LR"
//#define CHNAME_2 "H1:SUS-ITMX_COIL_UL"
//#define CHNAME_2 "H1:SUS-ITMX_COIL_UR"


//#define START 942449664
//#define START (942449664+8000)
#define START (942449664-500)

#define REBIN 5
//#define REBIN 1
#define NLOOP 100
//#define NLOOP 10
#define BATCH

#define BLC_THRESHOLD 0.5

#define ODIR_NAME "plots_bico_H1_SUS-ITMX_COIL_LL"

#define GRAPHID 0

#define SAVE

}
