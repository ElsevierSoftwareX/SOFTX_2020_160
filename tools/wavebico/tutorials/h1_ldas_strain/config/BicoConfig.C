/*
# Copyright (C) 2019 Gabriele Vedovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


//
// Config File for wavebico 
// Author : Gabriele Vedovato


{

#define SRATE 16384.		// sample rate of the first channel
#define BSIZE (16384*80)	// FFT size, define the final bi-spectrum resolution 
#define BCM_ORDER 2		// order = [1/2] -> [coherence/bicoherence]  

// line 120.0
#define BCM_X_NUM_SLICES 512	// number of slices used for channel 1 
                                // SLICE_WIDTH_X = SRATE/BCM_X_NUM_SLICES/2
#define BCM_Y_NUM_SLICES 2048   // number of slices used for channel 2 
                                // SLICE_WIDTH_Y = SRATE/BCM_X_NUM_SLICES/2
#define BCM_X_SLICE_INDEX 7     // channel 1 slice displayed in the final bi-spectrum 
                                // BEGIN_SLICE_FREQ_X = BCM_X_SLICE_INDEX*SLICE_WIDTH_X
                                // END_SLICE_FREQ_X   = BEGIN_SLICE_FREQ_X+DFREQ_X
#define BCM_Y_SLICE_INDEX 0     // channel 2 slice displayed in the final bi-spectrum 
                                // BEGIN_SLICE_FREQ_Y = BCM_Y_SLICE_INDEX*SLICE_WIDTH_Y
                                // END_SLICE_FREQ_Y   = BEGIN_SLICE_FREQ_Y+SLICE_WIDTH_Y

#define FRLIST_NAME_1 "input/H1_LDAS_C02_L2.frl"  // frame list of the first channel
#define FRLIST_NAME_2 "input/H1_RDS_R_L1.frl"	  // frame list of the second channel

#define CHNAME_1 "H1:LDAS-STRAIN"		// name of the first channel

#define CHNAME_2 "H1:SUS-ITMX_COIL_LL"		// name of the second channel
//#define CHNAME_2 "H1:SUS-ITMX_COIL_LR"
//#define CHNAME_2 "H1:SUS-ITMX_COIL_UL"
//#define CHNAME_2 "H1:SUS-ITMX_COIL_UR"

#define START (942449664+8000)			// gps start of bicoherence

//#define REBIN 1
#define REBIN 5					// final rebinning of bi-spectrum 	
                                                // WARNING !!! low value of REBIN can hidden bi-coherent pixels 
#define NLOOP 100				// number of FFT used for final bi-spectrum
                                                // to increase SNR one must increase NLOOP
                                                // WARNING !!! - noisy data period can destroy bi-coherence 
#define BATCH					// if defined the bi-spectrum is displayed  
#define GRAPHID 0				// if GRAPHID = [0/1] spectrum of channel 1/2 are displayed 
                                                // together with the bi-spectrum 

#define ODIR_NAME "output"			// output directory to save bi-spectra
#define SAVE					// if define the bi-spectrum is saved in the output dir

}
