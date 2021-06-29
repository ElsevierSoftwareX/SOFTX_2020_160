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


/***************************************************************************
                          CWB::Filter.hh  -  description
                             -------------------
    begin                : Jun 29 2011
    copyright            : (C) 2011 by Gabriele Vedovato
    email                : vedovato@lnl.infn.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CWBFILTER_H
#define CWBFILTER_H


/**
  *@author Gabriele Vedovato
  */

/* mkfilter -- given n, compute recurrence relation
   to implement Butterworth, Bessel or Chebyshev filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 http://www-users.cs.york.ac.uk/~fisher/mkfilter */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <ctype.h>
#include "mkfcomplex.hh"
#ifdef _USE_WAT
#include "wavearray.hh"
#endif

#define global
#define unless(x)   if(!(x))
#define until(x)    while(!(x))

#define VERSION	    "4.6"
#undef	PI
#define PI	    3.14159265358979323846  /* Microsoft C++ does not define M_PI ! */
#define TWOPI	    (2.0 * PI)
#define EPS	    1e-10
#define MAXORDER    10
#define MAXPZ	    512	    /* .ge. 2*MAXORDER, to allow for doubling of poles in BP filter; high values needed for FIR filters */
#define MAXSTRING   256

//typedef void (*proc)();
typedef unsigned int uint;
/*
extern "C"
  { double atof(const char*);
    int atoi(char*);
    void exit(int);
  };
*/
extern char *progname;
extern void readdata(char*, double&, int&, double*, int&, double*);

inline double sqr(double x)	    { return x*x;			       }
inline bool seq(char *s1, char *s2) { return strcmp(s1,s2) == 0;	       }
inline bool onebit(uint m)	    { return (m != 0) && ((m & (m-1)) == 0);   }

//inline double asinh(double x)       // ADA-0.61.32
//  { /* Microsoft C++ does not define */
//    return log(x + sqrt(1.0 + sqr(x)));
//  }

inline double fix(double x)
  { /* nearest integer */
    return (x >= 0.0) ? floor(0.5+x) : -floor(0.5-x);
  }

struct MkfPar {
  uint    index;
  uint    options;
  double  freq1;
  double  freq2;
  double  width1;
  double  width2;
  double  qfactor;
  int     order;
};

#define opt_be 0x00001	/* -Be		Bessel characteristic	       */
#define opt_bu 0x00002	/* -Bu		Butterworth characteristic     */
#define opt_ch 0x00004	/* -Ch		Chebyshev characteristic       */
#define opt_re 0x00008	/* -Re		Resonator		       */
#define opt_pi 0x00010	/* -Pi		proportional-integral	       */
#define opt_lr 0x80000	/* -Lz		Lorentzian Ratio	       */

#define opt_lp 0x00020	/* -Lp		lowpass			       */
#define opt_hp 0x00040	/* -Hp		highpass		       */
#define opt_bp 0x00080	/* -Bp		bandpass		       */
#define opt_bs 0x00100	/* -Bs		bandstop		       */
#define opt_ap 0x00200	/* -Ap		allpass			       */
#define opt_bc 0x40000	/* -Bc		bandcut			       */

#define opt_a  0x00400	/* -a		alpha value		       */
#define opt_l  0x00800	/* -l		just list filter parameters    */
#define opt_o  0x01000	/* -o		order of filter		       */
#define opt_p  0x02000	/* -p		specified poles only	       */
#define opt_w  0x04000	/* -w		don't pre-warp		       */
#define opt_z  0x08000	/* -z		use matched z-transform	       */
#define opt_Z  0x10000	/* -Z		additional zero		       */
#define opt_s  0x20000	/* -o		sample rate		       */

inline void CWBFilterHelp() {
    printf(" \n"
           "-------------------------------------------\n"
           " -Be	Bessel characteristic	       \n"
           " -Bu	Butterworth characteristic     \n"
           " -Ch	Chebyshev characteristic       \n"
           " -Re	Resonator		       \n"
           " -Pi	proportional-integral	       \n"
           " -Lz	Lorentzian Ratio	       \n"
           "-------------------------------------------\n"
           " -Lp	lowpass			       \n"
           " -Hp	highpass		       \n"
           " -Bp	bandpass		       \n"
           " -Bs	bandstop		       \n"
           " -Ap	allpass			       \n"
           " -Bc	bandcut			       \n"
           "-------------------------------------------\n"
           " -a 	alpha value		       \n"
           " -l 	just list filter parameters    \n"
           " -o 	order of filter		       \n"
           " -p 	specified poles only	       \n"
           " -w 	don't pre-warp		       \n"
           " -z 	use matched z-transform	       \n"
           " -Z 	additional zero		       \n"
           " -o 	sample rate		       \n"
           "-------------------------------------------\n"
           "\n");
    };

struct pzrep
  { _complex poles[MAXPZ], zeros[MAXPZ];
    int numpoles, numzeros;
  };


static c_complex bessel_poles[] =
  { /* table produced by /usr/fisher/bessel --	N.B. only one member of each C.Conj. pair is listed */
    { -1.00000000000e+00, 0.00000000000e+00}, { -1.10160133059e+00, 6.36009824757e-01},
    { -1.32267579991e+00, 0.00000000000e+00}, { -1.04740916101e+00, 9.99264436281e-01},
    { -1.37006783055e+00, 4.10249717494e-01}, { -9.95208764350e-01, 1.25710573945e+00},
    { -1.50231627145e+00, 0.00000000000e+00}, { -1.38087732586e+00, 7.17909587627e-01},
    { -9.57676548563e-01, 1.47112432073e+00}, { -1.57149040362e+00, 3.20896374221e-01},
    { -1.38185809760e+00, 9.71471890712e-01}, { -9.30656522947e-01, 1.66186326894e+00},
    { -1.68436817927e+00, 0.00000000000e+00}, { -1.61203876622e+00, 5.89244506931e-01},
    { -1.37890321680e+00, 1.19156677780e+00}, { -9.09867780623e-01, 1.83645135304e+00},
    { -1.75740840040e+00, 2.72867575103e-01}, { -1.63693941813e+00, 8.22795625139e-01},
    { -1.37384121764e+00, 1.38835657588e+00}, { -8.92869718847e-01, 1.99832584364e+00},
    { -1.85660050123e+00, 0.00000000000e+00}, { -1.80717053496e+00, 5.12383730575e-01},
    { -1.65239648458e+00, 1.03138956698e+00}, { -1.36758830979e+00, 1.56773371224e+00},
    { -8.78399276161e-01, 2.14980052431e+00}, { -1.92761969145e+00, 2.41623471082e-01},
    { -1.84219624443e+00, 7.27257597722e-01}, { -1.66181024140e+00, 1.22110021857e+00},
    { -1.36069227838e+00, 1.73350574267e+00}, { -8.65756901707e-01, 2.29260483098e+00},
  };

using namespace std;

namespace CWB {

class Filter {
public: 

  Filter(const int argc, const char *argv[]);
  Filter(char *instr);
  ~Filter();

  void GetZeros(int& numzeros, double *&zeros);
  void GetPoles(int& numpoles, double *&poles);
  double GetGain();
  double Arma(double value);
  void Reset();
  _complex GetFrequencyResponseFunction(double f);
  unsigned int GetOptions() {return options;}
  MkfPar  GetMkfParameters() {return mkf_par;}
  void usage();

  void Help() {
    printf(" \n"
           "-------------------------------------------\n"
           " -Be	Bessel characteristic	       \n"
           " -Bu	Butterworth characteristic     \n"
           " -Ch	Chebyshev characteristic       \n"
           " -Re	Resonator		       \n"
           " -Pi	proportional-integral	       \n"
           " -Lz	Lorentzian Ratio	       \n"
           "-------------------------------------------\n"
           " -Lp	lowpass			       \n"
           " -Hp	highpass		       \n"
           " -Bp	bandpass		       \n"
           " -Bs	bandstop		       \n"
           " -Ap	allpass			       \n"
           " -Bc	bandcut			       \n"
           "-------------------------------------------\n"
           " -a		alpha value		       \n"
           " -l		just list filter parameters    \n"
           " -o		order of filter		       \n"
           " -p		specified poles only	       \n"
           " -w		don't pre-warp		       \n"
           " -z		use matched z-transform	       \n"
           " -Z		additional zero		       \n"
           " -o		sample rate		       \n"
           "-------------------------------------------\n"
           "\n");
    }

#ifdef _USE_WAT
  wavearray<double> data;
#endif

private:

  void Init(const int argc, const char *argv[]);
  void readcmdline(char *argv[]);
  uint decodeoptions(char *s);
  uint optbit(char c);
  double getfarg(char *s);
  int getiarg(char *s);
  void checkoptions();
  void opterror(char *msg, int p1 = 0, int p2 = 0);
  void setdefaults();
  void compute_s(); /* compute S-plane poles for prototype LP filter */
  void choosepole(_complex z);
  void prewarp();
  void normalize();		/* called for trad, not for -Re or -Pi */
  void compute_z_blt(); /* given S-plane poles & zeros, compute Z-plane poles & zeros, by bilinear transform */
  _complex blt(_complex pz);
  void compute_z_mzt(); /* given S-plane poles & zeros, compute Z-plane poles & zeros, by matched z-transform */
  void compute_notch();
  void compute_apres();
  _complex reflect(_complex z);
  void compute_bpres();
  void add_extra_zero();
  void expandpoly(); /* given Z-plane poles & zeros, compute top & bot polynomials in Z, and then recurrence relation */
  void expand(_complex pz[], int npz, _complex coeffs[]);
  void multin(_complex w, int npz, _complex coeffs[]);
  void printresults(char *argv[]);
  void computegain(char *argv[]);
  void printcmdline(char *argv[]);
  void printcoeffs(char *pz, int npz, double coeffs[]);
  void printfilter();
  void printgain(char *str, _complex gain);
  void printrat_s();	/* print S-plane poles and zeros */
  void printrat_z();	/* print Z-plane poles and zeros */
  void printpz(_complex *pzvec, int num);
  void printrecurrence(); /* given (real) Z-plane poles & zeros, compute & print recurrence relation */
  void prcomplex(_complex z);

  pzrep splane, zplane;
  int order;
  double raw_alpha1, raw_alpha2, raw_alphaz;
  _complex dc_gain, fc_gain, hf_gain;
  uint options;
  double warped_alpha1, warped_alpha2, chebrip, qfactor, width1, width2;
  double sample_rate;
  bool infq;
  uint polemask;
  double xcoeffs[MAXPZ+1], ycoeffs[MAXPZ+1];
	double gain;
	double xbp[MAXPZ+1];
	double ybp[MAXPZ+1];
  MkfPar mkf_par;
  const char* mkf_argv[16];
};

} // end namespace

#ifdef _USE_WAT
// put operator
wavearray<double>& operator >> (CWB::Filter& filter, wavearray<double>& x);
// get operator
CWB::Filter& operator >> (wavearray<double>& x, CWB::Filter& filter);
#endif

#endif
