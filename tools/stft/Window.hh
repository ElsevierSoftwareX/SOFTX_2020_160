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
                          Window.hh  -  description
                             -------------------
    begin                : Wen Dec 21 2011
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

#ifndef CWBWINDOW_H
#define CWBWINDOW_H


/**
  *@author Gabriele Vedovato
  */

/* Reference:
 * Oppenheim, A. V. and Schafer, R. W., ``Discrete-Time Signal Processing.''
 * ISBN 0-13-216771-9, Prentice-Hall, Inc., 1989, pp. 447--448
 *
 * http://www.mathworks.com/access/helpdesk/help/toolbox/signal/signal.shtml
 */

/* this is a c++ porting of the window package in the "Sonic Flow" tool box
   The homepage of Sonic Flow is at http://sonicflow.sourceforge.net/.
 */

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

using namespace std;

namespace CWB {

class Window {
public: 
	Window(char* formula, unsigned n, double fParameter=0);
	~Window();

  TString GetVersion(char c='s') {
    if(c=='s') return TString("window101");
    if(c=='d') return TString("101");
    else       return TString("window-1.0.1");
  }

  double GetValue(unsigned i);
  double GetSize() {return size;}

private:

  double* window;
  unsigned size;
  void Normalize (double* out_window, unsigned n);
  double fParameter;

  /* Compute a modified Bartlett-Hann window.
   * Like Bartlett, Hann, and Hamming windows, this window has a
   * mainlobe at the origin and asymptotically decaying sidelobes on
   * both sides. It is a linear combination of weighted Bartlett and
   * Hann windows with near sidelobes lower than both Bartlett and Hann
   * and with far sidelobes lower than both Bartlett and Hamming
   * windows. The mainlobe width of the modified Bartlett-Hann window is
   * not increased relative to either Bartlett or Hann window mainlobes.
   */
  void
  barthann (double* out_window, unsigned n);

  /* Compute a Bartlett window.
   * The Bartlett window is very similar to a triangular window as
   * returned by the triangular function. The Bartlett
   * window always ends with zeros at samples 0 and n -1, however,
   * while the triangular window is nonzero at those points. For length
   * odd, the center n-2 points of bartlett(,n) are equivalent
   * to those of sf_window_triangular(,n-2).
   */
  void
  bartlett (double* out_window, unsigned n);

  /* Compute a Blackman window
   * Blackman windows have slightly wider central lobes and less sideband
   * leakage than equivalent length Hamming and Hann windows.
   */
  void
  blackman (double* out_window, unsigned n);

  /* Compute a minimum 4-term Blackman-harris window.
   * The window is minimum in the sense that its maximum sidelobes are
   * minimized.
   */
  void
  blackmanharris (double* out_window, unsigned n);

  /* Compute a Bohman window.
   * A Bohman window is the convolution of two half-duration cosine
   * lobes. In the time domain, it is the product of a triangular window
   * and a single cycle of a cosine with a term added to set the first
   * derivative to zero at the boundary. Bohman windows fall off as 1/w^4.
   */
  void
  bohman (double* out_window, unsigned n);

  /* [UNIMPLEMENTED] */
  /*  void */
  /*  chebyshev (double * out_window, unsigned n, */
  /*			 double  r = 100.0); */

  /* Compute a Flat Top weighted window
   * Flat Top windows have very low passband ripple (< 0.01 dB) and are
   * used primarily for calibration purposes. Their bandwidth is
   * approximately 2.5 times wider than a Hann window.
   */
  void
  flattop (double* out_window, unsigned n);

  /* Compute a Gaussian window where alpha >= 2 is the reciprocal of the
   * standard deviation. The width of the window is inversely related to
   * the value of alpha; a larger value of alpha produces a more narrow
   * window. Default value for alpha : 2.5
   */
  void
  gauss (double* out_window, unsigned n,
		     double  alpha);

  /* Compute a Hamming window
   */
  void
  hamming (double* out_window, unsigned n);

  /* Compute a Hann (Hanning) window
   */
  void
  hann (double* out_window, unsigned n);

  /* [UNIMPLEMENTED] */
  /*  void */
  /*  kaiser (double * out_window, unsigned n, */
  /*		      SF_Real beta = 0.5); */

  /* Compute a minimum 4-term Blackman-Harris window, as defined by Nuttall
   * The window is minimum in the sense that its maximum sidelobes are
   * minimized. The coefficients for this window differ from the
   * Blackman-Harris window coefficients computed with
   * blackmanharris and produce slightly lower sidelobes.
   */
  void
  nuttall (double* out_window, unsigned n);

  /* [UNIMPLEMENTED] */
  /*  void */
  /*  parzen (double* out_window, unsigned n); */

  /* Compute a rectangular window
   */
  void
  rectangular (double* out_window, unsigned n);

  /* [UNIMPLEMENTED] */
  /* void */
  /* saramaki (double* out_window, unsigned n); */

  /* [UNIMPLEMENTED] */
  /*  void */
  /*  transversal (double* out_window, unsigned n); */

  /* Compute a triangular window.
   */
  void
  triangular (double* out_window, unsigned n);

  /* Compute a Tukey window
   * Tukey windows are cosine-tapered windows. r is the ratio of taper to
   * constant sections and is between 0 and 1. r <= 0 is a
   * rectangular window and r>= 1 is a hann window.
   * default value for r : 0.5
   */
  void
  tuckey (double* out_window, unsigned n,
		     double  r);

  /* Compute a welch window.
   */
  void
  welch (double* out_window, unsigned n);
};

} // end namespace

#endif

