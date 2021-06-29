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
                          CWB::Window.cc  -  description
                             -------------------
    begin                : Sun May 2 2004
    copyright            : (C) 2004 by Gabriele Vedovato
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

#include "Window.hh"

CWB::Window::Window(char* formula, unsigned n, double fParameter) {

  if (n <= 0) {cout << "CWB::Window::Window - window size must be > 0" << endl;exit(1);}
  if (formula == NULL) {cout << "CWB::Window::Window - window formula not defined" << endl;exit(1);}

  window = new double[n];
  size = n;

  bool  iswindow = false;

  // initialize window

  if (strcmp(formula,"barthann") == 0)       {barthann (window, n);iswindow=true;}
  if (strcmp(formula,"bartlett") == 0)       {bartlett (window, n);iswindow=true;}
  if (strcmp(formula,"blackman") == 0)       {blackman (window, n);iswindow=true;}
  if (strcmp(formula,"blackmanharris") == 0) {blackmanharris (window, n);iswindow=true;}
  if (strcmp(formula,"bohman") == 0)         {bohman (window, n);iswindow=true;}
  if (strcmp(formula,"flattop") == 0)        {flattop (window, n);iswindow=true;}
  if (strcmp(formula,"gauss") == 0)          {gauss (window, n, fParameter);iswindow=true;}
  if (strcmp(formula,"hamming") == 0)        {hamming (window, n);iswindow=true;}
  if (strcmp(formula,"hann") == 0)           {hann (window, n);iswindow=true;}
  if (strcmp(formula,"nuttall") == 0)        {nuttall (window, n);iswindow=true;}
  if (strcmp(formula,"rectangular") == 0)    {rectangular (window, n);iswindow=true;}
  if (strcmp(formula,"triangular") == 0)     {triangular (window, n);iswindow=true;}
  if (strcmp(formula,"tuckey") == 0)         {tuckey (window, n, fParameter);iswindow=true;}
  if (strcmp(formula,"welch") == 0)          {welch (window, n);iswindow=true;}

  // user window
  if (iswindow==false) {cout << "CWB::Window::Window - window not defined" << endl;exit(1);}

  // normalize window
  Normalize(window, n);
}

CWB::Window::~Window() {
  delete [] window;
}

double
CWB::Window::GetValue(unsigned i) {
  if ((i<0) || (i>size)) {cout << "CWB::Window::GetValue - index not allowed" << endl;exit(1);}
  return window[i];
}

void
CWB::Window::Normalize (double* out_window, unsigned n) {
	double norm = 0;
	for (unsigned int i=0;i<n;i++) norm += pow(out_window[i],2);
	norm /= n;
	for (unsigned int i=0;i<n;i++) out_window[i] /= sqrt(norm);
}

void
CWB::Window::barthann (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++) {
	  double f = 0.0;
	  f = ((double) i) / ((double) (n - 1));
	  out_window[i] = 0.62 -0.48*(f - 0.5) +
	  0.38 * cos(2 * M_PI * (f - 0.5));
  }
}

void
CWB::Window::bartlett (double* out_window, unsigned n) {
  unsigned i = 0;
  unsigned odd = 0;

  odd = n % 2;

  for (i = 0; i < (n/2) +odd; i++)
    out_window[i] = 2.0 * ((double) i) / ((double) (n - 1));
  for (i = (n/2) + odd; i < n; i++)
    out_window[i] = 2.0 - 2.0 * ((double) i) / ((double) (n - 1));
}

void
CWB::Window::blackman (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    double f = 0.0;
    f = ((double) i) / ((double) (n - 1));
    out_window[i] = (0.42
		  - 0.5 * cos (2.0 * M_PI * f)
		  + 0.08 * cos (4.0 * M_PI * f));
  }
}

void
CWB::Window::blackmanharris (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    double f = 0.0;
    f = ((double) i) / ((double) (n - 1));
    out_window[i] = 0.35875 -
      0.48829 * cos(2.0 * M_PI * f) +
      0.14128 * cos(4.0 * M_PI * f) -
      0.01168 * cos(6.0 * M_PI * f);
  }
}

void
CWB::Window::bohman (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    double f = 0.0;
    f = (((double) i) -  ((double) (n / 2))) /
        ((double) (n / 2));
    out_window[i] = (1.0 - f) * cos(M_PI * f) +
      (1 / M_PI) * sin(M_PI * f);
  }
}

/* [UNIMPLEMENTED] chebyshev */

void
CWB::Window::flattop (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    double f = 0.0;
    f = ((double) i) / ((double) (n - 1));
    out_window[i] = 1.0 -
      1.93 * cos(2.0 * M_PI * f) +
      1.29 * cos(4.0 * M_PI * f) -
      0.388 * cos(6.0 * M_PI * f) +
      0.322 * cos(8.0 * M_PI * f);
  }
}

void
CWB::Window::gauss (double* out_window, unsigned n, double alpha) {
  unsigned i = 0;

  if (alpha < 2) alpha = 2.5;

  for (i = 0; i < n; i++)
  {
    double f = 0.0;
    f = (((double) i) -  ((double) (n / 2))) /
        ((double) (n / 2));
    out_window[i] = exp(-0.5 * (alpha * f) * (alpha * f));
  }
}

void
CWB::Window::hamming (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    double f = 0.0;
    f = ((double) i) / ((double) (n - 1));
    out_window[i] = 0.54 - 0.46 * cos (2.0 * M_PI * f);
  }
}

void
CWB::Window::hann (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    double f = 0.0;
    f = ((double) i) / ((double) (n - 1));
    out_window[i] = 0.5 - 0.5 * cos (2.0 * M_PI * f);
  }
}

/* [UNIMPLEMENTED] kaiser */

void
CWB::Window::nuttall (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    double f = 0.0;
    f = ((double) i) / ((double) (n - 1));
    out_window[i] = 0.3635819 -
      0.4891775 * cos(2.0 * M_PI * f) +
      0.1365995 * cos(4.0 * M_PI * f) -
      0.0106411 * cos(6.0 * M_PI * f);
  }
}

/* [UNIMPLEMENTED] parzen */
void
CWB::Window::rectangular (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    out_window[i] = 1.0;
  }
}

/* [UNIMPLEMENTED] saramaki */

/* [UNIMPLEMENTED] transversal */

void
CWB::Window::triangular (double* out_window, unsigned n) {
  unsigned i = 0;
  unsigned odd = 0;
  unsigned mirror = 0;

  /* helpful constants */
  odd = n % 2;
  mirror = (n + odd) / 2;

  /* fill the first half */
  for (i = 0; i < mirror; i++)
  {
    unsigned k = 0;
    k = 2 * (i + 1) + odd - 1;
    out_window[i] = (((double) k) / ((double) (n + odd)));
  }
  /* and mirror the other */
  for (i = mirror; i < n; i++)
  {
    out_window[i] = out_window[n - i - 1];
  }
}

void
CWB::Window::tuckey (double* out_window, unsigned n, double r) {
  unsigned i = 0;

  for (i = 0; i < (((double) n) / 2.0) * (1 + r); i++)
    out_window[i] = 1.0;
  for (; i < n; i++)
  {
    double f = 0.0;
    f = ((double) i) - (((double) n) / 2.0) * (1 + r);
    f = f / (((double) n) * (1 - r));
    out_window[i] = 0.5 * ( 1.0 + cos(M_PI * f));
  }
}

void
CWB::Window::welch (double* out_window, unsigned n) {
  unsigned i = 0;

  for (i = 0; i < n; i++)
  {
    out_window[i] = 1-pow(((double)i-(double)n/2.0)/((double)n/2.0),2);
  }
}
