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


/* mkfilter -- given n, compute recurrence relation
   to implement Butterworth, Bessel or Chebyshev filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 */

/* Routines for complex arithmetic */

#include <math.h>
#include "mkfcomplex.hh"

static _complex eval(_complex[], int, _complex);
static double Xsqrt(double);


_complex evaluate(_complex topco[], int nz, _complex botco[], int np, _complex z)
  { /* evaluate response, substituting for z */
    return eval(topco, nz, z) / eval(botco, np, z);
  }

_complex eval(_complex coeffs[], int npz, _complex z)
  { /* evaluate polynomial in z, substituting for z */
    _complex sum = _complex(0.0);
    for (int i = npz; i >= 0; i--) sum = (sum * z) + coeffs[i];
    return sum;
  }

_complex csqrt(_complex x)
  { double r = hypot(x);
    _complex z = _complex(Xsqrt(0.5 * (r + x.re)),
			Xsqrt(0.5 * (r - x.re)));
    if (x.im < 0.0) z.im = -z.im;
    return z;
  }

static double Xsqrt(double x)
  { /* because of deficiencies in hypot on Sparc, it's possible for arg of Xsqrt to be small and -ve,
       which logically it can't be (since r >= |x.re|).	 Take it as 0. */
    return (x >= 0.0) ? sqrt(x) : 0.0;
  }

_complex cexp(_complex z)
  { return exp(z.re) * expj(z.im);
  }

_complex expj(double theta)
  { return _complex(cos(theta), sin(theta));
  }

_complex operator * (_complex z1, _complex z2)
  { return _complex(z1.re*z2.re - z1.im*z2.im,
		   z1.re*z2.im + z1.im*z2.re);
  }

_complex operator / (_complex z1, _complex z2)
  { double mag = (z2.re * z2.re) + (z2.im * z2.im);
    return _complex (((z1.re * z2.re) + (z1.im * z2.im)) / mag,
		    ((z1.im * z2.re) - (z1.re * z2.im)) / mag);
  }

