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


struct c_complex
  { double re, im;
  };

struct _complex
  { double re, im;
    _complex(double r, double i = 0.0) { re = r; im = i; }
    _complex() { }					/* uninitialized complex */
    _complex(c_complex z) { re = z.re; im = z.im; }	/* init from denotation */
  };

// ADA 0.61.32
//extern _complex csqrt(_complex), cexp(_complex), expj(double);	    /* from complex.C */
//extern _complex evaluate(_complex[], int, _complex[], int, _complex);   /* from complex.C */
extern _complex evaluate(_complex*, int, _complex*, int, _complex);   /* from complex.C */
extern _complex csqrt(_complex);   /* from complex.C */
extern _complex cexp(_complex);    /* from complex.C */
extern _complex expj(double);	    /* from complex.C */

inline double hypot(_complex z) { return ::hypot(z.im, z.re); }
inline double atan2(_complex z) { return ::atan2(z.im, z.re); }

inline _complex cconj(_complex z)
  { z.im = -z.im;
    return z;
  }

inline _complex operator * (double a, _complex z)
  { z.re *= a; z.im *= a;
    return z;
  }

inline _complex operator / (_complex z, double a)
  { z.re /= a; z.im /= a;
    return z;
  }

inline void operator /= (_complex &z, double a)
  { z = z / a;
  }

extern _complex operator * (_complex, _complex);
extern _complex operator / (_complex, _complex);

inline _complex operator + (_complex z1, _complex z2)
  { z1.re += z2.re;
    z1.im += z2.im;
    return z1;
  }

inline _complex operator - (_complex z1, _complex z2)
  { z1.re -= z2.re;
    z1.im -= z2.im;
    return z1;
  }

inline _complex operator - (_complex z)
  { return 0.0 - z;
  }

inline bool operator == (_complex z1, _complex z2)
  { return (z1.re == z2.re) && (z1.im == z2.im);
  }

inline _complex sqr(_complex z)
  { return z*z;
  }

