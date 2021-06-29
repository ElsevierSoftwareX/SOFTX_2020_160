/*
# Copyright (C) 2019 Sergey Klimenko, Valentin Necula
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


#include "numpy.hh"

double interp(double v, double* x, double* y, int n)
{  int i;
   for(i=0; i<n; ++i)if(v<x[i])break;
   if(i==0)return y[0];
   if(i==n)return y[n-1];
   return y[i-1] + (v-x[i-1])/(x[i] - x[i-1]) *(y[i] - y[i-1]);
}

double interp(double v, double* x, double* y, int n, double left, double right)
{  int i;
   for(i=0; i<n; ++i)if(v<x[i])break;
   if(i==0)return left;
   if(i==n)return right;
   return y[i-1] + (v-x[i-1])/(x[i] - x[i-1]) *(y[i] - y[i-1]);
}

double* interp(double* v, int Nv, double* x, double* y, int n, double left, double right)
{  double* res = new double[Nv];
   for(int i=0; i<Nv; ++i)res[i] = interp(v[i], x, y, n, left, right);
   return res;
}

Complex operator+(const Complex& x, const Complex& y)
{  return Complex(x.Re() + y.Re(), x.Im() + y.Im());
}

Complex operator*(const Complex& x, const Complex& y)
{  return Complex(x.Re()*y.Re() - x.Im()*y.Im(), x.Re()*y.Im() + x.Im()*y.Re());
}

Complex Exp(double phase)
{  return Complex(cos(phase), sin(phase));
}

double fabs(const Complex& x)
{  return sqrt(x.Re()*x.Re() + x.Im()*x.Im());
}
