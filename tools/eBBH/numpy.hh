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


#ifndef NUMPY_HH
#define NUMPY_HH

#include "math.h"

double interp(double v, double* x, double* y, int n);
double interp(double v, double* x, double* y, int n, double left, double right);

double* interp(double* v, int Nv, double* x, double* y, int n, double left, double right);

class Complex{
public:
   Complex(double r=0., double i=0.){ a = r; b = i;}
   double Re() const{   return a;};
   double Im() const{   return b;};
private:
   double a, b;
};


Complex operator+(const Complex& x, const Complex& y);

Complex operator*(const Complex& x, const Complex& y);

Complex Exp(double phase);

double fabs(const Complex& x);

#endif
