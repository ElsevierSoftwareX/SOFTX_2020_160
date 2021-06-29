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


#ifndef COSFOURIER_HH
#define COSFOURIER_HH

#include "Function.hh"


// Basic Simpson integration; the implicit integration interval is 1,  
// so the result must be multiplied by the actual integration interval, if different 
double IntSimpson(double* v, int n);


// Performs a cosine expansion with a specialized method thought to reduce numerical errors
// function f defined on the interval [0,1]
// maxN : maximum number of coefficients to be computed in the Cosine expansion expansion
// n : number of intervals (#points - 1) within a half-period used for numerical integration
// integ : integrator	
double* CosFourier(const Function<double>& f, int maxN, int n=2000, 
	double (*integ)(double*, int) = IntSimpson); 

#endif
