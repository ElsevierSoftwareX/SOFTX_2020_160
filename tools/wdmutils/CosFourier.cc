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


#include "CosFourier.hh"

#include "Function.hh"
#include "StdFunction.hh"

#include "math.h"
#include "stdio.h"

static const double Pi = 3.14159265358979312;

double IntSimpson(double* v, int n)
{	double sum1 = 0, sum2 = 0;
	if(n%2 || n<2){
		printf("IntSimpson error: n must be even and non-null. return 0.\n"); 
		return 0;}
	
	for(int i=1; i<n-1; i+=2){
		sum1 += v[i]; 
		sum2 += v[i+1];
	}
	sum1 = v[0] + v[n] + 2*(2*(sum1+v[n-1]) + sum2) ;
	return sum1/(3*n); 
}

double* CosFourier(	const Function<double>& f, int maxN, int nSubInt, double (*integrator)(double*, int))
{	double* res = new double[maxN+1];   		   // expansion coeffs (result)
	double* fi = new double[nSubInt+1]; 		   // intermediate sums on f 
	double* Cos = new double[nSubInt+1];
	for(int i=0; i<=nSubInt; ++i)Cos[i] = cos(i*Pi/nSubInt);
	
	bool* done = new bool[maxN+1];
	for(int i=0; i<=maxN; ++i)done[i]=false;		// done[i] = true : coefficient i has been calculatedd
	
   while(1){
		int i;	                                 // a_i index
		for(i=maxN; i>0; --i)if(!done[i])break;	// found new coefficient to compute
		if(i==0)break;
	
		int nSI2 = nSubInt*i;							// total number of subintervals for [0,1]
		double aux = nSI2;
		 
		double* func = new double[nSI2+1];			// sampling the function
		for(int j=0; j<=nSI2; ++j) func[j] = f(j/aux);
		
		for(int j=i; j>0; --j)if(!done[j])if(i%j==0){	//divisor, a_j index to compute
			int mult = i/j;
			int HalfPeriodStep = nSubInt*mult;		// corresponding to j
			for(int k=0; k<=nSubInt; ++k){
				fi[k]=0;
				int i0 = k*mult;							// subsampling
				for(int l=0; l<j; ++l){ 				// loop over half-periods
					if(l&1)fi[k] -= func[i0];
					else fi[k] += func[i0];
					i0 += HalfPeriodStep;
				}
				fi[k] *= Cos[k];
			}
			res[j] = sqrt(2)*integrator(fi, nSubInt)/j; 
			// the range 1./j is not accounted for by the integrator
			
         done[j] = true;
		}
      
		if(i==maxN){
			for(int j=0, k=0; j<=nSI2; j+= maxN) fi[k++] = func[j];
			res[0] = integrator(fi, nSubInt);
			done[0] = true;
		}
		delete [] func;
	}
	
	delete [] fi;
	delete [] done;
	delete [] Cos;
	return res;
} 
