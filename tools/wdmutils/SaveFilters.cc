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


#include "StdFunction.hh"
#include "NuFunction.hh"
#include "CompositeFunction.hh"
#include "CosFourier.hh"

#include "math.h"
#include "stdio.h"


// number of coefficents needed in the Fourier Expansion of
// cos, cos2, sin(2x)
const int nFourierCoeff[8][3]= { 0, 0, 0, 
 	0, 0, 0, 	
	400, 100, 600,	// n = 2
	500, 60, 800, 
	150, 45, 150, 	// n = 4
	125, 35, 150,
	80, 35, 80,
	60, 30, 80
};


// number of subintervals used in the numerical integration 
// required for the cosine expansion  of
// cos, cos2, sin(2x)
const int nFourierIntegral[8][3] = {0,0,0,   
	0,0,0,   		  // n=1
	2000, 500, 4200, // n=2 
	2700, 200, 6000, // n=3
	400,  150, 800,  // n=4
	300,  100, 700,  // n=5
	200,  100, 300,  // n=6
	200,  100, 300   // n=7
};


// the functions that nneed to be expanded :
// (besides cosine, which is predefined)

double cos2(double x)
{	double c = cos(x);
	return c*c;
}

double doublesin(double x)
{	return sin(2*x);
}


int main()
{  double (*const Func[3])(double) = { cos, cos2, doublesin };
   const char* varName[3] = {"Cos", "Cos2", "SinCos"};
   int n;
   for(int i = 0; i<3; ++i){
      StdFunction stdF(Func[i]);
      for(int nu=2; nu<8; ++nu){
         printf("   %s[%d] = new double[%d];\n", varName[i], nu, n = nFourierCoeff[nu][i]+1);
         printf("   %sSize[%d] = %d;\n", varName[i], nu, n);
         if(i)printf("   for(int i=0; i<%d; ++i)%s[%d][i] = 0;\n", n, varName[i], nu);
         NuFunction nuF(nu);
         CompositeFunction<double> smoothEdge(stdF, nuF);
         double* v = CosFourier(smoothEdge, n, nFourierIntegral[nu][i]);
         if(i==1)v[0] = 0.5;
         for(int j = 0; j<n; ++j){
            if(i==1)if(j && j%2==0)continue;
            if(i==2)if(j&1)continue;
            printf("   %s[%d][%d] = %+.17le;\n", varName[i], nu, j, v[j]);
         }
         printf("\n");
      }
	}
}
