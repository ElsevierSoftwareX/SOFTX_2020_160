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


#include "NuFunction.hh"
#include "TMath.h"

#include "stdio.h"

NuFunction::NuFunction(int n)
{	this->n = n;
	if(n>0)return;
	printf("NuFunction::NuFunction : n value too small, reset to 1\n");
	n=1;
}

const Function<double>& NuFunction::Clone() const
{	return *new NuFunction(n);
}


double NuFunction::operator()(const double& x) const
{	static const double PiO2 = TMath::Pi()/2;
	if(n==1)return x*PiO2;
	double y=x*x;
   switch(n){ 
      case 2:  return (3-2*x)*y*PiO2;
      case 3:  return y*(10*x - (15 - 6*x)*y)*PiO2;      
      case 4:  return y*y*(35 - 84*x + (70 - 20*x)*y)*PiO2; 
      case 5:  return y*y*(126*x + y*( -420 + 540*x  + y*(-315 + 70*x) ))*PiO2; 
      default: return TMath::BetaDistI(x, n, n) *PiO2;
   }
}
	
