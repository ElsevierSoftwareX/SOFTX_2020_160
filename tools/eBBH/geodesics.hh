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


#ifndef GEODESICS_HH
#define GEODESICS_HH

#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"  /* serial N_Vector types, fcts., macros */

void printfutil(double *y);
void f_util(double* y, double Jspin, double mu);

class geodesic{
public:
   geodesic(double r0, double phi0, double E0, double L0, 
      int dir, double Jspin, double mu=0, double dt=1);
   ~geodesic();
   
   bool integrate(double D_t, int& N, double rmin=0, bool stop_alr=false);
   //returns "success" and the variables below:
   double* r, *pr, *phi, *tau, *t, *Et, *Lt, *hplus, *hcross;
      
private:
   N_Vector y;
   realtype udata[2], tret;      //udata[0] = Jspin, udata[1] = mu
   realtype* ydata;
   void* cvode_mem;

};

#endif
