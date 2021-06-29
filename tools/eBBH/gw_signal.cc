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


#include "egw_model.hh"
#include "geodesics.hh"
#include "merger.hh"
#include "numpy.hh"

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
    
inline double min(double x, double y)
{  return x<y? x: y;}


int main(int argc, char** argv)
{  
   if(argc<3){
      printf("Usage: %s <rmin>  <eccentricity> [<mass ratio>] [PLOT] [<TMAX>]\n",
         argv[0]);
      return -1;
   }
   
   double rmin0 = atof(argv[1]);
   if(rmin0<0){
      printf("Error! rmin must be >=0\n");
      return -2;
   }
   double e0 = atof(argv[2]);
   if(e0>1){
      printf("Error! Orbit must be initially parabolic or bound.\n");
      return -3;
   }
   double q = 0.25;
   
   double t_end = 10000000.0;
   bool PLOT = false;
   if(argc>3){
      q = atof(argv[3]);
      if(q<0 || q>1){
         printf("Error! mass ratio must be between 0 and 1\n");
         return -4;
      }
      if(argc>4){
         PLOT = true;
         if(argc>5)t_end = atof(argv[5]);
      }
   }
   
   double rmin = rmin0;
   double e = e0;
   double Mbh = 1/(1+q);
   double Mns = 1 - Mbh;
   printf("Mbh = %.1lf\nMns = %.1lf\n", Mbh, Mns);
   double mu = Mns*Mbh;
   
   double abh = 0.0;
   double dt = 1;
   double a_final = 0;
   double t_post_merger = 400;
   
   double Jspin = abh*Mbh*Mbh;
	double L0 = ang_mom_eff_geo(rmin0, e0, Jspin, mu);
	double a0 = Jspin+mu*L0;
	double E0 = energy_geo(rmin0, e0, a0);
	double r0 = 1000;
	if (e0<1 && r0 > rmin0*(1+e0)/(1-e0)) r0 = rmin0*(1+e0)/(1-e0);
   
	int dir = -1;  //inward
	int N_plot = t_end/dt; //printf("%d\n", N_plot);
	double phi0 = 0; //printf("%.6le %.6le %.6le %.6le %d %.6le %.6le\n", r0, phi0, E0, L0, dir, Jspin, mu);
	geodesic geo(r0, phi0, E0, L0, dir, Jspin, mu);
	
   bool success = geo.integrate(t_end, N_plot, 0, true);
   printf("%d\n", N_plot);
      
   FILE* f = fopen("path_cc.dat", "w");
   for(int i=0; i<N_plot; ++i)
   fprintf(f, "%1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n", 
      geo.t[i], geo.r[i], geo.pr[i], geo.phi[i],  geo.tau[i], geo.Et[i], geo.Lt[i],
      geo.hplus[i], geo.hcross[i]);
   fclose(f); 
   
   double tmerge;
   if (success){
		tmerge = -1;
		t_post_merger=0;
   }
	else{
		double* r_LR = new double[N_plot];
      for(int i=0; i<N_plot; ++i)
         r_LR[i] = 2*(1+cos(2*(acos(-min(mu*geo.Lt[i] + Jspin, 1.))/3))) - geo.r[i];  //Light ring 
		tmerge = interp(0, r_LR, geo.t, N_plot, 1., 1.);
      delete [] r_LR;
   }
   
   int Ntsig = (geo.t[N_plot-1]+t_post_merger)/dt;
 	double* tsig = new double[Ntsig];
   for(int i=0; i<Ntsig; ++i) tsig[i] = dt*i;
   	
   double* hr = interp(tsig, Ntsig, geo.t, geo.hplus, N_plot, 0, 0);
   double* hi = interp(tsig, Ntsig, geo.t, geo.hcross, N_plot, 0, 0);
   
	Complex* hsig = new Complex[Ntsig];
   for(int i=0; i<Ntsig; ++i)hsig[i] = Complex(hr[i], hi[i]);
   
   a_final = mu*geo.Lt[N_plot-1] + abh*Mbh*Mbh;  // Spin of final BH
	
	double tstart_merger = tmerge; 
   if (tmerge>0) irs_merger(dt, tmerge, tsig, hsig, Ntsig, a_final, tstart_merger);

   f=fopen("wave_cc.dat", "w");
   for(int i=0; i<Ntsig; ++i)
      fprintf(f, "%1.8e %1.8e %1.8e\n", tsig[i], hsig[i].Re(), hsig[i].Im());
   fclose(f);
   //hsigr = real(hsig)
   //hsigi = imag(hsig)
  
   //N.savetxt('t.dat', tsig, fmt='%1.8e', delimiter='\n')
   //N.savetxt('hp.dat', hsigr, fmt='%1.8e', delimiter='\n')
   //N.savetxt('hc.dat', hsigi, fmt='%1.8e', delimiter='\n')
}
