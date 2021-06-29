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


#include "eBBH.hh"

#include "egw_model.hh"
#include "geodesics.hh"
#include "merger.hh"
#include "numpy.hh"


#include "stdio.h"
#include "stdlib.h"
#include "math.h"
    
inline double min(double x, double y)
{  return x<y? x: y;}

int Sample(Complex* H, double* t, int N, double totalMass, double Rate, wavearray<double>& hps, wavearray<double>& hxs)
{  const double C = 299792458;
   const double G = 6.67428e-11;
   const double SolarMass = 1.98892e30;
   const double L = G/C/C*SolarMass*totalMass;
   const double T = L/C; //printf("T = %lf\n", T);
	
	double* tt = new double[N]; 
   
	for(int i=0; i<N; ++i){
		tt[i] = t[i]*T;
		if(i)tt[i] -= tt[0];
	}
	tt[0] = 0;
	//printf("T*16384 = %lf\n", T*16384);
	double tMax = tt[N-1]; //printf("tMax = %lf\n", tMax);
   int nSamples = 0;
   if(tMax*Rate<2.1e9)nSamples = tMax*Rate;
   //printf("nSamples = %d\n", nSamples);
  
	hps.resize(nSamples); hps.rate(Rate);
	hxs.resize(nSamples); hxs.rate(Rate);
	
	int j=1;
	for(int i=0; i<nSamples; ++i){
      double ts = i/Rate;
		while(ts>tt[j])++j;
      if(j>=N)printf("ERROR\n");
		double frac = (ts - tt[j-1])/(tt[j] - tt[j-1]);
		hps[i] = H[j-1].Re() + (H[j].Re() - H[j-1].Re())*frac;
		hxs[i] = H[j-1].Im() + (H[j].Im() - H[j-1].Im())*frac;
	}
   delete [] tt;

	return nSamples;
}

int getEBBH(double m1, double m2, double rmin0, double e0, wavearray<double>& Hp, wavearray<double>& Hx, double t_end)
{  if(m1<0 || m2<0){
      printf("Error! masses must be >=0\n");
      return -1;
   }
   
   if(rmin0<0){
      printf("Error! rmin must be >=0\n");
      return -2;
   }
   if(e0>1){
      printf("Error! Orbit must be initially parabolic or bound.\n");
      return -3;
   }
   
   //double q = m1*m2/(m1+m2)/(m1+m2);
   double q = m1/m2; 
   if(q>1)q=1./q;
   double rmin = rmin0;
   double e = e0;
   double Mbh = 1/(1+q);
   double Mns = 1 - Mbh;
   //printf("Mbh = %.3lf\nMns = %.3lf\n", Mbh, Mns);
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
	double phi0 = 0; 
   //printf("%.6le %.6le %.6le %.6le %d %.6le %.6le\n", r0, phi0, E0, L0, dir, Jspin, mu);
	geodesic geo(r0, phi0, E0, L0, dir, Jspin, mu);
	
   bool success = geo.integrate(t_end, N_plot, 0, true);    // memory LEAK!
   
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

   int retVal = Sample(hsig, tsig, Ntsig, m1+m2, 16384., Hp, Hx)==0;
   delete [] hsig;
   delete [] hr;
   delete [] hi;
   delete [] tsig;
   return retVal;
}
