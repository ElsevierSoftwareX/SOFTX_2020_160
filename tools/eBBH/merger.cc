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

static const double pi = 3.14159265358979312;


// Quasinormal ringdown of spin a BH (assuming M=1)
Complex ringdown(double a, double t)
{	double fQNR = (1.0-0.63*pow(1.0-a,0.3))/(2.0*pi);
	double Q = 2.0/pow(1.0-a,0.45);
	return exp(-pi*fQNR*t/Q)*Exp(2*pi*fQNR*t);
}

// IRS merger-ringdown, see http://arxiv.org/abs/0805.1428 and http://arxiv.org/abs/1107.1181
void irs_merger(double dt, double tmerge, double* tsig, Complex* hsig, 
int lentsig, double a, double tstart=-1)
{  if (tstart<0)tstart = tmerge;
	double hmerge = fabs(hsig[int(tmerge/dt)]);
   Complex& aux = hsig[int(tstart/dt)];
	double phiIRS = 0.5*atan2(aux.Im(), aux.Re());
	double omQNM = 1-0.63*pow(1-a, 0.3);
	double Q = 2.0/pow(1.0-a, 0.45);
	double b = 2.0*Q/omQNM;
	//#kappa = 0.426
	double kappa = 0.644;
	//#c = 0.252
	double c = 0.26;
	//#c = 2.0*(1.0-2.0/omQNM/omega0)/((1.0+1.0/kappa)**(1+kappa)-(1.0+1.0/kappa))
	double alpha = 72.3/(Q*Q);
	double fhat = c/2.*( pow(1 + 1./kappa, 1+kappa) - (1+1./kappa) );
   
	double OmIRS = omQNM/2.0*(1.0-fhat);
   double fdot = -c/b;
	double Amppeak = sqrt(fabs(fdot)/(1+alpha*(fhat*fhat- pow(fhat,4))))/(2.0*OmIRS);

	for(int i=tstart/dt; i<lentsig; ++i){
      double t = tsig[i]-tmerge;
	   fhat = c/2*pow(1+1./kappa,1+kappa)*( 1 - pow(1+exp(-2*t/b)/kappa, -kappa) );
      OmIRS = omQNM/2*(1-fhat);
      fdot = -c/b*pow(1+1./kappa,1+kappa)*pow(1+exp(-2*t/b)/kappa,-1.0-kappa)*exp(-2*t/b);
		double Amp = hmerge*sqrt(fabs(fdot)/(1+alpha*(fhat*fhat-pow(fhat, 4))))/(2*OmIRS*Amppeak);
      phiIRS += OmIRS*dt;
		hsig[i] = Amp*Exp(2*phiIRS);
   }
}
