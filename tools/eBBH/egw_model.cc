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

#include "stdio.h"

static const double pi = 3.14159265358979312;

double E_ecc(double e)
{  double e2 = e*e;
   return (1+73.0/24*e2+37.0/96*e2*e2)/pow(1.0+e,3.5);
}

double J_ecc(double e)
{  return (1.0+7.0/8*e*e)/(1.0+e)/(1+e);
}

double bisect(double (*f)(double), double a, double b)
{  double fa = f(a);
   double fb = f(b);
   if (fa == 0.) return a;     
   if (fb == 0.) return b;     
   
   double tol = 1.0e-14;
        
   while(fabs(b-a)>tol){
      double c = 0.5*(a+b);
      double fc = f(c);
      if (fc==0.)return c;
      if (fb*fc < 0.0){
         a = c;
         fa = fc;
      }
      else{
         b = c;
         fb = fc;
      }
   }
   return 0.5*(a+b);
}

double zw_func(double p, double e, double a)
{  //"Equals zero for zoom-whirl parameters"
   // Formulae from Glampedakis and Kennefick
   double Mbh = 1.0;
   double F = ( p*p*p - 2*Mbh*(3.0+e*e)*p*p + Mbh*Mbh*pow(3.0+e*e, 2)*p -
                4*Mbh*a*a*pow(1.0-e*e, 2) )/p/p/p;
   double N = (2.0/p)*(-Mbh*p*p+(Mbh*Mbh*(3.0+e*e)-a*a)*p-Mbh*a*a*(1+3*e*e) ); 
   double C = a*a-Mbh*p;
   C *= C;
   double xsq, des = fabs(N*N-4.0*F*C);
   if (a==0.)des = 0;
   if (a>0)xsq = (-N - sqrt(des))/2.0/F; //  Prograde motion
   else xsq = (-N + sqrt(des))/2.0/F;    //  Retrograde motion

   return p*p-xsq*(1.0+e)*(3.0-e);  
}

double rw_e, rw_a;

double rw_f(double p)
{  return zw_func(p, rw_e, rw_a);
}

double rzoom_whirl(double e, double a)
{  //"Value of rp for zoom-whirl orbit in Kerr spacetime"
   rw_e = e;
   rw_a = a;
   double plow, phigh;
   if (a>0){
      plow = 1.0+e;
      phigh =  6.0 + 2.0*e;
   }
   else if(a<0){
      plow = 6.0+2.0*e;
      phigh = 15;
   }
   else return (6.0+2.0*e)/(1.0+e);
   
   double p = bisect(rw_f, plow, phigh);
   return p/(1.0+e);
}

double scaling_exp_zw(double a, double r0)
{  // "Instability exponent as a function of BH spin and radius of unstable circular orbit"
   double r = r0;
   double m = 1;
   double w, R0 =sqrt(r*r + a*a*(1+2*m/r));
   if (a>0)w = m/(m*a + r*sqrt(m*r));
   else w = m/(m*a- r*sqrt(m*r));
   double Delta = r*r + a*a - 2*m*r;
   return  r*r/(2*pi)/sqrt( fabs( 3.0*r*r*Delta + 4.0*m/w/w*(r*R0*R0*w*w - 4*m*a*w - r + 2*m) ) );
}

double scaling_exp(double e, double a)
{  double gamma_parabolic = 0.19; //# Exponent for e=1 a=a_para
   double a_para  = 0.5; 
   return scaling_exp_zw(a, rzoom_whirl(e,a))/scaling_exp_zw(a_para, rzoom_whirl(1,a_para))*gamma_parabolic;
}


void close_encounter_pm(double e, double rmin, double mu, double& e_new, double& rmin_new)
{ //"Change in orbital parameters due to a close encounter according to Peters & Mathews"
	double a = mu/rmin; 
   double E = a*(e-1.0)/2.0;
   double const b = (64.0*3.14159265358979312)/5.0; 
	double L = mu*sqrt(rmin*(1.0+e));
	double deltaE = b*mu*mu*pow(1.0/rmin, 3.5)*E_ecc(e);
	double deltaL = b*a*a*J_ecc(e);
   
   E -= deltaE;
   L -= deltaL;
	e_new = sqrt(1.0+2*E*L*L/mu/mu/mu);
   rmin_new = (e_new-1.0)/2.0/(E/mu); 
}

void gwave_data(double e, double rmin, double mu, double& e_new, double& rmin_new)
{  //      "Change in orbital parameters due to g-wave results from simulations"
	double E = mu*(e-1.0)/2.0/rmin;
	double L = mu*sqrt(rmin*(1.0+e));
	double rp[7]  = {6.95, 7.22, 7.5, 8.75, 10, 12.5, 15};
   double E_data[7] = { 0.006753, 0.003573, 0.002420, 
                        0.000730, 0.000334, 0.000105, 3.081720e-05};
	double J_data[7] = { 0.071858, 0.044782, 0.035159, 
                        0.015815, 0.009692, 0.004601, 2.144661e-03};
         
	
	double deltaE = interp(rmin, rp, E_data, 7)*E_ecc(e)/E_ecc(1);
	double deltaL = interp(rmin, rp, J_data, 7)*J_ecc(e)/J_ecc(1);
   E -= deltaE;
   L -= deltaL;
   e_new = rmin_new = 0;
	if ((1.0+2*E*L*L/mu/mu/mu)<0)return;  
   e_new = sqrt(1.0+2*E*L*L/mu/mu/mu);
   rmin_new = (e_new-1.0)/2.0/(E/mu); 
}

double getDelta()
{  return 3;
}

double getE0(double e, double rmin, double mu, double rc)
{   return mu*(e-1.0)/2.0/rmin + mu/2.0/rc;
}

            
double newt_crit_radius(double e, double a)
{	//        "Effective periastron distance at large separation for the BH-NS system to merge function of eccentricity and BH spin"
   double rparabolic = 6.89;  // # Capture for e=1, a=0.5 orbit
   double a_para  = 0.5;
   double r = (rparabolic/rzoom_whirl(1.0,a_para))*rzoom_whirl(e,a);
   return r;
}         

double newt_getDeltaE(double e, double rmin, double mu, double a)
{  double rc = newt_crit_radius(e, a);
   double E0 = getE0(e,rmin,mu,rc);
   double Delta = getDelta();
   double gamma = scaling_exp(e,a);
   return  E0*( 1 - pow((rmin-rc)/Delta,gamma) );
}

double newt_getDeltaL(double e, double rmin, double mu, double a)
{  double rc = newt_crit_radius(e, a);
   double L0 = mu*(sqrt(rmin*(1.0+e)) - sqrt(rc));
	double Delta = getDelta();
   double gamma = scaling_exp(e,a);
   return L0*( 1 - pow((rmin-rc)/Delta, gamma) );
}
         
void newt_close_encounter(double e, double rmin, double mu, double a, double& e_new, double& rmin_new)
{  //"Change in orbital parameters due to a close encounter according to ZW model"
	double E = mu*(e-1.0)/2.0/rmin;
	double L = mu*sqrt(rmin*(1.0+e));
	double deltaE =  newt_getDeltaE(e, rmin, mu, a);
   double deltaL =  newt_getDeltaL(e, rmin, mu, a);
	double Delta = getDelta();
   double rc = newt_crit_radius(e, a);
   // Stitch on P&M for large rmin
   
   if ((rmin-rc)>=Delta or (64.0*pi)/5.0*mu*mu*pow(1.0/rmin,3.5)*E_ecc(e)>deltaE)
      return close_encounter_pm(e, rmin, mu, e_new, rmin_new);

   E -= deltaE;
	L -= deltaL;
   double f = 1.0+2*E*L*L/mu/mu/mu;
	if (f<0){
      printf("E=%e deltaE=%e rmin=%e e=%e\n", E, deltaE, rmin, e);
      f = 0;
   }
   e_new = sqrt(f);
   rmin_new = (e_new-1.0)/2.0/(E/mu);
}
           

double newt_amp_enhance(double e, double rmin, double mu, double a)
{  double gamma = scaling_exp(e,a);
	double Delta = getDelta();
	double rc = newt_crit_radius(e, a);
   double E0 = getE0(e,rmin,mu,rc);
   double deltaE = E0;
   if (rmin>rc)deltaE =  E0*(1- pow((rmin-rc)/Delta, gamma) );
                
   double deltaEpm = (64.0*pi)/5.0*mu*mu*pow(1.0/rmin, 3.5)*E_ecc(e);
   double ratio = deltaE/deltaEpm;
   if (fabs(rmin-rc)>Delta || ratio < 1)ratio = 1;
   return sqrt(ratio);
}


double energy_geo(double rp, double e, double a)
{
   double p = rp*(1+e);
   double Mbh = 1.0;
   double F = ( p*p*p - 2*Mbh*(3+e*e)*p*p + p*Mbh*Mbh*pow(3.0+e*e, 2)-4*Mbh*a*a*pow(1.0-e*e,2) )/p/p/p;
   double N = (2/p)*( -Mbh*p*p + p*( Mbh*Mbh*(3.0+e*e)-a*a) - Mbh*a*a*(1+3*e*e) );
   double C = a*a-Mbh*p; 
   C *=C;
   
   double xsq, des = fabs(N*N-4*F*C);
   if(a==0)des = 0;
   if (a>=0)xsq = (-N - sqrt(des))/2/F;   // Prograde motion
   else xsq = (-N + sqrt(des))/2/F;       // Retrograde motion
   return sqrt( fabs( 1 - (Mbh/p)*(1.0-e*e)*(1.0-xsq/p/p*(1.0-e*e))  ) );
}

double ang_mom_geo(double rp, double e, double a)
{  
   double p = rp*(1+e);
   double Mbh = 1.0;
   double F = ( p*p*p - 2*Mbh*(3+e*e)*p*p + p*Mbh*Mbh*pow(3.0+e*e, 2)-4*Mbh*a*a*pow(1.0-e*e,2) )/p/p/p;
   double N = (2/p)*( -Mbh*p*p + p*( Mbh*Mbh*(3.0+e*e)-a*a) - Mbh*a*a*(1+3*e*e) );
   double C = a*a-Mbh*p; 
   C *=C;
   
   double xsq, des = fabs(N*N-4*F*C);
   if(a==0)des = 0;
   if (a>=0)xsq = (-N - sqrt(des))/2/F;   // Prograde motion
   else xsq = (-N + sqrt(des))/2/F;       // Retrograde motion
   double E = sqrt( 1 - (Mbh/p)*(1.0-e*e)*(1.0-xsq/p/p*(1.0-e*e)) );
   return sqrt(fabs(xsq)) + a*E;
}

double ang_mom_eff_geo(double rp, double e, double Jspin, double mu)
{  double a = Jspin+mu*sqrt((1+e)*rp);
   double a1 = Jspin+mu*ang_mom_geo(rp, e, a);
   while(fabs(a-a1)>1.0e-5){
      a = a1;
      a1 = Jspin+mu*ang_mom_geo(rp,e,a);
   }         
   return ang_mom_geo(rp,e,a);
}

double opg_x;
double opg_E;
double opg_L;
double opg_a;

double opg_ecc(double p)
{  double x = opg_x;
   double E = opg_E;
   double des = pow(p/x,4) + 4*p*p*p/x/x*(E*E-1);
   if (des<0)des = 0;
   double y = 0.5*(pow(p/x, 2)-sqrt(des));
   return sqrt(fabs(1-y));                
}

double opg_f(double p)
{  double e = opg_ecc(p);
   return opg_E - energy_geo(p/(1+e),e,opg_a);
}

double opg_g(double p)
{  return ang_mom_geo(p/2.0,1.0, opg_a)- opg_L;
}
 
void orbital_param_geo(double E, double L, double a, double& rp, double& e)
{  opg_x = L-a*E;
   opg_E = E;
   opg_a = a;      
   opg_L = L;
   if (E>1)printf("Error in orbital_param_geo, E>0\n");
   double p;
   if (E==1.)p = bisect(opg_g,1.0,100);
   else p = bisect(opg_f, 1.0, 100);
   e = opg_ecc(p);
   rp = p/(1.0+e);
}

double geo_crit_radius(double e, double a)
{  return rzoom_whirl(e, a);
}
         
double geo_getDeltaE(double e, double rmin, double mu, double a)
{  double rc = geo_crit_radius(e, a);
   double E0 = getE0(e,rmin,mu,rc);
   double Delta = getDelta();
   double gamma = scaling_exp(e,a);
   return E0*(1.0 - pow ((rmin-rc)/Delta, gamma) );
}


double geo_getDeltaL(double e, double rmin, double mu, double a)
{  double rc = geo_crit_radius(e, a);
   double L0 = mu*(sqrt(rmin*(1.0+e)) - sqrt(rc));
   double Delta = getDelta();
   double gamma = scaling_exp(e,a);
   return L0*(1.0 - pow( (rmin-rc)/Delta, gamma ) );
}
                
void geo_close_encounter(double e, double rmin, double mu, double a, double& e_new, double& rmin_new)
{  double E = mu*energy_geo(rmin,e,a); 
	double L = mu*ang_mom_geo(rmin,e,a); 

	double deltaE =  geo_getDeltaE(e, rmin, mu, a);
   double deltaL =  geo_getDeltaL(e, rmin, mu, a);
	double Delta = getDelta();
   double rc = geo_crit_radius(e, a);
	
   if ((rmin-rc)>=Delta || (64.0*pi)/5.0*mu*mu*pow(1.0/rmin, 3.5)*E_ecc(e)>deltaE){
	   deltaE = (64.0*pi)/5.0*mu*mu*pow(1.0/rmin,3.5)*E_ecc(e);
		deltaL = (64.0*pi)/5.0*mu*mu/rmin/rmin*J_ecc(e);
   }
	E = E - deltaE;
   L = L - deltaL;
   orbital_param_geo(E/mu,L/mu,a,rmin_new, e_new);
	if (e_new==1)
			printf("%lf %lf %lf %lf\n", E/mu, L/mu,
            energy_geo(rmin_new,e_new,a), ang_mom_geo(rmin_new,e_new,a) ); 
}


double crit_radius(double e, double a)
{  // "Effective periastron distance at large separation for the BH-NS system to merge function of eccentricity and BH spin"
   double rparabolic = 6.89;  // Capture for e=1, a=0.5 orbit        
   double a_para  = 0.5; 
   return  ( rparabolic/rzoom_whirl(1, a_para) )*rzoom_whirl(e,a);
}


double geo_amp_enhance(double e, double rmin, double mu, double a)
{  return 1;}

double amp_enhance(double e, double rmin, double mu, double a)
{  double gamma = scaling_exp(e,a);
   double Delta = getDelta();
   double rc = crit_radius(e, a);
   double E0 = getE0(e,rmin,mu,rc); 
   double deltaE =  E0;
   if (rmin>rc)deltaE =  E0*(1.0-pow((rmin-rc)/Delta, gamma) );
     
   double deltaEpm = (64.0*pi)/5.0*mu*mu*pow(1.0/rmin,3.5)*E_ecc(e);     
   double ratio = deltaE/deltaEpm;
   if (fabs(rmin-rc)>Delta || ratio < 1)ratio = 1.0;
   return sqrt(ratio);
}


void orbital_param_newt_to_geo(double rmin,double e, double a, double& rp, double& ee)
{  double E = 1+(e-1.0)/2.0/rmin;
   double L = sqrt(rmin*(1.0+e));
   orbital_param_geo(E,L,a, rp, ee);
}
       
void orbital_param_geo_to_newt(double rmin, double e, double a, double& rmin_new, double& e_new)
{  double E = energy_geo(rmin,e,a) - 1.0;
   double L = ang_mom_geo(rmin,e,a);
   double f = 1.0+2*E*L*L;
   e_new = sqrt(f);
   if (E>0)rmin_new = (e_new-1.0)/2.0/E;
   else if (E==0.) rmin_new = L*L/2.0;
}


double a_eff(double e, double rp, double mu, double Jbh)
{        return 0.5*sqrt( (1+e)*rp/(2*6.89) )+Jbh; 
}       
