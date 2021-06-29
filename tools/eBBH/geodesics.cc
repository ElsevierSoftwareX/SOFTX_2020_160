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


#include "geodesics.hh"

#include "math.h"
#include "stdlib.h"

//return values of f_util:

static double r, pr, phi, tau, r2, a2, r3, Delta, R02, Q, A, B, r_dot, w, pr_dot, Edot, Ldot;

// Here we compute second and third time derivatives of r and phi
// for a geodesic for use in the quadrupole formula
// REMOVED first three arguments! (made global)


void geo_higher_deriv(double E, double L, double a, 
   double& r_dotdot, double& r_dotdotdot, double& wdot, double& wdotdot)
{  
   // r2=r*r; a2=a*a; r3=r2*r;
   // Delta=r2+a2-2*r; R02=r2+a2*(1+2/r); Q=Delta/(E*R02-2*a*L/r) // Q=tau_dot
   // w=(L*Q+2*a/r)/R02 # w=phi_dot ; A=(w**2*(r3-a2)+2*a*w-1); B=(a2-r)/r3
   // r_dot=Delta*Q*pr/r2
   // pr_dot=A/r2/Q + B*pr**2*Q
   
   double F = r3+a2*r-2*r2;
   double Fdot = (3*r2+a2-4*r)*r_dot;
   double G = E*(r3+a2*r+2*a2)-2*a*L;
   double Gdot = E*(3*r2+a2)*r_dot;
   double Qdot = (Fdot*G-F*Gdot)/G/G;
   double H = L*Q*r+2*a;
   double Hdot = L*(Q*r_dot+Qdot*r);
   double I = (r3+a2*r+2*a2);
   double Idot = (3*r2+a2)*r_dot;
   wdot = (Hdot*I-H*Idot)/I/I;
   double  Adot = 2*w*wdot*(r3-a2)+w*w*(3*r2*r_dot)+2*a*wdot;
   double Bdot = (-3*a2/r2/r2+2/r3)*r_dot;
   double C = r2+a2-2*r;
   C *= C;
   double Cdot = 4*(r-1)*(r2+a2-2*r)*r_dot;
   double D = E*R02*r2-2*a*L*r;
   double Ddot = E*(4*r3+2*a2*r)*r_dot+r_dot*(2*a2*E-2*a*L);
   double alpha = Delta*Q/r2;
   double alphadot = (Cdot*D-C*Ddot)/D/D;
   r_dotdot =  alphadot*pr + alpha*pr_dot;
   double Cdotdot = 4*(r-1)*(r2+a2-2*r)*r_dotdot +  4*(r_dot)*(r2+a2-2*r)*r_dot 
      + 8*(r-1)*(r-1)*r_dot*r_dot;
   double Ddotdot = E*(4*r3+2*a2*r)*r_dotdot+ E*(12*r2+2*a2)*r_dot*r_dot +
   r_dotdot*(2*a2*E-2*a*L);
   double alphadotdot = (Cdotdot*D-C*Ddotdot)/D/D - 2*(Cdot*D-C*Ddot)*Ddot/D/D/D;
   double pr_dotdot = Adot/r2/Q-A*(2*r*r_dot*Q+r2*Qdot)/(r2*Q)/(r2*Q)+Bdot*pr*pr*Q
      + 2*B*pr*pr_dot*Q+B*pr*pr*Qdot;
   r_dotdotdot = alphadotdot*pr + 2*alphadot*pr_dot + alpha*pr_dotdot;
   double Fdotdot = (3*r2+a2-4*r)*r_dotdot+(6*r-4)*r_dot*r_dot;
   double Gdotdot = E*(3*r2+a2)*r_dotdot + E*6*r*r_dot*r_dot;
   double Qdotdot = (Fdotdot*G-F*Gdotdot)/G/G-2*(Fdot*G-F*Gdot)*Gdot/G/G/G;
   double Hdotdot = L*(Q*r_dotdot+Qdotdot*r+2*Qdot*r_dot);
   double Idotdot = (3*r2+a2)*r_dotdot + 6*r*r_dot*r_dot;
   wdotdot = (Hdotdot*I-H*Idotdot)/I/I-2*(Hdot*I-H*Idot)*Idot/I/I/I;
}


void reduced_quad_dot(double r, double r_t, double r_tt, double r_ttt, double phi,
double phi_t, double phi_tt, double phi_ttt, double (&I_tt)[3][3], double (&I_ttt)[3][3])
{  double X = r*r;
   double X_t = 2*r*r_t;
   double X_tt = 2*r_t*r_t+2*r*r_tt;
   double X_ttt = 6*r_t*r_tt+2*r*r_ttt;

   double cos2phi = cos(2*phi);
   double sin2phi = sin(2*phi);
   
   double A = cos2phi;
   double A_t = -2*sin2phi*phi_t;
   double A_tt = -4*cos2phi*phi_t*phi_t-2*sin2phi*phi_tt;
   double A_ttt = 8*sin2phi*pow(phi_t,3)-8*cos2phi*phi_t*phi_tt
         -4*cos2phi*phi_t*phi_tt-2*sin2phi*phi_ttt;

   double B = sin2phi;
   double B_t = 2*cos2phi*phi_t;
   double B_tt = -4*sin2phi*phi_t*phi_t+2*cos2phi*phi_tt;
   double B_ttt = -8*cos2phi*pow(phi_t,3) - 8*sin2phi*phi_t*phi_tt
         -4*sin2phi*phi_t*phi_tt + 2*cos2phi*phi_ttt;

   double M[3][3] = {1+3*A,3*B,0,   3*B,1-3*A,0,  0,0,-2};
   double M_t[3][3] = {3*A_t, 3*B_t, 0,   3*B_t, -3*A_t, 0,   0,0,0};
   double M_tt[3][3] = {3*A_tt, 3*B_tt, 0,   3*B_tt, -3*A_tt, 0,   0,0,0};
   double M_ttt[3][3] ={3*A_ttt, 3*B_ttt, 0,   3*B_ttt, -3*A_ttt, 0,   0,0,0};

   for(int i=0; i<3; ++i)for(int j=0; j<3; ++j){
      I_tt[i][j] = (X_tt*M[i][j]+X*M_tt[i][j]+2*X_t*M_t[i][j])/6;
      I_ttt[i][j] = (X_ttt*M[i][j]+X*M_ttt[i][j]+3*X_tt*M_t[i][j]+3*X_t*M_tt[i][j])/6;
   }
}
  
void quad_loss(double r, double r_t, double r_tt, double r_ttt, double phi,
    double phi_t, double phi_tt, double phi_ttt, double& E_dot, double& Lz_dot)

/* Compute rate of energy and angular momentum loss according to quadrupole formula
   Specifically, for a two body system (in units of total mass, M:=m1+m2=1) 
   where object one is located at d1 = d*m2 and d2 =d*m1 where
   d = (r*cos(phi),r*sin(phi)) this returns (dE/dt)/mu**2
   and (dL/dt)/mu**2 where mu is the reduced mass
   The arguments are r and its first three time derivatives
   and phi and its first three time derivatives. */
   
{  double I_tt[3][3], I_ttt[3][3];
   reduced_quad_dot(r,r_t,r_tt,r_ttt,phi,phi_t,phi_tt,phi_ttt, I_tt, I_ttt);

   double sum = 0;
   for(int i=0; i<3; ++i)for(int j=0; j<3; ++j) sum +=I_ttt[i][j]*I_ttt[i][j];

   E_dot = -0.2*sum;

   // Assuming angular momentum only changes in z direction
   sum = 0;
   for(int i=0; i<3; ++i)sum +=  I_tt[0][i]*I_ttt[i][1]- I_tt[1][i]*I_ttt[i][0];
   Lz_dot = -0.4*sum;
}
  
   
void f_util(double* y, double Jspin, double mu)
{  r = y[0];
   pr = y[1];
   phi = y[2];
   tau = y[3];
   double& E = y[4];
   double& L = y[5];
   double a = Jspin +mu*L;
   
   r2 = r*r;
   a2 = a*a;
   r3 = r2*r;
   
   Delta = r2 + a2 - 2*r;
   R02 = r2 + a2*(1+2/r);
   Q=Delta/(E*R02-2*a*L/r);  // Q=tau_dot
   
   w=(L*Q+2*a/r)/R02;   // w=phi_dot
   A=(w*w*(r3-a2)+2*a*w-1);
   B=(a2-r)/r3;
   
   r_dot=Delta*Q*pr/r2;
   pr_dot=A/r2/Q + B*pr*pr*Q;
   
   double r_dotdot, r_dotdotdot, wdot, wdotdot;
   geo_higher_deriv(E, L, a, r_dotdot, r_dotdotdot, wdot, wdotdot);

   quad_loss(r, r_dot, r_dotdot, r_dotdotdot, phi, w, wdot, wdotdot, Edot, Ldot);
   Edot = mu*Edot;
   Ldot = mu*Ldot;
}


// Returns h plus and h cross (times r) of GW from quadrupole formula.
// See note above quad_loss
void h_quad(double r, double r_t, double  r_tt, double  r_ttt,
   double  phi, double  phi_t, double phi_tt, double phi_ttt,
   double& hplus,  double& hcross)
{
   double I_tt[3][3];
   double I_ttt[3][3];
   reduced_quad_dot(r, r_t, r_tt, r_ttt, phi, phi_t, phi_tt, phi_ttt, I_tt, I_ttt);

   hplus = 2*I_tt[0][0];
   hcross = 2*I_tt[0][1];
}


// ud has Jspin, mu
int f(realtype t, N_Vector y, N_Vector ydot, void* ud)
{  realtype *ydata, *ydotdata;
   realtype* udr = (realtype*)ud;
   ydata = NV_DATA_S(y);
   ydotdata = NV_DATA_S(ydot);
   
   f_util(ydata, udr[0], udr[1]);
   ydotdata[0] = r_dot;
   ydotdata[1] = pr_dot;
   ydotdata[2] = w;
   ydotdata[3] = Q;
   ydotdata[4] = Edot;
   ydotdata[5] = Ldot;
   return 0; 
}

geodesic::geodesic(double r0, double phi0, double E0, double L0, int dir, double Jspin,
double mu, double dt)
{  if(dir!=1 && dir!=-1){
      printf("Error... dir must be +-1\n"); 
      exit(1);
   }
   y = N_VNew_Serial(6);
   ydata = NV_DATA_S(y);
   ydata[0] = r0;
   ydata[1] = 0.0;
   ydata[2] = phi0;
   ydata[3] = 0.0;
   ydata[4] = E0;
   ydata[5] = L0;
   udata[0] = Jspin;
   udata[1] = mu;
   tret = 0;
   f_util(ydata, Jspin, mu);
   double Q2 = Q*Q;
   ydata[1] = dir*sqrt(fabs(r2/Delta/Q2*(1-2/::r-Q2-w*(R02*w-4*sqrt(a2)/::r))));
 
#ifdef SUNDIALS_PACKAGE_VERSION 
   cvode_mem =CVodeCreate(CV_ADAMS, CV_FUNCTIONAL); 			  // SUNDIALS_VERSION = 2.7.0
#else
#if SUNDIALS_VERSION_MAJOR >   4 || (SUNDIALS_VERSION_MAJOR ==  4 && \
   (SUNDIALS_VERSION_MINOR >   0 || (SUNDIALS_VERSION_MINOR ==  0 && \
    SUNDIALS_VERSION_PATCH >=  1                                  )))     // SUNDIALS_VERSION >= 4.0.1
   cvode_mem =CVodeCreate(CV_ADAMS);
#else
   cvode_mem =CVodeCreate(CV_ADAMS, CV_FUNCTIONAL); 
#endif
#endif
   if(!cvode_mem){
      printf("cvode_mem initialization failed\n");
      exit(1);
   }
   CVodeSetUserData(cvode_mem, udata);
   
   if(CVodeInit(cvode_mem, f, tret, y) != CV_SUCCESS){
      printf("cvodeInit failed\n");
      exit(1);
   }
        
   if(CVodeSStolerances(cvode_mem, 1e-4, 1e-7)!= CV_SUCCESS){
      printf("failed setting tolerances\n");
      exit(1);
   } 
   
   r = pr = phi = tau = t = Et = Lt = hplus = hcross =0;
}

geodesic::~geodesic()
{  delete [] r; 
   delete [] pr; 
   delete [] phi;
   delete [] tau;
   delete [] t;
   delete [] Et;
   delete [] Lt;
   delete [] hplus;
   delete [] hcross;
}


double min(double x, double y)
{  return x<y ? x : y;}

bool geodesic::integrate(double D_t, int& N, double rmin, bool stop_alr)
{  double dt = D_t/N;
   this->r = new double[N+1];
   this->pr = new double[N+1];
   this->phi = new double[N+1];
   this->tau = new double[N+1];
   this->Et = new double[N+1];
   this->Lt = new double[N+1];
   this->t = new double[N+1];
   
   this->r[0] = ydata[0];
   this->pr[0] = ydata[1];
   this->phi[0] = ydata[2];
   this->tau[0] = ydata[3];
   this->Et[0] = ydata[4];
   this->Lt[0] = ydata[5];
   this->t[0] = tret;
   
   double r_stop = 0;
   if(stop_alr) r_stop = 2*( 1+cos( 2*(acos(-min(udata[0] + udata[1]*this->Lt[0], 1))/3 )));
   if(rmin>r_stop)r_stop = rmin;
   int i=0;
   while(i<N && r[i]>r_stop){
      if(CVode(cvode_mem, tret+dt, y, &tret, CV_NORMAL)!=CV_SUCCESS)return false;
      ++i;
      this->r[i] = ydata[0];
      this->pr[i] = ydata[1];
      this->phi[i] = ydata[2];
      this->tau[i] = ydata[3];
      this->Et[i] = ydata[4];
      this->Lt[i] = ydata[5];
      this->t[i] = tret;
      if(stop_alr){
         r_stop = 2*( 1+cos( 2*(acos(-min(udata[0] + udata[1]*this->Lt[i], 1))/3 )));
         if(rmin>r_stop)r_stop = rmin;
      }
   }
   bool ret = (i==N);
   
   N=i+1;
   hplus = new double[N];
   hcross = new double[N];
   double r_tt, r_ttt, phi_tt, phi_ttt;
   double hp, hc;
   
   for(i = 0; i<N; ++i){
      ydata[0] = this->r[i];
      ydata[1] = this->pr[i];
      ydata[2] = this->phi[i];
      ydata[3] = this->tau[i];
      ydata[4] = this->Et[i];
      ydata[5] = this->Lt[i];
      if(i==2000){
         // ydata[0] = 7.256704e+00;
         // ydata[1] = -4.276535e-02;
         // ydata[2] = 5.112497e+01;
         // ydata[3] = 1.728791e+03;
         // ydata[4] = 9.438350e-01;
         // ydata[5] = 3.180806e+00;     
      }
      f_util(ydata, udata[0], udata[1]);
      geo_higher_deriv(ydata[4], ydata[5], sqrt(a2), r_tt, r_ttt, phi_tt, phi_ttt);
      
      h_quad(::r,r_dot, r_tt, r_ttt, ::phi, w, phi_tt, phi_ttt, hp, hc);
      hplus[i] = udata[1]*hp;
      hcross[i] = udata[1]*hc;
   }
   return ret;
}


void printfutil(double* y)
{  printf("%.6le %.6le %.6le %.6le %.6le %.6le\n", r, pr, phi, tau, r2, a2);
   printf("%.6le %.6le %.6le %.6le %.6le %.6le\n", r3, Delta, R02, Q, A, B);
   printf("%.6le %.6le %.6le %.6le %.6le\n",r_dot, w, pr_dot, Edot, Ldot);
   double r_tt, r_ttt, phi_tt, phi_ttt;
   double hp, hc;
   
   geo_higher_deriv(y[0], y[1], sqrt(a2), r_tt, r_ttt, phi_tt, phi_ttt);
   printf("%.6le %.6le %.6le %.6le %.6le\n", sqrt(a2), r_tt, r_ttt, phi_tt, phi_ttt);
   h_quad(::r,r_dot, r_tt, r_ttt, ::phi, w, phi_tt, phi_ttt, hp, hc);
   printf("%.6le %.6le\n", hp, hc);
}
