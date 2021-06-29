/*
# Copyright (C) 2019 Sergey Klimenko, Valentin Necula, Vaibhav Tiwari
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


#include "GNGen.hh"

#include "TF1.h"

#include "math.h"

static const double G  = 6.67384e-11;
static const double SM = 1.9891e30;
static const double PC = 3.08567758e16;
static const double C = 299792458;
static const double Pi = 3.14159265358979312;

static const double DistUnit = G*SM/C/C;
static const double TimeUnit = G*SM/C/C/C;

// work in units G=c=SM=1


GNGen::GNGen(double mSMBH, double mmin, double mmax, double beta): rnd(0)
{  
   this->beta = beta;
   smbhM = mSMBH;
   
   minM = mmin; 
   maxM = mmax; 
   freqCutoff = 20.;
   
   double normm = (beta-1)/(1./pow(minM, beta-1) - 1./pow(maxM, beta-1));
    
   int nBins = 100;
   double dm = (maxM - minM)/nBins;
         
   double rmin = 0.001;      // pc
   double rmax = 1;         // pc   
   double rrmin = rmin*PC/DistUnit;
   double rrmax = rmax*PC/DistUnit;
   
   double dr = (rmax - rmin)/nBins;
   double r = rmin + dr/2;
   
   double p0 = 0.5;
   double const bmaxc = pow(340*Pi/3, 1./7);
   double totrate = 0;
   dGammadmdMdr = new TH3F("d", "", nBins, mmin, mmax, nBins, mmin, mmax, nBins, rmin, rmax); 
   
   for(int i=0; i<nBins; ++i){
      double rr = r*PC/DistUnit;          // r in the chosen units
      double vofr = sqrt(mSMBH/rr);
      //printf("Sanity check: %lf vs %lf\n", vofr, sqrt(G*mSMBH*SM/r/PC)/C);
      
      for(int j=0; j<nBins-1; ++j){
         double m1 = mmin + (j+0.5)*dm;
         double tmp = 3./2 - p0*m1/maxM;
         double norm1 = tmp/( pow(rrmax, tmp) - pow(rrmin, tmp) )/4/Pi ; 
         
         for(int k=j+1; k<nBins; ++k){
            double m2 = mmin + (k+0.5)*dm;
            double M = m1 + m2;
            double eta = m1*m2/M/M;
                
            double bmax = bmaxc*M*pow(eta, 1./7)/pow(vofr, 9./7);    // eq. 17
            double sigmaC = Pi*bmax*bmax; 
            
            tmp = 3./2 - p0*m2/maxM;
            double norm2 = tmp/( pow(rrmax, tmp) - pow(rrmin, tmp) )/4/Pi ;
            double bhDensityProd = normm*normm*norm1*norm2/pow(rr, 3 + p0*M/maxM)/pow(m1*m2, beta);  // n1(r)*n2(r)
            
            double rate = bhDensityProd*sigmaC*vofr*4*Pi*rr*rr;
            
            totrate += rate;
            dGammadmdMdr->Fill(m1, m2, r, rate/1e-48);      //4Pi ignored
            //dGammadmdMdr->Fill(m1, m2, r, 1./pow(m1, beta)/pow(m2, beta));
         }
      }
      r+= dr;
   }  
   
   //printf("totrate (yr) = %le\n", totrate*dm*dm*dr*PC/DistUnit/TimeUnit*365.*24*3600);
}

GNGen::GNGen(const GNGen& x) : rnd(0)
{  minM = x.minM;
   maxM = x.minM;
   beta = x.beta;
   smbhM = x.smbhM;
   dGammadmdMdr = (TH3F*)x.dGammadmdMdr->Clone();
   freqCutoff = x.freqCutoff;
}

GNGen::~GNGen()
{  delete dGammadmdMdr;
}

void GNGen::setFreqCutoff(double f)
{  freqCutoff = f;
}

static double g(double e)   // used in a(e)
{  const double c1 = 12./19;
   const double c2 = 121./304;
   const double c3 = 870./2299;
   double e2 = e*e;
   return exp(c1*log(e)+c3*log(1+c2*e2))/(1-e2);
   //TF1 gG("gG4", "19./48*exp ( 4*0.631578947368421018*log(x) + 4*0.37842540234884730*log(1+0.398026315789473673*x*x) - 0.5*log(1-x*x))", 0, 1);
   
}

static double Duration(double m1, double m2, double rp, double ra)
{  
   static TF1 F("F", "exp(1.52631578947368429*log(x) + 0.513701609395389336*log(1+0.398026315789473673*x*x) - 1.5*log(1-x*x))", 0, 1);
   
   double tM = (m1+m2)*SM;
   double rM = m1*m2/(m1+m2)*SM;
   double e = (ra-rp)/(ra+rp);
   double a = (ra+rp)*(m1+m2)*DistUnit/2;
   
   
   double aux = C*a/g(e)/G/sqrt(tM);
   aux *= aux;
   aux *= aux;
   return F.Integral(0, e)*aux*15./304*C*G/rM;
 
}

double HT()
{  double tM = 2.8, rM = 0.7;
   double T0 = 7.75*3600;
   double e = 0.617; 
   double a0 = exp( log( G*tM*SM/(4*Pi*Pi/T0/T0)  )/3 );
   // e = 1- rp/a
   a0 /= tM*DistUnit;
   double rp = (1-e)*a0;
   double ra = 2*a0 - rp;
   return Duration(1.4, 1.4, rp, ra);
}

double DiffEq( double ee, double aa)
{
	double dade = ((12*aa*(1. + (73./24)*ee*ee + (37./96)*pow(ee,4.)))/(19.*ee*(1-ee*ee)*(1. + (121./304)*ee*ee)));
	return dade;
}

void GNGen::EvolveRa(double m1, double m2, double& rp, double& ra)
{  
   
   double a = (ra + rp)/2.;
   double e = (ra- rp)/(ra + rp);
   double totMass = m1+m2; 
   double da, de=.1;
   double k1, k2, k3, k4, k5, K=1.e10;
   
   while(fabs(1.-k1/K)>.01 || fabs(1.-k2/K)>.01 || fabs(1.-k3/K)>.01 || fabs(1-k4/K)>.01 || fabs(1.-k5/K)>.01) {
	   k1 = de * DiffEq(e, a);
	   k2 = de * DiffEq(e + de/3., a + k1/3.);
	   k3 = de * DiffEq(e + de/3., a + k1/6. + k2/6.);
	   k4 = de * DiffEq(e + de/2., a + k1/8. + 3.*k3/8.);
	   k5 = de * DiffEq(e + de, a - k1/2. -3*k3/2. + 2*k4);
	   da = (k1 + 4.*k4 + k5)/6.;
	   de /= 2.;
	   K = (k1 + k2 + k3 +k4 +k5)/5.;
   }
  
   while(rp>10.) {
	   k1 = de * DiffEq(e, a);
	   k2 = de * DiffEq(e + de/3., a + k1/3.);
	   k3 = de * DiffEq(e + de/3., a + k1/6. + k2/6.);
	   k4 = de * DiffEq(e + de/2., a + k1/8. + 3.*k3/8.);
	   k5 = de * DiffEq(e + de, a - k1/2. -3*k3/2. + 2*k4);
	   da = (k1 + 4.*k4 + k5)/6.;
       a -= da;
       e -= de;
       rp = (1. - e) * a;
       ra = (1. + e) * a;        
       if(da/a>.01) de /=2.;
       if(da/a<.001) de *=2.; 
      
       double freq = sqrt(0.5/rp)/(rp*Pi*totMass*TimeUnit);
       if(freq>freqCutoff) break;   
   }
}


double GNGen::generateEvent(double& m1, double& m2, double &rp, double& e)
{  const double dE_const = 85*Pi/12/sqrt(2);
   double r, ra, vra;
   double M , m, eta, rr, vofr, E0, rpm;
   
   do{
      dGammadmdMdr->GetRandom3(m1, m2, r);
   
      M = m1 + m2;
      m = m1*m2/M;
      eta = m/M;
   
      rr = r*PC/DistUnit;       // r in G=c=SM=1 units
      vofr = sqrt(smbhM/rr);    // v(r)
      E0 = m*vofr*vofr/2;
   
      rpm = pow(dE_const*eta*m/E0, 2./7)*M;
  }
  while(rpm<7.5*M);
   do{
      double rndnr = rnd.Uniform(); //printf("random nr = %lf\n", rndnr);
      rp = rndnr*rpm;
      //std::cout<<"m1: "<<m1<<" m2:"<<m2<<" rp/M:"<<rp/M<<" rndnr: "<<rndnr<<std::endl;
   }
   while(rp<7.5*M);
      
   double L0 = m*sqrt(2*rp*M + pow(vofr*rp, 2));      // bmw ; b = rp*sqrt(1 + 2GM/w^2rp)
   
   double deltaE = -dE_const*eta*m*pow(M/rp, 3.5);
   double deltaL = -6*Pi*pow(M*m/rp, 2);
   
   //accounting for both energy and angular momentum losses:
   double lom = (L0 + deltaL)/m;
   double aux = M/lom;
   double aux2 = sqrt( aux*aux + 2*(E0+deltaE)/m );
   
   ra = lom/(aux  - aux2);
   vra = lom/ra;
   rp = lom/(aux  + aux2);
   //double e0 = (ra - rp)(ra+rp);
   
   rp /= M;   // in units of total M
   ra /= M;   // in units of total M
   
   EvolveRa(m1, m2, rp, ra);
  
   //t = HT()/365.25/24./3600;
   e = (ra-rp)/(ra+rp);
   return Duration(m1, m2, rp, ra); 
}

void GNGen::generateEvents(int n, char* fn)
{  double m1, m2, rp, e;
   FILE* f = stdout;
   if(fn)f=fopen(fn, "w");
   for(int i=0; i<n; ++i){
      
      generateEvent(m1, m2, rp, e);
      //double q = m1/m2;
      //if(q>1) q = 1./q;
      
      fprintf(f, "%d %lf %lf %lf %.8lf\n", i, m1, m2, rp, e);
      //printf("%d %lf %lf %lf %lf %le %.8lf\n", i, m1, m2, ra, rp, vra, (ra-rp)/(ra+rp));
   }
   if(fn)fclose(f);
}
