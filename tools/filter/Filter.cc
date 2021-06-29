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


/***************************************************************************
                          CWBFilter.cpp  -  description
                             -------------------
    begin                : Jun 29 2011
    copyright            : (C) 2011 by Gabriele Vedovato
    email                : vedovato@lnl.infn.it
 ***************************************************************************/

/***************************************************************************

   mkfilter -- given n, compute recurrence relation
   to implement Butterworth, Bessel or Chebyshev filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 http://www-users.cs.york.ac.uk/~fisher/mkfilte

 ***************************************************************************/

#include "Filter.hh"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

CWB::Filter::Filter(char *instr) {

   if(instr==NULL) {cout << "CWB::Filter::Filter : NULL input params" << endl;exit(1);}

   int mkf_argc = 2;
   for(int i=0;i<16;i++) mkf_argv[i] = new char[32];
   sprintf(const_cast<char*>(mkf_argv[0]),"mkf");
   sprintf(const_cast<char*>(mkf_argv[1]),"-l");
   TObjArray* token = TString(instr).Tokenize(TString(' '));
   for(int i=0;i<token->GetEntries();i++) {
      TObjString* tok = (TObjString*)token->At(i);
      TString stok = tok->GetString(); 
      strcpy(const_cast<char*>(mkf_argv[mkf_argc++]),stok.Data());
      if (mkf_argc>=16) {cout << "CWB::Filter::Filter : Max Filters params = 14" << endl;exit(1);} 
   }
   mkf_argv[mkf_argc++]=NULL;  
   
   for(int i=0;i<mkf_argc;i++) {
     //cout << i << " mkf " << mkf_argv[i] << endl;
   }

   Init(mkf_argc,mkf_argv);
}

CWB::Filter::Filter(const int argc, const char *argv[]) {
   for(int i=0;i<16;i++) mkf_argv[i] = NULL;
   Init(argc,argv);
}

void
CWB::Filter::Init(const int argc, const char *argv[]) {

   readcmdline(const_cast<char**>(argv));
   checkoptions();
   setdefaults();

   if (options & opt_bc) return;	   /* bandcut (notch) */
   if (options & opt_lr) return;	   /* lorentzian ration */
   else
   { if (options & opt_re)
     { if (options & opt_bp) compute_bpres();	   /* bandpass resonator	 */
       if (options & opt_bs) compute_notch();	   /* bandstop resonator (notch) */
       if (options & opt_ap) compute_apres();	   /* allpass resonator		 */
     }
   else
   { if (options & opt_pi)
     { prewarp();
       splane.poles[0] = 0.0;
       splane.zeros[0] = -TWOPI * warped_alpha1;
       splane.numpoles = splane.numzeros = 1;
     }
     else
     { compute_s();
       prewarp();
       normalize();
     }
       if (options & opt_z) compute_z_mzt(); else compute_z_blt();
     }
   }
   if (options & opt_Z) add_extra_zero();
   expandpoly();
   computegain(const_cast<char**>(argv));
   //printresults(const_cast<char**>(argv));
   Reset();
}

CWB::Filter::~Filter(){
   for(int i=0;i<16;i++) if(mkf_argv[i] != NULL) delete [] mkf_argv[i];
}

void
CWB::Filter::Reset(){
  for (int i=0;i<MAXPZ+1;i++) xbp[i] = 0;	
  for (int i=0;i<MAXPZ+1;i++) ybp[i] = 0;
}

void
CWB::Filter::readcmdline(char *argv[])
  {
    options = order = polemask = 0;sample_rate = 0;width1 = 0;width2 = 0;
    int ap = 0;
    unless (argv[ap] == NULL) ap++; /* skip program name */
    until (argv[ap] == NULL)
    { 
        uint m = decodeoptions(argv[ap++]);
    	if (m & opt_ch) chebrip = getfarg(argv[ap++]);
    	if (m & opt_a)
  	  { raw_alpha1 = getfarg(argv[ap++]);
  	    raw_alpha2 = (argv[ap] != NULL && argv[ap][0] != '-') ? getfarg(argv[ap++]) : raw_alpha1;
  	  }
    	if (m & opt_Z) raw_alphaz = getfarg(argv[ap++]);
    	if (m & opt_o) order = getiarg(argv[ap++]);
    	if (m & opt_s) sample_rate = getfarg(argv[ap++]);
    	if (m & opt_p)
  	    { while (argv[ap] != NULL && argv[ap][0] >= '0' && argv[ap][0] <= '9')
  	      { int p = atoi(argv[ap++]);
  		      if (p < 0 || p > 31) p = 31; /* out-of-range value will be picked up later */
  		        polemask |= (1 << p);
  	      }
  	    }
  	  if (m & opt_re)
  	  { char *s = argv[ap++];
  	    if (s != NULL && seq(s,const_cast<char*>("Inf"))) infq = true;
  	    else { qfactor = getfarg(s); infq = false; }
  	  }
  	  if (m & opt_lr)
  	  { char *s = argv[ap++];
  	    if (s != NULL && seq(s,const_cast<char*>("Inf"))) infq = true;
  	    else { width1 = getfarg(s); infq = false; }
        char *a = argv[ap++];
  	    if (a != NULL && seq(a,const_cast<char*>("Inf"))) infq = true;
  	    else { width2 = getfarg(a); infq = false; }
  	  }
      options |= m;
    }
    if (sample_rate !=0) {raw_alpha1 /= sample_rate;raw_alpha2 /= sample_rate;}

    mkf_par.freq1=raw_alpha1*sample_rate;
    mkf_par.freq2=raw_alpha2*sample_rate;
    mkf_par.width1=width1;
    mkf_par.width2=width2;
    mkf_par.options=options;
    mkf_par.order=order;
    mkf_par.index=-1;
  }

uint
CWB::Filter::decodeoptions(char *s)
  { unless (*(s++) == '-') usage();
    uint m = 0;
    if (seq(s,const_cast<char*>("Be"))) m |= opt_be;
    else if (seq(s,const_cast<char*>("Bu"))) m |= opt_bu;
    else if (seq(s,const_cast<char*>("Ch"))) m |= opt_ch;
    else if (seq(s,const_cast<char*>("Re"))) m |= opt_re;
    else if (seq(s,const_cast<char*>("Lr"))) m |= opt_lr;
    else if (seq(s,const_cast<char*>("Pi"))) m |= opt_pi;
    else if (seq(s,const_cast<char*>("Lp"))) m |= opt_lp;
    else if (seq(s,const_cast<char*>("Hp"))) m |= opt_hp;
    else if (seq(s,const_cast<char*>("Bp"))) m |= opt_bp;
    else if (seq(s,const_cast<char*>("Bs"))) m |= opt_bs;
    else if (seq(s,const_cast<char*>("Ap"))) m |= opt_ap;
    else if (seq(s,const_cast<char*>("Bc"))) m |= opt_bc;
    else
      { until (*s == '\0')
	  { uint bit = optbit(*(s++));
	    if (bit == 0) usage();
	    m |= bit;
	  }
      }
    return m;
  }

uint
CWB::Filter::optbit(char c)
  { switch (c)
      { default:    return 0;
      	case 'a':   return opt_a;
      	case 'l':   return opt_l;
      	case 'o':   return opt_o;
      	case 's':   return opt_s;
      	case 'p':   return opt_p;
      	case 'w':   return opt_w;
      	case 'z':   return opt_z;
      	case 'Z':   return opt_Z;
      }
  }

double
CWB::Filter::getfarg(char *s)
  { if (s == NULL) usage();
    return atof(s);
  }

int
CWB::Filter::getiarg(char *s)
  { if (s == NULL) usage();
    return atoi(s);
  }

void
CWB::Filter::usage()
  { fprintf(stderr, "Mkfilter V.%s from <fisher@minster.york.ac.uk>\n", VERSION);
    fprintf(stderr, "Usage: mkfilter [-Be | -Bu | -Ch <r> | -Pi] [-Lp | -Hp | -Bp | -Bs] [-p <n1> <n2> ...] [-{lwz}] "
				     "[-Z <alphaz>] "
				     "-o <order> -a <alpha1> [ <alpha2> ]\n");
    fprintf(stderr, "   mkfilter -Re <q> [-Bp | -Bs | -Ap] [-l] -a <alpha>\n\n");
    fprintf(stderr, "  -Be, Bu             = Bessel, Butterworth\n");
    fprintf(stderr, "  -Ch <r>             = Chebyshev (r = dB ripple)\n");
    fprintf(stderr, "  -Pi                 = Proportional-Integral\n");
    fprintf(stderr, "  -Re <q>             = 2-pole resonator (q = Q-factor)\n");
    fprintf(stderr, "  -Lr <w1> <w2>       = Lorentzian ratio (w1 = Width1, w2 = Width2)\n");
    fprintf(stderr, "  -Lp, Hp, Bp, Bs, Ap, Bc = lowpass, highpass, bandpass, bandstop, allpass, bandcut\n");
    fprintf(stderr, "  -p                  = use listed poles only (ni = 0 .. order-1)\n");
    fprintf(stderr, "  -l                  = just list <order> parameters\n");
    fprintf(stderr, "  -w                  = don't pre-warp frequencies\n");
    fprintf(stderr, "  -z                  = use matched z-transform\n");
    fprintf(stderr, "  -Z                  = additional z-plane zero\n");
    fprintf(stderr, "  order = 1..%d;  alpha = f(corner)/f(sample)\n\n", MAXORDER);
    exit(1);
  }

bool optsok;

void
CWB::Filter::checkoptions()
  { optsok = true;
    unless (onebit(options & (opt_be | opt_bu | opt_ch | opt_re | opt_pi | opt_bc | opt_lr)))
      opterror(const_cast<char*>("must specify exactly one of -Be, -Bu, -Ch, -Re, -Pi -Bc -Lr"));
    if (options & opt_re)
      { unless (onebit(options & (opt_bp | opt_bs | opt_ap)))
	        opterror(const_cast<char*>("must specify exactly one of -Bp, -Bs, -Ap with -Re"));
	      if (options & (opt_lp | opt_hp | opt_o | opt_p | opt_w | opt_z))
	          opterror(const_cast<char*>("can't use -Lp, -Hp, -o, -p, -w, -z with -Re"));
      }
    else if (options & opt_lr)
      {
	      if (options & (opt_lp | opt_hp | opt_o | opt_p | opt_w | opt_z | opt_bp | opt_bs | opt_ap))
	          opterror(const_cast<char*>("can't use -Lp, -Hp, -Bp, -Bs, -Ap, -o, -p, -w, -z with -Lr"));
      }
    else if (options & opt_pi)
      { if (options & (opt_lp | opt_hp | opt_bp | opt_bs | opt_ap))
	          opterror(const_cast<char*>("-Lp, -Hp, -Bp, -Bs, -Ap illegal in conjunction with -Pi"));
	      unless ((options & opt_o) && (order == 1)) opterror(const_cast<char*>("-Pi implies -o 1"));
      }
    else if (options & opt_bc)
      { if (options & (opt_lp | opt_hp | opt_bp | opt_bs | opt_ap))
	          opterror(const_cast<char*>("-Lp, -Hp, -Bp, -Bs, -Ap illegal in conjunction with -Bc"));
	      if (options & opt_o) opterror(const_cast<char*>("-o illegal in conjunction with -Bc"));
      }
    else
      { unless (onebit(options & (opt_lp | opt_hp | opt_bp | opt_bs)))
	        opterror(const_cast<char*>("must specify exactly one of -Lp, -Hp, -Bp, -Bs"));
	      if (options & opt_ap) opterror(const_cast<char*>("-Ap implies -Re"));
	      if (options & opt_o)
	        { unless (order >= 1 && order <= MAXORDER) opterror(const_cast<char*>("order must be in range 1 .. %d"), MAXORDER);
	          if (options & opt_p)
	            { uint m = (1 << order) - 1; /* "order" bits set */
		            if ((polemask & ~m) != 0)
		              opterror(const_cast<char*>("order=%d, so args to -p must be in range 0 .. %d"), order, order-1);
	            }
	        }
	      else opterror(const_cast<char*>("must specify -o"));
      }
    unless (options & opt_a) opterror(const_cast<char*>("must specify -a"));
    unless (optsok) exit(1);
  }

void
CWB::Filter::opterror(char *msg, int p1, int p2)
  { fprintf(stderr, "mkfilter: "); fprintf(stderr, msg, p1, p2); putc('\n', stderr);
    optsok = false;
  }

void
CWB::Filter::setdefaults()
  { unless (options & opt_p) polemask = ~0; /* use all poles */
    unless (options & (opt_bp | opt_bs | opt_bc | opt_lr)) raw_alpha2 = raw_alpha1;
  }

void
CWB::Filter::compute_s() /* compute S-plane poles for prototype LP filter */
  { splane.numpoles = 0;
    if (options & opt_be)
      { /* Bessel filter */
	int p = (order*order)/4; /* ptr into table */
	if (order & 1) choosepole(bessel_poles[p++]);
	for (int i = 0; i < order/2; i++)
	  { choosepole(bessel_poles[p]);
	    choosepole(cconj(bessel_poles[p]));
	    p++;
	  }
      }
    if (options & (opt_bu | opt_ch))
      { /* Butterworth filter */
	for (int i = 0; i < 2*order; i++)
	  { double theta = (order & 1) ? (i*PI) / order : ((i+0.5)*PI) / order;
	    choosepole(expj(theta));
	  }
      }
    if (options & opt_ch)
      { /* modify for Chebyshev (p. 136 DeFatta et al.) */
	if (chebrip >= 0.0)
	  { fprintf(stderr, "mkfilter: Chebyshev ripple is %g dB; must be .lt. 0.0\n", chebrip);
	    exit(1);
	  }
	double rip = pow(10.0, -chebrip / 10.0);
	double eps = sqrt(rip - 1.0);
	double y = asinh(1.0 / eps) / (double) order;
	if (y <= 0.0)
	  { fprintf(stderr, "mkfilter: bug: Chebyshev y=%g; must be .gt. 0.0\n", y);
	    exit(1);
	  }
	for (int i = 0; i < splane.numpoles; i++)
	  { splane.poles[i].re *= sinh(y);
	    splane.poles[i].im *= cosh(y);
	  }
      }
  }

void
CWB::Filter::choosepole(_complex z)
  { if (z.re < 0.0)
      { if (polemask & 1) splane.poles[splane.numpoles++] = z;
	polemask >>= 1;
      }
  }

void
CWB::Filter::prewarp()
  { /* for bilinear transform, perform pre-warp on alpha values */
    if (options & (opt_w | opt_z))
      { warped_alpha1 = raw_alpha1;
	warped_alpha2 = raw_alpha2;
      }
    else
      { warped_alpha1 = tan(PI * raw_alpha1) / PI;
	warped_alpha2 = tan(PI * raw_alpha2) / PI;
      }
  }

void
CWB::Filter::normalize()		/* called for trad, not for -Re or -Pi */
  { double w1 = TWOPI * warped_alpha1;
    double w2 = TWOPI * warped_alpha2;
    /* transform prototype into appropriate filter type (lp/hp/bp/bs) */
    switch (options & (opt_lp | opt_hp | opt_bp| opt_bs))
      { case opt_lp:
	  { for (int i = 0; i < splane.numpoles; i++) splane.poles[i] = splane.poles[i] * w1;
	    splane.numzeros = 0;
	    break;
	  }

	case opt_hp:
	  { int i;
	    for (i=0; i < splane.numpoles; i++) splane.poles[i] = w1 / splane.poles[i];
	    for (i=0; i < splane.numpoles; i++) splane.zeros[i] = 0.0;	 /* also N zeros at (0,0) */
	    splane.numzeros = splane.numpoles;
	    break;
	  }

	case opt_bp:
	  { double w0 = sqrt(w1*w2), bw = w2-w1; int i;
	    for (i=0; i < splane.numpoles; i++)
	      { _complex hba = 0.5 * (splane.poles[i] * bw);
      		_complex temp = csqrt(1.0 - sqr(w0 / hba));
      		splane.poles[i] = hba * (1.0 + temp);
      		splane.poles[splane.numpoles+i] = hba * (1.0 - temp);
	      }
	    for (i=0; i < splane.numpoles; i++) splane.zeros[i] = 0.0;	 /* also N zeros at (0,0) */
	    splane.numzeros = splane.numpoles;
	    splane.numpoles *= 2;
	    break;
	  }

	case opt_bs:
	  { double w0 = sqrt(w1*w2), bw = w2-w1; int i;
	    for (i=0; i < splane.numpoles; i++)
	      { _complex hba = 0.5 * (bw / splane.poles[i]);
      		_complex temp = csqrt(1.0 - sqr(w0 / hba));
      		splane.poles[i] = hba * (1.0 + temp);
      		splane.poles[splane.numpoles+i] = hba * (1.0 - temp);
	      }
	    for (i=0; i < splane.numpoles; i++)	   /* also 2N zeros at (0, +-w0) */
	      { splane.zeros[i] = _complex(0.0, +w0);
		      splane.zeros[splane.numpoles+i] = _complex(0.0, -w0);
	      }
	    splane.numpoles *= 2;
	    splane.numzeros = splane.numpoles;
	    break;
	  }
      }
  }

void
CWB::Filter::compute_z_blt() /* given S-plane poles & zeros, compute Z-plane poles & zeros, by bilinear transform */
  { int i;
    zplane.numpoles = splane.numpoles;
    zplane.numzeros = splane.numzeros;
    for (i=0; i < zplane.numpoles; i++) zplane.poles[i] = blt(splane.poles[i]);
    for (i=0; i < zplane.numzeros; i++) zplane.zeros[i] = blt(splane.zeros[i]);
    while (zplane.numzeros < zplane.numpoles) zplane.zeros[zplane.numzeros++] = -1.0;
  }

_complex
CWB::Filter::blt(_complex pz)
  { return (2.0 + pz) / (2.0 - pz);
  }

void
CWB::Filter::compute_z_mzt() /* given S-plane poles & zeros, compute Z-plane poles & zeros, by matched z-transform */
  { int i;
    zplane.numpoles = splane.numpoles;
    zplane.numzeros = splane.numzeros;
    for (i=0; i < zplane.numpoles; i++) zplane.poles[i] = cexp(splane.poles[i]);
    for (i=0; i < zplane.numzeros; i++) zplane.zeros[i] = cexp(splane.zeros[i]);
  }

void
CWB::Filter::compute_notch()
  { /* compute Z-plane pole & zero positions for bandstop resonator (notch filter) */
    compute_bpres();		/* iterate to place poles */
    double theta = TWOPI * raw_alpha1;
    _complex zz = expj(theta);	/* place zeros exactly */
    zplane.zeros[0] = zz; zplane.zeros[1] = cconj(zz);
  }

void
CWB::Filter::compute_apres()
  { /* compute Z-plane pole & zero positions for allpass resonator */
    compute_bpres();		/* iterate to place poles */
    zplane.zeros[0] = reflect(zplane.poles[0]);
    zplane.zeros[1] = reflect(zplane.poles[1]);
  }

_complex
CWB::Filter::reflect(_complex z)
  { double r = hypot(z);
    return z / sqr(r);
  }

void
CWB::Filter::compute_bpres()
  { /* compute Z-plane pole & zero positions for bandpass resonator */
    zplane.numpoles = zplane.numzeros = 2;
    zplane.zeros[0] = 1.0; zplane.zeros[1] = -1.0;
    double theta = TWOPI * raw_alpha1; /* where we want the peak to be */
    if (infq)
      { /* oscillator */
      	_complex zp = expj(theta);
      	zplane.poles[0] = zp; zplane.poles[1] = cconj(zp);
      }
    else
      { /* must iterate to find exact pole positions */
      	_complex topcoeffs[MAXPZ+1]; expand(zplane.zeros, zplane.numzeros, topcoeffs);
      	double r = exp(-theta / (2.0 * qfactor));
      	double thm = theta, th1 = 0.0, th2 = PI;
      	bool cvg = false;
      	for (int i=0; i < 50 && !cvg; i++)
      	  { _complex zp = r * expj(thm);
      	    zplane.poles[0] = zp; zplane.poles[1] = cconj(zp);
      	    _complex botcoeffs[MAXPZ+1]; expand(zplane.poles, zplane.numpoles, botcoeffs);
      	    _complex g = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, expj(theta));
      	    double phi = g.im / g.re; /* approx to atan2 */
      	    if (phi > 0.0) th2 = thm; else th1 = thm;
      	    if (fabs(phi) < EPS) cvg = true;
      	    thm = 0.5 * (th1+th2);
      	  }
	      unless (cvg) fprintf(stderr, "mkfilter: warning: failed to converge\n");
      }
  }

void
CWB::Filter::add_extra_zero()
  { if (zplane.numzeros+2 > MAXPZ)
      { fprintf(stderr, "mkfilter: too many zeros; can't do -Z\n");
	      exit(1);
      }
    double theta = TWOPI * raw_alphaz;
    _complex zz = expj(theta);
    zplane.zeros[zplane.numzeros++] = zz;
    zplane.zeros[zplane.numzeros++] = cconj(zz);
    while (zplane.numpoles < zplane.numzeros) zplane.poles[zplane.numpoles++] = 0.0;	 /* ensure causality */
  }

void
CWB::Filter::expandpoly() /* given Z-plane poles & zeros, compute top & bot polynomials in Z, and then recurrence relation */
  { _complex topcoeffs[MAXPZ+1], botcoeffs[MAXPZ+1]; int i;
    expand(zplane.zeros, zplane.numzeros, topcoeffs);
    expand(zplane.poles, zplane.numpoles, botcoeffs);
    dc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, 1.0);
    double theta = TWOPI * 0.5 * (raw_alpha1 + raw_alpha2); /* "jwT" for centre freq. */
    fc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, expj(theta));
    hf_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, -1.0);
    for (i = 0; i <= zplane.numzeros; i++) xcoeffs[i] = +(topcoeffs[i].re / botcoeffs[zplane.numpoles].re);
    for (i = 0; i <= zplane.numpoles; i++) ycoeffs[i] = -(botcoeffs[i].re / botcoeffs[zplane.numpoles].re);
  }

void
CWB::Filter::expand(_complex pz[], int npz, _complex coeffs[])
  { /* compute product of poles or zeros as a polynomial of z */
    int i;
    coeffs[0] = 1.0;
    for (i=0; i < npz; i++) coeffs[i+1] = 0.0;
    for (i=0; i < npz; i++) multin(pz[i], npz, coeffs);
    /* check computed coeffs of z^k are all real */
    for (i=0; i < npz+1; i++)
      { if (fabs(coeffs[i].im) > EPS)
	  { fprintf(stderr, "mkfilter: coeff of z^%d is not real; poles/zeros are not complex conjugates\n", i);
	    exit(1);
	  }
      }
  }

void
CWB::Filter::multin(_complex w, int npz, _complex coeffs[])
  { /* multiply factor (z-w) into coeffs */
    _complex nw = -w;
    for (int i = npz; i >= 1; i--) coeffs[i] = (nw * coeffs[i]) + coeffs[i-1];
    coeffs[0] = nw * coeffs[0];
  }

void
CWB::Filter::GetZeros(int& numzeros, double *&zeros)
	{
	  if (options & (opt_bc|opt_lr)) return;	   /* bandcut|lorentzian_ratio (notch */
		numzeros = zplane.numzeros;
  	zeros = new double[numzeros+1];
   	for (int i=0;i<=numzeros;i++) zeros[i] = xcoeffs[i];
 	}

void
CWB::Filter::GetPoles(int& numpoles, double *&poles)
	{
	  if (options & (opt_bc|opt_lr)) return;	   /* bandcut|lorentzian_ratio (notch) */
		numpoles = zplane.numpoles;
  	poles = new double[numpoles+1];
   	for (int i=0;i<=numpoles;i++) poles[i] = ycoeffs[i];
 	}

double
CWB::Filter::GetGain()
	{
	  if (options & (opt_bc|opt_lr)) return 1;	   /* bandcut|lorentzian_ratio (notch) */
 		return gain;
  }

double
CWB::Filter::Arma(double value) {

  if (options & (opt_bc | opt_lr)) {	   /* bandcut/lorentzian_ratio (notch) */
    fprintf(stderr, "CWB::Filter::Arma Arma method not allowed with Band Cut option\n");
    exit(1);
  }

  for (int i=0;i<zplane.numzeros;i++) xbp[i] = xbp[i+1];	
  for (int i=0;i<zplane.numpoles;i++) ybp[i] = ybp[i+1];

  xbp[zplane.numzeros] = value / gain;
  ybp[zplane.numpoles] = 0;
  for (int i=0;i<zplane.numzeros;i++) ybp[zplane.numzeros] += xcoeffs[i]*xbp[i];
  ybp[zplane.numzeros] += xbp[zplane.numzeros];
  for (int i=0;i<zplane.numpoles;i++) ybp[zplane.numpoles] += ycoeffs[i]*ybp[i];

  return ybp[zplane.numpoles];

}

void
CWB::Filter::computegain(char *argv[])
  { if (options & opt_l)
      { /* just list parameters */
        //printcmdline(argv);
        _complex _gain = (options & opt_pi) ? hf_gain :
		       (options & opt_lp) ? dc_gain :
		       (options & opt_hp) ? hf_gain :
		       (options & (opt_bp | opt_ap)) ? fc_gain :
		       (options & opt_bs) ? csqrt(dc_gain * hf_gain) : _complex(1.0);
        gain = hypot(_gain);
      }
    else
      { printf("Command line: ");
      	printcmdline(argv);
      	printfilter();
      }
  }

void
CWB::Filter::printresults(char *argv[])
  { if (options & opt_l)
      { /* just list parameters */
        printcmdline(argv);
        _complex gain = (options & opt_pi) ? hf_gain :
		       (options & opt_lp) ? dc_gain :
		       (options & opt_hp) ? hf_gain :
		       (options & (opt_bp | opt_ap)) ? fc_gain :
		       (options & opt_bs) ? csqrt(dc_gain * hf_gain) : _complex(1.0);
       	printf("G  = %.10e\n", hypot(gain));
       	printcoeffs(const_cast<char*>("NZ"), zplane.numzeros, xcoeffs);
       	printcoeffs(const_cast<char*>("NP"), zplane.numpoles, ycoeffs);
      }
    else
      { printf("Command line: ");
      	printcmdline(argv);
      	printfilter();
      }
  }

void
CWB::Filter::printcmdline(char *argv[])
  { int k = 0;
    until (argv[k] == NULL)
      { if (k > 0) putchar(' ');
	      fputs(argv[k++], stdout);
      }
    putchar('\n');
 }

void
CWB::Filter::printcoeffs(char *pz, int npz, double coeffs[])
  { printf("%s = %d\n", pz, npz);
    for (int i = 0; i <= npz; i++) printf("%18.15e\n", coeffs[i]);
  }

void
CWB::Filter::printfilter()
  { printf("raw alpha1    = %14.10f\n", raw_alpha1);
    printf("raw alpha2    = %14.10f\n", raw_alpha2);
    unless (options & (opt_re | opt_w | opt_z | opt_lr))
      { printf("warped alpha1 = %14.15f\n", warped_alpha1);
	      printf("warped alpha2 = %14.15f\n", warped_alpha2);
      }
    printgain(const_cast<char*>("dc    "), dc_gain);
    printgain(const_cast<char*>("centre"), fc_gain);
    printgain(const_cast<char*>("hf    "), hf_gain);
    putchar('\n');
    unless (options & opt_re) printrat_s();
    printrat_z();
    printrecurrence();
  }

void
CWB::Filter::printgain(char *str, _complex gain)
  { double r = hypot(gain);
    printf("gain at %s:   mag = %15.9e", str, r);
    if (r > EPS) printf("   phase = %14.15f pi", atan2(gain) / PI);
    putchar('\n');
  }

void
CWB::Filter::printrat_s()	/* print S-plane poles and zeros */
  { printf("S-plane zeros:\n");
    printpz(splane.zeros, splane.numzeros);
    printf("S-plane poles:\n");
    printpz(splane.poles, splane.numpoles);
  }

void
CWB::Filter::printrat_z()	/* print Z-plane poles and zeros */
  { printf("Z-plane zeros:\n");
    printpz(zplane.zeros, zplane.numzeros);
    printf("Z-plane poles:\n");
    printpz(zplane.poles, zplane.numpoles);
  }

void
CWB::Filter::printpz(_complex *pzvec, int num)
  { int n1 = 0;
    while (n1 < num)
      { putchar('\t'); prcomplex(pzvec[n1]);
      	int n2 = n1+1;
      	while (n2 < num && pzvec[n2] == pzvec[n1]) n2++;
      	if (n2-n1 > 1) printf("\t%d times", n2-n1);
      	putchar('\n');
      	n1 = n2;
      }
    putchar('\n');
  }

void
CWB::Filter::printrecurrence() /* given (real) Z-plane poles & zeros, compute & print recurrence relation */
  { printf("Recurrence relation:\n");
    printf("y[n] = ");
    int i;
    for (i = 0; i < zplane.numzeros+1; i++)
      { if (i > 0) printf("     + ");
      	double x = xcoeffs[i];
      	double f = fmod(fabs(x), 1.0);
       	char *fmt = const_cast<char*>((f < EPS || f > 1.0-EPS) ? "%3g" : "%14.10f");
      	//char *fmt = (f < EPS || f > 1.0-EPS) ? "%3g" : "%14.10f";
      	putchar('('); printf(fmt, x); printf(" * x[n-%2d])\n", zplane.numzeros-i);
      }
    putchar('\n');
    for (i = 0; i < zplane.numpoles; i++)
      { printf("     + (%14.10f * y[n-%2d])\n", ycoeffs[i], zplane.numpoles-i);
      }
    putchar('\n');
  }

void
CWB::Filter::prcomplex(_complex z)
  { printf("%14.10f + j %14.15f", z.re, z.im);
  }


_complex 		
CWB::Filter::GetFrequencyResponseFunction(double f) {

	if (options & opt_bc) {	   /* bandcut (notch) */
    if ((f>=2*M_PI*raw_alpha1*sample_rate)&&(f<=2*M_PI*raw_alpha2*sample_rate)) {
      return _complex(0,0);
    } else {
      return _complex(1,0);
    }
  }
	if (options & opt_lr) {	   /* lorentzian_ratio (notch) */
    double rfreq1 = raw_alpha1*sample_rate;
    double rfreq2 = raw_alpha2*sample_rate;
    f/=2*M_PI;f=fabs(f);
    double value = sqrt((pow(f-rfreq1,2)+pow(width1,2))/(pow(f-rfreq2,2)+pow(width2,2)));
    if (value==0) value=1;
    //printf("%e %e\n",f,1/value);
    return _complex(value,0);
  }

	_complex F = _complex(0,-fabs(f/sample_rate));	
	_complex N = _complex(0,0);
	_complex D = _complex(0,0);
	_complex Y = _complex(0,0);
	
  for (int i=0;i<zplane.numzeros;i++) N = N+cexp(F*i)*xcoeffs[i];
  N = N+cexp(F*zplane.numzeros);
  for (int i=0;i<zplane.numpoles;i++) D = D+cexp(F*i)*ycoeffs[i];
  D = D-cexp(F*zplane.numpoles);

  if ((D.re != 0) || (D.im !=0)) Y=N/D;

  return cconj(Y/gain);
}

#ifdef _USE_WAT
wavearray<double>& operator >> (CWB::Filter& filter, wavearray<double>& x) {
  x=filter.data;
  return x;
}

CWB::Filter& operator >> (wavearray<double>& x, CWB::Filter& filter) {
  filter.data=x;
  for (int i=0;i<(int)x.size();i++) filter.data[i]=filter.Arma(x.data[i]);
  return filter;
}
#endif
