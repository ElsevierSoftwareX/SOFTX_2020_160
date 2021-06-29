
# INTRODUCTION

This package is a C++ version of the Python code developed by Sean McWilliams and others
( http://arxiv.org/abs/1212.0837) for eBBH simulations. To interface it with other code please use
the function

   getEBBH(double m1, double m2, double rmin0, double e0, wavearray<double>& Hp, wavearray<double>& Hx)

declared in 

   eBBH.hh 

The code requires the installation of an external package for ordinary differential equations integration.

NB. The distance from the source is always GM/c^2, so please rescale it accordingly.


# CVODE PACKAGE INSTALLATION

See [online manual](https://gwburst.gitlab.io/documentation/latest/html/install.html#cvode-installation)


# USAGE

External code needs to link/load the WAT (for wavearray), eBBH and two CVODE shared libraries.


# EVENT GENERATOR

There is also an event generator for Galactic Nuclei eccentric binary
formation as per http://arxiv.org/abs/0807.2638

First, create an object with the constructor:

   GNGen x(double mSMBH, double mmin, double mmax, double beta=2);
   ex.: GNGen x(3.5e6, 10, 25); 

mSMBH = mass of supermassive BH 
mmin, mmax = BH mass range
beta value is a model parameter, the paper explores values 2, 3
NB: masses are in solar mass units

Then, call

   x.generateEvent(double& m1, double& m2, double& rp, double& e);
   
to obtain one event. The return value is the merger time in seconds (leading order estimate). 

Alternatively, call

   x.generateEvents(int b, char* filename);
   
to generate a list of events to be printed on the screen (filename=0) or save in a file.

