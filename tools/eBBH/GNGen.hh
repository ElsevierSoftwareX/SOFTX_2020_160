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


#ifndef GNGEN_HH
#define GNGEN_HH

#include "TRandom3.h"
#include "TH3F.h"

class GNGen{
public:
   // all masses in units of solar mass
   GNGen(double mSMBH, double mmin, double mmax, double beta=2); 
   GNGen(const GNGen& x);
   ~GNGen();
   
   // 
   void setFreqCutoff(double f);
   
   // prints or saves in a file a list of m1, m2, rp, e values to be used by PN codes! 
   void generateEvents(int n, char* fn=0);  
   
   // generates m1, m2, rp, e (with evolution of orbit to higher frequency)
   // returns leading order merger time estimate (not very accurate)
   double generateEvent(double& m1, double& m2, double& rp, double& e);
   
private:
   double minM, maxM, beta, smbhM;
   TH3F* dGammadmdMdr;
   TRandom3 rnd;
   double freqCutoff;
   void EvolveRa(double m1, double m2, double& rp, double& ra);
};

#endif
