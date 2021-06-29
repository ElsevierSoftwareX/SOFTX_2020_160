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


#ifndef STDFUCTION_HH
#define STDFUCTION_HH

#include "Function.hh"


// this class abstracts an implemented C mathematical function 

class StdFunction : public Function<double>{
public:
	// constructor
   // arg1 : the C function abstracted by this class
   StdFunction(double (*)(double));
   
	// clone
   // return value: the cloned function
   const Function<double>& Clone() const;
	
   // function evaluation
   // arg1: function argument
   // return value : function value at the given argument
   double operator()(const double& x) const;
	
private:
	double (*f)(double);
};

#endif
