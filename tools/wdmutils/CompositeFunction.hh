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


#ifndef COMPOSITEFUNCTION_HH
#define COMPOSITEFUNCTION_HH

#include "Function.hh"

// This class implements the concept of a composite function
// as in  f(g(x)); ex. the Meyer taper function cos ( nu (x) )


template <class T>
class CompositeFunction : public Function<T>{
public:
	// constructor creating the composite function f(g(.))
   // arg1 : function f
   // arg2 : function g
   CompositeFunction(const Function<T>&, const Function<T>&);
	
   // copy construction
   // arg1 : the Composite function to be copied
   CompositeFunction(const CompositeFunction<T>&);
   
   // clone
   // returns the cloned function
	const Function<T>& Clone() const;
   
   // destruction
	~CompositeFunction();
   
   // function evaluation
   // arg1: function argument
   // return value : function value at the given argument
	T operator()(const T&) const;
	
private:
	const Function<T>& f;
	const Function<T>& g;
};


// creating a composite function by overriding 
// the binary operator |
template <class T>
const CompositeFunction<T> operator|(const Function<T>& , const Function<T>&);

#endif
