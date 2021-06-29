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


#ifndef FUNCTION_HH
#define FUNCTION_HH


// this class abstracts functions in a most limited way
// but it will make life easier later

template <class T>
class Function{
public:
	virtual const Function<T>& Clone() const  = 0;
	virtual ~Function();
   
   // function evaluation
   // arg1 : the point at which the function is evaluated
   // returns the function value at that point
	virtual T operator()(const T&) const = 0;
};


// samples a function
// arg1 : the function to be sampled
// arg2, arg3 : the sampling interval [a,b] 
// arg4 : how many subintervals are defined by this sampling ( #points - 1)
// return value : the sampled values
template <class T>
T* Sample(const Function<T>& f, T a, T b, int n);  //x[0] = a; x[n]=b


// samples a function
// arg1 : the function to be sampled
// arg2, arg3 : the sampling interval [a,b] 
// arg4 : how many subintervals are defined by this sampling ( #points - 1)
// arg5 : stores the sampled values; must be preallocated by the calling function
template <class T>
void Sample(const Function<T>& f, T a, T b, int n, T* res); // x[0] = a; x[n]=b

// same as above but the sampled values are added 
// to the values present in the vector given in the fifth argument
template <class T>
void AddSampling(const Function<T>& f, T a, T b, int n, T* res);


#endif
