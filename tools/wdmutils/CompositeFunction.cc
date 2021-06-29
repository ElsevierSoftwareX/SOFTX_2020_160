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


#include "CompositeFunction.hh"
//#include "stdio.h"

template <class T>
CompositeFunction<T>::
CompositeFunction(const Function<T>& a, const Function<T>& b) : 
f(a.Clone()), g(b.Clone())
{}

template <class T>
CompositeFunction<T>::CompositeFunction(const CompositeFunction<T>& x) :
f(x.f.Clone()), g(x.g.Clone())
{}


template <class T> 
const Function<T>&  CompositeFunction<T>::Clone() const
{	return *new CompositeFunction<T>(f, g);
}

	
template <class T>
CompositeFunction<T>::~CompositeFunction()
{	//printf("CompositeFunction object destroyed\n");
}

template <class T>
T CompositeFunction<T>::operator()(const T& x) const
{	return f(g(x));
}


template <class T>
const CompositeFunction<T> operator|(const Function<T>& f, const Function<T>& g)
{	return CompositeFunction<T>(f, g);
}


template class CompositeFunction<double>;

template const CompositeFunction<double> 
operator|(const Function<double>& , const Function<double>&);
