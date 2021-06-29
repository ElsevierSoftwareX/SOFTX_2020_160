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


#include "Function.hh"

template <class T>
Function<T>::~Function()
{}

template <class T>
T* Sample(const Function<T>& f, T a, T b, int n)
{	T* res = new T[n+1];
	T step = (b-a)/n;
	for(int i=0; i<=n; ++i) res[i] = f(step*i);
	return res;
}

template <class T>
void Sample(const Function<T>& f, T a, T b, int n, T* res)
{	T step = (b-a)/n;
	for(int i=0; i<=n; ++i) res[i] = f(step*i);
}


template <class T>
void AddSampling(const Function<T>& f, T a, T b, int n, T* res)
{	T step = (b-a)/n;
	for(int i=0; i<=n; ++i) res[i] += f(step*i);
}

template class Function<double>;

template double* Sample(const Function<double>&, double, double, int);
template void Sample(const Function<double>&, double, double, int, double*);
template void AddSampling(const Function<double>&, double, double, int, double*);
