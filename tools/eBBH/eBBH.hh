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


#ifndef EBBH
#define EBBH

#include "wavearray.hh"

int getEBBH(double m1, double m2, double rmin0, double e0, wavearray<double>& Hp, wavearray<double>& Hx, double t_end=10000000.);


#endif
