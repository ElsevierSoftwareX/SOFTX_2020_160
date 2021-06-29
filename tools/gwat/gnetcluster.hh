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


/**********************************************************
 * Package:      graphic wat Class Library
 * File name:    gnecluster.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef GNETCLUSTER_HH
#define GNETCLUSTER_HH

#include "gwavearray.hh"
#include "netcluster.hh"

class gnetcluster : public netcluster 
{

public:
 
  // Constructors
  gnetcluster() : netcluster() {Init();}
  gnetcluster(netcluster&  nc);
  virtual ~gnetcluster();

  // operators
  gnetcluster& operator  = (const gnetcluster& nc){netcluster::operator=(nc);return *this;}
  gnetcluster& operator  = (const  netcluster& nc){netcluster::operator=(nc);return *this;}

private:

  void Init();

  ClassDef(gnetcluster,1)
};  

#endif
