/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#include <vector>
#include <life/lifefem/elemOper_ext.hpp>
namespace LifeV
{
Real g1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 0.;
    break;
  }
  return 0.;
};

class Sphere: public Function {
 public:
  Real operator()(Real x, Real y, Real z) {
    Real s;
    s = sqrt( (x - 0.5) * (x - 0.5) + (y - 0.75) * (y - 0.75) + (z - 0.5) * (z - 0.5) ) - 0.20;

    return s;
  }
};
}

