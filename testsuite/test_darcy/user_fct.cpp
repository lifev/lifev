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

#include "user_fct.hpp"

namespace LifeV
{
double zero(const double& t, const double& x, const double& y, const double& z, const ID& i) {
  return 0.;
}

double g1(const double& t, const double& x, const double& y, const double& z, const ID& i) {
  switch(i){
  case 1:
    return 1.;
    break;
  }
  return 0;
}

double g2(const double& t, const double& x, const double& y, const double& z, const ID& i) {
  switch(i){
  case 1:
    return 2*z; //2*z;
    break;
  }
  return 0;
}

double g3(const double& t, const double& x, const double& y, const double& z, const ID& i) {
  switch(i){
  case 1:
    return -1.;
    break;
  }
  return 0;
}

double mixte_coeff(const double& t, const double& x, const double& y, const double& z, const ID& i) {
  return 1.;
}
//////////////////////////////////////////////////////////////////////

//! Analytical solution function (0 bc on the hexahedron [0,1]x[0,2]x[0,5])
double Vfct(const double& t, const double& x, const double& y, const double& z, const ID& i)
{
  switch(i){
  case 1:
    return x*(1-x) * y*(2-y) * z*(5-z);
    break;
  }
  return 0.;
}

//! First derivatives
double VfctDer1(const double& t, const double& x, const double& y, const double& z, const ID& i)
{
  switch(i){
  case 1:   //! Dx
    return (1-2*x) * y*(2-y) * z*(5-z);
    break;
  case 2:   //! Dy
    return x*(1-x) * (2-2*y) * z*(5-z);
    break;
  case 3:   //! Dz
    return x*(1-x) * y*(2-y) * (5-2*z);
    break;
  }
  return 0.;
}

//! Second derivatives
double VfctDer2(const double& t, const double& x, const double& y, const double& z, const ID& i)
{
  switch(i){
  case 1:  //! Dxx
    return  -2 * y*(2-y) * z*(5-z);
    break;
  case 2:  //! Dyy
    return  - x*(1-x)* 2 * z*(5-z);
    break;
  case 3:  //! Dzz
    return - x*(1-x) * y*(2-y) * 2;
    break;
  }
  return 0.;
}

//!< := -laplace(Vfct)
double minusLaplaceVfct(const double& t, const double& x, const double& y, const double& z, const ID& i)
{
  return - ( VfctDer2(t,x,y,z,1) + VfctDer2(t,x,y,z,2) + VfctDer2(t,x,y,z,3) );
}
}
