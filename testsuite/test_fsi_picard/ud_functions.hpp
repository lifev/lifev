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
namespace LifeV
{
Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    return 0.;
}

Real u1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}

Real fZero(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}

// Initial velocity
Real u0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}


Real u2(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  //  Real pi = 3.14159265358979;
  switch(i) {
  case 1:
  case 2:
    return 0.0;
    break;
  case 3:
    if (t <= 0.003 )
      //return 6500*( 1.0 - cos( pi*t / 2.5e-3 )  );
      return 1.3332e4;
    else
      return 0.0;
    break;
  }
  return 0;
}


// Initial displacement and velocity
Real d0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  switch(i) {
  case 1:
    //    return -z*(z-5)*x/50;
    return 0.;
    break;
  case 2:
    //    return -z*(z-5)*y/50;
    return 0.;
    break;
  case 3:
    return 0.;
    break;
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }
}

Real w0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{

  switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return 0.0;
    break;
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }
}
}
