/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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


  Real gu1 (const Real & t, const Real & x, const Real & y, const Real & z,
        const ID & i, const Real & u)
  {
    return u * (z + y);        //convection
  }
  Real gu3 (const Real & t, const Real & x, const Real & y, const Real & z,
        const ID & i, const Real & u)
  {
    return 0.0;
  }
  Real gu2 (const Real & t, const Real & x, const Real & y, const Real & z,
        const ID & i, const Real & u)
  {
    if (t < 1.0)
      return 10.0 * t;        //not a d dirac
    return 10.0;
  }
  Real nu (Real t, Real x, Real y, Real z, Real u)
  {
    return 1.0;
  }
  Real mu (Real t, Real x, Real y, Real z, Real u)
  {
    return 1.0*(1+u*u/10000.0);
  }
  Real sigma (Real t, Real x, Real y, Real z, Real u)
  {
    return 0.0;
  }
  Real fct (Real t, Real x, Real y, Real z, Real u)
  {
    return 0.0;
  }




}
