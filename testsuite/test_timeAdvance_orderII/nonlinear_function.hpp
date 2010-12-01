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

#ifndef NLF_HPP
#define NLF_HPP

#include <life/lifecore/life.hpp>

// ===================================================
//! User functions
// ===================================================

namespace LifeV
{

double alpha = 1;
Real Pi2 = Pi*Pi;

class AnalyticalSol
{
public:
    inline Real operator()(Real t, Real x,Real y,Real z, UInt /*ic*/=0) const
    {
        return alpha*exp(-sin(Pi/2*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
    }
    inline Real grad(UInt icoor, Real t, Real x,Real y,Real z, UInt /*ic*/=0) const
    {
        switch (icoor)
        {
        case 1: // der_x
            return -alpha *Pi*exp(-sin(Pi/2*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
        case 2: // der_y
            return -alpha*Pi*exp(-sin(Pi/2*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
        case 3: // der_z
            return -alpha*Pi*exp(-sin(Pi/2*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
        default:
            return 0;
        }
    }
};

//solution on the boundary
Real uexact( const Real&  t ,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID&  icomp)
{
    return  alpha*exp(-sin(Pi/2*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
}

// u = alpha*exp(-sin(Pi/2*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z
// v = -Pi/2*cos(Pi/2*t)*u;
// w =  Pi2/4*sin(Pi/2*t)*u+Pi2/4cos(Pi/2*t)*u

Real source_in( const Real&  t ,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  icomp)
{
    return (3 + 1./4.*( sin(Pi/2*t) + cos(Pi/2 * t) * cos(Pi/2 * t) ) )
           * Pi2 * exp(-sin(Pi/2*t)) *cos(Pi*x) *cos(Pi*y) *cos(Pi*z);
}

Real d0( const Real&  t ,
         const Real& x,
         const Real& y,
         const Real& z,
         const ID&  icomp)
{
    return alpha*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
}

Real v0( const Real&  t ,
         const Real& x,
         const Real& y,
         const Real& z,
         const ID&  icomp)
{
    return -alpha*Pi/2*cos(Pi*x)*cos(Pi*y)*cos(Pi*z) * cos(Pi/2*t) * exp(-sin(Pi/2*t));
}

Real a0( const Real&  t ,
         const Real& x,
         const Real& y,
         const Real& z,
         const ID&  icomp)
{
    return + 1./4.*( sin(Pi/2*t) + cos(Pi/2 * t) * cos(Pi/2 * t) )
           * Pi2 * exp(-sin(Pi/2*t)) *cos(Pi*x) *cos(Pi*y) *cos(Pi*z);
}

Real UZero( const Real& /* t */,
            const Real& ,
            const Real& ,
            const Real& ,
            const ID&   )
{
    return 0.;
}



}
#endif
