//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

#ifndef NLF_HPP
#define NLF_HPP

#include <lifev/core/LifeV.hpp>

// ===================================================
//! User functions
// ===================================================

namespace LifeV
{

Real Pi2 = Pi* Pi;

class AnalyticalSol
{
public:
    inline Real operator() (Real t, Real x, Real y, Real z, UInt /*ic*/ = 0) const
    {
        return std::exp (-std::sin (Pi / 2 * t) ) * std::cos (Pi * x) * std::cos (Pi * y) * std::cos (Pi * z);
    }
    inline Real grad (UInt icoor, Real t, Real x, Real y, Real z, UInt /*ic*/ = 0) const
    {
        switch (icoor)
        {
            case 1:
                return -Pi * std::exp (-std::sin (Pi / 2 * t) ) * std::sin (Pi * x) * std::cos (Pi * y) * std::cos (Pi * z);
            case 2:
                return -Pi * std::exp (-std::sin (Pi / 2 * t) ) * std::cos (Pi * x) * std::sin (Pi * y) * std::cos (Pi * z);
            case 3:
                return -Pi * std::exp (-std::sin (Pi / 2 * t) ) * std::cos (Pi * x) * std::cos (Pi * y) * std::sin (Pi * z);
            default:
                return 0;
        }
    }
};

//solution on the boundary
Real uexact ( const Real&  t ,
              const Real& x,
              const Real& y,
              const Real& z,
              const ID&  icomp)
{
    return std::exp (-std::sin (Pi / 2 * t) ) * std::cos (Pi * x) * std::cos (Pi * y) * std::cos (Pi * z);
}


Real source_in ( const Real&  t ,
                 const Real& x,
                 const Real& y,
                 const Real& z,
                 const ID&  icomp)
{
    return  (3 * Pi - 1. / 2.*std::cos (Pi / 2 * t) ) * Pi * std::exp (-std::sin (Pi / 2 * t) ) * std::cos (Pi * x) * std::cos (Pi * y) * std::cos (Pi * z);
}


Real d0 ( const Real&  t ,
          const Real& x,
          const Real& y,
          const Real& z,
          const ID&  icomp)
{
    return std::exp (-std::sin (Pi / 2.*t) ) * std::cos (Pi * x) * std::cos (Pi * y) * std::cos (Pi * z) ;
}

Real v0 ( const Real&  t ,
          const Real& x,
          const Real& y,
          const Real& z,
          const ID&  icomp)
{
    return Pi / 2.* (std::cos (Pi * x) * std::cos (Pi * y) * std::cos (Pi * z) ) * std::cos (Pi / 2.*t) * std::exp (-std::sin (Pi / 2.*t) );
}

Real a0 ( const Real&  t ,
          const Real& x,
          const Real& y,
          const Real& z,
          const ID&  icomp)
{
    return Pi2 / 4 * ( std::sin (Pi / 2 * t) + std::cos (Pi / 2 * t) * std::cos (Pi / 2 * t) ) * (std::cos (Pi * x) * std::cos (Pi * y) * std::cos (Pi * z) ) * std::exp (-std::sin (Pi / 2 * t) ) ;
}

Real UZero ( const Real& /* t */,
             const Real& ,
             const Real& ,
             const Real& ,
             const ID&   )
{
    return 0.;
}

}
#endif
