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

#ifndef LF_HPP
#define LF_HPP

#include <lifev/core/LifeV.hpp>

// ===================================================
//! User functions
// ===================================================

namespace LifeV
{
const Real Pi = 3.14159265358979323846264338328;

class AnalyticalSol
{
public:
    inline Real operator() (Real t, Real x, Real y, Real z, UInt /*ic*/ = 0) const
    {
        return exp (-sin (Pi / 2 * t) ) * (x + y + z);
    }
    inline Real grad (UInt icoor, Real t, Real /*x*/, Real /*y*/, Real /*z*/, UInt /*ic*/ = 0) const
    {
        switch (icoor)
        {
            case 1: // der_x
                return exp (-sin (Pi / 2 * t) );
            case 2: // der_y
                return exp (-sin (Pi / 2 * t) );
            case 3: // der_z
                return exp (-sin (Pi / 2 * t) );
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
              const ID&  /*icomp*/)
{
    return exp (-sin (Pi / 2 * t) ) * (x + y + z);
}


Real source_in ( const Real&  t ,
                 const Real& x,
                 const Real& y,
                 const Real& z,
                 const ID&  /*icomp*/)
{
    return  -Pi / 2 * cos (Pi / 2 * t) * exp (-sin (Pi / 2 * t) ) * (x + y + z) ;
    //  return (3*Pi*Pi-alpha)*exp(-alpha*t)*sin(Pi*x+Pi/2)*sin(Pi*y+Pi/2)*sin(Pi*z+Pi/2);
}


Real d0       ( const Real& ,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  )
{
    return x + y + z;
}

Real v0 ( const Real& t ,
          const Real& x,
          const Real& y,
          const Real& z,
          const ID&)
{
    return -Pi / 2 * cos (Pi / 2 * t) * exp (-sin (Pi / 2 * t) ) * (x + y + z);
}

Real a0 ( const Real&  ,
          const Real& ,
          const Real& ,
          const Real& ,
          const ID&  )
{
    return 0;
}

Real UOne ( const Real& /* t */,
            const Real& ,
            const Real& ,
            const Real& ,
            const ID&   )
{
    return 1.;
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
