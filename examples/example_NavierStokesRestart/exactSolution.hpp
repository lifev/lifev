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

#include <life/lifecore/LifeV.hpp>

// ===================================================
//! User functions
// ===================================================

namespace LifeV
{
const Real Pi = 3.14159265358979323846264338328;
const Real R = 0.5; //radius
const Real mu = 0.035; //Dynamic viscosity (The density = 1.0 )
const Real Re = 300; //Reynols
const Real L = 10; //Length
const Real Vavg = ( Re * mu ) / ( 2 * R ); //Characteristic velocity
const Real Vmax = 2 * Vavg; //Maximum velocity
const Real DeltaP = - ( 8 * Vavg * mu * L ) / ( R * R ); //Pressure drop
const Real Konstant = -(1/(4*mu))*(DeltaP/L);


class AnalyticalSolVelocity
{
public:
    inline Real operator()(Real t, Real x,Real y,Real z, UInt /*ic*/=0) const
    {

      return Konstant*( R*R- (x*x + y*y) );

    }
    inline Real grad(UInt icoor, Real t, Real x,Real y,Real z, UInt /*ic*/=0) const
    {

      switch (icoor)
        {
        case 0: // der_x
	  return - 2 * Konstant * x;
        case 1: // der_y
	  return - 2 * Konstant * y; 
        case 2: // der_z
   	    return 0;
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
             const ID&  /*icomp*/)
{

  return Konstant*( R*R- (x*x + y*y) );
}


class AnalyticalSolPressure
{
public:
    inline Real operator()(Real t, Real x,Real y,Real z, UInt /*ic*/=0) const
    {      
      return (DeltaP)*(z-L);
    }
    inline Real grad(UInt icoor, Real t, Real x,Real y,Real z, UInt /*ic*/=0) const
    {

      switch (icoor)
        {
        case 0: // der_x
	  return 0;
        case 1: // der_y
	  return 0;
        case 2: // der_z
	  return DeltaP;
        default:
	  return 0;
        }
    }
};

//solution on the boundary
Real pexact( const Real&  t ,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID&  /*icomp*/)
{

  return (DeltaP)*(z-L);
}


Real source_in( const Real&  t ,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  /*icomp*/)
{
    return  -Pi/2*cos(Pi/2*t)*exp(-sin(Pi/2*t))*(x+y+z) ;
//  return (3*Pi*Pi-alpha)*exp(-alpha*t)*sin(Pi*x+Pi/2)*sin(Pi*y+Pi/2)*sin(Pi*z+Pi/2);
}


Real d0       ( const Real& ,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  )
{
    return x+y+z;
}

Real v0( const Real& t ,
         const Real& x,
         const Real& y,
         const Real& z,
         const ID&)
{
    return -Pi/2*cos(Pi/2*t)*exp(-sin(Pi/2*t))*(x+y+z);
}

Real a0( const Real&  ,
         const Real& ,
         const Real& ,
         const Real& ,
         const ID&  )
{
    return 0;
}

Real UOne( const Real& /* t */,
           const Real& ,
           const Real& ,
           const Real& ,
           const ID&   )
{
    return 1.;
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
