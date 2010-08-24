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

Real Pi2 = Pi*Pi;

class AnalyticalSol
{
public:
	inline Real operator()(Real t, Real x,Real y,Real z, UInt /*ic*/=0) const
        {
          return exp(-sin(Pi/2*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);		
	}
	inline Real grad(UInt icoor, Real t, Real x,Real y,Real z, UInt /*ic*/=0) const {
		switch(icoor)
		{
  	  	case 1: 
		  return -Pi*exp(-sin(Pi/2*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
                case 2: 
                  return -Pi*exp(-sin(Pi/2*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
                case 3: 
                  return -Pi*exp(-sin(Pi/2*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
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
  return exp(-sin(Pi/2*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
}


Real source_in( const Real&  t ,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  icomp)
{
  return  (3 * Pi - 1./2.*cos(Pi/2*t) )*Pi*exp(-sin(Pi/2*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
}


Real d0 ( const Real&  t ,
	  const Real& x,
	  const Real& y,
	  const Real& z,
	  const ID&  icomp)
{
  return exp(-sin(Pi/2.*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z) ;
}

Real v0( const Real&  t ,
	 const Real& x,
	 const Real& y,
	 const Real& z,
	 const ID&  icomp)
{
  return Pi/2.*(cos(Pi*x)*cos(Pi*y)*cos(Pi*z) )* cos(Pi/2.*t) * exp(-sin(Pi/2.*t));
}

Real a0( const Real&  t ,
	 const Real& x,
	 const Real& y,
	 const Real& z,
	 const ID&  icomp)
{
  return Pi2 / 4*( sin(Pi/2*t)+cos(Pi/2*t)*cos(Pi/2*t) )*(cos(Pi*x)*cos(Pi*y)*cos(Pi*z) )*exp(-sin(Pi/2*t)) ;
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
