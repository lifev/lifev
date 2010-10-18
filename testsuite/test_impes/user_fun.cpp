/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
       Date: 2010-07-29

  Copyright (C) 2010 EPFL, Politecnico di Milano

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/

/**
   @file user_fun.cpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2010-07-29
*/

#include "user_fun.hpp"

// ===================================================
//! User functions for the pressure equation
// ===================================================

// Inverse of permeability matrix
/* In this case the permeability matrix is
K = [2 1 0
     1 1 0
     0 0 1]
*/
Matrix pressurePermeability( const Real& /*t*/,
                             const Real& x,
                             const Real& y,
                             const Real& z )
{
    Matrix inversePermeabilityMatrix( static_cast<UInt>(3), static_cast<UInt>(3) );

    // First row
    Real Entry00 = 1.;
    Real Entry01 = -1.;
    Real Entry02 = 0.;

    // Second row
    Real Entry11 = 2.;
    Real Entry12 = 0.;

    // Third row
    Real Entry22 = 1.;

    // Fill in of the inversePermeabilityMatrix
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(0) ) = Entry00;
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(1) ) = Entry01;
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(2) ) = Entry02;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(0) ) = Entry01;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(1) ) = Entry11;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(2) ) = Entry12;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(0) ) = Entry02;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(1) ) = Entry12;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(2) ) = Entry22;

    return inversePermeabilityMatrix;
}

// Source term
Real pressureSource( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& /*z*/,
                     const ID&  /*icomp*/)
{
    return -2.*x*x - 4.*y*y - 8.*x*y;
}

// Boundary condition of Dirichlet
Real pressureDirichlet( const Real& /* t */,
                        const Real& x,
                        const Real& y,
                        const Real& z,
                        const ID&   /*icomp*/)
{
    return x*x*y*y + 6.*x + 5.*z;
}

// Boundary condition of Neumann
Real pressureNeumann( const Real& /* t */,
                       const Real& x,
                       const Real& y,
                       const Real& z,
                       const ID&   icomp)
{
	switch(icomp){
  		case 1:   //! Dx
            return  -1.*(4.*x*y*y + 2.*x*x*y + 12.);
		break;
		case 2:   //! Dy
            return 0.;
	 	break;
  		case 3:   //! Dz
    		return 0.;
    	break;
  	}
	return 0.;
}

// Boundary condition of Robin
Real pressureMixte( const Real& /* t */,
                    const Real& x,
                    const Real& y,
                    const Real& z,
                    const ID&   /*icomp*/)
{
    return -2.*y*x*x - 2*x*y*y - 6. + x*x*y*y + 6.*x + 5.*z;
}



// ===================================================
//! User functions for the saturation equation
// ===================================================

// Inverse of permeability matrix
/* In this case the permeability matrix is
K = [p^2+2 1   0
     1     1   0
     0     0   2]
*/
Matrix saturationPermeability( const Real& u,
                               const Real& t,
                               const Real& x,
                               const Real& y,
                               const Real& z )
{
    Matrix inversePermeabilityMatrix( static_cast<UInt>(3), static_cast<UInt>(3) );

    // First row
    Real Entry00 = 1./(u*u + 1);
    Real Entry01 = -1./(u*u + 1);
    Real Entry02 = 0.;

    // Second row
    Real Entry11 = (u*u + 2)/(u*u + 1);
    Real Entry12 = 0.;

    // Third row
    Real Entry22 = 1./2.;

    // Fill in of the inversePermeabilityMatrix
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(0) ) = Entry00;
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(1) ) = Entry01;
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(2) ) = Entry02;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(0) ) = Entry01;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(1) ) = Entry11;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(2) ) = Entry12;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(0) ) = Entry02;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(1) ) = Entry12;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(2) ) = Entry22;

    return inversePermeabilityMatrix;
}

// Physical flux function
Vector saturationPhysicalFlux( const Real& /*t*/,
                               const Real& x,
                               const Real& y,
                               const Real& z,
                               const Real& u )
{
    Vector physicalFluxVector( static_cast<UInt>(3) );

    // First row
    Real Entry0 = u;

    // Second row
    Real Entry1 = 0.;

    // Third row
    Real Entry2 = 0.;

    physicalFluxVector( static_cast<UInt>(0) ) = Entry0;
    physicalFluxVector( static_cast<UInt>(1) ) = Entry1;
    physicalFluxVector( static_cast<UInt>(2) ) = Entry2;

    return physicalFluxVector;
}

// First derivative in u of the physical flux function
Vector saturationFirstDerivativePhysicalFlux( const Real& /*t*/,
                                              const Real& x,
                                              const Real& y,
                                              const Real& z,
                                              const Real& u )
{
    Vector firstDerivativePhysicalFluxVector( static_cast<UInt>(3) );

    // First row
    Real Entry0 = 1.;

    // Second row
    Real Entry1 = 0.;

    // Third row
    Real Entry2 = 0.;

    firstDerivativePhysicalFluxVector( static_cast<UInt>(0) ) = Entry0;
    firstDerivativePhysicalFluxVector( static_cast<UInt>(1) ) = Entry1;
    firstDerivativePhysicalFluxVector( static_cast<UInt>(2) ) = Entry2;

    return firstDerivativePhysicalFluxVector;
}

// Source term
Real saturationSource( const Real& /*t*/,
                       const Real& x,
                       const Real& y,
                       const Real& /*z*/,
                       const ID&  /*icomp*/)
{
    return -2.*x*x - 4.*y*y - 8.*x*y;
}

// Initial condition
Real saturationInitialCondition( const Real& /* t */,
                                 const Real& x,
                                 const Real& y,
                                 const Real& z,
                                 const ID&   /*icomp*/ )
{
    return 1.;
}

// Boundary condition of Dirichlet
Real saturationDirichlet( const Real& /* t */,
                          const Real& x,
                          const Real& y,
                          const Real& z,
                          const ID&   /*icomp*/)
{
    return x*x*y*y + 6.*x + 5.*z;
}

// Boundary condition of Neumann
Real saturationNeumann( const Real& /* t */,
                        const Real& x,
                        const Real& y,
                        const Real& z,
                        const ID&   icomp)
{
	switch(icomp){
  		case 1:   //! Dx
            return  -1.*(4.*x*y*y + 2.*x*x*y + 12.);
		break;
		case 2:   //! Dy
            return 0.;
	 	break;
  		case 3:   //! Dz
    		return 0.;
    	break;
  	}
	return 0.;
}

// Boundary condition of Robin
Real saturationMixte( const Real& /* t */,
                      const Real& x,
                      const Real& y,
                      const Real& z,
                      const ID&   /*icomp*/)
{
    return -2.*y*x*x - 2*x*y*y - 6. + x*x*y*y + 6.*x + 5.*z;
}
