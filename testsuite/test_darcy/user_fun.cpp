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

namespace dataProblem
{

// ===================================================
//!                    Problem data
// ===================================================

// Inverse of permeability matrix
/* In this case the permeability matrix is
K = [2 1 0
     1 1 0
     0 0 1]
*/
Matrix inversePermeability( const Real& /*t*/,
                            const Real& /*x*/,
                            const Real& /*y*/,
                            const Real& /*z*/, const std::vector<Real>& )
{
    Matrix inversePermeabilityMatrix( static_cast<UInt>(3), static_cast<UInt>(3) );

    /* Real Entry00, Entry01, Entry02, Entry11, Entry12, Entry22;

    if ( ((x > 1000) && (x < 1250) && (y > 0) &&  (y < 1250) )
         || ((x > 2250) && (x < 2500) && (y > 1000) &&  (y < 3000))
         || ((x > 3500) && (x < 3750) && (y > 500) &&  (y < 3000)) )
    {
        // First row
        Entry00 = 100.;
        Entry01 = 0.;
        Entry02 = 0.;

        // Second row
        Entry11 = 100.;
        Entry12 = 0.;

        // Third row
        Entry22 =  100.;

    }
    else
    { // First row
        Entry00 = 1.;
        Entry01 = 0.;
        Entry02 = 0.;

        // Second row
        Entry11 = 1.;
        Entry12 = 0.;

        // Third row
        Entry22 = 1.;
        }*/


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
Real source_in( const Real& /*t*/,
                const Real& x,
                const Real& y,
                const Real& /*z*/,
                const ID&  /*icomp*/)
{
    return -2.*x*x - 4.*y*y - 8.*x*y;//0.;
}

// Initial time primal variable for transient and non-linear transient solvers
Real initialPrimal( const Real& /*t*/,
                    const Real& x,
                    const Real& /*y*/,
                    const Real& /*z*/,
                    const ID& /*ic*/)
{
    return 6.*x;
}

// Zero iteration primal variable for non-linear solver
Real primalZeroIteration( const Real& /*t*/,
                          const Real& /*x*/,
                          const Real& /*y*/,
                          const Real& /*z*/,
                          const ID& /*ic*/)
{
    return 1.;
}

// Mass function for time dependent problem
Real mass( const Real& /*t*/,
           const Real& /*x*/,
           const Real& /*y*/,
           const Real& /*z*/,
           const ID&   /*ic*/)
{
    return 1.;
}


// ===================================================
//!                    Boundary data
// ===================================================

// Boundary condition of Dirichlet
Real dirichlet( const Real& /* t */,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&   /*icomp*/)
{
    return x*x*y*y + 6.*x + 5.*z;
}

// Boundary condition of Neumann
Real neumann1( const Real& /* t */,
               const Real& x,
               const Real& y,
               const Real& /*z*/,
               const ID&   icomp)
{
    switch (icomp)
    {
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

Real neumann2( const Real& /* t */,
               const Real& x,
               const Real& y,
               const Real& /*z*/,
               const ID&   icomp)
{
    switch (icomp)
    {
    case 1:   //! Dx
        return  4.*x*y*y + 2.*x*x*y + 12.;
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
Real robin( const Real& /* t */,
            const Real& x,
            const Real& y,
            const Real& z,
            const ID&   /*icomp*/)
{
    return -2.*y*x*x - 2*x*y*y - 6. + x*x*y*y + 6.*x + 5.*z;
}

// ===================================================
//!                 Analytical solution
// ===================================================

// Analytical solution
Real analyticalSolution( const Real&/*t*/,
                         const Real& x,
                         const Real& y,
                         const Real& z,
                         const ID& /*ic*/)
{

    return x*x*y*y + 6.*x + 5.*z;

}

// Gradient of the analytical solution
Real analyticalFlux( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& /*z*/,
                     const ID& icomp)
{

    switch ( icomp )
    {
    case 1:
        return -1. * (4.*x*y*y + 12. + 2.*x*x*y);

    case 2:
        return -1. * (2.*y*x*x + 2.*x*y*y + 6.);

    case 3:
        return -5.;

    default:
        return 0.;
    }

}

}
