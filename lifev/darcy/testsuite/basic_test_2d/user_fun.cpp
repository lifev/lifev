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

/**
    @file
    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author Anna Scotti <anna.scotti@mail.polimi.it>

    @date 2012-03-30
*/

#include "user_fun.hpp"

namespace dataProblem
{

// ===================================================
//!                    Problem data
// ===================================================

// Inverse of permeability matrix
/* In this case the permeability matrix is
K = [2 1
     1 2]
*/
Matrix inversePermeability ( const Real& /*t*/,
                             const Real& /*x*/,
                             const Real& /*y*/,
                             const Real& /*z*/, const std::vector<Real>& )
{
    Matrix invK ( static_cast<UInt>(2), static_cast<UInt>(2) );

    // First row
    const Real Entry00 = 2./3.;
    const Real Entry01 = -1./3.;

    // Second row
    const Real Entry11 = 2./3.;

    // Fill in of the inversePermeabilityMatrix
    invK ( static_cast<UInt>(0), static_cast<UInt>(0) ) = Entry00;
    invK ( static_cast<UInt>(0), static_cast<UInt>(1) ) = Entry01;
    invK ( static_cast<UInt>(1), static_cast<UInt>(0) ) = Entry01;
    invK ( static_cast<UInt>(1), static_cast<UInt>(1) ) = Entry11;

    return invK;

}

// Source term
Real source ( const Real& /*t*/,
              const Real& /*x*/,
              const Real& y,
              const Real& /*z*/,
              const ID&  /*icomp*/ )
{
    return -1. + 4. * y ;
}

// Vector source term
Vector vectorSource ( const Real& /*t*/,
                      const Real& x,
                      const Real& y,
                      const Real& /*z*/,
                      const ID& /*icomp*/ )
{
    Vector vectorSource( static_cast<UInt>(2) );

    const Real Entry0 = x - y;
    const Real Entry1 = y * y;

    vectorSource ( static_cast<UInt>(0) ) = Entry0;
    vectorSource ( static_cast<UInt>(1) ) = Entry1;

    return vectorSource;
}

// ===================================================
//!                    Boundary data
// ===================================================

// Boundary condition of Dirichlet
Real dirichlet ( const Real& /*t*/,
                 const Real& x,
                 const Real& y,
                 const Real& /*z*/,
                 const ID&   /*icomp*/ )
{
    return x * x + x * y - y * y;
}

// Boundary condition of Neumann
Real neumann1 ( const Real& /*t*/,
                const Real& x,
                const Real& y,
                const Real& /*z*/,
                const ID&   icomp )
{
    switch ( icomp )
    {
    case 0:   // Dx
        return  -3. * x - 2. * y + y * y;
        break;
    case 1:   // Dy
        return  0;
        break;
    }
    return 0.;
}

Real neumann2 ( const Real& /*t*/,
                const Real& x,
                const Real& y,
                const Real& /*z*/,
                const ID&   icomp )
{
    switch ( icomp )
    {
    case 0:   // Dx
        return -3. * x + 2. * y + 2. * y * y;
        break;
    case 1:   // Dy
        return 0.;
        break;
    }
    return 0.;
}

// Boundary condition of Robin
Real robin ( const Real& /* t */,
             const Real& x,
             const Real& y,
             const Real& /*z*/,
             const ID&   /*icomp*/ )
{
    return x * x + x * y - y * y - ( - 3. * x + 2. * y + 2. * y * y );
}

// ===================================================
//!                 Analytical solution
// ===================================================

// Analytical solution
Real analyticalSolution ( const Real& /*t*/,
                          const Real& x,
                          const Real& y,
                          const Real& /*z*/,
                          const ID& /*ic*/ )
{

    return x * x + x * y - y * y;

}

// Gradient of the analytical solution
Real analyticalFlux ( const Real& /*t*/,
                      const Real& x,
                      const Real& y,
                      const Real& /*z*/,
                      const ID& icomp )
{

    switch ( icomp )
    {
    case 0:
        return - 3. * x - 2. * y + y * y;

    case 1:
        return - 3. * x + 2. * y + 2. * y * y;

    default:
        return 0.;
    }

}

} // Namespace DataProblem
