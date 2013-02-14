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
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2010-07-29
*/
#include <iomanip>
#include "user_fun.hpp"

namespace dataProblem
{

// ===================================================
//!                    Problem data
// ===================================================

// Inverse of permeability matrix
/* In this case the permeability matrix is
K = [1 0
     0 1+p^2 ]
*/
Matrix inversePermeability::eval ( const UInt& iElem, const Vector3D& P, const Real& time ) const
{
    Matrix invK ( static_cast<UInt> (2), static_cast<UInt> (2) );

    const Real unkown_n = scalarField (0).eval ( iElem, P, time );

    // First row
    const Real Entry00 = 1.;
    const Real Entry01 = 0.;

    // Second row
    const Real Entry11 = 1. / ( unkown_n * unkown_n + 1. );

    // Fill in of the inversePermeabilityMatrix
    invK ( static_cast<UInt> (0), static_cast<UInt> (0) ) = Entry00;
    invK ( static_cast<UInt> (0), static_cast<UInt> (1) ) = Entry01;
    invK ( static_cast<UInt> (1), static_cast<UInt> (0) ) = Entry01;
    invK ( static_cast<UInt> (1), static_cast<UInt> (1) ) = Entry11;

    return invK;
}

// Reaction term
Real reactionTerm::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& /*time*/ ) const
{
    return 1;
}

// Scalar source term
Real scalarSource::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& time ) const
{
    const Real x (P[0]), y (P[1]), t (time);
    return  t * x * x - 2. * t * t - t -
            ( 6. * y + 2. * t * t * y ) * ( 1. + t * t * t * t * x * x * x * x +
                                            y * y * y * y * y * y +
                                            2. * t * t * x * x * y * y * y ) -
            ( 3. * y * y + t * t * y * y ) * ( 6. * y * y * y * y * y +
                                               6. * t * t * x * x * y * y ) +
            t * t * x * x + y * y * y;
}

// Vector source term
Vector vectorSource::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& time ) const
{
    const Real x (P[0]), y (P[1]), t (time);
    Vector source ( static_cast<UInt> (2) );

    const Real Entry0 = t * ( - x + y );
    const Real Entry1 = t * t * ( - y * y );

    source ( static_cast<UInt> (0) ) = Entry0;
    source ( static_cast<UInt> (1) ) = Entry1;

    return source;
}

// Initial time primal variable for transient and non-linear transient solvers
Real initialCondition::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& /*time*/ ) const
{
    const Real y (P[1]);
    return y * y * y;
}

// Mass function for time dependent problem
Real massFunction::eval ( const UInt& /*iElem*/, const Vector3D& /*P*/, const Real& /*time*/ ) const
{
    return 0.5;
}

// ===================================================
//!                    Boundary data
// ===================================================
void setBoundaryConditions ( bcHandlerPtr_Type& bcDarcy )
{

    BCFunctionBase dirichletBDfun;
    dirichletBDfun.setFunction ( dirichlet );

    bcDarcy->addBC ( "Top",    BCFlags::TOP,    Essential, Scalar, dirichletBDfun );
    bcDarcy->addBC ( "Bottom", BCFlags::BOTTOM, Essential, Scalar, dirichletBDfun );
    bcDarcy->addBC ( "Left",   BCFlags::LEFT,   Essential, Scalar, dirichletBDfun );
    bcDarcy->addBC ( "Right",  BCFlags::RIGHT,  Essential, Scalar, dirichletBDfun );

}

// Boundary condition of Dirichlet
Real dirichlet ( const Real& t,
                 const Real& x,
                 const Real& y,
                 const Real& z,
                 const ID&   icomp )
{
    return analyticalSolution ( t, x, y, z, icomp );
}

// ===================================================
//!                 Analytical solution
// ===================================================

// Analytical solution
Real analyticalSolution ( const Real& t,
                          const Real& x,
                          const Real& y,
                          const Real& z,
                          const ID& /*ic*/ )
{
    return  t * t * x * x + y * y * y;
}

// Gradient of the analytical solution
Real analyticalFlux ( const Real& t,
                      const Real& x,
                      const Real& y,
                      const Real& z,
                      const ID& icomp )
{
    switch ( icomp )
    {
        case 0:
            return -1. * ( 2. * t * t * x + t * ( x - y ) );
        case 1:
            return -1. * ( ( 3. * y * y + t * t * y * y ) *
                           ( 1. + t * t * t * t * x * x * x * x +
                             y * y * y * y * y * y +
                             2. * t * t * x * x * y * y * y ) );
        default:
            return 0.;
    }
}

} // Namespace DataProblem
