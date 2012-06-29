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
K = [2 1 0
     1 1 0
     0 0 1]
*/
Matrix inversePermeability::eval ( const UInt& iElem, const Vector3D& P, const Real& time ) const
{
    Matrix invK ( static_cast<UInt>(2), static_cast<UInt>(2) );

//    Real unkown_n = scalarField(0).eval( iElem, P, time );

    // First row
    const Real Entry00 = 1.;
    const Real Entry01 = 0.;

    // Second row
    const Real Entry11 = 1.;
    //const Real Entry11 = 1./( unkown_n * unkown_n + 1.);

    // Fill in of the inversePermeabilityMatrix
    invK ( static_cast<UInt>(0), static_cast<UInt>(0) ) = Entry00;
    invK ( static_cast<UInt>(0), static_cast<UInt>(1) ) = Entry01;
    invK ( static_cast<UInt>(1), static_cast<UInt>(0) ) = Entry01;
    invK ( static_cast<UInt>(1), static_cast<UInt>(1) ) = Entry11;

    return invK;
}

// Reaction term
Real reactionTerm::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& /*time*/ ) const
{
    return 0;//1;
}

// Scalar source term
Real scalarSource::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& t/*ime*/ ) const
{
    const Real x(P[0]), y(P[1]);
    return  x;
}

// Vector source term
Vector vectorSource::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& time ) const
{
    Vector source( static_cast<UInt>(2) );

    const Real Entry0 = time * ( - P[0] + P[1] );
    const Real Entry1 = time * time * ( - P[1] * P[1] );

    source ( static_cast<UInt>(0) ) = Entry0;
    source ( static_cast<UInt>(1) ) = Entry1;

    return 0*source;
}

// Initial time primal variable for transient and non-linear transient solvers
Real initialCondition::eval ( const UInt& /*iElem*/, const Vector3D& /*P*/, const Real& /*time*/ ) const
{
    return 0.;
}

// Zero iteration primal variable for non-linear solver
Real nonLinearStartValue::eval ( const UInt& /*iElem*/, const Vector3D& /*P*/, const Real& /*time*/ ) const
{
    return 1.;
}

// Mass function for time dependent problem
Real massFunction::eval ( const UInt& /*iElem*/, const Vector3D& /*P*/, const Real& /*time*/ ) const
{
    return 1.;
}

// ===================================================
//!                    Boundary data
// ===================================================
void setBoundaryConditions ( bcHandlerPtr_Type & bcDarcy )
{

    BCFunctionBase dirichletBDfun, neumannBDfun1, neumannBDfun2;
    BCFunctionRobin robinBDfun;

    dirichletBDfun.setFunction ( dirichlet );
    neumannBDfun1.setFunction  ( neumann1 );
    neumannBDfun2.setFunction  ( neumann2 );
    // dp/dn = first_parameter + second_parameter * p
    robinBDfun.setFunctions_Robin( robin, robinMass );

    bcDarcy->addBC( "Top",    BCFlags::TOP,    Essential, Scalar, dirichletBDfun  );
    bcDarcy->addBC( "Bottom", BCFlags::BOTTOM, Essential, Scalar, dirichletBDfun  );
    bcDarcy->addBC( "Left",   BCFlags::LEFT,   Essential, Scalar, dirichletBDfun  );
    bcDarcy->addBC( "Right",  BCFlags::RIGHT,  Essential, Scalar, dirichletBDfun  );

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

// Boundary condition of Neumann
Real neumann1 ( const Real& t,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&   icomp )
{
    return analyticalFlux ( t, x, y, z, 0 );
}

Real neumann2 ( const Real& t,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&   icomp )
{
    return -1. * analyticalFlux ( t, x, y, z, 1 );
}

// Boundary condition of Robin
Real robin ( const Real& t,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID&   icomp )
{
    return -1. * analyticalFlux ( t, x, y, z, 2 ) +
            analyticalSolution ( t, x, y, z, icomp );
}

Real robinMass( const Real& /*t*/,
                const Real& /*x*/,
                const Real& /*y*/,
                const Real& /*z*/,
                const ID&   /*icomp*/ )
{
    return 1.;
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
    return  t*x;
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
        return -t;
    case 1:
        return 0.;
    default:
        return 0.;
    }
}

} // Namespace DataProblem
