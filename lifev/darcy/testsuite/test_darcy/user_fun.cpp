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
Matrix inversePermeability::eval ( const UInt& /*iElem*/, const Vector3D& /*P*/, const Real& /*time*/ ) const
{
    Matrix invK ( static_cast<UInt>(3), static_cast<UInt>(3) );

    // First row
    const Real Entry00 = 1.;
    const Real Entry01 = -1.;
    const Real Entry02 = 0.;

    // Second row
    const Real Entry11 = 2.;
    const Real Entry12 = 0.;

    // Third row
    const Real Entry22 = 1.;

    // Fill in of the inversePermeabilityMatrix
    invK ( static_cast<UInt>(0), static_cast<UInt>(0) ) = Entry00;
    invK ( static_cast<UInt>(0), static_cast<UInt>(1) ) = Entry01;
    invK ( static_cast<UInt>(0), static_cast<UInt>(2) ) = Entry02;
    invK ( static_cast<UInt>(1), static_cast<UInt>(0) ) = Entry01;
    invK ( static_cast<UInt>(1), static_cast<UInt>(1) ) = Entry11;
    invK ( static_cast<UInt>(1), static_cast<UInt>(2) ) = Entry12;
    invK ( static_cast<UInt>(2), static_cast<UInt>(0) ) = Entry02;
    invK ( static_cast<UInt>(2), static_cast<UInt>(1) ) = Entry12;
    invK ( static_cast<UInt>(2), static_cast<UInt>(2) ) = Entry22;

    return invK;
}

// Reaction term
Real reactionTerm::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& /*time*/ ) const
{
    return (P[0] * P[1] - 0.5 * P[2]);
}

// Scalar source term
Real scalarSource::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& /*time*/ ) const
{
    return - 4. * P[1] * P[1] + 4. * P[0] * P[0] - 8. * P[0] * P[1] + 6. +
         ( P[0] * P[1] - 0.5 * P[2] )  * ( P[0] * P[0] * P[1] * P[1] + 6. * P[0] + 5. * P[2] );
}

// Vector source term
Vector vectorSource::eval ( const UInt& /*iElem*/, const Vector3D& P, const Real& /*time*/ ) const
{
    Vector source( static_cast<UInt>(3) );

    const Real Entry0 = std::pow( P[0], 3 );
    const Real Entry1 = 2. * P[1];
    const Real Entry2 = 4. * P[2];

    source ( static_cast<UInt>(0) ) = Entry0;
    source ( static_cast<UInt>(1) ) = Entry1;
    source ( static_cast<UInt>(2) ) = Entry2;

    return source;
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
enum BCNAME
{

    // Flags for structured meshes
    BACK   = 3,
    FRONT  = 1,
    LEFT   = 4,
    RIGHT  = 2,
    BOTTOM = 5,
    TOP    = 6
};

void setBoundaryConditions ( boost::shared_ptr<BCHandler> & bcDarcy )
{

    BCFunctionBase dirichletBDfun, neumannBDfun1, neumannBDfun2;
    BCFunctionRobin robinBDfun;

    dirichletBDfun.setFunction ( dirichlet );
    neumannBDfun1.setFunction  ( neumann1 );
    neumannBDfun2.setFunction  ( neumann2 );
    // dp/dn = first_parameter + second_parameter * p
    robinBDfun.setFunctions_Robin( robin, robinMass );

    bcDarcy->addBC( "Top",     TOP,     Natural,    Full,    neumannBDfun1, 1 );
    bcDarcy->addBC( "Bottom",  BOTTOM,  Robin,      Scalar,  robinBDfun      );
    //bcDarcy->addBC(   "Top",    TOP,    Essential,  Scalar,  dirichletBDfun  );
    //bcDarcy->addBC("Bottom", BOTTOM,    Essential,  Scalar,  dirichletBDfun  );
    bcDarcy->addBC(  "Left",   LEFT,    Essential,  Scalar,  dirichletBDfun  );
    bcDarcy->addBC( "Right",  RIGHT,    Essential,  Scalar,  dirichletBDfun  );
    bcDarcy->addBC( "Back",  BACK,    Essential,  Scalar,  dirichletBDfun  );
    //bcDarcy->addBC(  "Front",   FRONT,    Essential,  Scalar,  dirichletBDfun  );
    bcDarcy->addBC( "Front",    FRONT,    Natural,    Full,    neumannBDfun2, 1 );

}

// Boundary condition of Dirichlet
Real dirichlet ( const Real& /*t*/,
                 const Real& x,
                 const Real& y,
                 const Real& z,
                 const ID&   /*icomp*/ )
{
    return x * x * y * y + 6. * x + 5. * z;
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
    case 0:   //! Dx
        return  -1. * ( 4. * x * y * y + 2. * x * x * y + 12. - 2. * x * x * x - 2. * y );
        break;
    case 1:   //! Dy
        return 0.;
        break;
    case 2:   //! Dz
        return 0.;
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
    case 0:   //! Dx
        return   ( 4. * x * y * y + 2. * x * x * y + 12. - 2. * x * x * x - 2. * y );
        break;
    case 1:   //! Dy
        return 0.;
        break;
    case 2:   //! Dz
        return 0.;
        break;
    }
    return 0.;
}

// Boundary condition of Robin
Real robin ( const Real& /*t*/,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID&   /*icomp*/ )
{
    return -2. * y * x * x - 2. * x * y * y - 6. + 2. * y + x * x * x + x * x * y * y + 6. * x + 5. * z;
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
Real analyticalSolution ( const Real& /*t*/,
                          const Real& x,
                          const Real& y,
                          const Real& z,
                          const ID& /*ic*/ )
{
    return  x * x * y * y + 6. * x + 5. * z;
}

// Gradient of the analytical solution
Real analyticalFlux ( const Real& /*t*/,
                      const Real& x,
                      const Real& y,
                      const Real& z,
                      const ID& icomp )
{
    switch ( icomp )
    {
    case 0:
        return -1. * ( 4. * x * y * y + 2. * x * x * y + 12. - 2. * x * x * x - 2. * y );

    case 1:
        return -1. * ( 2. * y * x * x + 2. * x * y * y + 6. - 2. * y - x * x * x );

    case 2:
        return -1. * ( 5. - 4. * z );

    default:
        return 0.;
    }
}

} // Namespace DataProblem
