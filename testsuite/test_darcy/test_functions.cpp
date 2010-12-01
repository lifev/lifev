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

/*!
    @file
    @brief A short description of the file content

    @author Alessio Fumagalli <fumagalli@mox4.mate.polimi.it>
    @date 24 Nov 2010

    A more detailed description of the file (if necessary)
 */

/*
  This file contain a test case with analytic solution for each of the following problems:
  1 - Darcy solver;
  2 - Darcy transient solver;
  3 - Darcy solver with non linear permeability tensor;
  4 - Darcy transient solver with non linear permeability tensor.
  The computational domain is the unit cube [0,1]^3.
*/

// For all the test we assume the following boundary flags
enum BCNAME
{
    /*
      BACK   = 1,
      FRONT  = 2,
      LEFT   = 3,
      RIGHT  = 4,
      BOTTOM = 5,
      TOP    = 6
    */

    LEFT   = 4,
    RIGHT  = 2,
    FRONT  = 1,
    BACK   = 3,
    TOP    = 6,
    BOTTOM = 5
};

// ===================================================
//!                     Darcy Solver
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
                     const Real& z,
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

// Inverse of permeability matrix
/* In this case the permeability matrix is
K = [2 1 0
     1 1 0
     0 0 1]
*/
Matrix inversePermeability( const Real& /*t*/,
                            const Real& /*x*/,
                            const Real& /*y*/,
                            const Real& /*z*/,
                            const std::vector<Real>& /*u*/)
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
Real source_in( const Real& /*t*/,
                const Real& /*x*/,
                const Real& /*y*/,
                const Real& /*z*/,
                const ID&   /*icomp*/)
{
    return -2.*x*x - 4.*y*y - 8.*x*y;
}

// Boundary conditions

// dp/dn = first_parameter + second_parameter * p
mixteBDfun.setFunctions_Mixte( mixte,
                               Members->getUOne() );

BCHandler bcDarcy( 6 );

bcDarcy.addBC( "Top",     TOP,     Natural,    Full,    neumann1, 1 );
bcDarcy.addBC( "Bottom",  BOTTOM,  Mixte,      Scalar,  mixteBDfun  );
bcDarcy.addBC(  "Left",   LEFT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Right",  RIGHT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Front",  FRONT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Back",    BACK,    Natural,    Full,    neumann2, 1 );

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
               const Real& z,
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
               const Real& z,
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
Real mixte( const Real& /* t */,
            const Real& x,
            const Real& y,
            const Real& z,
            const ID&   /*icomp*/)
{
    return -2.*y*x*x - 2*x*y*y - 6. + x*x*y*y + 6.*x + 5.*z;
}

// ===================================================
//!                 Darcy Transient Solver
// ===================================================

// Analytical solution
Real analyticalSolution( const Real& t,
                         const Real& x,
                         const Real& y,
                         const Real& z,
                         const ID& /*ic*/)
{
    return x*x*y*y*t*t + 6.*x + 5.*z*t;
}

// Analytical flux
Real analyticalFlux( const Real& t,
                     const Real& x,
                     const Real& y,
                     const Real& z,
                     const ID& icomp)
{
    switch (icomp)
    {
    case 1: // \frac{\partial }{\partial x}
        return -1. * (4.*x*y*y*t*t + 12. + 2.*x*x*y*t*t);
    case 2: // \frac{\partial }{\partial y}
        return -1. * (2.*x*y*y*t*t + 6. + 2.*x*x*y*t*t);
    case 3: // \frac{\partial }{\partial z}
        return -5.*t;
    default:
        return 0.;
    }
}

// Inverse of permeability matrix
/* In this case the permeability matrix is
K = [2 1 0
     1 1 0
     0 0 1]
*/
Matrix inversePermeability( const Real& /*t*/,
                            const Real& /*x*/,
                            const Real& /*y*/,
                            const Real& /*z*/,
                            const std::vector<Real>& /*u*/)
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

// Initial time primal variable
Real initialPrimal( const Real& /*t*/,
                    const Real& x,
                    const Real& /*y*/,
                    const Real& /*z*/,
                    const ID&   /*ic*/)
{
    return 6.*x;
}

// Source term
Real source_in( const Real& t,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  /*icomp*/)
{
    return -(4.*y*y*t*t + 8.*x*y*t*t + 2.*x*x*t*t ) + 2.*x*x*y*y*t + 5*z;
}

// Boundary conditions

// dp/dn = first_parameter + second_parameter * p
mixteBDfun.setFunctions_Mixte( mixte,
                               Members->getUOne() );

BCHandler bcDarcy( 6 );

bcDarcy.addBC( "Top",     TOP,     Natural,    Full,    neumann1, 1 );
bcDarcy.addBC( "Bottom",  BOTTOM,  Mixte,      Scalar,  mixteBDfun  );
bcDarcy.addBC(  "Left",   LEFT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Right",  RIGHT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Front",  FRONT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Back",    BACK,    Natural,    Full,    neumann2, 1 );

// Boundary condition of Dirichlet
Real dirichlet( const Real& t,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&   /*icomp*/)
{
    return x*x*y*y*t*t + 6.*x + 5.*z*t;
}

// Boundary condition of Neumann
Real neumann1( const Real& t,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   icomp)
{
    switch (icomp)
    {
    case 1:   //! Dx
        return -1.*(4.*x*y*y*t*t + 12. + 2.*x*x*y*t*t);
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

// Boundary condition of Neumann
Real neumann2( const Real& t,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   icomp)
{
    switch (icomp)
    {
    case 1:   //! Dx
        return (4.*x*y*y*t*t + 12. + 2.*x*x*y*t*t);
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
Real mixte( const Real& t,
            const Real& x,
            const Real& y,
            const Real& z,
            const ID&   /*icomp*/)
{
    return -(2.*x*y*y*t*t + 6. + 2.*x*x*y*t*t) + x*x*y*y*t*t + 6.*x + 5*z*t;
}

// ===================================================
//!                Darcy Solver Non Linear
// ===================================================

// Analytical solution
Real analyticalSolution( const Real& /*t*/,
                         const Real& x,
                         const Real& y,
                         const Real& z,
                         const ID& /*ic*/)
{
    return x*x*y*y + 6.*x + 5.*z;
}

// Analytical flux
Real analyticalFlux( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& z,
                     const ID& icomp)
{
    switch (icomp)
    {
    case 1: // \frac{\partial }{\partial x}
        return -1.*(((x*x*y*y + 6.*x + 5.*z)*(x*x*y*y + 6.*x + 5.*z) + 1) * (2.*x*y*y + 6.) + 2.*x*x*y);
    case 2: // \frac{\partial }{\partial y}
        return -1.*(2.*x*y*y + 6. + 2.*x*x*y);
    case 3: // \frac{\partial }{\partial z}
        return -10.;
    default:
        return 0.;
    }
}

// Non-linear inverse of permeability matrix.
/* In this case the permeability matrix is
K = [p^2+2 1   0
     1     1   0
     0     0   2]
*/
Matrix inversePermeability( const Real& /*t*/,
                            const Real& /*x*/,
                            const Real& /*y*/,
                            const Real& /*z*/,
                            const std::vector<Real>& u )
{

    Matrix inversePermeabilityMatrix( static_cast<UInt>(3), static_cast<UInt>(3) );

    // First row
    Real Entry00 = 1. / ( u[0] * u[0] + 1 );
    Real Entry01 = -1. / ( u[0] * u[0] + 1);
    Real Entry02 = 0.;

    // Second row
    Real Entry11 = ( u[0] * u[0] + 2 ) / ( u[0] * u[0] + 1 );
    Real Entry12 = 0.;

    // Third row
    Real Entry22 = 1. / 2.;

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

// Zero iteration primal variable
Real zeroItarationPrimal( const Real& /*t*/,
                          const Real& /*x*/,
                          const Real& /*y*/,
                          const Real& /*z*/,
                          const ID&   /*ic*/)
{
    return 1;
}

// Source term
Real source_in( const Real& /*t*/,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  /*icomp*/)
{
    return -1.*(2*(x*x*y*y+6*x+5*z)*(2*x*y*y+6)*(2*x*y*y+6)+2*((x*x*y*y+6*x+5*z)*(x*x*y*y+6*x+5*z)+2)*y*y+8*x*y+2*x*x);
}

// Boundary conditions

// dp/dn = first_parameter + second_parameter * p
mixteBDfun.setFunctions_Mixte( mixte,
                               Members->getUOne() );

BCHandler bcDarcy( 6 );

bcDarcy.addBC( "Top",     TOP,     Natural,    Full,    neumann1, 1 );
bcDarcy.addBC( "Bottom",  BOTTOM,  Mixte,      Scalar,  mixteBDfun  );
bcDarcy.addBC(  "Left",   LEFT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Right",  RIGHT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Front",  FRONT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Back",    BACK,    Natural,    Full,    neumann2, 1 );

// Boundary condition of Dirichlet
Real dirichlet( const Real& /*t*/,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&   /*icomp*/)
{
    return x*x*y*y + 6.*x + 5.*z;
}

// Boundary condition of Neumann
Real neumann1( const Real& /*t*/,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   icomp)
{
    switch (icomp)
    {
    case 1:   //! Dx
        return -1.*( ((x*x*y*y+6*x+5*z)*(x*x*y*y+6*x+5*z)+2)*(2*x*y*y+6)+2*x*x*y );
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

// Boundary condition of Neumann
Real neumann2( const Real& /*t*/,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   icomp)
{
    switch (icomp)
    {
    case 1:   //! Dx
        return ((x*x*y*y+6*x+5*z)*(x*x*y*y+6*x+5*z)+2)*(2*x*y*y+6)+2*x*x*y;
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
Real mixte( const Real& /*t*/,
            const Real& x,
            const Real& y,
            const Real& z,
            const ID&   /*icomp*/)
{
    return -1.*(2*x*y*y+6+2*x*x*y) +  x*x*y*y + 6.*x + 5.*z;
}

// ===================================================
//!           Darcy Solver Transient Non Linear
// ===================================================

// Analytical solution
Real analyticalSolution( const Real& t,
                         const Real& x,
                         const Real& y,
                         const Real& z,
                         const ID& /*ic*/)
{
    return x*x*y*y*t*t + 6.*x + 5.*z*t;
}

// Analytical flux
Real analyticalFlux( const Real& t,
                     const Real& x,
                     const Real& y,
                     const Real& z,
                     const ID& icomp)
{
    switch (icomp)
    {
    case 1: // \frac{\partial }{\partial x}
        return -1.*(((x*x*y*y*t*t + 6.*x + 5.*z*t)*(x*x*y*y*t*t + 6.*x + 5.*z*t) + 1) * (2.*x*y*y*t*t + 6.) + 2.*x*x*y*t*t);
    case 2: // \frac{\partial }{\partial y}
        return -1.*(2.*x*y*y*t*t + 6. + 2.*x*x*y*t*t);
    case 3: // \frac{\partial }{\partial z}
        return -10.*t;
    default:
        return 0.;
    }
}

// Inverse of permeability matrix
/* In this case the permeability matrix is
K = [p^2+2 1   0
     1     1   0
     0     0   2]
*/
Matrix inversePermeability( const Real& /*t*/,
                            const Real& /*x*/,
                            const Real& /*y*/,
                            const Real& /*z*/,
                            const std::vector<Real> & u )
{
    Matrix inversePermeabilityMatrix( static_cast<UInt>(3), static_cast<UInt>(3) );

    // First row
    Real Entry00 = 1. / ( u[0] * u[0] + 1 );
    Real Entry01 = -1. / ( u[0] * u[0] + 1 );
    Real Entry02 = 0.;

    // Second row
    Real Entry11 = ( u[0] * u[0] + 2) / ( u[0] * u[0] + 1);
    Real Entry12 = 0.;

    // Third row
    Real Entry22 = 1. / 2.;

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

// Initial time primal variable
Real initialPrimal( const Real& /*t*/,
                    const Real& x,
                    const Real& /*y*/,
                    const Real& /*z*/,
                    const ID&   /*ic*/)
{
    return 6.*x;
}

// Boundary conditions

// dp/dn = first_parameter + second_parameter * p
mixteBDfun.setFunctions_Mixte( mixte,
                               Members->getUOne() );

BCHandler bcDarcy( 6 );

bcDarcy.addBC( "Top",     TOP,     Natural,    Full,    neumann1, 1 );
bcDarcy.addBC( "Bottom",  BOTTOM,  Mixte,      Scalar,  mixteBDfun  );
bcDarcy.addBC(  "Left",   LEFT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Right",  RIGHT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Front",  FRONT,    Essential,  Scalar,  dirichlet   );
bcDarcy.addBC( "Back",    BACK,    Natural,    Full,    neumann2, 1 );

// Boundary condition of Dirichlet
Real dirichlet( const Real& t,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&   /*icomp*/)
{
    return x*x*y*y*t*t + 6.*x + 5.*z*t;
}

// Boundary condition of Neumann
Real neumann1( const Real& t,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   icomp)
{
    switch (icomp)
    {
    case 1:   //! Dx
        return -1.*(((x*x*y*y*t*t+6*x+5*z*t)*(x*x*t*t*y*y*t*t+6*x+5*z*t)+2)*(2*x*y*y*t*t+6)+2*x*x*y*t*t);
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

// Boundary condition of Neumann
Real neumann2( const Real& t,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   icomp)
{
    switch (icomp)
    {
    case 1:   //! Dx
        return (((x*x*y*y*t*t+6*x+5*z*t)*(x*x*t*t*y*y*t*t+6*x+5*z*t)+2)*(2*x*y*y*t*t+6)+2*x*x*y*t*t);
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
Real mixte( const Real& t,
            const Real& x,
            const Real& y,
            const Real& z,
            const ID&   /*icomp*/)
{
    return -(2*x*y*y*t*t+6+2*x*x*y*t*t) + x*x*y*y*t*t + 6.*x + 5*z*t;
}

// Source term
Real source_in( const Real& t,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  /*icomp*/)
{
    return -(2*(x*x*y*y*t*t+6*x+5*z*t)*(2*x*y*y*t*t+6)*(2*x*y*y*t*t+6)+2*((x*x*y*y*t*t+6*x+5*z*t)*(x*x*y*y*t*t+6*x+5*z*t)+2)*y*y*t*t+8*x*y*t*t+2*x*x*t*t ) + 2.*x*x*y*y*t + 5*z;
}
