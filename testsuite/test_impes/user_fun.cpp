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
namespace dataProblem
{

namespace dataPhysical
{

// Porosity
const Real Phi ( const Real& x, const Real& y, const Real& z )
{
    return 0.4;
}

// Relative permeability for the wetting phase
const Real k_rw ( const Real& S_w )
{
    // Define the effective saturation
    const Real barS_w = (S_w - S_wr) / (1. - S_wr - S_nr);

    return pow( barS_w, (2. + 3.*lambda) / lambda );
}

// First derivative of the relative permeability for the wetting phase
const Real Dk_rw ( const Real& S_w )
{
    // Define the effective saturation
    const Real barS_w = (S_w - S_wr) / (1. - S_wr - S_nr);

    return (2. + 3.*lambda) / lambda * pow(barS_w, (2. + 3.*lambda) / lambda - 1.) / (1. - S_wr - S_nr);
}

// Relative permeability for the non-wetting phase
const Real k_rn ( const Real& S_n )
{
    // Define the effective saturation
    const Real barS_n = (S_n - S_nr) / (1. - S_wr - S_nr);

    return pow( barS_n, 2.) * (1 - pow( 1 - barS_n, (2. + lambda) / lambda));
}

// First derivative of the relative permeability for the non-wetting phase
const Real Dk_rn ( const Real& S_n )
{
    // Define the effective saturation
    const Real barS_n = (S_n - S_nr) / (1. - S_wr - S_nr);

    return ( 2.*barS_n * (1 - pow( 1 - barS_n, (2. + lambda) / lambda)) -
             pow( barS_n, 2.) * (2. + lambda) / lambda * pow(1 - barS_n, (2. + lambda) / lambda - 1.) ) /
           (1. - S_wr - S_nr);
}

// Capillary pressure
const Real pc ( const Real& S_w ) // [Pa]
{
    // Define the effective saturation
    const Real barS_w = (S_w - S_wr) / (1. - S_wr - S_nr);

    return pd * pow(barS_w, -1. / lambda);
}

// First derivative of the capillary pressure
const Real Dpc ( const Real& S_w ) // [Pa]
{
    // Define the effective saturation
    const Real barS_w = (S_w - S_wr) / (1. - S_wr - S_nr);

    return pd / lambda * pow(barS_w, -1. / lambda - 1.) / (1. - S_wr - S_nr);
}

// Absolute inverse permeability
const Matrix invK ( const Real& t, const Real& x, const Real& y, const Real& z ) // [m^2]
{
    Matrix inversePermeabilityMatrix( static_cast<UInt>(3), static_cast<UInt>(3) );

    const Real highInvPermeability = 1e5;
    const Real lowInvPermeability = 1.;

    Real Entry00, Entry01, Entry02, Entry11, Entry12, Entry22;

    if ( ((x > 1000) && (x < 1250) && (y > 0) &&  (y < 1250) )
            || ((x > 2250) && (x < 2500) && (y > 1000) &&  (y < 3000))
            || ((x > 3500) && (x < 3750) && (y > 500) &&  (y < 3000)) )
    {
        // First row
        Entry00 = highInvPermeability;
        Entry01 = 0.;
        Entry02 = 0.;

        // Second row
        Entry11 = highInvPermeability;
        Entry12 = 0.;

        // Third row
        Entry22 =  highInvPermeability;

    }
    else
    {
        // First row
        Entry00 = lowInvPermeability;
        Entry01 = 0.;
        Entry02 = 0.;

        // Second row
        Entry11 = lowInvPermeability;
        Entry12 = 0.;

        // Third row
        Entry22 = lowInvPermeability;
    }

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


}


// Inverse of permeability matrix
/* In this case the permeability matrix is
   K = [2 1 0
        1 1 0
        0 0 1]
*/
Matrix pressurePermeability( const Real& t,
                             const Real& x,
                             const Real& y,
                             const Real& z,
                             const std::vector<Real> & u )
{

    // Alias for the non-wetting phase saturation
    const Real& S_n = u[0];

    // Compute the phase mobility
    const Real lambda_w = dataPhysical::k_rw( 1. - S_n ) / dataPhysical::mu_w;
    const Real lambda_n = dataPhysical::k_rn( S_n ) / dataPhysical::mu_n;

    return - dataPhysical::invK(t, x, y, z) / ( lambda_n + lambda_w );
}

// Source term
Real pressureSource( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& /*z*/,
                     const ID&  /*icomp*/)
{
    return 0.;
}

// Boundary condition of Dirichlet
Real pressureDirichlet1( const Real& /* t */,
                         const Real& /* x */,
                         const Real& /* y */,
                         const Real& /* z */,
                         const ID&   /*icomp*/)
{
    return 11000000.;
}

// Boundary condition of Dirichlet
Real pressureDirichlet2( const Real& /* t */,
                         const Real& /* x */,
                         const Real& /* y */,
                         const Real& /* z */,
                         const ID&   /*icomp*/)
{
    return 10000000.;
}

// Boundary conditisaturationDirichletBDfunon of Dirichlet
Real pressureDirichlet3( const Real& /* t */,
                         const Real& /* x */,
                         const Real& /* y */,
                         const Real& /* z */,
                         const ID&   /*icomp*/)
{
    return 10500000.;
}

// Boundary condition of Neumann
Real pressureNeumann( const Real& /* t */,
                      const Real& x,
                      const Real& y,
                      const Real& z,
                      const ID&   icomp)
{
    return 0.;

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

// Boundary condition of Robin
Real pressureMixte( const Real& /* t */,
                    const Real& /* x */,
                    const Real& /* y */,
                    const Real& /* z */,
                    const ID&   /*icomp*/)
{
    return 0.;
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
Matrix saturationPermeability( const Real& t,
                               const Real& x,
                               const Real& y,
                               const Real& z,
                               const std::vector<Real> & u )
{
    // Alias for the non-wetting phase saturation
    const Real& S_n = u[0];

    // Compute the phase mobility
    const Real lambda_w = dataPhysical::k_rw( 1. - S_n ) / dataPhysical::mu_w;
    const Real lambda_n = dataPhysical::k_rn( S_n ) / dataPhysical::mu_n;

    // Compute the fractional flow
    const Real f_n = lambda_n / ( lambda_w + lambda_n );

    return - dataPhysical::invK(t, x, y, z) / ( lambda_w * f_n * dataPhysical::Dpc ( 1. - S_n ) ) ;
}

// Physical flux function
Vector saturationPhysicalFlux( const Real& t,
                               const Real& x,
                               const Real& y,
                               const Real& z,
                               const std::vector<Real>& u )
{
    Vector physicalFluxVector( static_cast<UInt>(3) );

    // Alias for the non-wetting phase saturation
    const Real& S_n = u[0];

    // Compute the phase mobility
    const Real lambda_w = dataPhysical::k_rw( 1. - S_n ) / dataPhysical::mu_w;
    const Real lambda_n = dataPhysical::k_rn( S_n ) / dataPhysical::mu_n;

    // Compute the fractional flow
    const Real f_n = lambda_n / ( lambda_w + lambda_n );

    // Compute the last column of the inverse of the inverse permeability
    const Matrix invK = dataPhysical::invK( t, x, y, z );

    // Compute the denominator of the last column of the inverse of the inverse permeability
    const Real denominator = invK(0, 0) * invK(1, 1) * invK(2, 2)
                             - invK(0, 0) * invK(1, 2) * invK(1, 2)
                             - invK(0, 1) * invK(0, 1) * invK(2, 2)
                             + 2. * invK(0, 1) * invK(0, 2) * invK(1, 2)
                             - invK(0, 2) * invK(0, 2) * invK(1, 1);

    // Compute the first component of the last column of the inverse of the inverse permeability
    const Real K02 = ( invK(0, 1) * invK(1, 2) - invK(0, 2) * invK(1, 1) ) / denominator;

    // Compute the second component of the last column of the inverse of the inverse permeability
    const Real K12 = ( invK(0, 0) * invK(1, 2) - invK(0, 1) * invK(0, 2) ) / denominator;

    // Compute the third component of the last column of the inverse of the inverse permeability
    const Real K22 = ( invK(0, 0) * invK(1, 1) - invK(0, 1) * invK(0, 1) ) / denominator;

    // Compute a common constant
    const Real lfrg = lambda_w * f_n * (dataPhysical::rho_w - dataPhysical::rho_n) * dataPhysical::g;

    // First row
    const Real Entry0 = f_n * u[1] - lfrg * K02;

    // Second row
    const Real Entry1 = f_n * u[2] - lfrg * K12;

    // Third row
    const Real Entry2 = f_n * u[3] - lfrg * K22;

    physicalFluxVector( static_cast<UInt>(0) ) = Entry0;
    physicalFluxVector( static_cast<UInt>(1) ) = Entry1;
    physicalFluxVector( static_cast<UInt>(2) ) = Entry2;

    return physicalFluxVector;
}

// First derivative in u of the physical flux function
Vector saturationFirstDerivativePhysicalFlux( const Real& t,
                                              const Real& x,
                                              const Real& y,
                                              const Real& z,
                                              const std::vector<Real>& u )
{
    Vector firstDerivativePhysicalFluxVector( static_cast<UInt>(3) );

    // Alias for the non-wetting phase saturation
    const Real& S_n = u[0];

    // Compute the phase mobility
    const Real lambda_w = dataPhysical::k_rw( 1. - S_n ) / dataPhysical::mu_w;
    const Real lambda_n = dataPhysical::k_rn( S_n ) / dataPhysical::mu_n;

    // Compute the first derivative of the phase mobility
    const Real Dlambda_w = dataPhysical::Dk_rw( 1. - S_n ) / dataPhysical::mu_w;
    const Real Dlambda_n = dataPhysical::Dk_rn( S_n ) / dataPhysical::mu_n;

    // Compute the fractional flow
    const Real f_n = lambda_n / ( lambda_w + lambda_n );

    // Compute the first derivative of the fractional flow
    const Real Df_n = (Dlambda_w * (lambda_w + lambda_n) - lambda_w * (Dlambda_w + Dlambda_n) ) /
                      pow(lambda_w + lambda_n, 2);

    // Compute the last column of the inverse of the inverse permeability
    const Matrix invK = dataPhysical::invK( t, x, y, z );

    // Compute the denominator of the last column of the inverse of the inverse permeability
    const Real denominator = invK(0, 0) * invK(1, 1) * invK(2, 2)
                             - invK(0, 0) * invK(1, 2) * invK(1, 2)
                             - invK(0, 1) * invK(0, 1) * invK(2, 2)
                             + 2. * invK(0, 1) * invK(0, 2) * invK(1, 2)
                             - invK(0, 2) * invK(0, 2) * invK(1, 1);

    // Compute the first component of the last column of the inverse of the inverse permeability
    const Real K02 = ( invK(0, 1) * invK(1, 2) - invK(0, 2) * invK(1, 1) ) / denominator;

    // Compute the second component of the last column of the inverse of the inverse permeability
    const Real K12 = ( invK(0, 0) * invK(1, 2) - invK(0, 1) * invK(0, 2) ) / denominator;

    // Compute the third component of the last column of the inverse of the inverse permeability
    const Real K22 = ( invK(0, 0) * invK(1, 1) - invK(0, 1) * invK(0, 1) ) / denominator;

    // Compute a common constant
    const Real lfrg = (Dlambda_w * f_n + lambda_w * Df_n) * (dataPhysical::rho_w - dataPhysical::rho_n) * dataPhysical::g;

    // First row
    Real Entry0 = u[1] * Df_n - lfrg * K02;

    // Second row
    Real Entry1 = u[2] * Df_n - lfrg * K12;

    // Third row
    Real Entry2 = u[3] * Df_n - lfrg * K22;

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
    return 0.;
}

// Initial condition
Real saturationInitialCondition( const Real& /* t */,
                                 const Real& x,
                                 const Real& y,
                                 const Real& z,
                                 const ID&   /*icomp*/ )
{
    return 0.1;
}

// Mass function
Real saturationMass( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& z,
                     const ID&   /*icomp*/ )
{
    return dataPhysical::Phi( x, y, z );
}

// Boundary condition of Dirichlet
Real saturationDirichlet1( const Real& t,
                           const Real& /* x */,
                           const Real& /* y */,
                           const Real& /* z */,
                           const ID&   /*icomp*/)
{
    if ( t < 5.)
        return 0.8;
    else
        return 0.1;
}

// Boundary condition of Dirichlet
Real saturationDirichlet2( const Real& /* t */,
                           const Real& /* x */,
                           const Real& /* y */,
                           const Real& /* z */,
                           const ID&   /*icomp*/)
{
    return 0.5;
}

// Boundary condition of Dirichlet
Real saturationDirichlet3( const Real& /* t */,
                           const Real& /* x */,
                           const Real& /* y */,
                           const Real& /* z */,
                           const ID&   /*icomp*/)
{
    return 0.01;
}

// Boundary condition of Neumann
Real saturationNeumann( const Real& /* t */,
                        const Real& x,
                        const Real& y,
                        const Real& z,
                        const ID&   icomp)
{
    return 0.;

    switch (icomp)
    {
    case 1:   //! Dx
        return 0.;
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
    return 0.;
}

}
