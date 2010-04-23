//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a class for non linear 1D model flux function.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#include "OneDimensionalModel_Flux_NonLinear.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_Flux_NonLinear::OneDimensionalModel_Flux_NonLinear() :
    super   ()
{}

OneDimensionalModel_Flux_NonLinear::OneDimensionalModel_Flux_NonLinear( const Physics_PtrType Physics ) :
    super   ( Physics )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_Flux_NonLinear::operator()( const Real& _A, const Real& _Q,
                                                const ID& ii, const UInt& indz ) const
{
    Real F2;
    Real Area0, alphaCor, beta0, beta1, rho;
    Real AoverA0, AoverA0POWbeta1;

    if( ii == 1 ) // F1
    {
        return _Q;
    }

    if( ii == 2 ) // F2
    {

        Area0 = M_Physics->Data()->Area0(indz);
        alphaCor = M_Physics->Data()->RobertsonCorrection() * M_Physics->Data()->AlphaCor(indz);
        beta0 = M_Physics->Data()->Beta0(indz);
        beta1 = M_Physics->Data()->Beta1(indz);
        rho   = M_Physics->Data()->DensityRho();

        AoverA0 = _A / Area0;
        AoverA0POWbeta1 = std::pow( AoverA0, beta1 );

        F2 = alphaCor * _Q * _Q / _A ;
        F2 += M_Physics->Data()->RobertsonCorrection() * beta0 * beta1 / ( rho * ( beta1 + 1) ) * _A * AoverA0POWbeta1;

        // added Robertson, Zakaria correction factor 3/4 (TP 05/07)
        //        return _M_Robertson_correction * F2;
        return F2;
    }
    ERROR_MSG("The flux function has only 2 components.");
    return -1.;
}

Real
OneDimensionalModel_Flux_NonLinear::diff( const Real& _A, const Real& _Q,
                                          const ID& ii,   const ID& jj, const UInt& indz ) const
{
    Real dF2dA;
    Real Area0, alphaCor, beta0, beta1, rho;
    Real AoverA0, AoverA0POWbeta1;

    if( ii == 1 && jj == 1 ) // dF1/dA
    {
        return 0.;
    }

    if( ii == 1 && jj == 2 )  // dF1/dQ
    {
        return 1.;
    }

    if( ii == 2 && jj == 1 )  // dF2/dA
    {
        Area0 = M_Physics->Data()->Area0(indz);
        alphaCor = M_Physics->Data()->RobertsonCorrection() * M_Physics->Data()->AlphaCor(indz);
        beta0 = M_Physics->Data()->Beta0(indz);
        beta1 = M_Physics->Data()->Beta1(indz);
        rho   = M_Physics->Data()->DensityRho();

        AoverA0 = _A / Area0;
        AoverA0POWbeta1 = std::pow( AoverA0, beta1 );

        dF2dA = - alphaCor * _Q * _Q / ( _A * _A );
        dF2dA += M_Physics->Data()->RobertsonCorrection() * beta0 * beta1 / rho * AoverA0POWbeta1;

        // added Robertson, Zakaria correction factor 3/4 (TP 05/07)
        return dF2dA;
    }
    if( ii == 2 && jj == 2 )  // dF2/dQ
    {
        alphaCor = M_Physics->Data()->RobertsonCorrection() * M_Physics->Data()->AlphaCor(indz);
        // added Robertson, Zakaria correction factor 3/4 (TP 05/07)
        return 2 * alphaCor * _Q / _A;
    }

    ERROR_MSG("Flux's differential function has only 4 components.");
    return -1.;
}

void
OneDimensionalModel_Flux_NonLinear::jacobian_EigenValues_Vectors( const Real& _A,
                                                                  const Real& _Q,
                                                                        Real& eig1,
                                                                        Real& eig2,
                                                                        Real& lefteigvec11,
                                                                        Real& lefteigvec12,
                                                                        Real& lefteigvec21,
                                                                        Real& lefteigvec22,
                                                                  const UInt& indz ) const
{
    Debug(6312) << "[OneDimensionalModel_Flux_NonLinear]::jabocian_EigenValues_Vectors\n";

    Real Area0, alphaCor, beta0, beta1, rho;
    Real AoverA0, AoverA0POWbeta1, QoverA;

    Real celeralpha;

    Area0           = M_Physics->Data()->Area0(indz - 1);
    alphaCor        = M_Physics->Data()->RobertsonCorrection() * M_Physics->Data()->AlphaCor(indz - 1);
    beta0           = M_Physics->Data()->Beta0(indz - 1);
    beta1           = M_Physics->Data()->Beta1(indz - 1);
    rho             = M_Physics->Data()->DensityRho();

    AoverA0         = _A / Area0;
    AoverA0POWbeta1 = std::pow( AoverA0, beta1 );

    QoverA          = _Q / _A;

    celeralpha      = alphaCor * ( alphaCor - 1) * QoverA * QoverA;
    celeralpha     += M_Physics->Data()->RobertsonCorrection() * beta0 * beta1 / rho * AoverA0POWbeta1;
    celeralpha      = std::sqrt( celeralpha );

    /*
      std::cout << "\nArea in compute eigenvalues/vectors " << _A
      << "\nFlux in compute eigenvalues/vectors " << _Q
      << "\nArea0 in compute eigenvalues/vectors " << Area0
      << "\nalphaCor in compute eigenvalues/vectors " << alphaCor
      << "\nbeta0 in compute eigenvalues/vectors " << beta0
      << "\nbeta1 in compute eigenvalues/vectors " << beta1
      << "\nrho in compute eigenvalues/vectors " << rho
      << "\nAoverA0 in compute eigenvalues/vectors " << AoverA0
      << "\nAoverA0POWbeta1 in compute eigenvalues/vectors " << AoverA0POWbeta1
      << "\nQoverA in compute eigenvalues/vectors " << QoverA
      << "\nceleralpha in compute eigenvalues/vectors " << celeralpha
      << "\nRobertson_correction in compute eigenvalues/vectors " << M_Physics->Data()->RobertsonCorrection()
      << std::endl;
    */

    // eigen values
    eig1 =   celeralpha + alphaCor * QoverA;
    eig2 = - celeralpha + alphaCor * QoverA;
    // eigen vectors
    lefteigvec11 = - eig2 / _A;
    lefteigvec12 = 1. / _A;
    lefteigvec21 = - eig1 / _A;
    lefteigvec22 = 1. / _A;
}

Real
OneDimensionalModel_Flux_NonLinear::diff2( const Real& _A, const Real& _Q,
                                           const ID& ii,   const ID& jj, const ID& kk,
                                           const UInt& indz ) const
{
    Real d2F2dA2;
    Real Area0, alphaCor, beta0, beta1, rho;
    Real AoverA0, AoverA0POWbeta1divA;

    // diff second of F1 is always 0.
    if( ii == 1 ) // d2F1/dUjdUk = 0.
    {
        if( ( jj == 1 || jj == 2 ) && ( kk == 1 || kk == 2 ) )
            return 0.;
    }
    if( ii == 2 )
    {
        if( jj == 1 && kk == 1 )  // d2F2/dA2
        {
            Area0 = M_Physics->Data()->Area0(indz);
            alphaCor = M_Physics->Data()->RobertsonCorrection() * M_Physics->Data()->AlphaCor(indz);
            beta0 = M_Physics->Data()->Beta0(indz);
            beta1 = M_Physics->Data()->Beta1(indz);
            rho   = M_Physics->Data()->DensityRho();

            AoverA0 = _A / Area0;
            AoverA0POWbeta1divA = std::pow( AoverA0, beta1 ) / _A;

            d2F2dA2 = 2 * alphaCor * _Q * _Q / ( _A * _A * _A );
            d2F2dA2 += M_Physics->Data()->RobertsonCorrection() * beta0 * beta1 * beta1 / rho * AoverA0POWbeta1divA;

            return d2F2dA2;
        }
        // cross terms (equal)
        if( (jj == 1 && kk == 2) || (jj == 2 && kk == 1) ) // d2F2/dAdQ=d2F2/dQdA
        {
            alphaCor = M_Physics->Data()->RobertsonCorrection() * M_Physics->Data()->AlphaCor(indz);
            return - 2 * alphaCor * _Q / ( _A * _A );
        }
        if( jj == 2 && kk == 2 ) // d2F2/dQ2
        {
            alphaCor = M_Physics->Data()->RobertsonCorrection() * M_Physics->Data()->AlphaCor(indz);
            return 2 * alphaCor / _A;
        }
    }
    ERROR_MSG("Flux's second differential function has only 8 components.");
    return -1.;
}

}
