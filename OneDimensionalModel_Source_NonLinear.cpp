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
 *  @brief
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date ???
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#include "OneDimensionalModel_Source_NonLinear.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_Source_NonLinear::OneDimensionalModel_Source_NonLinear() :
    super    ()
{}

OneDimensionalModel_Source_NonLinear::OneDimensionalModel_Source_NonLinear( const Physics_PtrType Physics ) :
    super    ( Physics )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_Source_NonLinear::operator()( const Real& _A, const Real& _Q,
                                                  const ID& ii,   const UInt& indz) const
{
    Real B2(0.);
    Real Area0, beta0, beta1, rho, Kr;
    Real beta1plus1, AoverA0, AoverA0POWbeta1plus1;
    Real tmp;


    Real dArea0dz = 0.;
    Real dbeta0dz = 0.;
    Real dbeta1dz = 0.;

    if( ii == 1 )  // B1
    {
        return 0.;
    }
    if( ii == 2 ) // B2
    {
        Area0 = M_Physics->Data()->Area0(indz);
        beta0 = M_Physics->Data()->Beta0(indz);
        beta1 = M_Physics->Data()->Beta1(indz);
        dArea0dz = M_Physics->Data()->dArea0dz(indz);
        dbeta0dz = M_Physics->Data()->dBeta0dz(indz);
        dbeta1dz = M_Physics->Data()->dBeta1dz(indz);
        Kr    = M_Physics->Data()->FrictionKr(indz);
        rho   = M_Physics->Data()->DensityRho();

        beta1plus1 = beta1 + 1;
        AoverA0 = _A / Area0;
        AoverA0POWbeta1plus1 = std::pow( AoverA0, beta1plus1 );

        tmp = beta0 / ( rho * beta1plus1 ) * AoverA0POWbeta1plus1;
        // friction term
        B2 = Kr * _Q /_A;

        // term with the derivative of A0 with respect to z
        B2 += - tmp * beta1 * dArea0dz;

        // term with the derivative of beta0 with respect to z
        B2 += Area0 / rho * ( AoverA0POWbeta1plus1 / beta1plus1  -  AoverA0 )
            * dbeta0dz;

        // term with the derivative of beta1 with respect to z
        B2 += tmp * Area0 * ( std::log( AoverA0 ) - 1. / beta1plus1 )
            * dbeta1dz;

        return M_Physics->Data()->RobertsonCorrection() * B2;
    }
    ERROR_MSG("The flux function has only 2 components.");
    return -1.;
}

Real
OneDimensionalModel_Source_NonLinear::diff( const Real& _A, const Real& _Q,
                                            const ID& ii,   const ID& jj,
                                            const UInt& indz) const
{
    Real dB2dA;
    Real Area0, beta0, beta1, rho, Kr;
    Real AoverA0, AoverA0POWbeta1;
    Real tmp;

    Real dArea0dz = 0.;
    Real dbeta0dz = 0.;
    Real dbeta1dz = 0.;

    if( ii == 1 ) // B1 = 0 so...
    {
        if( jj == 1 || jj == 2 ) // dB2/dUj = 0
        {
            return 0.;
        }
    }
    if( ii == 2 && jj == 1 ) // dB2/dA
    {
        Area0 = M_Physics->Data()->Area0(indz);
        beta0 = M_Physics->Data()->Beta0(indz);
        beta1 = M_Physics->Data()->Beta1(indz);
        dArea0dz = M_Physics->Data()->dArea0dz(indz);
        dbeta0dz = M_Physics->Data()->dBeta0dz(indz);
        dbeta1dz = M_Physics->Data()->dBeta1dz(indz);
        Kr    = M_Physics->Data()->FrictionKr(indz);
        rho   = M_Physics->Data()->DensityRho();

        AoverA0 = _A / Area0;
        AoverA0POWbeta1 = std::pow( AoverA0, beta1 );

        tmp = beta0 / rho * AoverA0POWbeta1;

        // friction term
        dB2dA = - Kr * _Q / ( _A * _A );

        // term with the derivative of A0 with respect to z
        dB2dA += - tmp * beta1 / Area0 * dArea0dz;

        // term with the derivative of beta0 with respect to z
        dB2dA += ( AoverA0POWbeta1 - 1. ) / rho * dbeta0dz;

        // term with the derivative of beta1 with respect to z
        dB2dA += tmp * std::log( AoverA0 ) * dbeta1dz;

        return M_Physics->Data()->RobertsonCorrection() * dB2dA;
    }
    if( ii == 2 && jj == 2 )  // dB2/dQ
    {
        Kr    = M_Physics->Data()->FrictionKr(indz);
        return M_Physics->Data()->RobertsonCorrection() * Kr / _A;
    }

    ERROR_MSG("Source's differential function has only 4 components.");
    return -1.;
}

// Second derivative tensor d2Bi/(dxj dxk)
Real
OneDimensionalModel_Source_NonLinear::diff2( const Real& _A, const Real& _Q,
                                             const ID& ii,   const ID& jj, const ID& kk,
                                             const UInt& indz ) const
{
    Real d2B2dA2;
    Real Area0, beta0, beta1, rho, Kr;
    Real AoverA0, AoverA0POWbeta1divA;
    Real tmp;

    Real dArea0dz = 0.;
    Real dbeta0dz = 0.;
    Real dbeta1dz = 0.;

    // B1 = 0 so ...
    if( ii == 1 ) // d2B1/dUjdUk = 0.
    {
        if( ( jj == 1 || jj == 2 ) && ( kk == 1 || kk == 2 ) )
        {
            return 0.;
        }
    }
    if( ii == 2 )
    {
        if( jj == 1 && kk == 1 ) // d2B2/dA2
        {
            // this term is not strictly necessary as it is always multiplied by 0.
            // but for the sake of generality...

            Area0 = M_Physics->Data()->Area0(indz);
            beta0 = M_Physics->Data()->Beta0(indz);
            beta1 = M_Physics->Data()->Beta1(indz);
            dArea0dz = M_Physics->Data()->dArea0dz(indz);
            dbeta0dz = M_Physics->Data()->dBeta0dz(indz);
            dbeta1dz = M_Physics->Data()->dBeta1dz(indz);
            Kr    = M_Physics->Data()->FrictionKr(indz);
            rho   = M_Physics->Data()->DensityRho();

            AoverA0 = _A / Area0;
            AoverA0POWbeta1divA = std::pow( AoverA0, beta1 ) / _A;

            tmp = beta0 / rho * AoverA0POWbeta1divA;

            // friction term
            d2B2dA2 = 2 * Kr * _Q / ( _A * _A * _A);

            // term with the derivative of A0 with respect to z
            d2B2dA2 += - tmp * beta1 * beta1 / Area0 * dArea0dz;

            // term with the derivative of beta0 with respect to z
            d2B2dA2 += beta1 / rho * AoverA0POWbeta1divA * dbeta0dz;

            // term with the derivative of beta1 with respect to z
            d2B2dA2 += tmp * ( beta1 * std::log( AoverA0 ) + 1. ) * dbeta1dz;

            return M_Physics->Data()->RobertsonCorrection() * d2B2dA2;
        }
        // cross terms (equal)
        if( (jj == 1 && kk == 2) || (jj == 2 && kk == 1) ) // d2B2/dAdQ=d2B2/dQdA
        {
            Kr    = M_Physics->Data()->FrictionKr(indz);
            return - M_Physics->Data()->RobertsonCorrection() * Kr / ( _A * _A );
        }

        if( jj == 2 && kk == 2 ) // d2B2/dQ2
        {
            return 0.;
        }
    }

    ERROR_MSG("Source's second differential function has only 8 components.");
    return -1.;
}

Real
OneDimensionalModel_Source_NonLinear::QuasiLinearSource( const Real& _A, const Real& _Q,
                                                         const ID& ii,   const UInt& indz ) const
{
    Real Sql2;
    Real Area0, beta0, beta1, rho, Kr;
    Real AoverA0, AoverA0POWbeta1timesA;
    Real tmp;

    Real dArea0dz = 0.;
    Real dbeta0dz = 0.;
    Real dbeta1dz = 0.;
    Real dalphadz = 0.;

    if( ii == 1 ) // Sql1
    {
        return 0.;
    }
    if( ii == 2 ) // Sql2
    {

        Area0 = M_Physics->Data()->Area0(indz);
        beta0 = M_Physics->Data()->Beta0(indz);
        beta1 = M_Physics->Data()->Beta1(indz);
        dArea0dz = M_Physics->Data()->dArea0dz(indz);
        dbeta0dz = M_Physics->Data()->dBeta0dz(indz);
        dbeta1dz = M_Physics->Data()->dBeta1dz(indz);
        dalphadz = M_Physics->Data()->dAlphadz(indz);
        Kr    = M_Physics->Data()->FrictionKr(indz);
        rho   = M_Physics->Data()->DensityRho();

        AoverA0               = _A/Area0;
        AoverA0POWbeta1timesA = std::pow( AoverA0, beta1 ) * _A;

        tmp  = beta0 / rho * AoverA0POWbeta1timesA;
        // friction term
        Sql2 = Kr * _Q/_A;
        // term with the derivative of A0 with respect to z
        Sql2 += - tmp * beta1 / Area0 * dArea0dz;
        // term with the derivative of beta0 with respect to z
        Sql2 += 1 / rho * ( AoverA0POWbeta1timesA  -  _A ) * dbeta0dz;
        // term with the derivative of beta1 with respect to z
        Sql2 += tmp * std::log( AoverA0 ) * dbeta1dz;
        // term with the derivative of alpha with respect to z
        Sql2 += _Q * _Q / _A * dalphadz;

        return M_Physics->Data()->RobertsonCorrection() * Sql2;
    }
    ERROR_MSG("The QL source function has only 2 components.");
    return -1.;
}

}
