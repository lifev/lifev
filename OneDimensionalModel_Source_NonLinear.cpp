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

namespace LifeV
{

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
OneDimensionalModel_Source_NonLinear::operator()( const Real& A, const Real& Q,
                                                  const ID& ii,  const UInt& i ) const
{
    if ( ii == 1 ) // B1
    {
        return 0.;
    }
    if ( ii == 2 ) // B2
    {
        Real beta1plus1( M_Physics->Data()->Beta1(i) + 1 );
        Real AoverA0( A / M_Physics->Data()->Area0(i) );
        Real C0( 1 / ( M_Physics->Data()->DensityRho() * beta1plus1 ) );
        Real C ( 1 / ( M_Physics->Data()->DensityRho() * beta1plus1 ) * std::pow(  AoverA0, beta1plus1 ) );

        return ( M_Physics->Data()->Friction() * Q / A
                 + Q * Q / A * M_Physics->Data()->dAlphadz(i)
                 + C  * ( M_Physics->Data()->Area0(i) * M_Physics->Data()->dBeta0dz(i)
                          - M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i) * M_Physics->Data()->dArea0dz(i)
                          + M_Physics->Data()->Area0(i) * M_Physics->Data()->Beta0(i) * ( std::log( AoverA0 ) - 1. / beta1plus1 ) * M_Physics->Data()->dBeta1dz(i)
                        )
                 - C0 * ( M_Physics->Data()->Area0(i) * M_Physics->Data()->dBeta0dz(i)
                          - M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i) * M_Physics->Data()->dArea0dz(i)
                          - M_Physics->Data()->Area0(i) * M_Physics->Data()->Beta0(i) / beta1plus1 * M_Physics->Data()->dBeta1dz(i)
                        )
                 + ( M_Physics->Data()->Area0(i) - A ) / M_Physics->Data()->DensityRho() * M_Physics->Data()->dBeta0dz(i)
               ) * M_Physics->Data()->RobertsonCorrection();
    }

    ERROR_MSG("The flux function has only 2 components.");

    return -1.;
}

Real
OneDimensionalModel_Source_NonLinear::diff( const Real& A, const Real& Q,
                                            const ID& ii,  const ID& jj,
                                            const UInt& i) const
{
    if ( ii == 1 ) // B1
    {
        if ( jj == 1 || jj == 2 ) // dB2/dUj = 0
        {
            return 0.;
        }
    }
    if ( ii == 2 ) // B2
    {
        if ( jj == 1 ) // dB2/dA
        {
            Real AoverA0( A / M_Physics->Data()->Area0(i) );
            Real C ( std::pow(  AoverA0, M_Physics->Data()->Beta1(i) ) / M_Physics->Data()->DensityRho() );

            return ( -M_Physics->Data()->Friction() * Q / A / A
                     - Q * Q / ( A * A ) * M_Physics->Data()->dAlphadz(i)
                     + C  * ( M_Physics->Data()->dBeta0dz(i)
                              - M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i) / M_Physics->Data()->Area0(i) * M_Physics->Data()->dArea0dz(i)
                              + M_Physics->Data()->Beta0(i) * std::log( AoverA0 ) * M_Physics->Data()->dBeta1dz(i)
                            )
                     - 1. / M_Physics->Data()->DensityRho() * M_Physics->Data()->dBeta0dz(i)
                   ) * M_Physics->Data()->RobertsonCorrection();
        }
        if ( jj == 2 ) // dB2/dQ
        {
            return M_Physics->Data()->RobertsonCorrection() * ( M_Physics->Data()->Friction() / A + 2 * Q / A * M_Physics->Data()->dAlphadz(i) );
        }
    }

    ERROR_MSG("Source's differential function has only 4 components.");

    return -1.;
}

// Second derivative tensor d2Bi/(dxj dxk)
//Real
//OneDimensionalModel_Source_NonLinear::diff2( const Real& A, const Real& Q,
//                                             const ID& ii,   const ID& jj, const ID& kk,
//                                             const UInt& i ) const
//{
//    Real d2B2dA2;
//    Real Area0, beta0, beta1, rho, Kr;
//    Real AoverA0, AoverA0POWbeta1divA;
//    Real tmp;
//
//    Real dArea0dz = 0.;
//    Real dbeta0dz = 0.;
//    Real dbeta1dz = 0.;
//
//    // B1 = 0 so ...
//    if( ii == 1 ) // d2B1/dUjdUk = 0.
//    {
//        if( ( jj == 1 || jj == 2 ) && ( kk == 1 || kk == 2 ) )
//        {
//            return 0.;
//        }
//    }
//    if( ii == 2 )
//    {
//        if( jj == 1 && kk == 1 ) // d2B2/dA2
//        {
//            // this term is not strictly necessary as it is always multiplied by 0.
//            // but for the sake of generality...
//
//            Area0 = M_Physics->Data()->Area0(i);
//            beta0 = M_Physics->Data()->Beta0(i);
//            beta1 = M_Physics->Data()->Beta1(i);
//            dArea0dz = M_Physics->Data()->dArea0dz(i);
//            dbeta0dz = M_Physics->Data()->dBeta0dz(i);
//            dbeta1dz = M_Physics->Data()->dBeta1dz(i);
//            Kr    = M_Physics->Data()->Friction();
//            rho   = M_Physics->Data()->DensityRho();
//
//            AoverA0 = A / Area0;
//            AoverA0POWbeta1divA = std::pow( AoverA0, beta1 ) / A;
//
//            tmp = beta0 / rho * AoverA0POWbeta1divA;
//
//            // friction term
//            d2B2dA2 = 2 * Kr * Q / ( A * A * A);
//
//            // term with the derivative of A0 with respect to z
//            d2B2dA2 += - tmp * beta1 * beta1 / Area0 * dArea0dz;
//
//            // term with the derivative of beta0 with respect to z
//            d2B2dA2 += beta1 / rho * AoverA0POWbeta1divA * dbeta0dz;
//
//            // term with the derivative of beta1 with respect to z
//            d2B2dA2 += tmp * ( beta1 * std::log( AoverA0 ) + 1. ) * dbeta1dz;
//
//            return M_Physics->Data()->RobertsonCorrection() * d2B2dA2;
//        }
//        // cross terms (equal)
//        if( (jj == 1 && kk == 2) || (jj == 2 && kk == 1) ) // d2B2/dAdQ=d2B2/dQdA
//        {
//            Kr    = M_Physics->Data()->Friction();
//            return - M_Physics->Data()->RobertsonCorrection() * Kr / ( A * A );
//        }
//
//        if( jj == 2 && kk == 2 ) // d2B2/dQ2
//        {
//            return 0.;
//        }
//    }
//
//    ERROR_MSG("Source's second differential function has only 8 components.");
//    return -1.;
//}

Real
OneDimensionalModel_Source_NonLinear::interpolatedQuasiLinearSource( const Real& A, const Real& Q,
                                                                     const ID& ii, const Container2D_Type& bcNodes, const Real& cfl ) const
{
    if ( ii == 1 ) // QLS1
    {
        return 0.;
    }
    if ( ii == 2 ) // QLS2
    {
        // Interpolate quantities
        Real Area0      = ( 1 - cfl ) * M_Physics->Data()->Area0(bcNodes[0])    + cfl * M_Physics->Data()->Area0(bcNodes[1]);
        Real Beta0      = ( 1 - cfl ) * M_Physics->Data()->Beta0(bcNodes[0])    + cfl * M_Physics->Data()->Beta0(bcNodes[1]);
        Real Beta1      = ( 1 - cfl ) * M_Physics->Data()->Beta1(bcNodes[0])    + cfl * M_Physics->Data()->Beta1(bcNodes[1]);
        Real dArea0dz   = ( 1 - cfl ) * M_Physics->Data()->dArea0dz(bcNodes[0]) + cfl * M_Physics->Data()->dArea0dz(bcNodes[1]);
        Real dBeta0dz   = ( 1 - cfl ) * M_Physics->Data()->dBeta0dz(bcNodes[0]) + cfl * M_Physics->Data()->dBeta0dz(bcNodes[1]);
        Real dBeta1dz   = ( 1 - cfl ) * M_Physics->Data()->dBeta1dz(bcNodes[0]) + cfl * M_Physics->Data()->dBeta1dz(bcNodes[1]);
        Real dAlphadz   = ( 1 - cfl ) * M_Physics->Data()->dAlphadz(bcNodes[0]) + cfl * M_Physics->Data()->dAlphadz(bcNodes[1]);

        Real AoverA0( A / Area0 );
        Real C( A / M_Physics->Data()->DensityRho() * std::pow(  AoverA0, Beta1 ) );

        return ( M_Physics->Data()->Friction() * Q / A
                 + Q * Q / A * dAlphadz
                 + C * ( dBeta0dz
                         - Beta0 * Beta1 / Area0 * dArea0dz
                         + Beta0 * std::log( AoverA0 ) * dBeta1dz
                       )
                 - A / M_Physics->Data()->DensityRho() * dBeta0dz
               ) * M_Physics->Data()->RobertsonCorrection();
    }

    ERROR_MSG("The QLS function has only 2 components.");

    return -1.;
}

}
