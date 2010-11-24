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
OneDimensionalModel_Flux_NonLinear::operator()( const Real& A, const Real& Q,
                                                const ID& ii,  const UInt& i ) const
{
    if( ii == 1 ) // F1
    {
        return Q;
    }

    if( ii == 2 ) // F2
    {
        return ( M_Physics->Data()->Alpha(i) * Q * Q / A
               + M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i) * M_Physics->Data()->Area0(i)
               / ( ( M_Physics->Data()->Beta1(i) + 1) * M_Physics->Data()->DensityRho() )
               * ( std::pow( A / M_Physics->Data()->Area0(i), M_Physics->Data()->Beta1(i) + 1 ) - 1 )
               ) * M_Physics->Data()->RobertsonCorrection();
    }

    ERROR_MSG("The flux function has only 2 components.");

    return -1.;
}

Real
OneDimensionalModel_Flux_NonLinear::diff( const Real& A, const Real& Q,
                                          const ID& ii,   const ID& jj, const UInt& i ) const
{
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
        return ( M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i) / M_Physics->Data()->DensityRho()
               * std::pow( A / M_Physics->Data()->Area0(i), M_Physics->Data()->Beta1(i) )
               - M_Physics->Data()->Alpha(i) * Q * Q / A / A
               ) * M_Physics->Data()->RobertsonCorrection();
    }
    if( ii == 2 && jj == 2 )  // dF2/dQ
    {
        return M_Physics->Data()->RobertsonCorrection() * 2 * M_Physics->Data()->Alpha(i) * Q / A;
    }

    ERROR_MSG("Flux's differential function has only 4 components.");

    return -1.;
}

//Real
//OneDimensionalModel_Flux_NonLinear::diff2( const Real& A, const Real& Q,
//                                           const ID& ii,   const ID& jj, const ID& kk,
//                                           const UInt& i ) const
//{
//    // diff second of F1 is always 0.
//    if( ii == 1 ) // d2F1/dUjdUk = 0.
//    {
//        if( ( jj == 1 || jj == 2 ) && ( kk == 1 || kk == 2 ) )
//            return 0.;
//    }
//    if( ii == 2 )
//    {
//        if( jj == 1 && kk == 1 )  // d2F2/dA2
//        {
//            return M_Physics->Data()->RobertsonCorrection()
//                       * M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i) * M_Physics->Data()->Beta1(i)
//                       / ( M_Physics->Data()->DensityRho() * M_Physics->Data()->Area0(i) )
//                       * std::pow( A / M_Physics->Data()->Area0(i), M_Physics->Data()->Beta1(i) - 1);
//        }
//        // cross terms (equal)
//        if( (jj == 1 && kk == 2) || (jj == 2 && kk == 1) ) // d2F2/dAdQ=d2F2/dQdA
//        {
//            return -M_Physics->Data()->RobertsonCorrection() * M_Physics->Data()->Alpha(i) * Q / ( A * A );
//        }
//        if( jj == 2 && kk == 2 ) // d2F2/dQ2
//        {
//            return M_Physics->Data()->RobertsonCorrection() * 2 * M_Physics->Data()->Alpha(i) / A;
//        }
//    }
//    ERROR_MSG("Flux's second differential function has only 8 components.");
//
//    return -1.;
//}

void
OneDimensionalModel_Flux_NonLinear::EigenValuesEigenVectors( const Real& A,
                                                             const Real& Q,
                                                                   Container2D_Type& eigenvalues,
                                                                   Container2D_Type& leftEigenvector1,
                                                                   Container2D_Type& leftEigenvector2,
                                                             const UInt& i ) const
{
    Debug(6312) << "[OneDimensionalModel_Flux_NonLinear]::jabocian_EigenValues_Vectors\n";

    Real celerity;
    celerity       = std::sqrt( M_Physics->Data()->Alpha(i) * ( M_Physics->Data()->Alpha(i) - 1) * Q * Q / ( A * A )
                              + M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i)
                              / M_Physics->Data()->DensityRho() * std::pow( A / M_Physics->Data()->Area0(i), M_Physics->Data()->Beta1(i) ) );

    eigenvalues[0] = M_Physics->Data()->Alpha(i) * Q / A + celerity;
    eigenvalues[1] = M_Physics->Data()->Alpha(i) * Q / A - celerity;

    leftEigenvector1[0] = - eigenvalues[1] / A;
    leftEigenvector1[1] = 1. / A;
    leftEigenvector2[0] = - eigenvalues[0] / A;
    leftEigenvector2[1] = 1. / A;
}

void
OneDimensionalModel_Flux_NonLinear::deltaEigenValuesEigenVectors( const Real& A,
                                                                  const Real& Q,
                                                                        Container2D_Type& deltaEigenvalues,
                                                                        Container2D_Type& deltaLeftEigenvector1,
                                                                        Container2D_Type& deltaLeftEigenvector2,
                                                                  const UInt& i ) const
{
    Real deltaCelerity;

    Real AoverA0( A / M_Physics->Data()->Area0(i) );
    Real C ( std::pow(  AoverA0, M_Physics->Data()->Beta1(i) ) / M_Physics->Data()->DensityRho() );

    deltaCelerity  = 0.5 / std::sqrt( M_Physics->Data()->Alpha(i) * ( M_Physics->Data()->Alpha(i) - 1) * Q * Q / ( A * A )
                                    + M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i) * C )
                         * ( C * (  M_Physics->Data()->Beta1(i) * M_Physics->Data()->dBeta0dz(i)
                                  - M_Physics->Data()->Beta0(i) * M_Physics->Data()->Beta1(i) * M_Physics->Data()->Beta1(i) /  M_Physics->Data()->Area0(i) * M_Physics->Data()->dArea0dz(i)
                                  + M_Physics->Data()->Beta0(i) * ( 1 + M_Physics->Data()->Beta0(i) * std::log( AoverA0 ) ) * M_Physics->Data()->dBeta1dz(i) )
                             + ( 2 * M_Physics->Data()->Alpha(i) - 1 ) * Q * Q / ( A * A ) * M_Physics->Data()->dAlphadz(i) );

    deltaEigenvalues[0] = M_Physics->Data()->dAlphadz(i) * Q / A + deltaCelerity;
    deltaEigenvalues[1] = M_Physics->Data()->dAlphadz(i) * Q / A - deltaCelerity;

    deltaLeftEigenvector1[0] = - deltaEigenvalues[1] / A;
    deltaLeftEigenvector1[1] = 0;
    deltaLeftEigenvector2[0] = - deltaEigenvalues[0] / A;
    deltaLeftEigenvector2[1] = 0;
}

}
