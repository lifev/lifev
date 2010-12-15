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
 *  @file
 *  @brief File containing a class providing non linear physical operations for the 1D model data.
 *
 *  @version 1.0
 *  @date 01-07-2004
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 13-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/OneDimensionalModel_Physics_NonLinear.hpp>

namespace LifeV
{
// ===================================================
// Conversion methods
// ===================================================
void
OneDimensionalModel_Physics_NonLinear::fromUToW( Real& W1, Real& W2, const Real& A,  const Real& Q, const UInt& indz ) const
{
    Real celerity( celerity0(indz) * std::sqrt( std::pow( A / M_data -> area0(indz), M_data -> beta1(indz) ) ) );

    Real add( std::sqrt( M_data -> robertsonCorrection() ) * ( celerity - celerity0(indz) ) * 2 / M_data -> beta1(indz) );

    Real QoverA  = Q / A;

    W1 = QoverA + add;
    W2 = QoverA - add;
}

void
OneDimensionalModel_Physics_NonLinear::fromWToU( Real& A, Real& Q, const Real& W1, const Real& W2, const UInt& indz ) const
{
    Real rhooverbeta0beta1 ( M_data -> densityRho() / ( M_data -> beta0(indz) * M_data -> beta1(indz) ) );

    Real beta1over4SQRTchi( M_data -> beta1(indz) / ( std::sqrt(M_data -> robertsonCorrection() ) * 4 ) );

    A = M_data -> area0(indz)
        * std::pow( rhooverbeta0beta1, (1/M_data -> beta1(indz)) )
        * std::pow( beta1over4SQRTchi * (W1 - W2) + celerity0(indz), (2/M_data -> beta1(indz)) );

    Q = A * ( W1 + W2 ) / 2;
}

Real
OneDimensionalModel_Physics_NonLinear::fromWToP( const Real& W1, const Real& W2, const UInt& indz ) const
{
    Real rhooverbeta0beta1 ( M_data -> densityRho() / ( M_data -> beta0(indz) * M_data -> beta1(indz) ) );

    Real beta1over4SQRTchi( M_data -> beta1(indz) / ( std::sqrt(M_data -> robertsonCorrection()) * 4 ) );

    return M_data -> beta0(indz) * ( rhooverbeta0beta1 * std::pow( beta1over4SQRTchi * (W1 - W2) + celerity0(indz), 2 ) - 1 );
}

Real
OneDimensionalModel_Physics_NonLinear::fromPToW( const Real& P, const Real& W, const ID& i, const UInt& indz ) const
{
    Real SQRTbeta0beta1overrho( M_data -> beta0(indz) * M_data -> beta1(indz) / M_data -> densityRho() );
    SQRTbeta0beta1overrho = std::sqrt( SQRTbeta0beta1overrho );

    Real SQRTchi4overbeta1( std::sqrt(M_data -> robertsonCorrection()) * 4 / M_data -> beta1(indz) );

    Real add( SQRTchi4overbeta1 * SQRTbeta0beta1overrho
              * ( pow( ( P / M_data -> beta0(indz) + 1 ), 0.5 ) - 1 ) );

#ifdef HAVE_LIFEV_DEBUG
    Debug(6320) << "[OneDimensionalModel_Physics_NonLinear::W_fromP] "
    << "SQRTchi4overbeta1 = " << SQRTchi4overbeta1
    << ", beta0beta1overrho = " << SQRTbeta0beta1overrho
    << ", pow( ( P / M_data -> beta0(indz) + 1 ), 0.5 ) = " << pow( ( P / M_data -> beta0(indz) + 1 ), 0.5 ) << "\n";
    Debug(6320) << "[OneDimensionalModel_Physics_NonLinear::W_fromP] add term = " << add << "\n";
#endif

    if ( i == 1 )
        return W - add;
    if ( i == 2 )
        return W + add;

    ERROR_MSG("You can only find W1 or W2 as function of P");
    return -1.;
}

Real
OneDimensionalModel_Physics_NonLinear::fromQToW( const Real& Q, const Real& W_n, const Real& W, const ID& i, const UInt& indz ) const
{
    Real K0( M_data -> beta1(indz) / ( std::sqrt(M_data -> robertsonCorrection()) * 4 ) );

    Real K1( (M_data -> area0(indz) / 2) );
    K1 *= pow( M_data -> densityRho() / (M_data -> beta0(indz) * M_data -> beta1(indz)), 1/M_data -> beta1(indz) );
    K1 *= pow( K0, 2/M_data -> beta1(indz) );

    Real w_k = W_n;
    Real f_k, df_k, tau_k(0);

    if ( i == 1 ) // W1 given
    {
        f_k = pow( W - w_k + celerity0(indz) / K0, 2/M_data -> beta1(indz) );
        tau_k = pow( W - w_k + celerity0(indz) / K0, 2/M_data -> beta1(indz) );
        df_k = (-2 / M_data -> beta1(indz)) * pow( W - w_k + celerity0(indz) / K0, 2/M_data -> beta1(indz) - 1 );
    }
    if ( i == 2 ) // W2 given
    {
        f_k = pow( w_k - W + celerity0(indz) / K0, 2/M_data -> beta1(indz) );
        tau_k = pow( w_k - W + celerity0(indz) / K0, 2/M_data -> beta1(indz) );
        df_k = (-2 / M_data -> beta1(indz)) * pow( w_k - W + celerity0(indz) / K0, 2/M_data -> beta1(indz) - 1 );
    }
    f_k *= (W + w_k);
    f_k += - Q / K1;
    df_k *= (W + w_k);
    df_k += f_k;

#ifdef HAVE_LIFEV_DEBUG
    Debug(6320) << "[OneDimensionalModel_Physics_NonLinear::W_fromQ] "
    << "K0 = " << K0
    << ", K1 = " << K1
    << ", tau_k = " << tau_k << "\n";
#endif

    Real w_kp1 = Q / (K1 * tau_k) - W;

    return w_kp1;
}

// ===================================================
// Derivatives methods
// ===================================================
Real
OneDimensionalModel_Physics_NonLinear::dPdW( const Real& W1, const Real& W2, const ID& i, const UInt& indz ) const
{
    Real rhoover2SQRTchi ( M_data -> densityRho() / ( std::sqrt(M_data -> robertsonCorrection()) * 2 ) );

    Real beta1over4SQRTchi( M_data -> beta1(indz) / ( std::sqrt(M_data -> robertsonCorrection()) * 4 ) );

    Real result( beta1over4SQRTchi * (W1 - W2)  );
    result += celerity0(indz);
    result *= rhoover2SQRTchi;

    if ( i == 1 ) //! dP/dW1
    {
        return result;
    }

    if ( i == 2 ) //! dP/dW2
    {
        return -result;
    }

    ERROR_MSG("P(W1,W2)'s differential function has only 2 components.");
    return -1.;
}

}
