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

#include <life/lifesolver/OneDimensionalPhysicsNonLinear.hpp>

namespace LifeV
{
// ===================================================
// Conversion methods
// ===================================================
void
OneDimensionalPhysicsNonLinear::fromUToW( Real& W1, Real& W2, const Real& A,  const Real& Q, const UInt& iNode ) const
{
    Real celerity( celerity0( iNode ) * std::sqrt( OneDimensional::pow05( A / M_data -> area0( iNode ), M_data -> beta1( iNode ) ) ) );

    Real add( std::sqrt( M_data -> robertsonCorrection() ) * ( celerity - celerity0( iNode ) ) * 2 / M_data -> beta1( iNode ) );

    Real QoverA  = Q / A;

    W1 = QoverA + add;
    W2 = QoverA - add;
}

void
OneDimensionalPhysicsNonLinear::fromWToU( Real& A, Real& Q, const Real& W1, const Real& W2, const UInt& iNode ) const
{
    Real rhooverbeta0beta1 ( M_data -> densityRho() / ( M_data -> beta0( iNode ) * M_data -> beta1( iNode ) ) );

    Real beta1over4SQRTchi( M_data -> beta1( iNode ) / ( std::sqrt(M_data -> robertsonCorrection() ) * 4 ) );

    A = M_data -> area0( iNode )
        * OneDimensional::pow20( rhooverbeta0beta1, 1 / M_data -> beta1( iNode ) )
        * OneDimensional::pow40( beta1over4SQRTchi * (W1 - W2) + celerity0( iNode ), 2 / M_data -> beta1( iNode ) );

    Q = A * ( W1 + W2 ) / 2;
}

Real
OneDimensionalPhysicsNonLinear::fromWToP( const Real& W1, const Real& W2, const UInt& iNode ) const
{
    Real rhooverbeta0beta1 ( M_data -> densityRho() / ( M_data -> beta0( iNode ) * M_data -> beta1( iNode ) ) );

    Real beta1over4SQRTchi( M_data -> beta1( iNode ) / ( std::sqrt(M_data -> robertsonCorrection()) * 4 ) );

    return M_data -> beta0( iNode ) * ( rhooverbeta0beta1 * ( beta1over4SQRTchi * (W1 - W2) + celerity0( iNode ) ) * ( beta1over4SQRTchi * (W1 - W2) + celerity0( iNode ) ) - 1 );
}

Real
OneDimensionalPhysicsNonLinear::fromPToW( const Real& P, const Real& W, const ID& iW, const UInt& iNode ) const
{
    Real SQRTbeta0beta1overrho( M_data -> beta0( iNode ) * M_data -> beta1( iNode ) / M_data -> densityRho() );
    SQRTbeta0beta1overrho = std::sqrt( SQRTbeta0beta1overrho );

    Real SQRTchi4overbeta1( std::sqrt(M_data -> robertsonCorrection()) * 4 / M_data -> beta1( iNode ) );

    Real add( SQRTchi4overbeta1 * SQRTbeta0beta1overrho
              * ( std::sqrt( P / M_data -> beta0( iNode ) + 1 ) - 1 ) );

#ifdef HAVE_LIFEV_DEBUG
    Debug(6320) << "[OneDimensionalModel_Physics_NonLinear::W_fromP] "
    << "SQRTchi4overbeta1 = " << SQRTchi4overbeta1
    << ", beta0beta1overrho = " << SQRTbeta0beta1overrho
    << ", pow( ( P / M_data -> beta0( iNode ) + 1 ), 0.5 ) = " << std::sqrt( ( P / M_data -> beta0( iNode ) + 1 ) ) << "\n";
    Debug(6320) << "[OneDimensionalModel_Physics_NonLinear::W_fromP] add term = " << add << "\n";
#endif

    if ( iW == 0 )
        return W - add;
    if ( iW == 1 )
        return W + add;

    ERROR_MSG("You can only find W1 or W2 as function of P");
    return -1.;
}

Real
OneDimensionalPhysicsNonLinear::fromQToW( const Real& Q, const Real& W_tn, const Real& W, const ID& iW, const UInt& iNode ) const
{
    Real K0( M_data -> beta1( iNode ) / ( std::sqrt(M_data -> robertsonCorrection()) * 4 ) );

    Real K1( (M_data -> area0( iNode ) / 2) );
    K1 *= OneDimensional::pow20( M_data -> densityRho() / (M_data -> beta0( iNode ) * M_data -> beta1( iNode )), 1 / M_data -> beta1( iNode ) );
    K1 *= OneDimensional::pow40( K0, 2/M_data -> beta1( iNode ) );

    Real f_k, df_k, tau_k(0);

    if ( iW == 0 ) // W1 given
    {
        f_k = OneDimensional::pow40( W - W_tn + celerity0( iNode ) / K0, 2/M_data -> beta1( iNode ) );
        tau_k = OneDimensional::pow40( W - W_tn + celerity0( iNode ) / K0, 2/M_data -> beta1( iNode ) );
        df_k = (-2 / M_data -> beta1( iNode )) * OneDimensional::pow30( W - W_tn + celerity0( iNode ) / K0, 2/M_data -> beta1( iNode ) - 1 );
    }
    if ( iW == 1 ) // W2 given
    {
        f_k = OneDimensional::pow40( W_tn - W + celerity0( iNode ) / K0, 2/M_data -> beta1( iNode ) );
        tau_k = OneDimensional::pow40( W_tn - W + celerity0( iNode ) / K0, 2/M_data -> beta1( iNode ) );
        df_k = (-2 / M_data -> beta1( iNode )) * OneDimensional::pow30( W_tn - W + celerity0( iNode ) / K0, 2/M_data -> beta1( iNode ) - 1 );
    }
    f_k *= (W + W_tn);
    f_k += - Q / K1;
    df_k *= (W + W_tn);
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
OneDimensionalPhysicsNonLinear::dPdW( const Real& W1, const Real& W2, const ID& iW, const UInt& iNode ) const
{
    Real rhoover2SQRTchi ( M_data -> densityRho() / ( std::sqrt(M_data -> robertsonCorrection()) * 2 ) );

    Real beta1over4SQRTchi( M_data -> beta1( iNode ) / ( std::sqrt(M_data -> robertsonCorrection()) * 4 ) );

    Real result( beta1over4SQRTchi * (W1 - W2)  );
    result += celerity0( iNode );
    result *= rhoover2SQRTchi;

    if ( iW == 0 ) //! dP/dW1
    {
        return result;
    }

    if ( iW == 1 ) //! dP/dW2
    {
        return -result;
    }

    ERROR_MSG("P(W1,W2)'s differential function has only 2 components.");
    return -1.;
}

}
