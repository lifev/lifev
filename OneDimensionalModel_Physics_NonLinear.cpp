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
    @brief File containing a class providing non linear physical operations for the 1D model data.

    @version 1.0
    @date 01-07-2004
    @author Vincent Martin

    @version 2.0
    @date 13-04-2010
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>

    @contributor Simone Rossi <simone.rossi@epfl.ch>

    @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */


#include <lifemc/lifesolver/OneDimensionalModel_Physics_NonLinear.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
OneDimensionalModel_Physics_NonLinear::OneDimensionalModel_Physics_NonLinear() :
    super   ()
{
}

OneDimensionalModel_Physics_NonLinear::OneDimensionalModel_Physics_NonLinear( const dataPtr_Type Data ) :
    super   ( Data )
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_Physics_NonLinear::W_from_U(       Real& W1,       Real& W2,
                                                       const Real& A,  const Real& Q, const UInt& indz ) const
{
    Real celerity( Celerity0(indz) * std::sqrt( std::pow( A / M_data -> Area0(indz), M_data -> Beta1(indz) ) ) );

    Real add( std::sqrt( M_data -> RobertsonCorrection() ) * ( celerity - Celerity0(indz) ) * 2 / M_data -> Beta1(indz) );

    Real QoverA  = Q / A;

    W1 = QoverA + add;
    W2 = QoverA - add;
}

void
OneDimensionalModel_Physics_NonLinear::U_from_W(       Real& A,        Real& Q,
                                                       const Real& W1, const Real& W2, const UInt& indz ) const
{
    Real rhooverbeta0beta1 ( M_data -> DensityRho() / ( M_data -> Beta0(indz) * M_data -> Beta1(indz) ) );

    Real beta1over4SQRTchi( M_data -> Beta1(indz) / ( std::sqrt(M_data -> RobertsonCorrection()) * 4 ) );

    A = M_data -> Area0(indz)
        * std::pow( rhooverbeta0beta1, (1/M_data -> Beta1(indz)) )
        * std::pow( beta1over4SQRTchi * (W1 - W2) + Celerity0(indz), (2/M_data -> Beta1(indz)) );

    Q = A * ( W1 + W2 ) / 2;
}

Real
OneDimensionalModel_Physics_NonLinear::pressure_W( const Real& W1, const Real& W2, const UInt& indz ) const
{
    Real rhooverbeta0beta1 ( M_data -> DensityRho() / ( M_data -> Beta0(indz) * M_data -> Beta1(indz) ) );

    Real beta1over4SQRTchi( M_data -> Beta1(indz) / ( std::sqrt(M_data -> RobertsonCorrection()) * 4 ) );

    return M_data -> Beta0(indz) * ( rhooverbeta0beta1 * std::pow( beta1over4SQRTchi * (W1 - W2) + Celerity0(indz), 2 ) - 1 );
}

Real
OneDimensionalModel_Physics_NonLinear::pressure_WDiff( const Real& W1, const Real& W2,
                                                       const ID& i,     const UInt& indz ) const
{
    Real rhoover2SQRTchi ( M_data -> DensityRho() / ( std::sqrt(M_data -> RobertsonCorrection()) * 2 ) );

    Real beta1over4SQRTchi( M_data -> Beta1(indz) / ( std::sqrt(M_data -> RobertsonCorrection()) * 4 ) );

    Real result( beta1over4SQRTchi * (W1 - W2)  );
    result += Celerity0(indz);
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

Real
OneDimensionalModel_Physics_NonLinear::W_from_P( const Real& _P, const Real& _W,
                                                 const ID& i, const UInt& indz ) const
{
    Real SQRTbeta0beta1overrho( M_data -> Beta0(indz) * M_data -> Beta1(indz) / M_data -> DensityRho() );
    SQRTbeta0beta1overrho = std::sqrt( SQRTbeta0beta1overrho );

    Real SQRTchi4overbeta1( std::sqrt(M_data -> RobertsonCorrection()) * 4 / M_data -> Beta1(indz) );

    Real add( SQRTchi4overbeta1 * SQRTbeta0beta1overrho
              * ( pow( ( _P / M_data -> Beta0(indz) + 1 ), 0.5 ) - 1 ) );

    Debug(6320) << "[OneDimensionalModel_Physics_NonLinear::W_from_P] "
    << "SQRTchi4overbeta1 = " << SQRTchi4overbeta1
    << ", beta0beta1overrho = " << SQRTbeta0beta1overrho
    << ", pow( ( _P / M_data -> Beta0(indz) + 1 ), 0.5 ) = " << pow( ( _P / M_data -> Beta0(indz) + 1 ), 0.5 ) << "\n";
    Debug(6320) << "[OneDimensionalModel_Physics_NonLinear::W_from_P] add term = " << add << "\n";

    if ( i == 1 )
        return _W - add;
    if ( i == 2 )
        return _W + add;

    ERROR_MSG("You can only find W1 or W2 as function of P");
    return -1.;
}

Real
OneDimensionalModel_Physics_NonLinear::W_from_Q( const Real& _Q, const Real& _W_n,
                                                 const Real& _W, const ID& i, const UInt& indz ) const
{
    Real K0( M_data -> Beta1(indz) / ( std::sqrt(M_data -> RobertsonCorrection()) * 4 ) );

    Real K1( (M_data -> Area0(indz) / 2) );
    K1 *= pow( M_data -> DensityRho() / (M_data -> Beta0(indz) * M_data -> Beta1(indz)), 1/M_data -> Beta1(indz) );
    K1 *= pow( K0, 2/M_data -> Beta1(indz) );

    Real w_k = _W_n;
    Real f_k, df_k, tau_k(0);

    if ( i == 1 ) // W1 given
    {
        f_k = pow( _W - w_k + Celerity0(indz) / K0, 2/M_data -> Beta1(indz) );
        tau_k = pow( _W - w_k + Celerity0(indz) / K0, 2/M_data -> Beta1(indz) );
        df_k = (-2 / M_data -> Beta1(indz)) * pow( _W - w_k + Celerity0(indz) / K0, 2/M_data -> Beta1(indz) - 1 );
    }
    if ( i == 2 ) // W2 given
    {
        f_k = pow( w_k - _W + Celerity0(indz) / K0, 2/M_data -> Beta1(indz) );
        tau_k = pow( w_k - _W + Celerity0(indz) / K0, 2/M_data -> Beta1(indz) );
        df_k = (-2 / M_data -> Beta1(indz)) * pow( w_k - _W + Celerity0(indz) / K0, 2/M_data -> Beta1(indz) - 1 );
    }
    f_k *= (_W + w_k);
    f_k += - _Q / K1;
    df_k *= (_W + w_k);
    df_k += f_k;

    Debug(6320) << "[OneDimensionalModel_Physics_NonLinear::W_from_Q] "
    << "K0 = " << K0
    << ", K1 = " << K1
    << ", tau_k = " << tau_k << "\n";

    Real w_kp1 = _Q / (K1 * tau_k) - _W;

    /* for debugging purposes
    std::ofstream ofile;
    ofile.open( "imposed_W_from_Q.m", std::ios::app );
    ofile << _Q << "\t"
      << w_kp1 << "\n";
    ofile.close();
    */

    return w_kp1;

    //    ERROR_MSG("You can only find W1 or W2 as function of P");
    //    return -1.;
}

}
