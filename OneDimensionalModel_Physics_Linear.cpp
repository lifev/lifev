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
 *  @brief File containing a class providing linear physical operations for the 1D model data.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date 01-07-2004
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 13-04-2010
 */

#include <lifemc/lifesolver/OneDimensionalModel_Physics_Linear.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
OneDimensionalModel_Physics_Linear::OneDimensionalModel_Physics_Linear()
{
}

OneDimensionalModel_Physics_Linear::OneDimensionalModel_Physics_Linear( const Data_PtrType Data ) :
        super   ( Data )
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_Physics_Linear::W_from_U(       Real& _W1,       Real& _W2,
                                                    const Real& _U1, const Real& _U2, const UInt& indz ) const
{
    _W1 = _U2 + Celerity0(indz) * ( _U1 - M_Data->Area0(indz) );

    _W2 = _U2 - Celerity0(indz) * ( _U1 - M_Data->Area0(indz) );

    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::W_from_U] Q " << _U2 << "\n";
    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::W_from_U] W1 " << _W1 << "\n";
    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::W_from_U] W2 " << _W2 << "\n";
    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::W_from_U] Celerity " << Celerity0(indz) << "\n";
    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::W_from_U] ( _U1 - Area0(indz) ) " << ( _U1 - M_Data->Area0(indz) ) << "\n";
}

void
OneDimensionalModel_Physics_Linear::U_from_W(       Real& _U1,       Real& _U2,
                                                    const Real& _W1, const Real& _W2, const UInt& indz ) const
{
    _U1 = M_Data->Area0(indz) + (_W1 - _W2) / ( 2 * Celerity0(indz) );

    _U2 = ( _W1 + _W2 ) / 2;
}

Real
OneDimensionalModel_Physics_Linear::pressure_W( const Real& _W1, const Real& _W2, const UInt& indz ) const
{
    return ( M_Data->Beta0(indz)
             * ( std::pow( 1 / M_Data->Area0(indz), M_Data->Beta1(indz) )
                 * std::pow( (_W1 - _W2) / ( 2*Celerity0(indz) ) + M_Data->Area0(indz), M_Data->Beta1(indz) )
                 - 1 )
           );
}

Real
OneDimensionalModel_Physics_Linear::pressure_WDiff( const Real& _W1, const Real& _W2,
                                                    const ID& i,     const UInt& indz ) const
{
    Real beta0beta1overA0beta1 ( M_Data->Beta0(indz) * M_Data->Beta1(indz) / std::pow( M_Data->Area0(indz), M_Data->Beta1(indz) ) );

    Real oneover2celerity( 1 / ( 2 * Celerity0(indz) ) );

    Real result( beta0beta1overA0beta1 * oneover2celerity );
    result *= ( ( _W1 - _W2 ) * oneover2celerity + M_Data->Area0(indz) );

    if ( i == 1 ) //! dP/dW1
        return result;

    if ( i == 2 ) //! dP/dW2
        return -result;

    ERROR_MSG("P(W1,W2)'s differential function has only 2 components.");
    return -1.;
}

Real
OneDimensionalModel_Physics_Linear::W_from_P( const Real& _P, const Real& _W, const ID& i, const UInt& indz ) const
{
    Real add( 2 * Celerity0(indz) * M_Data->Area0(indz)
              * ( pow( ( _P / M_Data->Beta0(indz) + 1 ), 1 / M_Data->Beta1(indz) ) - 1 ) );

    Debug(6320) << "[W_from_P] "
    << "2 * Celerity0(indz) * Area0(indz) = " << 2 * Celerity0(indz) * M_Data->Area0(indz)
    << ", pow( ( _P / Beta0(indz) + 1 ), 1 / Beta1(indz) ) = "
    << pow( ( _P / M_Data->Beta0(indz) + 1 ), 1 / M_Data->Beta1(indz)  ) << "\n";
    Debug(6320) << "[W_from_P] add term = " << add << "\n";

    if ( i == 1 )
        return _W - add;
    if ( i == 2 )
        return _W + add;

    ERROR_MSG("You can only find W1 or W2 as function of P");
    return -1.;
}

Real
OneDimensionalModel_Physics_Linear::W_from_Q( const Real& _Q, const Real& /*_W_n*/, const Real& _W, const ID& i, const UInt& /*indz*/ ) const
{
    Real add( 2 * _Q );

    if ( i == 1 ) // W1 given
        return add - _W;

    if ( i == 2 ) // W2 given
        return add - _W;

    ERROR_MSG("You can only find W1 or W2 as function of Q");
    return -1.;
}

}
