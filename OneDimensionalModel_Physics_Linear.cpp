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
    @brief File containing a class providing linear physical operations for the 1D model data.

    @version 1.0
    @date 01-07-2004
    @author Vincent Martin

    @version 2.0
    @date 13-04-2010
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>

    @contributor Simone Rossi <simone.rossi@epfl.ch>

    @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/OneDimensionalModel_Physics_Linear.hpp>

namespace LifeV
{
// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_Physics_Linear::fromUToW(       Real& W1,       Real& W2,
                                                    const Real& U1, const Real& U2, const UInt& indz ) const
{
    W1 = U2 + celerity0(indz) * ( U1 - M_data->area0(indz) );

    W2 = U2 - celerity0(indz) * ( U1 - M_data->area0(indz) );

    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::fromUToW] Q " << U2 << "\n";
    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::fromUToW] W1 " << W1 << "\n";
    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::fromUToW] W2 " << W2 << "\n";
    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::fromUToW] celerity " << celerity0(indz) << "\n";
    Debug( 6320 ) << "[OneDimensionalModel_Physics_Linear::fromUToW] ( _U1 - area0(indz) ) " << ( U1 - M_data->area0(indz) ) << "\n";
}

void
OneDimensionalModel_Physics_Linear::fromWToU(       Real& U1,       Real& U2,
                                                    const Real& W1, const Real& W2, const UInt& indz ) const
{
    U1 = M_data -> area0(indz) + ( W1 - W2) / ( 2 * celerity0(indz) );

    U2 = ( W1 + W2 ) / 2;
}

Real
OneDimensionalModel_Physics_Linear::fromWToP( const Real& W1, const Real& W2, const UInt& indz ) const
{
    return ( M_data -> beta0(indz)
             * ( std::pow( 1 / M_data->area0(indz), M_data -> beta1(indz) )
                 * std::pow( (W1 - W2 ) / ( 2 * celerity0(indz) ) + M_data -> area0(indz), M_data -> beta1(indz) )
                 - 1 )
           );
}

Real
OneDimensionalModel_Physics_Linear::dPdW( const Real& W1, const Real& W2,
                                                    const ID& i,     const UInt& indz ) const
{
    Real beta0beta1overA0beta1 ( M_data->beta0(indz) * M_data -> beta1(indz) / std::pow( M_data -> area0(indz), M_data -> beta1(indz) ) );

    Real oneover2celerity( 1 / ( 2 * celerity0(indz) ) );

    Real result( beta0beta1overA0beta1 * oneover2celerity );
    result *= ( ( W1 - W2 ) * oneover2celerity + M_data -> area0(indz) );

    if ( i == 1 ) //! dP/dW1
        return result;

    if ( i == 2 ) //! dP/dW2
        return -result;

    ERROR_MSG("P(W1,W2)'s differential function has only 2 components.");
    return -1.;
}

Real
OneDimensionalModel_Physics_Linear::fromPToW( const Real& P, const Real& W, const ID& i, const UInt& indz ) const
{
    Real add( 2 * celerity0(indz) * M_data -> area0(indz) * ( pow( ( P / M_data -> beta0(indz) + 1 ), 1 / M_data -> beta1(indz) ) - 1 ) );

    Debug(6320) << "[fromPToW] "
    << "2 * celerity0(indz) * area0(indz) = " << 2 * celerity0(indz) * M_data -> area0(indz)
    << ", pow( ( P / beta0(indz) + 1 ), 1 / beta1(indz) ) = "
    << pow( ( P / M_data -> beta0(indz) + 1 ), 1 / M_data -> beta1(indz) ) << "\n";
    Debug(6320) << "[fromPToW] add term = " << add << "\n";

    if ( i == 1 )
        return W - add;
    if ( i == 2 )
        return W + add;

    ERROR_MSG("You can only find W1 or W2 as function of P");
    return -1.;
}

Real
OneDimensionalModel_Physics_Linear::fromQToW( const Real& Q, const Real& /*W_n*/, const Real& W, const ID& i, const UInt& /*indz*/ ) const
{
    Real add( 2 * Q );

    if ( i == 1 ) // W1 given
        return add - W;

    if ( i == 2 ) // W2 given
        return add - W;

    ERROR_MSG("You can only find W1 or W2 as function of Q");
    return -1.;
}

}
