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
 *  @brief File containing a class for linear 1D model flux function.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#include "OneDimensionalModel_Flux_Linear.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_Flux_Linear::OneDimensionalModel_Flux_Linear() :
    super   ()
{}

OneDimensionalModel_Flux_Linear::OneDimensionalModel_Flux_Linear( const Physics_PtrType Physics ) :
    super   ( Physics )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_Flux_Linear::operator()( const Real& _U1, const Real& _U2,
                                             const ID& ii,    const UInt& indz ) const
{
    if( ii == 1 ) // F1
    {
        return M_Physics->Data()->Flux11( indz ) * _U1 + M_Physics->Data()->Flux12( indz ) * _U2;
    }
    if( ii == 2 ) // F2
    {
        return M_Physics->Data()->Flux21( indz ) * _U1 + M_Physics->Data()->Flux22( indz ) * _U2;
    }
    ERROR_MSG("The flux function has only 2 components.");
    return -1.;
}

Real
OneDimensionalModel_Flux_Linear::diff( const Real& /*_U1*/, const Real& /*_U2*/,
                                       const ID& ii, const ID& jj, const UInt& indz) const
{
    if( ii == 1 && jj == 1 ) // dF1/dU1
    {
        return M_Physics->Data()->Flux11( indz );
    }
    if( ii == 1 && jj == 2 ) // dF1/dU2
    {
        return M_Physics->Data()->Flux12( indz );
    }
    if( ii == 2 && jj == 1 ) // dF2/dU1
    {
        return M_Physics->Data()->Flux21( indz );
    }
    if( ii == 2 && jj == 2 ) // dF2/dU2
    {
        return M_Physics->Data()->Flux22( indz );
    }

    ERROR_MSG("Flux's differential function has only 4 components.");
    return -1.;
}

void
OneDimensionalModel_Flux_Linear::jacobian_EigenValues_Vectors( const Real& /*_U1*/, const Real& /*_U2*/,
                                                                     Real& eig1,          Real& eig2,
                                                                     Real& lefteigvec11,  Real& lefteigvec12,
                                                                     Real& lefteigvec21,  Real& lefteigvec22,
                                                               const UInt& indz ) const
{
    // eigen values
    eig1 = M_Physics->Data()->Celerity1( indz -1 );
    eig2 = M_Physics->Data()->Celerity2( indz -1 );

    // eigen vectors
    lefteigvec11 = M_Physics->Data()->LeftEigenVector11( indz -1 );
    lefteigvec12 = M_Physics->Data()->LeftEigenVector12( indz -1 );
    lefteigvec21 = M_Physics->Data()->LeftEigenVector21( indz -1 );
    lefteigvec22 = M_Physics->Data()->LeftEigenVector22( indz -1 );
}

Real
OneDimensionalModel_Flux_Linear::diff2( const Real& /*_U1*/, const Real& /*_U2*/,
                                        const ID& ii, const ID& jj, const ID& kk,
                                        const UInt& /*indz*/ ) const
{
    if( (0 < ii && ii < 3) && (0 < jj && jj < 3) && (0 < kk && kk < 3) )
    {
        return 0.;
    }
    ERROR_MSG("Flux's second differential function has only 8 components.");
    return -1.;
}

}
