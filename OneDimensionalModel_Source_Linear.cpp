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
 *  @brief File containing a class for the linear source function B of the 1D hyperbolic problem.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-04-2010
 */

#include "OneDimensionalModel_Source_Linear.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_Source_Linear::OneDimensionalModel_Source_Linear() :
    super    ()
{}

OneDimensionalModel_Source_Linear::OneDimensionalModel_Source_Linear( const Physics_PtrType Physics ) :
    super    ( Physics )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_Source_Linear::operator()( const Real& _U1, const Real& _U2,
                                               const ID& ii,    const UInt& indz ) const
{
    if( ii == 1 ) // S1
    {
        return M_Physics->Data()->Source10( indz ) +
               M_Physics->Data()->Source11( indz ) * _U1 +
               M_Physics->Data()->Source12( indz ) * _U2;
    }
    if( ii == 2 ) // S2
    {
        return M_Physics->Data()->Source20( indz ) +
               M_Physics->Data()->Source21( indz ) * _U1 +
               M_Physics->Data()->Source22( indz ) * _U2;
    }
    ERROR_MSG("The flux function has only 2 components.");
    return -1.;
}

Real
OneDimensionalModel_Source_Linear::diff( const Real& /*_U1*/, const Real& /*_U2*/,
                                         const ID& ii,        const ID& jj,
                                         const UInt& indz ) const
{
    if( ii == 1 && jj == 1) // dS1/dU1 = 0
    {
        return M_Physics->Data()->Source11( indz );
    }
    if( ii == 1 && jj == 2) // dS1/dU2 = 0
    {
        return M_Physics->Data()->Source12( indz );
    }
    if( ii == 2 && jj == 1 ) // dS2/dU1
    {
        return M_Physics->Data()->Source21( indz );
    }
    if( ii == 2 && jj == 2 ) // dS2/dU2
    {
        return M_Physics->Data()->Source22( indz );
    }
    ERROR_MSG("Source's differential function has only 4 components.");
    return -1.;
}

//Real
//OneDimensionalModel_Source_Linear::diff2( const Real& /*_U1*/, const Real& /*_U2*/,
//                                          const ID& ii,        const ID& jj, const ID& kk,
//                                          const UInt& /*indz*/ ) const
//{
//    if( (0 < ii && ii < 3) && (0 < jj && jj < 3) && (0 < kk && kk < 3) )
//    {
//        return 0.;
//    }
//    ERROR_MSG("Source's second differential function has only 8 components.");
//    return -1.;
//}

Real
OneDimensionalModel_Source_Linear::interpolatedQuasiLinearSource( const Real& _U1, const Real& _U2,
                                                                  const ID& ii,    const Container2D_Type& bcNodes, const Real& /*cfl*/ ) const
{
    //TODO Implement the interpolation as done for the non-linear case
    return this->operator()(_U1, _U2, ii, bcNodes[0]);
}

}
