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
 *  @brief File containing a base class for linear 1D model flux function.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include "OneDimensionalModel_Flux_Linear.hpp"

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalFluxLinear::flux( const Real& U1, const Real& U2, const ID& ii, const UInt& i ) const
{
    if ( ii == 1 ) // F1
    {
        return M_physics->data()->flux11( i ) * U1 + M_physics->data()->flux12( i ) * U2;
    }
    if ( ii == 2 ) // F2
    {
        return M_physics->data()->flux21( i ) * U1 + M_physics->data()->flux22( i ) * U2;
    }
    ERROR_MSG("The flux function has only 2 components.");
    return -1.;
}

Real
OneDimensionalFluxLinear::dFdU( const Real& /*U1*/, const Real& /*U2*/, const ID& ii, const ID& jj, const UInt& i) const
{
    if ( ii == 1 && jj == 1 ) // dF1/dU1
    {
        return M_physics->data()->flux11( i );
    }
    if ( ii == 1 && jj == 2 ) // dF1/dU2
    {
        return M_physics->data()->flux12( i );
    }
    if ( ii == 2 && jj == 1 ) // dF2/dU1
    {
        return M_physics->data()->flux21( i );
    }
    if ( ii == 2 && jj == 2 ) // dF2/dU2
    {
        return M_physics->data()->flux22( i );
    }

    ERROR_MSG("Flux's differential function has only 4 components.");
    return -1.;
}

//Real
//OneDimensionalFluxLinear::diff2( const Real& /*U1*/, const Real& /*U2*/,
//                                        const ID& ii, const ID& jj, const ID& kk,
//                                        const UInt& /*i*/ ) const
//{
//    if( (0 < ii && ii < 3) && (0 < jj && jj < 3) && (0 < kk && kk < 3) )
//    {
//        return 0.;
//    }
//    ERROR_MSG("Flux's second differential function has only 8 components.");
//    return -1.;
//}

void
OneDimensionalFluxLinear::eigenValuesEigenVectors( const Real& /*U1*/, const Real& /*U2*/,
                                                   container2D_Type& eigenvalues,
                                                   container2D_Type& leftEigenvector1,
                                                   container2D_Type& leftEigenvector2,
                                                   const UInt& i ) const
{
    eigenvalues[0] = M_physics->data()->celerity1( i );
    eigenvalues[1] = M_physics->data()->celerity2( i );

    leftEigenvector1[0] = M_physics->data()->leftEigenVector11( i );
    leftEigenvector1[1] = M_physics->data()->leftEigenVector12( i );
    leftEigenvector2[0] = M_physics->data()->leftEigenVector21( i );
    leftEigenvector2[1] = M_physics->data()->leftEigenVector22( i );
}

void
OneDimensionalFluxLinear::deltaEigenValuesEigenVectors( const Real& /*U1*/, const Real& /*U2*/,
                                                        container2D_Type& deltaEigenvalues,
                                                        container2D_Type& deltaLeftEigenvector1,
                                                        container2D_Type& deltaLeftEigenvector2,
                                                        const UInt& /*i*/ ) const
{
    deltaEigenvalues[0] = 0;
    deltaEigenvalues[1] = 0;

    deltaLeftEigenvector1[0] = 0;
    deltaLeftEigenvector1[1] = 0;
    deltaLeftEigenvector2[0] = 0;
    deltaLeftEigenvector2[1] = 0;
}

}
