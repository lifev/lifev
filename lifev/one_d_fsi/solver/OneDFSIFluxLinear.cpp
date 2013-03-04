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
 *  @maintainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/one_d_fsi/solver/OneDFSIFluxLinear.hpp>

namespace LifeV
{

// ===================================================
// Methods
// ===================================================
Real
OneDFSIFluxLinear::flux ( const Real& U1, const Real& U2, const ID& row, const UInt& iNode ) const
{
    if ( row == 0 ) // F1
    {
        return M_physicsPtr->data()->flux11 ( iNode ) * U1 + M_physicsPtr->data()->flux12 ( iNode ) * U2;
    }
    if ( row == 1 ) // F2
    {
        return M_physicsPtr->data()->flux21 ( iNode ) * U1 + M_physicsPtr->data()->flux22 ( iNode ) * U2;
    }
    ERROR_MSG ("The flux function has only 2 components.");
    return -1.;
}

Real
OneDFSIFluxLinear::dFdU ( const Real& /*U1*/, const Real& /*U2*/, const ID& row, const ID& column, const UInt& iNode) const
{
    if ( row == 0 && column == 0 ) // dF1/dU1
    {
        return M_physicsPtr->data()->flux11 ( iNode );
    }
    if ( row == 0 && column == 1 ) // dF1/dU2
    {
        return M_physicsPtr->data()->flux12 ( iNode );
    }
    if ( row == 1 && column == 0 ) // dF2/dU1
    {
        return M_physicsPtr->data()->flux21 ( iNode );
    }
    if ( row == 1 && column == 1 ) // dF2/dU2
    {
        return M_physicsPtr->data()->flux22 ( iNode );
    }

    ERROR_MSG ("Flux's differential function has only 4 components.");
    return -1.;
}

void
OneDFSIFluxLinear::eigenValuesEigenVectors ( const Real& /*U1*/, const Real& /*U2*/,
                                             container2D_Type& eigenvalues,
                                             container2D_Type& leftEigenvector1,
                                             container2D_Type& leftEigenvector2,
                                             const UInt& iNode ) const
{
    eigenvalues[0] = M_physicsPtr->data()->celerity1 ( iNode );
    eigenvalues[1] = M_physicsPtr->data()->celerity2 ( iNode );

    leftEigenvector1[0] = M_physicsPtr->data()->leftEigenVector11 ( iNode );
    leftEigenvector1[1] = M_physicsPtr->data()->leftEigenVector12 ( iNode );
    leftEigenvector2[0] = M_physicsPtr->data()->leftEigenVector21 ( iNode );
    leftEigenvector2[1] = M_physicsPtr->data()->leftEigenVector22 ( iNode );
}

void
OneDFSIFluxLinear::deltaEigenValuesEigenVectors ( const Real& /*U1*/, const Real& /*U2*/,
                                                  container2D_Type& deltaEigenvalues,
                                                  container2D_Type& deltaLeftEigenvector1,
                                                  container2D_Type& deltaLeftEigenvector2,
                                                  const UInt& /*iNode*/ ) const
{
    deltaEigenvalues[0] = 0;
    deltaEigenvalues[1] = 0;

    deltaLeftEigenvector1[0] = 0;
    deltaLeftEigenvector1[1] = 0;
    deltaLeftEigenvector2[0] = 0;
    deltaLeftEigenvector2[1] = 0;
}

}
