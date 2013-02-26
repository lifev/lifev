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
 *  @brief File containing a base class for non linear 1D model flux function.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 15-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/one_d_fsi/solver/OneDFSIFluxNonLinear.hpp>

namespace LifeV
{
// ===================================================
// Methods
// ===================================================
Real
OneDFSIFluxNonLinear::flux ( const Real& A, const Real& Q, const ID& row,  const UInt& iNode ) const
{
    if ( row == 0 ) // F1
    {
        return Q;
    }

    if ( row == 1 ) // F2
    {
        return ( M_physicsPtr->data()->alpha ( iNode ) * Q * Q / A +
                 M_physicsPtr->data()->beta0 ( iNode ) * M_physicsPtr->data()->beta1 ( iNode ) *
                 M_physicsPtr->data()->area0 ( iNode ) / ( ( M_physicsPtr->data()->beta1 ( iNode ) + 1 ) *
                                                           M_physicsPtr->data()->densityRho() ) * ( OneDFSI::pow15 ( A / M_physicsPtr->data()->area0 ( iNode ),
                                                                   M_physicsPtr->data()->beta1 ( iNode ) + 1 ) - 1 ) ) *
               M_physicsPtr->data()->robertsonCorrection();
    }

    ERROR_MSG ("The flux function has only 2 components.");

    return -1.;
}

Real
OneDFSIFluxNonLinear::dFdU ( const Real& A, const Real& Q, const ID& row, const ID& column, const UInt& iNode ) const
{
    if ( row == 0 && column == 0 ) // dF1/dA
    {
        return 0.;
    }

    if ( row == 0 && column == 1 ) // dF1/dQ
    {
        return 1.;
    }

    if ( row == 1 && column == 0 ) // dF2/dA
    {
        return ( M_physicsPtr->data()->beta0 ( iNode ) *
                 M_physicsPtr->data()->beta1 ( iNode ) /
                 M_physicsPtr->data()->densityRho() * OneDFSI::pow05 ( A / M_physicsPtr->data()->area0 ( iNode ),
                                                                       M_physicsPtr->data()->beta1 ( iNode ) ) -
                 M_physicsPtr->data()->alpha ( iNode ) * Q * Q / A / A ) *
               M_physicsPtr->data()->robertsonCorrection();
    }
    if ( row == 1 && column == 1 ) // dF2/dQ
    {
        return M_physicsPtr->data()->robertsonCorrection() * 2 * M_physicsPtr->data()->alpha ( iNode ) * Q / A;
    }

    ERROR_MSG ("Flux's differential function has only 4 components.");

    return -1.;
}

void
OneDFSIFluxNonLinear::eigenValuesEigenVectors ( const Real& A,
                                                const Real& Q,
                                                container2D_Type& eigenvalues,
                                                container2D_Type& leftEigenvector1,
                                                container2D_Type& leftEigenvector2,
                                                const UInt& iNode ) const
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream (6312) << "[OneDFSIModel_Flux_NonLinear]::jacobian_EigenValues_Vectors\n";
#endif

    Real celerity;
    celerity       = std::sqrt ( M_physicsPtr->data()->alpha ( iNode ) * (
                                     M_physicsPtr->data()->alpha ( iNode ) - 1) * Q * Q / ( A * A ) +
                                 M_physicsPtr->data()->beta0 ( iNode ) *
                                 M_physicsPtr->data()->beta1 ( iNode ) /
                                 M_physicsPtr->data()->densityRho() * OneDFSI::pow05 ( A / M_physicsPtr->data()->area0 ( iNode ),
                                         M_physicsPtr->data()->beta1 ( iNode ) ) );

    eigenvalues[0] = M_physicsPtr->data()->alpha ( iNode ) * Q / A + celerity;
    eigenvalues[1] = M_physicsPtr->data()->alpha ( iNode ) * Q / A - celerity;

    leftEigenvector1[0] = - eigenvalues[1] / A;
    leftEigenvector1[1] = 1. / A;
    leftEigenvector2[0] = - eigenvalues[0] / A;
    leftEigenvector2[1] = 1. / A;
}

void
OneDFSIFluxNonLinear::deltaEigenValuesEigenVectors ( const Real& A,
                                                     const Real& Q,
                                                     container2D_Type& deltaEigenvalues,
                                                     container2D_Type& deltaLeftEigenvector1,
                                                     container2D_Type& deltaLeftEigenvector2,
                                                     const UInt& iNode ) const
{
    Real deltaCelerity;

    Real AoverA0 ( A / M_physicsPtr->data()->area0 ( iNode ) );
    Real C ( OneDFSI::pow05 (  AoverA0, M_physicsPtr->data()->beta1 ( iNode ) ) / M_physicsPtr->data()->densityRho() );

    deltaCelerity  = 0.5 / std::sqrt ( M_physicsPtr->data()->alpha ( iNode ) * (
                                           M_physicsPtr->data()->alpha ( iNode ) - 1) * Q * Q / ( A * A ) +
                                       M_physicsPtr->data()->beta0 ( iNode ) *
                                       M_physicsPtr->data()->beta1 ( iNode ) * C ) * ( C * (
                                                   M_physicsPtr->data()->beta1 ( iNode ) *
                                                   M_physicsPtr->data()->dBeta0dz ( iNode ) -
                                                   M_physicsPtr->data()->beta0 ( iNode ) *
                                                   M_physicsPtr->data()->beta1 ( iNode ) *
                                                   M_physicsPtr->data()->beta1 ( iNode ) /
                                                   M_physicsPtr->data()->area0 ( iNode ) *
                                                   M_physicsPtr->data()->dArea0dz ( iNode ) +
                                                   M_physicsPtr->data()->beta0 ( iNode ) * ( 1 +
                                                           M_physicsPtr->data()->beta0 ( iNode ) * std::log ( AoverA0 ) ) *
                                                   M_physicsPtr->data()->dBeta1dz ( iNode ) ) + ( 2 *
                                                           M_physicsPtr->data()->alpha ( iNode ) - 1 ) * Q * Q / ( A * A ) *
                                               M_physicsPtr->data()->dAlphadz ( iNode ) );

    deltaEigenvalues[0] = M_physicsPtr->data()->dAlphadz ( iNode ) * Q / A + deltaCelerity;
    deltaEigenvalues[1] = M_physicsPtr->data()->dAlphadz ( iNode ) * Q / A - deltaCelerity;

    deltaLeftEigenvector1[0] = - deltaEigenvalues[1] / A;
    deltaLeftEigenvector1[1] = 0;
    deltaLeftEigenvector2[0] = - deltaEigenvalues[0] / A;
    deltaLeftEigenvector2[1] = 0;
}

}
