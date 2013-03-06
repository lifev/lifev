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
 *  @brief File containing a class for the non linear source function B of the 1D hyperbolic problem
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

#include <lifev/one_d_fsi/solver/OneDFSISourceNonLinear.hpp>

namespace LifeV
{
// ===================================================
// Methods
// ===================================================
Real
OneDFSISourceNonLinear::source ( const Real& A, const Real& Q, const ID& row, const UInt& iNode ) const
{
    if ( row == 0 ) // B1
    {
        return 0.;
    }
    if ( row == 1 ) // B2
    {
        Real beta1plus1 ( M_physicsPtr->data()->beta1 ( iNode ) + 1 );
        Real AoverA0 ( A / M_physicsPtr->data()->area0 ( iNode ) );
        Real C0 ( 1 / ( M_physicsPtr->data()->densityRho() * beta1plus1 ) );
        Real C ( 1 / ( M_physicsPtr->data()->densityRho() * beta1plus1 ) * OneDFSI::pow15 ( AoverA0, beta1plus1 ) );

        return ( M_physicsPtr->data()->friction() * Q / A
                 + Q * Q / A * M_physicsPtr->data()->dAlphadz ( iNode )
                 + C  * ( M_physicsPtr->data()->area0 ( iNode ) * M_physicsPtr->data()->dBeta0dz ( iNode )
                          - M_physicsPtr->data()->beta0 ( iNode ) * M_physicsPtr->data()->beta1 ( iNode ) * M_physicsPtr->data()->dArea0dz ( iNode )
                          + M_physicsPtr->data()->area0 ( iNode ) * M_physicsPtr->data()->beta0 ( iNode )
                          * ( std::log ( AoverA0 ) - 1. / beta1plus1 ) * M_physicsPtr->data()->dBeta1dz ( iNode )
                        )
                 - C0 * ( M_physicsPtr->data()->area0 ( iNode ) * M_physicsPtr->data()->dBeta0dz ( iNode )
                          - M_physicsPtr->data()->beta0 ( iNode ) * M_physicsPtr->data()->beta1 ( iNode ) * M_physicsPtr->data()->dArea0dz ( iNode )
                          - M_physicsPtr->data()->area0 ( iNode ) * M_physicsPtr->data()->beta0 ( iNode ) / beta1plus1 * M_physicsPtr->data()->dBeta1dz ( iNode )
                        )
                 + ( M_physicsPtr->data()->area0 ( iNode ) - A ) / M_physicsPtr->data()->densityRho() * M_physicsPtr->data()->dBeta0dz ( iNode )
               ) * M_physicsPtr->data()->robertsonCorrection();
    }

    ERROR_MSG ("The flux function has only 2 components.");

    return -1.;
}

Real
OneDFSISourceNonLinear::dSdU ( const Real& A, const Real& Q, const ID& row, const ID& column, const UInt& iNode) const
{
    if ( row == 0 ) // B1
    {
        if ( column == 0 || column == 1 ) // dB2/dUj = 0
        {
            return 0.;
        }
    }
    if ( row == 1 ) // B2
    {
        if ( column == 0 ) // dB2/dA
        {
            Real AoverA0 ( A / M_physicsPtr->data()->area0 ( iNode ) );
            Real C ( OneDFSI::pow05 ( AoverA0, M_physicsPtr->data()->beta1 ( iNode ) ) / M_physicsPtr->data()->densityRho() );

            return ( -M_physicsPtr->data()->friction() * Q / A / A
                     - Q * Q / ( A * A ) * M_physicsPtr->data()->dAlphadz ( iNode )
                     + C  * ( M_physicsPtr->data()->dBeta0dz ( iNode )
                              - M_physicsPtr->data()->beta0 ( iNode ) * M_physicsPtr->data()->beta1 ( iNode )
                              / M_physicsPtr->data()->area0 ( iNode ) * M_physicsPtr->data()->dArea0dz ( iNode )
                              + M_physicsPtr->data()->beta0 ( iNode ) * std::log ( AoverA0 ) * M_physicsPtr->data()->dBeta1dz ( iNode )
                            )
                     - 1. / M_physicsPtr->data()->densityRho() * M_physicsPtr->data()->dBeta0dz ( iNode )
                   ) * M_physicsPtr->data()->robertsonCorrection();
        }
        if ( column == 1 ) // dB2/dQ
        {
            return M_physicsPtr->data()->robertsonCorrection() * ( M_physicsPtr->data()->friction() / A +
                                                                   2 * Q / A * M_physicsPtr->data()->dAlphadz ( iNode ) );
        }
    }

    ERROR_MSG ("Source's differential function has only 4 components.");

    return -1.;
}

Real
OneDFSISourceNonLinear::interpolatedNonConservativeSource ( const Real& A, const Real& Q,
                                                            const ID& row, const container2D_Type& bcNodes, const Real& cfl ) const
{
    if ( row == 0 ) // QLS1
    {
        return 0.;
    }
    if ( row == 1 ) // QLS2
    {
        // Interpolate quantities
        Real area0      = ( 1 - cfl ) * M_physicsPtr->data()->area0 (bcNodes[0])    + cfl * M_physicsPtr->data()->area0 (bcNodes[1]);
        Real beta0      = ( 1 - cfl ) * M_physicsPtr->data()->beta0 (bcNodes[0])    + cfl * M_physicsPtr->data()->beta0 (bcNodes[1]);
        Real beta1      = ( 1 - cfl ) * M_physicsPtr->data()->beta1 (bcNodes[0])    + cfl * M_physicsPtr->data()->beta1 (bcNodes[1]);
        Real dArea0dz   = ( 1 - cfl ) * M_physicsPtr->data()->dArea0dz (bcNodes[0]) + cfl * M_physicsPtr->data()->dArea0dz (bcNodes[1]);
        Real dBeta0dz   = ( 1 - cfl ) * M_physicsPtr->data()->dBeta0dz (bcNodes[0]) + cfl * M_physicsPtr->data()->dBeta0dz (bcNodes[1]);
        Real dBeta1dz   = ( 1 - cfl ) * M_physicsPtr->data()->dBeta1dz (bcNodes[0]) + cfl * M_physicsPtr->data()->dBeta1dz (bcNodes[1]);
        Real dAlphadz   = ( 1 - cfl ) * M_physicsPtr->data()->dAlphadz (bcNodes[0]) + cfl * M_physicsPtr->data()->dAlphadz (bcNodes[1]);

        Real AoverA0 ( A / area0 );
        Real C ( A / M_physicsPtr->data()->densityRho() * OneDFSI::pow05 ( AoverA0, beta1 ) );

        return ( M_physicsPtr->data()->friction() * Q / A
                 + Q * Q / A * dAlphadz
                 + C * ( dBeta0dz
                         - beta0 * beta1 / area0 * dArea0dz
                         + beta0 * std::log ( AoverA0 ) * dBeta1dz
                       )
                 - A / M_physicsPtr->data()->densityRho() * dBeta0dz
               ) * M_physicsPtr->data()->robertsonCorrection();
    }

    ERROR_MSG ("The QLS function has only 2 components.");

    return -1.;
}

}
