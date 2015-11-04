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
    @brief This file contains the definition of the methods for gradient recovery procedures

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 10 Jan 2012
 */

#ifndef GRADIENT_RECOVERY_HPP
#define GRADIENT_RECOVERY_HPP 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/CurrentFE.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/core/fem/ReferenceElement.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

namespace GradientRecovery
{

/*! Gradient recovery procedure from Zienkiewicz and Zhu.

  @param fespace The finite element space describing the data
  @param inputData The vector of data (pass it as repeated if possible)
  @param dxi The component to be recovered
  @return recovered gradient (unique map!)

 */
template<typename FESpaceType, typename VectorType>
VectorType ZZGradient (std::shared_ptr<FESpaceType> fespace,
                       const VectorType& inputData,
                       const UInt& dxi)
{
    // Repeated vector is needed
    if (inputData.mapType() != Repeated)
    {
        return ZZGradient (fespace, VectorType (inputData, Repeated), dxi);
    };

    // Get the area of the reference element
    Real refElemArea (0);

    switch ( fespace->refFE().shape() )
    {
        case TETRA:
            refElemArea = 1.0 / 6.0;
            break;
        case PRISM:
            refElemArea = 1.0 / 2.0;
            break;
        case HEXA:
            refElemArea = 1.0;
            break;
        case QUAD:
            refElemArea = 1.0;
            break;
        case TRIANGLE:
            refElemArea = 1.0 / 2.0;
            break;
        case LINE:
            refElemArea = 1.0;
            break;
        case POINT:
            refElemArea = 1.0; // Makes things consistent afterwards
            break;
        default:
            std::cerr << "ZZ Gradient Recovery: unknown shape! Aborting. " << std::endl;
            std::abort();
    }

    // Define the specific QR to be used
    // so that values in the QR correspond to
    // the values in the nodes

    QuadratureRule interpQuad;
    interpQuad.setDimensionShape ( shapeDimension (fespace->refFE().shape() ) , fespace->refFE().shape() );

    Real wQuad (refElemArea / fespace->refFE().nbDof() );

    for (UInt iQuadPt (0); iQuadPt < fespace->refFE().nbDof(); ++ iQuadPt)
    {
        interpQuad.addPoint (QuadraturePoint ( fespace->refFE().xi (iQuadPt),
                                               fespace->refFE().eta (iQuadPt),
                                               fespace->refFE().zeta (iQuadPt),
                                               wQuad) );
    }

    // Initialization of the two vectors

    VectorType patchArea ( inputData, Repeated );
    patchArea *= 0.0;

    VectorType gradientSum ( inputData, Repeated );
    gradientSum *= 0.0;

    // Build the structure and the constants

    CurrentFE interpCFE ( fespace->refFE(), getGeometricMap (*fespace->mesh() ), interpQuad );

    const UInt nbElement (fespace->mesh()->numElements() );
    const UInt nbLocalDof (fespace->dof().numLocalDof() );

    // Now loop over the elements

    for (UInt iElement (0); iElement < nbElement; ++iElement)
    {
        interpCFE.update ( fespace->mesh()->element (iElement), UPDATE_DPHI | UPDATE_WDET);

        for (UInt iDof (0); iDof < nbLocalDof; ++iDof)
        {
            for (UInt iDim (0); iDim < fespace->fieldDim(); ++iDim)
            {
                UInt globaliDofID ( fespace->dof().localToGlobalMap (iElement, iDof)
                                    + iDim * fespace->dof().numTotalDof() );

                patchArea[globaliDofID] += interpCFE.measure();

                for (UInt jDof (0); jDof < nbLocalDof; ++jDof)
                {
                    UInt globaljDofID ( fespace->dof().localToGlobalMap (iElement, jDof)
                                        + iDim * fespace->dof().numTotalDof() );

                    gradientSum[globaliDofID] += interpCFE.measure()
                                                 * inputData[globaljDofID]
                                                 * interpCFE.dphi (jDof, dxi, iDof);
                }
            }
        }
    }

    // Assembly

    return VectorType (gradientSum, Unique, Add) / VectorType (patchArea, Unique, Add);
}


/*! Laplacian recovery following Zienkiewicz and Zhu.

  @param fespace The finite element space describing the data
  @param inputData The vector of data (pass it as repeated if possible)
  @return recovered laplacian (unique map!)

 */
template<typename FESpaceType, typename VectorType>
VectorType ZZLaplacian (std::shared_ptr<FESpaceType> fespace,
                        const VectorType& inputData)
{
    // Need a repeated input
    if (inputData.mapType() != Repeated)
    {
        return ZZLaplacian (fespace, VectorType (inputData, Repeated) );
    }

    VectorType laplacian (inputData, Unique);
    laplacian *= 0.0;

    // Loop over the components
    for (UInt iDim (0); iDim < fespace->fieldDim(); ++iDim)
    {
        laplacian += ZZGradient (fespace, ZZGradient (fespace, inputData, iDim) , iDim);
    }

    return laplacian;
}


} // Namespace GradientRecovery

} // Namespace LifeV

#endif /* GRADIENTRECOVERY_H */
