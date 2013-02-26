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
    @brief File containing the procedures for the local assembly of the differential operators

    @contributor Claudia Colciago <claudia.colciago@epfl.ch>
    @mantainer Claudia Colciago <claudia.colciago@epfl.ch>

 */


#ifndef _ELEMOPERBD_H_INCLUDED
#define _ELEMOPERBD_H_INCLUDED

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/CurrentBoundaryFE.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/fem/DOF.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{
//! @name Public typedefs
//@{
typedef boost::numeric::ublas::matrix<Real> Matrix;
typedef boost::numeric::ublas::vector<Real> Vector;
typedef boost::numeric::ublas::zero_matrix<Real> ZeroMatrix;
//@}

/*! /namespace AssemblyElementalBoundary

  This namespace is specially designed to contain the elementary
  operations (corresponding to differential operators) that build
  the local contributions to be used in the assembly procedures.

 */
template<typename GeoShapeType>
class AssemblyElementalBoundary
{

private:

    AssemblyElementalBoundary() {}

public:

    static AssemblyElementalBoundary& instance()
    {

        static AssemblyElementalBoundary singleton;

        return singleton;

    }

    ~AssemblyElementalBoundary() {}

    void laplaceBeltrami (MatrixElemental& localLB,
                          const CurrentFE& LBCFE,
                          const CurrentBoundaryFE& LBCBdFE,
                          const Real& coefficient,
                          const ID    LBCFEID,
                          const UInt& fieldDim);

    void stiffStrainBoundary (MatrixElemental& localLB,
                              const CurrentFE& LBCFE,
                              const CurrentBoundaryFE& LBCBdFE,
                              const Real& coefficient,
                              const ID    LBCFEID,
                              const UInt& fieldDim);


    void divDivBoundary (MatrixElemental& localLB,
                         const CurrentFE& LBCFE,
                         const CurrentBoundaryFE& LBCBdFE,
                         const Real& coefficient,
                         const ID    LBCFEID,
                         const UInt& fieldDim);

};


template<typename GeoShapeType>
void AssemblyElementalBoundary<GeoShapeType>::laplaceBeltrami (MatrixElemental& localLB,
                                                               const CurrentFE& LBCFE,
                                                               const CurrentBoundaryFE& LBCBdFE,
                                                               const Real& coefficient,
                                                               const ID    LBCFEID,
                                                               const UInt& fieldDim)
{

    typedef GeoShapeType geoShape_Type;

    const UInt nbFEDof (LBCBdFE.nbNode() );
    const UInt nbQuadPt (LBCBdFE.nbQuadPt() );
    Real localValue (0);
    Real normalGradj (0);
    Real normalGradi (0);

    UInt iDofE, jDofE;


    // Assemble the local diffusion
    for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
    {
        // Extract the view of the matrix
        MatrixElemental::matrix_view localView = localLB.block (iterFDim, iterFDim);

        // Loop over the basis functions
        for (UInt iDofFace (0); iDofFace < nbFEDof ; ++iDofFace)
        {
            iDofE = geoShape_Type::faceToPoint ( LBCFEID, iDofFace ); // local vertex number (in element)

            // Build the local matrix only where needed:
            // Lower triangular + diagonal parts
            for (UInt jDofFace (0); jDofFace <= iDofFace; ++jDofFace)
            {
                localValue = 0.0;

                jDofE = geoShape_Type::faceToPoint ( LBCFEID, jDofFace ); // local vertex number (in element)

                //Loop on the quadrature nodes
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                {

                    normalGradj = 0;
                    normalGradi = 0;

                    for ( int iDim = 0; iDim < 3; ++iDim )
                    {
                        normalGradi += LBCFE.dphi ( iDofE, iDim, iQuadPt ) * LBCBdFE.normal ( iDim, iQuadPt );
                        normalGradj += LBCFE.dphi ( jDofE, iDim, iQuadPt ) * LBCBdFE.normal ( iDim, iQuadPt );
                    }

                    for (UInt iDim (0); iDim < 3; ++iDim)
                    {
                        localValue += ( LBCFE.dphi (iDofE, iDim, iQuadPt) - normalGradi * LBCBdFE.normal ( iDim, iQuadPt) )
                                      * ( LBCFE.dphi (jDofE, iDim, iQuadPt) - normalGradj * LBCBdFE.normal ( iDim, iQuadPt) )
                                      * LBCBdFE.weightMeas (iQuadPt);
                    }
                }

                localValue *= coefficient;

                // Add on the local matrix
                localView (iDofE, jDofE) += localValue;

                if (iDofFace != jDofFace)
                {
                    localView (jDofE, iDofE) += localValue;
                }
            }
        }
    }
}


template<typename GeoShapeType>
void AssemblyElementalBoundary<GeoShapeType>::divDivBoundary (MatrixElemental& localLB,
                                                              const CurrentFE& LBCFE,
                                                              const CurrentBoundaryFE& LBCBdFE,
                                                              const Real& coefficient,
                                                              const ID    LBCFEID,
                                                              const UInt& fieldDim)
{

    typedef GeoShapeType geoShape_Type;

    const UInt nbFEDof (LBCBdFE.nbNode() );
    const UInt nbQuadPt (LBCBdFE.nbQuadPt() );
    Real localValue (0);
    Real normalGradj (0);
    Real normalGradi (0);

    UInt iDofE, jDofE;


    // // Assemble the local diffusion
    // for (UInt iterFDim(0); iterFDim<fieldDim; ++iterFDim)
    // {
    //     // Extract the view of the matrix
    //     MatrixElemental::matrix_view localView = localLB.block(iterFDim,iterFDim);

    //     // Loop over the basis functions
    //     for (UInt iDofFace(0); iDofFace < nbFEDof ; ++iDofFace)
    //     {
    //      iDofE = geoShape_Type::faceToPoint( LBCFEID, iDofFace ); // local vertex number (in element)

    //         // Build the local matrix only where needed:
    //         // Lower triangular + diagonal parts
    //         for (UInt jDofFace(0); jDofFace <= iDofFace; ++jDofFace)
    //         {
    //             localValue = 0.0;

    //      jDofE = geoShape_Type::faceToPoint( LBCFEID, jDofFace ); // local vertex number (in element)

    //             //Loop on the quadrature nodes
    //             for (UInt iQuadPt(0); iQuadPt < nbQuadPt; ++iQuadPt)
    //             {

    //          normalGradj = 0;
    //          normalGradi = 0;

    //          for ( int iDim = 0; iDim < 3; ++iDim )
    //          {
    //              normalGradi += LBCFE.dphi( iDofE, iDim, iQuadPt ) * LBCBdFE.normal( iDim, iQuadPt);
    //              normalGradj += LBCFE.dphi( jDofE, iDim, iQuadPt ) * LBCBdFE.normal( iDim, iQuadPt);
    //          }


    //          localValue += ( LBCFE.dphi(iDofE,iterFDim,iQuadPt) -
    //                             normalGradi * LBCBdFE.normal( iterFDim, iQuadPt) )
    //                     * ( LBCFE.dphi(jDofE,iterFDim,iQuadPt)-
    //                             normalGradj * LBCBdFE.normal( iterFDim, iQuadPt)  )
    //                                   * LBCBdFE.weightMeas(iQuadPt);
    //             }

    //             localValue*=coefficient;

    //             // Add on the local matrix
    //             localView(iDofE,jDofE)+=localValue;

    //             if (iDofFace != jDofFace)
    //             {
    //                 localView(jDofE,iDofE)+=localValue;
    //             }
    //         }
    //     }
    // }

    for ( UInt iFDim (0); iFDim < fieldDim; ++iFDim )
    {
        for ( UInt jFDim (0); jFDim < fieldDim; ++jFDim )
        {
            MatrixElemental::matrix_view localView = localLB.block ( iFDim, jFDim );

            for ( UInt iDofFace (0); iDofFace < nbFEDof; ++iDofFace )
            {

                iDofE = geoShape_Type::faceToPoint ( LBCFEID, iDofFace ); // local vertex number (in element)

                for ( UInt jDofFace (0); jDofFace < nbFEDof; ++jDofFace )
                {
                    jDofE = geoShape_Type::faceToPoint ( LBCFEID, jDofFace ); // local vertex number (in element)

                    localValue = 0.0;

                    for ( UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt )
                    {

                        normalGradj = 0;
                        normalGradi = 0;

                        for ( int iDim = 0; iDim < 3; ++iDim )
                        {
                            normalGradi += LBCFE.dphi ( iDofE, iDim, iQuadPt ) * LBCBdFE.normal ( iDim, iQuadPt );
                            normalGradj += LBCFE.dphi ( jDofE, iDim, iQuadPt ) * LBCBdFE.normal ( iDim, iQuadPt );
                        }

                        localValue += ( LBCFE.dphi (iDofE, iFDim, iQuadPt) -
                                        normalGradi * LBCBdFE.normal ( iFDim, iQuadPt) )
                                      * ( LBCFE.dphi (jDofE, jFDim, iQuadPt) -
                                          normalGradj * LBCBdFE.normal ( jFDim, iQuadPt)  )
                                      * LBCBdFE.weightMeas (iQuadPt);

                    }

                    localView ( iDofE, jDofE ) += coefficient * localValue;
                }
            }
        }
    }
}


template<typename GeoShapeType>
void AssemblyElementalBoundary<GeoShapeType>::stiffStrainBoundary (MatrixElemental& localLB,
                                                                   const CurrentFE& LBCFE,
                                                                   const CurrentBoundaryFE& LBCBdFE,
                                                                   const Real& coefficient,
                                                                   const ID    LBCFEID,
                                                                   const UInt& fieldDim)
{
    typedef GeoShapeType geoShape_Type;

    const UInt nbFEDof (LBCBdFE.nbNode() );
    const UInt nbQuadPt (LBCBdFE.nbQuadPt() );
    Real localValue (0);
    const Real newCoefficient (coefficient * 0.5);
    Real normalGradj (0);
    Real normalGradi (0);

    UInt iDofE, jDofE;

    // Assemble the local diffusion
    for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
    {
        // Extract the view of the matrix
        MatrixElemental::matrix_view localView = localLB.block (iterFDim, iterFDim);

        // Loop over the basis functions
        for (UInt iDofFace (0); iDofFace < nbFEDof ; ++iDofFace)
        {

            iDofE = geoShape_Type::faceToPoint ( LBCFEID, iDofFace ); // local vertex number (in element)

            // Build the local matrix only where needed:
            // Lower triangular + diagonal parts
            for (UInt jDofFace (0); jDofFace <= iDofFace; ++jDofFace)
            {

                jDofE = geoShape_Type::faceToPoint ( LBCFEID, jDofFace ); // local vertex number (in element)

                localValue = 0.0;

                //Loop on the quadrature nodes
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                {

                    normalGradj = 0;
                    normalGradi = 0;

                    for ( int iDim = 0; iDim < 3; ++iDim )
                    {
                        normalGradi += LBCFE.dphi ( iDofE, iDim, iQuadPt ) * LBCBdFE.normal ( iDim, iQuadPt );
                        normalGradj += LBCFE.dphi ( jDofE, iDim, iQuadPt ) * LBCBdFE.normal ( iDim, iQuadPt );
                    }

                    for (UInt iDim (0); iDim < 3; ++iDim)
                    {
                        localValue += ( LBCFE.dphi (iDofE, iDim, iQuadPt) - normalGradi * LBCBdFE.normal ( iDim, iQuadPt) )
                                      * ( LBCFE.dphi (jDofE, iDim, iQuadPt) - normalGradj * LBCBdFE.normal ( iDim, iQuadPt) )
                                      * LBCBdFE.weightMeas (iQuadPt);

                    }
                }

                localValue *= newCoefficient;

                // Add on the local matrix
                localView (iDofE, jDofE) += localValue;

                if (iDofFace != jDofFace)
                {
                    localView (jDofE, iDofE) += localValue;
                }
            }
        }
    }

    for ( UInt iFDim (0); iFDim < fieldDim; ++iFDim )
    {
        for ( UInt jFDim (0); jFDim < fieldDim; ++jFDim )
        {
            MatrixElemental::matrix_view localView = localLB.block ( iFDim, jFDim );

            for ( UInt iDofFace (0); iDofFace < nbFEDof; ++iDofFace )
            {

                iDofE = geoShape_Type::faceToPoint ( LBCFEID, iDofFace ); // local vertex number (in element)

                for ( UInt jDofFace (0); jDofFace < nbFEDof; ++jDofFace )
                {
                    jDofE = geoShape_Type::faceToPoint ( LBCFEID, jDofFace ); // local vertex number (in element)

                    localValue = 0.0;

                    for ( UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt )
                    {

                        normalGradj = 0;
                        normalGradi = 0;

                        for ( int iDim = 0; iDim < 3; ++iDim )
                        {
                            normalGradi += LBCFE.dphi ( jDofE, iDim, iQuadPt ) * LBCBdFE.normal ( iDim, iQuadPt );
                            normalGradj += LBCFE.dphi ( iDofE, iDim, iQuadPt ) * LBCBdFE.normal ( iDim, iQuadPt );
                        }

                        localValue += ( LBCFE.dphi (jDofE, iFDim, iQuadPt) - normalGradi * LBCBdFE.normal ( iFDim, iQuadPt) )
                                      * ( LBCFE.dphi (iDofE, jFDim, iQuadPt) - normalGradj * LBCBdFE.normal ( jFDim, iQuadPt) )
                                      * LBCBdFE.weightMeas (iQuadPt);
                    }

                    localView ( iDofE, jDofE ) += newCoefficient * localValue;
                }
            }
        }
    }
}


} // namespace LifeV
#endif
