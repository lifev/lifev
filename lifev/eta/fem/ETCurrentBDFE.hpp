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
    @brief A short description of the file content

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 21 Feb 2012

    A more detailed description of the file (if necessary)
 */

#ifndef ETCURRENTBDFE_H
#define ETCURRENTBDFE_H 1

#include <lifev/core/LifeV.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>

namespace LifeV
{

//! ETCurrentBDFE - Short description of the class

template<UInt spaceDim>
class ETCurrentBDFE
{
    //! @name Friends
    //@{

    //!Friend to allow direct access to the raw data
    template <UInt dim>
    friend class ExpressionAssembly::EvaluationHK;

    //!Friend to allow direct access to the raw data
    template <UInt dim>
    friend class ExpressionAssembly::EvaluationPosition;

    //@}

public:

    //! @name Constructors, destructor
    //@{

    //! Full constructor
    /*!
      @param refFE The reference element for the FE
      @param geoMap The geometric map from the reference element to the current element
      @param qr The quadrature rule
     */
    ETCurrentBDFE (const GeometricMap& geoMap, const QuadratureRule& qr)
        : M_geometricMap (&geoMap),
          M_quadrature (new QuadratureRule (qr) ),
          M_nbMapDof (geoMap.nbDof() ),
          M_nbQuadPt (qr.nbQuadPt() )
    {
        setupInternalConstants();
    };


    //! Constructor without quadrature rule
    /*!
      @param refFE The reference element for the FE
      @param geoMap The geometric map from the reference element to the current element
     */
    ETCurrentBDFE (const GeometricMap& geoMap)
        : M_geometricMap (&geoMap),
          M_nbMapDof (geoMap.nbDof() )
    {};

    //! Copy constructor
    /*!
      @param otherFE The currentFE to be copied
     */
    ETCurrentBDFE (const ETCurrentBDFE<spaceDim>& otherFE)
        : M_geometricMap (otherFE.M_geometricMap),
          M_quadrature ( new QuadratureRule (*otherFE.M_quadrature) ),

          M_nbMapDof ( otherFE.M_nbMapDof),
          M_nbQuadPt ( otherFE.M_nbQuadPt),

          M_currentID (otherFE.M_currentID),
          M_currentLocalID (otherFE.M_currentLocalID),

          M_phiMap (otherFE.M_phiMap),
          M_dphiGeometricMap (otherFE.M_dphiGeometricMap),
          M_refTangents ( otherFE.M_refTangents),

          M_cellNode (otherFE.M_cellNode),
          M_quadNode (otherFE.M_quadNode),
          M_jacobian (otherFE.M_jacobian),
          M_tangents (otherFE.M_tangents),

          M_metricTensor (otherFE.M_metricTensor),
          M_measure (otherFE.M_measure),

          M_wMeas (otherFE.M_wMeas),
          M_normal (otherFE.M_normal)
    {};


    //! Destructor
    virtual ~ETCurrentBDFE()
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
    }

    //@}

    void setQuadratureRule (const QuadratureRule& qr)
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
        M_quadrature = new QuadratureRule (qr);
        M_nbQuadPt = M_quadrature->nbQuadPt();

        setupInternalConstants();
    }

    void setRefTangents ( std::vector< VectorSmall< spaceDim> > tangents)
    {
        M_refTangents = tangents;
    }

    template< typename ElementType>
    void update ( const ElementType& cell)
    {
        M_currentID = cell.id();
        M_currentLocalID = cell.localId();

        // Cell nodes
        for (UInt i (0); i < M_nbMapDof; ++i)
        {
            for (UInt j (0); j < spaceDim; ++j)
            {
                M_cellNode[i][j] = cell.point (i).coordinate (j);
            }
        }

        for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        {
            for (UInt iDim (0); iDim < spaceDim; ++iDim)
            {
                M_quadNode[iq][iDim] = 0.0;

                for (UInt iDof (0); iDof < M_nbMapDof; ++iDof)
                {
                    M_quadNode[iq][iDim] += M_cellNode[iDof][iDim] * M_phiMap[iq][iDof];
                }
            }
        }

        // Update the jacobian
        Real partialSum (0.0);
        for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        {
            for (UInt iDim (0); iDim < spaceDim; ++iDim)
            {
                for (UInt jDim (0); jDim < spaceDim; ++jDim)
                {
                    partialSum = 0.0;
                    for (UInt iMap (0); iMap < M_nbMapDof; ++iMap)
                    {
                        partialSum += M_cellNode[iMap][iDim] * M_dphiGeometricMap[iq][iMap][jDim];
                    }
                    M_jacobian[iq][iDim][jDim] = partialSum;
                }
            }
        }

        // Update the tangents
        for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        {
            for (UInt iDim (0); iDim < spaceDim - 1; ++iDim)
            {
                for (UInt jDim (0); jDim < spaceDim; ++jDim)
                {
                    partialSum = 0.0;
                    for (UInt kDim (0); kDim < spaceDim; ++kDim)
                    {
                        partialSum += M_jacobian[iq][jDim][kDim] * M_refTangents[iDim][kDim]; // jdim<->kdim
                    }
                    M_tangents[iq][iDim][jDim] = partialSum;
                }
            }
        }

        for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        {
            for (UInt iDim (0); iDim < spaceDim - 1; ++iDim)
            {
                for (UInt jDim (0); jDim < spaceDim - 1; ++jDim)
                {
                    partialSum = 0.0;
                    for (UInt kDim (0); kDim < spaceDim; ++kDim)
                    {
                        partialSum += M_tangents[iq][iDim][kDim] * M_tangents[iq][jDim][kDim];
                    }
                    M_metricTensor[iq][iDim][jDim] = partialSum;
                }
            }
        }

        // SPECIALIZED FOR THE 3D
        for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        {
            M_measure[iq] = std::sqrt ( M_metricTensor[iq][0][0] * M_metricTensor[iq][1][1]
                                        - M_metricTensor[iq][0][1] * M_metricTensor[iq][1][0]);
        }

        for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        {
            M_wMeas[iq] = M_measure[iq] * M_quadrature->weight (iq);
        }

        // SPECIALIZED FOR THE 3D
        Real n0 (0.0);
        Real n1 (0.0);
        Real n2 (0.0);
        Real nnorm (0.0);
        for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        {
            n0 = M_tangents[iq][0][1] * M_tangents[iq][1][2]
                 - M_tangents[iq][0][2] * M_tangents[iq][1][1];
            n1 = M_tangents[iq][0][2] * M_tangents[iq][1][0]
                 - M_tangents[iq][0][0] * M_tangents[iq][1][2];
            n2 = M_tangents[iq][0][0] * M_tangents[iq][1][1]
                 - M_tangents[iq][0][1] * M_tangents[iq][1][0];
            nnorm = std::sqrt (n0 * n0 + n1 * n1 + n2 * n2);

            M_normal[iq][0] = -n0 / nnorm;
            M_normal[iq][1] = -n1 / nnorm;
            M_normal[iq][2] = -n2 / nnorm;
        }
    }

    //private:

    void setupInternalConstants()
    {
        M_phiMap.resize (M_nbQuadPt);
        for (UInt q (0); q < M_nbQuadPt; ++q)
        {
            M_phiMap[q].resize (M_nbMapDof);
            for (UInt i (0); i < M_nbMapDof; ++i)
            {
                M_phiMap[q][i] = M_geometricMap->phi (i, M_quadrature->quadPointCoor (q) );
            }
        }


        M_dphiGeometricMap.resize (M_nbQuadPt);
        for (UInt q (0); q < M_nbQuadPt; ++q)
        {
            M_dphiGeometricMap[q].resize (M_nbMapDof);
            for (UInt i (0); i < M_nbMapDof; ++i)
            {
                M_dphiGeometricMap[q][i].resize (spaceDim);
                for (UInt j (0); j < spaceDim; ++j)
                {
                    M_dphiGeometricMap[q][i][j] = M_geometricMap->dPhi (i, j, M_quadrature->quadPointCoor (q) );
                }
            }
        }

        M_cellNode.resize (M_nbMapDof);
        for (UInt i (0); i < M_nbMapDof; ++i)
        {
            M_cellNode[i].resize (spaceDim);
        }

        M_quadNode.resize (M_nbQuadPt);

        M_jacobian.resize (M_nbQuadPt);
        for (UInt iq (0); iq < M_nbQuadPt; ++iq)
        {
            M_jacobian[iq].resize (spaceDim);
            for (UInt iDim (0); iDim < spaceDim; ++iDim)
            {
                M_jacobian[iq][iDim].resize (spaceDim);
            }
        }

        M_refTangents.resize (spaceDim - 1);

        M_tangents.resize (M_nbQuadPt);
        for (UInt i (0); i < M_nbQuadPt; ++i)
        {
            M_tangents[i].resize (spaceDim - 1);
        }

        M_metricTensor.resize (M_nbQuadPt);
        for (UInt i (0); i < M_nbQuadPt; ++i)
        {
            M_metricTensor[i].resize (spaceDim - 1);
            for (UInt j (0); j < spaceDim - 1; ++j)
            {
                M_metricTensor[i][j].resize (spaceDim - 1);
            }
        }

        M_measure.resize (M_nbQuadPt);

        M_wMeas.resize (M_nbQuadPt);

        M_normal.resize (M_nbQuadPt);

    }


    const GeometricMap* M_geometricMap;
    const QuadratureRule* M_quadrature;

    UInt M_nbMapDof;
    UInt M_nbQuadPt;


    UInt M_currentID;
    UInt M_currentLocalID;

    std::vector< std::vector < Real > > M_phiMap;
    std::vector< std::vector < std::vector < Real> > > M_dphiGeometricMap;
    std::vector < VectorSmall<spaceDim> > M_refTangents; // ASSUME CONSTANT TANGENTS


    std::vector< std::vector < Real > > M_cellNode;
    std::vector< VectorSmall<spaceDim> > M_quadNode;
    std::vector< std::vector < std::vector < Real > > > M_jacobian;
    std::vector< std::vector < VectorSmall<spaceDim> > > M_tangents;
    std::vector< std::vector < std::vector <Real> > > M_metricTensor;
    std::vector< Real > M_measure;
    std::vector< Real > M_wMeas;

    std::vector< VectorSmall<spaceDim> > M_normal;

};


} // Namespace LifeV

#endif /* ETCURRENTBDFE_H */
