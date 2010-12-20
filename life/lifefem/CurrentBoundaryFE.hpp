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
    @brief CurrentFE living on the sides of the elements

    @author Jean-Frederic Gerbeau
    @date 00-09-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */

#ifndef _CURRENTBDFE_H
#define _CURRENTBDFE_H

#include <life/lifecore/life.hpp>

#include <life/lifefem/GeometricMap.hpp>
#include <life/lifefem/ReferenceFE.hpp>
#include <life/lifefem/staticBdFE.hpp>

namespace LifeV
{
/*!
  \class CurrentBdFE
  \brief The class for a boundary finite element
  \author J.-F. Gerbeau
  \date 09/2002

  This class is used for a current boundary elements, i.e. the boundary of
  a CurrentFE. As for the CurrentFE (and contrarily to StaticBdFE) on must
  update this element with a geometrical element before using it. If only
  the measure is needed - to perform for example surface integrations -
  call updateMeas(...), which by the way compute also the tangent.  If, the
  normal at the integration point is needed, call updateMeasNormal(...).
  See the description of the base class StaticBdFE for further details.
*/

class CurrentBdFE:
        public StaticBdFE
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Constructor with reference FE and geometric mapping
    CurrentBdFE( const RefFE& refFE, const GeoMap& geoMap );

    //! Constructor with reference FE, geometric mapping and quadrature rule
    CurrentBdFE( const RefFE& refFE, const GeoMap& geoMap, const QuadRule& qr );

    //! Destructor
    virtual ~CurrentBdFE();

    //@}


    //! @name Methods
    //@{

    /*!
      Compute only the coordinates of the nodes on the current boundary element
    */
    template <typename GeometricType>
    void update( const GeometricType& geometricEntity );

    /*!
      Compute the arrays meas, weightMeas, tangent
      on the current boundary element
    */
    template <typename GeometricType>
    void updateMeas( const GeometricType& geometricEntity );

    /*!
      Compute the arrays meas, weightMeas, tangent
      and quadrature points on the current boundary element
    */
    template <typename GeometricType>
    void updateMeasQuadPt( const GeometricType& geometricEntity );

    /*!
      Compute the arrays meas, weightMeas, tangent
      and normal on the current boundary element
    */
    template <typename GeometricType>
    void updateMeasNormal( const GeometricType& geometricEntity );

    /*!
      Compute the arrays meas, weightMeas, tangent,
      normal and quadrature points on the current boundary element
    */
    template <typename GeometricType>
    void updateMeasNormalQuadPt( const GeometricType& geometricEntity );

    //@}
};

// ===================================================
// Methods
// ===================================================

template <typename GeometricType>
void
CurrentBdFE::
update( const GeometricType& geometricEntity )
{
#ifdef TEST_PRE
    M_hasMeasure = false;
    M_hasTangent = false;
    M_hasNormal = false;
    M_hasQuadPtCoor = false;
    M_hasFirstDerivative = false;
#endif

    M_currentID = geometricEntity.id();
    // update the definition of the geo points
    for ( UInt i(0); i < M_nbGeoNode; i++ )
    {
        for (UInt icoor(0); icoor < nDimensions; icoor++)
        {
            M_point( i, icoor ) = geometricEntity.point( i + 1 ).coordinatesArray()[icoor];
        }
    }
}

template <typename GeometricType>
void
CurrentBdFE::
updateMeas( const GeometricType& geometricEntity )
{
#ifdef TEST_PRE
    M_hasMeasure = true;
    M_hasTangent = true;
    M_hasNormal = false;
    M_hasQuadPtCoor = false;
    M_hasFirstDerivative = false;
#endif

    M_currentID = geometricEntity.id();
    // update the definition of the geo points

    for ( UInt i = 0; i < M_nbGeoNode; i++ )
    {
        for (UInt icoor=0; icoor<nDimensions; icoor++)
        {
            M_point( i, icoor ) = geometricEntity.point( i + 1 ).coordinatesArray()[icoor];
        }
    }

    // compute the measure
    computeMeasure();
}

template <typename GeometricType>
void
CurrentBdFE::
updateMeasQuadPt( const GeometricType& geometricEntity )
{
#ifdef TEST_PRE
    M_hasMeasure = true;
    M_hasTangent = true;
    M_hasNormal = false;
    M_hasQuadPtCoor = true;
    M_hasFirstDerivative = false;
#endif

    M_currentID = geometricEntity.id();
    // update the definition of the geo points

    for ( UInt i = 0; i < M_nbGeoNode; i++ )
    {
        for (UInt icoor=0; icoor<nDimensions; icoor++)
        {
            M_point( i, icoor ) = geometricEntity.point( i + 1 ).coordinatesArray()[icoor];
        }
    }

    // compute the measure
    computeMeasure();
    // compute the coordinates of the quad points
    computeQuadPointCoordinate();
}

template <typename GeometricType>
void
CurrentBdFE::
updateMeasNormal( const GeometricType& geometricEntity )
{
#ifdef TEST_PRE
    M_hasMeasure = true;
    M_hasTangent = true;
    M_hasNormal = true;
    M_hasQuadPtCoor = false;
    M_hasFirstDerivative = false;
#endif

    M_currentID = geometricEntity.id();
    // update the definition of the geo points

    for ( UInt i = 0; i < M_nbGeoNode; i++ )
    {
        for (UInt icoor=0; icoor<nDimensions; icoor++)
        {
            M_point( i, icoor ) = geometricEntity.point( i + 1 ).coordinatesArray()[icoor];
        }
    }

    // compute the measure and the normal
    computeMeasureNormal();
}

template <typename GeometricType>
void
CurrentBdFE::
updateMeasNormalQuadPt( const GeometricType& geometricEntity )
{
#ifdef TEST_PRE
    M_hasMeasure = true;
    M_hasTangent = true;
    M_hasNormal = true;
    M_hasQuadPtCoor = true;
    M_hasFirstDerivative = false;
#endif

    M_currentID = geometricEntity.id();
    // update the definition of the geo points

    for ( UInt i = 0; i < M_nbGeoNode; i++ )
    {
        for (UInt icoor=0; icoor<nDimensions; icoor++)
        {
            M_point( i, icoor ) = geometricEntity.point( i + 1 ).coordinatesArray()[icoor];
        }
    }

    // compute the measure and the normal
    computeMeasureNormal();

    // compute the coordinates of the quad points
    computeQuadPointCoordinate();
}

}
#endif
