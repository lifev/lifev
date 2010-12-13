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
    @brief Extension of the mesh classes interfaces for ALE schemes

    @author Miguel Fernandez Ruiz <miguel.fernandezruiz@epfl.ch>
    @contributor Nur Aiman Fadel <nur.fadel@mail.polimi.it>
    @maintainer Nur Aiman Fadel <nur.fadel@mail.polimi.it>

    @date

    Introduces the RegionMesh3D class for ALE schemes.

    @todo RegionMesh3D_ALE( RegionMesh3D_ALE<GEOSHAPE, MC> const & m ) is not implemented.
 */

#ifndef _REGIONMESH3DALE_HH_
#define _REGIONMESH3DALE_HH_ 1

#include <life/lifemesh/regionMesh3D.hpp>

#include <life/lifefem/refFEScalar.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/quadRule.hpp>

namespace LifeV
{

//! RegionMesh3D_ALE - Extension of the mesh classes interfaces for ALE schemes
/*!
    @author Miguel Fernandez Ruiz
    @see

    This class is an extension of regionMesh3D.hpp in order to use ALE schemes.<br>
 */

template <typename GEOSHAPE, typename MC = DefMarkerCommon >
class RegionMesh3D_ALE:
        public RegionMesh3D<GEOSHAPE, MC>
{
public:

    //! @name Public Types
    //@{

    typedef typename RegionMesh3D<GEOSHAPE, MC>::ElementShape ElementShape;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    explicit RegionMesh3D_ALE();

    //! Constructor providing ID
    /*!
        It is a constructor which requires the identifier.
        @param id, the identifier.
    */
    explicit RegionMesh3D_ALE( ID id );

    //! Copy constructor
    /*!
        It is the copy constructor.
        @param m, the constructor to be copied.
        @todo It is not implemented.
    */
    explicit RegionMesh3D_ALE( RegionMesh3D_ALE<GEOSHAPE, MC> const & m );

    //! Destructor
    ~RegionMesh3D_ALE <GEOSHAPE, MC>();
    //@}

    //! @name Operators
    //@{

    //! operator equivalence
    RegionMesh3D_ALE<GEOSHAPE, MC> operator=( RegionMesh3D_ALE<GEOSHAPE, MC> const & m );
    //@}

    //! @name Get Methods
    //@{

    //! Get the reference RefFE object associated to this mesh.
    /*!
      Get the reference RefFE object associated to this mesh.
      It is necessary when implementing mesh movement routines based on
      harmonic reconstruction, in order to build the finite element discretisation.
    */
    const RefFE& getRefFE() const;

    //! Get the reference GeoMap object  associated to the mesh
    /*!
      Get the reference GeoMap object associated to this mesh.
      It is necessary when implementing mesh movement routines based on
      harmonic reconstruction, in order to build the finite element discretisation.
    */
    const GeoMap& getGeoMap() const;
}; // Class RegionMesh3D_ALE

// ===================================================
// Constructors & Destructor
// ===================================================
template <typename GEOSHAPE, typename MC>
RegionMesh3D_ALE<GEOSHAPE, MC>::RegionMesh3D_ALE() :
        RegionMesh3D<GEOSHAPE, MC>()
{}

template <typename GEOSHAPE, typename MC>
RegionMesh3D_ALE<GEOSHAPE, MC>::RegionMesh3D_ALE( ID id ) :
        RegionMesh3D<GEOSHAPE, MC>( id )
{}

template <typename GEOSHAPE, typename MC>
RegionMesh3D_ALE<GEOSHAPE, MC>::RegionMesh3D_ALE( RegionMesh3D_ALE<GEOSHAPE, MC> const & m ) :
        RegionMesh3D<GEOSHAPE, MC>( m )
{}

template <typename GEOSHAPE, typename MC>
RegionMesh3D_ALE<GEOSHAPE, MC>::~RegionMesh3D_ALE()
{}

// ===================================================
// Operators
// ===================================================

// operator equivalence
template <typename GEOSHAPE, typename MC>
RegionMesh3D_ALE<GEOSHAPE, MC>
RegionMesh3D_ALE<GEOSHAPE, MC>::operator=( RegionMesh3D_ALE<GEOSHAPE, MC> const & m )
{
    RegionMesh3D<GEOSHAPE, MC>::operator=( m );
    return *this;
}

// ===================================================
// Get Methods
// ===================================================

// Get the reference RefFE object associated to the mesh
template <typename GEOSHAPE, typename MC>
const RefFE& RegionMesh3D_ALE<GEOSHAPE, MC>::getRefFE() const
{
    switch ( ElementShape::Shape )
    {
    case HEXA:
        if ( ElementShape::numPoints == 8 )
            return feHexaQ1;
        else
            ERROR_MSG( "Finite Element not implemented for the mesh motion" );
        break;
    case TETRA:
        if ( ElementShape::numPoints == 4 )
            return feTetraP1;
        else
            ERROR_MSG( "Finite Element not implemented for the mesh motion" );
        break;
    default:
        ERROR_MSG( "Finite Element not implemented for the mesh motion" );
    }
} // Method RegionMesh3D_ALE<GEOSHAPE, MC>::getRefFE

// Get the reference GeoMap object  associated to the mesh
template <typename GEOSHAPE, typename MC>
const GeoMap& RegionMesh3D_ALE<GEOSHAPE, MC>::getGeoMap() const
{
    switch ( ElementShape::Shape )
    {
    case HEXA:
        if ( ElementShape::numPoints == 8 )
            return geoBilinearHexa;
        else
            ERROR_MSG( "Finite Element not implemented for the mesh motion" );
        break;
    case TETRA:
        if ( ElementShape::numPoints == 4 )
            return geoLinearTetra;
        else
            ERROR_MSG( "Finite Element not implemented for ALE" );
        break;
    default:
        ERROR_MSG( "Finite Element not implemented for ALE" );
    }
} // Method RegionMesh3D_ALE<GEOSHAPE, MC>::getGeoMap
} // Namespace LifeV
#endif  /* REGIONMESH3D_ALE_H */
