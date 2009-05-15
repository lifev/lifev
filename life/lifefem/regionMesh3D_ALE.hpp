/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*! file regionMesh3D_ALE.h
  \brief Extension of the  mesh classes interfaces for ALE schemes

  \version $Revision: 1.7 $ Miguel Fernandez

  Introduces the RegionMesh3D class
*/

#ifndef _REGIONMESH3DALE_HH_
#define _REGIONMESH3DALE_HH_
#include <life/lifecore/life.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifemesh/basisElSh.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/quadRule.hpp>

namespace LifeV
{
template <typename GEOSHAPE, typename MC = DefMarkerCommon >
class RegionMesh3D_ALE
            :
            public RegionMesh3D<GEOSHAPE, MC>
{
public:
    typedef typename RegionMesh3D<GEOSHAPE, MC>::ElementShape ElementShape;

    explicit RegionMesh3D_ALE();
    //! Constructor providing ID
    explicit RegionMesh3D_ALE( ID id );
    //! Copy contructor: not implemented yet!
    explicit RegionMesh3D_ALE( RegionMesh3D_ALE<GEOSHAPE, MC> const & m );

    RegionMesh3D_ALE<GEOSHAPE, MC> operator=( RegionMesh3D_ALE<GEOSHAPE, MC> const & m );

    ~RegionMesh3D_ALE <GEOSHAPE, MC>();
    //@{
    /*! Get the reference RefFE object associated to this mesh.
      It is necessary when implementing mesh movement routines based on
      harmonic reconstruction, in order to build the finite element discretisation.
    */
    const RefFE& getRefFE() const;
    /*! Get the reference GeoMap object associated to this mesh.
      It is necessary when implementing mesh movement routines based on
      harmonic reconstruction, in order to build the finite element discretisation.
    */
    const GeoMap& getGeoMap() const;
};


/* ---------------------------------------------------------------------
   RegionMesh3D_ALE Implementation
   -----------------------------------------------------------------------*/
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
RegionMesh3D_ALE<GEOSHAPE, MC>
RegionMesh3D_ALE<GEOSHAPE, MC>::operator=( RegionMesh3D_ALE<GEOSHAPE, MC> const & m )
{
    RegionMesh3D<GEOSHAPE, MC>::operator=( m );
    return *this;
}

template <typename GEOSHAPE, typename MC>
RegionMesh3D_ALE<GEOSHAPE, MC>::~RegionMesh3D_ALE()
{}

//! Modif Miguel:11/2002
//!< Get the reference RefFE object associated to the mes
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
}

//!< Get the reference GeoMap object  associated to the mesh
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
}
}
#endif
