/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
/*! file geoElement.h
\brief Geometric elements
\version $Revision: 1.10 $ Luca Formaggia

  Introduces all the geometric elements
*/

#ifndef _GEOELEMENT_HH_
#define _GEOELEMENT_HH_

#include <life/lifecore/life.hpp>
#include <life/lifemesh/markers.hpp>
#include <life/lifemesh/geoND.hpp>
#include <life/lifemesh/bareItems.hpp>

namespace LifeV
{
//using namespace std; //Pretty useless

//     *********** Geometrical Elements *****************
//! \defgroup GeoEle Geometry Element classes

/*! \warning It is UP to the user to parametrise the templates correctly. */

/*@{*/
//! Class for Points and Vertices
/*! MC is a markercommon which defines the markers for
  all geometric entities. */
template <typename MC = DefMarkerCommon>
class
            GeoElement0D: public Geo0D, public MC::PointMarker
{
public:
    GeoElement0D();

    //! Declares item id and if it is on boundary
    GeoElement0D( ID id, bool boundary = false );

    //! Declares item id and if it is on boundary, and provides coordinate
    //! data.
    GeoElement0D( ID id, Real x, Real y, Real z, bool boundary = false );

    GeoElement0D( GeoElement0D const & g );
    GeoElement0D( Geo0D const & g, MC const & m );
    GeoElement0D & operator = ( GeoElement0D const & g );

    typedef typename MC::PointMarker Marker;

};


//!  Class for edges.
/*! In the 2D case, we store the IDs of the adjacent
  2Delements and their relative position.*/
template
<typename GEOSHAPE, typename MC = DefMarkerCommon>
class GeoElement1D : public GeoND<GEOSHAPE, GeoElement0D<MC> >, public MC::EdgeMarker
{
public:
    GeoElement1D( ID id = 0 );

    typedef GEOSHAPE GeoShape;
    typedef typename MC::EdgeMarker Marker;
    typedef GeoElement0D<MC> GeoBElement;
    typedef GeoElement0D<MC> PointType;

#ifdef TWODIM
    //! ID of first adjacent 2D element
    ID ad_first() const
    {
        return efirst;
    };
    //! ID of second adjacent 2D element
    ID ad_second() const
    {
        return esecond;
    };
    //! ID of first adjacent 2D element
    ID & ad_first()
    {
        return efirst;
    }
    //! ID of second adjacent 2D element
    ID & ad_second()
    {
        return esecond;
    }
    //! Position of first adjacent 2D element
    ID pos_first() const
    {
        return posfirst;
    };
    //! Position of second adjacent 2D element
    ID pos_second() const
    {
        return possecond;
    };
    //! Position of fist adjacent 2D element
    ID & pos_first()
    {
        return posfirst;
    }
    //! Position of second adjacent 2D element
    ID & pos_second()
    {
        return possecond;
    }
private:
    ID efirst;
    ID esecond;
    ID posfirst;
    ID possecond;
#endif
};

//!  Class for 2D elements.

/*! In the 3D case, we store the IDs of the adjacent
  3Delements and their relative position.*/
template
<typename GEOSHAPE, typename MC = DefMarkerCommon>
class GeoElement2D
            :
            public GeoND<GEOSHAPE, GeoElement0D<MC> >,
            public MC::FaceMarker
{
public:

    GeoElement2D( ID id = 0 );
    //! Number of element edges
    static const UInt numLocalEdges = GeoND<GEOSHAPE, GeoElement0D<MC> >::numEdges; //For compatibility

    typedef GEOSHAPE GeoShape;
    typedef typename MC::FaceMarker Marker;
    typedef typename GEOSHAPE::GeoBShape EdgeShape;
    typedef GeoElement1D<EdgeShape, MC> EdgeType;
    typedef GeoElement0D<MC> PointType;
    typedef EdgeType GeoBElement;


#ifdef THREEDIM
    //! ID of first adjacent 2D element
    ID ad_first() const
    {
        return efirst;
    };
    //! ID of second adjacent 2D element
    ID ad_second() const
    {
        return esecond;
    };
    //! ID of first adjacent 2D element
    ID & ad_first()
    {
        return efirst;
    }
    //! ID of second adjacent 2D element
    ID & ad_second()
    {
        return esecond;
    }
    //! Position of first adjacent 2D element
    ID pos_first() const
    {
        return posfirst;
    };
    //! Position of second adjacent 2D element
    ID pos_second() const
    {
        return possecond;
    };
    //! Position of fist adjacent 2D element
    ID & pos_first()
    {
        return posfirst;
    }
    //! Position of second adjacent 2D element
    ID & pos_second()
    {
        return possecond;
    }
private:
    ID efirst;
    ID esecond;
    ID posfirst;
    ID possecond;
#endif
};


//! Class for 3D elements
template
<typename GEOSHAPE, typename MC = DefMarkerCommon>
class GeoElement3D
            :
            public GeoND<GEOSHAPE, GeoElement0D<MC> >,
            public MC::VolumeMarker
{
public:

    GeoElement3D( ID id = 0 );

    typedef GEOSHAPE GeoShape;
    typedef typename MC::VolumeMarker Marker;
    typedef typename GEOSHAPE::GeoBShape FaceShape;
    typedef typename FaceShape::GeoBShape EdgeShape;

    typedef GeoElement1D<EdgeShape, MC> EdgeType;
    typedef GeoElement2D<FaceShape, MC> FaceType;
    typedef GeoElement0D<MC> PointType;
    typedef FaceType GeoBElement;

    //! Number of local Vertices
    static const UInt numLocalVertices = GeoND<GEOSHAPE, GeoElement0D<MC> >::numVertices;
    //! Number of local Faces
    static const UInt numLocalFaces = GeoND<GEOSHAPE, GeoElement0D<MC> >::numFaces;
    //! Number of local Edges (using Euler Formula)
    static const UInt numLocalEdges = GeoND<GEOSHAPE, GeoElement0D<MC> >::numEdges;
};
/*@}*/


/*-------------------------------------------------------------------------
  GeoElement0D
  --------------------------------------------------------------------------*/
template <typename MC>
GeoElement0D<MC>::GeoElement0D() :
        Geo0D(), MC::PointMarker()
{}

template <typename MC>
GeoElement0D<MC>::GeoElement0D( ID id, bool boundary ) :
        Geo0D( id, boundary ), MC::PointMarker()
{}

template <typename MC>
GeoElement0D<MC>::GeoElement0D( ID id, Real x, Real y, Real z, bool boundary ) :
        Geo0D( id, x, y, z, boundary ), MC::PointMarker()
{}

template <typename MC>
GeoElement0D<MC>::GeoElement0D( GeoElement0D<MC> const & g ) :
        Geo0D( g ), MC::PointMarker( g )
{}

//! It calls operator= of base classes, just to be sure to do the right thing.
template <typename MC>
GeoElement0D<MC> &
GeoElement0D<MC>::operator = ( GeoElement0D<MC> const & g )
{
    if ( this != &g )
    {
        Geo0D::operator=( g );
        Marker::operator=( g );
    }
    return *this;
}

/*-------------------------------------------------------------------------
  GeoElement1D
  --------------------------------------------------------------------------*/

#ifdef TWODIM
template <typename GEOSHAPE, typename MC>
GeoElement1D<GEOSHAPE, MC>::GeoElement1D( ID id ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( id ),
        efirst( 0 ),
        esecond( 0 ),
        posfirst( 0 ),
        possecond( 0 )
#else
template <typename GEOSHAPE, typename MC>
GeoElement1D<GEOSHAPE, MC>::GeoElement1D( ID id ) : GeoND<GEOSHAPE, GeoElement0D<MC> >( id )
#endif
{
    ASSERT_PRE( GEOSHAPE::nDim == 1 , "geoElement2D with incorrect GeoSHape" ) ;
}

/*-------------------------------------------------------------------------
  GeoElement2D
  --------------------------------------------------------------------------*/
template <typename GEOSHAPE, typename MC>
const UInt GeoElement2D<GEOSHAPE, MC>::numLocalEdges;

#ifdef THREEDIM
template <typename GEOSHAPE, typename MC>
GeoElement2D<GEOSHAPE, MC>::GeoElement2D( ID id ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( id ),
        efirst( 0 ),
        esecond( 0 ),
        posfirst( 0 ),
        possecond( 0 )
#else
template <typename GEOSHAPE, typename MC>
GeoElement2D<GEOSHAPE, MC>::GeoElement2D( ID id ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( id )
#endif
{
    ASSERT_PRE( GEOSHAPE::nDim == 2 , "geoElement2D with incorrect GeoSHape" ) ;
}


/*-------------------------------------------------------------------------
                 GeoElement3D
 --------------------------------------------------------------------------*/
template <typename GEOSHAPE, typename MC>
const UInt GeoElement3D<GEOSHAPE, MC>::numLocalVertices;
template <typename GEOSHAPE, typename MC>
const UInt GeoElement3D<GEOSHAPE, MC>::numLocalFaces;
template <typename GEOSHAPE, typename MC>
const UInt GeoElement3D<GEOSHAPE, MC>::numLocalEdges;

template <typename GEOSHAPE, typename MC>
GeoElement3D<GEOSHAPE, MC>::GeoElement3D( ID id ) : GeoND<GEOSHAPE, GeoElement0D<MC> >( id )
{
    ASSERT_PRE( GEOSHAPE::nDim == 3 , "geoElement3D with incorrect GeoSHape" )
}
}
#endif

