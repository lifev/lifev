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
    @brief Geometric elements

    @author Luca Formaggia <luca.formaggia@polimi.it>
    @contributor Marta D'Elia <mdelia2@mathcs.emory.edu>
    @maintainer Marta D'Elia <mdelia2@mathcs.emory.edu>

    @date 00-00-0000

 */

#ifndef _GEOELEMENT_HH_
#define _GEOELEMENT_HH_

#include <life/lifecore/life.hpp>
#include <life/lifemesh/markers.hpp>
#include <life/lifemesh/geoND.hpp>
#include <life/lifemesh/bareItems.hpp>

namespace LifeV
{

//! GeoElement0D - Class for Points and Vertices
/*!
    @author Luca Formaggia

    @sa markers.h
 */
template <typename MC = DefMarkerCommon>
class GeoElement0D: public Geo0D, public MC::PointMarker
{
public:

    //! @name Public Types
    //@{
    typedef typename MC::PointMarker Marker; //to be removed
    typedef typename MC::PointMarker marker_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    GeoElement0D();

    //! Declares item identity and states if it is on boundary
    /*!
    	@param identity Element identity
        @param boundary True if the element is on boundary
     */
    GeoElement0D( ID identity, bool boundary = false );

    //! Declares item identity, provides coordinate and states if it is on boundary
    /*!
    	@param identity Element identity
    	@param x Element x coordinate
    	@param y Element y coordinate
    	@param z Element z coordinate
    	@param boundary True if the element is on boundary
     */
    GeoElement0D( ID identity, Real x, Real y, Real z, bool boundary = false );

    //! Copy constructor
    /*!
        @param Element GeoElement0D to be copied
     */
    GeoElement0D( GeoElement0D const & Element );

    //! Copy constructor
    /*!
        @param Element GeoElement0D to be copied
        @param Marker Markercommon
     */
    GeoElement0D( Geo0D const & Element, MC const & Marker );

    //@}

    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
        @param Element Equivalent GeoElement0D
        @return Reference to a new GeoElement0D with the same content of GeoElement0D Element
     */
    GeoElement0D & operator = ( const GeoElement0D  & Element );

    //@}

};

//! GeoElement1D - Class for Edges
/*!
    @author Luca Formaggia
	@warning In the 2D case, Identities of the adjacent 2Delements and their relative position are stored
 */
template
<typename GEOSHAPE, typename MC = DefMarkerCommon>
class GeoElement1D : public GeoND<GEOSHAPE, GeoElement0D<MC> >, public MC::EdgeMarker
{
public:

    //! @name Public Types
    //@{

    typedef GEOSHAPE geoShape_Type;
    typedef typename MC::EdgeMarker marker_Type;
    typedef GeoElement0D<MC> geoBElement_Type;
    typedef GeoElement0D<MC> point_Type;
    typedef GEOSHAPE GeoShape; // to be removed
    typedef typename MC::EdgeMarker Marker; // to be removed
    typedef GeoElement0D<MC> GeoBElement; // to be removed
    typedef GeoElement0D<MC> PointType; // to be removed
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit GeoElement1D( ID identity = 0 );

    //! Copy constructor
    /*!
        @param Element GeoElement1D to be copied
     */
    GeoElement1D( const GeoElement1D& Element);

    //@}

#ifdef TWODIM
    //! @name Get Methods
    //@{

    //! Returns the identity of the first adjacent element
    /*!
    	@return Identity of the first adjacent element
     */
    ID ad_first() const
    {
        return M_FirstAdjacentElementIdentity;
    };
    //! Returns the identity of the second adjacent element
    /*!
    	@return Identity of the second adjacent element
     */
    ID ad_second() const
    {
        return M_SecondAdjacentElementIdentity;
    };
    //! Returns the identity of the first adjacent element
    /*!
    	@return Identity of the first adjacent element
     */
    ID & ad_first()
    {
        return M_FirstAdjacentElementIdentity;
    };
    //! Returns the identity of the second adjacent element
    /*!
    	@return Identity of the second adjacent element
     */
    ID & ad_second()
    {
        return M_SecondAdjacentElementIdentity;
    };
    //! Returns the position of the first adjacent element
    /*!
    	@return Position of the first adjacent element
     */
    ID pos_first() const
    {
        return M_FirstAdjacentElementPosition;
    };
    //! Returns the position of the second adjacent element
    /*!
    	@return Position of the second adjacent element
     */
    ID pos_second() const
    {
        return M_SecondAdjacentElementPosition;
    };
    //! Returns the position of the first adjacent element
    /*!
    	@return Position of the first adjacent element
     */
    ID & pos_first()
    {
        return M_FirstAdjacentElementPosition;
    };
    //! Returns the position of the second adjacent element
    /*!
    	@return Position of the second adjacent element
     */
    ID & pos_second()
    {
        return M_SecondAdjacentElementPosition;
    };
    //@}

private:
    ID M_FirstAdjacentElementIdentity;
    ID M_SecondAdjacentElementIdentity;
    ID M_FirstAdjacentElementPosition;
    ID M_SecondAdjacentElementPosition;
#endif
};

//! GeoElement2D - Class for Faces
/*!
    @author Luca Formaggia
	@warning In the 3D case, Identities of the adjacent 3D elements and their relative position are stored
 */
template
<typename GEOSHAPE, typename MC = DefMarkerCommon>
class GeoElement2D: public GeoND<GEOSHAPE, GeoElement0D<MC> >, public MC::FaceMarker
{

public:

    //! @name Public Types
    //@{

    //! Number of element edges, for compatibility
    static const UInt numLocalEdges = GeoND<GEOSHAPE, GeoElement0D<MC> >::numEdges;

    typedef GEOSHAPE geoShape_Type;
    typedef typename MC::FaceMarker marker_Type;
    typedef typename GEOSHAPE::GeoBShape edgeShape_Type;
    typedef GeoElement1D<edgeShape_Type, MC> edge_Type;
    typedef GeoElement0D<MC> point_Type;
    typedef edge_Type geoBElement_Type;

    typedef GEOSHAPE GeoShape; // to be removed
    typedef typename MC::FaceMarker Marker; // to be removed
    typedef typename GEOSHAPE::GeoBShape EdgeShape; // to be removed
    typedef GeoElement1D<EdgeShape, MC> EdgeType; // to be removed
    typedef GeoElement0D<MC> PointType; // to be removed
    typedef EdgeType GeoBElement; // to be removed
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit GeoElement2D( ID identity = 0 );

    //! Copy constructor
    /*!
        @param Element GeoElement2D to be copied
     */
    GeoElement2D( const GeoElement2D<GEOSHAPE, MC>& Element);

    //@}

    //! @name Get Methods
    //@{

    //! Returns the identity of the first adjacent element
    /*!
    	@return Identity of the first adjacent element
     */
    ID ad_first() const
    {
        return M_FirstAdjacentElementIdentity;
    };
    //! Returns the identity of the second adjacent element
    /*!
    	@return Identity of the second adjacent element
     */
    ID ad_second() const
    {
        return M_SecondAdjacentElementIdentity;
    };
    //! Returns the identity of the first adjacent element
    /*!
    	@return Identity of the first adjacent element
     */
    ID & ad_first()
    {
        return M_FirstAdjacentElementIdentity;
    };
    //! Returns the identity of the second adjacent element
    /*!
    	@return Identity of the second adjacent element
     */
    ID & ad_second()
    {
        return M_SecondAdjacentElementIdentity;
    };
    //! Returns the position of the first adjacent element
    /*!
    	@return Position of the first adjacent element
     */
    ID pos_first() const
    {
        return M_FirstAdjacentElementPosition;
    };
    //! Returns the position of the second adjacent element
    /*!
    	@return Position of the second adjacent element
     */
    ID pos_second() const
    {
        return M_SecondAdjacentElementPosition;
    };
    //! Returns the position of the first adjacent element
    /*!
    	@return Position of the first adjacent element
     */
    ID & pos_first()
    {
        return M_FirstAdjacentElementPosition;
    };
    //! Returns the position of the second adjacent element
    /*!
    	@return Position of the second adjacent element
     */
    ID & pos_second()
    {
        return M_SecondAdjacentElementPosition;
    };
    //@}

private:
    ID M_FirstAdjacentElementIdentity;
    ID M_SecondAdjacentElementIdentity;
    ID M_FirstAdjacentElementPosition;
    ID M_SecondAdjacentElementPosition;
};


//! GeoElement3D - Class for Volumes
/*!
    @author Luca Formaggia
 */
template
<typename GEOSHAPE, typename MC = DefMarkerCommon>
class GeoElement3D: public GeoND<GEOSHAPE, GeoElement0D<MC> >, public MC::VolumeMarker
{
public:

    //! @name Public Types
    //@{

    //! Number of local Vertices
    static const UInt numLocalVertices = GeoND<GEOSHAPE, GeoElement0D<MC> >::numVertices;
    //! Number of local Faces
    static const UInt numLocalFaces = GeoND<GEOSHAPE, GeoElement0D<MC> >::numFaces;
    //! Number of local Edges (using Euler Formula)
    static const UInt numLocalEdges = GeoND<GEOSHAPE, GeoElement0D<MC> >::numEdges;

    typedef GEOSHAPE geoShape_Type;
    typedef typename MC::VolumeMarker marker_Type;
    typedef typename GEOSHAPE::GeoBShape faceShape_Type;
    typedef typename faceShape_Type::GeoBShape edgeShape_Type;

    typedef GeoElement1D<edgeShape_Type, MC> edge_Type;
    typedef GeoElement2D<faceShape_Type, MC> face_Type;
    typedef GeoElement0D<MC> point_Type;
    typedef face_Type geoBElement_Type;

    typedef GEOSHAPE GeoShape;
    typedef typename MC::VolumeMarker Marker;
    typedef typename GEOSHAPE::GeoBShape FaceShape;
    typedef typename FaceShape::GeoBShape EdgeShape;

    typedef GeoElement1D<EdgeShape, MC> EdgeType; // to be removed
    typedef GeoElement2D<FaceShape, MC> FaceType; // to be removed
    typedef GeoElement0D<MC> PointType; // to be removed
    typedef FaceType GeoBElement; // to be removed

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit GeoElement3D( ID identity = 0 );

    //! Copy constructor
    /*!
        @param Element GeoElement3D to be copied
     */
    GeoElement3D( const GeoElement3D<GEOSHAPE, MC>& Element );

    //@}
};


/*-------------------------------------------------------------------------
  GeoElement0D
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================
template <typename MC>
GeoElement0D<MC>::GeoElement0D() :
        Geo0D(), MC::PointMarker()
{}

template <typename MC>
GeoElement0D<MC>::GeoElement0D( ID identity, bool boundary ) :
        Geo0D( identity, boundary ), MC::PointMarker()
{}

template <typename MC>
GeoElement0D<MC>::GeoElement0D( ID identity, Real x, Real y, Real z, bool boundary ) :
        Geo0D( identity, x, y, z, boundary ), MC::PointMarker()
{}

template <typename MC>
GeoElement0D<MC>::GeoElement0D( GeoElement0D<MC> const & Element ) :
        Geo0D( Element ), MC::PointMarker( Element )
{}

// ==========================================
// Operators
// ==========================================
//! It calls operator= of base classes, just to be sure to do the right thing.
template <typename MC>
GeoElement0D<MC> &
GeoElement0D<MC>::operator = ( GeoElement0D<MC> const & Element )
{
    if ( this != &Element )
    {
        Geo0D::operator=( Element );
        Marker::operator=( Element );
    }
    return *this;
}

/*-------------------------------------------------------------------------
  GeoElement1D
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================
#ifdef TWODIM
template <typename GEOSHAPE, typename MC>
GeoElement1D<GEOSHAPE, MC>::GeoElement1D( ID Identity ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( Identity ),
        MC::EdgeMarker (),
        M_FirstAdjacentElementIdentity    ( 0 ),
        M_SecondAdjacentElementIdentity   ( 0 ),
        M_FirstAdjacentElementPosition    ( 0 ),
        M_SecondAdjacentElementPosition   ( 0 )
#else
template <typename GEOSHAPE, typename MC>
GeoElement1D<GEOSHAPE, MC>::GeoElement1D( ID identity ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( identity )
#endif
{
    ASSERT_PRE( GEOSHAPE::nDim == 1 , "geoElement2D with incorrect GeoShape" ) ;
}

#ifdef TWODIM
template <typename GEOSHAPE, typename MC>
GeoElement1D<GEOSHAPE, MC>::GeoElement1D( const GeoElement1D<GEOSHAPE, MC>& Element ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( Element ),
        MC::EdgeMarker                    ( Element ),
        M_FirstAdjacentElementIdentity    ( Element.M_FirstAdjacentElementIdentity),
        M_SecondAdjacentElementIdentity   ( Element.M_SecondAdjacentElementIdentity),
        M_FirstAdjacentElementPosition    ( Element.M_FirstAdjacentElementPosition ),
        M_SecondAdjacentElementPosition   ( Element.M_SecondAdjacentElementPosition )
#else
template <typename GEOSHAPE, typename MC>
GeoElement1D<GEOSHAPE, MC>::GeoElement1D( const GeoElement1D<GEOSHAPE, MC>& geoElement1D ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( geoElement1D ),
        MC::EdgeMarker                    ( geoElement1D )
#endif
{
    ASSERT_PRE( GEOSHAPE::nDim == 1 , "geoElement2D with incorrect GeoShape" ) ;
}



/*-------------------------------------------------------------------------
  GeoElement2D
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================
template <typename GEOSHAPE, typename MC>
const UInt GeoElement2D<GEOSHAPE, MC>::numLocalEdges;

#ifdef TWODIM
template <typename GEOSHAPE, typename MC>
GeoElement2D<GEOSHAPE, MC>::GeoElement2D( ID identity ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( identity )
#else
template <typename GEOSHAPE, typename MC>
GeoElement2D<GEOSHAPE, MC>::GeoElement2D( ID identity ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( identity ),
        M_FirstAdjacentElementIdentity   ( 0 ),
        M_SecondAdjacentElementIdentity  ( 0 ),
        M_FirstAdjacentElementPosition   ( 0 ),
        M_SecondAdjacentElementPosition  ( 0 )
#endif
{
    ASSERT_PRE( GEOSHAPE::nDim == 2 , "geoElement2D with incorrect GeoShape" ) ;
}

#ifdef TWODIM
template <typename GEOSHAPE, typename MC>
GeoElement2D<GEOSHAPE, MC>::GeoElement2D( const GeoElement2D<GEOSHAPE, MC>& Element ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( Element ),
        MC::FaceMarker                    ( Element )
#else
template <typename GEOSHAPE, typename MC>
GeoElement2D<GEOSHAPE, MC>::GeoElement2D( const GeoElement2D<GEOSHAPE, MC>& Element ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( Element ),
        MC::FaceMarker                    ( Element ),
        M_FirstAdjacentElementIdentity    ( Element.M_FirstAdjacentElementIdentity),
        M_SecondAdjacentElementIdentity   ( Element.M_SecondAdjacentElementIdentity),
        M_FirstAdjacentElementPosition    ( Element.M_FirstAdjacentElementPosition ),
        M_SecondAdjacentElementPosition   ( Element.M_SecondAdjacentElementPosition )
#endif
{
    ASSERT_PRE( GEOSHAPE::nDim == 2 , "geoElement2D with incorrect GeoShape" ) ;
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

// ==========================================
// Constructor & Destructor
// ==========================================
template <typename GEOSHAPE, typename MC>
GeoElement3D<GEOSHAPE, MC>::GeoElement3D( ID identity ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( identity )
{
    ASSERT_PRE( GEOSHAPE::nDim == 3 , "geoElement3D with incorrect GeoShape" )
}

template <typename GEOSHAPE, typename MC>
GeoElement3D<GEOSHAPE, MC>::GeoElement3D( const GeoElement3D<GEOSHAPE, MC>& Element ) :
        GeoND<GEOSHAPE, GeoElement0D<MC> >( Element ),
        MC::VolumeMarker                  ( Element )
{
    ASSERT_PRE( GEOSHAPE::nDim == 3 , "geoElement3D with incorrect GeoShape" )
}
}
#endif

