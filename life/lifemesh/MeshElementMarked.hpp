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

#ifndef MESHELEMENTMARKED_H
#define MESHELEMENTMARKED_H

#include <life/lifecore/LifeV.hpp>
#include <life/lifemesh/MarkerDefinitions.hpp>
#include <life/lifemesh/MeshElement.hpp>
#include <life/lifemesh/MeshElementBare.hpp>

namespace LifeV
{

//! MeshElementMarked0D - Class for Points and Vertices
/*!
    @author Luca Formaggia

    @sa markers.h
 */
template <typename MC, int geoDim>
class MeshElementMarked0D: public MeshVertex, public MC::pointMarker_Type
{
public:

    //! @name Public Types
    //@{
    typedef typename MC::pointMarker_Type marker_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MeshElementMarked0D();

    //! Declares item identity and states if it is on boundary
    /*!
    	@param identity Element identity
        @param boundary True if the element is on boundary
     */
    MeshElementMarked0D( ID identity, bool boundary = false );

    //! Declares item identity, provides coordinate and states if it is on boundary
    /*!
    	@param identity Element identity
    	@param x Element x coordinate
    	@param y Element y coordinate
    	@param z Element z coordinate
    	@param boundary True if the element is on boundary
     */
    MeshElementMarked0D( ID identity, Real x, Real y, Real z, bool boundary = false );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked0D to be copied
     */
    MeshElementMarked0D( MeshElementMarked0D const & Element );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked0D to be copied
        @param Marker Markercommon
     */
    MeshElementMarked0D( MeshVertex const & Element, MC const & Marker );

    //! Destructor
    virtual ~MeshElementMarked0D()
    {
        // nothing to be done
    }

    //@}

    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
        @param Element Equivalent MeshElementMarked0D
        @return Reference to a new MeshElementMarked0D with the same content of MeshElementMarked0D Element
     */
    MeshElementMarked0D & operator = ( const MeshElementMarked0D  & Element );

    //@}


};

//! MeshElementMarked0Din2D - Class for Points and Vertices
/*!
    @author Mauro Perego

    @sa markers.h
 */
template <typename MC>
class MeshElementMarked0D<MC, 2>: public MeshVertex, public MC::pointMarker_Type
{
public:

    //! @name Public Types
    //@{
    typedef typename MC::pointMarker_Type marker_Type;

    void setPoint( ID const identity, MeshElementMarked0D<MC ,2> const * point );

    void setPoint( ID const /*identity*/, MeshElementMarked0D<MC ,2> const & point );

    MeshElementMarked0D<MC, 2> const & point ( ID const identity ) const;


    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MeshElementMarked0D();

    //! Declares item identity and states if it is on boundary
    /*!
    	@param identity Element identity
        @param boundary True if the element is on boundary
     */
    MeshElementMarked0D( ID identity, bool boundary = false );

    //! Declares item identity, provides coordinate and states if it is on boundary
    /*!
    	@param identity Element identity
    	@param x Element x coordinate
    	@param y Element y coordinate
    	@param z Element z coordinate
    	@param boundary True if the element is on boundary
     */
    MeshElementMarked0D( ID identity, Real x, Real y, Real z, bool boundary = false );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked0Din2D to be copied
     */
    MeshElementMarked0D( MeshElementMarked0D<MC, 2> const & Element );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked0Din2D to be copied
        @param Marker Markercommon
     */
    MeshElementMarked0D( MeshVertex const & Element, MC const & Marker );

    //! Destructor
    virtual ~MeshElementMarked0D()
    {
        // nothing to be done
    }

    //@}

    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
        @param Element Equivalent MeshElementMarked0Din2D
        @return Reference to a new MeshElementMarked0Din2D with the same content of MeshElementMarked0Din2D Element
     */
    MeshElementMarked0D & operator = ( const MeshElementMarked0D  & Element );

    //@}


};

template
<typename GeoShape, typename MC, int geoDim>
class MeshElementMarked1D :
    public MeshElement<GeoShape, MeshElementMarked0D<MC, geoDim> >,
    public MC::edgeMarker_Type
{};

//! MeshElementMarked1D - Class for Edges
/*!
    @author Luca Formaggia
	@warning In the 2D case, Identities of the adjacent 2Delements and their relative position are stored
 */
template
<typename GeoShape, typename MC>
class MeshElementMarked1D<GeoShape, MC, 3> :
    public MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >,
    public MC::edgeMarker_Type
{
public:

    //! @name Public Types
    //@{

    typedef GeoShape geoShape_Type;
    typedef typename MC::edgeMarker_Type marker_Type;
    typedef MeshElementMarked0D<MC, 3> geoBElement_Type;
    typedef MeshElementMarked0D<MC, 3> point_Type;
   //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked1D( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked1D to be copied
     */
    MeshElementMarked1D( const MeshElementMarked1D& Element);

    //! Destructor
    virtual ~MeshElementMarked1D()
    {
        // nothing to be done
    }
};
    //@}

//! MeshElementMarked2D - Class for Faces
/*!
    @author Luca Formaggia
	@warning In the 3D case, Identities of the adjacent 3D elements and their relative position are stored
 */
template
<typename GeoShape, typename MC, int geoDim>
class MeshElementMarked2D:
    public MeshElement<GeoShape, MeshElementMarked0D<MC, geoDim> >,
    public MC::faceMarker_Type
{};


//! MeshElementMarked2D - Class for Faces
/*!
    @author Luca Formaggia
	@warning In the 3D case, Identities of the adjacent 3D elements and their relative position are stored
 */
template
<typename GeoShape, typename MC>
class MeshElementMarked2D<GeoShape, MC, 3>:
    public MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >,
    public MC::faceMarker_Type
{

public:

    //! @name Public Types
    //@{

    //! Number of element edges, for compatibility
    static const UInt S_numLocalEdges = MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >::S_numEdges;

    typedef GeoShape geoShape_Type;
    typedef typename MC::faceMarker_Type marker_Type;
    typedef typename GeoShape::GeoBShape edgeShape_Type;
    typedef MeshElementMarked1D<edgeShape_Type, MC, 3> edge_Type;
    typedef MeshElementMarked0D<MC, 3> point_Type;
    typedef edge_Type geoBElement_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked2D( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked2D to be copied
     */
    MeshElementMarked2D( const MeshElementMarked2D<GeoShape, MC, 3>& Element);

    //! Destructor
    virtual ~MeshElementMarked2D()
    {
        // nothing to be done
    }

    //@}

    //! @name Get Methods
    //@{

    //! Returns the identity of the first adjacent element
    /*!
    	@return Identity of the first adjacent element
     */
    ID firstAdjacentElementIdentity() const
    {

    	return M_firstAdjacentElementIdentity;
    };

    //! Returns the identity of the second adjacent element
    /*!
    	@return Identity of the second adjacent element
     */
    ID secondAdjacentElementIdentity() const
    {
        return M_secondAdjacentElementIdentity;
    };

    //! Returns the identity of the first adjacent element
    /*!
    	@return Identity of the first adjacent element
     */
    ID & firstAdjacentElementIdentity()
    {
        return M_firstAdjacentElementIdentity;
    };

    //! Returns the identity of the second adjacent element
    /*!
    	@return Identity of the second adjacent element
     */
    ID & secondAdjacentElementIdentity()
    {
        return M_secondAdjacentElementIdentity;
    };

    //! Returns the position of the first adjacent element
    /*!
    	@return Position of the first adjacent element
     */
    ID firstAdjacentElementPosition() const
    {
        return M_firstAdjacentElementPosition;
    };

    //! Returns the position of the second adjacent element
    /*!
    	@return Position of the second adjacent element
     */
    ID secondAdjacentElementPosition() const
    {
        return M_secondAdjacentElementPosition;
    };


    //! Returns the position of the first adjacent element
    /*!
    	@return Position of the first adjacent element
     */
    ID & firstAdjacentElementPosition()
    {
        return M_firstAdjacentElementPosition;
    };

    //! Returns the position of the second adjacent element
    /*!
    	@return Position of the second adjacent element
     */
    ID & secondAdjacentElementPosition()
    {
        return M_secondAdjacentElementPosition;
    };

    //@}

private:
    ID M_firstAdjacentElementIdentity;
    ID M_secondAdjacentElementIdentity;
    ID M_firstAdjacentElementPosition;
    ID M_secondAdjacentElementPosition;
};




//! MeshElementMarked1Din2DGeo - Class for Faces
/*!
    @author Luca Formaggia
	@warning In the 3D case, Identities of the adjacent 3D elements and their relative position are stored
 */
template
<typename GeoShape, typename MC>
class MeshElementMarked1D<GeoShape, MC, 2>: public MeshElement<GeoShape, MeshElementMarked0D<MC,2> >, public MC::faceMarker_Type
{

public:

    //! @name Public Types
    //@{

    //! Number of element edges, for compatibility
    static const UInt S_numLocalVertices = MeshElement<GeoShape, MeshElementMarked0D<MC, 2> >::S_numVertices;

    typedef GeoShape geoShape_Type;
    typedef typename MC::faceMarker_Type marker_Type;
    typedef typename GeoShape::GeoBShape edgeShape_Type;
    typedef MeshElementMarked1D<edgeShape_Type, MC, 2> edge_Type;
    typedef MeshElementMarked0D<MC, 2> point_Type;
    typedef edge_Type geoBElement_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked1D( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked1Din2DGeo to be copied
     */
    MeshElementMarked1D( const MeshElementMarked1D<GeoShape, MC, 2>& Element);

    //! Destructor
    virtual ~MeshElementMarked1D()
    {
        // nothing to be done
    }

    //@}

    //! @name Get Methods
    //@{

    //! Returns the identity of the first adjacent element
    /*!
    	@return Identity of the first adjacent element
     */
    ID firstAdjacentElementIdentity() const
    {

    	return M_firstAdjacentElementIdentity;
    };

    //! Returns the identity of the second adjacent element
    /*!
    	@return Identity of the second adjacent element
     */
    ID secondAdjacentElementIdentity() const
    {
        return M_secondAdjacentElementIdentity;
    };

    //! Returns the identity of the first adjacent element
    /*!
    	@return Identity of the first adjacent element
     */
    ID & firstAdjacentElementIdentity()
    {
        return M_firstAdjacentElementIdentity;
    };

    //! Returns the identity of the second adjacent element
    /*!
    	@return Identity of the second adjacent element
     */
    ID & secondAdjacentElementIdentity()
    {
        return M_secondAdjacentElementIdentity;
    };

    //! Returns the position of the first adjacent element
    /*!
    	@return Position of the first adjacent element
     */
    ID firstAdjacentElementPosition() const
    {
        return M_firstAdjacentElementPosition;
    };

    //! Returns the position of the second adjacent element
    /*!
    	@return Position of the second adjacent element
     */
    ID secondAdjacentElementPosition() const
    {
        return M_secondAdjacentElementPosition;
    };


    //! Returns the position of the first adjacent element
    /*!
    	@return Position of the first adjacent element
     */
    ID & firstAdjacentElementPosition()
    {
        return M_firstAdjacentElementPosition;
    };

    //! Returns the position of the second adjacent element
    /*!
    	@return Position of the second adjacent element
     */
    ID & secondAdjacentElementPosition()
    {
        return M_secondAdjacentElementPosition;
    };

    //@}

private:
    ID M_firstAdjacentElementIdentity;
    ID M_secondAdjacentElementIdentity;
    ID M_firstAdjacentElementPosition;
    ID M_secondAdjacentElementPosition;
};





//! MeshElementMarked3D - Class for Volumes
/*!
    @author Luca Formaggia
 */
template
<typename GeoShape, typename MC = defaultMarkerCommon_Type>
class MeshElementMarked3D:
    public MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >,
    public MC::volumeMarker_Type
{
public:

    //! @name Public Types
    //@{

    //! Number of local Vertices
    static const UInt S_numLocalVertices = MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >::S_numVertices;
    //! Number of local Faces
    static const UInt S_numLocalFaces = MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >::S_numFaces;
    //! Number of local Edges (using Euler Formula)
    static const UInt S_numLocalEdges = MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >::S_numEdges;

    typedef GeoShape geoShape_Type;
    typedef typename MC::volumeMarker_Type marker_Type;
    typedef typename GeoShape::GeoBShape faceShape_Type;
    typedef typename faceShape_Type::GeoBShape edgeShape_Type;

    typedef MeshElementMarked1D<edgeShape_Type, MC, 3> edge_Type;
    typedef MeshElementMarked2D<faceShape_Type, MC, 3> face_Type;
    typedef MeshElementMarked0D<MC, 3> point_Type;
    typedef face_Type geoBElement_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked3D( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked3D to be copied
     */
    MeshElementMarked3D( const MeshElementMarked3D<GeoShape, MC>& Element );

    //! Destructor
    virtual ~MeshElementMarked3D()
    {
        // nothing to be done
    }

    //@}
};





//! MeshElementMarked3D - Class for Volumes
/*!
    @author Luca Formaggia
 */
template
<typename GeoShape, typename MC>
class MeshElementMarked2D<GeoShape, MC, 2>: public MeshElement<GeoShape, MeshElementMarked0D<MC, 2> >, public MC::volumeMarker_Type
{
public:

    //! @name Public Types
    //@{

    //! Number of local Vertices
    static const UInt S_numLocalVertices = MeshElement<GeoShape, MeshElementMarked0D<MC, 2> >::S_numVertices;
    //! Number of local Faces
    static const UInt S_numLocalEdges = MeshElement<GeoShape, MeshElementMarked0D<MC, 2> >::S_numEdges;

    typedef GeoShape geoShape_Type;
    typedef typename MC::volumeMarker_Type marker_Type;
    typedef typename GeoShape::GeoBShape faceShape_Type;
    typedef typename faceShape_Type::GeoBShape edgeShape_Type;

    typedef MeshElementMarked1D<edgeShape_Type, MC, 2> edge_Type;
    typedef MeshElementMarked2D<faceShape_Type, MC, 2> face_Type;
    typedef MeshElementMarked0D<MC, 2> point_Type;
    typedef face_Type geoBElement_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked2D( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked3D to be copied
     */
    MeshElementMarked2D( const MeshElementMarked2D<GeoShape, MC, 2>& Element );

    //! Destructor
    virtual ~MeshElementMarked2D()
    {
        // nothing to be done
    }

    //@}
};


/*-------------------------------------------------------------------------
  MeshElementMarked0D
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================
template <typename MC, int geoDim>
MeshElementMarked0D<MC, geoDim>::MeshElementMarked0D() :
        MeshVertex(), MC::pointMarker_Type()
{}

template <typename MC, int geoDim>
MeshElementMarked0D<MC, geoDim>::MeshElementMarked0D( ID identity, bool boundary ) :
        MeshVertex( identity, boundary ), MC::pointMarker_Type()
{}

template <typename MC, int geoDim>
MeshElementMarked0D<MC, geoDim>::MeshElementMarked0D( ID identity, Real x, Real y, Real z, bool boundary ) :
        MeshVertex( identity, x, y, z, boundary ), MC::pointMarker_Type()
{}

template <typename MC, int geoDim>
MeshElementMarked0D<MC, geoDim>::MeshElementMarked0D( MeshElementMarked0D<MC, geoDim> const & Element ) :
        MeshVertex( Element ), MC::pointMarker_Type( Element )
{}

// ==========================================
// Operators
// ==========================================
//! It calls operator= of base classes, just to be sure to do the right thing.
template <typename MC, int geoDim>
MeshElementMarked0D<MC, geoDim> &
MeshElementMarked0D<MC, geoDim>::operator = ( MeshElementMarked0D<MC, geoDim> const & Element )
{
    if ( this != &Element )
    {
        MeshVertex::operator=( Element );
        marker_Type::operator=( Element );
    }
    return *this;
}

/*-------------------------------------------------------------------------
  MeshElementMarked0D
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================
template <typename MC>
MeshElementMarked0D<MC, 2>::MeshElementMarked0D() :
	MeshVertex(), MC::pointMarker_Type()
{}

template <typename MC>
MeshElementMarked0D<MC, 2>::MeshElementMarked0D( ID identity, bool boundary ) :
MeshVertex( identity, boundary ), MC::pointMarker_Type()
{}

template <typename MC>
MeshElementMarked0D<MC, 2>::MeshElementMarked0D( ID identity, Real x, Real y, Real z, bool boundary ) :
	MeshVertex( identity, x, y, z, boundary ), MC::pointMarker_Type()
{}

template <typename MC>
MeshElementMarked0D<MC, 2>::MeshElementMarked0D( MeshElementMarked0D<MC, 2> const & Element ) :
	MeshVertex( Element ), MC::pointMarker_Type( Element )
{}

// ==========================================
// Operators
// ==========================================
//! It calls operator= of base classes, just to be sure to do the right thing.
template <typename MC>
MeshElementMarked0D<MC, 2> &
MeshElementMarked0D<MC, 2>::operator = ( MeshElementMarked0D<MC, 2> const & Element )
{
    if ( this != &Element )
    {
        MeshVertex::operator=( Element );
        marker_Type::operator=( Element );
    }
    return *this;
}

template <typename MC>
MeshElementMarked0D<MC, 2> const &
MeshElementMarked0D<MC, 2>::point ( ID const /*identity*/ ) const{
	return *this;
}

template <typename MC>
void
MeshElementMarked0D<MC, 2>::setPoint( ID const /*identity*/, MeshElementMarked0D<MC, 2> const * point )
{
	if(this != point)
		this = point;
	/*->x() = point->x();
	this->y() = point->y();
	this->z() = point->z();
	this->id() = point->id();*/
}

template <typename MC>
void
MeshElementMarked0D<MC, 2>::setPoint( ID const /*identity*/, MeshElementMarked0D<MC, 2> const & point )
{
	if(this != &point)
		*this = point;
	/*this->x() = point.x();
	this->y() = point.y();
	this->z() = point.z();
	this->id() = point.id();*/
}


/*-------------------------------------------------------------------------
  MeshElementMarked1D
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================

template <typename GeoShape, typename MC>
MeshElementMarked1D<GeoShape, MC, 3>::MeshElementMarked1D( ID identity ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >( identity )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 1 , "geoElement2D with incorrect GeoShape" ) ;
}


template <typename GeoShape, typename MC>
MeshElementMarked1D<GeoShape, MC, 3>::MeshElementMarked1D( const MeshElementMarked1D<GeoShape, MC, 3>& Element ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >( Element ),
        MC::edgeMarker_Type                   ( Element )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 1 , "geoElement2D with incorrect GeoShape" ) ;
}



/*-------------------------------------------------------------------------
  MeshElementMarked2D
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================
template <typename GeoShape, typename MC>
const UInt MeshElementMarked2D<GeoShape, MC, 3>::S_numLocalEdges;

template <typename GeoShape, typename MC>
MeshElementMarked2D<GeoShape, MC, 3>::MeshElementMarked2D( ID identity ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >( identity ),
        M_firstAdjacentElementIdentity   ( NotAnId ),
        M_secondAdjacentElementIdentity  ( NotAnId ),
        M_firstAdjacentElementPosition   ( NotAnId ),
        M_secondAdjacentElementPosition  ( NotAnId )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 2 , "geoElement2D with incorrect GeoShape" ) ;
}

template <typename GeoShape, typename MC>
MeshElementMarked2D<GeoShape, MC, 3>::MeshElementMarked2D( const MeshElementMarked2D<GeoShape, MC, 3>& Element ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >( Element ),
        MC::faceMarker_Type               ( Element ),
        M_firstAdjacentElementIdentity    ( Element.M_firstAdjacentElementIdentity),
        M_secondAdjacentElementIdentity   ( Element.M_secondAdjacentElementIdentity),
        M_firstAdjacentElementPosition    ( Element.M_firstAdjacentElementPosition ),
        M_secondAdjacentElementPosition   ( Element.M_secondAdjacentElementPosition )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 2 , "geoElement2D with incorrect GeoShape" ) ;
}


template <typename GeoShape, typename MC>
const UInt MeshElementMarked1D<GeoShape, MC, 2>::S_numLocalVertices;

template <typename GeoShape, typename MC>
MeshElementMarked1D<GeoShape, MC, 2>::MeshElementMarked1D( ID identity ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC,2> >( identity ),
        M_firstAdjacentElementIdentity   ( NotAnId ),
        M_secondAdjacentElementIdentity  ( NotAnId ),
        M_firstAdjacentElementPosition   ( NotAnId ),
        M_secondAdjacentElementPosition  ( NotAnId )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 1 , "geoElement1D in 2D geometry with incorrect GeoShape" ) ;
}

template <typename GeoShape, typename MC>
MeshElementMarked1D<GeoShape, MC,2>::MeshElementMarked1D( const MeshElementMarked1D<GeoShape, MC, 2>& Element ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC,2> >( Element ),
        MC::faceMarker_Type               ( Element ),
        M_firstAdjacentElementIdentity    ( Element.M_firstAdjacentElementIdentity),
        M_secondAdjacentElementIdentity   ( Element.M_secondAdjacentElementIdentity),
        M_firstAdjacentElementPosition    ( Element.M_firstAdjacentElementPosition ),
        M_secondAdjacentElementPosition   ( Element.M_secondAdjacentElementPosition )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 1 , "geoElement2D in 2D Geometry with incorrect GeoShape" ) ;
}


/*-------------------------------------------------------------------------
                 MeshElementMarked3D
 --------------------------------------------------------------------------*/
template <typename GeoShape, typename MC>
const UInt MeshElementMarked3D<GeoShape, MC>::S_numLocalVertices;
template <typename GeoShape, typename MC>
const UInt MeshElementMarked3D<GeoShape, MC>::S_numLocalFaces;
template <typename GeoShape, typename MC>
const UInt MeshElementMarked3D<GeoShape, MC>::S_numLocalEdges;


template <typename GeoShape, typename MC>
const UInt MeshElementMarked2D<GeoShape, MC, 2>::S_numLocalVertices;
template <typename GeoShape, typename MC>
const UInt MeshElementMarked2D<GeoShape, MC, 2>::S_numLocalEdges;

// ==========================================
// Constructor & Destructor
// ==========================================
template <typename GeoShape, typename MC>
MeshElementMarked3D<GeoShape, MC>::MeshElementMarked3D( ID identity ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >( identity )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 3 , "geoElement3D with incorrect GeoShape" )
}

template <typename GeoShape, typename MC>
MeshElementMarked3D<GeoShape, MC>::MeshElementMarked3D( const MeshElementMarked3D<GeoShape, MC>& Element ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC, 3> >( Element ),
        MC::volumeMarker_Type                ( Element )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 3 , "geoElement3D with incorrect GeoShape" )
}


template <typename GeoShape, typename MC>
MeshElementMarked2D<GeoShape, MC, 2>::MeshElementMarked2D( ID identity ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC, 2> >( identity )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 2 , "geoElement2D in 2D geometry with incorrect GeoShape" )
}

template <typename GeoShape, typename MC>
MeshElementMarked2D<GeoShape, MC, 2>::MeshElementMarked2D( const MeshElementMarked2D<GeoShape, MC, 2>& Element ) :
        MeshElement<GeoShape, MeshElementMarked0D<MC, 2> >( Element ),
        MC::volumeMarker_Type                  ( Element )
{
    ASSERT_PRE( GeoShape::S_nDimensions == 2 , "geoElement2D in 2D geometry with incorrect GeoShape" )
}

}
#endif

