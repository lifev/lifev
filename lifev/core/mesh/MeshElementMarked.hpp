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
    @contributor Mauro Perego <mdelia2@mathcs.emory.edu>
    @maintainer Marta D'Elia <mdelia2@mathcs.emory.edu>

    @date 00-00-0000

 */

#ifndef MESHELEMENTMARKED_H
#define MESHELEMENTMARKED_H

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/MarkerDefinitions.hpp>
#include <lifev/core/mesh/MeshElement.hpp>
#include <lifev/core/mesh/MeshElementBare.hpp>

namespace LifeV
{

//! Class for describing a geometric Entity immersed in 1D, 2D or 3D Geometry.
/* it is implemented using traits.
 * Template parameters:
 * elemDim: is the dimension of the entity: 0 for vertices, 1 for lines, 2 for triangles, 3 for tetrahedra.
 * geoDim: is the dimension of the geometry: (1, for a 1D geometry, 2 for a 2D geometry, 3 for a 3D geometry).
 * GeoShape: is the shape type: (eg. linearTetra, linearTriangle ecc. )
             notice that the dimension of GeoShape is GeoDim (except for GeoShape)
 * MC: is the marker
 */
template <int elemDim, int geoDim, typename GeoShape, typename MC>
class MeshElementMarked: public MeshVertex
{
public:
    typedef nullShape geoShape_Type;
    MeshElementMarked();
};


//! specialization for 0D entities (points).
template <int geoDim, typename GeoShape, typename MC>
class MeshElementMarked<0, geoDim, GeoShape, MC>: public MeshVertex, public MC::pointMarker_Type
{
public:

    //! @name Public Types
    //@{
    typedef typename MC::pointMarker_Type marker_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MeshElementMarked();

    //! Declares item identity and states if it is on boundary
    /*!
        @param identity Element identity
        @param boundary True if the element is on boundary
     */
    explicit MeshElementMarked ( ID identity, bool boundary = false );

    //! Declares item identity, provides coordinate and states if it is on boundary
    /*!
        @param identity Element identity
        @param x Element x coordinate
        @param y Element y coordinate
        @param z Element z coordinate
        @param boundary True if the element is on boundary
     */
    MeshElementMarked ( ID identity, Real x, Real y, Real z, bool boundary = false );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked0D to be copied
     */
    MeshElementMarked ( MeshElementMarked<0, geoDim, GeoShape, MC> const& Element );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked0D to be copied
        @param Marker Markercommon
     */
    MeshElementMarked ( MeshVertex const& Element, MC const& Marker );

    //! Destructor
    virtual ~MeshElementMarked()
    {
        // nothing to be done
    }

    //@}

    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
        @param Element Equivalent MeshElementMarked
        @return Reference to a new MeshElementMarked with the same content of MeshElementMarked Element
     */
    MeshElementMarked& operator = ( const MeshElementMarked<0, geoDim, GeoShape, MC>&   Element );

    //@}

    void setPoint ( ID const identity, MeshElementMarked<0, geoDim, GeoShape, MC> const* point );

    void setPoint ( ID const /*identity*/, MeshElementMarked<0, geoDim, GeoShape, MC> const& point );

    MeshElementMarked<0, geoDim, GeoShape, MC> const& point ( ID const identity ) const;


};

//! specialization for 0D entities (points) in 1D Geometry.
template
<typename GeoShape, typename MC>
class MeshElementMarked<0, 1, GeoShape, MC>: public MeshVertex, public MC::pointMarker_Type
{
public:

    //! @name Public Types
    //@{
    typedef typename MC::pointMarker_Type marker_Type;
    static const UInt S_numPoints = 1;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MeshElementMarked();

    //! Declares item identity and states if it is on boundary
    /*!
        @param identity Element identity
        @param boundary True if the element is on boundary
     */
    explicit MeshElementMarked ( ID identity, bool boundary = false );

    //! Declares item identity, provides coordinate and states if it is on boundary
    /*!
        @param identity Element identity
        @param x Element x coordinate
        @param y Element y coordinate
        @param z Element z coordinate
        @param boundary True if the element is on boundary
     */
    MeshElementMarked ( ID identity, Real x, Real y, Real z, bool boundary = false );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked0D to be copied
     */
    MeshElementMarked ( MeshElementMarked<0, 1, GeoShape, MC> const& Element );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked0D to be copied
        @param Marker Markercommon
     */
    MeshElementMarked ( MeshVertex const& Element, MC const& Marker );

    //! Destructor
    virtual ~MeshElementMarked()
    {
        // nothing to be done
    }

    //@}

    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
        @param Element Equivalent MeshElementMarked
        @return Reference to a new MeshElementMarked with the same content of MeshElementMarked Element
     */
    MeshElementMarked& operator = ( const MeshElementMarked<0, 1, GeoShape, MC>&   Element );

    //@}

    void setPoint ( ID const identity, MeshElementMarked<0, 1, GeoShape, MC> const* point );

    void setPoint ( ID const /*identity*/, MeshElementMarked<0, 1, GeoShape, MC> const& point );

    MeshElementMarked<0, 1, GeoShape, MC> const& point ( ID const identity ) const;

    //! @name Get Methods
    //@{

    //! Returns the identity of the first adjacent element
    /*!
        @return Identity of the first adjacent element
     */
    ID firstAdjacentElementIdentity() const
    {

        return M_firstAdjacentElementIdentity;
    }

    //! Returns the identity of the second adjacent element
    /*!
        @return Identity of the second adjacent element
     */
    ID secondAdjacentElementIdentity() const
    {
        return M_secondAdjacentElementIdentity;
    }

    //! Returns the identity of the first adjacent element
    /*!
        @return Identity of the first adjacent element
     */
    ID& firstAdjacentElementIdentity()
    {
        return M_firstAdjacentElementIdentity;
    }

    //! Returns the identity of the second adjacent element
    /*!
        @return Identity of the second adjacent element
     */
    ID& secondAdjacentElementIdentity()
    {
        return M_secondAdjacentElementIdentity;
    }

    //! Returns the position of the first adjacent element
    /*!
        @return Position of the first adjacent element
     */
    ID firstAdjacentElementPosition() const
    {
        return M_firstAdjacentElementPosition;
    }

    //! Returns the position of the second adjacent element
    /*!
        @return Position of the second adjacent element
     */
    ID secondAdjacentElementPosition() const
    {
        return M_secondAdjacentElementPosition;
    }


    //! Returns the position of the first adjacent element
    /*!
        @return Position of the first adjacent element
     */
    ID& firstAdjacentElementPosition()
    {
        return M_firstAdjacentElementPosition;
    }

    //! Returns the position of the second adjacent element
    /*!
        @return Position of the second adjacent element
     */
    ID& secondAdjacentElementPosition()
    {
        return M_secondAdjacentElementPosition;
    }

    //@}

private:
    ID M_firstAdjacentElementIdentity;
    ID M_secondAdjacentElementIdentity;
    ID M_firstAdjacentElementPosition;
    ID M_secondAdjacentElementPosition;

};


//! specialization for 1D entities (edges) in 1D Geometry.
template
<typename GeoShape, typename MC>
class MeshElementMarked<1, 1, GeoShape, MC> :
    public MeshElement<GeoShape, MeshElementMarked<0, 1, GeoShape, MC> >,
    public MC::edgeMarker_Type
{
public:

    //! @name Public Types
    //@{

    typedef GeoShape geoShape_Type;
    typedef typename MC::edgeMarker_Type marker_Type;
    typedef MeshElementMarked<0, 1, GeoShape, MC> geoBElement_Type;
    typedef MeshElementMarked<0, 1, GeoShape, MC> point_Type;
    static const UInt S_numLocalVertices = MeshElement<GeoShape, MeshElementMarked<0, 1, GeoShape, MC> >::S_numVertices;
    static const UInt S_numLocalFacets = S_numLocalVertices;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked ( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked to be copied
     */
    MeshElementMarked ( const MeshElementMarked<1, 1, GeoShape, MC>& Element);

    //! Destructor
    virtual ~MeshElementMarked()
    {
        // nothing to be done
    }
    //@}
};


//! specialization for 1D entities (edges) in 3D geometry.
template
<typename GeoShape, typename MC>
class MeshElementMarked<1, 3, GeoShape, MC> :
    public MeshElement<typename GeoShape::GeoBShape::GeoBShape, MeshElementMarked<0, 3, GeoShape, MC> >,
    public MC::edgeMarker_Type
{
public:

    //! @name Public Types
    //@{

    typedef typename GeoShape::GeoBShape::GeoBShape geoShape_Type;
    typedef typename MC::edgeMarker_Type marker_Type;
    typedef MeshElementMarked<0, 3, GeoShape, MC> geoBElement_Type;
    typedef MeshElementMarked<0, 3, GeoShape, MC> point_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked ( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked to be copied
     */
    MeshElementMarked ( const MeshElementMarked<1, 3, GeoShape, MC>& Element);

    //! Destructor
    virtual ~MeshElementMarked()
    {
        // nothing to be done
    }
    //@}
};



//! specialization for 1D entities (edges) in a 2D geometry. Identities of the adjacent 2D elements and their relative position are stored.
template
<typename GeoShape, typename MC>
class MeshElementMarked<1, 2, GeoShape, MC>: public MeshElement<typename GeoShape::GeoBShape, MeshElementMarked<0, 2, GeoShape, MC> >, public MC::faceMarker_Type
{

public:

    //! @name Public Types::GeoBShape
    //@{

    //! Number of element edges, for compatibility
    static const UInt S_numLocalVertices = MeshElement<typename GeoShape::GeoBShape, MeshElementMarked<0, 2, GeoShape, MC> >::S_numVertices;

    typedef typename GeoShape::GeoBShape geoShape_Type;
    typedef typename MC::faceMarker_Type marker_Type;
    typedef typename geoShape_Type::GeoBShape edgeShape_Type;
    typedef MeshElementMarked<1, 2, GeoShape, MC> edge_Type;
    typedef MeshElementMarked<0, 2, GeoShape, MC> point_Type;
    typedef edge_Type geoBElement_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked ( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked to be copied
     */
    MeshElementMarked ( const MeshElementMarked<1, 2, GeoShape, MC>& Element);

    //! Destructor
    virtual ~MeshElementMarked()
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
    }

    //! Returns the identity of the second adjacent element
    /*!
        @return Identity of the second adjacent element
     */
    ID secondAdjacentElementIdentity() const
    {
        return M_secondAdjacentElementIdentity;
    }

    //! Returns the identity of the first adjacent element
    /*!
        @return Identity of the first adjacent element
     */
    ID& firstAdjacentElementIdentity()
    {
        return M_firstAdjacentElementIdentity;
    }

    //! Returns the identity of the second adjacent element
    /*!
        @return Identity of the second adjacent element
     */
    ID& secondAdjacentElementIdentity()
    {
        return M_secondAdjacentElementIdentity;
    }

    //! Returns the position of the first adjacent element
    /*!
        @return Position of the first adjacent element
     */
    ID firstAdjacentElementPosition() const
    {
        return M_firstAdjacentElementPosition;
    }

    //! Returns the position of the second adjacent element
    /*!
        @return Position of the second adjacent element
     */
    ID secondAdjacentElementPosition() const
    {
        return M_secondAdjacentElementPosition;
    }


    //! Returns the position of the first adjacent element
    /*!
        @return Position of the first adjacent element
     */
    ID& firstAdjacentElementPosition()
    {
        return M_firstAdjacentElementPosition;
    }

    //! Returns the position of the second adjacent element
    /*!
        @return Position of the second adjacent element
     */
    ID& secondAdjacentElementPosition()
    {
        return M_secondAdjacentElementPosition;
    }

    //@}

private:
    ID M_firstAdjacentElementIdentity;
    ID M_secondAdjacentElementIdentity;
    ID M_firstAdjacentElementPosition;
    ID M_secondAdjacentElementPosition;
};


//! specialization for 2D entities (faces) in a 3D geometry. Identities of the adjacent 3D elements and their relative position are stored.
template
<typename GeoShape, typename MC>
class MeshElementMarked<2, 3, GeoShape, MC>:
    public MeshElement<typename GeoShape::GeoBShape, MeshElementMarked<0, 3, GeoShape, MC> >,
    public MC::faceMarker_Type
{

public:

    //! @name Public Types
    //@{

    //! Number of element edges, for compatibility
    static const UInt S_numLocalEdges = MeshElement<typename GeoShape::GeoBShape, MeshElementMarked<0, 3, GeoShape, MC> >::S_numEdges;

    typedef typename GeoShape::GeoBShape geoShape_Type;
    typedef typename MC::faceMarker_Type marker_Type;
    typedef typename geoShape_Type::GeoBShape edgeShape_Type;
    typedef MeshElementMarked<1, 3, GeoShape, MC> edge_Type;
    typedef MeshElementMarked<0, 3, GeoShape, MC> point_Type;
    typedef edge_Type geoBElement_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked ( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked to be copied
     */
    MeshElementMarked ( const MeshElementMarked<2, 3, GeoShape, MC>& Element);

    //! Destructor
    virtual ~MeshElementMarked()
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
    }

    //! Returns the identity of the second adjacent element
    /*!
        @return Identity of the second adjacent element
     */
    ID secondAdjacentElementIdentity() const
    {
        return M_secondAdjacentElementIdentity;
    }

    //! Returns the identity of the first adjacent element
    /*!
        @return Identity of the first adjacent element
     */
    ID& firstAdjacentElementIdentity()
    {
        return M_firstAdjacentElementIdentity;
    }

    //! Returns the identity of the second adjacent element
    /*!
        @return Identity of the second adjacent element
     */
    ID& secondAdjacentElementIdentity()
    {
        return M_secondAdjacentElementIdentity;
    }

    //! Returns the position of the first adjacent element
    /*!
        @return Position of the first adjacent element
     */
    ID firstAdjacentElementPosition() const
    {
        return M_firstAdjacentElementPosition;
    }

    //! Returns the position of the second adjacent element
    /*!
        @return Position of the second adjacent element
     */
    ID secondAdjacentElementPosition() const
    {
        return M_secondAdjacentElementPosition;
    }


    //! Returns the position of the first adjacent element
    /*!
        @return Position of the first adjacent element
     */
    ID& firstAdjacentElementPosition()
    {
        return M_firstAdjacentElementPosition;
    }

    //! Returns the position of the second adjacent element
    /*!
        @return Position of the second adjacent element
     */
    ID& secondAdjacentElementPosition()
    {
        return M_secondAdjacentElementPosition;
    }

    //@}

private:
    ID M_firstAdjacentElementIdentity;
    ID M_secondAdjacentElementIdentity;
    ID M_firstAdjacentElementPosition;
    ID M_secondAdjacentElementPosition;
};




//! specialization for 2D entities (faces) in a 2D geometry.
template
<typename GeoShape, typename MC>
class MeshElementMarked<2, 2, GeoShape, MC>: public MeshElement<GeoShape, MeshElementMarked<0, 2, GeoShape, MC> >, public MC::volumeMarker_Type
{
public:

    //! @name Public Types
    //@{

    //! Number of local Vertices
    static const UInt S_numLocalVertices = MeshElement<GeoShape, MeshElementMarked<0, 2, GeoShape, MC> >::S_numVertices;
    //! Number of local Faces
    static const UInt S_numLocalEdges = MeshElement<GeoShape, MeshElementMarked<0, 2, GeoShape, MC> >::S_numEdges;
    static const UInt S_numLocalFacets = S_numLocalEdges;
    static const UInt S_numLocalRidges = S_numLocalVertices;

    typedef GeoShape geoShape_Type;
    typedef typename MC::volumeMarker_Type marker_Type;
    typedef typename GeoShape::GeoBShape faceShape_Type;
    typedef typename faceShape_Type::GeoBShape edgeShape_Type;

    typedef MeshElementMarked<1, 2, GeoShape, MC> edge_Type;
    typedef MeshElementMarked<2, 2, GeoShape, MC> face_Type;
    typedef MeshElementMarked<0, 2, GeoShape, MC> point_Type;
    typedef face_Type geoBElement_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked ( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked3D to be copied
     */
    MeshElementMarked ( const MeshElementMarked<2, 2, GeoShape, MC>& Element );

    //! Destructor
    virtual ~MeshElementMarked()
    {
        // nothing to be done
    }

    //@}
};

//! specialization for 3D entities (cells) in a 3D geometry.
template
<typename GeoShape, typename MC>
class MeshElementMarked<3, 3, GeoShape, MC>:
    public MeshElement<GeoShape, MeshElementMarked<0, 3, GeoShape, MC> >,
    public MC::volumeMarker_Type
{
public:

    //! @name Public Types
    //@{

    //! Number of local Vertices
    static const UInt S_numLocalVertices = MeshElement<GeoShape, MeshElementMarked<0, 3, GeoShape, MC> >::S_numVertices;
    //! Number of local Faces
    static const UInt S_numLocalFaces = MeshElement<GeoShape, MeshElementMarked<0, 3, GeoShape, MC> >::S_numFaces;
    //! Number of local Edges (using Euler Formula)
    static const UInt S_numLocalEdges = MeshElement<GeoShape, MeshElementMarked<0, 3, GeoShape, MC> >::S_numEdges;

    static const UInt S_numLocalFacets = S_numLocalFaces;
    static const UInt S_numLocalRidges = S_numLocalEdges;

    typedef GeoShape geoShape_Type;
    typedef typename MC::volumeMarker_Type marker_Type;
    typedef typename GeoShape::GeoBShape faceShape_Type;
    typedef typename faceShape_Type::GeoBShape edgeShape_Type;

    typedef MeshElementMarked<1, 3, GeoShape, MC> edge_Type;
    typedef MeshElementMarked<2, 3, GeoShape, MC> face_Type;
    typedef MeshElementMarked<0, 3, GeoShape, MC> point_Type;
    typedef face_Type geoBElement_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Declares element identity
    /*!
        @param identity Element identity
     */
    explicit MeshElementMarked ( ID identity = NotAnId );

    //! Copy constructor
    /*!
        @param Element MeshElementMarked to be copied
     */
    MeshElementMarked ( const MeshElementMarked<3, 3, GeoShape, MC>& Element );

    //! Destructor
    virtual ~MeshElementMarked()
    {
        // nothing to be done
    }

    //@}
};


/*-------------------------------------------------------------------------
  MeshElementMarked (0D)
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================
template <int geoDim, typename GeoShape, typename MC>
MeshElementMarked<0, geoDim, GeoShape, MC>::MeshElementMarked() :
    MeshVertex(), MC::pointMarker_Type()
{}

template <int geoDim, typename GeoShape, typename MC>
MeshElementMarked<0, geoDim, GeoShape, MC>::MeshElementMarked ( ID identity, bool boundary ) :
    MeshVertex ( identity, boundary ), MC::pointMarker_Type()
{}

template <int geoDim, typename GeoShape, typename MC>
MeshElementMarked<0, geoDim, GeoShape, MC>::MeshElementMarked ( ID identity, Real x, Real y, Real z, bool boundary ) :
    MeshVertex ( identity, x, y, z, boundary ), MC::pointMarker_Type()
{}

template <int geoDim, typename GeoShape, typename MC>
MeshElementMarked<0, geoDim, GeoShape, MC>::MeshElementMarked ( MeshElementMarked<0, geoDim, GeoShape, MC> const& Element ) :
    MeshVertex ( Element ), MC::pointMarker_Type ( Element )
{}


// ==========================================
// Operators
// ==========================================
//! It calls operator= of base classes, just to be sure to do the right thing.
template <int geoDim, typename GeoShape, typename MC>
MeshElementMarked<0, geoDim, GeoShape, MC>&
MeshElementMarked<0, geoDim, GeoShape, MC>::operator = ( MeshElementMarked<0, geoDim, GeoShape, MC> const& Element )
{
    if ( this != &Element )
    {
        MeshVertex::operator= ( Element );
        marker_Type::operator= ( Element );
    }
    return *this;
}


// ==========================================
// Methods
// ==========================================
template <int geoDim, typename GeoShape, typename MC>
MeshElementMarked<0, geoDim, GeoShape, MC> const&
MeshElementMarked<0, geoDim, GeoShape, MC>::point ( ID const /*identity*/ ) const
{
    return *this;
}

template <int geoDim, typename GeoShape, typename MC>
void
MeshElementMarked<0, geoDim, GeoShape, MC>::setPoint ( ID const /*identity*/, MeshElementMarked<0, geoDim, GeoShape, MC> const* point )
{
    if (this != point)
    {
        this = point;
    }
}

template <int geoDim, typename GeoShape, typename MC>
void
MeshElementMarked<0, geoDim, GeoShape, MC>::setPoint ( ID const /*identity*/, MeshElementMarked<0, geoDim, GeoShape, MC> const& point )
{
    if (this != &point)
    {
        *this = point;
    }
}

/*-------------------------------------------------------------------------
  MeshElementMarked 0D in 1D geometry
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================
template <typename GeoShape, typename MC>
MeshElementMarked<0, 1, GeoShape, MC>::MeshElementMarked() :
    MeshVertex(), MC::pointMarker_Type()
{}

template <typename GeoShape, typename MC>
MeshElementMarked<0, 1, GeoShape, MC>::MeshElementMarked ( ID identity, bool boundary ) :
    MeshVertex ( identity, boundary ), MC::pointMarker_Type()
{}

template <typename GeoShape, typename MC>
MeshElementMarked<0, 1, GeoShape, MC>::MeshElementMarked ( ID identity, Real x, Real y, Real z, bool boundary ) :
    MeshVertex ( identity, x, y, z, boundary ), MC::pointMarker_Type()
{}

template <typename GeoShape, typename MC>
MeshElementMarked<0, 1, GeoShape, MC>::MeshElementMarked ( MeshElementMarked<0, 1, GeoShape, MC> const& Element ) :
    MeshVertex ( Element ), MC::pointMarker_Type ( Element )
{}

// ==========================================
// Operators
// ==========================================
//! It calls operator= of base classes, just to be sure to do the right thing.
template <typename GeoShape, typename MC>
MeshElementMarked<0, 1, GeoShape, MC>&
MeshElementMarked<0, 1, GeoShape, MC>::operator = ( MeshElementMarked<0, 1, GeoShape, MC> const& Element )
{
    if ( this != &Element )
    {
        MeshVertex::operator= ( Element );
        marker_Type::operator= ( Element );
    }
    return *this;
}


// ==========================================
// Methods
// ==========================================
template <typename GeoShape, typename MC>
MeshElementMarked<0, 1, GeoShape, MC> const&
MeshElementMarked<0, 1, GeoShape, MC>::point ( ID const /*identity*/ ) const
{
    return *this;
}

template <typename GeoShape, typename MC>
void
MeshElementMarked<0, 1, GeoShape, MC>::setPoint ( ID const /*identity*/, MeshElementMarked<0, 1, GeoShape, MC> const* point )
{
    if (this != point)
    {
        this = point;
    }
}

template <typename GeoShape, typename MC>
void
MeshElementMarked<0, 1, GeoShape, MC>::setPoint ( ID const /*identity*/, MeshElementMarked<0, 1, GeoShape, MC> const& point )
{
    if (this != &point)
    {
        *this = point;
    }
}


/*-------------------------------------------------------------------------
  MeshElementMarked (1D)
  --------------------------------------------------------------------------*/
// ==========================================
// Constructor & Destructor
// ==========================================

template <typename GeoShape, typename MC>
MeshElementMarked<1, 1, GeoShape, MC>::MeshElementMarked ( ID identity ) :
    MeshElement<GeoShape, MeshElementMarked<0, 1, GeoShape, MC> > ( identity )
{
    ASSERT_PRE ( GeoShape::S_nDimensions == 1 , "geoElement2D with incorrect GeoShape" ) ;
}


template <typename GeoShape, typename MC>
MeshElementMarked<1, 1, GeoShape, MC>::MeshElementMarked ( const MeshElementMarked<1, 1, GeoShape, MC>& Element ) :
    MeshElement<GeoShape, MeshElementMarked<0, 1, GeoShape, MC> > ( Element ),
    MC::edgeMarker_Type                   ( Element )
{
    ASSERT_PRE ( GeoShape::S_nDimensions == 1 , "geoElement2D with incorrect GeoShape" ) ;
}

/*-------------------------------------------------------------------------
  MeshElementMarked (1D in 3D geometry)
  --------------------------------------------------------------------------*/

// ==========================================
// Constructor & Destructor
// ==========================================

template <typename GeoShape, typename MC>
MeshElementMarked<1, 3, GeoShape, MC>::MeshElementMarked ( ID identity ) :
    MeshElement<geoShape_Type, MeshElementMarked<0, 3, GeoShape, MC> > ( identity )
{
    ASSERT_PRE ( geoShape_Type::S_nDimensions == 1 , "geoElement1D with incorrect GeoShape" ) ;
}


template <typename GeoShape, typename MC>
MeshElementMarked<1, 3, GeoShape, MC>::MeshElementMarked ( const MeshElementMarked<1, 3, GeoShape, MC>& Element ) :
    MeshElement<geoShape_Type, MeshElementMarked<0, 3, GeoShape, MC> > ( Element ),
    MC::edgeMarker_Type                   ( Element )
{
    ASSERT_PRE ( geoShape_Type::S_nDimensions == 1 , "geoElement1D with incorrect GeoShape" ) ;
}



/*-------------------------------------------------------------------------
  MeshElementMarked (1D in 2D geometry)
  --------------------------------------------------------------------------*/

template <typename GeoShape, typename MC>
const UInt MeshElementMarked<1, 2, GeoShape, MC>::S_numLocalVertices;

// ==========================================
// Constructor & Destructor
// ==========================================

template <typename GeoShape, typename MC>
MeshElementMarked<1, 2, GeoShape, MC>::MeshElementMarked ( ID identity ) :
    MeshElement<geoShape_Type, MeshElementMarked<0, 2, GeoShape, MC> > ( identity ),
    M_firstAdjacentElementIdentity   ( NotAnId ),
    M_secondAdjacentElementIdentity  ( NotAnId ),
    M_firstAdjacentElementPosition   ( NotAnId ),
    M_secondAdjacentElementPosition  ( NotAnId )
{
    ASSERT_PRE ( geoShape_Type::S_nDimensions == 1 , "geoElement1D in 2D geometry with incorrect GeoShape" ) ;
}

template <typename GeoShape, typename MC>
MeshElementMarked<1, 2, GeoShape, MC>::MeshElementMarked ( const MeshElementMarked<1, 2, GeoShape, MC>& Element ) :
    MeshElement<geoShape_Type, MeshElementMarked<0, 2, GeoShape, MC> > ( Element ),
    MC::faceMarker_Type               ( Element ),
    M_firstAdjacentElementIdentity    ( Element.M_firstAdjacentElementIdentity),
    M_secondAdjacentElementIdentity   ( Element.M_secondAdjacentElementIdentity),
    M_firstAdjacentElementPosition    ( Element.M_firstAdjacentElementPosition ),
    M_secondAdjacentElementPosition   ( Element.M_secondAdjacentElementPosition )
{
    ASSERT_PRE ( geoShape_Type::S_nDimensions == 1 , "geoElement1D in 2D Geometry with incorrect GeoShape" ) ;
}


/*-------------------------------------------------------------------------
  MeshElementMarked (2D in 3D geometry)
  --------------------------------------------------------------------------*/

template <typename GeoShape, typename MC>
const UInt MeshElementMarked<2, 3, GeoShape, MC>::S_numLocalEdges;


// ==========================================
// Constructor & Destructor
// ==========================================

template <typename GeoShape, typename MC>
MeshElementMarked<2, 3, GeoShape, MC>::MeshElementMarked ( ID identity ) :
    MeshElement<geoShape_Type, MeshElementMarked<0, 3, GeoShape, MC> > ( identity ),
    M_firstAdjacentElementIdentity   ( NotAnId ),
    M_secondAdjacentElementIdentity  ( NotAnId ),
    M_firstAdjacentElementPosition   ( NotAnId ),
    M_secondAdjacentElementPosition  ( NotAnId )
{
    ASSERT_PRE ( geoShape_Type::S_nDimensions == 2 , "geoElement2D with incorrect GeoShape" ) ;
}

template <typename GeoShape, typename MC>
MeshElementMarked<2, 3, GeoShape, MC>::MeshElementMarked ( const MeshElementMarked<2, 3, GeoShape, MC>& Element ) :
    MeshElement<geoShape_Type, MeshElementMarked<0, 3, GeoShape, MC> > ( Element ),
    MC::faceMarker_Type               ( Element ),
    M_firstAdjacentElementIdentity    ( Element.M_firstAdjacentElementIdentity),
    M_secondAdjacentElementIdentity   ( Element.M_secondAdjacentElementIdentity),
    M_firstAdjacentElementPosition    ( Element.M_firstAdjacentElementPosition ),
    M_secondAdjacentElementPosition   ( Element.M_secondAdjacentElementPosition )
{
    ASSERT_PRE ( geoShape_Type::S_nDimensions == 2 , "geoElement2D with incorrect GeoShape" ) ;
}



/*-------------------------------------------------------------------------
  MeshElementMarked (2D in 2D geometry)
  --------------------------------------------------------------------------*/

template <typename GeoShape, typename MC>
const UInt MeshElementMarked<2, 2, GeoShape, MC>::S_numLocalVertices;
template <typename GeoShape, typename MC>
const UInt MeshElementMarked<2, 2, GeoShape, MC>::S_numLocalEdges;

// ==========================================
// Constructor & Destructor
// ==========================================

template <typename GeoShape, typename MC>
MeshElementMarked<2, 2, GeoShape, MC>::MeshElementMarked ( ID identity ) :
    MeshElement<GeoShape, MeshElementMarked<0, 2, GeoShape, MC> > ( identity )
{
    ASSERT_PRE ( GeoShape::S_nDimensions == 2 , "geoElement2D in 2D geometry with incorrect GeoShape" )
}

template <typename GeoShape, typename MC>
MeshElementMarked<2, 2, GeoShape, MC>::MeshElementMarked ( const MeshElementMarked<2, 2, GeoShape, MC>& Element ) :
    MeshElement<GeoShape, MeshElementMarked<0, 2, GeoShape, MC> > ( Element ),
    MC::volumeMarker_Type                  ( Element )
{
    ASSERT_PRE ( GeoShape::S_nDimensions == 2 , "geoElement2D in 2D geometry with incorrect GeoShape" )
}



/*-------------------------------------------------------------------------
  MeshElementMarked (3D in 3D geometry)
  --------------------------------------------------------------------------*/
template <typename GeoShape, typename MC>
const UInt MeshElementMarked<3, 3, GeoShape, MC>::S_numLocalVertices;
template <typename GeoShape, typename MC>
const UInt MeshElementMarked<3, 3, GeoShape, MC>::S_numLocalFaces;
template <typename GeoShape, typename MC>
const UInt MeshElementMarked<3, 3, GeoShape, MC>::S_numLocalEdges;

// ==========================================
// Constructor & Destructor
// ==========================================

template <int elemDim, int geoDim, typename GeoShape, typename MC>
MeshElementMarked<elemDim, geoDim, GeoShape, MC>::MeshElementMarked()
{
    ;
}


template <typename GeoShape, typename MC>
MeshElementMarked<3, 3, GeoShape, MC>::MeshElementMarked ( ID identity ) :
    MeshElement<GeoShape, MeshElementMarked<0, 3, GeoShape, MC> > ( identity )
{
    ASSERT_PRE ( GeoShape::S_nDimensions == 3 , "geoElement3D with incorrect GeoShape" )
}

template <typename GeoShape, typename MC>
MeshElementMarked<3, 3, GeoShape, MC>::MeshElementMarked ( const MeshElementMarked<3, 3, GeoShape, MC>& Element ) :
    MeshElement<GeoShape, MeshElementMarked<0, 3, GeoShape, MC> > ( Element ),
    MC::volumeMarker_Type                ( Element )
{
    ASSERT_PRE ( GeoShape::S_nDimensions == 3 , "geoElement3D with incorrect GeoShape" )
}




}
#endif

