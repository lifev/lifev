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
/*! file regionMesh2D.h
  \brief The 2D mesh classes interfaces
  \version $Revision: 1.7 $ Luca Formaggia

  Introduces the RegionMesh2D class
*/

#ifndef _REGIONMESH2D_HH_
#define _REGIONMESH2D_HH_

#include "life.hpp"
#include "geoElement.hpp"
#include "switch.hpp"

#include "bareItems.hpp"
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "SimpleVect.hpp"  /* stl wrap of vector class template.
It supports numbering from one* /

namespace LifeV
{

/* ---------------------------------------------------------------------
                            RegionMesh 2D
-----------------------------------------------------------------------*/
//! The Region Mesh Class for 2D elements
/*!
  This is the class that stores the mesh entities for a single 2D region In
  a region elements are all of the same type
*/
template <typename GEOSHAPE, typename MC = DefMarkerCommon >
class RegionMesh2D : public MeshEntity, public MC::RegionMarker
{

public:


    explicit RegionMesh2D( ID id = 0 );
    explicit RegionMesh2D( RegionMesh2D<GEOSHAPE, MC> const & m );


    RegionMesh2D<GEOSHAPE, MC> operator=( RegionMesh2D<GEOSHAPE, MC> const & m );

    //! \name Markers_Types
    /*! From Marker Common (MC) template parameter */
    //@{
    typedef typename MC::PointMarker PointMarker;
    typedef typename MC::EdgeMarker EdgeMarker;
    typedef typename MC::FaceMarker FaceMarker;
    typedef typename MC::RegionMarker RegionMarker;
    typedef typename MC::RegionMarker Marker;
    //@}

    //! \name Basic_Element_Shapes_Types
    //@{
    typedef GEOSHAPE FaceShape;
    typedef typename GEOSHAPE::GeoBShape EdgeShape;
    //@}

    //! \name Geometric_Element_Types
    //@{
    typedef GeoElement2D<FaceShape, MC> FaceType;
    typedef GeoElement1D<EdgeShape, MC> EdgeType;
    typedef GeoElement0D<MC> PointType;
    //@}

    //! \name GeoElement_Container_Types
    /*! Typedefs for STL compliant containers of mesh geometric entities
      I Use SimpleVect for addressing from 1. */
    //@{
    //! Points Container
    typedef SimpleVect<PointType> Points;
    //! Faces Container
    typedef SimpleVect<FaceType> Faces;
    //! Edges Container: at least boundary edges
    typedef SimpleVect<EdgeType> Edges;
    //@}
    /*! \name Generic_Types
     * Generic types for all regionmeshXX These are part
     * of the generic generic interface common for all RegionMeshes (3D --
     * 1D).
     */
    //@{
    typedef GEOSHAPE ElementShape;
    typedef typename GEOSHAPE::GeoBShape BElementShape;
    typedef GeoElement2D<GEOSHAPE, MC> ElementType;
    typedef GeoElement1D<EdgeShape, MC> BElementType;
    typedef SimpleVect<FaceType > Elements;
    typedef SimpleVect<EdgeType> BElements;
    //@}

    /*! \name Switches_Methods

    * Switches are used to store the status of the RegionMesh The switches
    * are used internally to control whether some data structures have been
    * set up.
    *
    *     The possible Switches are
    *
    *    \verbatim
    *
    *     HAS_ALL_EDGES            HAS_FACE_TO_EDGES
    *     HAS_BEEN_CHECKED
    *     HAS_BOUNDARY_EDGES       EDGES_HAVE_ADIACENCY
    *
    *     \endverbatim
    *
    */

    //@{
    //! Returns the number of switches which have been set
    const UInt numSwitches() const
    {
        return switches.size();
    };
    //! Interrogate Switch
    bool getLinkSwitch( std::string const & _s ) const;
    //! Set a switch
    void setLinkSwitch( std::string const & _s );
    //! uset a switch
    void unsetLinkSwitch( std::string const & _s );
    //@}
    UInt numLocalVertices() const; //!< Number of local vertices for each (2D) element
    UInt numLocalEdges() const;  //!< Number of local edges for each (2D) element

    /*! \name Generic_Methods Generic methods for all regionmeshXX
     * These are the generic methods to get information about the number of
     * elements.  It is a generic interface common for all RegionMeshes (3D --
     * 1D)
     */
    //@{
    UInt numElements() const; //!< Number of elements in mesh (alias to numFaces())
    UInt & numElements(); //!< Number of 3D elements
    UInt numBElements() const; //!< Number of Boundary faces
    UInt & numBElements(); //!< Number of boundary faces
    ElementType & element( ID const & i );
    ElementType const & element( ID const & i ) const;
    BElementType & bElement( ID const & i );
    BElementType const & bElement( ID const & i ) const;
    //@}

    /* ============================================
                      Face Related Methods
       ============================================*/
    //! \name Face_Methods All methods which operates on 2D elements
    //@{
    UInt numFaces() const;    /*!< Returns number of Face elements in the mesh
             as given by the internal counter.*/
    UInt & numFaces();    /*!< Access number of Face (internal counter)*/
    UInt storedFaces() const; //!< faces actully stored in list
    UInt maxNumFaces() const; /*!< Current capacity of Faces Container,
             i.e. how many elements may be stored */
    void setMaxNumFaces( UInt const n, bool const setcounter = false ); /*!< Changes Current capacity of Faces
             (Optionally sets internal counter.) */
    FaceType & addFace(); //!< Adds faces. Id computed automatically. Return ref to added V.
    FaceType & addFace( FaceType const & v );
    FaceType & setFace( FaceType const & v, ID const pos ); //!< Add face to a specified position
    void setFaceCounter(); //! set numFaces counter
    FaceType & lastFace(); //!< Reference to last face stored in list. Useful for mesh readers
    FaceType const & face( ID const i ) const; //!< ith mesh 2Delement
    FaceType & face( ID const i ); //!<ith mesh 2Delement
    //@}

    /* ============================================
                      Edge Related Methods
       ============================================*/

    //! \name  EdgeMethods Methods to access/create/modify edges data
    /*!  There are different point counters which may bi interrogated or set:

      <UL> <LI>Number of Edges: is the declared number of total edges in the
      mesh. \b Beware: A value different from zero does NOT imply that the
      edges are actually stored.</LI>

      <LI>Number of Boundary Edges: is the declared number of boundary edges
      in the mesh. \b Beware: A value different from zero does NOT imply that
      the boundary edges are actually stored.</LI>

      <LI>Number of Stored Edges: It is the number of Edges actually stored
      on the edge container.</LI>

      <LI>Maximum number of stored edges. The number of edges that may stored
      before the container is resized. \b Important: This parameter has to be
      set BEFORE inserting edges in the container if we want that pointer
      into the container maintains their validity. The container will also
      have a better performance.</LI> </UL> \note To have more information on
      the Edge methods look at the documentation of th eanalogous Faces
      methods */

    //@{
    UInt numEdges() const; //!<Number of total edges in mesh (uses counter).
    UInt & numEdges();    //!<Number of total edges (uses counter, may modify)
    UInt storedEdges() const; //!< Number of stored edges
    UInt maxNumEdges() const;  //!< Max number of edges that can be stored
    //! Set Maximum Number of Edges in the Edges container
    void setMaxNumEdges( UInt const n, bool const setcounter = false );
    //! Adds a edge to list. Returns reference to it
    EdgeType & addEdge( bool const boundary = false );
    //! Adds a edge (optionally a boundary edge) to the end of the list
    //! and adjourn its ID. Returns reference to the added edge
    EdgeType & addEdge( EdgeType const & f, bool const boundary = false );
    //! Adds a edge (optionally a boundary edge) and adjourn its ID
    EdgeType & setEdge( EdgeType const & f, ID position, bool const boundary = false );
    //! The last edge in the container
    EdgeType & lastEdge();
    EdgeType const & edge( ID const i ) const; //<! ith mesh edge
    EdgeType & edge( ID const i ); //!< ith mesh edge (may  modify!)
    EdgeType const & boundaryEdge( ID const i ) const; //!< ith boundary edge.
    EdgeType & boundaryEdge( ID const i ); //!< ith boundary  edge.
    void setNumBEdges( UInt const n ) ; //<! Set counter of boundary edges
    bool hasEdges() const;  //<! Do I store mesh edges?
    bool hasInternalEdges() const; //<! Do I store also internal edges?
    UInt numBEdges() const;  //<! Number of Boundary Edges
    bool isBoundaryEdge( EdgeType const & f ) const;  //<!Is this edge on boundary?
    bool isBoundaryEdge( ID const & id ) const;  //<!Is  edge whose ID is id on boundary?
    /*! Does this ID corresponds to a full edge? A FULL EDGE is a 2DElement  that is actually
    stored in the Edge container */
    bool isFullEdge( UInt const & id ) const;

    /*@{*/
    //!ID of the  Element adjacent to a Edge.
    /*!  Edge ID given.  Pos =1 or 2 indicates first or second element.  The
      first element is the one <em>ORIENTED coherently with the edge</em> (AS
      STORED in Edges). It means that the edge orientation is OUTWARD with
      respect to the element. The second element is either null (boundary
      edge) or indicates that the normal of the edge appears INWARD with
      respect to that element*/

    UInt edgeElement( ID const edgeId, UInt const Pos ) const;
    //! Element adjacent to a EDGE. Edge reference given.
    UInt edgeElement( EdgeType const & f, UInt const Pos ) const;
    /*@}*/
    //@}


    /* ============================================
                      Points/Vertices Related Methods
       ============================================*/
    //! \name PointMethods Methods to access/create/modify Points/Vertices
    //!data
    /*!  There are different Point counters which may be interrogated or set
         <UL> <LI>Number of Points: is the declared number of total points in
         the mesh.</LI>

      <LI>Number of Boundary Points: is the declared number of boundary
      points in the mesh.  does NOT imply that the boundary points are
      actually stored.</LI>

      <LI>Number of Stored Points: It is the number of Points actually stored
      on the points container.</LI>

      <LI>Maximum number of stored points. The number of points that may
      stored before the container is resized.  \b Very \b Important: This
      parameter has to be set \b BEFORE inserting points in the
      container. since Geoelements will contain POINTERS into the container!
      This is a debatable point! </LI> </UL>

      \note To have more information on the Point methods look at the
      documentation of the analogous Edges methods */

    //@{
    //! Returns number of points in the mesh (interrogate counter)
    UInt numPoints() const;
    //! Changes number of points in the mesh (interrogate counter)
    UInt & numPoints();
    //! Returns number of points currently stored in Points container
    /*! (interrogates container)*/
    UInt storedPoints() const;
    UInt storedBPoints() const;
    UInt maxNumPoints() const;
    void setMaxNumPoints( UInt const n, bool const setcounter = false );
    PointType & addPoint( bool const boundary = false, bool const vertices = false );
    PointType & addPoint( PointType const & p, bool const boundary = false, bool const vertices = false );
    PointType & setPoint( PointType const & p, ID const position, bool const boundary = false, bool const vertices = false ); //!< adds point
    PointType & lastPoint(); //!< last mesh point
    PointType const & point( ID const i ) const; //!< ith mesh point
    PointType & point( ID const i ); //!< ith mesh point
    PointType const & boundaryPoint( ID const i ) const; //!< ith b. point.
    PointType & boundaryPoint( ID const i ); //!< ith b. point.
    UInt numBPoints() const; //!< counter of boundary points
    void setNumBPoints( UInt const n ); //<! Sets counter of boundary points
    UInt numVertices() const; //!< Number of vertices in Region
    UInt & numVertices(); //!< Allows to change number of vertices in Region
    UInt numBVertices() const; //!< Number of boundary vertices in RegionMesh
    UInt & numBVertices();
    bool isVertex( ID const & id ) const;  //<!Is this point a Vertex?
    bool isVertex( PointType const & p ) const;  //<!Is this point a Vertex?
    bool isBoundaryPoint( PointType const & p ) const;  //<!Is this point on boundary?
    bool isBoundaryPoint( ID const & id ) const;  //<!Is this point on boundary?
    //@}

    /*! \name Element_Adiacency_Methods

     Methods to obtain the ID of Edges belonging to an element.
     Accessing this information requires that the appropriate data
     structures have been set up by using the updateElementEdges() method.

     Often, it is NOT required to have the full information about edges and faces:
     The ID of the Face and Edge entities may be calculated without
     contructing the corresponding Edge of Face Object. This saves memory.  */
    //@{

    //! Is the array for local Edges set up?
    /*! It does not use switches, but interrogates the container directly*/
    bool hasLocalEdges() const;  //!<Check wether edge info is available
    //!Edge n. locF around face. Returns edgeID
    UInt localEdgeId( const FaceType & iface, UInt const locE ) const;
    //!Edge n. locF around face. Returns edgeID
    UInt localEdgeId( ID const facId, UInt const locE ) const;
    void updateElementEdges(); //!<Build localEdgeId table
    //! Destroys element-to-face container
    void cleanElementEdges();
    //@}


    //! Prints some mesh info
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout );
    //! Basic tests for mesh consistency.
    /*! For more estensive test see \link mesh_util.h */
    int check( int level = 0, bool const fix = false, bool const verbose = true, std::ostream & out = std::cerr );

    /*! \name RegionContainers
     *

     * Stl compliant containers for basic structure I expose them since they
     * are standard containers
     */
    //@{
    Points pointList; //!< Container of mesh points/verices
    Faces faceList;  //!< Container of mesh  Faces
    Edges edgeList;  //!< Container of mesh Edges.
    SimpleVect<PointType * > _bPoints; //!< Boundary points list
    //@}

    //! Switches
    // \sa Switch
    Switch switches;

protected:
    /*! Arrays containing the ids of Edges and Faces of each element
    I use a Define to use localto global array or directly the
    bareedges */
#ifdef SAVEMEMORY

    BareItemsHandler<BareEdge> _FToE;
#else

    SimpleArray<UInt> _FToE;
#endif

#ifdef NOT_BDATA_FIRST

    SimpleVect<EdgeType * > _bEdges;
#endif
    // function prototipes
    template < typename T >
    UInt numItems( SimpleVect< T> const & list ) const;
    template < typename T >
    UInt maxNumItems( SimpleVect< T> const & list ) const;
    template < typename T >
    void setMaxNumItems( SimpleVect< T> & list, UInt n, char * title );


    // Internal counters
    UInt _numVertices;
    UInt _numBVertices;
    UInt _numPoints;
    UInt _numBPoints;
    UInt _numFaces;
    UInt _numBFaces;
    UInt _numEdges;
    UInt _numBEdges;
};

/* ---------------------------------------------------------------------
                            RegionMesh2D Implementations
-----------------------------------------------------------------------*/
void set_switches_for_regionmesh2( Switch & sw )
{
    sw.create( "HAS_ALL_EDGES" );
    sw.create( "HAS_BOUNDARY_EDGES" );
    sw.create( "HAS_FACE_TO_EDGES" );
    sw.create( "HAS_BEEN_CHECKED" );
    sw.create( "EDGES_HAVE_ADIACENCY" );
}


template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::RegionMesh2D( ID id ) :
        MeshEntity( id ),
        MC::RegionMarker(),
        switches(),
        _numVertices( 0 ),
        _numBVertices( 0 ),
        _numPoints( 0 ),
        _numBPoints( 0 ),
        _numFaces( 0 ),
        _numEdges( 0 ),
        _numBEdges( 0 )
{
    set_switches_for_regionmesh( switches );
}

template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::RegionMesh2D( RegionMesh2D<GEOSHAPE, MC> const & m )
{
    ASSERT( true, "Copy Costructor Not Yet Implemented for RegionMesh2D" ) ;
}

template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>
RegionMesh2D<GEOSHAPE, MC>::operator=( RegionMesh2D<GEOSHAPE, MC> const & m )
{
    ASSERT( true, "Assignement Operator  Yet Not Implemented for RegionMesh2D" ) ;
}

template <typename GEOSHAPE, typename MC>
INLINE
void
RegionMesh2D<GEOSHAPE, MC>::setLinkSwitch( std::string const & _s )
{
    ASSERT0( switches.set( _s ), "Switch named " << _s << " is not allowed" );
};

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::getLinkSwitch( std::string const & _s ) const
{
    switches.test( _s );
};

template <typename GEOSHAPE, typename MC>
INLINE
void
RegionMesh2D<GEOSHAPE, MC>::unsetLinkSwitch( std::string const & _s )
{
    ASSERT0( switches.unset( _s ), "Switch named " << _s << " is not allowed" );
};


template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numLocalVertices() const
{
    return FaceType::numLocalVertices;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numLocalEdges() const
{
    return FaceType::numLocalEdges;
}

// ************** Generic Methods
template <typename GEOSHAPE, typename MC>
UInt RegionMesh2D<GEOSHAPE, MC>::numElements() const
{
    return _numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt & RegionMesh2D<GEOSHAPE, MC>::numElements()
{
    return _numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt RegionMesh2D<GEOSHAPE, MC>::numBElements() const
{
    return _numBEdges;
}

template <typename GEOSHAPE, typename MC>
UInt &RegionMesh2D<GEOSHAPE, MC>::numBElements()
{
    return _numBEdges;
}

template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::ElementType &
RegionMesh2D<GEOSHAPE, MC>::element( ID const & i )
{
    return face( i );
}

template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::ElementType const &
RegionMesh2D<GEOSHAPE, MC>::element( ID const & i ) const
{
    return face( i );
}

template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::BElementType &
RegionMesh2D<GEOSHAPE, MC>::bElement( ID const & i )
{
    return boundaryEdge( i );
}

template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::BElementType const &
RegionMesh2D<GEOSHAPE, MC>::bElement( ID const & i ) const
{
    return boundaryEdge( i );
}


//-------------------------------------------------------
// Templates for interrogate and modify list capacities
// \todo make them global and friends to regionmeshXX.
//--------------------------------------------------------
template <typename GEOSHAPE, typename MC>
template <typename T>
UInt
RegionMesh2D<GEOSHAPE, MC>::numItems( SimpleVect< T> const & list ) const
{
    return list.size();
}

template <typename GEOSHAPE, typename MC>
template <typename T>
UInt
RegionMesh2D<GEOSHAPE, MC>::maxNumItems( SimpleVect< T> const & list ) const
{
    return list.capacity();
}

template <typename GEOSHAPE, typename MC>
template <typename T>
void
RegionMesh2D<GEOSHAPE, MC>::
setMaxNumItems( SimpleVect< T> & list, UInt n, char * title )
{
    if ( list.capacity() == 0 )
    {
        list.reserve( n );
    }
    else if ( list.capacity() == n )
    {
#ifdef VERBOSE
        std::cerr << "WARNING: Capacity of " << title << "list already set to " << n << std::endl;
#endif

    }
    else
    {
#ifdef VERBOSE
        std::cerr << "WARNING: Resetting " << title << " list size to " << n << std::endl;
        std::cerr << "         ALL PREVIOUS POINTERS TO THE LIST (IF ANY) ARE NOW INVALID" << std::endl;
#endif

        list.reserve( n );
    }
}
//-------------------------------------------------------------------

// ***************************** FACES
template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numFaces() const
{
    return _numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numFaces()
{
    return _numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::maxNumFaces() const
{
    return maxNumItems( faceList );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::storedFaces() const
{
    return numItems( faceList );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumFaces( UInt const n, bool const setcounter )
{
    setMaxNumItems( faceList, n, "Face" );
    if ( setcounter )
        _numFaces = n;
}

// \todo use addItem
template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::addFace()
{
    return addFace( FaceType() );
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::addFace( FaceType const & v )
{
    ASSERT_PRE( faceList.size() < faceList.capacity() , "Face list size exceeded" <<
                faceList.size() + 1 << " " << faceList.capacity() ) ;
    faceList.push_back( v );
    ( faceList.back() ).id() = faceList.size();
    return faceList.back();
}
// \todo Use setItem

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::setFace( FaceType const & v, ID const pos )
{
    ASSERT_PRE( pos <= faceList.capacity() , "position requested exceed capacity" <<
                pos << " " << faceList.capacity() ) ;
    faceList( pos ) = v;
    faceList( pos ).id() = pos;
    return faceList( pos );
}

template <typename GEOSHAPE, typename MC>
INLINE
void
RegionMesh2D<GEOSHAPE, MC>::setFaceCounter()
{
    _numFaces = faceList.size();
}


template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::lastFace()
{
    return faceList.back();
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::FaceType const &
RegionMesh2D<GEOSHAPE, MC>::face( ID const i ) const
{
    ASSERT_BD( i > 0 && i <= faceList.size() ) ;
    return faceList( i );
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::face( ID const i )
{
    ASSERT_BD( i > 0 && i <= faceList.size() ) ;
    return faceList( i );
}

// ************************* EDGES ******************************
template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numEdges() const
{
    return _numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numEdges()
{
    return _numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::storedEdges() const
{
    return numItems( edgeList );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::maxNumEdges() const
{
    return maxNumItems( edgeList );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumEdges( UInt const n, bool const setcounter )
{
    setMaxNumItems( edgeList, n, "Edge" );
    if ( setcounter )
        _numEdges = n;
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::addEdge( bool const boundary )
{
    return addEdge( EdgeType(), boundary );
}


template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::addEdge( EdgeType const & f, bool const boundary )
{
    ASSERT_PRE( edgeList.size() < edgeList.capacity(), "Edge list size exceeded" <<
                edgeList.size() + 1 << " " << edgeList.capacity() ) ;
    edgeList.push_back( f );
    ( edgeList.back() ).id() = edgeList.size();
    if ( boundary )
    {
#ifdef NOT_BDATA_FIRST
        ASSERT_PRE( _bEdges.size() < _bEdges.capacity(), "Boundary Edge list size exceeded" <<
                    _bEdges.size() + 1 << " " << bEdges.capacity() ) ;
        _bEdges.push_back( &edgeList.back() );
#endif

    }
    return edgeList.back();
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::setEdge( EdgeType const & f, ID position, bool const boundary )
{
    ASSERT_PRE( position <= edgeList.capacity(), "Edge list size exceeded" <<
                position << " " << edgeList.capacity() ) ;
    edgeList( position ) = f;
    edgeList( position ).id() = position;
#ifdef NOT_BDATA_FIRST

    if ( boundary )
    {
        ASSERT_PRE( position <= _bEdges.capacity(), "Boundary Edge list size exceeded" <<
                    _bEdges.size() << " " << bEdges.capacity() ) ;
        _bEdges.push_back( &( edgeList( position ) ) );
    }
#endif
    return edgeList( position );
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::lastEdge()
{
    return edgeList.back();
}


template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::EdgeType const &
RegionMesh2D<GEOSHAPE, MC>::edge( ID const i ) const
{
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::edge( ID const i )
{
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
}


template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::EdgeType const &
RegionMesh2D<GEOSHAPE, MC>::boundaryEdge( ID const i ) const
{
#ifdef NOT_BDATA_FIRST
    ASSERT_PRE( _bEdges.size() != 0 " Boundary Edges not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bEdges.size() ) ;
    return *( _bEdges( i ) );
#else

    ASSERT_PRE( edgeList.size() != 0, "Boundary Edges not stored" ) ;
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
#endif
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::boundaryEdge( ID const i )
{
#ifdef NOT_BDATA_FIRST
    ASSERT_PRE( _bEdges.size() != 0 " Boundary Edges not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bEdges.size() ) ;
    return *( _bEdges( i ) );
#else

    ASSERT_PRE( edgeList.size() != 0, "Boundary Edges not stored" ) ;
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
#endif
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setNumBEdges( UInt const n )
{
    _numBEdges = n;
#ifdef NOT_BDATA_FIRST

    _bEdges.reserve( n );
#endif
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::hasEdges() const
{
    return ! edgeList.empty();
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::hasInternalEdges() const
{
    return edgeList.size() > _numBEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numBEdges() const
{
    return _numBEdges;
}




template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryEdge( EdgeType const & e ) const
{
#ifdef NOT_BDATA_FIRST
    //ASSERT(false,"In this version Boundary edges must be stored first");
    bool isboundary = true;
    for ( UInt k = 1;k <= EdgeType::numVertices;++k )
    {
        isboundary = isboundary & e.point( k ).boundary();
    }
    return isboundary;
#else

    return e.id() <= _numBEdges;
#endif
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryEdge( ID const & id ) const
{
    return isBoundaryEdge( edge( id ) );
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isFullEdge( ID const & id ) const
{
    return edgeList.size() >= id;
}

template <typename GEOSHAPE, typename MC>
INLINE
UInt
RegionMesh2D<GEOSHAPE, MC>::edgeElement( ID const i, UInt const Pos ) const
{
    ASSERT_PRE( i <= edgeList.size(), "Not enough faces stored" ) ;
    ASSERT_BD( i > 0 ) ;
    return edgeElement( edge( i ), Pos );
};

template <typename GEOSHAPE, typename MC>
INLINE
UInt
RegionMesh2D<GEOSHAPE, MC>::edgeElement( EdgeType const & f, UInt const Pos ) const
{
    ASSERT_BD( ! edgeList.empty() ) ;
    ASSERT_PRE( Pos == 1 || Pos == 2 , "Wrong position (1 or 2)" ) ;
    ASSERT_BD( i > 0 ) ;
    if ( Pos == 1 )
    {
        return f.ad_first();
    }
    else
    {
        return f.ad_second();
    }
};

// ************************ Points/Vertices

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numPoints() const
{
    return _numPoints;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numPoints()
{
    return _numPoints;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::storedPoints() const
{
    return numItems( pointList );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::storedBPoints() const
{
    return _bPoints.size();
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::maxNumPoints() const
{
    return maxNumItems( pointList );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumPoints( UInt const n, bool const setcounter )
{
    setMaxNumItems( pointList, n, "Point" );
    if ( setcounter )
        _numPoints = n;
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::addPoint( bool const boundary, bool const vertex )
{
    return addPoint( PointType(), boundary, vertex );
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::addPoint
( PointType const & p, bool const boundary, bool const vertex )
{
    ASSERT_PRE( pointList.size() < pointList.capacity(), "Point list size exceeded" <<
                pointList.size() + 1 << " " << pointList.capacity() ) ;
    pointList.push_back( p );
    PointType * pp = & pointList.back();
    pp->id() = pointList.size();
    if ( boundary )
    {
        ASSERT_PRE( _bPoints.size() < _bPoints.capacity(), "Boundary Point list size exceeded" <<
                    _bPoints.size() + 1 << " " << _bPoints.capacity() ) ;
        _bPoints.push_back( pp );
        pp->boundary() = true;
    }
    return pointList.back();
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::setPoint
( PointType const & p, ID position, bool const boundary, bool const vertex )
{
    ASSERT_PRE( position <= pointList.capacity(), "Position  exceed lpoint list capacity" <<
                position << " " << pointList.capacity() ) ;
    bool found( false );
    pointList( position ) = p;
    PointType * pp = & pointList( position );
    pp->id() = position;
    if ( boundary )
    {
        pp->boundary() = true;
        // This is rather complex, since I do not know a priori
        // if point was already stored in the list!
        // No way to avoid it, sorry

        for ( SimpleVect<PointType *>::iterator bp = _bPoints.begin(); bp != _bPoints.end(); ++bp )
        {
            if ( ( *bp ) ->id() == position )
            {
                found = true;
                break;
            }
        }
        if ( ! found )
            _bPoints.push_back( pp );
    }
    return *pp;
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::lastPoint()
{
    return pointList.back();
}


template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::PointType const &
RegionMesh2D<GEOSHAPE, MC>::point( UInt const i ) const
{
    ASSERT_BD( i > 0 && i <= pointList.size() ) ;
    return pointList( i );
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::point( UInt const i )
{
    ASSERT_BD( i > 0 && i <= pointList.size() ) ;
    return pointList( i );
}


template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::PointType const &
RegionMesh2D<GEOSHAPE, MC>::boundaryPoint( ID const i ) const
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
INLINE
RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::boundaryPoint( ID const i )
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numBPoints() const
{
    return _numBPoints;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setNumBPoints( UInt const n )
{
    _numBPoints = n;
    _bPoints.reserve( n );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numVertices() const
{
    return _numVertices;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numVertices()
{
    return _numVertices;
}


template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numBVertices() const
{
    return _numBVertices;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numBVertices()
{
    return _numBVertices;
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isVertex( PointType const & p ) const
{
    return p.id() <= _numVertices;
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isVertex( ID const & id ) const
{
    return id <= _numVertices;
}


template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryPoint( ID const & id ) const
{
    return point( id ).boundary();
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryPoint( PointType const & p ) const
{
    return p.boundary();
}


/********************************************************************************
                        Element Adiacency Methods
*******************************************************************************/

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::hasLocalEdges() const
{
    return ! _FToE.empty();
}

#ifdef SAVEMEMORY
class BareEdge;

template <typename GEOSHAPE, typename MC>
INLINE
ID
RegionMesh3D<GEOSHAPE, MC>::localEdgeId( const FaceType & ifac, ID const locE ) const
{
    ASSERT_PRE( !_FToE.empty(), "Face to Edges array not  set" );
    ASSERT_BD( locE > 0 && locE <= FaceType::numLocalEdges );
    pair<BareEdge, bool> it;
    ID i1, i2;
    i1 = GEOSHAPE::eToP( locE, 1 );
    i2 = GEOSHAPE::eToP( locE, 2 );
    i1 = ( ifac.point( i1 ) ).id();
    i2 = ( ifac.point( i2 ) ).id();
    it = makeBareEdge( i1, i2 );
    return _FToE.id( it.first );
}

template <typename GEOSHAPE, typename MC>
INLINE
ID
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( ID const facId, ID const locE ) const
{
    ASSERT_BD( facId > 0 && facId <= _numFaces );
    return localEdgeId( face( facId ), locE );
}

#else

template <typename GEOSHAPE, typename MC>
INLINE
ID
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( const FaceType & ifac, ID const locE )
const
{
    return _FToE( ifac.id(), locE );
}

template <typename GEOSHAPE, typename MC>
INLINE
ID
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( ID const facId, ID const locE )
const
{
    ASSERT_PRE( !_FToE.empty(), "Face to Edges array not  set" );
    ASSERT_BD( facId > 0 && facId <= _numFaces );
    ASSERT_BD( locE > 0 && locE <= FaceType::numLocalEdges );
    return _FToE( locE, facId );
}

#endif

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::updateElementEdges()
{
#ifndef SAVEMEMORY
    BareItemsHandler<BareEdge> _be;
    pair<UInt, bool> e;
    _FToE.reshape( numLocalEdges(), numFaces() ); // DIMENSION ARRAY
#endif

    UInt vid, i1, i2;
    pair<BareEdge, bool> _edge;
    GEOSHAPE ele;
    // First We check if we have already Edges stored
    if ( ! edgeList.empty() )
    {
        // dump first edges, to maintain the correct numbering
        // if everything is correct the numbering in the bareedge
        // structure will reflect the actual edge numbering
        pair<UInt, bool> _check;
        for ( Edges::iterator j = edgeList.begin(); j != edgeList.end();++j )
        {
            i1 = ( j->point( 1 ) ).id();
            i2 = ( j->point( 2 ) ).id();
            _edge = makeBareEdge( i1, i2 );
#ifdef SAVEMEMORY

            _check = _FToE.addIfNotThere( _edge.first );
#else

            _check = _be.addIfNotThere( _edge.first );
#endif
#ifdef TEST_PRE
            // This precondition is hard to test otherwise!
            if ( !_check.second )
            {
                std::cerr << "Two identical Edges stored in EdgeList" << std::endl;
                std::cerr << "Vertex ids: " << i1 << ", " << i2 << std::endl;
                abort();
            }
            if ( _check.first != edgeList[ j ].id() )
            {
                std::cerr << "Edges in EdgeList have no correct id" << std::endl;
                std::cerr << "Stored Id: " << edgeList[ j ].id() << ", Found " << _check.first << std::endl;
                abort();
            }
#endif

        }
    }
    for ( Faces::iterator ifac = faceList.begin();
            ifac != faceList.end(); ++ifac )
    {
        vid = ifac->id();
        // REMEMBER: numbering from 1
        for ( UInt j = 1;j <= numLocalEdges();j++ )
        {
            i1 = ele.eToP( j, 1 );
            i2 = ele.eToP( j, 2 );
            // go to global
            i1 = ( ifac->point( i1 ) ).id();
            i2 = ( ifac->point( i2 ) ).id();
            _edge = makeBareEdge( i1, i2 );
#ifdef SAVEMEMORY

            _FToE.addIfNotThere( _edge.first );
#else

            e = _be.addIfNotThere( _edge.first );
            _FToE( j, vid ) = e.first;
#endif

        }
    }
#ifdef SAVEMEMORY
    UInt n = _FToE.maxId();
#else

    UInt n = _be.maxId();
#endif

    if ( _numEdges == 0 || _numEdges == _numBEdges )
        _numEdges = n;
    ASSERT_POS( n == _numEdges , "#Edges found not equal that in RegionMesh" << n << " " << _numEdges ) ;
    setLinkSwitch( std::string( "HAS_VOLUME_TO_EDGES" ) );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::cleanElementEdges()
{
    _FToE.clear();
    unsetLinkSwitch( "HAS_VOLUME_TO_EDGES" );
}


// *************** GENERAL INFO *******************

template <typename GEOSHAPE, typename MC>
std::ostream & RegionMesh2D<GEOSHAPE, MC>::showMe( bool verbose, std::ostream & out )
{
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << "                      RegionMesh2D                " << std::endl;
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << " ID: " << _id << std::endl;
    out << "Edges local to  faces stored: " << hasLocalEdges() << std::endl;
    //out <<"Edges local to  faces   stored:"<<switches.test("FACEtoEDGE")<<std::endl;
    //out <<"Faces adjacent to Edges stored: "<<switches.test("EDGEtoFACE")<<std::endl<<std::endl;
    out << "Edges Stored: " << hasEdges() << " Internal: "
    << hasInternalEdges() << std::endl;
    out << "**************************************************" << std::endl;
    out << "NumPoints=" << numPoints() << "  " << "numBPoints=" << numBPoints() << std::endl;
    out << "NumVertices=" << numVertices() << "  " << "numBVerices=" << numBVertices() << std::endl;
    out << "numFaces=" << numFaces() << std::endl;
    out << "numEdges=" << numEdges() << "NumBEdges=" << numBEdges() << std::endl;
    out << "**************************************************" << std::endl;
    switches.showMe( verbose, out );
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    if ( verbose )
    {
        std::cout << "Verbose version not implemented yet";
    }
    return out;

}
template <typename GEOSHAPE, typename MC>
int
RegionMesh2D<GEOSHAPE, MC>::check( int level, bool const fix, bool const verb, std::ostream & out )
{
    if ( verb )
    {
        out << "**************************************************" << std::endl;
        out << "         Checkin  RegionMesh2D                " << std::endl;
        out << " ID: " << _id << std::endl;
        out << "**************************************************" << std::endl;
    }
    if ( pointList.size() != _numPoints )
    {
        out << " Point list size " << pointList.size() << " not equal to internal counter value "
        << _numPoints << std::endl;
        if ( fix )
        {
            _numPoints = pointList.size();
            out << "Fixed";
            out.flush();
        }
    }
    if ( edgeList.size() == 0 )
        if ( verb )
            out << "Warning: No Edges Stored" << std::endl;

    if ( faceList.size() == 0 )
        if ( verb )
            out << "Warning: No Faces Stored" << std::endl;

    UInt count = 0;
    for ( Points::iterator i = pointList.begin(); i != pointList.end(); ++i )

        if ( i->boundary() )
            ++count;
    if ( count != _numBPoints )
    {
        out << " Num Boundary points " << count << " not equal to internal counter value "
        << _numBPoints << std::endl;
        if ( fix )
        {
            _numBPoints = count;
            out << "Fixed";
            out.flush();
        }
    }
    UInt badid = 0;
    for ( UInt i = 1; i <= storedPoints(); ++i )
        //for (UInt i=1; i<= _numVertices; ++i)// PERICOLOSISSIMO DOMANDARE ALEX
        if ( point( i ).id() != i )
            ++badid;
    if ( badid != 0 )
        out << " SEVERITY ERROR:" << badid << "Points ids are wrong";

    badid = 0;
    for ( UInt i = 1; i <= storedEdges(); ++i )
        if ( edge( i ).id() != i )
            ++badid;
    if ( badid != 0 )
        out << " SEVERITY ERROR:" << badid << "Edges ids are wrong";

    badid = 0;
    for ( UInt i = 1; i <= storedFaces(); ++i )
        if ( face( i ).id() != i )
            ++badid;
    if ( badid != 0 )
        out << " SEVERITY ERROR:" << badid << "Faces ids are wrong";

    badid = 0;

    if ( _numVertices == 0 )
        out << " SEVERITY ERROR: internal Vertices Counter unset";
    if ( _numPoints == 0 )
        out << " SEVERITY ERROR: internal Points Counter unset";
    if ( _numPoints == 0 )
        out << " SEVERITY ERROR: internal Points Counter unset";
    if ( _numBPoints == 0 )
        out << " SEVERITY ERROR: boundary Points Counter unset";
    if ( _numBVertices == 0 )
        out << " SEVERITY ERROR: boundary Vertices Counter unset";

    if ( verb )
        out << "   Check Finished              " << std::endl <<
        "***********************************************" << std::endl;
}
}
#endif
