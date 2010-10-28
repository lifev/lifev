//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
               2006-2010 EPFL, Politecnico di Milano

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing 1D mesh classes
 *
 *  @version 1.0
 *  @author Luca Formaggia
 *  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *  @author Tiziano Passerini
 *  @date
 *
 *  @version 1.1
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 27-04-2010
 */

#ifndef REGIONMESH1D_H
#define REGIONMESH1D_H

#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/switch.hpp>

#include <life/lifemesh/geoElement.hpp>
#include <life/lifemesh/bareItems.hpp>
#include <life/lifemesh/basisElSh.hpp>
#include <life/lifearray/SimpleVect.hpp>

namespace LifeV {
//! RegionMesh1D
/*!
 * @author Luca Formaggia
 * @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 * @author Tiziano Passerini
 *
 * This is the class that stores the mesh entities for a single 1D region In
 * a region elements are all of the same type
 */
template <typename GEOSHAPE, typename MC = DefMarkerCommon >
class RegionMesh1D : public MeshEntity,
                     public MC::RegionMarker
{
public:

    //! \name Markers_Types
    //@{

    typedef MC MarkerCommon;
    typedef typename MC::PointMarker     PointMarker;
    typedef typename MC::EdgeMarker      EdgeMarker;
    typedef typename MC::RegionMarker    RegionMarker;
    typedef typename MC::RegionMarker    Marker;

    //@}

    //! \name Basic_Element_Shapes_Types
    //@{

    typedef GEOSHAPE                     VolumeShape;
    typedef GEOSHAPE                     EdgeShape;
    typedef typename GEOSHAPE::GeoBShape PointShape;

    //@}

    //! \name Geometric_Element_Types
    //@{

    typedef GeoElement1D<EdgeShape, MC>  VolumeType;
    typedef GeoElement1D<EdgeShape, MC>  FaceType;
    typedef GeoElement1D<EdgeShape, MC>  EdgeType;
    typedef GeoElement0D<MC>             PointType;

    //@}

    //! \name GeoElement_Container_Types
    /*!
     * Typedefs for STL compliant containers of mesh geometric entities
     * I Use SimpleVect for addressing from 1.
     */
    //@{

    //! Points Container
    typedef SimpleVect<PointType>  Points;
    //! Elements Container - compatibility
    typedef SimpleVect<VolumeType> Volumes;
    //! Faces Container - compatibility
    typedef SimpleVect<FaceType>   Faces;
    //! Edges Container: at least boundary edges
    typedef SimpleVect<EdgeType>   Edges;

    //@}

    //! \name Generic_Types
    //@{

    typedef GEOSHAPE                     ElementShape;
    typedef typename GEOSHAPE::GeoBShape BElementShape;

    typedef GeoElement1D<GEOSHAPE, MC>   ElementType;
    typedef GeoElement0D<MC>             BElementType;

    typedef SimpleVect<EdgeType >        Elements;
    typedef SimpleVect<PointType>        BElements;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Default constructor
    /*!
     * @param id marker of the RegionMesh1D
     */
    explicit RegionMesh1D( ID id = 0 );

    //! Copy constructor
    /*!
     * @param mesh a RegionMesh1D
     */
    explicit RegionMesh1D( const RegionMesh1D<GEOSHAPE, MC>& mesh );

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
     * @param mesh a RegionMesh1D
     * @return the newly copied RegionMesh1D
     */
    RegionMesh1D<GEOSHAPE, MC> operator=( RegionMesh1D<GEOSHAPE, MC> const & mesh );

    //@}


    //! @name Methods
    //@{

    //! Setup mesh
    /*!
     * Mesh construction by hand.
     */
    void setup( const Real& Length, const UInt& NumberOfElements );

    //! Transform the mesh using boost::numeric::ublas; \author Cristiano Malossi:27/04/2010
    /*! Scale, rotate and translate the mesh (operations performed in this order!).
     *  NOTES:
     *    -  Rotation follows Paraview conventions: first rotate around z-axis,
     *       then around y-axis and finally around x-axis;
     *    -  All the vectors must allow the command: operator[];
     * @param scale        - vector of three components for (x,y,z) scaling of the mesh
     * @param rotate       - vector of three components (radiants) for rotating the mesh
     * @param translate    - vector of three components for (x,y,z) translation the mesh
     */
    template <typename VECTOR>
    void transformMesh( const VECTOR& scale, const VECTOR& rotate, const VECTOR& translate );

    //! Get the maximum H over all the edges of the mesh; @author Cristiano Malossi: 27/04/2010
    /*!
     * @return maximum H
     */
    Real maxH() const;


    //! Get the minumum H over all the edges of the mesh; @author Cristiano Malossi: 27/04/2010
    /*!
     * @return minumum H
     */
    Real minH() const;

    //! Get the mean H over all the edges of the mesh; @author Cristiano Malossi: 27/04/2010
    /*!
     * @return maximum H
     */
    Real meanH() const;

    //@}


    //! @name Switches Methods
    /*!
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

    //! numSwitches
    /*!
     * @return the number of switches which have been set
     */
    const UInt& numSwitches() const
    {
        return switches.size();
    }

    void set_switches_for_regionmesh1D( Switch & sw );

    /**
     * Interrogate Switch
     * @param _s name of the switch
     * @return true if the name is in the switch, false otherwise
     */
    bool getLinkSwitch( std::string const & _s ) const;

    /**
     *  Set a switch
     * @param _s
     */
    void setLinkSwitch( std::string const & _s );

    /**
     * unset a switch
     * @param _s
     */
    void unsetLinkSwitch( std::string const & _s );

    //@}


    //! @name Generic Methods
    /*!
     * These are the generic methods to get information about the number of
     * elements.  It is a generic interface common for all RegionMeshes (3D --
     * 1D)
     */
    //@{

    /**
     * Number of local vertices for each (1D) element
     * @return Number of local vertices for each (1D) element
     */
    UInt numLocalVertices() const;

    /**
     * Number of local edges for each (1D) element
     * @return Number of local edges for each (1D) element
     */
    //UInt numLocalEdges() const;

    /**
     * Number of elements in mesh (alias to numEdges())
     * @return Number of elements in mesh (alias to numEdges())
     */
    UInt numElements() const;

    /**
     * Number of 3D elements
     * @return Number of 3D elements
     */
    UInt & numElements();

    /**
     * Number of Global elements in mesh (alias to numEdges())
     * @return Number of elements in mesh (alias to numEdges())
     */
    UInt numGlobalElements() const;

    /**
     * Number of Global 3D elements
     * @return Number of 3D elements
     */
    UInt & numGlobalElements();

    /**
     * Number of Boundary edges
     * @return Number of Boundary edges
     */
    UInt numBElements() const;

    /**
     * Number of boundary edges
     * @return Number of boundary edges
     */
    UInt& numBElements();

    /**
     * get element at the i-th index
     * @param i index of the element
     * @return element at index i
     */
    ElementType& element( const ID& i );

    /**
     * get element at the i-th index
     * @param i index of the element
     * @return element at index i
     */
    const ElementType& element( const ID& i ) const;

    /**
     * get boundary element at the i-th index
     * @param i index of the boundary element
     * @return boundary element at the i-th index
     */
    BElementType& bElement( const ID& i );

    /**
     * get boundary element at the i-th index
     * @param i index of the boundary element
     * @return boundary element at the i-th index
     */
    const BElementType& bElement( const ID& i ) const;

    //@}


    //! \name Volume_Methods - For compatibility reasons this methods do the same as edge methods
    //@{

    UInt numVolumes()       const                { return 0; }
    UInt numGlobalVolumes() const                { return 0; }
    const VolumeType& volume( const ID i ) const { ASSERT_BD( i > 0 && i <= edgeList.size() ); return edgeList( i ); }
    VolumeType& volume( const ID i )             { ASSERT_BD( i > 0 && i <= edgeList.size() ); return edgeList( i ); }

    //@}


    //! \name Faces Methods - Not used in 1D
    //@{

    UInt numFaces()         const {return 0;}
    UInt numGlobalFaces()   const {return 0;}
    bool hasLocalFaces() const { return true; }
    ID localFaceId( const ID /*volId*/, const ID /*locF*/ ) const { return 0; }
    ID localFaceId( const VolumeType& /*iv*/, const ID /*locF*/ ) const {return 0;}
    void updateElementFaces( bool createFaces = false, UInt estimateFaceNumber = 0 );
    void cleanElementFaces() {}

    //@}

    /* ============================================
                      Edge Related Methods
       ============================================*/
    //! \name Edge_Methods All methods which operates on 1D elements
    //@{
    /**
     * Returns number of Edge elements in the mesh
     * as given by the internal counter
     * @return
     */
    UInt numEdges() const;
    UInt numGlobalEdges() const;


    /**
     * Access number of Edge (internal counter)
     * @return Access number of Edge (internal counter)
     */
    UInt & numEdges();

    /**
     * edges actually stored in list
     * @return number of edges actually stored in list
     */
    UInt storedEdges() const;

    /**
     *
     * @return
     */
    UInt maxNumEdges() const; /*!< Current capacity of Edges Container,
             i.e. how many elements may be stored */

    /**
     * Changes Current capacity of Edges
     * (Optionally sets internal counter.)
     * @param n maximum number of edges
     * @param setcounter true to set the counter, false otherwise
     */
    void setMaxNumEdges( UInt const n, bool const setcounter = false );

    /**
     * Adds edges. Id computed automatically.
     * @return ref to added edge
     */
    EdgeType & addEdge();

    /**
     * Adds edges. Id computed automatically.
     * @param v edge to add
     * @return reference to the newly added edge
     */
    EdgeType & addEdge( EdgeType const & v );

    /**
     * Add edge to a specified position
     * @param v edge to add
     * @param pos position of the edge
     * @return reference to the newly added edge
     */
    EdgeType & setEdge( EdgeType const & v, ID const pos );

    /**
     * set numEdges counter
     */
    void setEdgeCounter();

    /**
     * Reference to last edge stored in list.
     * Useful for mesh readers
     * @return reference of the last edge in the list
     */
    EdgeType & lastEdge();

    /**
     * ith mesh 1Delement
     * @param i index of the mesh 1Delement
     * @return the i-th edge
     */
    EdgeType const & edge( ID const i ) const;

    /**
     * ith mesh 1Delement
     * @param i index of the mesh edge
     * @return reference to the ith mesh edge
     */
    EdgeType & edge( ID const i ); //!<

    /**
     * ith mesh 1Delement length
     * @param i index of the mesh 1Delement
     * @return length of the i-th edge
     */
    const Real edgeLength( const ID& i ) const;

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
      the Edge methods look at the documentation of th eanalogous Edges
      methods */

    //@{
//     UInt numEdges() const; //!<Number of total edges in mesh (uses counter).
//     UInt & numEdges();    //!<Number of total edges (uses counter, may modify)
//     UInt storedEdges() const; //!< Number of stored edges
//     UInt maxNumEdges() const;  //!< Max number of edges that can be stored
     //! Set Maximum Number of Edges in the Edges container
//     void setMaxNumEdges( UInt const n, bool const setcounter = false );
    //! Adds a edge to list. Returns reference to it
    EdgeType & addEdge( bool const boundary = false );
    //! Adds a edge (optionally a boundary edge) to the end of the list
    //! and adjourn its ID. Returns reference to the added edge
    EdgeType & addEdge( EdgeType const & f, bool const boundary = false );
    //! Adds a edge (optionally a boundary edge) and adjourn its ID
    EdgeType & setEdge( EdgeType const & f, ID position, bool const boundary = false );
    //! The last edge in the container
//     EdgeType & lastEdge();
//     EdgeType const & edge( ID const i ) const; //<! ith mesh edge
//     EdgeType & edge( ID const i ); //!< ith mesh edge (may  modify!)
    EdgeType const & boundaryEdge( ID const i ) const; //!< ith boundary edge.
    EdgeType & boundaryEdge( ID const i ); //!< ith boundary  edge.
    void setNumBEdges( UInt const n ) ; //<! Set counter of boundary edges
    bool hasEdges() const;  //<! Do I store mesh edges?
    bool hasInternalEdges() const; //<! Do I store also internal edges?
    UInt numBEdges() const;  //<! Number of Boundary Edges
    bool isBoundaryEdge( EdgeType const & f ) const;  //<!Is this edge on boundary?
    bool isBoundaryEdge( ID const & id ) const;  //<!Is  edge whose ID is id on boundary?
    /*! Does this ID corresponds to a full edge? A FULL EDGE is a 1DElement  that is actually
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
    PointType & firstPoint(); //!< first mesh point
    PointType & lastPoint(); //!< last mesh point
    PointType const & point( ID const i ) const; //!< ith mesh point
    PointType & point( ID const i ); //!< ith mesh point
    PointType const & boundaryPoint( ID const i ) const; //!< ith b. point.
    PointType & boundaryPoint( ID const i ); //!< ith b. point.
    UInt numBPoints() const; //!< counter of boundary points
    void setNumBPoints( UInt const n ); //<! Sets counter of boundary points

    // Vertices

    UInt  numVertices () const; //!< Number of vertices in Region
    UInt& numVertices (); //!< Allows to change number of vertices in Region
    UInt  numBVertices() const; //!< Number of boundary vertices in RegionMesh
    UInt& numBVertices();

    UInt numGlobalVertices() const;
    void setNumGlobalVertices( UInt const n ){ M_numGlobalVertices = n; }

    bool isVertex        ( ID const & id )       const;  //<!Is this point a Vertex?
    bool isVertex        ( PointType const & p ) const;  //<!Is this point a Vertex?
    bool isBoundaryPoint ( PointType const & p ) const;  //<!Is this point on boundary?
    bool isBoundaryPoint ( ID const & id )       const;  //<!Is this point on boundary?
    //@}


    //! Prints some mesh info
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout );
    //! Basic tests for mesh consistency.
    /*! For more estensive test see \link mesh_util.h */
    int check( int level = 0, bool const fix = false, bool const verbose = true, std::ostream & out = std::cerr );


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
    bool hasLocalEdges() const {return true;}  //!<Check wether edge info is available
    //!Edge n. locF around face. Returns edgeID
    UInt localEdgeId( const FaceType& /*iface*/, const UInt /*locE*/ ) const {return 0;}
    //!Edge n. locF around face. Returns edgeID
    UInt localEdgeId( const ID /*facId*/, const UInt /*locE*/ ) const {return 0;}
    //
    void updateElementEdges(bool ce=false, UInt ee=0 ); //!<Build localEdgeId table
    //! Destroys element-to-face container
    void cleanElementEdges() {return;}

    //
    void setNumGlobalEdges( UInt const n ){ M_numGlobalEdges = n; }

    /*! \name RegionContainers
     *

     * Stl compliant containers for basic structure I expose them since they
     * are standard containers
     */
    //@{
    Volumes volumeList; // NOT USED IN 1D
    Points  pointList;  //!< Container of mesh points/verices
    Faces   faceList;   // NOT USED IN 1D
    Edges   edgeList;   //!< Container of mesh Edges.

    SimpleVect<PointType * > _bPoints; //!< Boundary points list
    //@}

    //! Switches
    Switch switches;

protected:

    template < typename T >
    UInt numItems( SimpleVect< T> const & list ) const;

    template < typename T >
    UInt maxNumItems( SimpleVect< T> const & list ) const;

    template < typename T >
    void setMaxNumItems( SimpleVect< T> & list, UInt n, std::string title );

    /*! Arrays containing the ids of Edges and Edges of each element
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

    // Internal counters
    UInt _numVertices;
    UInt _numBVertices;
    UInt _numPoints;
    UInt _numBPoints;
    UInt _numEdges;
    UInt _numBEdges;

    UInt M_numGlobalVertices;
    UInt M_numGlobalPoints;
    UInt M_numGlobalEdges;
};

// ===================================================
// Constructors
// ===================================================
template <typename GEOSHAPE, typename MC>
RegionMesh1D<GEOSHAPE, MC>::RegionMesh1D( ID id ) :
        MeshEntity          ( id ),
        MC::RegionMarker    (),
        switches            (),
        _numVertices        ( 0 ),
        _numBVertices       ( 0 ),
        _numPoints          ( 0 ),
        _numBPoints         ( 0 ),
        _numEdges           ( 0 ),
        _numBEdges          ( 0 ),
        M_numGlobalVertices (),
        M_numGlobalPoints   (),
        M_numGlobalEdges    ()
{
    set_switches_for_regionmesh1D( switches );
}

template <typename GEOSHAPE, typename MC>
RegionMesh1D<GEOSHAPE, MC>::RegionMesh1D( RegionMesh1D<GEOSHAPE, MC> const & m )
{
    ASSERT( true, "Copy Costructor Not Yet Implemented for RegionMesh1D" );
}

template <typename GEOSHAPE, typename MC>
RegionMesh1D<GEOSHAPE, MC>
RegionMesh1D<GEOSHAPE, MC>::operator=( RegionMesh1D<GEOSHAPE, MC> const & m )
{
    ASSERT( true, "Assignement Operator  Yet Not Implemented for RegionMesh1D" );
}

// ===================================================
// Methods
// ===================================================
template <typename GEOSHAPE, typename MC>
void
RegionMesh1D<GEOSHAPE, MC>::setup( const Real& Length, const UInt& NumberOfElements )
{
    ASSERT_PRE( NumberOfElements > 0, "BasicOneDMesh::BasicOneDMesh problems! The number of elements must be positive!");

    Real x_current = 0;

    setMaxNumPoints(NumberOfElements + 1, true);
    setNumBPoints(2);

    Real deltax = Length / NumberOfElements;
    ASSERT_PRE( deltax > 0 ,
                "BasicOneDMesh::BasicOneDMesh problems! The left point is on the right..." );

    PointType * pp = 0;

    for (UInt it = 0; it < NumberOfElements + 1; it++)
    {
        // insert a new Point1D in point list
        pp = &addPoint( (it == NumberOfElements) || ( it == 0) );
        pp->x() = x_current;
        pp->setId(it + 1);
        pp->setLocalId(it + 1);
        // move to the next point location
        x_current += deltax;
        //std::cout << "point " << pointList.back().id() << " " << pointList.back().x() << std::endl;
    }

    ASSERT_PRE( std::abs(Length - x_current + deltax ) < 1e-10 * deltax ,
                "BasicOneDMesh::BasicOneDMesh problems! Check xleft<xright?" );

    setMaxNumEdges(NumberOfElements, true);

    setNumGlobalVertices(pointList.size());
    _numVertices = pointList.size();

    EdgeType* pe = 0;

    for (UInt it = 0; it < NumberOfElements; it++)
    {
        pe = &addEdge( false );
        pe->setPoint(1, pointList(it + 1));
        pe->setPoint(2, pointList(it + 2));
        pe->setId(it + 1);
        pe->setLocalId(it + 1);
    }

    setEdgeCounter();
    setNumGlobalEdges(edgeList.size());
}

template <typename GEOSHAPE, typename MC>
template <typename VECTOR>
void RegionMesh1D<GEOSHAPE, MC>::transformMesh( const VECTOR& scale, const VECTOR& rotate, const VECTOR& translate )
{
    //Create the 3 planar rotation matrix and the scale matrix
    boost::numeric::ublas::matrix<Real> R(3,3), R1(3,3), R2(3,3), R3(3,3), S(3,3);

    R1(0,0) =  1.;             R1(0,1) =  0.;             R1(0,2) =  0.;
    R1(1,0) =  0.;             R1(1,1) =  cos(rotate[0]); R1(1,2) = -sin(rotate[0]);
    R1(2,0) =  0.;             R1(2,1) =  sin(rotate[0]); R1(2,2) =  cos(rotate[0]);

    R2(0,0) =  cos(rotate[1]); R2(0,1) =  0.;             R2(0,2) =  sin(rotate[1]);
    R2(1,0) =  0.;             R2(1,1) =  1.;             R2(1,2) = 0.;
    R2(2,0) = -sin(rotate[1]); R2(2,1) =  0.;             R2(2,2) =  cos(rotate[1]);

    R3(0,0) =  cos(rotate[2]); R3(0,1) = -sin(rotate[2]); R3(0,2) = 0.;
    R3(1,0) =  sin(rotate[2]); R3(1,1) =  cos(rotate[2]); R3(1,2) = 0.;
    R3(2,0) =  0;              R3(2,1) =  0.;             R3(2,2) = 1.;

     S(0,0) = scale[0];         S(0,1) = 0.;               S(0,2) = 0.;
     S(1,0) = 0.;               S(1,1) = scale[1];         S(1,2) = 0.;
     S(2,0) = 0.;               S(2,1) = 0.;               S(2,2) = scale[2];

    //The total rotation is: R = R1*R2*R3 (as in Paraview we rotate first around z, then around y, and finally around x).
    //We also post-multiply by S to apply the scale before the rotation.
    R = prod( R3, S );
    R = prod( R2, R );
    R = prod( R1, R );

    //Create the 3D translate vector
    boost::numeric::ublas::vector<Real> P(3), T(3);
    T(0) = translate[0]; T(1) = translate[1];  T(2) = translate[2];

    //Apply the transformation
    for ( UInt i(0); i < pointList.size(); ++i )
    {
        //P = pointList[ i ].coordinate(); // Try to avoid double copy if possible

        P( 0 ) = pointList[ i ].coordinate( 1 );
        P( 1 ) = pointList[ i ].coordinate( 2 );
        P( 2 ) = pointList[ i ].coordinate( 3 );

        P = T + prod( R, P );

        pointList[ i ].coordinate( 1 ) = P( 0 );
        pointList[ i ].coordinate( 2 ) = P( 1 );
        pointList[ i ].coordinate( 3 ) = P( 2 );
    }
}

template <typename GEOSHAPE, typename MC>
Real RegionMesh1D<GEOSHAPE, MC>::maxH() const
{
    Real MaxH(0);
    Real deltaX(0), deltaY(0), deltaZ(0);

    for ( UInt i(0); i < static_cast<UInt> ( edgeList.size() ); ++i )
    {
        deltaX = ( edgeList( i+1 ).point( 2 ) ).x() - ( edgeList( i+1 ).point( 1 ) ).x();
        deltaY = ( edgeList( i+1 ).point( 2 ) ).y() - ( edgeList( i+1 ).point( 1 ) ).y();
        deltaZ = ( edgeList( i+1 ).point( 2 ) ).z() - ( edgeList( i+1 ).point( 1 ) ).z();

        deltaX *= deltaX;
        deltaY *= deltaY;
        deltaZ *= deltaZ;

        MaxH = std::max( MaxH, deltaX+deltaY+deltaZ );
    }

    return std::sqrt( MaxH );
}

template <typename GEOSHAPE, typename MC>
Real RegionMesh1D<GEOSHAPE, MC>::minH() const
{
    Real MinH( 1E10 );
    Real deltaX(0), deltaY(0), deltaZ(0);

    for ( UInt i(0); i < static_cast<UInt> ( edgeList.size() ); ++i )
    {
        deltaX = ( edgeList( i+1 ).point( 2 ) ).x() - ( edgeList( i+1 ).point( 1 ) ).x();
        deltaY = ( edgeList( i+1 ).point( 2 ) ).y() - ( edgeList( i+1 ).point( 1 ) ).y();
        deltaZ = ( edgeList( i+1 ).point( 2 ) ).z() - ( edgeList( i+1 ).point( 1 ) ).z();

        deltaX *= deltaX;
        deltaY *= deltaY;
        deltaZ *= deltaZ;

        MinH = std::min( MinH, deltaX+deltaY+deltaZ );
    }

    return std::sqrt( MinH );
}

template <typename GEOSHAPE, typename MC>
Real RegionMesh1D<GEOSHAPE, MC>::meanH() const
{
    Real MeanH = 0;
    Real deltaX(0), deltaY(0), deltaZ(0);

    for ( UInt i(0); i < static_cast<UInt> ( edgeList.size() ); ++i )
    {
        deltaX = ( edgeList( i+1 ).point( 2 ) ).x() - ( edgeList( i+1 ).point( 1 ) ).x();
        deltaY = ( edgeList( i+1 ).point( 2 ) ).y() - ( edgeList( i+1 ).point( 1 ) ).y();
        deltaZ = ( edgeList( i+1 ).point( 2 ) ).z() - ( edgeList( i+1 ).point( 1 ) ).z();

        deltaX *= deltaX;
        deltaY *= deltaY;
        deltaZ *= deltaZ;

        MeanH += deltaX+deltaY+deltaZ;
    }

    return std::sqrt( MeanH / static_cast<Real> ( edgeList.size() ) );
}

// ===================================================
// Switch Methods
// ===================================================
template <typename GEOSHAPE, typename MC>
void
RegionMesh1D<GEOSHAPE, MC>::set_switches_for_regionmesh1D( Switch & sw )
{
    sw.create( "HAS_ALL_EDGES" );
    sw.create( "HAS_BOUNDARY_EDGES" );
    sw.create( "HAS_FACE_TO_EDGES" );
    sw.create( "HAS_BEEN_CHECKED" );
    sw.create( "EDGES_HAVE_ADIACENCY" );
}

template <typename GEOSHAPE, typename MC>
INLINE
void
RegionMesh1D<GEOSHAPE, MC>::setLinkSwitch( std::string const & _s )
{
    std::ostringstream _err_msg;
	_err_msg << "Switch named " << _s << " is not allowed";
    ASSERT0( switches.set( _s ), _err_msg.str().c_str() );
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::getLinkSwitch( std::string const & _s ) const
{
    return switches.test( _s );
}

template <typename GEOSHAPE, typename MC>
INLINE
void
RegionMesh1D<GEOSHAPE, MC>::unsetLinkSwitch( std::string const & _s )
{
    std::ostringstream _err_msg;
	_err_msg << "Switch named " << _s << " is not allowed";
    ASSERT0( switches.unset( _s ), _err_msg.str().c_str() );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::numLocalVertices() const
{
    return EdgeType::numLocalVertices;
}

// template <typename GEOSHAPE, typename MC>
// UInt
// RegionMesh1D<GEOSHAPE, MC>::numLocalEdges() const
// {
//     return EdgeType::numLocalEdges;
// }

// ************** Generic Methods
template <typename GEOSHAPE, typename MC>
UInt RegionMesh1D<GEOSHAPE, MC>::numElements() const
{
    return _numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt & RegionMesh1D<GEOSHAPE, MC>::numElements()
{
    return _numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt RegionMesh1D<GEOSHAPE, MC>::numGlobalElements() const
{
    return _numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt & RegionMesh1D<GEOSHAPE, MC>::numGlobalElements()
{
    return _numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt RegionMesh1D<GEOSHAPE, MC>::numBElements() const
{
    return _numBEdges;
}

template <typename GEOSHAPE, typename MC>
UInt &RegionMesh1D<GEOSHAPE, MC>::numBElements()
{
    return _numBEdges;
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh1D<GEOSHAPE, MC>::ElementType &
RegionMesh1D<GEOSHAPE, MC>::element( ID const & i )
{
    return edge( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh1D<GEOSHAPE, MC>::ElementType const &
RegionMesh1D<GEOSHAPE, MC>::element( ID const & i ) const
{
    return edge( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh1D<GEOSHAPE, MC>::BElementType &
RegionMesh1D<GEOSHAPE, MC>::bElement( ID const & i )
{
    return boundaryEdge( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh1D<GEOSHAPE, MC>::BElementType const &
RegionMesh1D<GEOSHAPE, MC>::bElement( ID const & i ) const
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
RegionMesh1D<GEOSHAPE, MC>::numItems( SimpleVect< T> const & list ) const
{
    return list.size();
}

template <typename GEOSHAPE, typename MC>
template <typename T>
UInt
RegionMesh1D<GEOSHAPE, MC>::maxNumItems( SimpleVect< T> const & list ) const
{
    return list.capacity();
}

template <typename GEOSHAPE, typename MC>
template <typename T>
void
RegionMesh1D<GEOSHAPE, MC>::setMaxNumItems( SimpleVect< T> & list, UInt n, std::string title )
{
    if ( list.capacity() == 0 )
    {
        list.reserve( n );
    }
    else if ( list.capacity() == n )
    {
        Warning() << "Capacity of " << title << "list already set to " << n << "\n";
    }
    else
    {
        Warning() << "Resetting " << title << " list size to " << n << "\n";
        Warning() << "ALL PREVIOUS POINTERS TO THE LIST (IF ANY) ARE NOW INVALID\n";

        list.reserve( n );
    }
}
//-------------------------------------------------------------------

// ***************************** EDGES



// \todo use addItem
template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::addEdge()
{
    return addEdge( EdgeType() );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::addEdge( EdgeType const & v )
{
    ASSERT_PRE( edgeList.size() < edgeList.capacity() , "Edge list size exceeded" <<
                edgeList.size() + 1 << " " << edgeList.capacity() ) ;
    edgeList.push_back( v );
    ( edgeList.back() ).id() = edgeList.size();
    return edgeList.back();
}
// \todo Use setItem

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::setEdge( EdgeType const & v, ID const pos )
{
    ASSERT_PRE( pos <= edgeList.capacity() , "position requested exceed capacity" <<
                pos << " " << edgeList.capacity() ) ;
    edgeList( pos ) = v;
    edgeList( pos ).id() = pos;
    return edgeList( pos );
}

template <typename GEOSHAPE, typename MC>
INLINE
void
RegionMesh1D<GEOSHAPE, MC>::setEdgeCounter()
{
    _numEdges = edgeList.size();
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType const &
RegionMesh1D<GEOSHAPE, MC>::edge( ID const i ) const
{
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::edge( ID const i )
{
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
}

template <typename GEOSHAPE, typename MC>
const Real
RegionMesh1D<GEOSHAPE, MC>::edgeLength( const ID& i ) const
{
    ASSERT_BD( i > 0 && i <= edgeList.size() );

    Real deltaX, deltaY, deltaZ;

    deltaX = ( edgeList( i+1 ).point( 2 ) ).x() - ( edgeList( i+1 ).point( 1 ) ).x();
    deltaY = ( edgeList( i+1 ).point( 2 ) ).y() - ( edgeList( i+1 ).point( 1 ) ).y();
    deltaZ = ( edgeList( i+1 ).point( 2 ) ).z() - ( edgeList( i+1 ).point( 1 ) ).z();

    deltaX *= deltaX;
    deltaY *= deltaY;
    deltaZ *= deltaZ;

    return std::sqrt( deltaX+deltaY+deltaZ );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::numEdges() const
{
    return _numEdges;
}


template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::numGlobalEdges() const
{
    return M_numGlobalEdges;
}



template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh1D<GEOSHAPE, MC>::numEdges()
{
    return _numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::storedEdges() const
{
    return numItems( edgeList );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::maxNumEdges() const
{
    return maxNumItems( edgeList );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh1D<GEOSHAPE, MC>::setMaxNumEdges( UInt const n, bool const setcounter )
{
    setMaxNumItems( edgeList, n, "Edge" );
    if ( setcounter )
        _numEdges = n;
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::addEdge( bool const boundary )
{
    return addEdge( EdgeType(), boundary );
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::addEdge( EdgeType const & f, bool const boundary )
{
    ASSERT_PRE( edgeList.size() < edgeList.capacity(), "Edge list size exceeded" <<
                edgeList.size() + 1 << " " << edgeList.capacity() ) ;

    edgeList.push_back( f );

    (edgeList.back()).setId(edgeList.size());
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
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::setEdge( EdgeType const & f, ID position, bool const boundary )
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
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::lastEdge()
{
    return edgeList.back();
}




template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType const &
RegionMesh1D<GEOSHAPE, MC>::boundaryEdge( ID const i ) const
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
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::boundaryEdge( ID const i )
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
RegionMesh1D<GEOSHAPE, MC>::setNumBEdges( UInt const n )
{
    _numBEdges = n;
#ifdef NOT_BDATA_FIRST

    _bEdges.reserve( n );
#endif
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::hasEdges() const
{
    return ! edgeList.empty();
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::hasInternalEdges() const
{
    return edgeList.size() > _numBEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::numBEdges() const
{
    return _numBEdges;
}




template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::isBoundaryEdge( EdgeType const & e ) const
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
RegionMesh1D<GEOSHAPE, MC>::isBoundaryEdge( ID const & id ) const
{
    return isBoundaryEdge( edge( id ) );
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::isFullEdge( ID const & id ) const
{
    return edgeList.size() >= id;
}

template <typename GEOSHAPE, typename MC>
INLINE
UInt
RegionMesh1D<GEOSHAPE, MC>::edgeElement( ID const i, UInt const Pos ) const
{
    ASSERT_PRE( i <= edgeList.size(), "Not enough edges stored" ) ;
    ASSERT_BD( i > 0 ) ;
    return edgeElement( edge( i ), Pos );
};

template <typename GEOSHAPE, typename MC>
INLINE
UInt
RegionMesh1D<GEOSHAPE, MC>::edgeElement( EdgeType const & f, UInt const Pos ) const
{
    ASSERT_BD( ! edgeList.empty() ) ;
    ASSERT_PRE( Pos == 1 || Pos == 2 , "Wrong position (1 or 2)" ) ;
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
RegionMesh1D<GEOSHAPE, MC>::numPoints() const
{
    return _numPoints;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh1D<GEOSHAPE, MC>::numPoints()
{
    return _numPoints;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::storedPoints() const
{
    return numItems( pointList );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::storedBPoints() const
{
    return _bPoints.size();
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::maxNumPoints() const
{
    return maxNumItems( pointList );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh1D<GEOSHAPE, MC>::setMaxNumPoints( UInt const n, bool const setcounter )
{
    setMaxNumItems( pointList, n, "Point" );
    if ( setcounter )
        _numPoints = n;
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::PointType &
RegionMesh1D<GEOSHAPE, MC>::addPoint( bool const boundary, bool const vertex )
{
    return addPoint( PointType(), boundary, vertex );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::PointType &
RegionMesh1D<GEOSHAPE, MC>::addPoint( PointType const & p, bool const boundary, bool const /*vertex*/ )
{
    ASSERT_PRE( pointList.size() < pointList.capacity(), "Point list size exceeded" <<
                pointList.size() + 1 << " " << pointList.capacity() ) ;
    pointList.push_back( p );
    PointType * pp = & pointList.back();
    pp->setId(pointList.size());
    if ( boundary )
    {
        ASSERT_PRE( _bPoints.size() < _bPoints.capacity(), "Boundary Point list size exceeded" <<
                    _bPoints.size() + 1 << " " << _bPoints.capacity() ) ;
        _bPoints.push_back( pp );
        pp->setBoundary(true);
    }
    return pointList.back();
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::PointType &
RegionMesh1D<GEOSHAPE, MC>::setPoint
( PointType const & p, ID position, bool const boundary, bool const vertex )
{
    ASSERT_PRE( position <= pointList.capacity(), "Position  exceed lpoint list capacity" <<
                position << " " << pointList.capacity() ) ;
    bool found( false );
    pointList( position ) = p;
    PointType * pp = & pointList( position );
    pp->setId(position);
    if ( boundary )
    {
        pp->boundary() = true;
        // This is rather complex, since I do not know a priori
        // if point was already stored in the list!
        // No way to avoid it, sorry

        for ( typename SimpleVect<PointType *>::iterator bp = _bPoints.begin(); bp != _bPoints.end(); ++bp )
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
typename RegionMesh1D<GEOSHAPE, MC>::PointType &
RegionMesh1D<GEOSHAPE, MC>::firstPoint()
{
    return pointList.front();
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::PointType &
RegionMesh1D<GEOSHAPE, MC>::lastPoint()
{
    return pointList.back();
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::PointType const &
RegionMesh1D<GEOSHAPE, MC>::point( UInt const i ) const
{
    ASSERT_BD( i > 0 && i <= pointList.size() ) ;
    return pointList( i );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::PointType &
RegionMesh1D<GEOSHAPE, MC>::point( UInt const i )
{
    ASSERT_BD( i > 0 && i <= pointList.size() ) ;
    return pointList( i );
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::PointType const &
RegionMesh1D<GEOSHAPE, MC>::boundaryPoint( ID const i ) const
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh1D<GEOSHAPE, MC>::PointType &
RegionMesh1D<GEOSHAPE, MC>::boundaryPoint( ID const i )
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::numBPoints() const
{
    return _numBPoints;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh1D<GEOSHAPE, MC>::setNumBPoints( UInt const n )
{
    _numBPoints = n;
    _bPoints.reserve( n );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::numVertices() const
{
    return _numVertices;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh1D<GEOSHAPE, MC>::numVertices()
{
    return _numVertices;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::numGlobalVertices() const
{
    return M_numGlobalVertices;
}


template <typename GEOSHAPE, typename MC>
UInt
RegionMesh1D<GEOSHAPE, MC>::numBVertices() const
{
    return _numBVertices;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh1D<GEOSHAPE, MC>::numBVertices()
{
    return _numBVertices;
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::isVertex( PointType const & p ) const
{
    return p.id() <= _numVertices;
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::isVertex( ID const & id ) const
{
    return id <= _numVertices;
}


template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::isBoundaryPoint( ID const & id ) const
{
    return point( id ).boundary();
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh1D<GEOSHAPE, MC>::isBoundaryPoint( PointType const & p ) const
{
    return p.boundary();
}


//


template <typename GEOSHAPE, typename MC>
void
RegionMesh1D<GEOSHAPE, MC>::updateElementEdges( bool /*ce*/, UInt /*ee*/ )
{
    return;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh1D<GEOSHAPE, MC>::updateElementFaces( bool /*cf*/, UInt /*ef*/ )
{
    return;
}

// *************** GENERAL INFO *******************

template <typename GEOSHAPE, typename MC>
std::ostream & RegionMesh1D<GEOSHAPE, MC>::showMe( bool verbose, std::ostream & out )
{
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << "                      RegionMesh1D                " << std::endl;
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << " ID: " << this->id() << std::endl;    //out <<"Edges local to  edges   stored:"<<switches.test("EDGEtoEDGE")<<std::endl;
    //out <<"Edges adjacent to Edges stored: "<<switches.test("EDGEtoEDGE")<<std::endl<<std::endl;
    out << "Edges Stored: " << hasEdges() << " Internal: "
    << hasInternalEdges() << std::endl;
    out << "**************************************************" << std::endl;
    out << "NumPoints=" << numPoints() << "  " << "numBPoints=" << numBPoints() << std::endl;
    out << "NumVertices=" << numVertices() << "  " << "numBVerices=" << numBVertices() << std::endl;
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
RegionMesh1D<GEOSHAPE, MC>::check( int level, bool const fix, bool const verb, std::ostream & out )
{
    int severity = 0;

    if ( verb )
    {
        out << "**************************************************" << std::endl;
        out << "         Checkin  RegionMesh1D                " << std::endl;
        out << " ID: " << this->id() << std::endl;
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
    {
        if ( verb )
            out << "Warning: No Edges Stored" << std::endl;
        if ( verb )
            out << "MAY NOT WORK IF WE HAVE DOF ON EDGES AND ESSENTIAL BC!!!!" << std::endl;
        severity = -1;
        unsetLinkSwitch( "HAS_ALL_EDGES" );
        unsetLinkSwitch( "HAS_BOUNDARY_EDGES" );
    }

    UInt count = 0;
    for ( typename Points::iterator i = pointList.begin(); i != pointList.end(); ++i )

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
    for ( UInt i = 1; i <= storedEdges(); ++i )
        if ( edge( i ).id() != i )
            ++badid;
    if ( badid != 0 )
        out << " SEVERITY ERROR:" << badid << "Edges ids are wrong";

    badid = 0;

    if ( _numVertices == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: internal Vertices Counter unset";
    }
    if ( _numPoints == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: internal Points Counter unset";
    }
    if ( _numPoints == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: internal Points Counter unset";
    }
    if ( _numBPoints == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: boundary Points Counter unset";
    }
    if ( _numBVertices == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: boundary Vertices Counter unset";
    }
    if ( verb )
        out << "   Check Finished              " << std::endl <<
        "***********************************************" << std::endl;

    return severity;
}

}

#endif //REGIONMESH1D_H
