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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with LifeV. If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing 2D Mesh Classes
 *
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *
 *  @contributor Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 *  @mantainer Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 */

#ifndef _REGIONMESH2D_HH_
#define _REGIONMESH2D_HH_

#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <life/lifecore/life.hpp>
#include <life/lifecore/switch.hpp>

#include <life/lifemesh/geoElement.hpp>
#include <life/lifemesh/bareItems.hpp>
#include <life/lifemesh/basisElSh.hpp>
#include <life/lifearray/SimpleVect.hpp>

namespace LifeV
{

/**
 *  @class RegionMesh2D
 *  @brief Class for 2D Mesh
 * 
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *
 *  This is the class that stores the mesh entities for a single 2D region.
 *
 *  In a region elements are all of the same type.
 */
template <typename GEOSHAPE, typename MC = DefMarkerCommon >
class RegionMesh2D
        :
        public MeshEntity,
        public MC::RegionMarker
{

public:
    /** @name Marker Types
     *  @ingroup public_types
     *  Markers for Point, Edge, Face and Region.
     * 
     *  @{
     */
    
    //! Point Marker
    typedef typename MC::PointMarker PointMarker;
    //! Edge Marker
    typedef typename MC::EdgeMarker EdgeMarker;
    //! Face Marker
    typedef typename MC::FaceMarker FaceMarker;
    //! Region Marker
    typedef typename MC::RegionMarker RegionMarker;
    //! Region Marker
    typedef typename MC::RegionMarker Marker;

    /** @} */ // End of group Marker Types

    /** @name Basic Element Shape Types
     *  @ingroup public_types
     *  Face and Edge geometric shapes.
     *  @{
     */

    //! Face Shape.
    typedef GEOSHAPE FaceShape;
    //! Edge Shape (Boundary Element).
    typedef typename GEOSHAPE::GeoBShape EdgeShape;

    /** @} */ // End of group Basic Element Shape Types

    /** @name Geometric Element Types
     *  @ingroup public_types
     *  Faces, Edges and Points.
     *  @{
     */

    //! Face Element (2D)
    typedef GeoElement2D<GEOSHAPE, MC> FaceType;
    //! Edge Element (1D)
    typedef GeoElement1D<EdgeShape, MC> EdgeType;
    //! Point Element (0D)
    typedef GeoElement0D<MC> PointType;

    /** @} */ // End of group Geometric Element Types

    /** @name Geometric Element Container Types
     *  @ingroup public_types
     *  Typedefs for STL compliant containers of mesh geometric entities.
     *
     *  I Use SimpleVect container for addressing from 1.
     *  @{
     */

    //! Points Container
    typedef SimpleVect<PointType> Points;
    //! Faces Container
    typedef SimpleVect<FaceType> Faces;
    //! Edges Container: at least boundary edges
    typedef SimpleVect<EdgeType> Edges;

    /** @} */ // End of group Geometric Element Container Types


    /** @name Generic Types
     *  @ingroup public_types
     * 
     *  Generic types for all regionMeshes.
     * 
     *  These are part of the generic generic interface common for all RegionMeshes (3D -- 1D).
     * 
     *  @{
     */
    
    //! Element Geometric Shape
    typedef GEOSHAPE                     ElementShape;
    //! Boundary Element Geometric Shape
    typedef typename GEOSHAPE::GeoBShape BElementShape;

    //! Element Geometric Type
    typedef GeoElement2D<GEOSHAPE, MC>   ElementType;
    //! Boundary Element Geometric Type
    typedef GeoElement1D<EdgeShape, MC>  BElementType;

    //! Element Geometric Shape Container Type
    typedef SimpleVect<FaceType >        Elements;
    //! Boundary Element Geometric Shape Container Type
    typedef SimpleVect<EdgeType>         BElements;

    /** @} */ // End of group Generic Types


    /** @name Constructors & Destructor
     *  Default and Copy Constructor for the class.
     *  @{
     */

    //! Default constructor
    explicit RegionMesh2D();
    
    //! Default constructor
    /**
     *  @param id marker of the RegionMesh2D
     */
    explicit RegionMesh2D( UInt id );

    //! Copy constructor
    /**
     *  @param m a RegionMesh2D
     */
    explicit RegionMesh2D( RegionMesh2D<GEOSHAPE, MC> const & m );

    //! Destructor
    ~RegionMesh2D<GEOSHAPE, MC>();
     
    /** @} */ // End of group Constructors & Destructor


    /** @name Operators
     *  Public Operator Methods
     *
     *  @{
     */

    //! Assignment operator
    /**
     *  @param m a RegionMesh2D
     *  @return the newly copied RegionMesh2D
     */
    RegionMesh2D<GEOSHAPE, MC> operator=( RegionMesh2D<GEOSHAPE, MC> const & m );

    /** @} */ // End of group Operators


    /** @defgroup public_methods Public Methods
     *
     */

    /** @name Debugging Methods
     *  @ingroup public_methods
     *  Debugging methods
     *
     *  @{
     */

    //! Display general information about the content of the class.
    /**
     *  Pretty output of information about the mesh.
     *
     *  @param verbose If true output is verbose, false otherwise (default);
     *  @param out Output stream (std::cout default);
     *  @return Output stream (for concatenation).
     */
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout );

    //! Basic tests for mesh consistency.
    /**
     *  Check consistency of the mesh and fix errors.
     *
     *  @param level How much deep in the check we should go (not used yet).
     *  @param fix If true errors will be fixed, not otherwise (default).
     *  @param verbose If true output is verbose, false otherwise (default);
     *  @param out Output stream (std::cerr default);
     *  @return Severity of the errors.
     */
    int check( int level = 0, bool const fix = false, bool const verbose = true, std::ostream & out = std::cerr );

    /** @} */ // End of group Debugging Methods


    /** @name Switches Methods
     *  @ingroup public_methods
     *  Switches are used to store the status of the RegionMesh The switches
     *  are used internally to control whether some data structures have been
     *  set up.
     *
     *  The possible Switches are:
     *  - \c HAS_ALL_EDGES
     *  - \c HAS_FACE_TO_EDGES
     *  - \c HAS_BEEN_CHECKED
     *  - \c HAS_BOUNDARY_EDGES
     *  - \c EDGES_HAVE_ADIACENCY
     *
     *  @{
     */

    //! Get the number of switch which have been set.
    /**
     *  @return The number of switches which have been set.
     */
    const UInt& numSwitches() const
    {
        return switches.size();
    }

    //! Interrogate Switch.
    /**
     *  @param _s Name of the switch.
     *  @return true if the name is in the switch, false otherwise.
     */
    bool getLinkSwitch( std::string const & _s ) const;

    //! Set a switch using a given name
    /**
     *  @param _s Name of the switch.
     */
    void setLinkSwitch( std::string const & _s );

    //! Unset a switch.
    /**
     *  @param _s Name of the switch.
     */
    void unsetLinkSwitch( std::string const & _s );

    /** @} */ // End of group Switch Methods


    /** @name Generic Methods
     *  @ingroup public_methods
     *
     *  These are the generic methods to get information about the number of
     *  elements.
     *
     *  It is a generic interface common for all RegionMeshes (3D -- 1D).
     *
     *  @{
     */

    //! Number of elements in mesh.
    /**
     *  @return Number of elements in mesh.
     *  @sa numFaces
     *  @note Alias to numFaces()
     */
    UInt numElements() const;
    
    //! Number of elements in mesh.
    /**
     *  @return Number of elements in mesh.
     *  @sa numFaces
     *  @note Alias to numFaces()
     */
    UInt & numElements();
    
    //! Number of Boundary faces.
    /**
     *  @return Number of Boundary faces.
     */
    UInt numBElements() const;

    //! Number of boundary faces.
    /**
     *  @return Number of boundary faces.
     */
    UInt& numBElements();

    //! Get element at the i-th index.
    /**
     * @param i Index of the element
     * @return Element at index i
     */
    ElementType& element( const UInt& i );

    //! Get element at the i-th index.
    /**
     * @param i Index of the element
     * @return Element at index i
     */
    const ElementType& element( const UInt& i ) const;

    //! Get boundary element at the i-th index.
    /**
     * @param i Index of the element
     * @return Boundary element at index i
     */
    BElementType& bElement( const UInt& i );

    //! Get boundary element at the i-th index.
    /**
     * @param i Index of the element
     * @return Boundary element at index i
     */
    const BElementType& bElement( const UInt& i ) const;

    /** @} */ // End of group Generic Methods


    /** @name Volume Methods
     *  @ingroup public_methods
     *
     *  These are generic methods related to volumes.
     *
     *  @{
     */

    //! Returns Number of Volumes (zeros in 2D)
    UInt numVolumes()       const                { return 0; }
    //! Returns Global Number of Volumes (zeros in 2D)
    UInt numGlobalVolumes() const                { return 0; }

    /** @} */ // End of group Volume Methods


    /** @name Faces Methods
     *  @ingroup public_methods
     *
     *  All methods which operates on 2D elements.
     *
     *  @{
     */

    //! Returns Number of Faces
    /**
     *  Returns number of Face elements in the mesh
     *  as given by the internal counter.
     * 
     *  @return Number of Faces.
     */
    UInt numFaces() const;
    
    //! Returns Global Number of Faces
    /**
     *  Returns global number of Face elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Global Number of Faces.
     */
    UInt numGlobalFaces() const;

    //! Access to the Number of Faces
    /**
     *  Access number of Face (internal counter).
     *
     *  @return Access number of Face (internal counter).
     */
    UInt & numFaces();

    //! Number of Faces actually stored.
    /**
     *  Faces actually stored in list
     *  @return Number of faces actually stored in list.
     */
    UInt storedFaces() const;

    //! Current capacity of Faces Container.
    /**
     * How many elements may be stored.
     * @return Number of maximum element that can be stored.
     */
    UInt maxNumFaces() const;

    //! Changes Current capacity of Faces.
    /**
     *  Changes Current capacity of Faces (Optionally sets internal counter).
     * 
     *  @param n Maximum number of faces.
     *  @param setcounter true to set the counter, false otherwise.
     */
    void setMaxNumFaces( UInt const n, bool const setcounter = false );
    
    //! Changes Current capacity of Global Faces.
    /**
     *  Changes Current capacity of Global Faces (Optionally sets internal counter).
     *
     *  @param n maximum number of global faces.
     */
    void setMaxNumGlobalFaces( UInt const n );

    //! Adds faces.
    /**
     *  Adds faces. Id computed automatically.
     *  @return Reference to added face.
     */
    FaceType & addFace();

    //! Adds faces.
    /**
     *  Adds faces. Id computed automatically.
     *  @param v Face to add.
     *  @return Reference to the newly added face.
     */
    FaceType & addFace( FaceType const & v );

    //! Add face in a certain position.
    /**
     *  Add face to a specified position.
     *  @param v Face to add.
     *  @param pos Position of the face.
     *  @return Reference to the newly added face.
     */
    FaceType & setFace( FaceType const & v, UInt const pos );

    //! set numFaces counter.
    void setFaceCounter();

    //! Reference to last face stored in list.
    /**
     *  Reference to last face stored in list.
     *  Useful for mesh readers.
     *  @return reference of the last face in the list.
     */
    FaceType & lastFace();

    //! i-th mesh 2Delement.
    /**
     *  @param i index of the mesh 2Delement.
     *  @return the i-th face.
     */
    FaceType const & face( UInt const i ) const;

    //! i-th mesh 2Delement.
    /**
     *  @param i index of the mesh face.
     *  @return reference to the ith mesh face.
     */
    FaceType & face( UInt const i );

    /** @} */ // End of group Faces Methods


    /** @name Edge Methods
     *  @ingroup public_methods
     *  All methods which operates on 1D elements.
     *
     *  There are different point counters which may bi interrogated or set:
     *  - <strong>Number of Edges</strong>: is the declared number of total edges in the mesh.\n
     *   <strong>Beware</strong>: A value different from zero does NOT imply that the
     *   edges are actually stored.
     *  - <strong>Number of Boundary Edges</strong>: is the declared number of boundary edges
     *   in the mesh.\n
     *   <strong>Beware</strong>: A value different from zero does NOT imply that
     *   the boundary edges are actually stored.
     *  - <strong>Number of Stored Edges</strong>: It is the number of Edges actually stored
     *   on the edge container.
     *  - <strong>Maximum number of stored edges</strong>: The number of edges that may stored
     *   before the container is resized.\n
     *   <strong>Important</strong>: This parameter has to be
     *   set BEFORE inserting edges in the container if we want that pointer
     *   into the container maintains their validity. The container will also
     *   have a better performance.
     *
     *  @note To have more information on the Edge methods
     *  look at the documentation of th eanalogous Edges methods.
     *
     *  @{
     */

    //! Number of Edges.
    /**
     *  Returns number of Edge elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Number of Edges.
     */
    UInt numEdges() const;

    //! Global number of Edges.
    /**
     *  Returns global number of Edge elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Global number of Edges.
     */
    UInt numGlobalEdges() const;

    //! Number of local edges for each (2D) element.
    /**
     * @return Number of local edges for each (2D) element
     */
    UInt numLocalEdges() const;

    //! Access number of Edges.
    /**
     *  Access to the internal counter.
     *
     *  @return Access number of Edges.
     */
    UInt & numEdges();

    //! Number of stored Edges.
    /**
     *  Edges actually stored in list.
     *
     *  @return number of edges actually stored in list.
     */
    UInt storedEdges() const;

    //! Capacity of Edge Container.
    /**
     *  Current capacity of Edges Container,
     *  i.e. how many elements may be stored.
     *
     *  @return Edges Container Capacity.
     */
    UInt maxNumEdges() const;

    //! Changes Current capacity of Edges.
    /**
     *  Optionally sets internal counter.
     *
     *  @param n Maximum number of edges.
     *  @param setcounter true to set the counter, false otherwise.
     */
    void setMaxNumEdges( UInt const n, bool const setcounter = false );

    //! Changes Current capacity of Global Edges.
    /**
     *  Optionally sets internal counter.
     *
     *  @param n Maximum number of global edges.
     */
    void setMaxNumGlobalEdges( UInt const n );

    //! Adds Edges.
    /**
     *  Adds edges. Id computed automatically.
     *
     *  @param boundary true if is on boundary.
     *  @return Reference to added edge.
     */
    EdgeType & addEdge( bool const boundary = false );

    //! Adds Edges.
    /**
     *  Adds a edge (optionally a boundary edge) to the end of the list
     *  and adjourn its ID.
     *
     *  @param f Edge to add.
     *  @param boundary true if is on boundary.
     *  @return Reference to added edge.
     */
    EdgeType & addEdge( EdgeType const & f, bool const boundary = false );

    //! Adds Edges to specified position.
    /**
     *  Adds a edge (optionally a boundary edge) and adjourn its Id.
     *
     *  @param f Edge to add.
     *  @param position Position of the edge.
     *  @param boundary true if is on boundary.
     *  @return Reference to added edge.
     */
    EdgeType & setEdge( EdgeType const & f, UInt position, bool const boundary = false );

    //! Reference to last edge stored in list.
    /**
     *  Useful for mesh readers
     *  @return reference of the last edge in the list.
     */
    EdgeType & lastEdge();

    //! i-th mesh Edge.
    /**
     *  Returns the i-th Edge.
     *
     *  @param i Index of the mesh Edge.
     *  @return The i-th Edge.
     */
    EdgeType const & edge( UInt const i ) const;

    //! i-th mesh 1D Edge reference.
    /**
     *  Returns a reference to the i-th mesh Edge.
     *
     *  @param i Index of the mesh 1D Edge.
     *  @return Reference to the i-th Edge.
     */
    EdgeType & edge( UInt const i );

    //! i-th mesh 1D Boundary Edge.
    /**
     *  Returns the i-th mesh Boundary Edge.
     *
     *  @param i Index of the mesh 1D Boundary Edge.
     *  @return i-th Boundary Edge.
     */
    EdgeType const & boundaryEdge( UInt const i ) const;

    //! i-th mesh 1D Boundary Edge Reference.
    /**
     *  Returns a reference to the i-th mesh Boundary Edge.
     *
     *  @param i Index of the mesh 1D Boundary Edge.
     *  @return Reference to the i-th Boundary edge.
     */
    EdgeType & boundaryEdge( UInt const i );

    //! Set boundary Edge counter.
    /**
     *  Set the Boundary Edges counter to a given number.
     *
     *  @param n Count of Boundary Edge.
     */
    void setNumBEdges( UInt const n );

    //! Do I store mesh edges?
    /**
     *  Returns true if edges are stored.
     *
     *  @return true if edge are stored, false otherwise.
     */
    bool hasEdges() const;

    //! Do I store mesh internal edges?
    /**
     *  Returns true if internal edges are stored.
     *
     *  @return true if internal edge are stored, false otherwise.
     */
    bool hasInternalEdges() const;

    //! Number of Boundary Edges.
    /**
     *  Returns the number of boundary edges.
     *
     *  @return Number of boundary Edges.
     */
    UInt numBEdges() const;

    //! Edge on boundary check.
    /**
     *  Is this edge on boundary?
     *
     *  @param f The Edge.
     *  @return true if the edge is on the boundary, false otherwise.
     */
    bool isBoundaryEdge( EdgeType const & f ) const;

    //! Edge on boundary check by id.
    /**
     *  Is this edge, of given id, on boundary?
     *
     *  @param id Id of the edge.
     *  @return true if the edge is on the boundary, false otherwise.
     */
    bool isBoundaryEdge( UInt const & id ) const;

    //! Full Edge check by id.
    /**
     *  Does this ID corresponds to a full edge?
     *
     *  A FULL EDGE is a 1D Element that is actually stored in the Edge container.
     *
     *  @param id Id of the edge.
     *  @return true if the edge is actually stored, false otherwise.
     */
    bool isFullEdge( UInt const & id ) const;

    /** @} */ // End of group Edges Methods
    

    /** @name Points Methods
     *  @ingroup public_methods
     *  All methods which operates on Point Elements.
     *
     *  There are different Point counters which may be interrogated or set:
     *  - <strong>Number of Points</strong>: is the declared number of total points in
     *   the mesh.
     *  - <strong>Number of Boundary Points</strong>: is the declared number of boundary
     *   points in the mesh.\n
     *   Does NOT imply that the boundary points are actually stored.
     *  - <strong>Number of Stored Points</strong>: It is the number of Points actually stored
     *   on the points container.
     *  - <strong>Maximum number of stored points</strong>: The number of points that may
     *   stored before the container is resized.\n
     *   <strong>Very Important</strong>: This parameter has to be set <strong>BEFORE</strong>
     *   inserting points in the container. Since GeoElements will contain POINTERS
     *   into the container!\n
     *   This is a debatable point!
     *
     *   @note To have more information on the Point methods look at the
     *   documentation of the analogous Edges methods.
     *
     *  @{
     */

    //! Returns number of points in the mesh.
    /**
     *  Interrogation to the counter of points in the mesh.
     *
     *  @return Number of points in the mesh.
     */
    UInt numPoints() const;

    //! Returns a reference to number of points in the mesh.
    /**
     *  Returns a reference to the internal counter of the number of points in the mesh.
     *
     *  @return Reference to the internal counter of number of points in the mesh.
     */
    UInt & numPoints();

    //! Returns number of points in the mesh actually stored.
    /**
     *  Returns number of points in the mesh actually stored interrogating the internal counter.
     *
     *  @return Number of stored points in the mesh.
     */
    UInt storedPoints() const;

    //! Returns number of boundary points in the mesh actually stored.
    /**
     *  Returns number of boundary points in the mesh actually stored interrogating the internal counter.
     *
     *  @return Number of stored boundary points in the mesh.
     */
    UInt storedBPoints() const;

    //! Returns the number of storable points in the mesh.
    /**
     *  Returns the number of storable points in the mesh interrogating the internal counter.
     *
     *  @return The number of storable points.
     */
    UInt maxNumPoints() const;

    //! Set the number of storable points in the mesh.
    /**
     *  Set the internal counter of storable points in the mesh.
     *
     *  @param n Maximum number of storable points in the mesh.
     *  @param setcounter If true, it sets the internal counter, otherwise not (default).
     */
    void setMaxNumPoints( UInt const n, bool const setcounter = false );

    //! Set the number of storable global points in the mesh.
    /**
     *  Set the internal counter of storable global points in the mesh.
     *
     *  @param n Maximum number of storable global points in the mesh.
     */
    void setMaxNumGlobalPoints( UInt const n);
    
    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    PointType & addPoint( bool const boundary = false, bool const vertices = false );

    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param p Point to be added.
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    PointType & addPoint( PointType const & p, bool const boundary = false, bool const vertices = false );

    //! Adds a Point in the mesh giving an id.
    /**
     *  Adds a Point inside the mesh giving an id, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param p Point to be added.
     *  @param position Desired id.
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    PointType & setPoint( PointType const & p, UInt const position, bool const boundary = false, bool const vertices = false );

    //! Returns the last mesh Point.
    /**
     *  Returns the last Point in the mesh.
     *
     *  @return Reference to the last mesh Point.
     */
    PointType & lastPoint();
    
    //! Returns the i-th mesh Point.
    /**
     *  Returns the i-th Point in the mesh.
     *
     *  @param i Id of the Point.
     *  @return i-th mesh Point.
     */
    PointType const & point( UInt const i ) const;

    //! Returns a reference to the i-th mesh Point.
    /**
     *  Returns the i-th Point in the mesh.
     *
     *  @param i Id of the Point.
     *  @return Reference i-th mesh Point.
     */
    PointType & point( UInt const i );

    //! Returns a reference to the i-th mesh Boundary Point.
    /**
     *  Returns the i-th Boundary Point in the mesh.
     *
     *  @param i Id of the Boundary Point.
     *  @return Reference i-th mesh Boundary Point.
     */
    PointType const & boundaryPoint( UInt const i ) const;

    //! Returns a reference to the i-th mesh Boundary Point.
    /**
     *  Returns the i-th Boundary Point in the mesh.
     *
     *  @param i Id of the Boundary Point.
     *  @return Reference i-th mesh Boundary Point.
     */
    PointType & boundaryPoint( UInt const i );

    //! Returns the number of Boundary Points.
    /**
     *  Returns the counter of the number of Boundary Points.
     *
     *  @return Number Boundary Points.
     */
    UInt numBPoints() const;

    //! Sets the number of Boundary Points.
    /**
     *  Sets the counter of the number of Boundary Points.
     *
     *  @param n Number Boundary Points.
     */
    void setNumBPoints( UInt const n );
    
    /** @} */ // End of group Points Methods


    /** @name Vertices Methods
     *  @ingroup public_methods
     *  All methods which operates on Vertex Elements.
     *
     *  @{
     */

    //! Number of Vertices in Region.
    /**
     *  Returns the number of Vertices in Region.
     *
     *  @return Number of Vertices in Region.
     */
    UInt  numVertices () const;

    //! Reference to the counter of Vertices in Region.
    /**
     *  Allows to change number of Vertices in Region.
     *
     *  @return Reference to the counter of Vertices in Region.
     */
    UInt& numVertices ();

    //! Number of local vertices for each (2D) element.
    /**
     *  @return Number of local vertices for each (2D) element.
     */
    UInt numLocalVertices() const;

    //! Returns the global number of vertices in the mesh.
    /**
     *  Interrogation to the counter of vertices in the mesh.
     *
     *  @return Number of vertices in the mesh.
     */
    UInt numGlobalVertices() const;

    //! Number of Boundary Vertices in Region.
    /**
     *  Returns the number of Boundary Vertices in Region.
     *
     *  @return Number of Vertices in Region.
     */
    UInt numBVertices() const;

    //! Reference to the counter of Boundary Vertices in Region.
    /**
     *  Allows to change number of Boundary Vertices in Region.
     *
     *  @return Reference to the counter of Boundary Vertices in Region.
     */
    UInt& numBVertices();

    //! Vertex check.
    /**
     *  Is this Point a Vertex?
     *
     *  @param id Point's id.
     *  @return true if the Point is a Vertex, false otherwise.
     */
    bool isVertex ( UInt const & id ) const;

    //! Vertex check.
    /**
     *  Is this Point a Vertex?
     *
     *  @param p A Point.
     *  @return true if the Point is a Vertex, false otherwise.
     */
    bool isVertex ( PointType const & p ) const;

    //! Vertex on boundary check.
    /**
     *  Is this Point on boundary?
     *
     *  @param id Point's id.
     *  @return true if the Point is on boundary, false otherwise.
     */
    bool isBoundaryPoint ( UInt const & id )       const;

    //! Vertex on boundary check.
    /**
     *  Is this Point on boundary?
     *
     *  @param p A Point.
     *  @return true if the Point is on boundary, false otherwise.
     */
    bool isBoundaryPoint ( PointType const & p ) const;

    
    //! Changes number of Vertices.
    /**
     *  Allows to change number of vertices in Region.
     *
     *  @param n Number of vertices in the mesh.
     */
    void setNumVertices(UInt const n);
    
    //! Set the number of vertices in the mesh.
    /**
     *  Set the internal counter of vertices points in the mesh.
     *
     *  @param n Number of vertices in the mesh.
     */
    void setNumGlobalVertices( UInt const n ) { M_numGlobalVertices = n; }
    
    //! Changes number of Boundary Vertices.
    /**
     *  Allows to change number of boundary vertices in Region.
     *
     *  @param n Number of boundary vertices in the mesh.
     */
    void setNumBVertices(UInt const n);
    
    /** @} */ // End of group Vertices Methods


    /** @name Element Adjacency Methods
     *  @ingroup public_methods
     *  Methods to obtain the ID of Edges belonging to an element.
     *
     *  Accessing this information requires that the appropriate data
     *  structures have been set up by using the updateElementEdges() method.
     *
     *  Often, it is NOT required to have the full information about edges and faces:
     *  The ID of the Face and Edge entities may be calculated without
     *  contructing the corresponding Edge of Face Object. This saves memory.
     *
     *  @{
     */

    //! Is the array for local Edges set up?
    /**
     *  It does not use switches, but interrogates the container directly.
     *
     *  @return true if edges information are available, false otherwise.
     */
    bool hasLocalEdges() const;

    //! Edge Id of a certain edge number around a Face.
    /**
     *  @param iface Reference to Face.
     *  @param locE Position of the surrounding edges.
     *  @return Edge Id.
     */
    UInt localEdgeId( const FaceType& iface, const UInt locE ) const;

    //! Edge Id of a certain edge number around a Face.
    /**
     *  @param facId Face Id.
     *  @param locE Position of the surrounding edges.
     *  @return Edge Id.
     */
    UInt localEdgeId( const UInt facId, const UInt locE ) const;

    //! Element adjacent to a Edge.
    /**
     *  Returns ID of adjacent Element to a given Edge.
     *
     *  The first element is the one <strong>ORIENTED coherently with the edge</strong> (AS
     *  STORED in Edges). It means that the edge orientation is OUTWARD with
     *  respect to the element.\n
     *  The second element is either null (boundary edge) or indicates that
     *  the normal of the edge appears INWARD with respect to that element.
     *
     *  @param edgeId Id of a given edge.
     *  @param Pos Position equal to 1 or 2 indicates first or second adjacent Element.
     *  @return Edge ID given.
     */
    UInt edgeElement( UInt const edgeId, UInt const Pos ) const;

    //! Element adjacent to a Edge by reference.
    /**
     *  Returns ID of adjacent Element to a given Edge.
     *
     *  The first element is the one <strong>ORIENTED coherently with the edge</strong> (AS
     *  STORED in Edges). It means that the edge orientation is OUTWARD with
     *  respect to the element.\n
     *  The second element is either null (boundary edge) or indicates that
     *  the normal of the edge appears INWARD with respect to that element.
     *
     *  @param f Reference of a given edge.
     *  @param Pos Position equal to 1 or 2 indicates first or second adjacent Element.
     *  @return Edge ID given.
     */
    UInt edgeElement( EdgeType const & f, UInt const Pos ) const;

    //! Builds Edge-To-Face lookup table.
    /**
     *  Builds the lookup table for adjacency queries.
     *
     *  @param ce true if edge Counter is set, false otherwise (default).
     *  @param ee Edge count for memory reservation (0 default).
     */
    void updateElementEdges(bool ce=false, UInt ee=0 );

    //! Destroys Edge-To-Face lookup table.
    void cleanElementEdges();

    //! Returns Global-to-Local Node map.
    /**
     *  @return Global-to-Local map.
     */
    std::map<int,int> & globalToLocalNode() {return M_globalToLocalNode;}
    //! Returns Local-to-Global Node map.
    /**
     *  @return Local-to-Global map.
     */
    std::map<int,int> & localToGlobalNode() {return M_localToGlobalNode;}

    /** @} */ // End of group Element Adjacency Methods


    /** @defgroup public_attributes Public Attributes
     *
     */

    /** @name Region Containers
     *  @ingroup public_attributes
     *
     *  STL compliant containers for basic structure.
     *
     *  I expose them since they are standard containers.
     *
     *  @{
     */

    //! Container of mesh Points/Vertices.
    Points pointList;
    //! Container of mesh Faces.
    Faces faceList;
    //! Container of mesh Edges.
    Edges edgeList;
    //! Boundary points list.
    SimpleVect<PointType * > _bPoints;

    /** @} */ // End of group Region Containers


    /** @name Switches
     *  @ingroup public_attributes
     *
     *  @{
     */

    //! Switches
    Switch switches;

    /** @} */ // End of group Switches


protected:

    /** @defgroup protected_methods Protected Methods
     */

    //! Number of Elements in a list.
    /**
     *  Returns the number of elements in a given list.
     *
     *  @param list SimpleVect list of elements.
     *  @return Number of elements in list.
     */
    template < typename T >
    UInt numItems( SimpleVect< T> const & list ) const;

    //! Maximum number of Elements in a list.
    /**
     *  Returns maximum number of elements in a given list.
     *
     *  @param list SimpleVect list of elements.
     *  @return Maximum number of elements in list.
     */
    template < typename T >
    UInt maxNumItems( SimpleVect< T> const & list ) const;

    //! Set maximum number of Elements in a list.
    /**
     *  Set maximum number of elements in a given list.
     *
     *  @param list SimpleVect list of elements.
     *  @param n Number of elements.
     *  @param title Title for verbose output.
     */
    template < typename T >
    void setMaxNumItems( SimpleVect< T> & list, UInt n, std::string title );

    /** @defgroup protected_attributes Protected Attributes
     */

    /** @name Face-To-Edge and Boundary Edges Containers
     *  @ingroup protected_attributes
     *
     *  Arrays containing the ids of Edges and Edges of each element
     *  I use a Define to use localto global array or directly the
     *  bareedges.
     *
     *  - \c SAVEMEMORY defined: Try to save some memory.
     *  - \c NOT_BDATA_FIRST defined: boundary edge are stored apart.
     *
     *  @{
     */
#ifdef SAVEMEMORY
    //! Face-To-Edge Container.
    BareItemsHandler<BareEdge> _FToE;
#else
    //! Face-To-Edge Container.
    SimpleArray<UInt> _FToE;
#endif

#ifdef NOT_BDATA_FIRST
    //! Boundary Edges Container
    SimpleVect<EdgeType * > _bEdges;
#endif

    /** @} */ // End of group Face-To-Edge and Boundary Containers

    /** @name Internal Counters
     *  @ingroup protected_attributes
     *
     *  @{
     */

    //! Vertices Number
    UInt M_numVertices;
    //! Boundary Vertices Number
    UInt M_numBVertices;
    //! Points Number
    UInt M_numPoints;
    //! Boundary Points Number
    UInt M_numBPoints;
    //! Edges Number
    UInt M_numEdges;
    //! Boundary Edges Number
    UInt M_numBEdges;
    //! Faces Number
    UInt M_numFaces;
    //! Boundary Faces Number
    UInt M_numBFaces;

    //! Global Vertices Number
    UInt M_numGlobalVertices;
    //! Global Points Number
    UInt M_numGlobalPoints;
    //! Global Edges Number
    UInt M_numGlobalEdges;
    //! Global Faces Number
    UInt M_numGlobalFaces;

    //! Nodes Global to Local map
    std::map<int, int>      M_globalToLocalNode;
    //! Nodes Local to Global map
    std::map<int, int>      M_localToGlobalNode;
    //! Edges Global to Local map
    std::map<int, int>      M_globalToLocalEdge;
    //! Faces Global to Local map
    std::map<int, int>      M_globalToLocalFace;
    //! Volumes Global to Local map
    std::map<int, int>      M_globalToLocalVolume;

    /** @} */ // End of group Internal Counters

}; // End of class RegionMesh2D


// =================================================== //
// =================================================== //
//                    IMPLEMENTATION                   //
// =================================================== //
// =================================================== //


template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::RegionMesh2D() :
        MeshEntity(),
        MC::RegionMarker(),
        switches(),
        M_numVertices( 0 ),
        M_numBVertices( 0 ),
        M_numPoints( 0 ),
        M_numBPoints( 0 ),
        M_numFaces( 0 ),
        M_numEdges( 0 ),
        M_numBEdges( 0 ),
        M_globalToLocalNode(),
        M_localToGlobalNode(),
        M_globalToLocalEdge(),
        M_globalToLocalFace(),
        M_globalToLocalVolume()
{ }


template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::RegionMesh2D( UInt id ) :
        MeshEntity( id ),
        MC::RegionMarker(),
        switches(),
        M_numVertices( 0 ),
        M_numBVertices( 0 ),
        M_numPoints( 0 ),
        M_numBPoints( 0 ),
        M_numFaces( 0 ),
        M_numEdges( 0 ),
        M_numBEdges( 0 )
{ }

template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::~RegionMesh2D() {}

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
    std::ostringstream _err_msg;
    _err_msg << "Switch named " << _s << " is not allowed";
    ASSERT0( switches.set( _s ), _err_msg.str().c_str() );
};

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::getLinkSwitch( std::string const & _s ) const
{
    return switches.test( _s );
};

template <typename GEOSHAPE, typename MC>
INLINE
void
RegionMesh2D<GEOSHAPE, MC>::unsetLinkSwitch( std::string const & _s )
{
    std::ostringstream _err_msg;
    _err_msg << "Switch named " << _s << " is not allowed";
    ASSERT0( switches.unset( _s ), _err_msg.str().c_str() );
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
    return M_numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt & RegionMesh2D<GEOSHAPE, MC>::numElements()
{
    return M_numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt RegionMesh2D<GEOSHAPE, MC>::numBElements() const
{
    return M_numBEdges;
}

template <typename GEOSHAPE, typename MC>
UInt &RegionMesh2D<GEOSHAPE, MC>::numBElements()
{
    return M_numBEdges;
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh2D<GEOSHAPE, MC>::ElementType &
RegionMesh2D<GEOSHAPE, MC>::element( UInt const & i )
{
    return face( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh2D<GEOSHAPE, MC>::ElementType const &
RegionMesh2D<GEOSHAPE, MC>::element( UInt const & i ) const
{
    return face( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh2D<GEOSHAPE, MC>::BElementType &
RegionMesh2D<GEOSHAPE, MC>::bElement( UInt const & i )
{
    return boundaryEdge( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh2D<GEOSHAPE, MC>::BElementType const &
RegionMesh2D<GEOSHAPE, MC>::bElement( UInt const & i ) const
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
RegionMesh2D<GEOSHAPE, MC>::setMaxNumItems( SimpleVect< T> & list, UInt n, std::string title )
{
    if ( list.capacity() == 0 )
    {
        list.reserve( n );
    }
    else if ( list.capacity() < n )
    {
        Debug(4000) << "Resetting " << title << " list size to " << n << "\n";
        Debug(4000) << "ALL PREVIOUS POINTERS TO THE LIST (IF ANY) ARE NOW INVALID\n";

        list.reserve( n );
    }
}
//-------------------------------------------------------------------

// ***************************** FACES
template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numFaces() const
{
    return M_numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numGlobalFaces() const
{
    return M_numGlobalFaces;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numFaces()
{
    return M_numFaces;
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
        M_numFaces = n;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumGlobalFaces( UInt const n)
{
    M_numGlobalFaces = n;
}

// \todo use addItem
template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::addFace()
{
    return addFace( FaceType() );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::addFace( FaceType const & v )
{
    ASSERT_PRE( faceList.size() < faceList.capacity() , "Face list size exceeded" <<
                faceList.size() + 1 << " " << faceList.capacity() ) ;
    faceList.push_back( v );
    ( faceList.back() ).setId(faceList.size());
    return faceList.back();
}
// \todo Use setItem

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::setFace( FaceType const & v, UInt const pos )
{
    ASSERT_PRE( pos <= faceList.capacity() , "position requested exceed capacity" <<
                pos << " " << faceList.capacity() ) ;
    faceList( pos ) = v;
    faceList( pos ).setId(pos);
    return faceList( pos );
}

template <typename GEOSHAPE, typename MC>
INLINE
void
RegionMesh2D<GEOSHAPE, MC>::setFaceCounter()
{
    M_numFaces = faceList.size();
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::lastFace()
{
    return faceList.back();
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::FaceType const &
RegionMesh2D<GEOSHAPE, MC>::face( UInt const i ) const
{
    ASSERT_BD( i > 0 && i <= faceList.size() ) ;
    return faceList( i );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::FaceType &
RegionMesh2D<GEOSHAPE, MC>::face( UInt const i )
{
    ASSERT_BD( i > 0 && i <= faceList.size() ) ;
    return faceList( i );
}

// ************************* EDGES ******************************
template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numEdges() const
{
    return M_numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numGlobalEdges() const
{
    return M_numGlobalEdges;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numEdges()
{
    return M_numEdges;
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
        M_numEdges = n;
}


template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumGlobalEdges( UInt const n)
{
    M_numGlobalEdges = n;
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::addEdge( bool const boundary )
{
    return addEdge( EdgeType(), boundary );
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::addEdge( EdgeType const & f, bool const boundary )
{
    ASSERT_PRE( edgeList.size() < edgeList.capacity(), "Edge list size exceeded" <<
                edgeList.size() + 1 << " " << edgeList.capacity() ) ;
    edgeList.push_back( f );
    ( edgeList.back() ).setId( edgeList.size() );
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
typename RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::setEdge( EdgeType const & f, UInt position, bool const boundary )
{
    ASSERT_PRE( position <= edgeList.capacity(), "Edge list size exceeded" <<
                position << " " << edgeList.capacity() ) ;
    edgeList( position ) = f;
    edgeList( position ).setId( position );
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
typename RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::lastEdge()
{
    return edgeList.back();
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::EdgeType const &
RegionMesh2D<GEOSHAPE, MC>::edge( UInt const i ) const
{
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::edge( UInt const i )
{
    if ( i <= 0 || i > edgeList.size() ) std::cout<< "i: " << i << "edgeList.size(): " <<edgeList.size()<<std::endl;
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::EdgeType const &
RegionMesh2D<GEOSHAPE, MC>::boundaryEdge( UInt const i ) const
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
typename RegionMesh2D<GEOSHAPE, MC>::EdgeType &
RegionMesh2D<GEOSHAPE, MC>::boundaryEdge( UInt const i )
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
    M_numBEdges = n;
#ifdef NOT_BDATA_FIRST

    M_bEdges.reserve( n );
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
    return edgeList.size() > M_numBEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numBEdges() const
{
    return M_numBEdges;
}




template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryEdge( EdgeType const & e ) const
{
#ifdef NOT_BDATA_FIRST
    //ASSERT(false,"In this version Boundary edges must be stored first");
    bool isboundary = true;
    for ( UInt k = 1; k <= EdgeType::numVertices; ++k )
    {
        isboundary = isboundary & e.point( k ).boundary();
    }
    return isboundary;
#else

    return e.id() <= M_numBEdges;
#endif
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryEdge( UInt const & id ) const
{
    return isBoundaryEdge( edge( id ) );
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isFullEdge( UInt const & id ) const
{
    return edgeList.size() >= id;
}

template <typename GEOSHAPE, typename MC>
INLINE
UInt
RegionMesh2D<GEOSHAPE, MC>::edgeElement( UInt const i, UInt const Pos ) const
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
    return M_numPoints;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numPoints()
{
    return M_numPoints;
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
        M_numPoints = n;
}


template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumGlobalPoints( UInt const n )
{
    M_numGlobalPoints = n;
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::addPoint( bool const boundary, bool const vertex )
{
    return addPoint( PointType(), boundary, vertex );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::addPoint( PointType const & p, bool const boundary, bool const /*vertex*/ )
{
    ASSERT_PRE( pointList.size() < pointList.capacity(), "Point list size exceeded" <<
                pointList.size() + 1 << " " << pointList.capacity() ) ;
    pointList.push_back( p );
    PointType * pp = & pointList.back();
    pp->setId( pointList.size() );
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
typename RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::setPoint
( PointType const & p, UInt position, bool const boundary, bool const vertex )
{
    ASSERT_PRE( position <= pointList.capacity(), "Position  exceed lpoint list capacity" <<
                position << " " << pointList.capacity() ) ;
    bool found( false );
    pointList( position ) = p;
    PointType * pp = & pointList( position );
    pp->setId( position );
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
typename RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::lastPoint()
{
    return pointList.back();
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::PointType const &
RegionMesh2D<GEOSHAPE, MC>::point( UInt const i ) const
{
    ASSERT_BD( i > 0 && i <= pointList.size() ) ;
    return pointList( i );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::point( UInt const i )
{
    ASSERT_BD( i > 0 && i <= pointList.size() ) ;
    return pointList( i );
}


template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::PointType const &
RegionMesh2D<GEOSHAPE, MC>::boundaryPoint( UInt const i ) const
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
INLINE
typename RegionMesh2D<GEOSHAPE, MC>::PointType &
RegionMesh2D<GEOSHAPE, MC>::boundaryPoint( UInt const i )
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numBPoints() const
{
    return M_numBPoints;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setNumBPoints( UInt const n )
{
    M_numBPoints = n;
    _bPoints.reserve( n );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numVertices() const
{
    return M_numVertices;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numVertices()
{
    return M_numVertices;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numGlobalVertices() const
{
    return M_numGlobalVertices;
}

template <typename GEOSHAPE, typename MC>
void RegionMesh2D<GEOSHAPE, MC>::setNumVertices(UInt const n)
{
    M_numVertices=n;
}


template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numBVertices() const
{
    return M_numBVertices;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numBVertices()
{
    return M_numBVertices;
}

template <typename GEOSHAPE, typename MC>
void RegionMesh2D<GEOSHAPE, MC>::setNumBVertices(UInt const n)
{
    M_numBVertices=n;
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isVertex( PointType const & p ) const
{
    return p.id() <= M_numVertices;
}

template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isVertex( UInt const & id ) const
{
    return id <= M_numVertices;
}


template <typename GEOSHAPE, typename MC>
INLINE
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryPoint( UInt const & id ) const
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
UInt
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( const FaceType & ifac, UInt const locE ) const
{
    ASSERT_PRE( !_FToE.empty(), "Face to Edges array not  set" );
    ASSERT_BD( locE > 0 && locE <= FaceType::numLocalEdges );
    std::pair<BareEdge, bool> it;
    UInt i1, i2;
    i1 = GEOSHAPE::eToP( locE, 1 );
    i2 = GEOSHAPE::eToP( locE, 2 );
    i1 = ( ifac.point( i1 ) ).id();
    i2 = ( ifac.point( i2 ) ).id();
    it = makeBareEdge( i1, i2 );
    return _FToE.id( it.first );
}

template <typename GEOSHAPE, typename MC>
INLINE
UInt
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( UInt const facId, UInt const locE ) const
{
    ASSERT_BD( facId > 0 && facId <= M_numFaces );
    return localEdgeId( face( facId ), locE );
}

#else

template <typename GEOSHAPE, typename MC>
INLINE
UInt
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( const FaceType & ifac, UInt const locE )
const
{
    return _FToE( ifac.id(), locE );
}

template <typename GEOSHAPE, typename MC>
INLINE
UInt
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( UInt const facId, UInt const locE )
const
{
    ASSERT_PRE( !_FToE.empty(), "Face to Edges array not  set" );
    ASSERT_BD( facId > 0 && facId <= M_numFaces );
    ASSERT_BD( locE > 0 && locE <= FaceType::numLocalEdges );
    return _FToE( locE, facId );
}

#endif

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::updateElementEdges( bool ce, UInt ee )
{
//  std::cout<< "locEdge: " << localEdgeId( 1, 3 )<< std::endl;
    std::cout << "     Updating element edges ... " << std::flush;

    ASSERT0( ! ce || M_numBEdges > 0, std::stringstream( std::string("Boundary Edges Must have been set") +
                                                         std::string("in order to call updateElementEdges with createEdges=true") +
                                                         std::string("\nUse buildBoundaryEdges(..) from mesh_util.h") ).str().c_str() );
    // If the counter is set we trust it! Otherwise we use Euler formula

    if ( ce && ee == 0 )
        //check!!!
        ee = M_numEdges > M_numBEdges ? M_numEdges : M_numFaces + M_numVertices -1;
    ASSERT( ce || numEdges() > 0 , "Mesh is not properly set!" );

    if ( ce )
        edgeList.reserve( ee );

    EdgeType edg;


    BareItemsHandler<BareEdge> _be;
    std::pair<UInt, bool> e;
    _FToE.reshape( numLocalEdges(), numFaces() ); // DIMENSION ARRAY

    UInt fid, i1, i2;
    std::pair<BareEdge, bool>_edge;
    GEOSHAPE ele;
    // If we have all edges and the edges store all adjacency info
    // everything is easier
    if ( (edgeList.size() == numEdges()) & getLinkSwitch( "EDGES_HAVE_ADIACENCY" ) & getLinkSwitch( "HAS_ALL_EDGES" ) )
    {
        for ( typename Edges::iterator ite = edgeList.begin(); ite != edgeList.end(); ++ite )
        {
            if ( ite->pos_first() != 0 )
                _FToE( ite->pos_first() , ite->ad_first() ) = ite->localId();
            if ( ite->pos_second() != 0 )
                _FToE( ite->pos_second(), ite->ad_second() ) = ite->localId();
        }
        // we finish here
        setLinkSwitch( "HAS_FACE_TO_EDGES" );
        std::cout << " done." << std::endl;
        return ;
    }

    // If I have only boundary faces I need to process them first to keep the correct numbering

    // First We check if we have already Faces stored
    if ( ! edgeList.empty() )
    {
        // dump all faces in the container, to maintain the correct numbering
        // if everything is correct the numbering in the bareface structure
        // will reflect the actual face numbering However, if I want to create
        // the internal faces I need to make sure that I am processing only the
        // boundary ones. So I resize the container!
        if ( ce )
            edgeList.resize( M_numBEdges );

        std::pair<UInt, bool> _check;
        for ( UInt j = 0; j < edgeList.size(); ++j )
        {
            i1 = ( edgeList[ j ].point( 1 ) ).localId();
            i2 = ( edgeList[ j ].point( 2 ) ).localId();
            _edge = makeBareEdge( i1, i2);
            _check = _be.addIfNotThere( _edge.first );
        }
    }

    for ( typename Faces::iterator iface = faceList.begin();
            iface != faceList.end(); ++iface )
    {
        fid = iface->localId();
        for ( UInt j = 1; j <= numLocalEdges(); j++ )
        {
            i1 = ele.eToP( j, 1 );
            i2 = ele.eToP( j, 2 );
            i1 = ( iface->point( i1 ) ).localId();
            i2 = ( iface->point( i2 ) ).localId();
            _edge = makeBareEdge( i1, i2);
            e = _be.addIfNotThere( _edge.first );
            _FToE( j, fid ) = e.first;
            if ( ce )
            {
                if ( e.second )
                {
                    // a new face It must be internal.
                    for ( UInt k = 1; k <= EdgeType::numPoints; ++k )//
                        //iface->point( ele.eToP( j, k ) );
                        edg.setPoint( k, iface->point( ele.eToP( j, k ) ) );
                    edg.ad_first()  = fid;
                    edg.pos_first() = j;

                    // gets the marker from the RegionMesh

                    edg.setMarker( this->marker() );
                    //        inheritWeakerMarker( edg );
                    addEdge( edg, false ); //The id should be correct
                }
                else
                {
                    // We assume that BEdges have been already set so we have to do
                    // nothing if the edge is on the boundary
                    if ( e.first > M_numBEdges )
                    {
                        edgeList( e.first ).ad_second() = fid;
                        edgeList( e.first ).pos_second() = j;

                    }
                }
            }
        }
    }

    UInt n = _be.maxId();
    // Fix _numfaces if it was not set or set to just the # of Bfaces
    if (!ce)
    {
        if ( M_numEdges == 0 || M_numEdges == M_numBEdges )
            M_numEdges = n;
    }
    else
    {
        M_numGlobalEdges = n;
    }

    std::cout << n << " edges ";
    ASSERT_POS( n == M_numEdges , "#Edges found inconsistent with that stored in RegionMesh" ) ;
    setLinkSwitch( "HAS_FACE_TO_EDGES" );
    if ( ce )
        setLinkSwitch( "HAS_ALL_EDGES" );
    if ( ce )
        setLinkSwitch( "EDGES_HAVE_ADIACENCY" );

    std::cout << " done." << std::endl;
}




template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::cleanElementEdges()
{
    _FToE.clear();
    unsetLinkSwitch( "HAS_FACE_TO_EDGES" );
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
    out << " ID: " << this->id() << std::endl;
    out << "Edges local to  faces stored: " << hasLocalEdges() << std::endl;
    //out <<"Edges local to  faces   stored:"<<switches.test("FACEtoEDGE")<<std::endl;
    //out <<"Faces adjacent to Edges stored: "<<switches.test("EDGEtoFACE")<<std::endl<<std::endl;
    out << "Edges Stored: " << hasEdges() << " Internal: "
    << hasInternalEdges() << std::endl;
    out << "**************************************************" << std::endl;
    out << "numPoints=" << numPoints() << "  " << "numBPoints=" << numBPoints() << std::endl;
    out << "numVertices=" << numVertices() << "  " << "numBVerices=" << numBVertices() << std::endl;
    out << "numFaces=" << numFaces() << std::endl;
    out << "numEdges=" << numEdges() << "  " << "numBEdges=" << numBEdges() << std::endl;
    out << "**************************************************" << std::endl;
    switches.showMe( verbose, out );
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    if ( verbose )
    {
        std::cout << "Verbose version not implemented yet" << std::endl;
    }
    return out;

}
template <typename GEOSHAPE, typename MC>
int
RegionMesh2D<GEOSHAPE, MC>::check( int /*level*/, bool const fix, bool const verb, std::ostream & out )
{
    int severity = 0;
    if ( verb )
    {
        out << "**************************************************" << std::endl;
        out << "         Checkin  RegionMesh2D                " << std::endl;
        out << " ID: " << this->id() << std::endl;
        out << "**************************************************" << std::endl;
    }
    if ( pointList.size() != M_numPoints )
    {
        out << " Point list size " << pointList.size() << " not equal to internal counter value "
        << M_numPoints << std::endl;
        if ( fix )
        {
            M_numPoints = pointList.size();
            out << "Fixed";
            out.flush();
        }
    }
    if ( edgeList.size() == 0 )
    {
        if ( verb )
            out << "Warning: No Edges Stored" << std::endl;
        severity = -1;
    }
    if ( faceList.size() == 0 )
    {
        if ( verb )
            out << "Warning: No Faces Stored" << std::endl;
        severity = 1;
    }
    UInt count = 0;
    for ( typename Points::iterator i = pointList.begin(); i != pointList.end(); ++i )

        if ( i->boundary() )
            ++count;
    if ( !count ) severity = 4;
    if ( count != M_numBPoints )
    {
        out << " Num Boundary points " << count << " not equal to internal counter value "
        << M_numBPoints << std::endl;
        if ( fix )
        {
            M_numBPoints = count;
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
    {
        out << " SEVERITY ERROR:" << badid << "Points ids are wrong";
        severity = 5;
    }

    badid = 0;
    for ( UInt i = 1; i <= storedEdges(); ++i )
        if ( edge( i ).id() != i )
            ++badid;
    if ( badid != 0 )
    {
        out << " SEVERITY ERROR:" << badid << "Edges ids are wrong";
        severity = 5;
    }

    badid = 0;
    for ( UInt i = 1; i <= storedFaces(); ++i )
        if ( face( i ).id() != i )
            ++badid;
    if ( badid != 0 )
    {
        out << " SEVERITY ERROR:" << badid << "Faces ids are wrong";
        severity = 5;
    }

    badid = 0;

    if ( M_numVertices == 0 )
    {
        out << " SEVERITY ERROR: internal Vertices Counter unset";
        severity = 6;
    }
    if ( M_numPoints == 0 )
    {
        out << " SEVERITY ERROR: internal Points Counter unset";
        severity = 6;
    }
    if ( M_numPoints == 0 )
    {
        out << " SEVERITY ERROR: internal Points Counter unset";
        severity = 6;
    }
    if ( M_numBPoints == 0 )
    {
        out << " SEVERITY ERROR: boundary Points Counter unset";
        severity = 6;
    }
    if ( M_numBVertices == 0 )
    {
        out << " SEVERITY ERROR: boundary Vertices Counter unset";
        severity = 6;
    }

    if ( verb )
        out << "   Check Finished              " << std::endl <<
        "***********************************************" << std::endl;

    return severity;
}

} // End of namespace LifeV

#endif //REGIONMESH2D_H
