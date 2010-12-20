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
 *  @brief File containing 1D Mesh Classes
 *
 *  @date 27-04-2010
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *  @author Tiziano Passerini <tiziano.passerini@gmail.com>
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 *  @mantainer Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 */

#ifndef REGIONMESH1D_H
#define REGIONMESH1D_H

#include <life/lifecore/life.hpp>
#include <life/lifecore/LifeDebug.hpp>
#include <life/lifecore/Switch.hpp>

#include <life/lifemesh/geoElement.hpp>
#include <life/lifemesh/bareItems.hpp>
#include <life/lifemesh/basisElSh.hpp>
#include <life/lifearray/SimpleVect.hpp>

#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace LifeV
{
/** @class RegionMesh1D
 *  @brief Class for 1D Mesh
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *  @author Tiziano Passerini <tiziano.passerini@gmail.com>
 *
 *  This is the class that stores the mesh entities for a single 1D region.
 *
 *  In a region elements are all of the same type.
 */
template <typename GEOSHAPE, typename MC = DefMarkerCommon >
class RegionMesh1D : public MeshEntity,
        public MC::RegionMarker
{
public:
    /** @defgroup public_types Public Types
     *
     */

    /** @name Marker Types
     *  @ingroup public_types
     *  Markers for Point, Edge and Region
     *  @{
     */

    //! Common Marker Class
    typedef MC MarkerCommon;
    //! Point Marker
    typedef typename MC::PointMarker     PointMarker;
    //! Edge Marker
    typedef typename MC::EdgeMarker      EdgeMarker;
    //! Region Marker
    typedef typename MC::RegionMarker    RegionMarker;

    /** @} */ // End of group Marker Types

    /** @name Basic Element Shape Types
     *  @ingroup public_types
     *  Volume, Face and Edge geometric shapes.
     *  @{
     */

    //! Volume Shape.
    typedef GEOSHAPE VolumeShape;
    //! Face Shape (Boundary Element).
    typedef GEOSHAPE EdgeShape;
    //! Edge Shape (Boundary of Boundary Element)
    typedef typename GEOSHAPE::GeoBShape PointShape;

    /** @} */ // End of group Basic Element Shape Types

    /** @name Geometric Element Types
     *  @ingroup public_types
     *  Volumes, Faces, Edges and Points.
     *  @{
     */

    //! Volume Element (1D)
    typedef GeoElement1D<EdgeShape, MC>  VolumeType;
    //! Face Element (1D)
    typedef GeoElement1D<EdgeShape, MC>  FaceType;
    //! Edge Element (1D)
    typedef GeoElement1D<EdgeShape, MC>  EdgeType;
    //! Point Element (0D)
    typedef GeoElement0D<MC>             point_Type;

    /** @} */ // End of group Geometric Element Types

    /** @name Geometric Element Container Types
     *  @ingroup public_types
     *  Typedefs for STL compliant containers of mesh geometric entities.
     *
     *  I Use SimpleVect container for addressing from 1.
     *  @{
     */

    //! Points Container
    typedef SimpleVect<point_Type>  Points;
    //! Elements Container - compatibility
    typedef SimpleVect<VolumeType> Volumes;
    //! Faces Container - compatibility
    typedef SimpleVect<FaceType>   Faces;
    //! Edges Container: at least boundary edges
    typedef SimpleVect<EdgeType>   Edges;

    /** @} */ // End of group Geometric Element Container Types

    /** @name Generic Types
     *  @ingroup public_types
     *  Elements and Boundary Elements Types.
     *  @{
     */

    //! Element Geometric Shape
    typedef GEOSHAPE                     ElementShape;
    //! Boundary Element Geometric Shape
    typedef typename GEOSHAPE::GeoBShape BElementShape;

    //! Element Geometric Type
    typedef GeoElement1D<GEOSHAPE, MC>   ElementType;
    //! Boundary Element Geometric Type
    typedef GeoElement0D<MC>             BElementType;

    //! Element Geometric Shape Container Type
    typedef SimpleVect<EdgeType>         Elements;
    //! Boundary Element Geometric Shape Container Type
    typedef SimpleVect<point_Type>        BElements;

    /** @} */ // End of group Generic Types

    /** @name Constructors & Destructor
     *  Default and Copy Constructor for the class.
     *  @{
     */

    //! Default constructor
    /**
     *  @param id marker of the RegionMesh1D
     */
    explicit RegionMesh1D( UInt id = 0 );

    //! Copy constructor
    /**
     *  @param mesh a RegionMesh1D
     */
    explicit RegionMesh1D( const RegionMesh1D<GEOSHAPE, MC>& mesh );

    /** @} */ // End of group Constructors & Destructor

    /** @name Operators
     *  Public Operator Methods
     *
     *  @{
     */

    //! Assignment operator
    /**
     *  @param mesh a RegionMesh1D
     *  @return the newly copied RegionMesh1D
     */
    RegionMesh1D<GEOSHAPE, MC> operator=( RegionMesh1D<GEOSHAPE, MC> const & mesh );

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


    /** @name Geometric Methods
     *  @ingroup public_methods
     *  Geometric operations on mesh.
     *
     *  @{
     */

    //! Setup mesh
    /**
     *  @param Length Total length of the mesh.
     *  @param NumberOfElements Number of elements inside the mesh.
     *
     *  Mesh construction by hand.
     */
    void setup( const Real& Length, const UInt& NumberOfElements );

    //! Transform the mesh using boost::numeric::ublas.
    /** Scale, rotate and translate the mesh (operations performed in this order).
     *  @date   27/04/2010
     *  @author Cristiano Malossi
     *  @note - Rotation follows Paraview conventions: first rotate around z-axis,
     *        then around y-axis and finally around x-axis;
     *  @note - All the vectors must allow the command: operator[];
     *  @param scale        vector of three components for (x,y,z) scaling of the mesh
     *  @param rotate       vector of three components (radiants) for rotating the mesh
     *  @param translate    vector of three components for (x,y,z) translation the mesh
     *
     */
    template <typename VECTOR>
    void transformMesh( const VECTOR& scale, const VECTOR& rotate, const VECTOR& translate );

    //! Get the maximum H over all the edges of the mesh.
    /**
     *  @date 27/04/2010
     *  @author Cristiano Malossi
     *  @return Maximum element length H
     */
    Real maxH() const;

    //! Get the minumum H over all the edges of the mesh.
    /**
     *  @date 27/04/2010
     *  @author Cristiano Malossi
     *  @return Minimum element length H
     */
    Real minH() const;

    //! Get the mean H over all the edges of the mesh.
    /**
     *  @date 27/04/2010
     *  @author Cristiano Malossi
     *  @return Average element length H
     */
    Real meanH() const;

    /** @} */ // End of group Geometric Methods


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

    //! Set a given switch as switch inside the class.
    /**
     *  @param sw Switch to be set.
     */
    void setSwitch( Switch & sw );

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

    //! Number of local vertices for each (1D) element.
    /**
     *  @return Number of local vertices for each (1D) element.
     */
    UInt numLocalVertices() const;

    //! Number of elements in mesh.
    /**
     *  @return Number of elements in mesh.
     *  @sa numEdges
     *  @note Alias to numEdges()
     */
    UInt numElements() const;

    //! Number of 3D elements.
    /**
     *  @return Number of 3D elements.
     */
    UInt & numElements();

    //! Number of Global elements in mesh.
    /**
     *  @return Number of elements in mesh.
     *  @sa numEdges
     *  @note Alias to numEdges()
     */
    UInt numGlobalElements() const;

    //! Number of Global 3D elements.
    /**
     *  @return Number of 3D elements.
     */
    UInt & numGlobalElements();

    //! Number of Boundary edges.
    /**
     * @return Number of Boundary edges.
     */
    UInt numBElements() const;

    //! Number of boundary edges
    /**
     *  @return Number of boundary edges
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
     *  @note For compatibility reasons this methods do the same as edge methods.
     *
     *  @{
     */

    //! Returns Number of Volumes (not used in 1D)
    UInt numVolumes()       const                { return 0; }
    //! Returns Global Number of Volumes (not used in 1D)
    UInt numGlobalVolumes() const                { return 0; }
    //! Returns Volume for a given Id
    /**
     *  For compatibility reasons, it returns an edge (i.e. Volume for 1D meshes)
     *  @param i Volume id
     *  @return Reference to the Volume
     */
    const VolumeType& volume( const UInt i ) const { ASSERT_BD( i > 0 && i <= edgeList.size() ); return edgeList( i ); }
    //! Returns Volume for a given Id
    /** @copydoc volume */
    VolumeType& volume( const UInt i )             { ASSERT_BD( i > 0 && i <= edgeList.size() ); return edgeList( i ); }

    /** @} */ // End of group Volume Methods


    /** @name Faces Methods
     *  @ingroup public_methods
     *
     *  These are generic methods related to faces.
     *
     *  @note These methods are not used in 1D.
     *
     *  @{
     */

    //! Returns Number of Faces (not used in 1D)
    UInt numFaces()         const { return 0; }
    //! Returns Global Number of Faces (not used in 1D)
    UInt numGlobalFaces()   const { return 0; }
    //! Returns true if there is local Faces table (not used in 1D)
    bool hasLocalFaces() const { return true; }
    //! Returns local Face Id for a given Volume (not used in 1D)
    UInt localFaceId( const UInt /*volId*/, const UInt /*locF*/ ) const { return 0; }
    //! Returns local Face Id for a given Volume (not used in 1D)
    UInt localFaceId( const VolumeType& /*iv*/, const UInt /*locF*/ ) const {return 0;}
    //! Update Element-To-Faces lookup table (not implemented)
    void updateElementFaces( bool createFaces = false, UInt estimateFaceNumber = 0 );
    //! Clean Element-To-Faces lookup table (not implemented)
    void cleanElementFaces() {}

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

    //! Adds Edges.
    /**
     *  Adds edges. Id computed automatically.
     *
     *  @return Reference to added edge.
     */
    EdgeType & addEdge();

    //! Adds Edges.
    /**
     *  Adds edges. Id computed automatically.
     *
     *  @param v Edge to add.
     *  @return Reference to added edge.
     */
    EdgeType & addEdge( EdgeType const & v );

    //! Adds Edges to specified position.
    /**
     *  Adds edges. Id computed automatically.
     *
     *  @param v Edge to add.
     *  @param pos Position of the edge.
     *  @return Reference to added edge.
     */
    EdgeType & setEdge( EdgeType const & v, UInt const pos );

    //! Set numEdges counter.
    /**
     *  Set internal numEdges counter.
     */
    void setEdgeCounter();

    //! Reference to last edge stored in list.
    /**
     *  Useful for mesh readers
     *  @return reference of the last edge in the list.
     */
    EdgeType & lastEdge();

    //! i-th mesh 1D Element.
    /**
     *  Returns the i-th mesh Element.
     *
     *  @param i Index of the mesh 1D Element.
     *  @return The i-th edge.
     */
    EdgeType const & edge( UInt const i ) const;

    //! i-th mesh 1D Element reference.
    /**
     *  Returns a reference to the i-th mesh Element.
     *
     *  @param i Index of the mesh 1D Element.
     *  @return Reference to the i-th edge.
     */
    EdgeType & edge( UInt const i );

    //! i-th mesh 1D Element length.
    /**
     *  Returns the i-th mesh 1D Element length.
     *
     *  @param i Index of the mesh 1D Element.
     *  @return Length of the i-th Edge.
     */
    Real edgeLength( const UInt& i ) const;

    //! Adds a edge to list.
    /**
     *  @param boundary true if it's a boundary edge.
     *  @return Reference to the newly added edge.
     */
    EdgeType & addEdge( bool const boundary = false );

    //! Adds a edge to the end of the list and adjourn its ID.
    /**
     *  @param f Reference of edge to be added.
     *  @param boundary true if it's a boundary edge.
     *  @return Reference to the newly added edge.
     */
    EdgeType & addEdge( EdgeType const & f, bool const boundary = false );

    //! Adds a edge to the end of the list and adjourn its ID.
    /**
     *  @param f Reference of edge to be added.
     *  @param position Desired id for the edge.
     *  @param boundary true if it's a boundary edge.
     *  @return Reference to the newly added edge.
     */
    EdgeType & setEdge( EdgeType const & f, UInt position, bool const boundary = false );

    //! i-th mesh 1D Boundary Element.
    /**
     *  Returns the i-th mesh Boundary Element.
     *
     *  @param i Index of the mesh 1D Boundary Element.
     *  @return i-th boundary edge.
     */
    EdgeType const & boundaryEdge( UInt const i ) const;

    //! i-th mesh 1D Boundary Element Reference.
    /**
     *  Returns a reference to the i-th mesh Boundary Element.
     *
     *  @param i Index of the mesh 1D Boundary Element.
     *  @return Reference to the i-th boundary edge.
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

    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    point_Type & addPoint( bool const boundary = false, bool const vertices = false );

    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param p Point to be added.
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    point_Type & addPoint( point_Type const & p, bool const boundary = false, bool const vertices = false );

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
    point_Type & setPoint( point_Type const & p, UInt const position, bool const boundary = false, bool const vertices = false );

    //! Returns the first mesh Point.
    /**
     *  Returns the first Point in the mesh.
     *
     *  @return Reference to the first mesh Point.
     */
    point_Type & firstPoint();

    //! Returns the last mesh Point.
    /**
     *  Returns the last Point in the mesh.
     *
     *  @return Reference to the last mesh Point.
     */
    point_Type & lastPoint();

    //! Returns the i-th mesh Point.
    /**
     *  Returns the i-th Point in the mesh.
     *
     *  @param i Id of the Point.
     *  @return i-th mesh Point.
     */
    point_Type const & point( UInt const i ) const;

    //! Returns a reference to the i-th mesh Point.
    /**
     *  Returns the i-th Point in the mesh.
     *
     *  @param i Id of the Point.
     *  @return Reference i-th mesh Point.
     */
    point_Type & point( UInt const i );

    //! Returns a reference to the i-th mesh Boundary Point.
    /**
     *  Returns the i-th Boundary Point in the mesh.
     *
     *  @param i Id of the Boundary Point.
     *  @return Reference i-th mesh Boundary Point.
     */
    point_Type const & boundaryPoint( UInt const i ) const;

    //! Returns a reference to the i-th mesh Boundary Point.
    /**
     *  Returns the i-th Boundary Point in the mesh.
     *
     *  @param i Id of the Boundary Point.
     *  @return Reference i-th mesh Boundary Point.
     */
    point_Type & boundaryPoint( UInt const i );

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

    //! Number of Boundary Vertices in Region.
    /**
     *  Returns the number of Boundary Vertices in Region.
     *
     *  @return Number of Vertices in Region.
     */
    UInt  numBVertices() const;

    //! Reference to the counter of Boundary Vertices in Region.
    /**
     *  Allows to change number of Boundary Vertices in Region.
     *
     *  @return Reference to the counter of Boundary Vertices in Region.
     */
    UInt& numBVertices();

    //! Returns the global number of vertices in the mesh.
    /**
     *  Interrogation to the counter of vertices in the mesh.
     *
     *  @return Number of vertices in the mesh.
     */
    UInt numGlobalVertices() const;

    //! Set the number of vertices in the mesh.
    /**
     *  Set the internal counter of vertices points in the mesh.
     *
     *  @param n Number of vertices in the mesh.
     */
    void setNumGlobalVertices( UInt const n ) { M_numGlobalVertices = n; }

    //! Vertex check.
    /**
     *  Is this Point a Vertex?
     *
     *  @param id Point's id.
     *  @return true if the Point is a Vertex, false otherwise.
     */
    bool isVertex        ( UInt const & id )       const;

    //! Vertex check.
    /**
     *  Is this Point a Vertex?
     *
     *  @param p A Point.
     *  @return true if the Point is a Vertex, false otherwise.
     */
    bool isVertex        ( point_Type const & p ) const;

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
    bool isBoundaryPoint ( point_Type const & p ) const;

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
    bool hasLocalEdges() const { return true; }

    //! Edge Id of a certain edge number around a Face.
    /**
     *  Returns always zero. Parameters not used.
     *
     *  @param iface Reference to Face.
     *  @param locE Position of the surrounding edges.
     *  @return Edge Id.
     */
    UInt localEdgeId( const FaceType& iface, const UInt locE ) const { return 0; }

    //! Edge Id of a certain edge number around a Face.
    /**
     *  Returns always zero. Parameters not used.
     *
     *  @param facId Face Id.
     *  @param locE Position of the surrounding edges.
     *  @return Edge Id.
     */
    UInt localEdgeId( const UInt facId, const UInt locE ) const { return 0; }

    //! Builds Edge-To-Face lookup table.
    /**
     *  Builds the lookup table for adjacency queries (not used).
     *
     *  It does nothing in 1D mesh.
     *
     *  @param ce true if edge Counter is set, false otherwise (default).
     *  @param ee Edge count for memory reservation (0 default).
     */
    void updateElementEdges(bool ce=false, UInt ee=0 );

    //! Destroys Edge-To-Face lookup table.
    /**
     *  Destroys the lookup table for adjacency queries (not used).
     *
     *  It does nothing in 1D mesh.
     */
    void cleanElementEdges() {return;}

    //! Global number of Edges.
    /**
     *  Set the internal counter of global number of Edge elements in the mesh.
     *
     *  @param n Number of Edges.
     */
    void setNumGlobalEdges( UInt const n ) { M_numGlobalEdges = n; }

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

    //! Container of mesh Volumes (not used).
    Volumes volumeList;
    //! Container of mesh Points/Vertices.
    Points  pointList;
    //! Container of mesh Faces (not used).
    Faces   faceList;
    //! Container of mesh Edges.
    Edges   edgeList;
    //! Boundary points list.
    SimpleVect<point_Type * > _bPoints;

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
    UInt _numVertices;
    //! Boundary Vertices Number
    UInt _numBVertices;
    //! Points Number
    UInt _numPoints;
    //! Boundary Points Number
    UInt _numBPoints;
    //! Edges Number
    UInt _numEdges;
    //! Boundary Edges Number
    UInt _numBEdges;

    //! Global Vertices Number
    UInt M_numGlobalVertices;
    //! Global Points Number
    UInt M_numGlobalPoints;
    //! Global Edges Number
    UInt M_numGlobalEdges;

    /** @} */ // End of group Internal Counters

}; // End of class RegionMesh1D


// =================================================== //
// =================================================== //
//                    IMPLEMENTATION                   //
// =================================================== //
// =================================================== //


// ===================================================
// Constructors
// ===================================================
template <typename GEOSHAPE, typename MC>
RegionMesh1D<GEOSHAPE, MC>::RegionMesh1D( UInt id ) :
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
    setSwitch( switches );
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

    point_Type * pp = 0;

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

    R1(0,0) =  1.;
    R1(0,1) =  0.;
    R1(0,2) =  0.;
    R1(1,0) =  0.;
    R1(1,1) =  cos(rotate[0]);
    R1(1,2) = -sin(rotate[0]);
    R1(2,0) =  0.;
    R1(2,1) =  sin(rotate[0]);
    R1(2,2) =  cos(rotate[0]);

    R2(0,0) =  cos(rotate[1]);
    R2(0,1) =  0.;
    R2(0,2) =  sin(rotate[1]);
    R2(1,0) =  0.;
    R2(1,1) =  1.;
    R2(1,2) = 0.;
    R2(2,0) = -sin(rotate[1]);
    R2(2,1) =  0.;
    R2(2,2) =  cos(rotate[1]);

    R3(0,0) =  cos(rotate[2]);
    R3(0,1) = -sin(rotate[2]);
    R3(0,2) = 0.;
    R3(1,0) =  sin(rotate[2]);
    R3(1,1) =  cos(rotate[2]);
    R3(1,2) = 0.;
    R3(2,0) =  0;
    R3(2,1) =  0.;
    R3(2,2) = 1.;

    S(0,0) = scale[0];
    S(0,1) = 0.;
    S(0,2) = 0.;
    S(1,0) = 0.;
    S(1,1) = scale[1];
    S(1,2) = 0.;
    S(2,0) = 0.;
    S(2,1) = 0.;
    S(2,2) = scale[2];

    //The total rotation is: R = R1*R2*R3 (as in Paraview we rotate first around z, then around y, and finally around x).
    //We also post-multiply by S to apply the scale before the rotation.
    R = prod( R3, S );
    R = prod( R2, R );
    R = prod( R1, R );

    //Create the 3D translate vector
    boost::numeric::ublas::vector<Real> P(3), T(3);
    T(0) = translate[0];
    T(1) = translate[1];
    T(2) = translate[2];

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
RegionMesh1D<GEOSHAPE, MC>::setSwitch( Switch & sw )
{
    sw.create( "HAS_ALL_EDGES" );
    sw.create( "HAS_BOUNDARY_EDGES" );
    sw.create( "HAS_FACE_TO_EDGES" );
    sw.create( "HAS_BEEN_CHECKED" );
    sw.create( "EDGES_HAVE_ADIACENCY" );
}

template <typename GEOSHAPE, typename MC>
inline
void
RegionMesh1D<GEOSHAPE, MC>::setLinkSwitch( std::string const & _s )
{
    std::ostringstream _err_msg;
    _err_msg << "Switch named " << _s << " is not allowed";
    ASSERT0( switches.set( _s ), _err_msg.str().c_str() );
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh1D<GEOSHAPE, MC>::getLinkSwitch( std::string const & _s ) const
{
    return switches.test( _s );
}

template <typename GEOSHAPE, typename MC>
inline
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
    return EdgeType::S_numLocalVertices;
}

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
RegionMesh1D<GEOSHAPE, MC>::element( UInt const & i )
{
    return edge( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh1D<GEOSHAPE, MC>::ElementType const &
RegionMesh1D<GEOSHAPE, MC>::element( UInt const & i ) const
{
    return edge( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh1D<GEOSHAPE, MC>::BElementType &
RegionMesh1D<GEOSHAPE, MC>::bElement( UInt const & i )
{
    return boundaryEdge( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh1D<GEOSHAPE, MC>::BElementType const &
RegionMesh1D<GEOSHAPE, MC>::bElement( UInt const & i ) const
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
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::addEdge()
{
    return addEdge( EdgeType() );
}

template <typename GEOSHAPE, typename MC>
inline
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
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::setEdge( EdgeType const & v, UInt const pos )
{
    ASSERT_PRE( pos <= edgeList.capacity() , "position requested exceed capacity" <<
                pos << " " << edgeList.capacity() ) ;
    edgeList( pos ) = v;
    edgeList( pos ).id() = pos;
    return edgeList( pos );
}

template <typename GEOSHAPE, typename MC>
inline
void
RegionMesh1D<GEOSHAPE, MC>::setEdgeCounter()
{
    _numEdges = edgeList.size();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType const &
RegionMesh1D<GEOSHAPE, MC>::edge( UInt const i ) const
{
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::edge( UInt const i )
{
    ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
    return edgeList( i );
}

template <typename GEOSHAPE, typename MC>
Real
RegionMesh1D<GEOSHAPE, MC>::edgeLength( const UInt& i ) const
{
    ASSERT_BD( i >= 0 && i < edgeList.size() );

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
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::addEdge( bool const boundary )
{
    return addEdge( EdgeType(), boundary );
}


template <typename GEOSHAPE, typename MC>
inline
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
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::setEdge( EdgeType const & f, UInt position, bool const boundary )
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
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::lastEdge()
{
    return edgeList.back();
}




template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType const &
RegionMesh1D<GEOSHAPE, MC>::boundaryEdge( UInt const i ) const
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
inline
typename RegionMesh1D<GEOSHAPE, MC>::EdgeType &
RegionMesh1D<GEOSHAPE, MC>::boundaryEdge( UInt const i )
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
inline
bool
RegionMesh1D<GEOSHAPE, MC>::hasEdges() const
{
    return ! edgeList.empty();
}

template <typename GEOSHAPE, typename MC>
inline
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
inline
bool
RegionMesh1D<GEOSHAPE, MC>::isBoundaryEdge( EdgeType const & e ) const
{
#ifdef NOT_BDATA_FIRST
    //ASSERT(false,"In this version Boundary edges must be stored first");
    bool isboundary = true;
    for ( UInt k = 1; k <= EdgeType::S_numVertices; ++k )
    {
        isboundary = isboundary & e.point( k ).boundary();
    }
    return isboundary;
#else

    return e.id() <= _numBEdges;
#endif
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh1D<GEOSHAPE, MC>::isBoundaryEdge( UInt const & id ) const
{
    return isBoundaryEdge( edge( id ) );
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh1D<GEOSHAPE, MC>::isFullEdge( const UInt & id ) const
{
    return edgeList.size() >= id;
}

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh1D<GEOSHAPE, MC>::edgeElement( UInt const i, UInt const Pos ) const
{
    ASSERT_PRE( i <= edgeList.size(), "Not enough edges stored" ) ;
    ASSERT_BD( i > 0 ) ;
    return edgeElement( edge( i ), Pos );
};

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh1D<GEOSHAPE, MC>::edgeElement( EdgeType const & f, UInt const Pos ) const
{
    ASSERT_BD( ! edgeList.empty() ) ;
    ASSERT_PRE( Pos == 1 || Pos == 2 , "Wrong position (1 or 2)" ) ;
    if ( Pos == 1 )
    {
        return f.firstAdjacentElementIdentity();
    }
    else
    {
        return f.secondAdjacentElementIdentity();
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
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type &
RegionMesh1D<GEOSHAPE, MC>::addPoint( bool const boundary, bool const vertex )
{
    return addPoint( point_Type(), boundary, vertex );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type &
RegionMesh1D<GEOSHAPE, MC>::addPoint( point_Type const & p, bool const boundary, bool const /*vertex*/ )
{
    ASSERT_PRE( pointList.size() < pointList.capacity(), "Point list size exceeded" <<
                pointList.size() + 1 << " " << pointList.capacity() ) ;
    pointList.push_back( p );
    point_Type * pp = & pointList.back();
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
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type &
RegionMesh1D<GEOSHAPE, MC>::setPoint
( point_Type const & p, UInt position, bool const boundary, bool const vertex )
{
    ASSERT_PRE( position <= pointList.capacity(), "Position  exceed lpoint list capacity" <<
                position << " " << pointList.capacity() ) ;
    bool found( false );
    pointList( position ) = p;
    point_Type * pp = & pointList( position );
    pp->setId(position);
    if ( boundary )
    {
        pp->boundary() = true;
        // This is rather complex, since I do not know a priori
        // if point was already stored in the list!
        // No way to avoid it, sorry

        for ( typename SimpleVect<point_Type *>::iterator bp = _bPoints.begin(); bp != _bPoints.end(); ++bp )
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
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type &
RegionMesh1D<GEOSHAPE, MC>::firstPoint()
{
    return pointList.front();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type &
RegionMesh1D<GEOSHAPE, MC>::lastPoint()
{
    return pointList.back();
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type const &
RegionMesh1D<GEOSHAPE, MC>::point( UInt const i ) const
{
    ASSERT_BD( i > 0 && i <= pointList.size() ) ;
    return pointList( i );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type &
RegionMesh1D<GEOSHAPE, MC>::point( UInt const i )
{
    ASSERT_BD( i > 0 && i <= pointList.size() ) ;
    return pointList( i );
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type const &
RegionMesh1D<GEOSHAPE, MC>::boundaryPoint( UInt const i ) const
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh1D<GEOSHAPE, MC>::point_Type &
RegionMesh1D<GEOSHAPE, MC>::boundaryPoint( UInt const i )
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
inline
bool
RegionMesh1D<GEOSHAPE, MC>::isVertex( point_Type const & p ) const
{
    return p.id() <= _numVertices;
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh1D<GEOSHAPE, MC>::isVertex( UInt const & id ) const
{
    return id <= _numVertices;
}


template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh1D<GEOSHAPE, MC>::isBoundaryPoint( UInt const & id ) const
{
    return point( id ).boundary();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh1D<GEOSHAPE, MC>::isBoundaryPoint( point_Type const & p ) const
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

} // End of namespace LifeV

#endif //REGIONMESH1D_H
