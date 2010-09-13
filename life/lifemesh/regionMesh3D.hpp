/*
  This file is part of the LifeV library

  Authors: Luca Formaggia
  Miguel Fernandez

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
/*!
  \file regionMesh3D.hpp
  \brief The mesh classes interfaces
  Introduces the RegionMesh3D class
*/

#ifndef _REGIONMESH3D_HH_
#define _REGIONMESH3D_HH_

#include <cstdlib>
#include <iomanip>
#include <fstream>

#include <life/lifecore/life.hpp>
#include <life/lifemesh/geoElement.hpp>
#include <life/lifecore/switch.hpp>
#include <life/lifemesh/bareItems.hpp>

#include <life/lifearray/vecUnknown.hpp>
#include <life/lifearray/SimpleVect.hpp>

#include <life/lifemesh/basisElSh.hpp>
#include <life/lifemesh/mesh_util_base.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef HAVE_MPI
//headers useful only for reordering:
#include "mpi.h"
#include <parmetis.h>
#endif

namespace LifeV
{
    /**
       \class RegionMesh3D
       \brief The Region Mesh Class

       This is the class that stores the mesh entities for a single region In
       a region elements are all of the same type

       \author Luca Formaggia
       \author Miguel Fernandez
    */
    template <typename GEOSHAPE, typename MC = DefMarkerCommon >
    class RegionMesh3D
        :
        public MeshEntity,
        public MC::RegionMarker
    {

    public:
        explicit RegionMesh3D();
        //! Constructor providing ID
        explicit RegionMesh3D( ID id );
        ~RegionMesh3D<GEOSHAPE, MC>();

        //! \name Markers Types
        /*! From Marker Common (MC) template parameter */
        //@{
        typedef MC MarkerCommon;
        typedef typename MC::PointMarker PointMarker;
        typedef typename MC::EdgeMarker EdgeMarker;
        typedef typename MC::FaceMarker FaceMarker;
        typedef typename MC::VolumeMarker VolumeMarker;
        typedef typename MC::RegionMarker RegionMarker;
        typedef typename MC::RegionMarker Marker;
        //@}

        //! \name Basic Element Shapes
        //@{
        typedef GEOSHAPE VolumeShape;
        typedef typename GEOSHAPE::GeoBShape FaceShape;
        typedef typename FaceShape::GeoBShape EdgeShape;
        //@}

        //! \name Geometric Entities
        /*! Are the ones actually stored in the containers
         */
        //@{
        typedef GeoElement3D<GEOSHAPE, MC>  VolumeType;
        typedef GeoElement2D<FaceShape, MC> FaceType;
        typedef GeoElement1D<EdgeShape, MC> EdgeType;
        typedef GeoElement0D<MC>            PointType;
        //@}

        //! \name GeoElemen Containers
        /*! Typedefs for STL compliant containers of mesh geometric entities
          I Use SimpleVect<> for addressing from 1. See SimpleVect.h
        */
        //!@{
        //! Points Container
        typedef SimpleVect<PointType>   Points;
        //! Elements Container
        typedef SimpleVect<VolumeType > Volumes;
        //! Faces Container: it may contain only Boundary faces
        typedef SimpleVect<FaceType>    Faces;
        //! Edges Container: it may be empty
        typedef SimpleVect<EdgeType>    Edges;
        //!@}

        /*! \name Generic_Types
         * Generic types for all regionmeshXX
         * These are part of the generic generic interface common
         * for all RegionMeshes (3D -- 1D).
         * The idea is that in this way it is simpler to write functions operating on 2D
         * and 3D meshes
         */
        //@{
        typedef GEOSHAPE ElementShape;
        typedef typename GEOSHAPE::GeoBShape BElementShape;
        typedef GeoElement3D<GEOSHAPE, MC> ElementType;
        typedef GeoElement2D<FaceShape, MC> BElementType;
        typedef SimpleVect<VolumeType > Elements;
        //! Faces Container: it may contain only Boundary faces
        typedef SimpleVect<FaceType> BElements;
        //@}

        /*! \name Switch  Methods
         *   A switch is used to store the status of the RegionMesh
         *   It is  used internally to control whether
         *   some data structures have been set up.
         *
         *   The possible values for the "switches" are
         *
         *   \verbatim

         HAS_ALL_FACES            HAS_VOLUME_TO_FACES
         HAS_ALL_EDGES            HAS_VOLUME_TO_EDGES
         HAS_BOUNDARY_FACES       HAS_BEEN_CHECKED
         HAS_BOUNDARY_EDGES       FACES_HAVE_ADIACENCY

         \endverbatim
         *
         */
        //@{
        //! Returns the number of switches which have been set
        /*const*/ UInt numSwitches() const
            {
                return switches.size();
            };
        //! Interrogate Switch.
        /*! Example: getLinkSwitch("HAS_ALL_FACES")*/
        bool getLinkSwitch( std::string const & _s ) const;
        //! Set a switch
        /*! Example: setLinkSwitch("HAS_ALL_FACES")*/
        void setLinkSwitch( std::string const & _s );
        //! uset a switch
        /*! Example: unsetLinkSwitch("HAS_ALL_FACES")*/
        void unsetLinkSwitch( std::string const & _s );
        //! output switches contents
        std::ostream & showLinkSwitch( bool verbose = false, std::ostream & out = std::cout )
            {
                return switches.showMe( verbose, out );
            }
        //@}



        /** @name Mesh transformation
         */
        //@{

    	//! Utilities to implement mesh movement \author Miguel Fenrandez:11/2002
        /*! Move the mesh according to a given displacement stored in disp.
          Disp is a 3*numpoints() vector which stores the x-displacement first, then the y-displacements etc.
          The VECTOR object must have a size() and a standard [] addressing operator.
         */
        template <typename VECTOR>
        void moveMesh( const VECTOR & disp, int dim );

    	//! Transform the mesh using boost::numeric::ublas; \author Cristiano Malossi:14/09/2009
        /*! Scale, rotate and translate the mesh (operations performed in this order!).
         *  NOTES:
         *    -  Rotation follows Paraview conventions: first rotate around z-axis,
         *       then around y-axis and finally around x-axis;
         *    -  All the vectors must allow the command: operator[];
         * \param scale        - vector of three components for (x,y,z) scaling of the mesh
         * \param rotate       - vector of three components (radiants) for rotating the mesh
         * \param translate    - vector of three components for (x,y,z) translation the mesh
         */
        template <typename VECTOR>
        void transformMesh( const VECTOR& scale, const VECTOR& rotate, const VECTOR& translate );

        //@}

        //! Get the maximum H over all the edges of the mesh; @author Cristiano Malossi: 22/02/2010
        /*!
         * \return maximum H
         */
        Real maxH() const;

        //! Get the minumum H over all the edges of the mesh; @author Cristiano Malossi: 27/04/2010
        /*!
         * @return minumum H
         */
        Real minH() const;

        //! Get the mean H over all the edges of the mesh; @author Cristiano Malossi: 22/02/2010
        /*!
         * \return maximum H
         */
        Real meanH() const;

        UInt numLocalVertices() const; //!< Number of local faces for each 3Delement
        UInt numLocalFaces() const; //!< Number of local faces for each 3Delement
        UInt numLocalEdges() const;  //!< Number of local edges for each 3DElement
        UInt numLocalEdgesOfFace() const; //!< Number of edges on each face

        /*! \name Generic Methods Generic methods for all regionmeshXX
         * These are the generic methods to get information about the number of elements.
         * It is a generic interface common for all RegionMeshes (3D -- 1D).
         */
        //@{
        UInt numGlobalElements() const; //!< Number of global elements (volumes)
        UInt & numGlobalElements(); //!< Number of global elements (volumes)
        UInt numElements() const; //!< Number of elements (volumes)
        UInt & numElements(); //!< Number of elements (volumes)
        UInt numBElements() const; //!< Number of Boundary elements (faces)
        UInt & numBElements(); //!< Number of boundary elements (faces)
        ElementType & element( ID const & i );
        ElementType const & element( ID const & i ) const;
        BElementType & bElement( ID const & i );
        BElementType const & bElement( ID const & i ) const;
        //@}
        /* ============================================
           Volume Related Methods
           ============================================*/
        /*! \name Volume Methods

        \anchor volume_methods
        Methods which operates on 3D elements

        There are different way of counting gemetry entities

        <UL>

        <LI><b>Counter</b> A counter stored the number of entities (element, faces
        etc) in the mesh. It does NOT necessarily correspond to the numer of
        entities actually stored. Indeed, I may not stor edges, yet having a counter with
        the number of edges in the mesh.</LI>

        <li><b>Stored. The number of stored entities corresponds to the size()
        of the container.  A container may store all mesh entities (this is the
        default for the Points and the Volumes), only the boundary entities
        (this is the default with the Faces) or even none (this is the default
        with the edges).

        <LI><b>Capacity</b>. The maximum number of entities that may stored
        before the container is resized.

        <\UL>
        */
        //@{
        /*! Returns number of Volume elements in the mesh as given by the internal counter.*/
        UInt /*const*/ numVolumes() const;
        UInt /*const*/ numGlobalVolumes() const;
//        UInt & numVolumes();    /*!< Access number of Volumes (internal counter)*/
        UInt storedVolumes() const; //!< volumes actully stored in list
        /*! Current capacity of Volumes  Container, i.e. how many elements may be stored */
        UInt maxNumVolumes() const;
        /*! Changes Current capacity of Volumes  (Optionally sets internal counter.) */
        void setMaxNumVolumes      ( UInt const n, bool const setcounter = false );
        void setMaxNumGlobalVolumes( UInt const n );

        void setNumVolumes      ( UInt const n );
        /*! Adds volumes.
          Volume Id is computed automatically. Return ref to added Volume
        */
        VolumeType & addVolume();
        VolumeType & addVolume( VolumeType const & v );
        VolumeType & setVolume( VolumeType const & v, ID const pos );
        void setVolumeCounter(); //<! set numVolumes counter
        VolumeType & lastVolume(); //!< Reference to last volume stored in list. Useful for mesh readers

        VolumeType const & volume( ID const i ) const; //!< ith mesh 3Delement
        VolumeType & volume( ID const i ); //!<ith mesh 3Delement
        //@}

        /*! \name Element Adiacency Methods

        Methods to obtain the ID of Face and Edge belonging to an element.
        Accessing this information requires that the appropriate data
        structures have been set by using the updateElementEdges() or
        updateElementFaces() methods.

        It is NOT required to have the full information about edges and faces:
        The ID of the Face and Edge entities may be calculated withoud
        contructing the corresponding Edge of Face Object. This saves
        memory. However the  methods has the capability of
        building the actual internal Faces (the boundary GeoFaces must
        already exist!) of the Edges. In this case the MarkerFlag is inherited using the rules stated in the
        marker_traits class.
        */
        //@{
        //! Is the array for local Faces set up?
        /*! It does not use switches, but interrogates the container directly*/
        bool hasLocalFaces() const;
        //! Build localFaceId table and optionally fills the list of Faces
        /*!  \param createFaces is set true if we want also to create the actual list
          of internal faces. There is another utility (in mesh_util.h), which
          might be used for the same purpose if we want just to create the faces
          and not also the LocaFaceID table.

          \param estimateFaceNumber is a guess provided by the user of the total
          number of faces. It is relevant only when createFaces=true. Setting it
          to a proper value helps in reducing time and memory.

          \pre The routine assumes that the boundary faces are properly set, if not use the
          methods in mesh_util.h.
        */
        void updateElementFaces( bool createFaces = false, const bool verbose = false, UInt estimateFaceNumber = 0 );
        //! Destroys element-to-face container. Useful to save memory!
        void cleanElementFaces();
        //!local Face
        /*! \param volId Id of volume (element)
          \param locF local face number 1<LocF<=numLocalFaces()
          \return ID of the face
        */
        ID localFaceId( ID const volId, ID const locF ) const;
        //!local Face
        /*! \param iv Reference to a volume (element)
          \param locF local face number 1<LocF<=numLocalFaces()
          \return ID of the face
        */
        ID localFaceId( const VolumeType & iv, ID const locF ) const;
        bool hasLocalEdges() const;  //!<Same for Edges
        //! Build localEdgeId table and optionally fills the list of Edges
        /*!  \param createEdges is set true if we want also to create the actual list
          of edges. There is another utility (in mesh_util.h), which
          might be used for the same purpose if we want just to create the faces
          and not also the LocalEdgeID table.

          \param estimateEdgeNumber is a guess provided by the user of the total
          number of edges. It is relevant only when createFaces=true. Setting it
          to a proper value helps in reducing time and memory.

          \note This method does not assume that bounrary edges are stores, since
          this condition is NOT a a paradigm for a RegionMesh3D.
        */
        void updateElementEdges( bool createEdges = false, const bool verbose = false, UInt estimateEdgeNumber = 0 );
        //! Destroys element-to-face container
        void cleanElementEdges();
        //!local Edge
        /*! \param volId Id of volume (element)
          \param locE local edge number 1<LocE<=numLocalEdges()
          \return ID of the edge
        */
        UInt localEdgeId( UInt const volId, UInt const locE ) const;
        UInt localEdgeId( const VolumeType & iv, UInt const locE ) const;
        //@}

        /*! \name Face Methods

        Methods to access/create/modify faces data  \anchor face_methods

        There are different way of counting Faces:

        <UL>

        <LI>Number of Faces: is the declared number of total faces in the mesh.
        It uses an internal counter.

        <b>WARNING></b>: A value different from zero does NOT imply that the
        faces are actually stored. This counter declares the number of Faces
        that the mash have not those that are actually stored in the list of
        faces!</LI>

        <LI>Number of Boundary Faces: is the declared number of boundary faces
        in the mesh. It uses an internal counter.

        <b>WARNING</b> A value different from zero does NOT imply that the
        boundary faces are actually stored.</LI>

        <LI>Number of Stored Faces: It is the number of Faces actually stored
        on the face container.</LI>

        <LI>Maximum number of stored faces. The number of faces that may stored
        before the container is resized.

        <b>Important:</b> This parameter has to be set BEFORE inserting faces
        in the container if we want that pointer into the container maintains
        their validity. The container will also have a better performance.</LI>
        <\UL>
        See also \ref volume_methods
        */

        //@{
        UInt /*const*/ numFaces() const; //!<Number of total faces in mesh (uses counter).
        UInt /*const*/ numGlobalFaces() const; //!<Number of total faces in mesh (uses counter).
//        UInt & numFaces();    //!<Number of total faces (may modify)
        UInt storedFaces() const; //!< Number of stored faces
        UInt maxNumFaces() const;  //!< Max number of faces that can be stored
        //! Set Maximum Number of Faces in the Faces container
        void setMaxNumFaces( UInt const n, bool const setcounter = false );
        void setMaxNumGlobalFaces( UInt const n );
        //! Adds a face to list (optionally a boundary face)
        FaceType & addFace( bool const boundary = false );
        //! Adds a face (optionally a boundary face) and adjourn its ID
        FaceType & addFace( FaceType const & f, bool const boundary = false );
        //! Adds a face (optionally a boundary face) and adjourn its ID
        FaceType & setFace( FaceType const & f, ID position, bool const boundary = false );
        //! The last face in the container
        FaceType & lastFace();
        FaceType const & face( ID const i ) const; //<! ith mesh face
        FaceType & face( ID const i ); //!< ith mesh face (may  modify!)
        FaceType const & boundaryFace( ID const i ) const; //!< ith boundary face.
        FaceType & boundaryFace( ID const i ); //!< ith boundary  face.
        void setNumFaces( UInt const n ) ; //<! Set counter of boundary faces
        void setNumBFaces( UInt const n ) ; //<! Set counter of boundary faces
        bool hasFaces() const;  //<! Do I store mesh faces?
        bool hasInternalFaces() const; //<! Do I store also internal faces?
        UInt numBFaces() const;  //<! Number of Boundary Faces
        bool isBoundaryFace( FaceType const & f ) const;  //<!Is this face on boundary?
        bool isBoundaryFace( ID const & id ) const;  //<!Is  face whose ID is id on boundary?

        /*! Does this ID corresponds to a full face? A FULL FACE is a 2DElement  that is actually
          stored in the Face container.
          \return false if id does not corresponfd to a boundary face and internal faces are not stored
        */
        bool isFullFace( UInt const & id ) const;
        //@}

        /*@{*/
        //!ID of the Volume Element adjacent to a Face.
        /*!  The first element is the one <em>ORIENTED coherently with the
          face</em> (AS STORED in Faces). It means that the face orientation is
          OUTWARD with respect to the element. The second element is either null
          (boundary face) or indicates that the normal of the face appears INWARD
          with respect to that element

          \param Facid ID of the face
          \param Pos is equal to 1 or 2 and indicates   first or second element.
          \return ID of adjacent volume element or 0 if none.
        */
        ID faceElement( ID const faceId, UInt const Pos ) const;
        //! 3DElement adjacent to a FACE. Face reference given.
        ID faceElement( FaceType const & f, UInt const Pos ) const;
        /*@}*/

        /*! \name  Edge Methods

        Methods to access/create/modify edges data

        There are different point counters which may bi interrogated or set:

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
        have a better performance.</LI> </UL>

        \note To have more information on
        the Edge methods look at the documentation of th eanalogous Faces
        methods \ref face_methods and uin \ref volume_methods */

        /*                          EDGES METHODS                           */
        //@{
        UInt /*const*/ numEdges() const;
        UInt /*const*/ numGlobalEdges() const;
//        UInt & numEdges();
        UInt storedEdges() const;
        UInt maxNumEdges() const;
        void setMaxNumEdges      ( UInt const n, bool const setcounter = false );
        void setMaxNumGlobalEdges( UInt const n );
        //
        EdgeType & addEdge( bool const boundary = false );
        EdgeType & addEdge( EdgeType const & f, bool const boundary = false );
        EdgeType & setEdge( EdgeType const & f, ID position, bool const boundary = false );
        EdgeType & lastEdge();
        EdgeType const & edge( ID const i ) const; //!< ith mesh edge
        EdgeType & edge( ID const i ); // ith mesh edge
        EdgeType const & boundaryEdge( ID const i ) const; //!< ith b. edge.
        EdgeType & boundaryEdge( ID const i ); //!< ith b. edge.
        UInt numBEdges() const;
        void setNumEdges ( UInt const n);
        void setNumBEdges( ID const n ) ; //<! Set counter of boundary edges
        bool hasEdges() const;         //<! Do I store edges?
        bool hasInternalEdges() const;         //<! Do I store internal edges?
        bool isBoundaryEdge( EdgeType const & e ) const;  //<!Is this edge on boundary?
        bool isBoundaryEdge( ID const & id ) const;  //<!Is  edge with ID id on boundary?
        bool isFullEdge( ID const & id ) const;  //<!Is edge "id" a full edge?
        //@}

        /*                           POINTS/VERTICES METHODS */
        /*! \name PointMethods

        Methods to access/create/modify Points/Vertices data

        There are different Point counters which may be interrogated or set
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
        container. since all mesh entities (a part the Point)
        will contain POINTERS into the Points container!
        </LI> </UL>

        \note To have more information on the Point methods look at the
        documentation of the analogous Faces method (\ref face_methods).
        */
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
        void setMaxNumPoints      ( UInt const n, bool const setcounter = false );
        void setMaxNumGlobalPoints( UInt const n );
        //
        void setNumVertices       ( UInt const n );
        void setNumGlobalVertices ( UInt const n );
        void setNumBVertices      ( UInt const n );
        /*! Add a point
          Point ID is computed automatically. We may speciy isf the point is on the boundary or
          if it is a vertex */
        PointType & addPoint( bool const boundary = false, bool const vertices = false );
        /*! Add a point
          This version taked a reference to a point
        */
        PointType & addPoint( PointType const & p, bool const boundary = false, bool const vertices = false );
        /*! Add a point
          This method is for advanced use only
        */
        PointType & setPoint( PointType const & p, ID const position, bool const boundary = false, bool const vertices = false );
        /*! Add a point
          This method is for advanced use only
        */
        PointType & setPoint( ID const & position, bool const boundary = false, bool const vertices = false );

        //!< adds point
        UInt addPoint( ID const iden, bool const boundary = false, UInt const start = 1 );
        PointType & lastPoint(); //!< last mesh point
        PointType const & point( ID const i ) const; //!< ith mesh point
        PointType & point( ID const i ); //!< ith mesh point
        PointType const & boundaryPoint( ID const i ) const; //!< ith b. point.
        PointType & boundaryPoint( ID const i ); //!< ith b. point.
        UInt numBPoints() const; //!< counter of boundary points
        void setNumBPoints( UInt const n ); //<! Sets counter of boundary points
        UInt /*const*/ numVertices() const; //!< Number of vertices in Region
        UInt /*const*/ numGlobalVertices() const; //!< Number of vertices in Region
//        UInt & numVertices(); //!< Allows to change number of vertices in Region
        UInt numBVertices() const; //!< Number of boundary vertices in RegionMesh
        UInt & numBVertices();
        bool isVertex( ID const & id ) const;  //<!Is this point a Vertex?
        bool isVertex( PointType const & p ) const;  //<!Is this point a Vertex?
        bool isBoundaryPoint( PointType const & p ) const;  //<!Is this point on boundary?
        bool isBoundaryPoint( ID const & id ) const;  //<!Is this point on boundary?
        //@}


        /*! \name Utilities
          Useful methods
        */
        //@{
        //! Extracts from the mesh a list of entities matchin a marker flag
        /*!
          It adds the IDs of the geometric entities matching an EntityFlag
          to vector of IDs.
          The entity to extract is defined through the ReferenceGeometry enum:
          {VERTEX=0, EDGE = 1, FACE = 2, VOLUME = 3};
        */
        void extractEntityList(std::vector<ID> &, ReferenceGeometry const & r, EntityFlag const &) const;

        //! Prints some mesh info
        std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const;
        //! Basic tests for mesh consistency.
        /*! For more estensive test see \link mesh_util.h */
        int check(int level = 0, bool const fix = false, bool const verbose = true,
                  std::ostream & out = std::cerr );
        //@}

        /*! \name RegionContainers
          Stl compliant containers for basic structure I
          expose them since they are standard containers
        */
        //@{
        Points pointList; //!< Container of mesh points/vertices
        /*! Container of mesh points/vertices (mesh movement)

        Used only by mesh node movement routines. It contains the
        mesh nodes with the previous value
        */
        Points _pointList;
        Volumes volumeList; //!< Container of mesh 3DElements
        Faces faceList;  //!< Container of mesh  Faces
        Edges edgeList;  //!< Container of mesh Edges.
        SimpleVect<PointType * > _bPoints; //!< Boundary points list
        //@}

        //! Switches
        // \sa Switch
        Switch switches;
        void list_of_points( bool ( *fct ) ( double, double, double ), std::vector<UInt>& list_pts );

        std::map<int,int> & globalToLocalNode(){return M_globalToLocalNode;}
        std::map<int,int> & localToGlobalNode(){return M_localToGlobalNode;}
#ifdef HAVE_MPI
        void orderMesh(MPI_Comm comm);
#endif
        void printLtGMap(std::ostream & os);
        void edgeMarkers(std::map<ID, ID> const& locDof, UInt newMarker);

    private:

        /*! Arrays containing the ids of Edges and Faces of each element
          I use a Define to use localto global array or directly the
          bareedges */
        SimpleArray<UInt> _VToF;
        SimpleArray<UInt> _VToE;

        /*! \name help-functions
          function prototipes for handling lists \todo
          make them global and friend to RegionmeshXX , add addItem, setItem and
          lastItem method
        */
        //@{
        template < typename T >
        UInt numItems( SimpleVect< T> const & list ) const;
        template < typename T >
        UInt maxNumItems( SimpleVect< T> const & list ) const;
        template < typename T >
        void setMaxNumItems( SimpleVect< T> & list, UInt n, std::string title );
        //@}

        // Internal counters
        UInt M_numVolumes;

        UInt M_numVertices;
        UInt M_numBVertices;

        UInt M_numPoints;
        UInt M_numBPoints;

        UInt M_numFaces;
        UInt M_numBFaces;

        UInt M_numEdges;
        UInt M_numBEdges;

        UInt M_numGlobalPoints;
        UInt M_numGlobalVertices;
        UInt M_numGlobalEdges;
        UInt M_numGlobalFaces;
        UInt M_numGlobalVolumes;


        // Modif Miguel:11/2002
        bool M_moved; // True whenever the mesh has been moved

        //
    private:
        //! Copy contructor: not implemented yet!
        explicit RegionMesh3D( RegionMesh3D<GEOSHAPE, MC> const & m );

        RegionMesh3D<GEOSHAPE, MC> operator=( RegionMesh3D<GEOSHAPE, MC> const & m );

        std::map<int, int>      M_globalToLocalNode;
        std::map<int, int>      M_localToGlobalNode;

        std::map<int, int>      M_globalToLocalEdge;

        std::map<int, int>      M_globalToLocalFace;

        std::map<int, int>      M_globalToLocalVolume;

    };


    /* ---------------------------------------------------------------------
       RegionMesh3D Implementations
       -----------------------------------------------------------------------*/
    void set_switches_for_regionmesh( Switch & sw );

    template <typename GEOSHAPE, typename MC>
    RegionMesh3D<GEOSHAPE, MC>::RegionMesh3D() :
        MeshEntity(),
        MC::RegionMarker(),
        switches(),
        M_numVolumes( 0 ),
        M_numVertices( 0 ),
        M_numBVertices( 0 ),
        M_numPoints( 0 ),
        M_numBPoints( 0 ),
        M_numFaces( 0 ),
        M_numBFaces( 0 ),
        M_numEdges( 0 ),
        M_numBEdges( 0 ),
        M_moved( 0 ),
        M_globalToLocalNode(),
        M_globalToLocalEdge(),
        M_globalToLocalFace(),
        M_globalToLocalVolume()
    { //Modif Miguel:11/2002
        set_switches_for_regionmesh( switches );
    }


    template <typename GEOSHAPE, typename MC>
    RegionMesh3D<GEOSHAPE, MC>::RegionMesh3D( ID id ) :
        MeshEntity( id ),
        MC::RegionMarker(),
        switches(),
        M_numVolumes( 0 ),
        M_numVertices( 0 ),
        M_numBVertices( 0 ),
        M_numPoints( 0 ),
        M_numBPoints( 0 ),
        M_numFaces( 0 ),
        M_numBFaces( 0 ),
        M_numEdges( 0 ),
        M_numBEdges( 0 ),
        M_moved( 0 )
    { //Modif Miguel:11/2002
        set_switches_for_regionmesh( switches );
    }


    template <typename GEOSHAPE, typename MC>
    RegionMesh3D<GEOSHAPE, MC>::RegionMesh3D( RegionMesh3D<GEOSHAPE, MC> const & m ):
        MeshEntity      (m),
        MC::RegionMarker(m),
        pointList       (m.pointList),
        _pointList      (m._pointList),
        volumeList      (m.volumeList),
        faceList        (m.faceList),
        edgeList        (m.edgeList),
        _bPoints        (m._bPoints),
        switches        (m.switches),
        _VToF           (m._VToF),
        _VToE           (m._VToE),
        // Internal counters
        M_numVolumes     (m.M_numVolumes),
        M_numVertices    (m.M_numVertices),
        M_numBVertices   (m.M_numBVertices),
        M_numPoints      (m.M_numPoints),
        M_numBPoints     (m.M_numBPoints),
        M_numFaces       (m.M_numFaces),
        M_numBFaces      (m.M_numBFaces),
        M_numEdges       (m.M_numEdges),
        M_numBEdges      (m.M_numBEdges),
        M_moved          (m.M_moved)
    {
        ASSERT0( true, "Copy Costructor Not Yet Implemented for RegionMesh3D" ) ;
    }

    template <typename GEOSHAPE, typename MC>
    RegionMesh3D<GEOSHAPE, MC>
    RegionMesh3D<GEOSHAPE, MC>::operator=( RegionMesh3D<GEOSHAPE, MC> const & m )
    {
    }

    template <typename GEOSHAPE, typename MC>
    RegionMesh3D<GEOSHAPE, MC>::~RegionMesh3D()
    {
    }



    //!< Move the mesh from a given displacement
    template <typename GEOSHAPE, typename MC>
    template <typename VECTOR>
    void RegionMesh3D<GEOSHAPE, MC>::moveMesh( const VECTOR & disp, int dim )
    {

//         if ( disp.size() != nDimensions * M_numPoints )
//             ERROR_MSG( "We can not move the mesh with this displacement" );

        if ( !M_moved )
        { // We store the reference position of the mesh
            M_moved = 1;
            _pointList.reserve( M_numPoints );
            for ( typename Points::const_iterator ip = pointList.begin();ip != pointList.end();++ip )
                _pointList.push_back( *ip );
        }

        for ( unsigned int i = 0; i < pointList.size();++i )
        {
            for ( ID j = 1; j <= nDimensions; ++j )
            {
                int id = pointList[i].id();
                if ( disp.BlockMap().LID(id + dim*(j - 1)) >= 0 )
                {
                    pointList[ i ].coordinate( j ) = _pointList[ i ].coordinate( j ) + disp[ ( j - 1 ) * dim + id ];
                }
            }
        }
    }

    template <typename GEOSHAPE, typename MC>
    template <typename VECTOR>
    void RegionMesh3D<GEOSHAPE, MC>::transformMesh( const VECTOR& scale, const VECTOR& rotate, const VECTOR& translate )
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
    Real RegionMesh3D<GEOSHAPE, MC>::maxH() const
    {
        Real MaxH(0);
        Real deltaX(0), deltaY(0), deltaZ(0);

        for ( UInt i(0); i < static_cast<UInt> ( edgeList.size() ); ++i )
        {
            deltaX = ( edgeList[ i ].point( 2 ) ).x() - ( edgeList[ i ].point( 1 ) ).x();
            deltaY = ( edgeList[ i ].point( 2 ) ).y() - ( edgeList[ i ].point( 1 ) ).y();
            deltaZ = ( edgeList[ i ].point( 2 ) ).z() - ( edgeList[ i ].point( 1 ) ).z();

            deltaX *= deltaX;
            deltaY *= deltaY;
            deltaZ *= deltaZ;

            MaxH = std::max( MaxH, deltaX+deltaY+deltaZ );
        }

        return std::sqrt( MaxH );
    }

    template <typename GEOSHAPE, typename MC>
    Real RegionMesh3D<GEOSHAPE, MC>::minH() const
    {
        Real MinH( 1E10 );
        Real deltaX(0), deltaY(0), deltaZ(0);

        for ( UInt i(0); i < static_cast<UInt> ( edgeList.size() ); ++i )
        {
            deltaX = ( edgeList[ i ].point( 2 ) ).x() - ( edgeList[ i ].point( 1 ) ).x();
            deltaY = ( edgeList[ i ].point( 2 ) ).y() - ( edgeList[ i ].point( 1 ) ).y();
            deltaZ = ( edgeList[ i ].point( 2 ) ).z() - ( edgeList[ i ].point( 1 ) ).z();

            deltaX *= deltaX;
            deltaY *= deltaY;
            deltaZ *= deltaZ;

            MinH = std::min( MinH, deltaX+deltaY+deltaZ );
        }

        return std::sqrt( MinH );
    }

    template <typename GEOSHAPE, typename MC>
    Real RegionMesh3D<GEOSHAPE, MC>::meanH() const
    {
        Real MeanH = 0;
        Real deltaX(0), deltaY(0), deltaZ(0);

        for ( UInt i(0); i < static_cast<UInt> ( edgeList.size() ); ++i )
        {
            deltaX = ( edgeList[ i ].point( 2 ) ).x() - ( edgeList[ i ].point( 1 ) ).x();
            deltaY = ( edgeList[ i ].point( 2 ) ).y() - ( edgeList[ i ].point( 1 ) ).y();
            deltaZ = ( edgeList[ i ].point( 2 ) ).z() - ( edgeList[ i ].point( 1 ) ).z();

            deltaX *= deltaX;
            deltaY *= deltaY;
            deltaZ *= deltaZ;

            MeanH += deltaX+deltaY+deltaZ;
        }

        return std::sqrt( MeanH / static_cast<Real> ( edgeList.size() ) );
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numLocalVertices() const
    {
        return VolumeType::numLocalVertices;
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numLocalFaces() const
    {
        return VolumeType::numLocalFaces;
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numLocalEdges() const
    {
        return VolumeType::numLocalEdges;
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numLocalEdgesOfFace() const
    {
        return FaceType::numLocalEdges;
    }

    // Generic Methods
    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numElements() const
    {
        return M_numVolumes;
    }

    template <typename GEOSHAPE, typename MC>
    UInt & RegionMesh3D<GEOSHAPE, MC>::numElements()
    {
        return M_numVolumes;
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numGlobalElements() const
    {
        return M_numGlobalVolumes;
    }

    template <typename GEOSHAPE, typename MC>
    UInt & RegionMesh3D<GEOSHAPE, MC>::numGlobalElements()
    {
        return M_numGlobalVolumes;
    }

    template <typename GEOSHAPE, typename MC>
    UInt RegionMesh3D<GEOSHAPE, MC>::numBElements() const
    {
        return M_numBFaces;
    }

    template <typename GEOSHAPE, typename MC>
    UInt & RegionMesh3D<GEOSHAPE, MC>::numBElements()
    {
        return M_numBFaces;
    }

    template <typename GEOSHAPE, typename MC>
    typename RegionMesh3D<GEOSHAPE, MC>::ElementType &
    RegionMesh3D<GEOSHAPE, MC>:: element( ID const & i )
    {
        return volume( i );
    }

    template <typename GEOSHAPE, typename MC>
    typename RegionMesh3D<GEOSHAPE, MC>::ElementType const &
    RegionMesh3D<GEOSHAPE, MC>:: element( ID const & i ) const
    {
        return volume( i );
    }

    template <typename GEOSHAPE, typename MC>
    typename RegionMesh3D<GEOSHAPE, MC>::BElementType &
    RegionMesh3D<GEOSHAPE, MC>:: bElement( ID const & i )
    {
        return boundaryFace( i );
    }

    template <typename GEOSHAPE, typename MC>
    typename RegionMesh3D<GEOSHAPE, MC>::BElementType const &
    RegionMesh3D<GEOSHAPE, MC>:: bElement( ID const & i ) const
    {
        return boundaryFace( i );
    }


    // Templates for interrogate and modify list capacities

    template <typename GEOSHAPE, typename MC>
    template <typename T>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numItems( SimpleVect< T> const & list ) const
    {
        return list.size();
    }

    template <typename GEOSHAPE, typename MC>
    template <typename T>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::maxNumItems( SimpleVect< T> const & list ) const
    {
        return list.capacity();
    }

    template <typename GEOSHAPE, typename MC>
    template <typename T>
    void
    RegionMesh3D<GEOSHAPE, MC>::
    setMaxNumItems( SimpleVect< T> & list, UInt n, std::string title )
    {
        if ( list.capacity() == 0 )
        {
            list.reserve( n );
        }
        else if ( list.capacity() < n )
        {
            Debug(4000) << "WARNING: Resetting " << title << " list size to " << n << "\n";
            Debug(4000) << "         ALL PREVIOUS POINTERS TO THE LIST (IF ANY) ARE NOW INVALID\n";
            list.reserve( n );
        }
    }


    // ***************************** VOLUMES

    template <typename GEOSHAPE, typename MC>
    UInt /*const*/
    RegionMesh3D<GEOSHAPE, MC>::numVolumes() const
    {
        return M_numVolumes;
    }

    template <typename GEOSHAPE, typename MC>
    UInt /*const*/
    RegionMesh3D<GEOSHAPE, MC>::numGlobalVolumes() const
    {
        return M_numGlobalVolumes;
    }

//     template <typename GEOSHAPE, typename MC>
//     UInt &
//     RegionMesh3D<GEOSHAPE, MC>::numVolumes()
//     {
//         return M_numVolumes;
//     }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::maxNumVolumes() const
    {
        return maxNumItems( volumeList );
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::storedVolumes() const
    {
        return numItems( volumeList );
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setMaxNumVolumes( UInt const n, bool const setcounter )
    {
        setMaxNumItems( volumeList, n, "Volume" );
        if ( setcounter )
            M_numVolumes = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumVolumes( UInt const n )
    {
        M_numVolumes = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setMaxNumGlobalVolumes( UInt const n )
    {
        M_numGlobalVolumes = n;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::VolumeType &
    RegionMesh3D<GEOSHAPE, MC>::addVolume()
    {
        return addVolume( VolumeType() );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::VolumeType &
    RegionMesh3D<GEOSHAPE, MC>::addVolume( VolumeType const & v )
    {
        ASSERT_PRE( volumeList.size() < volumeList.capacity() , "Volume list size exceeded" <<
                    volumeList.size() + 1 << " " << volumeList.capacity() ) ;
        volumeList.push_back( v );

        ( volumeList.back() ).setId( volumeList.size() );
        return volumeList.back();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::VolumeType &
    RegionMesh3D<GEOSHAPE, MC>::setVolume( VolumeType const & v, ID const pos )
    {
        ASSERT_PRE( pos <= volumeList.capacity() , "position requested exceed capacity" <<
                    pos << " " << volumeList.capacity() ) ;
        volumeList( pos ) = v;
        volumeList( pos ).setId( pos );
        return volumeList( pos );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    void
    RegionMesh3D<GEOSHAPE, MC>::setVolumeCounter()
    {
        M_numVolumes = volumeList.size();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::VolumeType &
    RegionMesh3D<GEOSHAPE, MC>::lastVolume()
    {
        return volumeList.back();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::VolumeType const &
    RegionMesh3D<GEOSHAPE, MC>::volume( ID const i ) const
    {
        ASSERT_BD( i > 0 && i <= volumeList.size() ) ;
        return volumeList( i );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::VolumeType &
    RegionMesh3D<GEOSHAPE, MC>::volume( ID const i )
    {
        ASSERT_BD( i > 0 && i <= volumeList.size() ) ;
        return volumeList( i );
    }

    // ************************* FACES ******************************
    template <typename GEOSHAPE, typename MC>
    UInt /*const*/
    RegionMesh3D<GEOSHAPE, MC>::numFaces() const
    {
        return M_numFaces;
    }

    template <typename GEOSHAPE, typename MC>
    UInt /*const*/
    RegionMesh3D<GEOSHAPE, MC>::numGlobalFaces() const
    {
        return M_numGlobalFaces;
    }

//     template <typename GEOSHAPE, typename MC>
//     UInt &
//     RegionMesh3D<GEOSHAPE, MC>::numFaces()
//     {
//         return M_numFaces;
//     }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::storedFaces() const
    {
        return numItems( faceList );
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::maxNumFaces() const
    {
        return maxNumItems( faceList );
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setMaxNumFaces( UInt const n, bool const setcounter )
    {
        setMaxNumItems( faceList, n, "Face" );
        if ( setcounter )
            M_numFaces = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setMaxNumGlobalFaces( UInt const n )
    {
        M_numGlobalFaces = n;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::FaceType &
    RegionMesh3D<GEOSHAPE, MC>::addFace( bool const boundary )
    {
        return addFace( FaceType(), boundary );
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::FaceType &
    RegionMesh3D<GEOSHAPE, MC>::addFace( FaceType const & f, bool const /*boundary*/ )
    {
        ASSERT_PRE( faceList.size() < faceList.capacity(), "Face list size exceeded" <<
                    faceList.size() + 1 << " " << faceList.capacity() ) ;
        faceList.push_back( f );
        ( faceList.back() ).setId( faceList.size() );

        return faceList.back();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::FaceType &
    RegionMesh3D<GEOSHAPE, MC>::setFace( FaceType const & f, ID position, bool const boundary )
    {
        ASSERT_PRE( position <= faceList.capacity(), "Face list size exceeded" <<
                    position << " " << faceList.capacity() ) ;
        faceList( position ) = f;
        faceList( position ).setId( position );
#ifdef NOT_BDATA_FIRST

        if ( boundary )
        {
            ASSERT_PRE( position <= _bFaces.capacity(), "Boundary Face list size exceeded" <<
                        _bFaces.size() << " " << bFaces.capacity() ) ;
            _bFaces.push_back( &( faceList( position ) ) );
        }
#endif
        return faceList( position );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::FaceType &
    RegionMesh3D<GEOSHAPE, MC>::lastFace()
    {
        return faceList.back();
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::FaceType const &
    RegionMesh3D<GEOSHAPE, MC>::face( ID const i ) const
    {
        ASSERT_BD( i > 0 && i <= faceList.size() ) ;
        return faceList( i );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::FaceType &
    RegionMesh3D<GEOSHAPE, MC>::face( ID const i )
    {
        ASSERT_BD( i > 0 && i <= faceList.size() ) ;
        return faceList( i );
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::FaceType const &
    RegionMesh3D<GEOSHAPE, MC>::boundaryFace( ID const i ) const
    {
        ASSERT_PRE( faceList.size() != 0, "Boundary Faces not stored" ) ;
        ASSERT_BD( i > 0 && i <= faceList.size() ) ;
        return faceList( i );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::FaceType &
    RegionMesh3D<GEOSHAPE, MC>::boundaryFace( ID const i )
    {
        ASSERT_PRE( faceList.size() != 0, "Boundary Faces not stored" ) ;
        ASSERT_BD( i > 0 && i <= faceList.size() ) ;
        return faceList( i );
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numBFaces() const
    {
        return M_numBFaces;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumFaces( UInt const n )
    {
        M_numFaces = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumBFaces( UInt const n )
    {
        M_numBFaces = n;
    }

    // ************************* EDGES ******************************

    template <typename GEOSHAPE, typename MC>
    UInt /*const*/
    RegionMesh3D<GEOSHAPE, MC>::numEdges() const
    {
        return M_numEdges;
    }

    template <typename GEOSHAPE, typename MC>
    UInt /*const*/
    RegionMesh3D<GEOSHAPE, MC>::numGlobalEdges() const
    {
        return M_numGlobalEdges;
    }

//     template <typename GEOSHAPE, typename MC>
//     UInt &
//     RegionMesh3D<GEOSHAPE, MC>::numEdges()
//     {
//         return M_numEdges;
//     }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::storedEdges() const
    {
        return numItems( edgeList );
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::maxNumEdges() const
    {
        return maxNumItems( edgeList );
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setMaxNumEdges( UInt const n, bool const setcounter )
    {
        setMaxNumItems( edgeList, n, "Edge" );
        if ( setcounter )
            M_numEdges = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setMaxNumGlobalEdges( UInt const n )
    {
        M_numGlobalEdges = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumEdges( UInt const n)
    {
        M_numEdges = n;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::EdgeType &
    RegionMesh3D<GEOSHAPE, MC>::addEdge( bool const boundary )
    {
        return addEdge( EdgeType(), boundary );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::EdgeType &
    RegionMesh3D<GEOSHAPE, MC>::addEdge( EdgeType const & f, bool const /*boundary*/ )
    {
        ASSERT_PRE( edgeList.size() < edgeList.capacity(), "Edge list size exceeded" <<
                    edgeList.size() + 1 << " " << edgeList.capacity() ) ;

        edgeList.push_back( f );
        ( edgeList.back() ).setId(edgeList.size());

        return edgeList.back();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::EdgeType &
    RegionMesh3D<GEOSHAPE, MC>::setEdge( EdgeType const & f, ID position, bool const boundary )
    {
        ASSERT_PRE( position <= edgeList.capacity(), "Edge list size exceeded" <<
                    position << " " << edgeList.capacity() ) ;
        edgeList( position ) = f;
        edgeList( position ).setId( position );
        return edgeList( position );
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::EdgeType &
    RegionMesh3D<GEOSHAPE, MC>::lastEdge()
    {
        return edgeList.back();
    }



    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::EdgeType const &
    RegionMesh3D<GEOSHAPE, MC>::edge( ID const i ) const
    {
        ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
        return edgeList( i );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::EdgeType &
    RegionMesh3D<GEOSHAPE, MC>::edge( ID const i )
    {
        ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
        return edgeList( i );
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::EdgeType const &
    RegionMesh3D<GEOSHAPE, MC>::boundaryEdge( ID const i ) const
    {
        ASSERT_PRE( edgeList.size() != 0, "Boundary Edges not stored" ) ;
        ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
        return edgeList( i );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::EdgeType &
    RegionMesh3D<GEOSHAPE, MC>::boundaryEdge( ID const i )
    {
        ASSERT_PRE( edgeList.size() != 0, "Boundary Edges not stored" ) ;
        ASSERT_BD( i > 0 && i <= edgeList.size() ) ;
        return edgeList( i );
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numBEdges() const
    {
        return M_numBEdges;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumBEdges( ID const n )
    {
        M_numBEdges = n;
    }

    // ************************ Points/Vertices

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numPoints() const
    {
        return M_numPoints;
    }

    template <typename GEOSHAPE, typename MC>
    UInt &
    RegionMesh3D<GEOSHAPE, MC>::numPoints()
    {
        return M_numPoints;
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::storedPoints() const
    {
        return numItems( pointList );
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::storedBPoints() const
    {
        return _bPoints.size();
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::maxNumPoints() const
    {
        return maxNumItems( pointList );
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setMaxNumPoints( UInt const n, bool const setcounter )
    {
        setMaxNumItems( pointList, n, "Point" );
        if ( setcounter )
            M_numPoints = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumVertices( UInt const n )
    {
            M_numVertices = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumGlobalVertices( UInt const n )
    {
            M_numGlobalVertices = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumBVertices( UInt const n )
    {
            M_numBVertices = n;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setMaxNumGlobalPoints( UInt const n )
    {
        M_numGlobalPoints = n;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::PointType &
    RegionMesh3D<GEOSHAPE, MC>::addPoint( bool const boundary, bool const vertex )
    {
        return addPoint( PointType(), boundary, vertex );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::PointType &
    RegionMesh3D<GEOSHAPE, MC>::addPoint( PointType const & p, bool const boundary, bool const /*vertex*/ )
    {
        ASSERT_PRE( pointList.size() < pointList.capacity(), "Point list size exceeded" <<
                    pointList.size() + 1 << " " << pointList.capacity() ) ;

        pointList.push_back( p );

        PointType * pp = & pointList.back();
//        pp->id() = pointList.size();

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
    typename RegionMesh3D<GEOSHAPE, MC>::PointType &
    RegionMesh3D<GEOSHAPE, MC>::setPoint( PointType const & p, ID position, bool const boundary, bool const vertex)
    {
        ASSERT_PRE( position <= pointList.capacity(), "Position  exceed lpoint list capacity" <<
                    position << " " << pointList.capacity() ) ;

        bool found( false );


        pointList.push_back( p );
        PointType * pp = & pointList.back();


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
    typename RegionMesh3D<GEOSHAPE, MC>::PointType &
    RegionMesh3D<GEOSHAPE, MC>::
    setPoint(ID const & position, bool const boundary, bool const vertex)
    {
        ASSERT_PRE( position <= pointList.capacity(), "Position  exceed lpoint list capacity" <<
                    position << " " << pointList.capacity() ) ;
        bool found( false );
        PointType * pp = & pointList( position );
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
    //RegionMesh3D<GEOSHAPE,MC>::PointType &
    UInt
    RegionMesh3D<GEOSHAPE, MC>::addPoint( ID const iden, bool const boundary, UInt const start )
    {
        ASSERT_PRE( pointList.size() <= pointList.capacity(), "Point list size exceeded" <<
                    pointList.size() << " " << pointList.capacity() ) ;
        for ( UInt i = start;i <= pointList.size();++i )
            if ( pointList( i ).id() == iden )
                return i;

        pointList.push_back( PointType() );
        PointType * pp = & pointList.back();
        pp->setId( iden );
        if ( boundary )
        {
            ASSERT_PRE( _bPoints.size() <= _bPoints.capacity(), "Boundary Point list size exceeded" <<
                        _bPoints.size() << " " << _bPoints.capacity() ) ;
            _bPoints.push_back( pp );
            pp->boundary() = true;
        }
        return pointList.size();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::PointType &
    RegionMesh3D<GEOSHAPE, MC>::lastPoint()
    {
        return pointList.back();
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::PointType const &
    RegionMesh3D<GEOSHAPE, MC>::point( ID const i ) const
    {
        ASSERT_BD( i > 0 && i <= pointList.size() ) ;

        return pointList( i );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::PointType &
    RegionMesh3D<GEOSHAPE, MC>::point( ID const i )
    {
        ASSERT_BD( i > 0 && i <= pointList.size() ) ;

        return pointList( i );
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::PointType const &
    RegionMesh3D<GEOSHAPE, MC>::boundaryPoint( ID const i ) const
    {
        ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
        ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
        return *( _bPoints( i ) );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    typename RegionMesh3D<GEOSHAPE, MC>::PointType &
    RegionMesh3D<GEOSHAPE, MC>::boundaryPoint( ID const i )
    {
        ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
        ASSERT_BD( i > 0 && i <= _bPoints.size() ) ;
        return *( _bPoints( i ) );
    }

    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numBPoints() const
    {
        return M_numBPoints;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::setNumBPoints( UInt const n )
    {
        M_numBPoints = n;
        _bPoints.reserve( n );
    }

    template <typename GEOSHAPE, typename MC>
    UInt /*const*/
    RegionMesh3D<GEOSHAPE, MC>::numVertices() const
    {
        return M_numVertices;
    }

    template <typename GEOSHAPE, typename MC>
    UInt /*const*/
    RegionMesh3D<GEOSHAPE, MC>::numGlobalVertices() const
    {
        return M_numGlobalVertices;
    }

//     template <typename GEOSHAPE, typename MC>
//     UInt &
//     RegionMesh3D<GEOSHAPE, MC>::numVertices()
//     {
//         return M_numVertices;
//     }


    template <typename GEOSHAPE, typename MC>
    UInt
    RegionMesh3D<GEOSHAPE, MC>::numBVertices() const
    {
        return M_numBVertices;
    }

    template <typename GEOSHAPE, typename MC>
    UInt &
    RegionMesh3D<GEOSHAPE, MC>::numBVertices()
    {
        return M_numBVertices;
    }

    // *************** GENERAL *******************

    template <typename GEOSHAPE, typename MC>
    std::ostream &
    RegionMesh3D<GEOSHAPE, MC>::showMe( bool verbose, std::ostream & out ) const
    {
        out << "**************************************************" << std::endl;
        out << "**************************************************" << std::endl;
        out << "                      RegionMesh3D                " << std::endl;
        out << "**************************************************" << std::endl;
        out << "**************************************************" << std::endl;
        out << " ID: " << this->id() << " Marker Flag:";
        this->printFlag( out );
        out << std::endl;
        //  out <<"Faces local to  volumes stored: "<<hasLocalFaces()<<std::endl;
        //out <<"Edges local to  volumes stored: "<<hasLocalEdges()<<std::endl;
        //out <<"Edges local to  faces   stored:"<<switches.test("FACEtoEDGE")<<std::endl;
        //out <<"Volumes adjacent to Faces stored: "<<switches.test("FACEtoVOLUME")<<std::endl;
        //out <<"Faces adjacent to Edges stored: "<<switches.test("EDGEtoFACE")<<std::endl<<std::endl;
        //out <<"Internal Faces Stored: " << hasInternalFaces()<<std::endl;
        //out <<"Edges Stored: " << hasEdges()<<" Internal: "
        //    << hasInternalEdges()<<std::endl;
        out << "**************** COUNTERS ********************" << std::endl;
        out << "NumPoints=" << numPoints() << "  " << "numBPoints=" << numBPoints() << std::endl;
        out << "NumVertices=" << numVertices() << "  " << "numBVerices=" << numBVertices() << std::endl;
        out << "NumVolumes=" << numVolumes() << "  " << "numFaces=" << numFaces() << std::endl;
        out << "NumBFaces=" << numBFaces() << "  " << "numEdges=" << numEdges() << std::endl;
        out << "NumBEdges=" << numBEdges() << std::endl;
        out << "**************************************************" << std::endl;
        out << "************ACTUALLY STORED  ********************" << std::endl;
        out << "Points=" << pointList.size() << "  " << "Edges=  " << edgeList.size() << std::endl;
        out << "Faces= " << faceList.size() << "  " << "Volumes=" << volumeList.size() << std::endl;
        out << "**************************************************" << std::endl;
        out << "*******         STATUS OF SWITCHES    ************" << std::endl;
        switches.showMe( verbose, out );
        out << "**************************************************" << std::endl;
        out << "**************************************************" << std::endl;
        if ( verbose )
        {
            std::cout << "Verbose version not implemented yet";
        }
        return out;

    }
    /*! Check mesh

    \param level Indicates level. <ul><li>0 => minimal level, just sets the switches to the appropriate value</li>
    <li>1 => launch an external function from mesh_util.h which performs a series of
    sophisticated tests.
    </ul>
    \param fix If true and level=1 it fixes the mesh
    \param verb Verbosity level
    \param out Strem where output goes (normally std::cerr)

    It returns a severity level. If different from 0 the mesh has problem. If positive the problems are such that the mesh may not
    work in any case. If less than zero it may work in some cases.
    */

    template <class RegionMesh3D>
    extern
    bool checkMesh3D( RegionMesh3D & mesh, Switch & sw, bool fix, bool verbose,
                      std::ostream & out, std::ostream & err, std::ostream & clog );


    template <typename GEOSHAPE, typename MC>
    int
    RegionMesh3D<GEOSHAPE, MC>::check( int level, bool const fix, bool const verb, std::ostream & out )
    {
        int severity = 0;
        Switch testsw;
        if ( verb )
        {
            out << "**************************************************" << std::endl;
            out << "         Checkin  RegionMesh3D                " << std::endl;
            out << " ID: " << this->id() << std::endl;
            out << "**************************************************" << std::endl;
        }
        if ( level == 1 )
        {
            checkMesh3D( *this, testsw, fix, verb, out, out, out );
            if ( verb )
            {
                out << "**********************************************" << std::endl <<
                    "DETAILS OF EXTENDED  CHECK:" << std::endl;
                testsw.showMe( true, out );
                out << "**********************************************" << std::endl << std::endl;
                if ( testsw.test( "ABORT_CONDITION" ) )
                    return 1;

            }
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
                out << "WARNING: No Edges Stored" << std::endl;
            if ( verb )
                out << "MAY NOT WORK IF WE HAVE DOF ON EDGES AND ESSENTIAL BC!!!!" << std::endl;
            severity = -1;
            unsetLinkSwitch( "HAS_ALL_EDGES" );
            unsetLinkSwitch( "HAS_BOUNDARY_EDGES" );
        }
        else if ( edgeList.size() == numBEdges() )
        {

            setLinkSwitch( "HAS_BOUNDARY_EDGES" );
            unsetLinkSwitch( "HAS_ALL_EDGES" );
            if ( verb )
                out << "INFORMATION: Only Boundary Edges Stored" << std::endl;
        }
        else if ( edgeList.size() == numEdges() )
        {
            setLinkSwitch( "HAS_BOUNDARY_EDGES" );
            setLinkSwitch( "HAS_ALL_EDGES" );
        }
        else
        {
            if ( verb )
            {
                out << "ERROR: Inconsistent numbering of edges" << std::endl;
                out << "Counters: Boundary=" << numBEdges() << " Total=" << numEdges() << std::endl;
                out << "Stored =" << edgeList.size() << std::endl;
                severity = 1;
            }

        }

        if ( faceList.size() == 0 )
        {
            if ( verb )
                out << "ERROR: No Faces Stored: at least boundary faces are needed" << std::endl;
            if ( verb )
                out << "" << std::endl;
            severity = 1;
            unsetLinkSwitch( "HAS_ALL_FACES" );
            unsetLinkSwitch( "HAS_BOUNDARY_FACES" );
        }
        else if ( faceList.size() == numBFaces() )
        {
            setLinkSwitch( "HAS_BOUNDARY_FACES" );
            unsetLinkSwitch( "HAS_ALL_FACES" );
            if ( verb )
                out << "INFORMATION: Only Boundary Faces Stored" << std::endl;
        }
        else if ( faceList.size() == numFaces() )
        {
            setLinkSwitch( "HAS_BOUNDARY_FACES" );
            setLinkSwitch( "HAS_ALL_FACES" );
        }
        else
        {
            if ( verb )
            {
                out << "ERROR: Inconsistent numbering of faces" << std::endl;
                out << "Counters: Boundary=" << numBFaces() << " Total=" << numFaces() << std::endl;
                out << "Stored =" << faceList.size() << std::endl;
                severity = 1;
            }
        }


        UInt count = 0;
        for ( typename Points::iterator i = pointList.begin(); i != pointList.end(); ++i )
            if ( i->boundary() )
                ++count;
        if ( count == 0 )
            severity = 4;
        if ( count != M_numBPoints )
        {
            out << "Num Boundary points " << count << " not equal to internal counter value "
                << M_numBPoints << std::endl;
            if ( count != 0 & fix )
            {
                M_numBPoints = count;
                out << "Fixed Counter";
                out.flush();
            }
        }

#if 0
        UInt badid = 0;
        //  for (UInt i=1; i<= storedPoints(); ++i)
        for ( UInt i = 1; i <= M_numVertices; ++i ) // PERICOLOSISSIMO DOMANDARE ALEX
            if ( point( i ).id() != i )
                ++badid;
        if ( badid != 0 )
        {
            severity = 5;
            out << " ERROR:" << badid << "Points ids are wrong";
        }


        badid = 0;
        for ( UInt i = 1; i <= storedEdges(); ++i )
            if ( edge( i ).id() != i )
                ++badid;
        if ( badid != 0 )
        {
            severity = 5;

            out << " SEVERITY ERROR:" << badid << "Edges ids are wrong";
        }


        badid = 0;
        unsigned int adib( 0 );
        unsigned int adisb( 0 );
        unsigned int adi( 0 );
        for ( UInt i = 1; i <= storedFaces(); ++i )
        {

            if ( face( i ).id() != i )
                ++badid;
            if ( i <= numBFaces() )
            {
                if ( face( i ).ad_first() == 0 )
                    ++adib;
                if ( face( i ).ad_second() != 0 )
                    ++adisb;
            }
            else
            {

                if ( face( i ).ad_first() == 0 )
                    ++adi;
                if ( face( i ).ad_second() == 0 )
                    ++adi;
            }
        }
        if ( adib != 0 )
        {
            if ( verb )
            {
                out << "WARNING: " << adib << " Boundary faces have the adjacency Volume ID unset" << std::endl;
                out << " The construction of boundary DOF may be badly affected" << std::endl;
                severity = -2;
            }
        }
        if ( adisb != 0 )
        {
            if ( verb )
            {
                out << "WARNING: " << adib << " Boundary faces have the adjacency Volume ID set as second" << std::endl;
                out << " adiacent element. This is inconsistent!. The construction of boundary DOF may be badly affected" << std::endl;
                severity = -2;
            }
        }
        if ( adi != 0 )
        {
            if ( verb )
            {
                out << "WARNING: " << adib << " internal faces have the adjacency Volume ID unset" << std::endl;
                out << "This in general does not cause any problem" << std::endl;
                severity = -1;
            }
        }
        if ( adib == 0 & adi == 0 )
            setLinkSwitch( "FACES_HAVE_ADIACENCY" );


        if ( badid != 0 )
        {
            severity = 5;
            out << " SEVERITY ERROR:" << badid << "Faces ids are wrong";
        }

        badid = 0;
        for ( UInt i = 1; i <= storedVolumes(); ++i )
            if ( volume( i ).id() != i )
                ++badid;
        if ( badid != 0 )
        {
            severity = 5;
            out << " SEVERITY ERROR:" << badid << "Volumes ids are wrong";
        }
#endif

        if ( M_numVertices == 0 )
        {
            severity = 6;

            out << " SEVERITY ERROR: internal Vertices Counter unset";
        }

        if ( M_numPoints == 0 )
        {
            severity = 6;
            out << " SEVERITY ERROR: internal Points Counter unset";
        }
        if ( M_numPoints == 0 )
        {
            severity = 6;
            out << " SEVERITY ERROR: internal Points Counter unset";
        }
        if ( M_numBPoints == 0 )
        {
            severity = 6;
            out << " SEVERITY ERROR: boundary Points Counter unset";
        }
        if ( M_numBVertices == 0 )
        {
            severity = 6;
            out << " SEVERITY ERROR: boundary Vertices Counter unset";
        }

        if ( verb )
            out << "   Check Finished              " << std::endl <<
                "***********************************************" << std::endl;
        setLinkSwitch( "HAS_BEEN_CHECKED" );

        return severity;

    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::extractEntityList
    (std::vector<ID> & list, ReferenceGeometry const & geometry, EntityFlag const & flag) const
    {
        switch(geometry){
            case VERTEX:
                for (typename Points::const_iterator p=pointList.begin();p!=pointList.end();++p){
                    if (p->hasEqualEntityFlag(flag))list.push_back(p->id());
                }
                break;
            case EDGE:
                for (typename Edges::const_iterator p=edgeList.begin();p!=edgeList.end();++p){
                    if (p->hasEqualEntityFlag(flag))list.push_back(p->id());
                }
                break;
            case FACE:
                for (typename Faces::const_iterator p=faceList.begin();p!=faceList.end();++p){
                    if (p->hasEqualEntityFlag(flag))list.push_back(p->id());
                }
                break;
            case VOLUME:
                for (typename Volumes::const_iterator p=volumeList.begin();p!=volumeList.end();++p){
                    if (p->hasEqualEntityFlag(flag))list.push_back(p->id());
                }
                break;
            default:
                std::cerr<<"Something weird in ExtractEntityList ABORTING"<<endl;
        }
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    ID
    RegionMesh3D<GEOSHAPE, MC>::faceElement( ID const i, UInt const Pos ) const
    {
        ASSERT_PRE( i <= faceList.size(), "Not enough faces stored" ) ;
        ASSERT_PRE( Pos == 1 || Pos == 2 , "Wrong position (1 or 2)" ) ;
        ASSERT_BD( i > 0 ) ;
        if ( Pos == 1 )
        {
            return ( faceList[ i -1 ] ).ad_first();
        }
        else
        {
            return ( faceList[ i -1 ] ).ad_second();
        }
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    ID
    RegionMesh3D<GEOSHAPE, MC>::faceElement( FaceType const & f, UInt const Pos ) const
    {
        ASSERT_BD( ! faceList.empty() ) ;
        ASSERT_PRE( Pos == 1 || Pos == 2 , "Wrong position (1 or 2)" ) ;
        //ASSERT_BD( i >0 ) ;
        if ( Pos == 1 )
        {
            return f.ad_first();
        }
        else
        {
            return f.ad_second();
        }
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    void
    RegionMesh3D<GEOSHAPE, MC>::setLinkSwitch( std::string const & _s )
    {
        ASSERT0( switches.set( _s ), std::stringstream( "Switch named " + _s + " is not allowed" ).str().c_str() );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    void
    RegionMesh3D<GEOSHAPE, MC>::unsetLinkSwitch( std::string const & _s )
    {
        ASSERT0( switches.unset( _s ), std::stringstream( "Switch named " + _s + " is not allowed" ).str().c_str() );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::getLinkSwitch( std::string const & _s ) const
    {
        return switches.test( _s );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::hasInternalFaces() const
    {
        return faceList.size() > M_numBFaces;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::hasInternalEdges() const
    {
        return edgeList.size() > M_numBEdges;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::hasEdges() const
    {
        return ! edgeList.empty();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::hasFaces() const
    {
        return ! faceList.empty();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isVertex( PointType const & p ) const
    {
        return p.id() <= M_numVertices;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isVertex( ID const & id ) const
    {
        return id <= M_numVertices;
    }


    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isBoundaryPoint( ID const & id ) const
    {
        return point( id ).boundary();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isBoundaryPoint( PointType const & p ) const
    {
        return p.boundary();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isBoundaryEdge( EdgeType const & e ) const
    {
        return e.id() <= M_numBEdges;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isBoundaryEdge( ID const & id ) const
    {
        return id <= M_numBEdges;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isBoundaryFace( FaceType const & f ) const
    {
        return f.id() <= M_numBFaces;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isBoundaryFace( ID const & id ) const
    {
        return id <= M_numBFaces;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::hasLocalFaces() const
    {
        return ! _VToF.empty();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::hasLocalEdges() const
    {
        return ! _VToE.empty();
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    ID
    RegionMesh3D<GEOSHAPE, MC>::localFaceId( ID const volId, ID const locF ) const
    {
        ASSERT_PRE( !_VToF.empty(), "Volume to Face array not  set" );
        ASSERT_BD( volId > 0 && volId <= M_numVolumes );
        ASSERT_BD( locF > 0 && locF <= VolumeType::numLocalFaces );
        return _VToF.operator() ( locF, volId );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    UInt
    RegionMesh3D<GEOSHAPE, MC>::localEdgeId( UInt const volId, UInt const locE )
        const
    {
        ASSERT_PRE( !_VToE.empty(), "Volume to Edges array not  set" );
        ASSERT_BD( volId > 0 && volId <= M_numVolumes );
        ASSERT_BD( locE > 0 && locE <= VolumeType::numLocalEdges );
        return _VToE( locE, volId );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    ID
    RegionMesh3D<GEOSHAPE, MC>::localFaceId( const VolumeType & iv, ID const locF ) const
    {
        ASSERT_PRE( !_VToF.empty(), "Volume to Face array not  set" );
        ASSERT_BD( locF > 0 && locF <= VolumeType::numLocalFaces );
        return _VToF( locF, iv.id() );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    UInt
    RegionMesh3D<GEOSHAPE, MC>::localEdgeId( const VolumeType & iv, UInt const locE )
        const
    {
        ASSERT_PRE( !_VToE.empty(), "Volume to Edges array not  set" );
        ASSERT_BD( locE > 0 && locE <= VolumeType::numLocalEdges );
        return _VToE.operator() ( locE, iv.id() );
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isFullEdge( ID const & id ) const
    {
        return edgeList.size() >= id;
    }

    template <typename GEOSHAPE, typename MC>
    INLINE
    bool
    RegionMesh3D<GEOSHAPE, MC>::isFullFace( UInt const & id ) const
    {
        return faceList.size() >= id;
    }


    /********************************************************************************
                         ELEMENT3D:GLOBAL FACES/EDGES
    *******************************************************************************/
    // Forward Declarations
    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::updateElementEdges( bool ce, bool verbose, UInt ee )
    {
        // If the counter is set we trust it! Otherwise we use Euler formula
        // this is ok for domains with at most 1 hole!

        if (verbose)
            std::cout << "     Updating element edges ... " << std::flush;

        if ( ce && ee == 0 )
            ee = M_numEdges > M_numBEdges ? M_numEdges : ( GEOSHAPE::numFaces / 2 - 1 ) * numVolumes() + M_numBFaces / 2 + numVertices();


        if ( ce )
        {
            // We want to create the edges, yet we need to clear existing edges, since we start from scratch!
            edgeList.reserve( ee );
            edgeList.resize( 0 );
        }
        BareItemsHandler<BareEdge> _be;
        pair<UInt, bool> e;
        _VToE.reshape( numLocalEdges(), numVolumes() ); // DIMENSION ARRAY

        UInt vid, i1, i2;
        pair<BareEdge, bool> _edge;
        GEOSHAPE ele;
        FaceShape bele;
        // First We check if we have already Edges stored
        if ( ! edgeList.empty() )
        {
            // dump first the existing edges, to maintain the correct numbering
            // if everything is correct the numbering in the bareedge
            // structure will reflect the actual edge numbering
            pair<UInt, bool> _check;
            for ( UInt j = 0; j < edgeList.size();++j )
            {
                i1 = ( edgeList[ j ].point( 1 ) ).id();
                i2 = ( edgeList[ j ].point( 2 ) ).id();

                _edge  = makeBareEdge( i1, i2 );
                _check = _be.addIfNotThere( _edge.first );
            }
        }

        EdgeType edg;

        if ( edgeList.empty() )
        {
            // We want that the first edges be those on the boundary, in order to obay the paradigm for
            // a RegionMesh3D

            for ( typename Faces::iterator ifa = faceList.begin();
                  ifa != faceList.begin() + M_numBFaces; ++ifa )
            {


                for ( UInt j = 1; j <= numLocalEdgesOfFace(); j++ )
                {
                    i1 = bele.eToP( j, 1 );
                    i2 = bele.eToP( j, 2 );
                    // go to global
                    i1 = ( ifa->point( i1 ) ).id();
                    i2 = ( ifa->point( i2 ) ).id();

                    _edge = makeBareEdge( i1, i2 );

                    e = _be.addIfNotThere( _edge.first );

                    if ( ce && e.second )
                    {
                        for ( UInt k = 1; k <= 2 + FaceShape::nbPtsPerEdge; k++ )
                        {
                            UInt inode = bele.eToP(j, k);
                            edg.setPoint( k, ifa->point( inode ) );
                        }
                        inheritWeakerMarker( edg );
                        addEdge( edg, true );
                    }
                }

            }

        }

        if ( ce )
        {
            M_numBEdges = edgeList.size();
            setLinkSwitch( "HAS_BOUNDARY_EDGES" );
        }

        for ( typename Volumes::iterator iv = volumeList.begin();
              iv != volumeList.end(); ++iv )
        {
            vid = iv->localId();
            // REMEMBER: numbering from 1
            for ( UInt j = 1;j <= numLocalEdges();j++ )
            {
                i1 = ele.eToP( j, 1 );
                i2 = ele.eToP( j, 2 );
                // go to global
                i1 = ( iv->point( i1 ) ).id();
                i2 = ( iv->point( i2 ) ).id();
                _edge = makeBareEdge( i1, i2 );

                e = _be.addIfNotThere( _edge.first );
                _VToE.operator() ( j, vid ) = e.first;
                if ( ce && e.second )
                {
                    for ( UInt k = 1; k <= 2 + VolumeShape::nbPtsPerEdge;k++ )
                        {
                            UInt inode = ele.eToP(j, k);
                            edg.setPoint( k, iv->point( inode ) );
                        }
                    inheritWeakerMarker( edg );
                    addEdge( edg, false );
                }
            }
        }

        if ( ce )
        {
            M_numEdges = edgeList.size();
            setLinkSwitch( "HAS_ALL_EDGES" );
        }

        UInt n = _be.maxId();

        if (!ce)
            {
                if ( M_numEdges == 0 || M_numEdges == M_numBEdges )
                    M_numEdges = n;
            }
        else
            {
                    M_numGlobalEdges = n;
            }

        if (verbose)
            std::cout << n << " edges ";
        ASSERT_POS( n == M_numEdges , "#Edges found is not equal to that in RegionMesh" << n << " " << M_numEdges ) ;
        setLinkSwitch( std::string( "HAS_VOLUME_TO_EDGES" ) );

        if (verbose)
            std::cout << " done." << std::endl;
    }


//
// Update element faces
//

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::updateElementFaces( bool cf, const bool verbose, UInt ef )
    {

        if (verbose)
            std::cout << "     Updating element faces ... " << std::flush;

        ASSERT0( ! cf || M_numBFaces > 0, std::stringstream( std::string("Boundary Faces Must have been set") +
                                                             std::string("in order to call updateElementFaces with createFaces=true") +
                                                             std::string("\nUse buildBoundaryFaces(..) from mesh_util.h") ).str().c_str() );
        // If the counter is set we trust it! Otherwise we use Euler formula

        if ( cf && ef == 0 )
            ef = M_numFaces > M_numBFaces ? M_numFaces : ( GEOSHAPE::numFaces * numVolumes() + M_numBFaces ) / 2;

        ASSERT( cf || numFaces() > 0 , "Mesh is not properly set!" );

        if ( cf )
            faceList.reserve( ef );

        FaceType face;



        BareItemsHandler<BareFace> _be;
        pair<UInt, bool> e;
        _VToF.reshape( numLocalFaces(), numVolumes() ); // DIMENSION ARRAY

        UInt vid, i1, i2, i3, i4;
        pair<BareFace, bool>_face;
        GEOSHAPE ele;
        // If we have all faces and the faces store all adjacency info
        // everything is easier
        if ( (faceList.size() == numFaces()) & getLinkSwitch( "FACES_HAVE_ADIACENCY" ) & getLinkSwitch( "HAS_ALL_FACES" ) )
        {
            for ( typename Faces::iterator itf = faceList.begin();itf != faceList.end();++itf )
            {
                if ( itf->pos_first() != 0 )
                    _VToF( itf->pos_first() , itf->ad_first() ) = itf->localId();
                if ( itf->pos_second() != 0 )
                    _VToF( itf->pos_second(), itf->ad_second() ) = itf->localId();
            }
            // we finish here
            setLinkSwitch( "HAS_VOLUME_TO_FACES" );
            if (verbose)
                std::cout << " done." << std::endl;

            return ;
        }

        // If I have only boundary faces I need to process them first to keep the correct numbering

        // First We check if we have already Faces stored
        if ( ! faceList.empty() )
        {
            // dump all faces in the container, to maintain the correct numbering
            // if everything is correct the numbering in the bareface structure
            // will reflect the actual face numbering However, if I want to create
            // the internal faces I need to make sure that I am processing only the
            // boundary ones. So I resize the container!
            if ( cf )
                faceList.resize( M_numBFaces );

            pair<UInt, bool> _check;
            for ( UInt j = 0; j < faceList.size();++j )
            {
                i1 = ( faceList[ j ].point( 1 ) ).localId();
                i2 = ( faceList[ j ].point( 2 ) ).localId();
                i3 = ( faceList[ j ].point( 3 ) ).localId();
                if ( FaceShape::numVertices == 4 )
                {
                    i4 = ( faceList[ j ].point( 4 ) ).localId();
                    _face = makeBareFace( i1, i2, i3, i4 );
                }
                else
                {
                    _face = makeBareFace( i1, i2, i3 );
                }
                _check = _be.addIfNotThere( _face.first );
            }
        }

        for ( typename Volumes::iterator iv = volumeList.begin();
              iv != volumeList.end(); ++iv )
        {
            vid = iv->localId();
            for ( UInt j = 1;j <= numLocalFaces();j++ )
            {
                i1 = ele.fToP( j, 1 );
                i2 = ele.fToP( j, 2 );
                i3 = ele.fToP( j, 3 );
                i1 = ( iv->point( i1 ) ).localId();
                i2 = ( iv->point( i2 ) ).localId();
                i3 = ( iv->point( i3 ) ).localId();
                if ( FaceShape::numVertices == 4 )
                {
                    i4 = ele.fToP( j, 4 );
                    i4 = ( iv->point( i4 ) ).localId();
                    _face = makeBareFace( i1, i2, i3, i4 );
                }
                else
                {
                    _face = makeBareFace( i1, i2, i3 );
                }

                e = _be.addIfNotThere( _face.first );
                _VToF( j, vid ) = e.first;
                if ( cf )
                {
                    if ( e.second )
                    {
                        // a new face It must be internal.
                        for ( UInt k = 1; k <= FaceType::numPoints; ++k )
                            face.setPoint( k, iv->point( ele.fToP( j, k ) ) );

                        face.ad_first()  = vid;
                        face.pos_first() = j;

                        // gets the marker from the RegionMesh

                        face.setMarker( this->marker() );
                        addFace( face, false ); //The id should be correct
                    }
                    else
                    {
                        // We assume that BFaces have been already set so we have to do
                        // nothing if the face is on the boundary
                        if ( e.first > M_numBFaces )
                        {
                            faceList( e.first ).ad_second() = vid;
                            faceList( e.first ).pos_second() = j;

                        }
                        else
                            {
                            }
                    }
                }
            }
        }

        UInt n = _be.maxId();
        // Fix _numfaces if it was not set or set to just the # of Bfaces
        if (!cf)
            {
                if ( M_numFaces == 0 || M_numFaces == M_numBFaces )
                    M_numFaces = n;
            }
        else
            {
                    M_numGlobalFaces = n;
            }


        if (verbose) std::cout << n << " faces ";
        ASSERT_POS( n == M_numFaces , "#Faces found inconsistent with that stored in RegionMesh" ) ;
        setLinkSwitch( "HAS_VOLUME_TO_FACES" );
        if ( cf )
            setLinkSwitch( "HAS_ALL_FACES" );
        if ( cf )
            setLinkSwitch( "FACES_HAVE_ADIACENCY" );

        if (verbose)
            std::cout << " done." << std::endl;
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::cleanElementFaces()
    {
        _VToF.clean();

        unsetLinkSwitch( "HAS_VOLUME_TO_FACES" );
    }

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::cleanElementEdges()
    {
        _VToE.clean();
        unsetLinkSwitch( "HAS_VOLUME_TO_EDGES" );
    }


    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::
    list_of_points( bool ( *fct ) ( double, double, double ), std::vector<UInt>& list_pts )
    {
        for ( UInt i = 1;i <= M_numPoints;i++ )
        {
            Geo0D& pt = pointList( i );
            if ( fct( pt.x(), pt.y(), pt.z() ) )
            {
                list_pts.push_back( i );
            }
        }
    }
#ifdef HAVE_MPI
template <typename GEOSHAPE, typename MC>
void
RegionMesh3D<GEOSHAPE, MC>::
orderMesh(MPI_Comm comm) // serial reordering:
    //a fake communicator must be passed to the function if we intend to use the serial reordering.
    // For a parallel reordering we should pass a partitioned mesh and the real communicator (to be implemented).
{

    std::vector<int> vertexDist(2);
    std::vector<int> adjncy;
    std::vector<int> xadj/*(this->numVolumes()+1)*/;//dimension # elements * # nodes per element
    std::vector<int> order(this->M_numPoints);//(this->numVolumes());
    std::vector<int> size(2);

    UInt elementNodes;
    UInt elementFaces;

    typedef VolumeShape ElementShape;

    switch ( ElementShape::Shape )
    {
        case HEXA:
            elementNodes = 8;
            break;
        case TETRA:
            elementNodes = 4;
            break;
        default:
            ERROR_MSG( "Element shape not implement in partitionMesh" );
    }


    switch ( FaceShape::Shape )
    {
        case QUAD:
            elementFaces = 6;
            break;
        case TRIANGLE:
            elementFaces = 4;
            break;
        default:
            ERROR_MSG( "Face Shape not implemented in partitionMesh" );
    }


    adjncy.resize(0);
    xadj.resize(0);
    xadj.push_back(0);

    UInt sum = 0;

    bool flag=false;
    for ( UInt iNode = 0; iNode < M_numPoints; ++iNode )
    {
        flag=false;
                for( UInt ed = 0; ed < edgeList.size(); ++ed )
        {
                    if(pointList[iNode].id()==edgeList[ed].point(1).id())
                      {
                          adjncy.push_back(edgeList[ed].point(2).id()-1);
                          ++sum;
                          flag=true;
                      }
                    else if(pointList[iNode].id()==edgeList[ed].point(2).id())
                        {
                            adjncy.push_back(edgeList[ed].point(1).id()-1);
                            ++sum;
                            flag=true;

                        }
        }
                ASSERT(flag==true,"vertex not inserted");
                xadj.push_back(sum);
    }


    //std::vector<int> adjncy2(adjncy);
    //std::vector<int> xadj2(xadj)/*(this->numVolumes()+1)*/;//dimension # elements * # nodes per element

     int numflag=0;
     vertexDist[0]=0;
     //     vertexDist[1]=this->numVolumes()*elementNodes;
     vertexDist[1]=M_numPoints;

     ParMETIS_V3_NodeND((int*) &vertexDist[0],
                        (int*) &xadj[0],
                        (int*) &adjncy[0],
                        &numflag,
                        (int*) new int(0),//options[0],
                        &order[0],// output vector, size: nb of local nodes (equal to the vector part).
                        //&edgecut, &part[0] // output of ParMETIS_V3_PartKway
                        &size[0],// output vector, size: 2*nb of processors.
                        &comm );//dynamic_cast<Epetra_MpiComm*>(&this->Comm())->Comm());


    std::vector<Real> ics(pointList.size());
    std::vector<Real> ipsilon(pointList.size());
    std::vector<Real> zeta(pointList.size());
    std::vector<ID> mk(pointList.size());

    for ( UInt iv = 0; iv < pointList.size(); ++iv )
        {
            ics[order[iv]]=pointList[iv].x();
            ipsilon[order[iv]]=pointList[iv].y();
            zeta[order[iv]]=pointList[iv].z();
            mk[order[iv]]=EntityFlag( pointList[iv].marker() );
        }

    for ( UInt iv = 0; iv < pointList.size(); ++iv )
         {
             pointList[iv].x()=ics[iv];
             pointList[iv].y()=ipsilon[iv];
             pointList[iv].z()=zeta[iv];
             pointList[iv].setMarker( mk[iv] );

             pointList[iv].setId(order[iv]+1);
             pointList[iv].setLocalId(order[iv]+1);
        }
    /*        for(UInt i=0; i<numVertices(); ++i)
            {
               std::cout<<"pointList "<<i<<"= "<<pointList[i].id()<<" oldPointList "<<i<<"= "<<newList[i].id()<<std::endl;
               }*/

    }
template <typename GEOSHAPE, typename MC>
void
RegionMesh3D<GEOSHAPE, MC>::
edgeMarkers(std::map<ID, ID> const& locDof, UInt newMarker)
{
    std::map<ID, ID>::const_iterator IT;
    for(IT=locDof.begin(); IT!=locDof.end(); ++IT)
        {
            pointList(IT->second).setMarker(newMarker);
        }
}

#endif

    template <typename GEOSHAPE, typename MC>
    void
    RegionMesh3D<GEOSHAPE, MC>::
    printLtGMap(std::ostream & os)
    {
    	std::map<int,int>::iterator iter;

    	os << "[RegionMesh3D] Local to Global Map" << std::endl;
    	os << "Number of Local Points\t=\t" << M_numPoints << std::endl;
    	os << "Local ID\t\t\tGlobal ID" << std::endl;

    	for( iter = M_localToGlobalNode.begin(); iter != M_localToGlobalNode.end(); ++iter )
    	{
    	    os << iter->first << "\t\t\t" << iter->second << std::endl;
    	}
    }

}
#endif
