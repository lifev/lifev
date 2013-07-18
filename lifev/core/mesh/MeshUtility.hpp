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
    @brief Base utilities operating on meshes

    @contributor Tiziano Passerini <tiziano@mathcs.emory.edu>
    @maintainer Tiziano Passerini <tiziano@mathcs.emory.edu>

    This file contains a set of base utilities used to test mesh entities or
    operate on them
 */

#ifndef MESHUTILITY_H
#define MESHUTILITY_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <algorithm>
#include <iterator>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/Switch.hpp>
#include <lifev/core/mesh/MeshElementBare.hpp>
#include <lifev/core/mesh/MarkerDefinitions.hpp>
#include <lifev/core/mesh/MeshEntity.hpp>
#include <lifev/core/mesh/MeshEntityContainer.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/VectorSmall.hpp>

namespace LifeV
{

namespace MeshUtility
{
/*
   todo change the class names:
       "EnquireBEntity --> EnquireBoundaryEntity" and so on
       "GetCoordComponent --> GetCoordinatesComponent" and so on
 */


//! A locally used structure, not meant for general use
typedef std::map < BareFace, std::pair<ID, ID >,
        cmpBareItem<BareFace> > temporaryFaceContainer_Type;

//! A locally used structure, not meant for general use
typedef std::map < BareEdge, std::pair<ID, ID>,
        cmpBareItem<BareEdge> > temporaryEdgeContainer_Type;

/*
 *******************************************************************************
                            FUNCTORS
 *******************************************************************************
 */
//! @defgroup Predicates Some useful functors to be used to test mesh entities

//! Functor to check if a Face is on the boundary
/*!
    @ingroup Predicates

    This object uses the information contained in a FaceContainer produced
    by findBoundaryFaces(). It does not use the information contained in the
    mesh PointList, so it may be used to set up a proper mesh

    @pre boundaryFaceContainer has been previously set by a call to findBoundaryFaces()
 */
template <typename MeshType>
class EnquireBFace
{
public:

    //! @name Public Types
    //@{
    typedef MeshType                            mesh_Type;
    typedef mesh_Type const*                    meshPtr_Type;
    typedef typename mesh_Type::face_Type        face_Type;
    typedef typename mesh_Type::facetShape_Type       faceShape_Type;
    typedef temporaryFaceContainer_Type const* temporaryFaceContainerPtr_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructor taking a mesh object and a face Entity
    /*!
        @param boundaryFaceContainer a container of boundary faces
        @pre boundaryFaceContainer has been previously set by a call to findBoundaryFaces()
     */
    EnquireBFace ( temporaryFaceContainer_Type const& boundaryFaceContainer ) :
        boundaryFaceContainerPtr ( &boundaryFaceContainer )
    {}

    //! Virtual Destructor
    virtual ~EnquireBFace()
    {}
    //@}

    //! @name Operators
    //@{

    //! The function call operator
    /*!
        @param face a face entity in the mesh_Type
        @return true if the face is on the boundary, false otherwise
     */
    bool operator() ( const face_Type& face ) const
    {
        ID point1Id, point2Id, point3Id, point4Id;
        BareFace bareFace;

        point1Id = face.point ( 0 ).localId();
        point2Id = face.point ( 1 ).localId();
        point3Id = face.point ( 2 ).localId();
        if ( faceShape_Type::S_numVertices == 4 )
        {
            point4Id = face.point ( 3 ).localId();
            bareFace = ( makeBareFace ( point1Id, point2Id, point3Id, point4Id ) ).first;
        }
        else
        {
            bareFace = ( makeBareFace ( point1Id, point2Id, point3Id ) ).first;
        }
        return ( boundaryFaceContainerPtr->find ( bareFace ) != boundaryFaceContainerPtr->end() );
    }

    //@}
private:
    //! @name Private Methods
    //@{

    //! Empty Constructor
    EnquireBFace()
    {}
    //@}
    temporaryFaceContainerPtr_Type boundaryFaceContainerPtr;
};


/*! Functor to check if an edge is on the boundary

    @ingroup Predicates

    This object uses the information contained in an EdgeContainer produced
    by findBoundaryEdges(). It does not use the information contained in the mesh
    PointList, so it differs from EnquireBEntity.

    @pre boundaryEdgeContainer has been previously set by a call to findBoundaryEdges()

 */
template <typename MeshType>
class EnquireBEdge
{
public:

    //! @name Public Types
    //@{
    typedef MeshType                              mesh_Type;
    typedef mesh_Type const*                      meshPtr_Type;
    typedef typename mesh_Type::edge_Type         edge_Type;
    typedef temporaryEdgeContainer_Type const*    temporaryEdgeContainerPtr_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructor taking a mesh object and an edge container
    /*!
        @param mesh a mesh object
        @param boundaryEdgeContainer a container of boundary edges
     */
    EnquireBEdge ( temporaryEdgeContainer_Type const& boundaryEdgeContainer ) :
        boundaryEdgeContainerPtr ( &boundaryEdgeContainer )
    {}

    //! Copy Constructor
    EnquireBEdge ( EnquireBEdge const& enquireBoundaryEdge ) :
        boundaryEdgeContainerPtr ( enquireBoundaryEdge.boundaryEdgeContainerPtr )
    {}

    //! Virtual Destructor
    virtual ~EnquireBEdge()
    {}
    //@}

    //! @name Operators
    //@{

    //! The function call operator
    /*!
        @param edge an edge entity in the mesh_Type
        @return true if the edge is on the boundary, false otherwise
     */
    bool operator() ( const edge_Type& edge ) const
    {
        ID point1Id, point2Id;
        BareEdge bareEdge;

        point1Id = edge.point ( 0 ).localId();
        point2Id = edge.point ( 1 ).localId();
        bareEdge = ( makeBareEdge ( point1Id, point2Id ) ).first;
        return boundaryEdgeContainerPtr->find ( bareEdge ) != boundaryEdgeContainerPtr->end();
    }

    //@}
private:
    //! @name Private Methods
    //@{

    //! Empty Constructor
    EnquireBEdge()
    {}
    //@}

    temporaryEdgeContainerPtr_Type boundaryEdgeContainerPtr;
};


/*! Functor to check if a mesh entity with boundary indicator is on the boundary

    @ingroup Predicates

    This objects works on mesh entities with boundary indicator (for instance a GeoPoint)
    by enquiring its boundary flag. To be used only for checks.

    @warning It assumes that boundary indicators are correctly set.
 */
template <typename MeshType>
class EnquireBPoint
{
public:

    //! @name Public Types
    //@{
    typedef MeshType                      mesh_Type;
    typedef mesh_Type*                    meshPtr_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructor taking a mesh object
    /*!
        @param mesh a mesh object
     */
    explicit EnquireBPoint ( mesh_Type& mesh ) : meshPtr ( &mesh )
    {}

    //! Copy Constructor
    EnquireBPoint ( EnquireBPoint const& enquireBoundaryPoint ) :
        meshPtr ( enquireBoundaryPoint.meshPtr )
    {}

    //! Virtual Destructor
    virtual ~EnquireBPoint()
    {}
    //@}

    //! @name Operators
    //@{

    //! The function call operator
    /*!
        @param MeshEntity a mesh entity with boundary indicator
        @return true if the entity is on the boundary, false otherwise
     */
    bool operator() ( const MeshEntity& meshEntity ) const
    {
        return meshEntity.boundary();
    }
    //@}

private:
    //! @name Public Types
    //@{

    //! Empty Constructor
    EnquireBPoint()
    {}
    //@}

    meshPtr_Type meshPtr;
};


/*! @ingroup Test_Functors

    @brief This functor is used to do some geometry checks.

    It is instantiated with the index of the component to be extracted from a vector
    of geometric coordinates.
 */
class GetCoordComponent
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    GetCoordComponent();

    //! Constructor taking the component
    /*!
        @param i the index of the component to be extracted from a vector of geometric coordinates
     */
    explicit GetCoordComponent ( Int i );

    //! Copy constructor
    GetCoordComponent ( const GetCoordComponent& getCoordComponent ) :
        componentIndex ( getCoordComponent.componentIndex )
    {}

    //! Virtual Destructor
    virtual ~GetCoordComponent()
    {}
    //@}


    //! @name Operators
    //@{

    //! The function call operator
    /*!
        If componentIndex is valid (0 <= 2), the function call operator will operate on a vector
        by returning a new vector of geometric coordinates, with the single
        non null component corresponding to the componentIndex.

        Otherwise, the vector given as an argument will be returned unchanged.

        @param[in] x first component of the input vector
        @param[in] y second component of the input vector
        @param[in] z third component of the input vector
        @param[out] ret output vector
     */
    void operator() ( Real const& x, Real const& y,
                      Real const& z, Real ret[ 3 ] ) const;

    //@}
private:
    Int componentIndex;
};


/*! @ingroup Test_Functors

    This functor is used to do some geometry checks.
    It returns a vector of ones
 */
class GetOnes
{
public:
    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    GetOnes()
    {}

    //! Copy Constructor
    GetOnes ( const GetOnes& /* getOnes */ )
    {}

    //! Virtual Destructor
    virtual ~GetOnes()
    {}
    //@}


    //! @name Operators
    //@{

    //! The function call operator
    /*!
        @param[in] x first component of the input vector
        @param[in] y second component of the input vector
        @param[in] z third component of the input vector
        @param[out] ret output vector (always a vector of ones)
     */
    void operator() ( Real const& x, Real const& y,
                      Real const& z, Real ret[ 3 ] ) const;
    //@}
};


/*
 *******************************************************************************
                            EDGES/FACES FINDERS
 *******************************************************************************
 */

//! Finds mesh faces
/*!
    A low level routine, not meant to be called directly. It creates a
    container with all the information needed to set up properly the boundary
    faces connectivities.

    @param mesh A 3D mesh.

    @param boundaryFaceContainer[out] This container will eventually contain a map whose key are
    the BareFaces corresponding to the boundary faces; the map data are pairs of
    IDs: the ID of the adjacent element and the relative position of the face
    in the element. The search of boundary faces does not rely on the proper setting of boundary nodes.

    @param NumInternalFaces[out] A reference to an integer returning the number of internal faces found.

    @param[out] internalFaces A container that will possibly contain a map whose keys are
    the BareFaces corresponding to internal faces and whose data are pairs of IDs:
    the ID of the two elements adjacent to the face.

    @param buildAllFaces When this bool is set true the function will also construct the set
    of internal faces, stored in internalFaces.

    @return Number of boundary faces found

    @note this method is intended to work on 3D meshes
 */
template <typename MeshType>
UInt findFaces ( const MeshType& mesh, temporaryFaceContainer_Type& boundaryFaceContainer,
                 UInt& numInternalFaces, temporaryFaceContainer_Type& internalFaces,
                 bool buildAllFaces = false )
{
    UInt                                  point1Id, point2Id, point3Id, point4Id;
    BareFace                              bareFace;
    typename MeshType::elementShape_Type   volumeShape;
    typedef typename MeshType::volumes_Type    volumeContainer_Type;
    temporaryFaceContainer_Type::iterator faceContainerIterator;

    // clean first in case it has been already used
    boundaryFaceContainer.clear();
    if ( buildAllFaces )
    {
        internalFaces.clear();
    }
    numInternalFaces = 0;

    for ( typename volumeContainer_Type::const_iterator volumeContainerIterator = mesh.volumeList.begin();
            volumeContainerIterator != mesh.volumeList.end(); ++volumeContainerIterator )
    {
        for ( ID jFaceLocalId = 0; jFaceLocalId < mesh.numLocalFaces(); ++jFaceLocalId )
        {
            point1Id = volumeShape.faceToPoint ( jFaceLocalId, 0 );
            point2Id = volumeShape.faceToPoint ( jFaceLocalId, 1 );
            point3Id = volumeShape.faceToPoint ( jFaceLocalId, 2 );
            // go to global
            point1Id = ( volumeContainerIterator->point ( point1Id ) ).localId();
            point2Id = ( volumeContainerIterator->point ( point2Id ) ).localId();
            point3Id = ( volumeContainerIterator->point ( point3Id ) ).localId();
            if ( MeshType::facetShape_Type::S_numVertices == 4 )
            {
                point4Id = volumeShape.faceToPoint ( jFaceLocalId, 3 );
                point4Id = ( volumeContainerIterator->point ( point4Id ) ).localId();
                bareFace = ( makeBareFace ( point1Id, point2Id, point3Id, point4Id ) ).first;
            }
            else
            {
                bareFace = ( makeBareFace ( point1Id, point2Id, point3Id ) ).first;
            }

            if ( ( faceContainerIterator = boundaryFaceContainer.find ( bareFace ) ) == boundaryFaceContainer.end() )
            {
                boundaryFaceContainer.insert (
                    std::make_pair ( bareFace, std::make_pair ( volumeContainerIterator->localId(), jFaceLocalId ) ) );
            }
            else
            {
                if ( buildAllFaces && point1Id > point2Id )
                {
                    internalFaces.insert (
                        ( std::make_pair ( bareFace, std::make_pair ( volumeContainerIterator->localId(), jFaceLocalId ) ) ) );
                }
                boundaryFaceContainer.erase ( faceContainerIterator ); // counted twice: internal face
                ++numInternalFaces;
            }
        }
    }
    return boundaryFaceContainer.size();
}


//! Finds boundary faces
/*!
    A low level routine, not meant to be called directly. It creates a
    container with all the information needed to set up properly the boundary
    face connectivities.

    @param mesh A 3D mesh.

    @param boundaryFaceContainer[out] This container will eventually contain a map whose key are
    the BareFaces corresponding to the boundary faces; the map data are pairs of
    IDs: the ID of the adjacent element and the relative position of the face
    in the element. The search of boundary faces does not rely on the proper setting of boundary nodes.

    @param numInternalFaces[out] A reference to an integer returning the number of internal faces found.

    @return Number of boundary faces found.

    @note this method is intended to work on 3D meshes
 */
template <typename MeshType>
UInt findBoundaryFaces ( const MeshType& mesh,
                         temporaryFaceContainer_Type& boundaryFaceContainer,
                         UInt& numInternalFaces )
{
    temporaryFaceContainer_Type dummy;
    return findFaces ( mesh, boundaryFaceContainer, numInternalFaces, dummy, false );
}


//! Finds boundary edges
/*!
    A low level routine, not meant to be called directly. It creates a
    container with all the information needed to set up properly the boundary
    edges connectivities.

    @param mesh A mesh.

    @param boundaryEdgeContainer[out] This container will eventually contain a map whose key are
    the BareEdges corresponding to the boundary edges; the map data are pairs of
    IDs: the ID of the adjacent face and the relative position of the edge
    in the element. The search of boundary edges does not rely on the proper setting of boundary nodes.

    @return Number of boundary edges found.

    @pre The list of boundary faces in the mesh must be correctly set.
 */
template <typename MeshType>
UInt findBoundaryEdges ( const MeshType& mesh, temporaryEdgeContainer_Type& boundaryEdgeContainer )
{
    UInt                                 point1Id, point2Id;
    BareEdge                             bareEdge;
    typedef typename MeshType::facetShape_Type facetShape_Type;
    typedef typename MeshType::faces_Type     faceContainer_Type;


    if ( ! mesh.hasFaces() )
    {
        return 0;
    }

    // clean first in case it has been already used
    boundaryEdgeContainer.clear();

    // the following cycle assumes to visit only the boundary faces in mesh.faceList()
    for ( typename faceContainer_Type::const_iterator faceContainerIterator = mesh.faceList.begin();
            faceContainerIterator != mesh.faceList.begin() + mesh.numBFaces(); ++faceContainerIterator )
    {
        for ( ID jEdgeLocalId = 0; jEdgeLocalId < mesh.numLocalEdgesOfFace(); ++jEdgeLocalId )
        {
            point1Id = facetShape_Type::edgeToPoint ( jEdgeLocalId, 0 );
            point2Id = facetShape_Type::edgeToPoint ( jEdgeLocalId, 1 );
            // go to global
            point1Id = ( faceContainerIterator->point ( point1Id ) ).localId();
            point2Id = ( faceContainerIterator->point ( point2Id ) ).localId();
            bareEdge = ( makeBareEdge ( point1Id, point2Id ) ).first;
            boundaryEdgeContainer.insert (
                std::make_pair ( bareEdge, std::make_pair ( faceContainerIterator->localId(), jEdgeLocalId ) ) );
        }
    }
    return boundaryEdgeContainer.size();
}


//! Finds internal edges.
/*!
    A low level routine, not meant to be called directly. It creates a
    container with all the information needed to set up properly the edge
    connectivities.

    @param mesh A 3D mesh.

    @param boundaryEdgeContainer[in] This container contains a map whose key are
    the BareEdges corresponding to the boundary edges; the map data are pairs of
    IDs (if generated from findBoundaryEdges(), those pairs will contain be the ID
    of the adjacent mesh element and the relative position of the edge in the element).

    @param internalEdgeContainer[out] This container contains a map whose key are
    the BareEdges corresponding to the boundary edges; the map data are pairs of
    IDs: the ID of the adjacent element and the relative position of the edge
    in the element.

    @return Number of edges found.

 */
template <typename MeshType>
UInt findInternalEdges ( const MeshType& mesh,
                         const temporaryEdgeContainer_Type& boundaryEdgeContainer,
                         temporaryEdgeContainer_Type& internalEdgeContainer )
{
    UInt                                   point1Id, point2Id;
    BareEdge                               bareEdge;
    typedef typename MeshType::elementShape_Type volumeShape_Type;
    typedef typename MeshType::volumes_Type     volumeContainer_Type;
    temporaryEdgeContainer_Type            temporaryEdgeContainer;

    ASSERT0 ( mesh.numVolumes() > 0, "We must have some 3D elements stored n the mesh to use this function!" );

    internalEdgeContainer.clear();
    internalEdgeContainer.swap (temporaryEdgeContainer);

    for ( typename volumeContainer_Type::const_iterator volumeContainerIterator = mesh.volumeList.begin();
            volumeContainerIterator != mesh.volumeList.end(); ++volumeContainerIterator )
    {
        for ( ID jEdgeLocalId = 0; jEdgeLocalId < mesh.numLocalEdges(); ++jEdgeLocalId )
        {
            point1Id = volumeShape_Type::edgeToPoint ( jEdgeLocalId, 0 );
            point2Id = volumeShape_Type::edgeToPoint ( jEdgeLocalId, 1 );
            // go to global
            point1Id = ( volumeContainerIterator->point ( point1Id ) ).localId();
            point2Id = ( volumeContainerIterator->point ( point2Id ) ).localId();
            bareEdge = ( makeBareEdge ( point1Id, point2Id ) ).first;
            if ( boundaryEdgeContainer.find ( bareEdge ) == boundaryEdgeContainer.end() )
                internalEdgeContainer.insert
                ( std::make_pair ( bareEdge, std::make_pair ( volumeContainerIterator->localId(), jEdgeLocalId ) ) );
        }
    }
    return internalEdgeContainer.size();
}


/*
 *******************************************************************************
                            MARKERS HANDLERS
 *******************************************************************************
 */
//! @defgroup marker_policies Used to manage missing handlers


/*! Sets the marker ID of a MeshElementMarked according to the policy of the marker

    @ingroup marker_policies

    It gets the stronger marker of the MeshElementMarked points. The marker
    hierarchy is defined in the MarkerDefinitions.hpp file. It returns the new
    marker id for the MeshElementMarked. If any of the vertices has an unset marker
    the result is an unset marked ID for the MeshElementMarked.

    @sa MarkerDefinitions.hpp
    @warning It overrides the original marker ID.
    @return the new marker ID for geoElement

    @todo LF: It should be made a functor so to give the user a easier way
              to change the policy if needed
 */

template <typename MeshElementMarkedType>
markerID_Type inheritPointsStrongerMarker ( MeshElementMarkedType& geoElement )
{
    ASSERT_PRE ( MeshElementMarkedType::S_nDimensions > 0,
                 "A MeshElementMarked with ndim == 0 cannot inherit marker IDs" );

    geoElement.setMarkerID ( geoElement.point ( 0 ).markerID() );
    for ( ID jPointId = 1; jPointId < MeshElementMarkedType::S_numVertices; ++jPointId )
    {
        geoElement.setStrongerMarker ( geoElement.point ( jPointId ).markerID() );
    }
    return geoElement.markerID();

}


/*! @ingroup marker_policies

//! @brief Sets the marker ID of a MeshElementMarked of dimension greater one

    It gets the weaker marker of the MeshElementMarked points. The marker
    hierarchy is defined in the MarkerDefinitions.hpp file. It returns the new
    marker  id for the MeshElementMarked. If any of the vertices has an unset marker
    the result is an unset marker ID for the MeshElementMarked.

    @sa MarkerDefinitions.hpp
    @warning It overrides the original marker ID.
    @return the new marker ID for geoElement
 */
template <typename MeshElementMarkedType>
markerID_Type inheritPointsWeakerMarker ( MeshElementMarkedType& geoElement )
{
    ASSERT_PRE ( MeshElementMarkedType::S_nDimensions > 0,
                 "A MeshElementMarked with ndim == 0 cannot inherit marker IDs" );

    geoElement.setMarkerID ( geoElement.point ( 0 ).markerID() );
    for ( ID jPointId = 1; jPointId < MeshElementMarkedType::S_numVertices; ++jPointId )
    {
        geoElement.setWeakerMarkerID ( geoElement.point ( jPointId ).markerID() );
    }
    return geoElement.markerID();

}


/*!
    This routine tests if the topological description of boundary face is sane.
    In particular all boundary edges must be adjacent to only 2 surface elements
    and the orientation must be correct.

    @param mesh a mesh
    @param numBoundaryEdges The function also returns the number of boundary edges in numBoundaryEdges.

    @return It it returns 0 if the test has been passed. If not it returns the number of wrong boundary edges.
    @warning numBoundaryEdges is properly set only if the test has been passed.
 */
template <typename MeshType>
UInt testDomainTopology ( MeshType const& mesh, UInt& numBoundaryEdges )
{

    typedef std::set <BareEdge, cmpBareItem<BareEdge> > localTemporaryEdgeContainer_Type;
    localTemporaryEdgeContainer_Type                    localTemporaryEdgeContainer;
    UInt                                                point1Id, point2Id;
    BareEdge                                            bareEdge;
    typename MeshType::facetShape_Type                       faceShape;
    typedef typename MeshType::faces_Type                    faceContainer_Type;
    typedef typename MeshType::face_Type               face_Type;
    localTemporaryEdgeContainer_Type::iterator          edgeContainerIterator;

    typename faceContainer_Type::const_iterator faceContainerIterator = mesh.faceList.begin();

    for ( UInt kFaceId = 0; kFaceId < mesh.numBFaces(); ++kFaceId )
    {
        std::ostringstream errorStream;
        errorStream << " Trying to get not existing face"
                    << kFaceId << " " << mesh.numBFaces();
        ASSERT ( faceContainerIterator != mesh.faceList.end(), errorStream.str().c_str() );

        for ( ID jEdgeLocalId = 0; jEdgeLocalId < face_Type::S_numEdges; ++jEdgeLocalId )
        {
            point1Id = faceShape.edgeToPoint ( jEdgeLocalId, 0 );
            point2Id = faceShape.edgeToPoint ( jEdgeLocalId, 1 );
            // go to global
            point1Id = ( faceContainerIterator->point ( point1Id ) ).localId();
            point2Id = ( faceContainerIterator->point ( point2Id ) ).localId();
            bareEdge = ( makeBareEdge ( point1Id, point2Id ) ).first;

            if ( ( edgeContainerIterator = localTemporaryEdgeContainer.find ( bareEdge ) )
                    == localTemporaryEdgeContainer.end() )
            {
                localTemporaryEdgeContainer.insert ( bareEdge );
                ++numBoundaryEdges;
            }
            else
            {
                localTemporaryEdgeContainer.erase ( edgeContainerIterator );
            }
        }
        ++faceContainerIterator;
    }
    return localTemporaryEdgeContainer.size();
}


/*
 *******************************************************************************
                                MARKERS FIXING
 *******************************************************************************
 */

//! Check whether all markers of a the geometry entities stored in a list are set
template <typename MeshEntityListType>
bool checkIsMarkerSet ( const MeshEntityListType& meshEntityList )
{
    typedef typename MeshEntityListType::const_iterator MeshEntityListTypeConstIterator_Type;
    bool ok ( true );
    for ( MeshEntityListTypeConstIterator_Type meshEntityListIterator = meshEntityList.begin();
            meshEntityListIterator != meshEntityList.end(); ++meshEntityListIterator )
    {
        ok = ( ok & meshEntityListIterator->isMarkerSet() );
    }
    return ok;
}


//! Sets the marker id for all boundary edges by inheriting them from boundary points.
/*!
    The paradigm is that an edge <B>WHOSE MARKER HAS NOT ALREADY BEEN
    SET</B> will get the WEAKER marker ID among its VERTICES. For instance
    if a vertex is assigned to an Essential BC and the other to a Natural
    BC the edge will get the marker ID related to the Natural B.C.

    What is a weaker marker is set in the MarkerPolicy passed through the markers.

    @param mesh A mesh
    @param logStream stream to which a map edgeId -> TimeAdvanceNewmarker will be output
    @param errorStream stream to which error messages will be sent
    @param verbose if false, no messages will be sent to the logStream

    @todo it should take the way to handle missing marker ids as policy
    @todo errorStream is unused
 */
template <typename MeshType>
void
setBoundaryEdgesMarker ( MeshType& mesh, std::ostream& logStream = std::cout,
                         std::ostream& /*errorStream*/ = std::cerr, bool verbose = true )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );

    typename MeshType::edge_Type* edgePtr = 0;
    UInt                  counter ( 0 );

    if ( verbose )
        logStream << "NEW EDGE MARKER MAP" << std::endl
                  << " ID->New Marker" << std::endl;

    for ( ID kEdgeId = 0; kEdgeId < mesh.numBEdges(); ++kEdgeId )
    {
        edgePtr = & ( mesh.edge ( kEdgeId ) );
        if ( edgePtr->isMarkerUnset() )
        {
            inheritPointsWeakerMarker ( *edgePtr );
            if ( verbose )
            {
                logStream << edgePtr->localId() << " -> " << edgePtr->markerID();
                logStream << " ";
                if ( ++counter % 3 == 0 )
                {
                    logStream << std::endl;
                }
            }
        }
    }
    if ( verbose )
    {
        logStream << std::endl;
    }
}


//! Sets the marker ID for all boundary faces by inheriting them from boundary points.
/*!
    The paradigm is that a face <B>WHOSE MARKER ID HAS NOT ALREADY BEEN SET</B> will
    get the WEAKER marker ID among its VERTICES. For instance if a vertex
    is assigned to a Natural BC and the others to a Natural BC the face
    will get the marker ID related to the Natural BC

    @param mesh A mesh
    @param logStream stream to which a map faceId -> TimeAdvanceNewmarker will be output
    @param errorStream stream to which error messages will be sent
    @param verbose if false, no messages will be sent to the logStream

    @todo the way to handle missing IDs should be passed as a policy
 */
template <typename MeshType>
void
setBoundaryFacesMarker ( MeshType& mesh, std::ostream& logStream = std::cout,
                         std::ostream& /*errorStream*/ = std::cerr, bool verbose = true )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );

    typename MeshType::face_Type* facePtr = 0;
    UInt                  counter ( 0 );

    if ( verbose )
    {
        logStream << "**** NEW FACE MARKER MAP **************" << std::endl;
        logStream << " Face ID -> New Marker\tFace ID -> New Marker\tFace ID -> New Marker" << std::endl;
    }

    for ( UInt kFaceId = 0; kFaceId < mesh.numBFaces(); ++kFaceId )
    {
        facePtr = & ( mesh.face ( kFaceId ) );
        if ( facePtr->isMarkerUnset() )
        {
            inheritPointsWeakerMarker ( *facePtr );
            if ( verbose )
            {
                logStream << facePtr->localId() << " -> " << facePtr->markerID();
                logStream << "\t";
                if ( ++counter % 3 == 0 )
                {
                    logStream << std::endl;
                }
            }
        }
    }
    if ( verbose )
    {
        logStream << std::endl;
    }
}


//! It sets the marker ID for Points, by inheriting it from facets.
/*!
    The paradigm is that a point whose marker ID is unset will inherit
    the strongest marker ID of the surrounding facets, with the
    convention that if the marker ID of one of the surrounding facets is null,
    it is ignored.

    @param mesh A mesh
    @param logStream stream to which a map PointId -> MarkerId will be output
    @param errorStream stream to which error messages will be sent
    @param verbose if false, no messages will be sent to the logStream
    @pre The boundary faces must be correctly set and the points boundary flags as well
    @note it does not touch _bPoints. It must be set otherwhise.
 */
template <typename MeshType>
void
setBoundaryPointsMarker ( MeshType& mesh, std::ostream& logStream = std::cout,
                          std::ostream& /*errorStream*/ = std::cerr, bool verbose = false )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );

    // First looks at points whose marker has already been set
    std::vector<bool> isDefinedPointMarker ( mesh.storedPoints(), false );

    typedef typename MeshType::points_Type::iterator pointContainerIterator_Type;
    typedef typename MeshType::facetShape_Type    faceShape_Type;

    std::vector<bool>::iterator isDefinedPointMarkerIterator = isDefinedPointMarker.begin();

    for ( pointContainerIterator_Type pointContainerIterator = mesh.pointList.begin();
            pointContainerIterator != mesh.pointList.end(); ++pointContainerIterator )
    {
        * ( isDefinedPointMarkerIterator++ ) = pointContainerIterator->isMarkerSet();
    }

    typename MeshType::face_Type* facetPtr = 0;
    for ( UInt kFacetId = 0; kFacetId < mesh.numBFaces(); ++kFacetId )
    {
        facetPtr = & ( mesh.boundaryFacet ( kFacetId ) );
        if ( facetPtr->isMarkerSet() )
        {
            for ( UInt jPointId = 0; jPointId < faceShape_Type::S_numPoints; ++jPointId )
            {
                if ( !isDefinedPointMarker[ facetPtr->point ( jPointId ).localId() ] )
                    // A bit involved but it works
                    //todo operate directly on point using setStrongerMarkerID
                {
                    facetPtr->setStrongerMarkerIDAtPoint ( jPointId, facetPtr->markerID() );
                }
            }
        }
    }
    // now the internal
    for ( UInt i = 0; i < mesh.storedPoints(); ++i )
        if (!mesh.point (i).boundary() && !isDefinedPointMarker[i])
        {
            mesh.point (i).setMarkerID (mesh.markerID() );
        }
    UInt counter ( 0 );

    if ( verbose )
    {
        logStream << "**** NEW POINT MARKER MAP **************" << std::endl;
        logStream << " Point ID -> New Marker\tPoint ID -> New Marker\tPoint ID -> New Marker" << std::endl;
        isDefinedPointMarkerIterator = isDefinedPointMarker.begin();
        for ( pointContainerIterator_Type pointContainerIterator = mesh.pointList.begin();
                pointContainerIterator != mesh.pointList.end(); ++pointContainerIterator )
        {
            if ( *isDefinedPointMarkerIterator++ )
            {
                logStream << pointContainerIterator->localId() << " -> " <<
                          pointContainerIterator->markerID();
                logStream << "\t";
                if ( ++counter % 3 )
                {
                    logStream << std::endl;
                }
            }
        }
        logStream << std::endl;
    }
}


/*
 *******************************************************************************
                                FIXING ID AND COUNTERS
 *******************************************************************************
 */
//! @brief Verifies if a list of mesh entities have the ID properly set.
/* More precisely, the localId() must correspond to the position of the entity
   in the list.

   @pre The template argument MeshEntityListType must be a stl
   compliant container and its elements must have the method localId().
 */
template <typename MeshEntityListType>
bool checkId ( const MeshEntityListType& meshEntityList )
{
    typedef typename MeshEntityListType::const_iterator MeshEntityListTypeConstIterator_Type;
    bool ok ( true );
    UInt counter ( 0 );
    for ( MeshEntityListTypeConstIterator_Type meshEntityListIterator = meshEntityList.begin();
            meshEntityListIterator != meshEntityList.end() && ok; ++meshEntityListIterator, ++counter )
    {
        ok = ok && ( meshEntityListIterator->localId() == counter );
    }
    return ok;
}


//! @brief Fixes a a list of mesh entities so that the ID is properly set.
/* @post  The localId will correspond to the position of the entity
   in the list.

   @pre The template argument MeshEntityListType must be a stl
   compliant container and its elements must have the method setLocalId().
 */
template <typename MeshEntityListType>
void fixId ( MeshEntityListType& meshEntityList )
{
    UInt counter ( 0 );
    typedef typename MeshEntityListType::iterator Iter;
    for ( Iter meshEntityListIterator = meshEntityList.begin();
            meshEntityListIterator != meshEntityList.end(); ++meshEntityListIterator )
    {
        meshEntityListIterator->setLocalId ( counter++ );
    }
}


/*! @brief Fixes boundary points counter
  It fix the boundary points counter by counting
  how many points have the boundary flag set.
  It also resets the boundary points list.

  @pre It assumes that the points have the boundary flag correctly set
 */
template <typename MeshType>
void
setBoundaryPointsCounters ( MeshType& mesh )
{

    UInt boundaryPointCounter ( 0 );
    UInt boundaryVertexCounter ( 0 );

    mesh._bPoints.clear();

    for ( UInt kPointId = 0; kPointId < mesh.numVertices(); ++kPointId )
    {
        if ( mesh.isBoundaryPoint ( kPointId ) )
        {
            ++boundaryPointCounter;
            ++boundaryVertexCounter;
        }
    }

    for ( UInt kPointId = mesh.numVertices(); kPointId < mesh.storedPoints(); ++kPointId )
    {
        if ( mesh.isBoundaryPoint ( kPointId ) )
        {
            ++boundaryPointCounter;
        }
    }

    mesh.numBVertices() = boundaryVertexCounter;
    mesh.setNumBPoints ( boundaryPointCounter );
    mesh._bPoints.clear();
    mesh._bPoints.reserve ( boundaryPointCounter );

    for ( UInt kPointId = 0; kPointId < mesh.storedPoints(); ++kPointId )
    {
        if ( mesh.isBoundaryPoint ( kPointId ) )
        {
            mesh._bPoints.push_back ( &mesh.point ( kPointId ) );
        }
    }
}


/*
 *******************************************************************************
BOUNDARY INDICATOR FIXING
 *******************************************************************************
 */
//! It fixes boundary flag on points laying on boundary faces.
/*!
  @param mesh a mesh
  @param logStream logging stream
  @param errorStream error stream
  @param verbose If true you have a verbose output

  @pre mesh point list must exists and boundary face list must have been set properly.
 */
template <typename MeshType>
void
fixBoundaryPoints ( MeshType& mesh, std::ostream& logStream = std::cout,
                    std::ostream& /* errorStream */ = std::cerr, bool verbose = true )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );

    ASSERT_PRE ( mesh.numPoints() > 0, "The point list should not be empty" );
    ASSERT_PRE ( mesh.numBFaces() > 0,
                 "The boundary faces list should not be empty" );

    typedef typename MeshType::faces_Type     faceContainer_Type;
    typedef typename MeshType::facetShape_Type faceShape_Type;

    if ( verbose )
    {
        logStream << "Fixing BPoints" << std::endl;
    }
    std::vector<bool> boundaryPoints (mesh.numPoints(), false);
    // I may have launched the program for a P2 mesh
    // yet not all the points are there
    UInt numitems;
    if (mesh.storedPoints() == mesh.numVertices() )
    {
        numitems = faceShape_Type::S_numVertices;
    }
    else
    {
        numitems = faceShape_Type::S_numPoints;
    }

    for ( UInt kFaceId = 0; kFaceId < mesh.numBFaces(); ++kFaceId )
        for ( UInt jPointId = 0; jPointId < numitems; ++jPointId )
        {
            boundaryPoints[mesh.boundaryFacet (kFaceId).point (jPointId).localId()] = true;
        }
    for (ID  kPointId = 0; kPointId < mesh.storedPoints() ; ++kPointId )
        if (boundaryPoints[kPointId])
        {
            mesh.point (kPointId).setFlag (EntityFlags::PHYSICAL_BOUNDARY);
        }
        else
        {
            mesh.point (kPointId).unSetFlag (EntityFlags::PHYSICAL_BOUNDARY);
        }
    // anihilate
    boundaryPoints.clear();
    std::vector<bool>().swap (boundaryPoints);
    // Fix now the number of vertices/points and reset _bpoints list in the mesh
    setBoundaryPointsCounters ( mesh );
}




/*
 *******************************************************************************
 UTILITIES TO VERIFY/CREATE FACES/EDGES
 *******************************************************************************
*/
/**
   @brief It rearranges the faces stored in the mesh.

   It makes sure that

   -# Faces have the boundary flag correctly set
   -# Boundary faces are stored first

   It uses the information given by findBoundaryFaces, so it relies ONLY on the mesh topology

    It works on a mesh where faces have already been found! So it is not
    meant to be used for finding the boundary faces. Use the general utility BuildFaces
    for that purpose. Its main role, as the name says, is to make sure that faces are well ordered.
    It does not verify the consitency of the face-to-adjacentVolume information. Use fixBoundaryFaces for this purpose

    @param[out] mesh a mesh

    @param[out] logStream stream that will receive all information regarding the markers

    @param[out] errorStream stream for error messages

    @param[out] sw A switch that will contain information on what has been done
    Possible values are
    <ol>
    <li>FACES_REORDERED</li>
    </ol>

    @param numFaces[out] It returns the number of faces found by the function

    @param numBoundaryFaces[out] It returns the number of boundary faces found by the function

    @param[in] verbose if false nothing is written to logStream

    @param[out] externalFaceContainer. If not NULL it is a pointer to an external map of boundary faces, already
      produced by a call to findBoundaryFaces(). This parameter may be used to save a lot of computational work, since
      findBoundaryFaces() is rather expensive.

    @pre  Mesh must contain the faces, at least the boundary ones
 */


template <class MeshType>
bool rearrangeFaces ( MeshType& mesh,
                      std::ostream& logStream,
                      std::ostream& errorStream,
                      Switch& sw,
                      UInt& numFaces,
                      UInt& numBoundaryFaces,
                      bool verbose = false,
                      temporaryFaceContainer_Type* externalFaceContainer = 0 )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );
    typedef typename MeshType::faces_Type faceContainer_Type;
    typedef typename MeshType::face_Type face_Type;

    UInt                                  point1Id, point2Id, point3Id, point4Id;
    BareFace                              bareFace;
    typename faceContainer_Type::iterator faceContainerIterator;
    temporaryFaceContainer_Type*          boundaryFaceContainerPtr;
    temporaryFaceContainer_Type::iterator boundaryFaceContainerIterator;
    UInt                                  numInternalFaces;
    bool                                  externalContainerIsProvided ( false );

    if ( (externalContainerIsProvided = ( externalFaceContainer != 0 ) ) )
    {
        boundaryFaceContainerPtr = externalFaceContainer;
        numBoundaryFaces = boundaryFaceContainerPtr->size();
    }
    else
    {
        boundaryFaceContainerPtr = new temporaryFaceContainer_Type;
        numBoundaryFaces = findBoundaryFaces ( mesh, *boundaryFaceContainerPtr, numInternalFaces );
        numFaces = numBoundaryFaces + numInternalFaces;
    }


    bool notEnough = mesh.storedFaces() < numBoundaryFaces;



    if ( notEnough )
    {
        errorStream << "WARNING: number of B. Faces stored smaller" << std::endl;
        if (verbose)
        {
            logStream << "WARNING: number of B. Faces stored smaller" << std::endl;
        }
        errorStream << "         than the number of boundaryFaces found  and build is not set"
                    << std::endl;
        errorStream << "ABORT condition in rearrangeFaces" << std::endl;
        sw.create ( "BFACE_STORED_MISMATCH", true );
        return false;
    }

    if ( mesh.numBFaces() == 0 )
    {
        errorStream << "ERROR: Boundary Element counter was not set" << std::endl;
        if ( verbose )
        {
            logStream << "ERROR: Boundary Element counter was not set" << std::endl;
        }
        errorStream << "I Cannot proceed because the situation is ambiguous"
                    << std::endl;
        errorStream << "Please check and eventually either: (a) call buildFaces()" << std::endl;
        errorStream << "or (b) set the correct number of boundaryFaces in the mesh using mesh.numBFaces()" << std::endl;
        errorStream << "ABORT" << std::endl;
        sw.create ( "BELEMENT_COUNTER_UNSET", true );
        return false;
    }

    if ( mesh.numBFaces() != numBoundaryFaces )
    {
        errorStream << "WARNING: Boundary face counter in mesh is set to "
                    << mesh.numBFaces() << std::endl;
        errorStream << "         while I have found " << numBoundaryFaces
                    << " boundary elements in mesh." << std::endl;
        errorStream << "         Please check... I continue anyway" << std::endl;
        sw.create ( "BFACE_COUNTER_MISMATCH", true );
    }


    for ( faceContainerIterator = mesh.faceList.begin(); faceContainerIterator != mesh.faceList.end();
            ++faceContainerIterator )
    {
        point1Id = ( faceContainerIterator->point ( 0 ) ).localId();
        point2Id = ( faceContainerIterator->point ( 1 ) ).localId();
        point3Id = ( faceContainerIterator->point ( 2 ) ).localId();

        if ( MeshType::facetShape_Type::S_numVertices == 4 )
        {
            point4Id = ( faceContainerIterator->point ( 3 ) ).localId();
            bareFace = ( makeBareFace ( point1Id, point2Id, point3Id, point4Id ) ).first;
        }
        else
        {
            bareFace = ( makeBareFace ( point1Id, point2Id, point3Id ) ).first;
        }
        boundaryFaceContainerIterator = boundaryFaceContainerPtr->find ( bareFace );
        if ( boundaryFaceContainerIterator == boundaryFaceContainerPtr->end() )
        {
            faceContainerIterator->unSetFlag (EntityFlags::PHYSICAL_BOUNDARY);
        }
        else
        {
            faceContainerIterator->setFlag (EntityFlags::PHYSICAL_BOUNDARY);
        }
    }
    mesh.faceList.reorderAccordingToFlag (EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet);

    return true;
}

//! It fixes boundary faces so that they are consistently numbered with volumes.

/*! An important step for building degrees of freedom on faces.  It also
    fixes other face related data.
    It works on a mesh where boundary faces have already been found! So it is not
    meant to be used for finding the boundary faces. Use the general utility BuildFaces
    for that purpose. Its main role, as the name says, is to fix a partially broken mesh.
    In particular, it assures that the boundary faces are correctly set w.r.t. the adjacent volumes


    @param[out] mesh a mesh

    @param[out] logStream stream that will receive all information regarding the markers

    @param[out] errorStream stream for error messages

    @param[out] sw A switch that will contain information on what has been done
    Possible values are
    <ol>
    <li>NUM_FACES_MISMATCH</li>
    <li>FIXED_FACE_COUNTER</li>
    <li>BFACE_MISSING</li>
    <li>BFACE_STORED_MISMATCH</li>
    <li>BFACE_COUNTER_UNSET</li>
    <li>BFACE_STORED_MISMATCH</li>
    <li>FIXED_MAX_NUM_FACES</li>
    </ol>

    @param numFaces[out] It returns the number of faces found by the function

    @param numBoundaryFaces[out] It returns the number of boundary faces found by the function

    @param fixMarker[in] If set to the true value, all faces without a markerFlag set will inherit it from the points.
        todo remove this parameter (unused)

    @param[in] verbose if false nothing is written to logStream

    @param[out] externalFaceContainer. If not NULL it is a pointer to an external map of boundary faces, already
      produced by a call to findBoundaryFaces(). This parameter may be used to save a lot of computational work, since
      findBoundaryFaces() is rather expensive.

    @pre Boundary faces list must be properly set.
    @todo The policy to treat missing markers should be passed in the argument, so to allow changes
 */

template <class MeshType>
bool fixBoundaryFaces ( MeshType& mesh,
                        std::ostream& logStream,
                        std::ostream& errorStream,
                        Switch& sw,
                        UInt& numFaces,
                        UInt& numBoundaryFaces,
                        bool /* fixMarker */ = false,
                        bool verbose = false,
                        temporaryFaceContainer_Type* externalFaceContainer = 0 )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );
    typedef typename MeshType::volumes_Type volumeContainer_Type;
    typedef typename MeshType::volume_Type volume_Type;
    typedef typename MeshType::faces_Type faceContainer_Type;
    typedef typename MeshType::face_Type face_Type;

    UInt                                  point1Id, point2Id, point3Id, point4Id;
    BareFace                              bareFace;
    volume_Type*                           volumePtr;
    typename faceContainer_Type::iterator faceContainerIterator;
    typename MeshType::elementShape_Type        volumeShape;
    temporaryFaceContainer_Type*          boundaryFaceContainerPtr;
    temporaryFaceContainer_Type::iterator boundaryFaceContainerIterator;
    std::pair<ID, ID>                     volumeIdToLocalFaceIdPair;
    ID                                    jFaceLocalId;
    ID                                    volumeId;
    UInt                                  numInternalFaces;
    bool                                  notfound ( false );
    bool                                  externalContainerIsProvided ( false );

    if ( (externalContainerIsProvided = ( externalFaceContainer != 0 ) ) )
    {
        boundaryFaceContainerPtr = externalFaceContainer;
        numBoundaryFaces = boundaryFaceContainerPtr->size();
    }
    else
    {
        boundaryFaceContainerPtr = new temporaryFaceContainer_Type;
        numBoundaryFaces = findBoundaryFaces ( mesh, *boundaryFaceContainerPtr, numInternalFaces );
        numFaces = numBoundaryFaces + numInternalFaces;
    }


    bool notEnough = mesh.storedFaces() < numBoundaryFaces;



    if ( notEnough )
    {
        errorStream << "WARNING: number of B. Faces stored smaller" << std::endl;
        errorStream << "         than the number of boundaryFaces found  and build is not set"
                    << std::endl;
        errorStream << "POSSIBLE ERROR" << std::endl;
        sw.create ( "BFACE_STORED_MISMATCH", true );
    }

    if ( mesh.numBFaces() == 0 )
    {
        errorStream << "ERROR: Boundary Element counter was not set" << std::endl;
        errorStream << "I Cannot proceed because the situation is ambiguous"
                    << std::endl;
        errorStream << "Please check and eventually either: (a) call buildBoundaryFaces()" << std::endl;
        errorStream << "or (b) set the correct number of boundaryFaces in the mesh using mesh.numBFaces()" << std::endl;
        errorStream << "ABORT" << std::endl;
        sw.create ( "BFACE_COUNTER_UNSET", true );
    }

    if ( mesh.numBFaces() != numBoundaryFaces )
    {
        errorStream << "WARNING: Boundary face counter in mesh is set to "
                    << mesh.numBFaces() << std::endl;
        errorStream << "         while I have found " << numBoundaryFaces
                    << " boundary elements in mesh." << std::endl;
        errorStream << "         Please check... I continue anyway" << std::endl;
        sw.create ( "BFACE_COUNTER_MISMATCH", true );
    }

    if ( verbose )
    {
        logStream << "**** Fixed Marker Flags for Boundary Faces ***" << std::endl;
        logStream << " (it only contains those that were fixed because unset !)"
                  << std::endl;
        logStream << "id->marker   id->marker  id->marker" << std::endl;
    }

    UInt counter ( 0 );

    faceContainerIterator = mesh.faceList.begin();
    for ( UInt facid = 0; facid < mesh.numBFaces(); ++facid )
    {
        point1Id = ( faceContainerIterator->point ( 0 ) ).localId();
        point2Id = ( faceContainerIterator->point ( 1 ) ).localId();
        point3Id = ( faceContainerIterator->point ( 2 ) ).localId();
        if ( MeshType::facetShape_Type::S_numVertices == 4 )
        {
            point4Id = ( faceContainerIterator->point ( 3 ) ).localId();
            bareFace = ( makeBareFace ( point1Id, point2Id, point3Id, point4Id ) ).first;
        }
        else
        {
            bareFace = ( makeBareFace ( point1Id, point2Id, point3Id ) ).first;
        }
        boundaryFaceContainerIterator = boundaryFaceContainerPtr->find ( bareFace );
        if ( boundaryFaceContainerIterator == boundaryFaceContainerPtr->end() )
        {
            if (verbose)
            {
                if ( MeshType::facetShape_Type::S_numVertices == 3 )
                {
                    errorStream << "Face " << point1Id << " " << point2Id << " " << point3Id;
                }
                else
                {
                    errorStream << "Face " << point1Id << " " << point2Id << " " << point3Id << " " << point4Id;
                }
                errorStream << " stored as boundary face, it's not!" << std::endl;
            }
            notfound = true;
        }
        else
        {
            volumeIdToLocalFaceIdPair = boundaryFaceContainerIterator->second;
            volumeId = volumeIdToLocalFaceIdPair.first; // Element ID
            volumePtr = &mesh.volume ( volumeId ); // Element
            jFaceLocalId = volumeIdToLocalFaceIdPair.second;       // The local ID of face on element
            // Reset face point definition to be consistent with face.
            for ( UInt kPointId = 0; kPointId < face_Type::S_numPoints; ++kPointId )
            {
                faceContainerIterator->setPoint ( kPointId, volumePtr->point ( volumeShape.faceToPoint ( jFaceLocalId, kPointId ) ) );
            }
            // Correct extra info
            faceContainerIterator->firstAdjacentElementIdentity() = volumeId;
            faceContainerIterator->firstAdjacentElementPosition() = jFaceLocalId;
            faceContainerIterator->secondAdjacentElementIdentity() = NotAnId;
            faceContainerIterator->secondAdjacentElementPosition() = NotAnId;

            if ( faceContainerIterator->isMarkerUnset() )
            {
                inheritPointsWeakerMarker ( *faceContainerIterator );
                if ( verbose )
                {
                    logStream << faceContainerIterator->localId() << " -> " <<
                              faceContainerIterator->markerID();
                    logStream << " ";
                    if ( ++counter % 3 == 0 )
                    {
                        logStream << std::endl;
                    }
                }
            }
            // Take out face from temporary container
            boundaryFaceContainerPtr->erase ( boundaryFaceContainerIterator );
        }
        ++faceContainerIterator;
    }

    if ( !externalContainerIsProvided )
    {
        delete boundaryFaceContainerPtr;
    }

    if ( notfound )
    {
        errorStream << "WARNING: At least one boundary face has not been found on the list stored in MeshType\n";
        sw.create ( "BFACE_MISSING", true );
    }

    if ( verbose )
    {
        logStream << std::endl << "  *****  END OF LIST ****" << std::endl;
    }

    // Here I calculate the number of faces,

    if ( mesh.numFaces() != numFaces )
    {
        errorStream << "WARNING: faces counter in mesh should be " << numFaces
                    << std::endl;
        errorStream << "         (boundaryFaceContainerPtr->size()+numInternalFaces)" << std::endl;
        errorStream << "         it is instead " << mesh.numFaces() << std::endl;
        sw.create ( "NUM_FACES_MISMATCH", true );
    }
    mesh.setLinkSwitch ( std::string ( "HAS_BOUNDARY_FACETS" ) );

    return true;
}



//! Builds faces
/*! This is a major function.
  It may be used to build or partially build faces. So it may operate also on a mesh where
  not all boundary faces are set.
  Function may alternatively be used to build the compulsory boundary
  faces, all the mesh faces, or just add to an existing list of just boundary
  faces the internal ones. It will not destroy the basic info (marker id, etc) contained in the
  faces list already stored in the meash. So if you want to build everything from scratch you need
  to clear it first. It (re)build the adjacent volume info, sets the boundary flags and ensures that
  boundary faces are stored first.


  @param mesh A mesh

  @param logStream Log stream for information on the newly created markers

  @param errorStream Error stream

  @param numBoundaryFaces It returns the number of boundary faces

  @param numInternalFaces It returns the number of internal faces (only if externalFaceContainer is not provided!)

  @param buildBoundaryFaces if true the function builds boundary faces

  @param buildInternalFaces if true the function builds internal faces

  @param verbose If true markerFlags info is written on logStream.

  @param externalFaceContainer. If not NULL it is a pointer to an external map of boundary faces, already
  produced by a call to findBoundaryFaces(). This parameter may be used to save a lot of computational work, since
  findBoundaryFaces() is rather expensive.

  @pre If buildInternalFaces=true and buildBoundaryFaces=false the mesh must contain a proper list
  of boundary faces

  @note By setting buildInternalFaces=true and buildBoundaryFaces=true the function just fixes the counters
  with the number of faces in the mesh

  @return true if successful
 */

#ifndef TWODIM
template <class MeshType>
bool buildFaces ( MeshType& mesh,
                  std::ostream& logStream,
                  std::ostream& errorStream,
                  UInt& numBoundaryFaces,
                  UInt& numInternalFaces,
                  bool buildBoundaryFaces = true,
                  bool buildInternalFaces = false,
                  bool verbose = false,
                  temporaryFaceContainer_Type* externalFaceContainer = 0 )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );
    UInt                                  point1Id, point2Id, point3Id, point4Id;
    typename MeshType::elementShape_Type   volumeShape;
    typedef typename MeshType::volumes_Type    volumeContainer_Type;
    typedef typename MeshType::volume_Type volume_Type;
    typedef typename MeshType::faces_Type      faceContainer_Type;
    typedef typename MeshType::face_Type   face_Type;
    volume_Type*                          volumePtr;
    temporaryFaceContainer_Type*          boundaryFaceContainerPtr;
    temporaryFaceContainer_Type::iterator boundaryFaceContainerIterator;
    bool                                  externalContainerIsProvided ( false );

    std::pair<ID, ID>                     volumeIdToLocalFaceIdPair;
    ID                                    jFaceLocalId, newFaceId;
    ID                                    volumeId;
    std::map<BareFace, ID>                 existingFacesMap;
    std::map<BareFace, ID>::iterator       existingFacesMap_It;
    std::pair<std::map<BareFace, ID>::iterator, bool>   existingFacesMap_insert;
    bool                                  faceExists (false);
    // Handle boundary face container
    if ( (externalContainerIsProvided = ( externalFaceContainer != 0 ) ) )
    {
        boundaryFaceContainerPtr = externalFaceContainer;
        numBoundaryFaces = boundaryFaceContainerPtr->size();
    }
    else
    {
        boundaryFaceContainerPtr = new temporaryFaceContainer_Type;
        numBoundaryFaces = findBoundaryFaces ( mesh, *boundaryFaceContainerPtr, numInternalFaces );
    }
    // Maybe we have already faces stored, save them!
    for ( UInt jFaceId = 0; jFaceId < mesh.faceList.size(); ++jFaceId )
    {
        point1Id = ( mesh.faceList[ jFaceId ].point ( 0 ) ).localId();
        point2Id = ( mesh.faceList[ jFaceId ].point ( 1 ) ).localId();
        point3Id = ( mesh.faceList[ jFaceId ].point ( 2 ) ).localId();
        if ( MeshType::facetShape_Type::S_numVertices == 4 )
        {
            point4Id = ( mesh.faceList[ jFaceId ].point ( 3 ) ).localId();
            existingFacesMap_insert = existingFacesMap.insert (
                                          std::make_pair ( makeBareFace ( point1Id, point2Id, point3Id, point4Id).first, jFaceId )
                                      );

        }
        else
        {
            existingFacesMap_insert = existingFacesMap.insert (
                                          std::make_pair ( makeBareFace ( point1Id, point2Id, point3Id).first, jFaceId )
                                      );
        }
        if (! existingFacesMap_insert.second)
        {
            errorStream << "ERROR in BuildFaces. Mesh stores two identical faces" << std::endl;
            if ( !externalContainerIsProvided )
            {
                delete boundaryFaceContainerPtr;
            }
            return false;
        }
    }


    if ( buildBoundaryFaces )
    {
        mesh.setNumBFaces ( numBoundaryFaces );
    }
    if ( !buildInternalFaces )
    {
        mesh.setMaxNumFaces ( std::max (numBoundaryFaces, static_cast<UInt> (mesh.faceList.size() ) ), false );
        mesh.setNumFaces ( numInternalFaces + numBoundaryFaces );
    }
    else
    {
        mesh.setMaxNumFaces ( numInternalFaces + numBoundaryFaces, true );
    }

    face_Type face;

    if ( buildBoundaryFaces )
    {

        if ( verbose )
        {
            logStream << "**** Marker Flags for Newly Created Boundary Faces ***"
                      << std::endl;
            logStream << "id->marker   id->marker  id->marker" << std::endl;
        }

        for ( boundaryFaceContainerIterator = boundaryFaceContainerPtr->begin();
                boundaryFaceContainerIterator != boundaryFaceContainerPtr->end(); ++boundaryFaceContainerIterator )
        {
            existingFacesMap_It = existingFacesMap.find (boundaryFaceContainerIterator->first);
            if (existingFacesMap_It != existingFacesMap.end() )
            {
                faceExists = true;
                face = mesh.faceList[existingFacesMap_It->second];
                existingFacesMap.erase (existingFacesMap_It);
            }
            else
            {
                faceExists = false;
                face = face_Type();
                face.setId ( mesh.faceList.size() );
            }
            volumeIdToLocalFaceIdPair = boundaryFaceContainerIterator->second;
            volumeId = volumeIdToLocalFaceIdPair.first; // Element ID
            volumePtr = &mesh.volume ( volumeId ); // Element
            jFaceLocalId = volumeIdToLocalFaceIdPair.second;       // The local ID of face on element

            for ( UInt kPointId = 0; kPointId < face_Type::S_numPoints; ++kPointId )
            {
                face.setPoint ( kPointId, volumePtr->point ( volumeShape.faceToPoint ( jFaceLocalId, kPointId ) ) );
            }
            // Add extra info
            face.firstAdjacentElementIdentity() = volumeId;
            face.firstAdjacentElementPosition() = jFaceLocalId;
            face.secondAdjacentElementIdentity() = NotAnId;
            face.secondAdjacentElementPosition() = NotAnId;
            // Get marker value
            if ( face.isMarkerUnset() )
            {
                inheritPointsWeakerMarker ( face );
            }
            face.setBoundary (true);
            if (faceExists)
            {
                // reset the existing face with new info
                newFaceId                = face.localId();
                mesh.setFace (face, newFaceId);
            }
            else
            {
                // The face is new, add it to the mesh
                newFaceId = mesh.addFace ( face).localId();
            }
            if ( verbose )
            {
                if ( newFaceId % 3 == 0 )
                {
                    logStream << std::endl;
                }
                logStream << newFaceId << " -> " << face.markerID();
                logStream << " ";
            }
        }
        mesh.setLinkSwitch ( std::string ( "HAS_BOUNDARY_FACETS" ) );
        if ( ! buildInternalFaces )
        {
            mesh.unsetLinkSwitch ( std::string ( "HAS_ALL_FACETS" ) );
        }
        mesh.setLinkSwitch ( std::string ( "FACETS_HAVE_ADIACENCY" ) );
    }

    if ( !externalContainerIsProvided )
    {
        delete boundaryFaceContainerPtr;
    }
    // All possibly remaining faces are necessarly internal
    for (existingFacesMap_It = existingFacesMap.begin(); existingFacesMap_It != existingFacesMap.end();
            ++existingFacesMap_It)
    {
        mesh.faceList[existingFacesMap_It->second].setBoundary (false);
    }

    // If there where faces stored originally I need to be sure that bfaces go first!
    // I need to do it now because of the tests I do later
    if (!existingFacesMap.empty() )
    {
        mesh.faceList.reorderAccordingToFlag (EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet);
    }

    if ( ! buildInternalFaces )
    {
        // I am done
        // Trim BFaces, memory does not come for free!
        mesh.faceList.trim();

        return true;
    }



    if ( !buildBoundaryFaces )
    {
        if ( mesh.storedFaces() < mesh.numBFaces() )
        {
            errorStream << "ERROR: mesh has not boundary faces, cannot just create internal ones!!!" << std::endl;
            errorStream << "ABORT CONDITION" << std::endl;
            return false;
        }
    }


    /*
      I may get rid of the boundaryFaces container. Unfortunately now I need a more
      complex structure, a MeshElementBareHandler, in order to generate the internal
      faces id. An alternative would be to use the point data to identify
      boundary faces as the ones with all point on the boundary. Yet in this
      function we do not want to use a priori information, so that it might
      work even if the points boundary flag is not properly set.
     */

    MeshElementBareHandler<BareFace> bareFaceHandler;
    std::pair<UInt, bool> faceIdToBoolPair;
    std::pair<BareFace, bool> _face;
    existingFacesMap.clear();
    for ( UInt jFaceId = 0; jFaceId < mesh.faceList.size(); ++jFaceId )
    {
        point1Id = ( mesh.faceList[ jFaceId ].point ( 0 ) ).localId();
        point2Id = ( mesh.faceList[ jFaceId ].point ( 1 ) ).localId();
        point3Id = ( mesh.faceList[ jFaceId ].point ( 2 ) ).localId();
        if ( MeshType::facetShape_Type::S_numVertices == 4 )
        {
            point4Id = ( mesh.faceList[ jFaceId ].point ( 3 ) ).localId();
            _face = makeBareFace ( point1Id, point2Id, point3Id, point4Id );
        }
        else
        {
            _face = makeBareFace ( point1Id, point2Id, point3Id );
        }
        // Store only bfaces by now so if I not find the face is
        // certainly an internal face
        if (mesh.faceList[ jFaceId ].boundary() )
        {
            bareFaceHandler.addIfNotThere ( _face.first );
        }
        else
            // I need to track the numbering
        {
            existingFacesMap.insert (std::make_pair (_face.first, jFaceId) );
        }
    }
    if (bareFaceHandler.howMany() > numBoundaryFaces)
    {
        errorStream << "ERROR in BuildFaces. Not all boundary faces found, very strange" << std::endl;
        errorStream << "ABORT CONDITION" << std::endl;
        return false;
    }

    for ( typename volumeContainer_Type::iterator volumeContainerIterator = mesh.volumeList.begin();
            volumeContainerIterator != mesh.volumeList.end(); ++volumeContainerIterator )
    {
        volumeId = volumeContainerIterator->localId();
        for ( UInt jFaceLocalId = 0; jFaceLocalId < mesh.numLocalFaces(); jFaceLocalId++ )
        {
            point1Id = volumeShape.faceToPoint ( jFaceLocalId, 0 );
            point2Id = volumeShape.faceToPoint ( jFaceLocalId, 1 );
            point3Id = volumeShape.faceToPoint ( jFaceLocalId, 2 );
            // go to global
            point1Id = ( volumeContainerIterator->point ( point1Id ) ).localId();
            point2Id = ( volumeContainerIterator->point ( point2Id ) ).localId();
            point3Id = ( volumeContainerIterator->point ( point3Id ) ).localId();
            if ( MeshType::facetShape_Type::S_numVertices == 4 )
            {
                point4Id = volumeShape.faceToPoint ( jFaceLocalId, 3 );
                point4Id = ( volumeContainerIterator->point ( point4Id ) ).localId();
                _face = makeBareFace ( point1Id, point2Id, point3Id, point4Id );
            }
            else
            {
                _face = makeBareFace ( point1Id, point2Id, point3Id );
            }
            faceIdToBoolPair = bareFaceHandler.addIfNotThere ( _face.first );
            if ( faceIdToBoolPair.second )
            {
                // a new face It must be internal.
                existingFacesMap_It = existingFacesMap.find (_face.first);
                if ( existingFacesMap_It != existingFacesMap.end() )
                {
                    faceExists = true;
                    face = mesh.faceList[existingFacesMap_It->second];
                }
                else
                {
                    faceExists = false;
                    face = face_Type();
                    face.setId ( mesh.faceList.size() );
                }

                for ( UInt kPointId = 0; kPointId < face_Type::S_numPoints; ++kPointId )
                {
                    face.setPoint ( kPointId, volumeContainerIterator->point ( volumeShape.faceToPoint ( jFaceLocalId, kPointId ) ) );
                }
                face.firstAdjacentElementIdentity() = volumeId;
                face.firstAdjacentElementPosition() = jFaceLocalId;
                // Marker is unset
                if (!faceExists)
                {
                    face.setMarkerID (face.nullMarkerID() );
                }
                face.setBoundary (false);
                if (faceExists)
                {
                    mesh.setFace (face, face.localId() );
                    // Add it so we can recover the numbering later on!
                    existingFacesMap.insert ( std::make_pair (_face.first, face.localId() ) );
                }
                else
                {
                    mesh.addFace ( face);
                    // Add it so we can recover the numbering
                    existingFacesMap.insert ( std::make_pair (_face.first, mesh.lastFace().localId() ) );
                }
            }
            else
            {
                if ( faceIdToBoolPair.first > numBoundaryFaces )  // internal
                {
                    existingFacesMap_It = existingFacesMap.find (_face.first);
                    mesh.faceList ( existingFacesMap_It->second).secondAdjacentElementIdentity() = volumeId;
                    mesh.faceList ( existingFacesMap_It->second).secondAdjacentElementPosition() = jFaceLocalId;
                }
            }
        }
    }
    mesh.setLinkSwitch ( std::string ( "HAS_ALL_FACETS" ) );
    return true;
}
#endif


//! It builds edges.
/*!
    This function may alternatively be used to build the boundary edges,
    all the mesh edges, or just add the internal edges to an existing list of
    just boundary edges.

    @param mesh A mesh

    @param logStream Log stream for information on the newly created markers for boundary edges

    @param errorStream Error stream

    @param numBoundaryEdgesFound Returns the number of boundary edges

    @param numInternalEdgesFound Returns the number of internal edges

    @param buildBoundaryEdges if true the function builds boundary edges

    @param buildInternalEdges if true the function builds internal edges

    @param verbose. If true markerFlags info is written on logStream.

    @param externalEdgeContainer. If not NULL it is a pointer to an external map of bondary edges, already
    produced by a call to findBoundaryEdges(). This parameter may be used to save al lot of computational work, since
    findBoundaryEdges() is rather expensive.

    @return true if successful

    @pre If buildInternalEdges=true and buildBoundaryEdges=false the mesh must contain a proper list
    of boundary edges

    @pre The mesh must contain a proper list of boundary faces

    @note By setting buildInternalEdges=true and buildBoundaryEdges=true the function just fixes the counters
    with the number of edges in the mesh
 */

template <typename MeshType>
bool buildEdges ( MeshType& mesh,
                  std::ostream& logStream,
                  std::ostream& errorStream,
                  UInt& numBoundaryEdgesFound,
                  UInt& numInternalEdgesFound,
                  bool buildBoundaryEdges = true,
                  bool buildInternalEdges = false,
                  bool verbose = false,
                  temporaryEdgeContainer_Type* externalEdgeContainer = 0 )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );
    typedef typename MeshType::volumes_Type volumeContainer_Type;
    typedef typename MeshType::faces_Type faceContainer_Type;
    typedef typename MeshType::volume_Type volume_Type;
    typedef typename MeshType::elementShape_Type volumeShape_Type;
    typedef typename MeshType::edges_Type edges_Type;
    typedef typename MeshType::edge_Type edge_Type;
    typedef typename MeshType::face_Type face_Type;
    typedef typename MeshType::facetShape_Type faceShape_Type;
    typedef typename MeshType::edges_Type::iterator Edges_Iterator;
    typename MeshType::face_Type* facePtr;


    std::map<BareEdge, ID> existingEdges;
    typedef std::map<BareEdge, ID>::iterator ExistingEdges_Iterator;
    ExistingEdges_Iterator existingEdges_It;
    bool edgeExists (false);

    temporaryEdgeContainer_Type* temporaryEdgeContainer;
    temporaryEdgeContainer_Type edgeContainer;
    std::pair<ID, ID> faceIdToLocalEdgeIdPair;
    ID jEdgeLocalId, newEdgeId;
    ID faceId;
    BareEdge bareEdge;

    bool externalContainerIsProvided ( false );


    if ( (externalContainerIsProvided = ( externalEdgeContainer != 0 ) ) )
    {
        temporaryEdgeContainer = externalEdgeContainer;
        numBoundaryEdgesFound = temporaryEdgeContainer->size();
    }
    else
    {
        temporaryEdgeContainer = new temporaryEdgeContainer_Type;
        numBoundaryEdgesFound = findBoundaryEdges ( mesh, *temporaryEdgeContainer );
    }

    numInternalEdgesFound = findInternalEdges ( mesh, *temporaryEdgeContainer, edgeContainer );
    // free some memory if not needed!
    // Dump exisitng edges
    ID point1Id;
    ID point2Id;
    for (Edges_Iterator it = mesh.edgeList.begin(); it < mesh.edgeList.end(); ++it)
    {
        point1Id = it->point (0).localId();
        point2Id = it->point (1).localId();
        existingEdges.insert (std::make_pair ( makeBareEdge ( point1Id, point2Id ).first, it->localId() ) );
    }


    if ( !buildBoundaryEdges && buildInternalEdges )
    {
        if ( mesh.storedEdges() < numBoundaryEdgesFound )
        {
            errorStream << "ERROR in buildEdges(): mesh does not contain boundary edges" << std::endl;
            errorStream << "  I need to set buildBoundaryEdges=true" << std::endl;
            errorStream << "  ABORT CONDITION" << std::endl;
            return false;
        }
        else if ( mesh.storedEdges() > numBoundaryEdgesFound )
        {
            mesh.edgeList.resize ( numBoundaryEdgesFound );
        }
    }
    mesh.setNumBEdges ( numBoundaryEdgesFound );
    mesh.setNumEdges ( numBoundaryEdgesFound + numInternalEdgesFound );

    if ( buildBoundaryEdges && ! buildInternalEdges )
    {
        mesh.setMaxNumEdges ( numBoundaryEdgesFound, false );
    }
    if ( buildInternalEdges )
    {
        mesh.setMaxNumEdges ( numBoundaryEdgesFound + numInternalEdgesFound, true );
    }

    if (verbose)
    {
        errorStream << "Building edges" << std::endl;
    }

    edge_Type edge;

    if ( buildBoundaryEdges )
    {

        if ( verbose )
        {
            logStream << "**** Marker Flags for Newly Created Boundary Edges ***"
                      << std::endl;
            logStream << "Edgeid->marker" << std::endl;
        }

        // First boundary.
        for ( temporaryEdgeContainer_Type::iterator edgeContainerIterator = temporaryEdgeContainer->begin();
                edgeContainerIterator != temporaryEdgeContainer->end(); ++edgeContainerIterator )
        {
            faceIdToLocalEdgeIdPair = edgeContainerIterator->second;
            faceId = faceIdToLocalEdgeIdPair.first; // Face ID
            facePtr = &mesh.face ( faceId ); // Face
            jEdgeLocalId = faceIdToLocalEdgeIdPair.second;       // The local ID of edge on face
            point1Id = facePtr->point ( faceShape_Type::edgeToPoint ( jEdgeLocalId, 0) ).localId();
            point2Id = facePtr->point ( faceShape_Type::edgeToPoint ( jEdgeLocalId, 1) ).localId();
            existingEdges_It = existingEdges.find ( ( makeBareEdge ( point1Id, point2Id ) ).first);
            if (existingEdges_It != existingEdges.end() )
            {
                edge = mesh.edge (existingEdges_It->second);
                edgeExists = true;

            }
            else
            {
                edge = edge_Type();
                edge.setId ( mesh.edgeList.size() );
                edgeExists = false;
            }
            for ( UInt kPointId = 0; kPointId < edge_Type::S_numPoints; ++kPointId )
            {
                edge.setPoint ( kPointId, facePtr->point ( faceShape_Type::edgeToPoint ( jEdgeLocalId, kPointId ) ) );
            }

            // Get marker value inheriting from points if necessary
            if ( edge.isMarkerUnset() )
            {
                inheritPointsWeakerMarker ( edge );
            }
            edge.setBoundary (true);
            if (edgeExists)
            {
                newEdgeId = edge.localId();
                mesh.setEdge (edge, newEdgeId);
            }
            else
            {
                newEdgeId = mesh.addEdge ( edge).localId();
            }

            if ( verbose )
            {
                if ( newEdgeId % 6 == 0 )
                {
                    logStream << std::endl;
                }
                logStream << newEdgeId << " -> " << edge.markerID();
                logStream << " ";
            }
        }

        if ( verbose )
            logStream << std::endl << "  *****  END OF LIST OF BOUNDARY EDGES ****"
                      << std::endl;

        mesh.setLinkSwitch ( std::string ( "HAS_BOUNDARY_RIDGES" ) );
    }

    if ( !externalContainerIsProvided )
    {
        delete temporaryEdgeContainer;
    }

    if ( !buildInternalEdges )
    {
        mesh.unsetLinkSwitch ( std::string ( "HAS_ALL_RIDGES" ) );
        return true;
    }



    // Now internal edges
    // free some memory
    volume_Type* volumePtr;
    for ( temporaryEdgeContainer_Type::iterator edgeContainerIterator = edgeContainer.begin();
            edgeContainerIterator != edgeContainer.end(); ++edgeContainerIterator )
    {
        faceIdToLocalEdgeIdPair = edgeContainerIterator->second;
        faceId = faceIdToLocalEdgeIdPair.first; // Volume ID
        volumePtr = &mesh.volume ( faceId ); // Volume that generated the edge
        jEdgeLocalId = faceIdToLocalEdgeIdPair.second;       // The local ID of edge on volume
        point1Id = volumePtr->point ( volumeShape_Type::edgeToPoint ( jEdgeLocalId, 0) ).localId();
        point2Id = volumePtr->point ( volumeShape_Type::edgeToPoint ( jEdgeLocalId, 1) ).localId();
        existingEdges_It = existingEdges.find ( ( makeBareEdge ( point1Id, point2Id ) ).first);
        if (existingEdges_It != existingEdges.end() )
        {
            edge = mesh.edge (existingEdges_It->second);
            edgeExists = true;

        }
        else
        {
            edge = edge_Type();
            edge.setId ( mesh.edgeList.size() );
            edgeExists = false;
        }
        for ( UInt kPointId = 0; kPointId < edge_Type::S_numPoints; ++kPointId )
        {
            edge.setPoint ( kPointId, volumePtr->point ( volumeShape_Type::edgeToPoint ( jEdgeLocalId, kPointId ) ) );
        }

        edge.setMarkerID (edge.nullMarkerID() ); // Set marker to null
        edge.setBoundary (false);
        if (edgeExists)
        {
            mesh.setEdge (edge, edge.localId() );
        }
        else
        {
            mesh.addEdge ( edge);
        }
    }

    mesh.setLinkSwitch ( std::string ( "HAS_ALL_RIDGES" ) );

    return true;
}


/*
 *******************************************************************************
 UTILITIES TO TRANSFORM A MESH
 *******************************************************************************
 */
//! It builds a P2 mesh from P1 data.
/*!
    @author L.Formaggia.
    @version Version 1.0
    @pre All compulsory structures in mesh must have been already set: volumes and boundary faces.
    @pre Points list MUST have been dimensioned correctly!!!
    @note the function takes advantage of the fact that vertex are stored first
    @param mesh[out] A mesh
    @param logStream[out] Log stream for information on the newly created markers for boundary edges
 */
template <typename MeshType>
void
p2MeshFromP1Data ( MeshType& mesh, std::ostream& logStream = std::cout )
{
    bool verbose ( mesh.comm()->MyPID() == 0 );

    typedef typename MeshType::elementShape_Type  GeoShape;
    typedef typename MeshType::facetShape_Type GeoBShape;
    ASSERT_PRE ( GeoShape::S_numPoints > 4, "p2MeshFromP1Data ERROR: we need a P2 mesh" );

    if ( verbose ) logStream << "Building P2 mesh points and connectivities from P1 data"
                                 << std::endl;


    typename MeshType::point_Type*        pointPtr = 0;
    typename MeshType::edge_Type*         edgePtr = 0;
    typename MeshType::volume_Type*      elementPtr = 0;
    typename MeshType::face_Type*       facePtr = 0;
    typedef typename MeshType::volumes_Type  elementContainer_Type;
    typedef typename MeshType::faces_Type faceContainer_Type;

    MeshElementBareHandler<BareEdge>           bareEdgeHandler;
    std::pair<UInt, bool>                edgeIdToBoolPair;
    UInt                                 point1Id, point2Id, edgeId;
    std::pair<BareEdge, bool>            bareEdgeToBoolPair;
    typename MeshType::elementShape_Type      elementShape;

    if ( verbose )
    {
        logStream << "Processing " << mesh.storedEdges() << " P1 Edges" << std::endl;
    }
    UInt numBoundaryEdges = mesh.numBEdges();
    for ( UInt jEdgeId = 0; jEdgeId < mesh.storedEdges(); ++jEdgeId )
    {
        edgePtr = & mesh.edge ( jEdgeId );
        point1Id = ( edgePtr->point ( 0 ) ).localId();
        point2Id = ( edgePtr->point ( 1 ) ).localId();
        pointPtr = & mesh.addPoint ( jEdgeId < numBoundaryEdges, false ); // true for boundary points
        pointPtr->setId ( mesh.pointList.size() - 1 );
        pointPtr->x() = ( ( edgePtr->point ( 0 ) ).x() +
                          ( edgePtr->point ( 1 ) ).x() ) * .5;
        pointPtr->y() = ( ( edgePtr->point ( 0 ) ).y() +
                          ( edgePtr->point ( 1 ) ).y() ) * .5;
        pointPtr->z() = ( ( edgePtr->point ( 0 ) ).z() +
                          ( edgePtr->point ( 1 ) ).z() ) * .5;

        /*
          If we have set a marker for the boundary edge, that marker is
          inherited by the new created point. Otherwise the edge (and the new
          created point) gets the WORST marker among the two end Vertices
         */
        if ( edgePtr->isMarkerUnset() )
        {
            inheritPointsWeakerMarker ( *edgePtr );
        }
        pointPtr->setMarkerID ( edgePtr->markerID() );
        // todo check that the localId() of the new point is correctly set
        edgePtr->setPoint ( 3, pointPtr ); //use overloaded version that takes a pointer
        bareEdgeToBoolPair = makeBareEdge ( point1Id, point2Id );
        edgeIdToBoolPair = bareEdgeHandler.addIfNotThere ( bareEdgeToBoolPair.first, pointPtr->localId() );
    }
    // Now the other edges, of which I do NOT build the global stuff
    // (I would need to check the switch but I will do that part later on)
    if ( GeoShape::S_nDimensions == 3 )
    {
        UInt numBoundaryFaces = mesh.numBFaces();

        if ( verbose ) logStream << "Processing " << mesh.storedFaces() << " Face Edges"
                                     << std::endl;
        for ( UInt kFaceId = 0; kFaceId < mesh.storedFaces(); ++kFaceId )
        {
            facePtr = &mesh.face ( kFaceId );
            for ( UInt jEdgeLocalId = 0; jEdgeLocalId < mesh.numLocalEdgesOfFace(); jEdgeLocalId++ )
            {
                point1Id = GeoBShape::edgeToPoint ( jEdgeLocalId, 0 );
                point2Id = GeoBShape::edgeToPoint ( jEdgeLocalId, 1 );
                point1Id = ( facePtr->point ( point1Id ) ).localId();
                point2Id = ( facePtr->point ( point2Id ) ).localId();
                bareEdgeToBoolPair = makeBareEdge ( point1Id, point2Id );
                edgeId = bareEdgeHandler.id ( bareEdgeToBoolPair.first );
                if ( edgeId != NotAnId )
                {
                    pointPtr = &mesh.point ( edgeId );
                }
                else
                {
                    // new edge -> new Point
                    pointPtr = &mesh.addPoint ( kFaceId < numBoundaryFaces, false ); // true for boundary points
                    edgeIdToBoolPair = bareEdgeHandler.addIfNotThere ( bareEdgeToBoolPair.first, pointPtr->localId() );
                    pointPtr->setId ( mesh.pointList.size() - 1 );
                    pointPtr->x() = ( mesh.point ( point1Id ).x() +
                                      mesh.point ( point2Id ).x() ) * .5;
                    pointPtr->y() = ( mesh.point ( point1Id ).y() +
                                      mesh.point ( point2Id ).y() ) * .5;
                    pointPtr->z() = ( mesh.point ( point1Id ).z() +
                                      mesh.point ( point2Id ).z() ) * .5;
                    // If we have set a marker for the face, that marker is
                    // inherited by the new created point
                    pointPtr->setMarkerID ( facePtr->markerID() );
                }
                facePtr->setPoint ( GeoBShape::S_numVertices + jEdgeLocalId, pointPtr );
            }
        }
    }

    if ( verbose ) logStream << "Processing " << mesh.numElements() << " Mesh Elements"
                                 << std::endl;
    UInt nev = GeoShape::S_numVertices;
    for ( UInt kElementId = 0; kElementId < mesh.numElements(); ++kElementId )
    {
        elementPtr = &mesh.element ( kElementId );
        for ( UInt jEdgeLocalId = 0; jEdgeLocalId < mesh.numLocalEdges(); jEdgeLocalId++ )
        {
            point1Id = elementShape.edgeToPoint ( jEdgeLocalId, 0 );
            point2Id = elementShape.edgeToPoint ( jEdgeLocalId, 1 );
            point1Id = ( elementPtr->point ( point1Id ) ).localId();
            point2Id = ( elementPtr->point ( point2Id ) ).localId();
            bareEdgeToBoolPair = makeBareEdge ( point1Id, point2Id );
            edgeId = bareEdgeHandler.id ( bareEdgeToBoolPair.first );
            if ( edgeId != NotAnId )
            {
                pointPtr = &mesh.point ( edgeId );
            }
            else
            {
                // cannot be on boundary if the mesh is proper!
                pointPtr = &mesh.addPoint ( false, false );
                edgeIdToBoolPair = bareEdgeHandler.addIfNotThere ( bareEdgeToBoolPair.first, pointPtr->localId() );
                pointPtr->setId ( mesh.pointList.size() - 1 );
                pointPtr->x() = ( mesh.point ( point1Id ).x() +
                                  mesh.point ( point2Id ).x() ) * .5;
                pointPtr->y() = ( mesh.point ( point1Id ).y() +
                                  mesh.point ( point2Id ).y() ) * .5;
                pointPtr->z() = ( mesh.point ( point1Id ).z() +
                                  mesh.point ( point2Id ).z() ) * .5;
                pointPtr->setMarkerID ( edgePtr->markerID() );
            }
            elementPtr->setPoint ( nev + jEdgeLocalId, pointPtr );
        }
    }
    /*=============================*/
    if ( verbose ) logStream << " ******* Done Construction of P2 Mesh *******"
                                 << std::endl << std::endl;
}

/*! Fix mesh switches
 * Using some heuristics it tries to fix mesh switches
 */

// template<typename MeshType>
// void
// fixSwitches(MeshType ^ mesh, std::ostream & logStream=std::cout, bool verbose=false)
// {

//   logStream<<" ************** FIXING MESH SWITCHES **********************"<<std::endl;
//   logStream<<"            Mesh switches Status before fixing"<<std::endl;
//   mesh.showLinkSwitch(verbose,logStream);
//   if (mesh.storedFaces()> mesh.numBFaces()){
//     mesh.setLinkSwitch("HAS_ALL_FACES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_ALL_FACES");
//   }
//   if (mesh.numBFaces>0){
//     mesh.setLinkSwitch("HAS_BOUNDARY_FACES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDARY_FACES");
//   }
//   if (mesh.storedEdges()> mesh.numBEdges()){
//     mesh.setLinkSwitch("HAS_ALL_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_ALL_EDGES");
//   }
//   if (mesh.numBEdges()> 0){
//     mesh.setLinkSwitch("HAS_BOUNDARY_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDAY_EDGES");
//   }
//   if (mesh.storedFaces()>0){
//     if(mesh.face(0).firstAdjacentElementIdentity(
//     mesh.setLinkSwitch("HAS_BOUNDARY_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDAY_EDGES");
//   }

/** Class to transform a mesh.
 * A class that implements methods to transform a mesh without changing
 * mesh connectivities. It has a constructor that takes the mesh to be transformed
 * @note The Template RMTYPE is used to compile with IBM AIX compilers
 * @author Luca Formaggia
 * @date 2 August 2011
 */
//
template <typename REGIONMESH, typename RMTYPE = typename REGIONMESH::markerCommon_Type >
class MeshTransformer
{
public:
    /** the constructor may take a reference to the mesh to be manipulated */
    MeshTransformer (REGIONMESH& m);
    /** Move the mesh according to a given displacement.
    *
    *  It moves the mesh from the last position saved with savePoints()
    *  For backward compatibility, if it is called before without calling
    *  savePoints(), the first time it is called it will save the current mesh point and then
    *  apply the movement.
    *
    *  Displacement is a 3*numpoints() VECTOR which stores the x-displacement first,
    *  then the y-displacements etc.
    *
    *  The VECTOR object must comply with lifeV distributed vector concept EpetraVector
    *  in particular it must have the methods isGlobalIDPresent(Uint i).
    *
    *  @author Miguel Fernandez
    *  @date 11/2002
    *
    *  @param disp Displacement vector. In this version it must be an EpetraVector
    *  @param dim  Length of vector disp.
    */
    template <typename VECTOR>
    void moveMesh ( const VECTOR& disp, UInt dim);
    /** Transform the mesh. It uses  boost::numeric::ublas (3,3) matrices
     *  scale, rotate and translate to perform the mesh movement
     *  (operations performed in this order).
     *  @date   14/09/2009
     *  @author Cristiano Malossi
     *  @note - Rotation follows Paraview conventions: first rotate around z-axis,
     *          then around y-axis and finally around x-axis;
     *  @note - All the vectors must allow the command: operator[];
     *  @param scale        vector of three components for (x,y,z) scaling of the mesh
     *  @param rotate       vector of three components (radiants) for rotating the mesh
     *  @param translate    vector of three components for (x,y,z) translation the mesh
     *
     */
    template <typename VECTOR>
    void transformMesh ( const VECTOR& scale, const VECTOR& rotate, const VECTOR& translate );

    /** Transform the mesh according to a given mapping.
     *  Transform the mesh according to a given meshMapping(Real& x, Real& y, Real& z).
     *  @date   12/2010
     *  @author Mauro Perego
     *  @param meshMapping   function void meshMmapping(Real& x, Real& y, Real& z) which receive
     *                   x, y, z, and transform them according to a certain mapping
     */
    template <typename function>
    void transformMesh ( const function& meshMapping);

    //! Tells if we store old points
    /**
     * If true than we can interrogate the old point position
     */
    bool hasOldPoint() const
    {
        return ! (this->M_pointList.empty() );
    }
    //! Saves the mesh points
    /**
     * Useful for algorithms which require to remember the position of the mesh
     * before the movement
     */
    void savePoints();
    /** Resets movement. Next step is like the mesh has never moved
     */
    void resetMovement()
    {
        this->M_pointList.clear();
    }

    /** Returns the i-th mesh Point before the last movement.

     *
     *  If the mesh points have not been saved with a previous call to
     *  savePoints() the method returns the mesh point
     *
     *  @note Avoid extensive use: it is inefficient. use pointListInitial() instead
     *  @param i Id of the Point.
     *  @return i-th mesh Point before the last movement.
     */
    typename REGIONMESH::point_Type const& pointInitial ( ID const i ) const;
    /** Returns a constant reference to the list of Points before the last movement.
      *
      *  If the mesh points have not been saved with a previous call to
      *  savePoints() returns the current mesh Points
      *
      *  @return The list mesh Point before the last movement.
      */
    typename REGIONMESH::points_Type const& pointListInitial() const;
private:
    /** Appropriately sets internal switches
     *
     *  It must be called by any mesh transformation method
     *  to ensure that the handling of (possibly) stored points
     *  works;
     */
    REGIONMESH& M_mesh;
    typename REGIONMESH::points_Type M_pointList;
};
/** Mesh statistics.
 *  Namespace that groups functions which operate on a mesh to extract statistics.
 *  The functions do not modify mesh content
 *  @author Luca Formaggia
 *  @date 3 August 2011
 */
namespace MeshStatistics
{

/** It holds statistics on mesh size.
 *  Mesh spacings:
 *  meshSize.minH  Min h
 *  meshsize.maxH  Max h
 *  meshsize.meanH Average h
 */
struct meshSize
{
    Real maxH;
    Real minH;
    Real meanH;
};
//! Computes mesh sizes
template<typename REGIONMESH>
meshSize computeSize (const REGIONMESH&);
}// namespace MeshStatistics

// *****   IMPLEMENTATIONS ****
// The Template RMTYPE is used to compile with IBM compilers
template <typename REGIONMESH, typename RMTYPE >
MeshTransformer<REGIONMESH, RMTYPE >::MeshTransformer (REGIONMESH& m) : M_mesh (m), M_pointList() {}
/**
 * @todo this method should be changed to make sure not to generate invalid elements
 */
template <typename REGIONMESH, typename RMTYPE >
template <typename VECTOR>
void MeshTransformer<REGIONMESH, RMTYPE >::moveMesh ( const VECTOR& disp, UInt dim )
{
    // the method must be called with a Repeated vector
    if ( disp.mapType() == Unique )
    {
#ifdef HAVE_LIFEV_DEBUG
        std::cerr << "Info: moveMesh() requires a Repeated vector, a copy of the passed Unique vector will be created.\n"
                  << "To optimize your code, you should pass a repeated vector to avoid the conversion." << std::endl;
#endif
        this->moveMesh ( VECTOR ( disp, Repeated ), dim );
        return;
    }

    if ( !this->hasOldPoint() )
    {
        this->savePoints();
    }

    typedef typename REGIONMESH::points_Type points_Type;
    points_Type& pointList ( M_mesh.pointList );
    for ( UInt i = 0; i < M_mesh.pointList.size(); ++i )
    {
        for ( UInt j = 0; j < nDimensions; ++j )
        {
            int globalId = pointList[i].id();
            ASSERT ( disp.isGlobalIDPresent ( globalId + dim * j ), "global ID missing" );
            pointList[ i ].coordinate ( j ) = M_pointList[ i ].coordinate ( j ) + disp[ j * dim + globalId ];
        }
    }
}

template<typename REGIONMESH, typename RMTYPE >
void MeshTransformer<REGIONMESH, RMTYPE >::savePoints()
{
    if (M_pointList.capacity() < M_mesh.pointList.size() )
    {
        // Create space and add
        M_pointList.clear();
        M_pointList.reserve ( M_mesh.numPoints() );
        std::copy (M_mesh.pointList.begin(), M_mesh.pointList.end(),
                   std::back_inserter (M_pointList) );
    }
    else
    {
        // Overwrite
        std::copy (M_mesh.pointList.begin(), M_mesh.pointList.end(), M_pointList.begin() );
    }
}
//  The Template RMTYPE is used to compile with IBM compilers
template <typename REGIONMESH, typename RMTYPE >
const typename REGIONMESH::point_Type&
MeshTransformer<REGIONMESH, RMTYPE >::pointInitial ( ID const i ) const
{
    ASSERT_BD ( i < M_mesh.pointList.size() );
    return M_pointList.empty() ? M_mesh.pointList[i] : this->M_pointList[i];
}

template <typename REGIONMESH, typename RMTYPE >
const typename REGIONMESH::points_Type&
MeshTransformer<REGIONMESH, RMTYPE >::pointListInitial() const
{
    return M_pointList.empty() ? M_mesh.points_Type : M_pointList;
}
//  The Template RMTYPE is used to compile with IBM compilers
//! @todo Change using homogeneous coordinates to make it more efficient.
template <typename REGIONMESH, typename RMTYPE >
template <typename VECTOR>
void MeshTransformer<REGIONMESH, RMTYPE >::transformMesh ( const VECTOR& scale, const VECTOR& rotate, const VECTOR& translate )
{
    // Make life easier
    typename REGIONMESH::points_Type& pointList (M_mesh.pointList);

    //Create the 3 planar rotation matrix and the scale matrix
    boost::numeric::ublas::matrix<Real> R (3, 3), R1 (3, 3), R2 (3, 3), R3 (3, 3), S (3, 3);

    R1 (0, 0) =  1.;
    R1 (0, 1) =  0.;
    R1 (0, 2) =  0.;
    R1 (1, 0) =  0.;
    R1 (1, 1) =  std::cos (rotate[0]);
    R1 (1, 2) = -std::sin (rotate[0]);
    R1 (2, 0) =  0.;
    R1 (2, 1) =  std::sin (rotate[0]);
    R1 (2, 2) =  std::cos (rotate[0]);

    R2 (0, 0) =  std::cos (rotate[1]);
    R2 (0, 1) =  0.;
    R2 (0, 2) =  std::sin (rotate[1]);
    R2 (1, 0) =  0.;
    R2 (1, 1) =  1.;
    R2 (1, 2) = 0.;
    R2 (2, 0) = -std::sin (rotate[1]);
    R2 (2, 1) =  0.;
    R2 (2, 2) =  std::cos (rotate[1]);

    R3 (0, 0) =  std::cos (rotate[2]);
    R3 (0, 1) = -std::sin (rotate[2]);
    R3 (0, 2) = 0.;
    R3 (1, 0) =  std::sin (rotate[2]);
    R3 (1, 1) =  std::cos (rotate[2]);
    R3 (1, 2) = 0.;
    R3 (2, 0) =  0;
    R3 (2, 1) =  0.;
    R3 (2, 2) = 1.;

    S (0, 0) = scale[0];
    S (0, 1) = 0.;
    S (0, 2) = 0.;
    S (1, 0) = 0.;
    S (1, 1) = scale[1];
    S (1, 2) = 0.;
    S (2, 0) = 0.;
    S (2, 1) = 0.;
    S (2, 2) = scale[2];

    //The total rotation is: R = R1*R2*R3 (as in Paraview we rotate first around z, then around y, and finally around x).
    //We also post-multiply by S to apply the scale before the rotation.
    R = prod ( R3, S );
    R = prod ( R2, R );
    R = prod ( R1, R );

    //Create the 3D translate vector
    boost::numeric::ublas::vector<Real> P (3), T (3);
    T (0) = translate[0];
    T (1) = translate[1];
    T (2) = translate[2];

    //Apply the transformation
    for ( UInt i (0); i < pointList.size(); ++i )
    {
        //P = pointList[ i ].coordinate(); // Try to avoid double copy if possible

        P ( 0 ) = pointList[ i ].coordinate ( 0 );
        P ( 1 ) = pointList[ i ].coordinate ( 1 );
        P ( 2 ) = pointList[ i ].coordinate ( 2 );

        P = T + prod ( R, P );

        pointList[ i ].coordinate ( 0 ) = P ( 0 );
        pointList[ i ].coordinate ( 1 ) = P ( 1 );
        pointList[ i ].coordinate ( 2 ) = P ( 2 );
    }
}
//  The Template RMTYPE is used to compile with IBM compilers
template <typename REGIONMESH, typename RMTYPE >
template <typename function>
void MeshTransformer<REGIONMESH, RMTYPE >::transformMesh ( const function& meshMapping)
{
    // Make life easier
    typename REGIONMESH::points_Type& pointList (M_mesh.pointList);

    for ( UInt i = 0; i < pointList.size(); ++i )
    {
        typename REGIONMESH::point_Type& p = pointList[ i ];
        meshMapping (p.coordinate (0), p.coordinate (1), p.coordinate (2) );
    }
}

template <typename REGIONMESH>
MeshStatistics::meshSize MeshStatistics::computeSize (REGIONMESH const& mesh)
{
    const Real bignumber = std::numeric_limits<Real>::max();
    Real MaxH (0), MinH (bignumber), MeanH (0);
    Real deltaX (0), deltaY (0), deltaZ (0), sum (0);
    typedef typename REGIONMESH::edges_Type edges_Type;
    edges_Type const& edgeList (mesh.edgeList);

    ASSERT0 (edgeList.size() > 0, "computeSize requires edges!");

    for (typename edges_Type::const_iterator i = edgeList.begin(); i < edgeList.end() ; ++i )
    {
        deltaX =  i->point ( 1 ).x() - i->point ( 0 ).x();
        deltaY =  i->point ( 1 ).y() - i->point ( 0 ).y();
        deltaZ =  i->point ( 1 ).z() - i->point ( 0 ).z();

        deltaX *= deltaX;
        deltaY *= deltaY;
        deltaZ *= deltaZ;
        sum     = deltaX + deltaY + deltaZ;
        MaxH = std::max ( MaxH, sum);
        MinH = std::min ( MinH, sum);
        MeanH += sum;
    }

    MeshStatistics::meshSize out;
    out.minH  = std::sqrt ( MinH );
    out.meanH = std::sqrt ( MeanH / static_cast<Real> ( edgeList.size() ) );
    out.maxH  = std::sqrt ( MaxH );
    return out;
}

template <typename RegionMeshType, typename RegionFunctorType>
void assignRegionMarkerID ( RegionMeshType& mesh, const RegionFunctorType& fun )
{

    // Extract the element list.
    typename RegionMeshType::elements_Type& elementList = mesh.elementList ();
    const UInt elementListSize = elementList.size();

    // Loop on the elements and decide the flag.
    for ( UInt i = 0; i < elementListSize; ++i )
    {
        // Computes the barycentre of the element
        Vector3D barycentre;

        for ( UInt k = 0; k < RegionMeshType::element_Type::S_numPoints; k++ )
        {
            barycentre += elementList[i].point ( k ).coordinates();
        }
        barycentre /= RegionMeshType::element_Type::S_numPoints;

        // Set the marker Id.
        elementList[i].setMarkerID ( fun ( barycentre ) );

    }

} // assignRegionMarkerID

} // namespace MeshUtility

} // namespace LifeV

#endif /* MESHUTILITY_H */

