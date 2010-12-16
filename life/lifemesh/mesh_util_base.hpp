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

#ifndef MESHUTILITYBASE_H
#define MESHUTILITYBASE_H 1

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include <life/lifecore/life.hpp>
#include <life/lifecore/switch.hpp>
#include <life/lifemesh/bareItems.hpp>
#include <life/lifemesh/markers.hpp>
#include <life/lifemesh/meshEntity.hpp>

namespace LifeV
{

// namespace MeshUtilities
// {
/*
   todo change the class names:
       "EnquireBEntity --> EnquireBoundaryEntity" and so on
       "GetCoordComponent --> GetCoordinatesComponent" and so on
 */


//! A locally used structure, not meant for general use
typedef std::map<BareFace, std::pair<ID, ID >,
cmpBareItem<BareFace> > temporaryFaceContainer_Type;
// deprecated
typedef temporaryFaceContainer_Type __attribute__ (( deprecated )) TempFaceContainer;

//! A locally used structure, not meant for general use
typedef std::map<BareEdge, std::pair<ID, ID>,
cmpBareItem<BareEdge> > temporaryEdgeContainer_Type;
// deprecated
typedef temporaryEdgeContainer_Type __attribute__ (( deprecated )) TempEdgeContainer;

/*
*******************************************************************************
                            FUNCTORS
*******************************************************************************
*/
//! @defgroup Test_Functors Some useful functors to be used to test mesh entities

//! Functor to check if a Point, Face or Edge is on the boundary
/*!
    @ingroup Test_Functors

    A geometric entity is on the boundary if all its vertices are boundary vertices.

    @pre It assumes that boundary points in the mesh are correctly set.
    @pre The MeshType must export the typenames PointType, FaceType and EdgeType.
*/
template <typename MeshType>
class EnquireBEntity
{
public:

    //! @name Public Types
    //@{
    // todo use boost::shared_ptr
    typedef MeshType                      mesh_Type;
    typedef mesh_Type const *             meshPtr_Type;
    typedef typename mesh_Type::FaceType  face_Type;
    typedef typename mesh_Type::EdgeType  edge_Type;
    typedef typename mesh_Type::PointType point_Type;

    // The following should be removed
    // typedef typename mesh_Type::FaceType  FaceType;
    // typedef typename mesh_Type::EdgeType  EdgeType;
    // typedef typename mesh_Type::PointType PointType;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructor taking a mesh object
    /*!
        @param mesh a mesh object
     */
    explicit EnquireBEntity( mesh_Type const & mesh ) : meshPtr( &mesh )
    {}

    //! Copy Constructor
    EnquireBEntity( EnquireBEntity const & enquireBoundaryEntity ) :
            meshPtr( enquireBoundaryEntity.meshPtr )
    {}

    //! Virtual Destructor
    virtual ~EnquireBEntity()
    {}
    //@}

    //! @name Operators
    //@{

    //! The function call operator
    /*!
        @param face a face entity in the mesh_Type
        @return true if the face is on the boundary, false otherwise
     */
    bool operator() ( const face_Type & face ) const
    {
        bool isBoundary = true;
        for ( UInt kPointId = 1; kPointId <= face_Type::numVertices; ++kPointId )
        {
            isBoundary = isBoundary & face.point( kPointId ).boundary();
        }
        return isBoundary;
    }

    //! The function call operator
    /*!
        @param edge an edge entity in the mesh_Type
        @return true if the edge is on the boundary, false otherwise
     */
    bool operator() ( const edge_Type & edge ) const
    {
        bool isBoundary = true;
        for ( UInt kPointId = 1; kPointId <= edge_Type::numVertices; ++kPointId )
        {
            isBoundary = isBoundary & edge.point( kPointId ).boundary();
        }
        return isBoundary;
    }

    //! The function call operator
    /*!
        @param point a point entity in the mesh_Type
        @return true if the point is on the boundary, false otherwise
     */
    inline bool operator() ( const point_Type & point ) const
    {
        return point.boundary();
    }

    //@}

private:
    //! @name Private Types
    //@{

    //! Empty Constructor
    EnquireBEntity() {}

    //@}

    meshPtr_Type meshPtr;
};


//! Functor to check if a Face is on the boundary
/*!
    @ingroup Test_Functors

    This object uses the information contained in a FaceContainer produced
    by findBoundaryFaces(). It does not use the information contained in the
    mesh PointList, so it differs from EnquireBEntity.

    @pre boundaryFaceContainer has been previously set by a call to findBoundaryFaces()
*/
template <typename MeshType>
class EnquireBFace
{
public:

    //! @name Public Types
    //@{
    typedef MeshType                            mesh_Type;
    typedef mesh_Type const *                   meshPtr_Type;
    typedef typename mesh_Type::FaceType        face_Type;
    typedef typename mesh_Type::FaceShape       faceShape_Type;
    typedef temporaryFaceContainer_Type const * temporaryFaceContainerPtr_Type;
    // The following should be removed
    // typedef typename mesh_Type::FaceType  FaceType;
    // typedef typename mesh_Type::FaceShape FaceShape;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructor taking a mesh object and a face container
    /*!
        @param mesh a mesh object
        @param boundaryFaceContainer a container of boundary faces
        @pre boundaryFaceContainer has been previously set by a call to findBoundaryFaces()
     */
    EnquireBFace( mesh_Type const & mesh, temporaryFaceContainer_Type const & boundaryFaceContainer ) :
            meshPtr                 ( &mesh ),
            boundaryFaceContainerPtr( &boundaryFaceContainer )
    {}

    //! Copy Constructor
    EnquireBFace( EnquireBFace const & enquireBoundaryFace ) :
            meshPtr                 ( enquireBoundaryFace.meshPtr ),
            boundaryFaceContainerPtr( enquireBoundaryFace.boundaryFaceContainerPtr )
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
    bool operator() ( const face_Type & face ) const
    {
        ID point1Id, point2Id, point3Id, point4Id;
        BareFace bareFace;

        point1Id = face.point( 1 ).id();
        point2Id = face.point( 2 ).id();
        point3Id = face.point( 3 ).id();
        if ( faceShape_Type::numVertices == 4 )
        {
            point4Id = face.point( 4 ).id();
            bareFace = ( makeBareFace( point1Id, point2Id, point3Id, point4Id ) ).first;
        }
        else
        {
            bareFace = ( makeBareFace( point1Id, point2Id, point3Id ) ).first;
        }
        return ( boundaryFaceContainerPtr->find( bareFace ) != boundaryFaceContainerPtr->end() );
    }

    //@}
private:
    //! @name Private Methods
    //@{

    //! Empty Constructor
    EnquireBFace()
    {}
    //@}

    meshPtr_Type meshPtr;
    temporaryFaceContainerPtr_Type boundaryFaceContainerPtr;
};


/*! Functor to check if an edge is on the boundary

    @ingroup Test_Functors

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
    typedef MeshType                            mesh_Type;
    typedef mesh_Type const *                   meshPtr_Type;
    typedef typename mesh_Type::EdgeType        edge_Type;
    typedef typename mesh_Type::EdgeShape       edgeShape_Type;
    typedef temporaryEdgeContainer_Type const * temporaryEdgeContainerPtr_Type;
    // The following should be removed
    typedef typename mesh_Type::EdgeType  EdgeType;
    typedef typename mesh_Type::EdgeShape EdgeShape;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructor taking a mesh object and an edge container
    /*!
        @param mesh a mesh object
        @param boundaryEdgeContainer a container of boundary edges
     */
    EnquireBEdge( mesh_Type const & mesh, temporaryEdgeContainer_Type const & boundaryEdgeContainer ) :
            meshPtr                 ( &mesh ),
            boundaryEdgeContainerPtr( &boundaryEdgeContainer )
    {}

    //! Copy Constructor
    EnquireBEdge( EnquireBEdge const & enquireBoundaryEdge ) :
            meshPtr                 ( enquireBoundaryEdge.meshPtr ),
            boundaryEdgeContainerPtr( enquireBoundaryEdge.boundaryEdgeContainerPtr )
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
    bool operator() ( const edge_Type & edge ) const
    {
        ID point1Id, point2Id;
        BareEdge bareEdge;

        point1Id = edge.point( 1 ).id();
        point2Id = edge.point( 2 ).id();
        bareEdge = ( makeBareEdge( point1Id, point2Id ) ).first;
        return boundaryEdgeContainerPtr->find( bareEdge ) != boundaryEdgeContainerPtr->end();
    }

    //@}
private:
    //! @name Private Methods
    //@{

    //! Empty Constructor
    EnquireBEdge()
    {}
    //@}

    meshPtr_Type                   meshPtr;
    temporaryEdgeContainerPtr_Type boundaryEdgeContainerPtr;
};


/*! Functor to check if a mesh entity with boundary indicator is on the boundary

    @ingroup Test_Functors

    This objects works on mesh entities with boundary indicator (for instance a GeoPoint)
    by enquiring its boundary flag.

    @warning It assumes that boundary indicators are correctly set.
*/
template <typename MeshType>
class EnquireBPoint
{
public:

    //! @name Public Types
    //@{
    typedef MeshType                      mesh_Type;
    typedef mesh_Type *                   meshPtr_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructor taking a mesh object
    /*!
        @param mesh a mesh object
     */
    explicit EnquireBPoint( mesh_Type & mesh ) : meshPtr( &mesh )
    {}

    //! Copy Constructor
    EnquireBPoint( EnquireBPoint const & enquireBoundaryPoint ) :
            meshPtr( enquireBoundaryPoint.meshPtr )
    {}

    //! Virtual Destructor
    virtual ~EnquireBPoint()
    {}
    //@}

    //! @name Operators
    //@{

    //! The function call operator
    /*!
        @param meshEntityWithBoundary a mesh entity with boundary indicator
        @return true if the entity is on the boundary, false otherwise
     */
    bool operator() ( const MeshEntityWithBoundary & meshEntityWithBoundary ) const
    {
        return meshEntityWithBoundary.boundary();
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
    explicit GetCoordComponent( Int i );

    //! Copy constructor
    GetCoordComponent( const GetCoordComponent& getCoordComponent ) :
            componentIndex( getCoordComponent.componentIndex )
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
    GetOnes( const GetOnes& /* getOnes */ )
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
UInt findFaces( const MeshType & mesh, temporaryFaceContainer_Type & boundaryFaceContainer,
                UInt & numInternalFaces, temporaryFaceContainer_Type & internalFaces,
                bool buildAllFaces = false )
{
    UInt                                  point1Id, point2Id, point3Id, point4Id;
    BareFace                              bareFace;
    typename MeshType::VolumeShape        volumeShape;
    typedef typename MeshType::Volumes    volumeContainer_Type;
    temporaryFaceContainer_Type::iterator faceContainerIterator;

    // clean first in case it has been already used
    boundaryFaceContainer.clear();
    if ( buildAllFaces )
        internalFaces.clear();
    numInternalFaces = 0;

    for ( typename volumeContainer_Type::const_iterator volumeContainerIterator = mesh.volumeList.begin();
            volumeContainerIterator != mesh.volumeList.end(); ++volumeContainerIterator )
    {
        for ( ID jFaceLocalId = 1; jFaceLocalId <= mesh.numLocalFaces(); ++jFaceLocalId )
        {
            point1Id = volumeShape.fToP( jFaceLocalId, 1 );
            point2Id = volumeShape.fToP( jFaceLocalId, 2 );
            point3Id = volumeShape.fToP( jFaceLocalId, 3 );
            // go to global
            point1Id = ( volumeContainerIterator->point( point1Id ) ).id();
            point2Id = ( volumeContainerIterator->point( point2Id ) ).id();
            point3Id = ( volumeContainerIterator->point( point3Id ) ).id();
            if ( MeshType::FaceShape::numVertices == 4 )
            {
                point4Id = volumeShape.fToP( jFaceLocalId, 4 );
                point4Id = ( volumeContainerIterator->point( point4Id ) ).id();
                bareFace = ( makeBareFace( point1Id, point2Id, point3Id, point4Id ) ).first;
            }
            else
            {
                bareFace = ( makeBareFace( point1Id, point2Id, point3Id ) ).first;
            }

            if ( ( faceContainerIterator = boundaryFaceContainer.find( bareFace ) ) == boundaryFaceContainer.end() )
            {
                boundaryFaceContainer.insert(
                        std::make_pair( bareFace, std::make_pair( volumeContainerIterator->id(), jFaceLocalId ) ) );
            }
            else
            {
                if ( buildAllFaces && point1Id > point2Id )
                {
                    internalFaces.insert(
                            ( std::make_pair( bareFace, std::make_pair( volumeContainerIterator->id(), jFaceLocalId ) ) ) );
                }
                boundaryFaceContainer.erase( faceContainerIterator ); // counted twice: internal face
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
UInt findBoundaryFaces( const MeshType & mesh,
                        temporaryFaceContainer_Type & boundaryFaceContainer,
                        UInt & numInternalFaces )
{
    temporaryFaceContainer_Type dummy;
    return findFaces( mesh, boundaryFaceContainer, numInternalFaces, dummy, false );
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
UInt findBoundaryEdges( const MeshType & mesh, temporaryEdgeContainer_Type & boundaryEdgeContainer )
{
    UInt                                 point1Id, point2Id;
    BareEdge                             bareEdge;
    typedef typename MeshType::FaceShape faceShape_Type;
    typedef typename MeshType::Faces     faceContainer_Type;


    if ( ! mesh.hasFaces() )
        return 0;

    // clean first in case it has been already used
    boundaryEdgeContainer.clear();

    // the following cycle assumes to visit only the boundary faces in mesh.faceList()
    for ( typename faceContainer_Type::const_iterator faceContainerIterator = mesh.faceList.begin();
            faceContainerIterator != mesh.faceList.begin() + mesh.numBFaces(); ++faceContainerIterator )
    {
        for ( ID jEdgeLocalId = 1; jEdgeLocalId <= mesh.numLocalEdgesOfFace(); ++jEdgeLocalId )
        {
            point1Id = faceShape_Type::eToP( jEdgeLocalId, 1 );
            point2Id = faceShape_Type::eToP( jEdgeLocalId, 2 );
            // go to global
            point1Id = ( faceContainerIterator->point( point1Id ) ).id();
            point2Id = ( faceContainerIterator->point( point2Id ) ).id();
            bareEdge = ( makeBareEdge( point1Id, point2Id ) ).first;
            boundaryEdgeContainer.insert(
                    std::make_pair( bareEdge, std::make_pair( faceContainerIterator->id(), jEdgeLocalId ) ) );
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
UInt findInternalEdges( const MeshType & mesh,
                        const temporaryEdgeContainer_Type & boundaryEdgeContainer,
                        temporaryEdgeContainer_Type & internalEdgeContainer )
{
    UInt                                   point1Id, point2Id;
    BareEdge                               bareEdge;
    typedef typename MeshType::VolumeShape volumeShape_Type;
    typedef typename MeshType::Volumes     volumeContainer_Type;
    temporaryEdgeContainer_Type            temporaryEdgeContainer;

    ASSERT0( mesh.numVolumes() > 0, "We must have some 3D elements stored n the mesh to use this function!" );

    internalEdgeContainer.clear();
    internalEdgeContainer.swap(temporaryEdgeContainer);

    for ( typename volumeContainer_Type::const_iterator volumeContainerIterator = mesh.volumeList.begin();
            volumeContainerIterator != mesh.volumeList.end(); ++volumeContainerIterator )
    {
        for ( ID jEdgeLocalId = 1; jEdgeLocalId <= mesh.numLocalEdges(); ++jEdgeLocalId )
        {
            point1Id = volumeShape_Type::eToP( jEdgeLocalId, 1 );
            point2Id = volumeShape_Type::eToP( jEdgeLocalId, 2 );
            // go to global
            point1Id = ( volumeContainerIterator->point( point1Id ) ).id();
            point2Id = ( volumeContainerIterator->point( point2Id ) ).id();
            bareEdge = ( makeBareEdge( point1Id, point2Id ) ).first;
            if ( boundaryEdgeContainer.find( bareEdge ) == boundaryEdgeContainer.end() )
                internalEdgeContainer.insert
                ( std::make_pair( bareEdge, std::make_pair( volumeContainerIterator->id(), jEdgeLocalId ) ) );
        }
    }
    return internalEdgeContainer.size();
}


/*
*******************************************************************************
                            MARKERS HANDLERS
*******************************************************************************
*/
//! @defgroup marker_handlers Used to manage missing handlers

/*! Sets the marker flag of a GeoElement of dimension greater than one

    @ingroup marker_handlers

    It gets the stronger marker of the GeoElement points. The marker
    hierarchy is defined in the markers.hpp file. It returns the new
    flag for the GeoElement. If any of the vertices has an unset marker
    the result is an unset flag for the GeoElement.

    @sa markers.hpp
    @warning It overrides the original marker flag.
    @return the new flag for geoElement
*/
template <typename GeoElementType>
EntityFlag __attribute__ (( deprecated ))
inheritStrongerMarker( GeoElementType & geoElement )
{
	return inheritPointsStrongerMarker( geoElement );
}
template <typename GeoElementType>
EntityFlag inheritPointsStrongerMarker( GeoElementType & geoElement )
{
    ASSERT_PRE( GeoElementType::nDim > 0,
                "A GeoElement with ndim<1 cannot inherit marker flags" );

    geoElement.setMarker( geoElement.point( 1 ).marker() );
    for ( ID jPointId = 2; jPointId <= GeoElementType::numVertices; ++jPointId )
        geoElement.setStrongerMarker( geoElement.point( jPointId ).marker() );
    return geoElement.marker();

}


/*! @ingroup marker_handlers

//! @brief Sets the marker flag of a GeoElement of dimension greater one

    It gets the weaker marker of the GeoElement points. The marker
    hierarchy is defined in the markers.hpp file. It returns the new
    flag for the GeoElement. If any of the vertices has an unset marker
    the result is an unset flag for the GeoElement.

    @sa markers.hpp
    @warning It overrides the original marker flag.
    @return the new flag for geoElement
*/
template <typename GeoElementType>
EntityFlag __attribute__ (( deprecated ))
inheritWeakerMarker( GeoElementType & geoElement )
{
	return inheritPointsWeakerMarker( geoElement );
}
template <typename GeoElementType>
EntityFlag inheritPointsWeakerMarker( GeoElementType & geoElement )
{
    ASSERT_PRE( GeoElementType::nDim > 0,
                "A GeoElement with ndim<1 cannot inherit marker flags" );

    geoElement.setMarker( geoElement.point( 1 ).marker() );
    for ( ID jPointId = 2; jPointId <= GeoElementType::numVertices; ++jPointId )
        geoElement.setWeakerMarker( geoElement.point( jPointId ).marker() );
    return geoElement.marker();

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
UInt __attribute__ (( deprecated ))
testClosedDomain_Top( MeshType const & mesh, UInt & numBoundaryEdges )
{
	return testClosedDomain( mesh, numBoundaryEdges );
}
template <typename MeshType>
UInt testClosedDomain( MeshType const & mesh, UInt & numBoundaryEdges )
{

    typedef std::set <BareEdge, cmpBareItem<BareEdge> > localTemporaryEdgeContainer_Type;
    localTemporaryEdgeContainer_Type                    localTemporaryEdgeContainer;
    UInt                                                point1Id, point2Id;
    BareEdge                                            bareEdge;
    typename MeshType::BElementShape                    facetShape;
    typedef typename MeshType::Faces                    faceContainer_Type;
    typedef typename MeshType::FaceType                 face_Type;
    localTemporaryEdgeContainer_Type::iterator          edgeContainerIterator;

    typename faceContainer_Type::const_iterator faceContainerIterator = mesh.faceList.begin();

    for ( UInt kFaceId = 0; kFaceId < mesh.numBFaces(); ++kFaceId )
    {
        std::ostringstream errorStream;
        errorStream << " Trying to get not existing face"
        << kFaceId << " " << mesh.numBFaces();
        ASSERT( faceContainerIterator != mesh.faceList.end(), errorStream.str().c_str() );

        for ( ID jEdgeLocalId = 1; jEdgeLocalId <= face_Type::numEdges; ++jEdgeLocalId )
        {
            point1Id = facetShape.eToP( jEdgeLocalId, 1 );
            point2Id = facetShape.eToP( jEdgeLocalId, 2 );
            // go to global
            point1Id = ( faceContainerIterator->point( point1Id ) ).id();
            point2Id = ( faceContainerIterator->point( point2Id ) ).id();
            bareEdge = ( makeBareEdge( point1Id, point2Id ) ).first;

            if ( ( edgeContainerIterator = localTemporaryEdgeContainer.find( bareEdge ) )
                    == localTemporaryEdgeContainer.end() )
            {
                localTemporaryEdgeContainer.insert( bareEdge );
                ++numBoundaryEdges;
            }
            else
            {
                localTemporaryEdgeContainer.erase( edgeContainerIterator );
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
bool __attribute__ (( deprecated ))
checkMarkerSet( const MeshEntityListType & meshEntityList )
{
	return checkIsMarkerSetInEntityList( meshEntityList );
}
template <typename MeshEntityListType>
bool checkIsMarkerSetInEntityList( const MeshEntityListType & meshEntityList )
{
    typedef typename MeshEntityListType::const_iterator MeshEntityListTypeConstIterator_Type;
    bool ok( true );
    for ( MeshEntityListTypeConstIterator_Type meshEntityListIterator = meshEntityList.begin();
            meshEntityListIterator != meshEntityList.end(); ++meshEntityListIterator )
        ok = ( ok & meshEntityListIterator->isMarkerSet() );
    return ok;
}


//! Sets the marker flag for all boundary edges by inheriting them from boundary points.
/*!
    The paradigm is that an edge <B>WHOSE MARKER HAS NOT ALREADY BEEN
    SET</B> will get the WEAKER marker flag among its VERTICES. For instance
    if a vertex is assigned to an Essential BC and the other to a Natural
    BC the edge will get the flag related to the Natural B.C.

    @param mesh A mesh
    @param logStream stream to which a map edgeId -> NewMarker will be output
    @param errorStream stream to which error messages will be sent
    @param verbose if false, no messages will be sent to the logStream

    @todo better handling of flags: all function handling flags should be
    wrapped into a class
    @todo errorStream is unused
*/
template <typename MeshType>
void __attribute__ (( deprecated ))
setBEdgesMarker( MeshType & mesh, std::ostream & logStream = std::cout,
                 std::ostream & errorStream = std::cerr, bool verbose = true )
{
	setBoundaryEdgesMarker( mesh, logStream, errorStream, verbose );

}
template <typename MeshType>
void
setBoundaryEdgesMarker( MeshType & mesh, std::ostream & logStream = std::cout,
                 std::ostream & /*errorStream*/ = std::cerr, bool verbose = true )
{
    typename MeshType::EdgeType * edgePtr = 0;
    UInt                  counter( 0 );

    if ( verbose )
        logStream << "NEW EDGE MARKER MAP" << std::endl
        << " ID->New Marker" << std::endl;

    for ( ID kEdgeId = 1; kEdgeId <= mesh.numBEdges(); ++kEdgeId )
    {
        edgePtr = &( mesh.edge( kEdgeId ) );
        if ( edgePtr->isMarkerUnset() )
        {
        	inheritPointsWeakerMarker( *edgePtr );
            if ( verbose )
            {
                logStream << edgePtr->id() << " -> ";
                edgePtr->printFlag( logStream );
                logStream << " ";
                if ( ++counter % 3 == 0 )
                    logStream << std::endl;
            }
        }
    }
    if ( verbose )
        logStream << std::endl;
}


//! Sets the marker flag for all boundary faces by inheriting them from boundary points.
/*!
    The paradigm is that a face <B>WHOSE MARKER HAS NOT ALREADY BEEN SET</B> will
    get the WEAKER marker flag among its VERTICES. For instance if a vertex
    is assigned to a Natural BC and the others to a Natural BC the face
    will get the flag related to the Natural BC

    @param mesh A mesh
    @param logStream stream to which a map faceId -> NewMarker will be output
    @param errorStream stream to which error messages will be sent
    @param verbose if false, no messages will be sent to the logStream

    @todo better handling of flags: all function handling flags should be
    wrapped into a class
*/
template <typename MeshType>
void __attribute__ (( deprecated ))
setBFacesMarker( MeshType & mesh, std::ostream & logStream = std::cout,
                 std::ostream & errorStream = std::cerr, bool verbose = true )
{
	setBoundaryFacesMarker( mesh, logStream, errorStream, verbose );

}
template <typename MeshType>
void
setBoundaryFacesMarker( MeshType & mesh, std::ostream & logStream = std::cout,
                 std::ostream & /*errorStream*/ = std::cerr, bool verbose = true )
{
    typename MeshType::FaceType * facePtr = 0;
    UInt                  counter( 0 );

    if ( verbose )
    {
        logStream << "**** NEW FACE MARKER MAP **************" << std::endl;
        logStream << " Face ID -> New Marker\tFace ID -> New Marker\tFace ID -> New Marker" << std::endl;
    }

    for ( UInt kFaceId = 1; kFaceId <= mesh.numBFaces(); ++kFaceId )
    {
        facePtr = &( mesh.face( kFaceId ) );
        if ( facePtr->isMarkerUnset() )
        {
        	inheritPointsWeakerMarker( *facePtr );
            if ( verbose )
            {
                logStream << facePtr->id() << " -> ";
                facePtr->printFlag( logStream );
                logStream << "\t";
                if ( ++counter % 3 == 0 )
                    logStream << std::endl;
            }
        }
    }
    if ( verbose )
        logStream << std::endl;
}


//! It sets the marker flag of boundary points, by inheriting it from facets.
/*!
    The paradigm is that a point whose marker flag is unset will inherit
    the strongest marker flag of the surrounding facets, with the
    convention that if the marker flag of one of the surrounding facets is null,
    it is ignored.

    @param mesh A mesh
    @param logStream stream to which a map edgeId -> NewMarker will be output
    @param errorStream stream to which error messages will be sent
    @param verbose if false, no messages will be sent to the logStream

*/
template <typename MeshType>
void __attribute__ (( deprecated ))
setBPointsMarker( MeshType & mesh, std::ostream & logStream = std::cout,
                 std::ostream & errorStream = std::cerr, bool verbose = false )
{
	setBoundaryPointsMarker( mesh, logStream, errorStream, verbose );

}
template <typename MeshType>
void
setBoundaryPointsMarker( MeshType & mesh, std::ostream & logStream = std::cout,
                  std::ostream& /*errorStream*/ = std::cerr, bool verbose = false )
{
    // First looks at points whose marker has already been set
    std::vector<bool> isDefinedPointMarker( mesh.storedPoints(), false );

    typedef typename MeshType::Points::iterator pointContainerIterator_Type;
    typedef typename MeshType::BElementShape    facetShape_Type;

    std::vector<bool>::iterator isDefinedPointMarkerIterator = isDefinedPointMarker.begin();

    for ( pointContainerIterator_Type pointContainerIterator = mesh.pointList.begin();
            pointContainerIterator != mesh.pointList.end(); ++pointContainerIterator )
        *( isDefinedPointMarkerIterator++ ) = pointContainerIterator->isMarkerSet();

    typename MeshType::BElementType * facetPtr = 0;
    for ( UInt kFacetId = 1; kFacetId <= mesh.numBElements(); ++kFacetId )
    {
        facetPtr = &( mesh.bElement( kFacetId ) );
        if ( facetPtr->isMarkerSet() )
        {
            for ( UInt jPointId = 1; jPointId <= facetShape_Type::numPoints; ++jPointId )
            {
                if ( !isDefinedPointMarker[ ( facetPtr->point( jPointId ).id() ) - 1 ] )
                    facetPtr->setStrongerMarkerAtPoint( jPointId, facetPtr->marker() );
            }
        }
    }
    UInt counter( 0 );

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
                logStream << pointContainerIterator->id() << " -> ";
                pointContainerIterator->printFlag( logStream );
                logStream << "\t";
                if ( ++counter % 3 )
                    logStream << std::endl;
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
/* More precisely, the id() must correspond to the position of the entity
   in the list (starting from 1, since id=0 is reserved for unset entities.

   @pre The template argument MeshEntityListType must be a stl
   compliant container and its elements must have the method id().
*/
template <typename MeshEntityListType>
bool __attribute__ (( deprecated ))
checkIdnumber( const MeshEntityListType & meshEntityList )
{
	return checkIdInEntityList( meshEntityList );
}
template <typename MeshEntityListType>
bool checkIdInEntityList( const MeshEntityListType & meshEntityList )
{
    typedef typename MeshEntityListType::const_iterator MeshEntityListTypeConstIterator_Type;
    bool ok( true );
    UInt counter( 1 );
    for ( MeshEntityListTypeConstIterator_Type meshEntityListIterator = meshEntityList.begin();
    		meshEntityListIterator != meshEntityList.end() && ok; ++meshEntityListIterator, ++counter )
        ok = ( meshEntityListIterator->id() == counter );
    return ok;
}


//! @brief Fixes a a list of mesh entities so that the ID is properly set.
/* @post  The id will correspond to the position of the entity
   in the list (starting from 1, since id=0 is reserved for unset entities.

   @pre The template argument MeshEntityListType must be a stl
   compliant container and its elements must have the method UInt &id().
*/
template <typename MeshEntityListType>
void __attribute__ (( deprecated ))
fixIdnumber( MeshEntityListType & meshEntityList )
{
	fixIdInEntityList( meshEntityList );
}
template <typename MeshEntityListType>
void fixIdInEntityList( MeshEntityListType & meshEntityList )
{
    UInt counter( 0 );
    typedef typename MeshEntityListType::iterator Iter;
    for ( Iter meshEntityListIterator = meshEntityList.begin();
    		meshEntityListIterator != meshEntityList.end(); ++meshEntityListIterator )
        meshEntityListIterator->setId( ++counter );
}


/*! @brief Fixes boundary points counter
  It fix the boundary points counter by counting
  how many points have the boundary flag set.
  It also resets the boundary points list.

  @pre It assumes that the points have the boundary flag correctly set
*/

template <typename MeshType>
void __attribute__ (( deprecated ))
setBPointsCounters( MeshType & mesh )
{
	setBoundaryPointsCounters( mesh );
}
template <typename MeshType>
void
setBoundaryPointsCounters( MeshType & mesh )
{

    UInt boundaryPointCounter( 0 );
    UInt boundaryVertexCounter( 0 );

    mesh._bPoints.clear();

    for ( UInt kPointId = 1; kPointId <= mesh.numVertices(); ++kPointId )
    {
        if ( mesh.isBoundaryPoint( kPointId ) )
        {
            ++boundaryPointCounter;
            ++boundaryVertexCounter;
        }
    }

    for ( UInt kPointId = mesh.numVertices() + 1; kPointId <= mesh.storedPoints(); ++kPointId )
    {
        if ( mesh.isBoundaryPoint( kPointId ) )
        {
            ++boundaryPointCounter;
        }
    }

    mesh.numBVertices() = boundaryVertexCounter;
    mesh.setNumBPoints( boundaryPointCounter );
    mesh._bPoints.clear();
    mesh._bPoints.reserve( boundaryPointCounter );

    for ( UInt kPointId = 1; kPointId <= mesh.storedPoints(); ++kPointId )
    {
        if ( mesh.isBoundaryPoint( kPointId ) )
            mesh._bPoints.push_back( &mesh.point( kPointId ) );
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
void __attribute__ (( deprecated ))
fixBPoints( MeshType & mesh, std::ostream & logStream = std::cout,
            std::ostream & errorStream = std::cerr, bool verbose = true )
{
	fixBoundaryPoints( mesh, logStream, errorStream, verbose );
}
template <typename MeshType>
void
fixBoundaryPoints( MeshType & mesh, std::ostream & logStream = std::cout,
            std::ostream & /* errorStream */ = std::cerr, bool verbose = true )
{
    ASSERT_PRE( mesh.numPoints() > 0, "The point list should not be empty" );
    ASSERT_PRE( mesh.numBElements() > 0,
                "The BElements list should not be empty" );

    typedef typename MeshType::BElements     facetContainer_Type;
    typedef typename MeshType::BElementShape facetShape_Type;

    if ( verbose ) logStream << "Fixing BPoints" << std::endl;
    std::vector<bool> boundaryPoints(mesh.numPoints());
    // I may have launched the program for a P2 mesh
    // yet not all the points are there
    UInt numitems;
    if (mesh.storedPoints()==mesh.numVertices())
    {
        numitems=facetShape_Type::numVertices;
    }
    else
    {
        numitems=facetShape_Type::numPoints;
    }

    for ( UInt kFacetId = 1; kFacetId <= mesh.numBElements(); ++kFacetId )
        for ( UInt jPointId = 1; jPointId <= numitems; ++jPointId )
            boundaryPoints[mesh.bElement(kFacetId).point(jPointId).id()-1]=true;
    for (ID  kPointId = 1; kPointId <= mesh.storedPoints() ; ++kPointId )
        mesh.point(kPointId).setBoundary(boundaryPoints[kPointId-1]);
    boundaryPoints.clear();
    std::vector<bool> temp;
    boundaryPoints.swap(temp);
    // Fix now the number of vertices/points
    setBoundaryPointsCounters( mesh );
}


//!It makes sure that boundary edges are stored first
/*!
    Calls fixIdInEntityList (@sa fixIdInEntityList)
    @pre It assumes that boundary points are properly stored in the mesh
*/
template <typename MeshType>
void __attribute__ (( deprecated ))
setBoundaryEdgesFirst( MeshType & mesh )
{
	correctEdgesStoringOrder( mesh );
}
template <typename MeshType>
void
correctEdgesStoringOrder( MeshType & mesh )
{

    typedef typename MeshType::Edges Edges;
    // set the functor
    EnquireBEntity<MeshType > enquireBoundaryEdge( mesh );

    std::partition( mesh.edgeList.begin(), mesh.edgeList.end(), enquireBoundaryEdge );
    fixIdInEntityList( mesh.edgeList );
}


//!It makes sure that boundary faces are stored first
/*!
    Calls fixIdInEntityList (@sa fixIdInEntityList)
    @pre It assumes that boundary points are properly stored in the mesh
*/
template <typename MeshType>
void __attribute__ (( deprecated ))
setBoundaryFacesFirst( MeshType & mesh )
{
	correctFacesStoringOrder( mesh );
}
template <typename MeshType>
void
correctFacesStoringOrder( MeshType & mesh )
{

    typedef typename MeshType::Faces faceContainer_Type;
    // set the functor
    EnquireBEntity<MeshType> enquireBoundaryFace( mesh );

    std::partition( mesh.faceList.begin(), mesh.faceList.end(), enquireBoundaryFace );
    fixIdInEntityList( mesh.faceList );
}


//! Tests if boundary faces are stored first
/*! @return true if boundary faces are indeed stored first
  @pre It assumes that boundary points are set */
template <typename MeshType>
bool __attribute__ (( deprecated ))
checkBoundaryFacesFirst( const MeshType & mesh )
{
	return checkFacesStoringOrder( mesh );
}
template <typename MeshType>
bool checkFacesStoringOrder( const MeshType & mesh )
{

    typedef typename MeshType::Faces faceContainer_Type;

    // set the functor
    EnquireBEntity<MeshType> enquireBoundaryFace( mesh );
    bool ok( true );

    for ( UInt kFacetId = 1; kFacetId <= mesh.numFacets(); ++kFacetId )
        ok = ok && enquireBoundaryFace( mesh.boundaryFace( kFacetId ) );
    for ( UInt kFacetId = mesh.numBElements() + 1; kFacetId <= mesh.storedFaces(); ++kFacetId )
        ok = ok && ! enquireBoundaryFace( mesh.face( kFacetId ) );

    return ok;
}


//! Tests if boundary edges are stored first
/*! @return true if boundary edges are indeed stored first
  @pre It assumes that boundary points are set */
template <typename MeshType>
bool __attribute__ (( deprecated ))
checkBoundaryEdgesFirst( const MeshType & mesh )
{
	return checkEdgesStoringOrder( mesh );
}
template <typename MeshType>
bool checkEdgesStoringOrder( const MeshType & mesh )
{

    typedef typename MeshType::Edges Edges;

    // set the functor
    EnquireBEntity<MeshType> enquireBoundaryEdge( mesh );
    bool ok( true );

    for ( UInt kBEdgeId = 1; kBEdgeId <= mesh.numBEdges(); ++kBEdgeId )
        ok = ok && enquireBoundaryEdge( mesh.bareEdge( kBEdgeId ) );
    for ( UInt kBEdgeId = mesh.numBEdges() + 1; kBEdgeId <= mesh.storedEdges(); ++kBEdgeId )
        ok = ok && ! enquireBoundaryEdge( mesh.edge( kBEdgeId ) );
    return ok;
}


/*
*******************************************************************************
 UTILITIES TO VERIFY/CREATE FACES/EDGES
*******************************************************************************
*/
//! It fixes boundary faces so that they are consistently numbered with volumes.

/*! An important step for building degrees of freedom on faces.  It also
    fixes other face related data.
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
	<li>BELEMENT_COUNTER_UNSET</li>
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
*/


#ifndef TWODIM
template <class MeshType>
bool fixBoundaryFaces( MeshType & mesh,
                       std::ostream & logStream,
                       std::ostream &errorStream,
                       Switch & sw,
                       UInt & numFaces,
                       UInt & numBoundaryFaces,
                       bool /* fixMarker */ = false,
                       bool verbose = false,
                       temporaryFaceContainer_Type * externalFaceContainer = 0 )
{

    typedef typename MeshType::Volumes volumeContainer_Type;
    typedef typename MeshType::VolumeType VolumeType;
    typedef typename MeshType::Faces faceContainer_Type;
    typedef typename MeshType::FaceType face_Type;

    UInt                                  point1Id, point2Id, point3Id, point4Id;
    BareFace                              bareFace;
    VolumeType *                          volumePtr;
    typename faceContainer_Type::iterator faceContainerIterator;
    typename MeshType::VolumeShape        volumeShape;
    temporaryFaceContainer_Type *         boundaryFaceContainerPtr;
    temporaryFaceContainer_Type::iterator boundaryFaceContainerIterator;
    std::pair<ID, ID>                     volumeIdToLocalFaceIdPair;
    ID                                    jFaceLocalId;
    ID                                    volumeId;
    UInt                                  numInternalFaces;
    bool                                  notfound( false );
    bool                                  externalContainerIsProvided( false );

    if ( (externalContainerIsProvided = ( externalFaceContainer != 0 )) )
    {
        boundaryFaceContainerPtr = externalFaceContainer;
        numBoundaryFaces = boundaryFaceContainerPtr->size();
    }
    else
    {
        boundaryFaceContainerPtr = new temporaryFaceContainer_Type;
        numBoundaryFaces = findBoundaryFaces( mesh, *boundaryFaceContainerPtr, numInternalFaces );
        numFaces = numBoundaryFaces + numInternalFaces;
    }


    bool notEnough = mesh.storedFaces() < numBoundaryFaces;



    if ( notEnough )
    {
        errorStream << "WARNING: number of B. Faces stored smaller" << std::endl;
        errorStream << "         than the number of boundaryFaces found  and build is not set"
        << std::endl;
        errorStream << "POSSIBLE ERROR" << std::endl;
        sw.create( "BFACE_STORED_MISMATCH", true );
    }

    if ( mesh.numBElements() == 0 )
    {
        errorStream << "ERROR: Boundary Element counter was not set" << std::endl;
        errorStream << "I Cannot proceed because the situation is ambiguous"
        << std::endl;
        errorStream << "Please check and eventually either: (a) call buildBoundaryFaces()" << std::endl;
        errorStream << "or (b) set the correct number of boundaryFaces in the mesh using mesh.numBElements()" << std::endl;
        errorStream << "ABORT" << std::endl;
        sw.create( "BELEMENT_COUNTER_UNSET", true );
    }

    if ( mesh.numBFaces() != numBoundaryFaces )
    {
        errorStream << "WARNING: Boundary face counter in mesh is set to "
        << mesh.numBFaces() << std::endl;
        errorStream << "         while I have found " << numBoundaryFaces
        << " boundary elements in mesh." << std::endl;
        errorStream << "         Please check... I continue anyway" << std::endl;
        sw.create( "BFACE_COUNTER_MISMATCH", true );
    }

    if ( verbose )
    {
        logStream << "**** Fixed Marker Flags for Boundary Faces ***" << std::endl;
        logStream << " (it only contains those that were fixed because unset !)"
        << std::endl;
        logStream << "id->marker   id->marker  id->marker" << std::endl;
    }

    UInt counter( 0 );

    faceContainerIterator = mesh.faceList.begin();
    for ( UInt facid = 0; facid < mesh.numBElements(); ++facid )
    {
        point1Id = ( faceContainerIterator->point( 1 ) ).id();
        point2Id = ( faceContainerIterator->point( 2 ) ).id();
        point3Id = ( faceContainerIterator->point( 3 ) ).id();
        if ( MeshType::FaceShape::numVertices == 4 )
        {
            point4Id = ( faceContainerIterator->point( 4 ) ).id();
            bareFace = ( makeBareFace( point1Id, point2Id, point3Id, point4Id ) ).first;
        }
        else
        {
            bareFace = ( makeBareFace( point1Id, point2Id, point3Id ) ).first;
        }
        boundaryFaceContainerIterator = boundaryFaceContainerPtr->find( bareFace );
        if ( boundaryFaceContainerIterator == boundaryFaceContainerPtr->end() )
        {
            if (verbose)
            {
                if ( MeshType::FaceShape::numVertices == 3 )
                {
                    errorStream<<"Face "<<point1Id<<" "<<point2Id<<" "<<point3Id;
                }
                else
                {
                    errorStream<<"Face "<<point1Id<<" "<<point2Id<<" "<<point3Id<<" " <<point4Id;
                }
                errorStream<<" stored as boundary face, it's not!"<< std::endl;
            }
            notfound = true;
        }
        else
        {
            volumeIdToLocalFaceIdPair = boundaryFaceContainerIterator->second;
            volumeId = volumeIdToLocalFaceIdPair.first; // Element ID
            volumePtr = &mesh.volume( volumeId ); // Element
            jFaceLocalId = volumeIdToLocalFaceIdPair.second;       // The local ID of face on element
            // Reset face point definition to be consistent with face.
            for ( UInt kPointId = 1; kPointId <= face_Type::numPoints; ++kPointId )
            {
                faceContainerIterator->setPoint( kPointId, volumePtr->point( volumeShape.fToP( jFaceLocalId, kPointId ) ) );
            }
            // Correct extra info
            faceContainerIterator->ad_first() = volumeId;
            faceContainerIterator->pos_first() = jFaceLocalId;
            if ( faceContainerIterator->isMarkerUnset() )
            {
                inheritPointsWeakerMarker( *faceContainerIterator );
                if ( verbose )
                {
                    logStream << faceContainerIterator->id() << " -> ";
                    faceContainerIterator->printFlag( logStream );
                    logStream << " ";
                    if ( ++counter % 3 == 0 )
                        logStream << std::endl;
                }
            }
            // Take out face from temporary container
            boundaryFaceContainerPtr->erase( boundaryFaceContainerIterator );
        }
        ++faceContainerIterator;
    }

    if ( !externalContainerIsProvided )
        delete boundaryFaceContainerPtr;

    if ( notfound )
    {
        errorStream << "WARNING: At least one boundary face has not been found on the list stored in MeshType\n";
        sw.create( "BFACE_MISSING", true );
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
        sw.create( "NUM_FACES_MISMATCH", true );
    }
    mesh.setLinkSwitch( std::string( "HAS_BOUNDARY_FACES" ) );

    return true;
}
#endif


//! Builds faces
/*! This function may alternatively be used to build the compulsory boundary
  faces, all the mesh faces, or just add to an existing list of just boundary
  faces the internal ones.

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
bool buildFaces( MeshType & mesh,
                 std::ostream & logStream,
                 std::ostream & errorStream,
                 UInt & numBoundaryFaces,
                 UInt & numInternalFaces,
                 bool buildBoundaryFaces = true,
                 bool buildInternalFaces = false,
                 bool verbose = false,
                 temporaryFaceContainer_Type * externalFaceContainer = 0 )
{
    UInt                                  point1Id, point2Id, point3Id, point4Id;
    typename MeshType::VolumeShape        volumeShape;
    typedef typename MeshType::Volumes    volumeContainer_Type;
    typedef typename MeshType::VolumeType volume_Type;
    typedef typename MeshType::Faces      faceContainer_Type;
    typedef typename MeshType::FaceType   face_Type;
    volume_Type *                         volumePtr;
    temporaryFaceContainer_Type*          boundaryFaceContainerPtr;
    temporaryFaceContainer_Type::iterator boundaryFaceContainerIterator;
    bool                                  externalContainerIsProvided( false );

    std::pair<ID, ID>                     volumeIdToLocalFaceIdPair;
    ID                                    jFaceLocalId, newFaceId;
    ID                                    volumeId;

    if ( (externalContainerIsProvided = ( externalFaceContainer != 0 )) )
    {
        boundaryFaceContainerPtr = externalFaceContainer;
        numBoundaryFaces = boundaryFaceContainerPtr->size();
    }
    else
    {
        boundaryFaceContainerPtr = new temporaryFaceContainer_Type;
        numBoundaryFaces = findBoundaryFaces( mesh, *boundaryFaceContainerPtr, numInternalFaces );
    }

    if ( buildBoundaryFaces )
        mesh.faceList.clear();
    mesh.setNumBFaces( numBoundaryFaces );
    if ( !buildInternalFaces )
    {
        mesh.setMaxNumFaces( numBoundaryFaces, false );
        mesh.setNumFaces( numInternalFaces + numBoundaryFaces );
    }
    else
    {
        mesh.setMaxNumFaces( numInternalFaces + numBoundaryFaces, true );
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
            volumeIdToLocalFaceIdPair = boundaryFaceContainerIterator->second;
            volumeId = volumeIdToLocalFaceIdPair.first; // Element ID
            volumePtr = &mesh.volume( volumeId ); // Element
            jFaceLocalId = volumeIdToLocalFaceIdPair.second;       // The local ID of face on element

            for ( UInt kPointId = 1; kPointId <= face_Type::numPoints; ++kPointId )
                face.setPoint( kPointId, volumePtr->point( volumeShape.fToP( jFaceLocalId, kPointId ) ) );
            // Add extra info
            face.ad_first() = volumeId;
            face.pos_first() = jFaceLocalId;
            // Get marker value
            inheritPointsWeakerMarker( face );
            newFaceId = mesh.addFace( face, true ).id();
            if ( verbose )
            {
                if ( newFaceId % 3 == 0 )
                    logStream << std::endl;
                logStream << newFaceId << " -> ";
                face.printFlag( logStream );
                logStream << " ";
            }
        }
        mesh.setLinkSwitch( std::string( "HAS_BOUNDARY_FACES" ) );
        if ( ! buildInternalFaces )
            mesh.unsetLinkSwitch( std::string( "HAS_ALL_FACES" ) );
        mesh.setLinkSwitch( std::string( "FACES_HAVE_ADIACENCY" ) );
    }

    if ( !externalContainerIsProvided )
        delete boundaryFaceContainerPtr;

    if ( ! buildInternalFaces )
        return true;


    if ( !buildBoundaryFaces )
    {
        if ( mesh.storedFaces() < mesh.numBFaces() )
        {
            errorStream << "ERROR: mesh has not boundary faces, cannot just create internal ones!!!" << std::endl;
            errorStream << "ABORT CONDITION" << std::endl;
            return false;
        }
        else if ( mesh.storedFaces() > mesh.numBFaces() )
        {
            mesh.faceList.resize( mesh.numBFaces() );
        }
    }


    /*
      I may get rid of the boundaryFaces container. Unfortunately now I need a more
      complex structure, a BareItemsHandler, in order to generate the internal
      faces id. An alternative would be to use the point data to identify
      boundary faces as the ones with all point on the boundary. Yet in this
      function we do not want to use a priori information, so that it might
      work even if the points boundary flag is not properly set.
    */

    BareItemsHandler<BareFace> bareFaceHandler;
    std::pair<UInt, bool> faceIdToBoolPair;
    std::pair<BareFace, bool> _face;

    for ( UInt jFaceId = 0; jFaceId < mesh.faceList.size(); ++jFaceId )
    {
        point1Id = ( mesh.faceList[ jFaceId ].point( 1 ) ).id();
        point2Id = ( mesh.faceList[ jFaceId ].point( 2 ) ).id();
        point3Id = ( mesh.faceList[ jFaceId ].point( 3 ) ).id();
        if ( MeshType::FaceShape::numVertices == 4 )
        {
            point4Id = ( mesh.faceList[ jFaceId ].point( 4 ) ).id();
            _face = makeBareFace( point1Id, point2Id, point3Id, point4Id );
        }
        else
        {
            _face = makeBareFace( point1Id, point2Id, point3Id );
        }
        bareFaceHandler.addIfNotThere( _face.first );
    }

    EntityFlag meshMarker( mesh.marker() );
    // ID volumeId;

    for ( typename volumeContainer_Type::iterator volumeContainerIterator = mesh.volumeList.begin();
            volumeContainerIterator != mesh.volumeList.end(); ++volumeContainerIterator )
    {
        volumeId = volumeContainerIterator->id();
        // REMEMBER: numbering from 1
        for ( UInt jFaceLocalId = 1; jFaceLocalId <= mesh.numLocalFaces(); jFaceLocalId++ )
        {
            point1Id = volumeShape.fToP( jFaceLocalId, 1 );
            point2Id = volumeShape.fToP( jFaceLocalId, 2 );
            point3Id = volumeShape.fToP( jFaceLocalId, 3 );
            // go to global
            point1Id = ( volumeContainerIterator->point( point1Id ) ).id();
            point2Id = ( volumeContainerIterator->point( point2Id ) ).id();
            point3Id = ( volumeContainerIterator->point( point3Id ) ).id();
            if ( MeshType::FaceShape::numVertices == 4 )
            {
                point4Id = volumeShape.fToP( jFaceLocalId, 4 );
                point4Id = ( volumeContainerIterator->point( point4Id ) ).id();
                _face = makeBareFace( point1Id, point2Id, point3Id, point4Id );
            }
            else
            {
                _face = makeBareFace( point1Id, point2Id, point3Id );
            }
            faceIdToBoolPair = bareFaceHandler.addIfNotThere( _face.first );
            if ( faceIdToBoolPair.second )
            {
                // a new face It must be internal.
                for ( UInt kPointId = 1; kPointId <= face_Type::numPoints; ++kPointId )
                    face.setPoint( kPointId, volumeContainerIterator->point( volumeShape.fToP( jFaceLocalId, kPointId ) ) );
                face.ad_first() = volumeId;
                face.pos_first() = jFaceLocalId;
                // gets the marker from the MeshType
                face.setMarker( meshMarker );
                mesh.addFace( face, false ); //The id should be correct
            }
            else
            {
                if ( faceIdToBoolPair.first > numBoundaryFaces )  // internal
                {
                    mesh.faceList( faceIdToBoolPair.first ).ad_second() = volumeId;
                    mesh.faceList( faceIdToBoolPair.first ).pos_second() = jFaceLocalId;
                }
            }
        }
    }
    mesh.setLinkSwitch( std::string( "HAS_ALL_FACES" ) );
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
bool buildEdges( MeshType & mesh,
                 std::ostream & logStream,
                 std::ostream & errorStream,
                 UInt & numBoundaryEdgesFound,
                 UInt & numInternalEdgesFound,
                 bool buildBoundaryEdges = true,
                 bool buildInternalEdges = false,
                 bool verbose = false,
                 temporaryEdgeContainer_Type * externalEdgeContainer = 0 )
{
    typedef typename MeshType::Volumes volumeContainer_Type;
    typedef typename MeshType::Faces faceContainer_Type;
    typedef typename MeshType::VolumeType volume_Type;
    typedef typename MeshType::VolumeShape volumeShape_Type;
    typedef typename MeshType::Edges Edges;
    typedef typename MeshType::EdgeType edge_Type;
    typedef typename MeshType::FaceType face_Type;
    typedef typename MeshType::FaceShape faceShape_Type;

    face_Type * facePtr;
    volume_Type * volumePtr;

    temporaryEdgeContainer_Type * temporaryEdgeContainer;
    temporaryEdgeContainer_Type edgeContainer;
    std::pair<ID, ID> faceIdToLocalEdgeIdPair;
    ID jEdgeLocalId, newEdgeId;
    ID faceId;


    bool externalContainerIsProvided( false );


    if ( (externalContainerIsProvided = ( externalEdgeContainer != 0 )) )
    {
        temporaryEdgeContainer = externalEdgeContainer;
        numBoundaryEdgesFound = temporaryEdgeContainer->size();
    }
    else
    {
        temporaryEdgeContainer = new temporaryEdgeContainer_Type;
        numBoundaryEdgesFound = findBoundaryEdges( mesh, *temporaryEdgeContainer );
    }

    numInternalEdgesFound = findInternalEdges( mesh, *temporaryEdgeContainer, edgeContainer );
    // free some memory if not needed!
    if ( !buildInternalEdges )
        edgeContainer.clear();
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
            mesh.edgeList.resize( numBoundaryEdgesFound );
        }
    }
    mesh.setNumBEdges( numBoundaryEdgesFound );
    mesh.setNumEdges( numBoundaryEdgesFound + numInternalEdgesFound );
    if ( buildBoundaryEdges )
        mesh.edgeList.clear();
    if ( buildBoundaryEdges && ! buildInternalEdges )
        mesh.setMaxNumEdges( numBoundaryEdgesFound, false );
    if ( buildInternalEdges )
        mesh.setMaxNumEdges( numBoundaryEdgesFound + numInternalEdgesFound, true );

    if (verbose)
        errorStream << "Building edges from scratch" << std::endl;

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
            facePtr = &mesh.face( faceId ); // Face
            jEdgeLocalId = faceIdToLocalEdgeIdPair.second;       // The local ID of edge on face
            for ( UInt kPointId = 1; kPointId <= edge_Type::numPoints; ++kPointId )
            {
                edge.setPoint( kPointId, facePtr->point( faceShape_Type::eToP( jEdgeLocalId, kPointId ) ) );
            }

            // Get marker value inheriting from points
            inheritPointsWeakerMarker( edge );

            newEdgeId = mesh.addEdge( edge, true ).id();
            if ( verbose )
            {
                if ( newEdgeId % 6 == 0 )
                    logStream << std::endl;
                logStream << newEdgeId << " -> ";
                edge.printFlag( logStream );
                logStream << " ";
            }
        }

        if ( verbose )
            logStream << std::endl << "  *****  END OF LIST OF BOUNDARY EDGES ****"
            << std::endl;

        mesh.setLinkSwitch( std::string( "HAS_BOUNDARY_EDGES" ) );
    }

    if ( !externalContainerIsProvided )
        delete temporaryEdgeContainer;

    if ( !buildInternalEdges )
    {
        mesh.unsetLinkSwitch( std::string( "HAS_ALL_EDGES" ) );
        return true;
    }



    // Now internal edges
    // free some memory

    for ( temporaryEdgeContainer_Type::iterator edgeContainerIterator = edgeContainer.begin();
            edgeContainerIterator != edgeContainer.end(); ++edgeContainerIterator )
    {
        faceIdToLocalEdgeIdPair = edgeContainerIterator->second;
        faceId = faceIdToLocalEdgeIdPair.first; // Volume ID
        volumePtr = &mesh.volume( faceId ); // Volume that generated the edge
        jEdgeLocalId = faceIdToLocalEdgeIdPair.second;       // The local ID of edge on volume
        for ( UInt kPointId = 1; kPointId <= edge_Type::numPoints; ++kPointId )
            edge.setPoint( kPointId, volumePtr->point( volumeShape_Type::eToP( jEdgeLocalId, kPointId ) ) );
        edge.setMarker( mesh.marker() ); // Get marker value: that of the mesh
        mesh.addEdge( edge, false );
    }

    mesh.setLinkSwitch( std::string( "HAS_ALL_EDGES" ) );

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
	@note the function takes advantage of the fact that
    @param mesh[out] A mesh
	@param logStream[out] Log stream for information on the newly created markers for boundary edges
*/
template <typename MeshType>
void __attribute__ (( deprecated ))
p1top2( MeshType & mesh, std::ostream & logStream = std::cout )
{
	p2MeshFromP1Data( mesh, logStream );
}
template <typename MeshType>
void
p2MeshFromP1Data( MeshType & mesh, std::ostream & logStream = std::cout )
{

    typedef typename MeshType::ElementShape  GeoShape;
    typedef typename MeshType::BElementShape GeoBShape;
    ASSERT_PRE( GeoShape::numPoints > 4, "p2MeshFromP1Data ERROR: we need a P2 mesh" );

    logStream << "Building P2 mesh points and connectivities from P1 data"
    << std::endl;


    typename MeshType::PointType *       pointPtr = 0;
    typename MeshType::EdgeType *        edgePtr = 0;
    typename MeshType::ElementType *     elementPtr = 0;
    typename MeshType::BElementType *    facetPtr = 0;
    typedef typename MeshType::Elements  elementContainer_Type;
    typedef typename MeshType::BElements facetContainer_Type;

    BareItemsHandler<BareEdge>           bareEdgeHandler;
    std::pair<UInt, bool>                edgeIdToBoolPair;
    UInt                                 point1Id, point2Id, edgeId;
    std::pair<BareEdge, bool>            bareEdgeToBoolPair;
    typename MeshType::ElementShape      elementShape;

    logStream << "Processing " << mesh.storedEdges() << " P1 Edges" << std::endl;
    UInt numBoundaryEdges = mesh.numBEdges();
    for ( UInt jEdgeId = 1; jEdgeId <= mesh.storedEdges(); ++jEdgeId )
    {
        edgePtr = & mesh.edge( jEdgeId );
        point1Id = ( edgePtr->point( 1 ) ).id();
        point2Id = ( edgePtr->point( 2 ) ).id();
        pointPtr = & mesh.addPoint( jEdgeId <= numBoundaryEdges ); // true for boundary points
        pointPtr->x() = ( ( edgePtr->point( 1 ) ).x() +
                    ( edgePtr->point( 2 ) ).x() ) * .5;
        pointPtr->y() = ( ( edgePtr->point( 1 ) ).y() +
                    ( edgePtr->point( 2 ) ).y() ) * .5;
        pointPtr->z() = ( ( edgePtr->point( 1 ) ).z() +
                    ( edgePtr->point( 2 ) ).z() ) * .5;

        /*
          If we have set a marker for the boundary edge, that marker is
          inherited by the new created point. Otherwise the edge (and the new
          created point) gets the WORST marker among the two end Vertices
        */
        if ( edgePtr->isMarkerUnset() )
            inheritPointsWeakerMarker( *edgePtr );
        pointPtr->setMarker( edgePtr->marker() );
        // todo check that the id() of the new point is correctly set
        edgePtr->setPoint( 3, pointPtr ); //use overloaded version that takes a pointer
        bareEdgeToBoolPair = makeBareEdge( point1Id, point2Id );
        edgeIdToBoolPair = bareEdgeHandler.addIfNotThere( bareEdgeToBoolPair.first, pointPtr->id() );
    }
    // Now the other edges, of which I do NOT build the global stuff
    // (I would need to check the switch but I will do that part later on)
    if ( GeoShape::nDim == 3 )
    {
        UInt numBoundaryFaces = mesh.numBFaces();

        logStream << "Processing " << mesh.storedFaces() << " Face Edges"
        << std::endl;
        for ( UInt kFaceId = 1; kFaceId <= mesh.storedFaces(); ++kFaceId )
        {
            facetPtr = &mesh.face( kFaceId );
            for ( UInt jEdgeLocalId = 1; jEdgeLocalId <= mesh.numLocalEdgesOfFace(); jEdgeLocalId++ )
            {
                point1Id = GeoBShape::eToP( jEdgeLocalId, 1 );
                point2Id = GeoBShape::eToP( jEdgeLocalId, 2 );
                point1Id = ( facetPtr->point( point1Id ) ).id();
                point2Id = ( facetPtr->point( point2Id ) ).id();
                bareEdgeToBoolPair = makeBareEdge( point1Id, point2Id );
                edgeId = bareEdgeHandler.id( bareEdgeToBoolPair.first );
                if ( edgeId != 0 )
                {
                    pointPtr = &mesh.point( edgeId );
                }
                else
                {
                    // new edge -> new Point
                    pointPtr = &mesh.addPoint( kFaceId <= numBoundaryFaces );// true for boundary points
                    edgeIdToBoolPair = bareEdgeHandler.addIfNotThere( bareEdgeToBoolPair.first, pointPtr->id() );
                    pointPtr->x() = ( mesh.point( point1Id ).x() +
                                mesh.point( point2Id ).x() ) * .5;
                    pointPtr->y() = ( mesh.point( point1Id ).y() +
                                mesh.point( point2Id ).y() ) * .5;
                    pointPtr->z() = ( mesh.point( point1Id ).z() +
                                mesh.point( point2Id ).z() ) * .5;
                    // If we have set a marker for the face, that marker is
                    // inherited by the new created point
                    pointPtr->setMarker( facetPtr->marker() );
                }
                facetPtr->setPoint( GeoBShape::numVertices + jEdgeLocalId, pointPtr );
            }
        }
    }

    logStream << "Processing " << mesh.numElements() << " Mesh Elements"
    << std::endl;
    UInt nev = GeoShape::numVertices;
    for ( UInt kElementId = 1; kElementId <= mesh.numElements(); ++kElementId )
    {
        elementPtr = &mesh.element( kElementId );
        for ( UInt jEdgeLocalId = 1; jEdgeLocalId <= mesh.numLocalEdges(); jEdgeLocalId++ )
        {
            point1Id = elementShape.eToP( jEdgeLocalId, 1 );
            point2Id = elementShape.eToP( jEdgeLocalId, 2 );
            point1Id = ( elementPtr->point( point1Id ) ).id();
            point2Id = ( elementPtr->point( point2Id ) ).id();
            bareEdgeToBoolPair = makeBareEdge( point1Id, point2Id );
            edgeId = bareEdgeHandler.id( bareEdgeToBoolPair.first );
            if ( edgeId != 0 )
            {
                pointPtr = &mesh.point( edgeId );
            }
            else
            {
                // cannot be on boundary if the mesh is proper!
                pointPtr = &mesh.addPoint( false );
                edgeIdToBoolPair = bareEdgeHandler.addIfNotThere( bareEdgeToBoolPair.first, pointPtr->id() );
                pointPtr->x() = ( mesh.point( point1Id ).x() +
                            mesh.point( point2Id ).x() ) * .5;
                pointPtr->y() = ( mesh.point( point1Id ).y() +
                            mesh.point( point2Id ).y() ) * .5;
                pointPtr->z() = ( mesh.point( point1Id ).z() +
                            mesh.point( point2Id ).z() ) * .5;
                pointPtr->setMarker( edgePtr->marker() );
            }
            elementPtr->setPoint( nev + jEdgeLocalId, pointPtr );
        }
    }
    /*=============================*/
    logStream << " ******* Done Construction of P2 Mesh *******"
    << std::endl << std::endl;
}

//! Fix mesh switches
/*!
  Using some heuristics it tries to fix mesh switches
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
//     if(mesh.face(1).ad_first(
//     mesh.setLinkSwitch("HAS_BOUNDARY_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDAY_EDGES");
//   }


// } // namespace MeshUtilities

} // namespace LifeV
#endif
